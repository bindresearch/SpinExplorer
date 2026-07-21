import numpy as np
import nmrglue as ng
from numpy.typing import NDArray
from typing import Union, Optional
import pyfftw
import pyfftw.interfaces.numpy_fft as fft
from SpinExplorer.SpinProcess.Processing.IST.sampling_utils import apply_sampling_schedule_nd, apply_sampling_schedule_to_2D_signal
from joblib import Parallel, delayed


pyfftw.interfaces.cache.enable()
pyfftw.interfaces.cache.set_keepalive_time(30)  # keep plans cached for 30s


def write_as_nmrpipe(array: NDArray, spec_dic: dict, outfile: str):
    ng.pipe.write(outfile, spec_dic, array, overwrite=True)


def get_thresh_signal_3d(signal_real: NDArray,
                        signal_imag: NDArray,
                        threshold: float,
                        window_dim2: NDArray,window_dim3:NDArray):

    def soft_thresh(x, threshold):
        mag = np.abs(x)
        cutoff = threshold * np.max(mag)
        scale = np.maximum(0.0, mag - cutoff) / (mag + 1e-12)
        return x * scale

    thresh_real = soft_thresh(signal_real, threshold)
    thresh_imag = soft_thresh(signal_imag, threshold)

    thresh_fid = retrieve_signal_ist_3d(thresh_real, thresh_imag, window_dim2, window_dim3)

    thresh_real, thresh_imag = pack_signal_ist_3d(thresh_fid)


    return thresh_real, thresh_imag


def ist_iteration_3d(data: NDArray,
                    threshold: float,
                    sampling_schedule: Union[list[int], NDArray],
                    window_dim2:NDArray,
                    window_dim3:NDArray) -> tuple[NDArray, NDArray, NDArray, np.floating]:

    # FT along faster indirect (axis=1): collapse cos/sin pairs → complex, FT, take real

    data_real, data_imag = pack_signal_ist_3d(data)

    thresh_sig_real, thresh_sig_imag = get_thresh_signal_3d(data_real, data_imag, threshold, window_dim2, window_dim3)
    leftover_sig_real = data_real - thresh_sig_real
    leftover_sig_imag = data_imag - thresh_sig_imag
    l2_norm = np.linalg.norm(leftover_sig_real)

    # Inverse: re-expand both axes back to interleaved form
    leftover_sig = retrieve_signal_ist_3d(leftover_sig_real, leftover_sig_imag)  # (n3*2, n2*2)
    
    leftover_sig = apply_sampling_schedule_nd(leftover_sig, sampling_schedule, (1, 0))
  
    return leftover_sig, thresh_sig_real, thresh_sig_imag, l2_norm


def pack_signal_ist_3d(data: NDArray):

    data = data[:, 0::2] + 1.0j * data[:, 1::2]       
    data = fft.fft(data, axis=-1)
    data_real = np.real(data)                                 
    data_imag = np.imag(data)
    
    data_real = np.transpose(data_real)
    data_imag = np.transpose(data_imag)

    data_real = data_real[:,0::2] + 1.0j*data_real[:,1::2]
    data_imag = data_imag[:,0::2] + 1.0j*data_imag[:,1::2]
    
        
    stacked = np.stack([data_real, data_imag])  # (2, ...)
    stacked = fft.fft(stacked, axis=-1)
    data_real, data_imag = stacked[0], stacked[1]
    return data_real, data_imag

def retrieve_signal_ist_3d(leftover_sig_real: NDArray, leftover_sig_imag: NDArray, window_dim2: Optional[NDArray] = None, window_dim3: Optional[NDArray] = None):

    
    
    stacked = np.stack([leftover_sig_real, leftover_sig_imag])
    stacked = fft.ifft(stacked, axis=-1)
    leftover_sig_real, leftover_sig_imag = stacked[0], stacked[1]
    
    if window_dim2 is None:
        window_dim2 = np.ones(leftover_sig_real.shape[-1])

    leftover_sig_real*=window_dim2
    leftover_sig_imag*=window_dim2

    def unpack_complex(x: NDArray):
        out = np.empty(x.shape[:-1] + (x.shape[-1]*2,), dtype=np.float32)
        out[..., 0::2] = x.real
        out[..., 1::2] = x.imag
        return out

    leftover_sig_real = unpack_complex(leftover_sig_real)
    leftover_sig_imag = unpack_complex(leftover_sig_imag)

    
    leftover_sig_real = np.ascontiguousarray(leftover_sig_real.T)
    leftover_sig_imag = np.ascontiguousarray(leftover_sig_imag.T)

    data = leftover_sig_real + 1j*leftover_sig_imag

    data = fft.ifft(data, axis = -1)

    if window_dim3 is None:
        window_dim3 = np.ones(data.shape[-1])

    data *= window_dim3
    data = unpack_complex(data)
    

    return data


def ist_3d(input_spec: NDArray,
           sampling_schedule: Union[list[int], NDArray],
           max_time_dim2: float,
           max_time_dim3: float,
           r2_dim2: Optional[float]= None,
           r2_dim3: Optional[float]= None,
           threshold: float = 0.9,
           terminate: float = 1e-4,
           convergence_tol: float = 1e-8,
           max_iter: int = 800,
           mode: int = 1,
           sched_ord: int = 0,
           verb: bool = False) -> NDArray:
    
    """
    IST reconstruction of 3D NUS data
    """
    if sched_ord == 1:
        sampling_schedule = np.asarray(sampling_schedule)
        sampling_schedule = sampling_schedule[:, ::-1]

    if r2_dim2 is None:
        r2_dim2 = 1.5/(max_time_dim2)
    
    if r2_dim3 is None:
        r2_dim3 = 1.5/(max_time_dim3)

    window_vals_dim3 = np.arange(input_spec.shape[-1]/2)/input_spec.shape[-1]/2
    window_dim3 = np.exp(-window_vals_dim3*r2_dim3*max_time_dim3)


    window_vals_dim2 = np.arange(input_spec.shape[1]/2)/input_spec.shape[1]/2
    window_dim2 = np.exp(-window_vals_dim2*r2_dim2*max_time_dim2)

    def _reconstruct_until_convergence(nus_fid: NDArray) -> tuple[NDArray,NDArray]:
        reconstructed_r = None
        
        reconstructed_i = None
        prev_norm = 0.0


        for iteration in range(1, max_iter + 1):
            nus_fid, threshold_sig_real, threshold_sig_imag, _ = ist_iteration_3d(nus_fid, threshold, sampling_schedule, window_dim2, window_dim3)
            if reconstructed_r is None:
                reconstructed_r = np.zeros_like(threshold_sig_real)
                reconstructed_i = np.zeros_like(threshold_sig_imag)
                
            reconstructed_r += threshold_sig_real 
            reconstructed_i += threshold_sig_imag

            curr_norm = np.sqrt(np.vdot(reconstructed_r, reconstructed_r).real)
            relative_change = abs(curr_norm - prev_norm) / (curr_norm + 1e-10)
            prev_norm = curr_norm

            if relative_change < convergence_tol:
                print(f"  converged at iteration {iteration} — "
                      f"relative change: {relative_change:.2e}")
                break
            if iteration == max_iter:
                print(f"  reached max iterations of {iteration} — "
                      f"relative change: {relative_change:.2e}")

        return reconstructed_r, reconstructed_i

    def _reconstruct_until_l2(nus_fid: NDArray) -> NDArray:
        reconstructed = np.zeros((nus_fid.shape[1]//2,nus_fid.shape[0]//2))
        l2_norm = terminate * 1000.0

        for iteration in range(1, max_iter + 1):
            nus_fid, threshold_signal_real, threshold_signal_imag, l2_norm = ist_iteration_3d(nus_fid, threshold, sampling_schedule)
            reconstructed += (threshold_signal_real+1j*threshold_signal_imag)
            if l2_norm <= terminate:
                if verb:
                    print(f"  converged at iteration {iteration} — "
                        f"l2 norm: {l2_norm:.2e}")
                break
            if iteration == max_iter:
                if verb:
                    print(f"  reached max iterations — L2 norm: {l2_norm:.2e}")

        return reconstructed

    recon_spec = np.zeros_like(input_spec)
    reconstruct = _reconstruct_until_convergence if mode == 1 else _reconstruct_until_l2
    
    def _process_slice_3d(i):
        if verb:
            print(f"Doing IST slice {i}")
        
        slice_data = input_spec[i].copy()
        recon_real, recon_imag = reconstruct(slice_data)
        recon_slice = retrieve_signal_ist_3d(recon_real, recon_imag)
        return i, recon_slice

    results = Parallel(n_jobs = -1)(delayed(_process_slice_3d)(i) for i in range(input_spec.shape[0]))
    for i, result in results:
        recon_spec[i] = result


    return recon_spec



def get_thresh_signal(signal_ft: NDArray,
                        threshold: float,
                        window: NDArray):

    def soft_thresh(x, threshold):
        mag = np.abs(x)
        cutoff = threshold * np.max(mag)
        scale = np.maximum(0.0, mag - cutoff) / (mag + 1e-12)
        return x * scale

    thresh_real = soft_thresh(signal_ft, threshold)
    thresh_fid = fft.ifft(thresh_real)
    thresh_fid*=window

    thresh_real = fft.fft(thresh_fid)

    return thresh_real


def ist_iteration_2d(nus_fid:NDArray, 
                  threshold:float, 
                  sampling_schedule:Union[list[int],NDArray],
                  window: NDArray)->tuple[NDArray,NDArray,np.floating]:
    
    # consider whether we want to put zero-filling into this
    


    signal_ft = fft.fft(nus_fid)
    threshold_sig = get_thresh_signal(signal_ft, threshold, window)
    leftover_sig = signal_ft-threshold_sig
    l2_norm = l2_norm = np.sqrt(np.vdot(leftover_sig, leftover_sig).real)

    leftover_fid = fft.ifft(leftover_sig)
    leftover_fid = apply_sampling_schedule_nd(leftover_fid, sampling_schedule, (0,))

    return leftover_fid, threshold_sig, l2_norm


def ist_2d(input_spec: NDArray,
           sampling_schedule: Union[list[int], NDArray],
           max_time: float,
           r2: Optional[float] = None,
           threshold: float = 0.9,
           terminate: float = 0.001,
           convergence_tol: float = 1e-6,
           max_iter: int = 4000,
           mode: int = 1,
           verb: bool = False) -> NDArray:
    """
    IST reconstruction of 2D NUS data. Direct dimension is assumed to be
    already Fourier transformed. IST is applied column by column along
    the indirect dimension.

    Parameters
    ----------
    input_spec      : input array, shape (n_direct, n_indirect)
    sampling_schedule: indices of sampled indirect points (0-based)
    threshold       : IST soft threshold fraction
    terminate       : L2 norm termination criterion (mode 2 only)
    convergence_tol : relative change termination criterion (mode 1 only)
    max_iter        : maximum number of iterations
    mode            : 1 = converge on relative change, 2 = converge on L2 norm
    """

    if r2 is None:
        r2 = 1.2/(max_time)+1.0

    window_vals = np.arange(input_spec.shape[-1])/input_spec.shape[-1]
    window = np.exp(-window_vals*r2*max_time)


    def _reconstruct_until_convergence(nus_fid: NDArray) -> NDArray:
        reconstructed = np.zeros_like(nus_fid)
        prev_reconstructed = np.zeros_like(nus_fid)
        prev_norm = 0.0

        for iteration in range(1, max_iter + 1):
            nus_fid, threshold_signal, _ = ist_iteration_2d(nus_fid, threshold, sampling_schedule, window)
            reconstructed += threshold_signal

            curr_norm = np.sqrt(np.vdot(reconstructed, reconstructed).real)
            relative_change = abs(curr_norm - prev_norm) / (curr_norm + 1e-10)
            prev_norm = curr_norm

            if relative_change < convergence_tol:
                if verb:
                    print(f"  converged at iteration {iteration} — "
                        f"relative change: {relative_change:.2e}")
                break
            if iteration == max_iter:
                if verb:
                    print(f"  reached max iterations of {iteration+1} "
                        f"relative change: {relative_change:.2e}")

        return reconstructed

    def _reconstruct_until_l2(nus_fid: NDArray) -> NDArray:
        reconstructed = np.zeros_like(nus_fid)

        for iteration in range(1, max_iter + 1):
            nus_fid, threshold_signal, l2_norm = ist_iteration_2d(nus_fid, threshold, sampling_schedule, window)
            reconstructed += threshold_signal

            if l2_norm <= terminate:
                break
            if iteration == max_iter:
                if verb:
                    print(f"  reached max iterations — L2 norm: {l2_norm:.2e}")

        return reconstructed
    
    def _process_slice(i):
        if verb:
            print(f"Doing IST slice {i}")
        reconstructed_slice = reconstruct(input_spec[i].copy())
        return i, reconstructed_slice


    recon_spec = np.zeros_like(input_spec)
    reconstruct = _reconstruct_until_convergence if mode == 1 else _reconstruct_until_l2
    recon_buffer = np.zeros_like(input_spec)

    # for i in range(input_spec.shape[0]):
    #     print(f"IST slice {i + 1} / {input_spec.shape[0]}")
    #     recon_buffer[i] = reconstruct(input_spec[i].copy())
    
    results = Parallel(n_jobs = -1)(delayed(_process_slice)(i) for i in range(input_spec.shape[0]))
    for i, result in results:
        recon_buffer[i] = result

    recon_spec = fft.ifft(recon_buffer, axis = -1)

    return recon_spec