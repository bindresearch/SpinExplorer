import numpy as np
import nmrglue as ng
from numpy.typing import NDArray
from typing import Union
import pyfftw
import pyfftw.interfaces.numpy_fft as fft
from SpinExplorer.SpinProcess.Processing.IST.sampling_utils import apply_sampling_schedule_nd, apply_sampling_schedule_to_2D_signal
from joblib import Parallel, delayed


pyfftw.interfaces.cache.enable()
pyfftw.interfaces.cache.set_keepalive_time(30)  # keep plans cached for 30s


def write_as_nmrpipe(array: NDArray, spec_dic: dict, outfile: str):
    #print(array.shape)
    ng.pipe.write(outfile, spec_dic, array, overwrite=True)

# def get_thresh_signal_3d(signal_real: NDArray, signal_imaginary: NDArray, threshold: float) -> tuple[NDArray, NDArray]:
    
#     max_val = np.max(signal_real)
#     signal_cutoff = threshold * max_val
    
#     # Compute soft threshold on real part
#     thresh_real = signal_real - signal_cutoff
#     thresh_real[thresh_real < 0.0] = 0.0
    
#     # Build a scaling mask from real part: where real was below cutoff, scale is 0;
#     # elsewhere, scale is (thresholded / original) to apply proportional attenuation
#     scale = np.where(signal_real > signal_cutoff, thresh_real / signal_real, 0.0)
    
#     # Apply same scaling to imaginary part
#     thresh_imag = signal_imaginary * scale
    
#     return thresh_real, thresh_imag




def get_thresh_signal_3d(signal_real: NDArray,
                        signal_imag: NDArray,
                        threshold: float):

    def soft_thresh(x, threshold):
        mag = np.abs(x)
        cutoff = threshold * np.max(mag)
        scale = np.maximum(0.0, mag - cutoff) / (mag + 1e-12)
        return x * scale

    thresh_real = soft_thresh(signal_real, threshold)
    thresh_imag = soft_thresh(signal_imag, threshold)

    return thresh_real, thresh_imag

"""
# protect the original functions that are working
def ist_iteration_3d(data: NDArray,
                    threshold: float,
                    sampling_schedule: Union[list[int], NDArray]) -> tuple[NDArray, NDArray, NDArray, np.floating]:

    # FT along faster indirect (axis=1): collapse cos/sin pairs → complex, FT, take real

    data = data[:, 0::2] + 1.0j * data[:, 1::2]       
    #data = fft.fftshift(fft.fft(data, axis=-1), axes=-1)
    data = ng.proc_base.fft_positive(data)
    data_real = np.real(data)                                 
    data_imag = np.imag(data)
    
    data_real = np.transpose(data_real)
    data_imag = np.transpose(data_imag)

    data_real = data_real[:,0::2] + 1.0j*data_real[:,1::2]
    data_imag = data_imag[:,0::2] + 1.0j*data_imag[:,1::2]
    

    data_real = ng.proc_base.fft_positive(data_real)
    data_imag = ng.proc_base.fft_positive(data_imag)
    
    thresh_sig_real, thresh_sig_imag = get_thresh_signal_3d(data_real, data_imag, threshold)
    leftover_sig_real = data_real - thresh_sig_real
    leftover_sig_imag = data_imag - thresh_sig_imag
    l2_norm = np.linalg.norm(leftover_sig_real)
    thresh_to_write = retrieve_signal_ist_3d(thresh_sig_real, thresh_sig_imag)


    # Inverse: re-expand both axes back to interleaved form
    leftover_sig = retrieve_signal_ist_3d(leftover_sig_real, leftover_sig_imag)  # (n3*2, n2*2)
    
    leftover_sig = apply_sampling_schedule_nd(leftover_sig, sampling_schedule, (1, 0))
  
    return leftover_sig, thresh_sig_real, thresh_sig_imag, l2_norm


def retrieve_signal_ist_3d(leftover_sig_real, leftover_sig_imag):
    
    leftover_sig_real = ng.proc_base.ifft_positive(leftover_sig_real)
    leftover_sig_imag = ng.proc_base.ifft_positive(leftover_sig_imag)
    

    def unpack_complex(x):
        out = np.zeros(x.shape[:-1] + (x.shape[-1]*2,), dtype=np.float32)
        out[..., 0::2] = np.real(x)
        out[..., 1::2] = np.imag(x)
        return out

    leftover_sig_real = unpack_complex(leftover_sig_real)
    leftover_sig_imag = unpack_complex(leftover_sig_imag)

    leftover_sig_real = np.transpose(leftover_sig_real)
    leftover_sig_imag = np.transpose(leftover_sig_imag)

    data = leftover_sig_real + 1j*leftover_sig_imag


    data = ng.proc_base.ifft_positive(data)
    data = unpack_complex(data)

    return data
"""

def ist_iteration_3d(data: NDArray,
                    threshold: float,
                    sampling_schedule: Union[list[int], NDArray]) -> tuple[NDArray, NDArray, NDArray, np.floating]:

    # FT along faster indirect (axis=1): collapse cos/sin pairs → complex, FT, take real

    data = data[:, 0::2] + 1.0j * data[:, 1::2]    
    data = fft.fftshift(fft.fft(data, axis=-1), axes=-1)
    #data = ng.proc_base.fft_positive(data)
    data_real = np.real(data)                                 
    data_imag = np.imag(data)
    
    data_real = np.transpose(data_real)
    data_imag = np.transpose(data_imag)

    data_real = data_real[:,0::2] + 1.0j*data_real[:,1::2]
    data_imag = data_imag[:,0::2] + 1.0j*data_imag[:,1::2]    

    # data_real = ng.proc_base.fft_positive(data_real)
    # data_imag = ng.proc_base.fft_positive(data_imag)

    # data_real = fft.fftshift(fft.fft(data_real, axis=-1), axes=-1)
    # data_imag = fft.fftshift(fft.fft(data_imag, axis=-1), axes=-1)
        
    stacked = np.stack([data_real, data_imag])  # (2, ...)
    stacked = fft.fftshift(fft.fft(stacked, axis=-1), axes=-1)
    data_real, data_imag = stacked[0], stacked[1]

    thresh_sig_real, thresh_sig_imag = get_thresh_signal_3d(data_real, data_imag, threshold)
    leftover_sig_real = data_real - thresh_sig_real
    leftover_sig_imag = data_imag - thresh_sig_imag
    l2_norm = np.linalg.norm(leftover_sig_real)


    # Inverse: re-expand both axes back to interleaved form
    leftover_sig = retrieve_signal_ist_3d(leftover_sig_real, leftover_sig_imag)  # (n3*2, n2*2)
    
    leftover_sig = apply_sampling_schedule_nd(leftover_sig, sampling_schedule, (1, 0))
  
    return leftover_sig, thresh_sig_real, thresh_sig_imag, l2_norm


def retrieve_signal_ist_3d(leftover_sig_real, leftover_sig_imag):
    
    # leftover_sig_real = ng.proc_base.ifft_positive(leftover_sig_real)
    # leftover_sig_imag = ng.proc_base.ifft_positive(leftover_sig_imag)

    #leftover_sig_real = fft.ifft(fft.ifftshift(leftover_sig_real, axes = -1), axis=-1)
    #leftover_sig_imag = fft.ifft(fft.ifftshift(leftover_sig_imag, axes = -1), axis=-1)
    
    stacked = np.stack([leftover_sig_real, leftover_sig_imag])
    stacked = fft.ifft(fft.ifftshift(stacked, axes=-1), axis=-1)
    leftover_sig_real, leftover_sig_imag = stacked[0], stacked[1]

    def unpack_complex(x):
        out = np.empty(x.shape[:-1] + (x.shape[-1]*2,), dtype=np.float32)
        out[..., 0::2] = x.real
        out[..., 1::2] = x.imag
        return out

    leftover_sig_real = unpack_complex(leftover_sig_real)
    leftover_sig_imag = unpack_complex(leftover_sig_imag)

    
    #leftover_sig_real = np.transpose(leftover_sig_real)
    #leftover_sig_imag = np.transpose(leftover_sig_imag)
    
    leftover_sig_real = np.ascontiguousarray(leftover_sig_real.T)
    leftover_sig_imag = np.ascontiguousarray(leftover_sig_imag.T)

    data = leftover_sig_real + 1j*leftover_sig_imag


    # data = ng.proc_base.ifft_positive(data)
    data = fft.ifft(fft.ifftshift(data, axes = -1), axis=-1)
    data = unpack_complex(data)

    return data


def ist_3d(input_spec: NDArray,
           sampling_schedule: Union[list[int], NDArray],
           threshold: float = 0.9,
           terminate: float = 1e-4,
           convergence_tol: float = 1e-8,
           max_iter: int = 800,
           mode: int = 1,
           sched_ord: int = 0,
           verb: bool = False,
           ist_callback = None) -> NDArray:
    
    """
    IST reconstruction of 3D NUS data
    """
    if sched_ord == 1:
        sampling_schedule = np.asarray(sampling_schedule)
        sampling_schedule = sampling_schedule[:, ::-1]


    def _reconstruct_until_convergence(nus_fid: NDArray) -> tuple[NDArray,NDArray]:
        reconstructed_r = None
        
        reconstructed_i = None
        prev_norm = 0.0


        for iteration in range(1, max_iter + 1):
            nus_fid, threshold_sig_real, threshold_sig_imag, _ = ist_iteration_3d(nus_fid, threshold, sampling_schedule)
            if reconstructed_r is None:
                reconstructed_r = np.zeros_like(threshold_sig_real)
                reconstructed_i = np.zeros_like(threshold_sig_imag)
                
            reconstructed_r += threshold_sig_real 
            reconstructed_i += threshold_sig_imag

            curr_norm = np.sqrt(np.vdot(reconstructed_r, reconstructed_r).real)
            relative_change = abs(curr_norm - prev_norm) / (curr_norm + 1e-10)
            prev_norm = curr_norm

            if relative_change < convergence_tol:
                if(verb):
                    print(f"  converged at iteration {iteration} — "
                        f"relative change: {relative_change:.2e}")
                break
            if iteration == max_iter:
                if(verb):
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


    results = Parallel(n_jobs = -1, return_as="generator")(delayed(_process_slice_3d)(i) for i in range(input_spec.shape[0]))
    for i, result in results:
        if(ist_callback!=None):
            ist_callback()
        recon_spec[i] = result

    # for i in range(input_spec.shape[0]):
    # #for i in range(1):
    #     print(f"IST slice {i + 1} / {input_spec.shape[0]}")
    #     slice_data = input_spec[i].copy()
    #     recon_real, recon_imag = reconstruct(slice_data)
    #     recon_slice = retrieve_signal_ist_3d(recon_real, recon_imag)
    #     recon_spec[i] = recon_slice

    if verb:
        print(f"Finished IST slice {i}")
    return recon_spec


# def get_thresh_signal(signal_ft: NDArray, threshold: float)->NDArray:
    
#     max = np.max(signal_ft)
#     signal_cutoff = threshold*max
#     thresh_sig = signal_ft - signal_cutoff
#     thresh_sig[thresh_sig<0.0] = 0.0
#     return thresh_sig


def get_thresh_signal(signal_ft: NDArray,
                        threshold: float):

    def soft_thresh(x, threshold):
        mag = np.abs(x)
        cutoff = threshold * np.max(mag)
        scale = np.maximum(0.0, mag - cutoff) / (mag + 1e-12)
        return x * scale


    thresh_real = soft_thresh(signal_ft, threshold)

    return thresh_real


def ist_iteration_2d(nus_fid:NDArray, 
                  threshold:float, 
                  sampling_schedule:Union[list[int],NDArray])->tuple[NDArray,NDArray,np.floating]:
    # consider whether we want to put zero-filling into this
    
    """
    # this is the original function
    signal_ft = np.fft.fft(nus_fid)
    threshold_sig = get_thresh_signal(signal_ft, threshold)
    leftover_sig = signal_ft-threshold_sig
    l2_norm = np.linalg.norm(leftover_sig)
    leftover_fid = np.fft.ifft(leftover_sig)
    leftover_fid = apply_sampling_schedule_nd(leftover_fid, sampling_schedule, (0,))
    """

    # with faster fft 
    signal_ft = fft.fft(nus_fid)
    #signal_ft[...,0] = signal_ft[...,0]/0.5 
    threshold_sig = get_thresh_signal(signal_ft, threshold)
    leftover_sig = signal_ft-threshold_sig
    l2_norm = l2_norm = np.sqrt(np.vdot(leftover_sig, leftover_sig).real)

    leftover_fid = fft.ifft(leftover_sig)
    #leftover_fid[...,0] = leftover_fid[...,0]*2.0
    leftover_fid = apply_sampling_schedule_nd(leftover_fid, sampling_schedule, (0,))

    return leftover_fid, threshold_sig, l2_norm


def ist_2d(input_spec: NDArray,
           sampling_schedule: Union[list[int], NDArray],
           threshold: float = 0.9,
           terminate: float = 0.001,
           convergence_tol: float = 1e-10,
           max_iter: int = 4000,
           mode: int = 1,
           verb: bool = False,
           ist_callback = None) -> NDArray:
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
    ist_callback    : A function to update the SpinExplorer counter so that the % through
                      reconstruction can be reported.
    """

    def _reconstruct_until_convergence(nus_fid: NDArray) -> NDArray:
        reconstructed = np.zeros_like(nus_fid)
        prev_reconstructed = np.zeros_like(nus_fid)
        prev_norm = 0.0

        nus_fid[...,0] = nus_fid[...,0]/2.0
        for iteration in range(1, max_iter + 1):
            nus_fid, threshold_signal, _ = ist_iteration_2d(nus_fid, threshold, sampling_schedule)
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
            nus_fid, threshold_signal, l2_norm = ist_iteration_2d(nus_fid, threshold, sampling_schedule)
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
    
    results = Parallel(n_jobs = -1, return_as="generator")(delayed(_process_slice)(i) for i in range(input_spec.shape[0]))
    for i, result in results:
        if(ist_callback!=None):
            ist_callback()
        recon_buffer[i] = result


    recon_spec = fft.ifft(recon_buffer, axis = -1)
    recon_spec[...,0] = recon_spec[...,0]*2.0

    return recon_spec