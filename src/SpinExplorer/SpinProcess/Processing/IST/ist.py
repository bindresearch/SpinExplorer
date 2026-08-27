import numpy as np
import nmrglue as ng # type: ignore
from numpy.typing import NDArray
from typing import Union, Optional
import pyfftw # type: ignore
import pyfftw.interfaces.numpy_fft as fft # type: ignore
from SpinExplorer.SpinProcess.Processing.IST.sampling_utils import apply_sampling_schedule_nd, apply_sampling_schedule_to_2D_signal
from joblib import Parallel, delayed # type: ignore
import copy

pyfftw.interfaces.cache.enable()
pyfftw.interfaces.cache.set_keepalive_time(30)  # keep plans cached for 30s


def write_as_nmrpipe(array: NDArray, spec_dic: dict, outfile: str):
    ng.pipe.write(outfile, spec_dic, array, overwrite=True)


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

    leftover_real = signal_real - thresh_real
    leftover_max_val = np.max(np.abs(leftover_real))

    thresh_fid = retrieve_signal_ist_3d(thresh_real, thresh_imag)

    thresh_real, thresh_imag = pack_signal_ist_3d(thresh_fid)


    return thresh_real, thresh_imag, leftover_max_val


def ist_iteration_3d(data: NDArray,
                    threshold: float,
                    sampling_schedule: Union[list[int], NDArray]) -> tuple[NDArray, NDArray, NDArray, np.floating, np.floating]:

    # FT along faster indirect (axis=1): collapse cos/sin pairs → complex, FT, take real

    data_real, data_imag = pack_signal_ist_3d(data)

    thresh_sig_real, thresh_sig_imag, max_val = get_thresh_signal_3d(data_real, data_imag, threshold)

    leftover_sig_real = data_real - thresh_sig_real
    leftover_sig_imag = data_imag - thresh_sig_imag
    l2_norm = np.linalg.norm(leftover_sig_real)

    # Inverse: re-expand both axes back to interleaved form
    leftover_sig = retrieve_signal_ist_3d(leftover_sig_real, leftover_sig_imag)  # (n3*2, n2*2)
    
    leftover_sig = apply_sampling_schedule_nd(leftover_sig, sampling_schedule, (1, 0))

    return leftover_sig, thresh_sig_real, thresh_sig_imag, l2_norm, leftover_max_val


def pack_signal_ist_3d(data: NDArray):

    data = data[:, 0::2] + 1.0j * data[:, 1::2]
    data = np.pad(data, pad_width=[(0,0),(0,data.shape[-1])])       
    data = fft.fft(data, axis=-1)
    data_real = np.real(data)                                 
    data_imag = np.imag(data)
    
    data_real = np.transpose(data_real)
    data_imag = np.transpose(data_imag)

    data_real = data_real[:,0::2] + 1.0j*data_real[:,1::2]
    data_imag = data_imag[:,0::2] + 1.0j*data_imag[:,1::2]

    data_real = np.pad(data_real, pad_width=[(0,0),(0,data_real.shape[-1])])
    data_imag = np.pad(data_imag, pad_width=[(0,0),(0,data_imag.shape[-1])])
        
    data_real = fft.fft(data_real, axis=-1)
    data_imag = fft.fft(data_imag, axis=-1)
    return data_real, data_imag

def retrieve_signal_ist_3d(leftover_sig_real: NDArray, leftover_sig_imag: NDArray):

    
    leftover_sig_real, leftover_sig_imag = fid_from_absorption(leftover_sig_real), fid_from_absorption(leftover_sig_imag)

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
    
    data = fid_from_absorption(data)
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
           ist_callback = None) -> tuple[NDArray,int]:
    
    """
    IST reconstruction of 3D NUS data
    """
    converged_results = 0

    if sched_ord == 1:
        sampling_schedule = np.asarray(sampling_schedule)
        sampling_schedule = sampling_schedule[:, ::-1]


    cancelled_button_pressed = False


    def _reconstruct_until_convergence(nus_fid: NDArray) -> tuple[NDArray,NDArray,bool]:

        prev_norm = 0.0

        converged = False


        for iteration in range(1, max_iter + 1):

            # Check to see if a user has cancelled the IST reconstruction

            nus_fid, threshold_sig_real, threshold_sig_imag, _, max_val = ist_iteration_3d(nus_fid, threshold, sampling_schedule)

            if(iteration==1):
                reconstructed_r = np.zeros_like(threshold_sig_real)
                reconstructed_i = np.zeros_like(threshold_sig_imag)

                thresh_signal_real_max = max_val

                
            reconstructed_r += threshold_sig_real 
            reconstructed_i += threshold_sig_imag

            curr_norm = np.sqrt(np.vdot(reconstructed_r, reconstructed_r).real)
            relative_change = abs(curr_norm - prev_norm) / (curr_norm + 1e-10)
            prev_norm = curr_norm

            if(leftover_max_val/max_val < convergence_tol):
                if(verb):
                    print(f"  converged at iteration {iteration} — "
                        f"relative change: {relative_change:.2e}")
                converged = True
                break
            if iteration == max_iter:
                if(verb):
                    print(f"  reached max iterations of {iteration} — "
                        f"relative change: {relative_change:.2e}")


        # leftover_real, leftover_imag = pack_signal_ist_3d(nus_fid)

        return reconstructed_r, reconstructed_i, converged

    def _reconstruct_until_l2(nus_fid: NDArray) -> tuple[NDArray,NDArray,bool]:
        l2_norm = terminate * 1000.0

        converged = False

        for iteration in range(1, max_iter + 1):
            nus_fid, threshold_sig_real, threshold_sig_imag, l2_norm, _ = ist_iteration_3d(nus_fid, threshold, sampling_schedule)

            if iteration==1:
                reconstructed_r = np.zeros_like(threshold_sig_real)
                reconstructed_i = np.zeros_like(threshold_sig_imag)

                
            reconstructed_r += threshold_sig_real 
            reconstructed_i += threshold_sig_imag

            if l2_norm <= terminate:
                if verb:
                    print(f"  converged at iteration {iteration} — "
                        f"l2 norm: {l2_norm:.2e}")
                converged = True
                break
            if iteration == max_iter:
                if verb:
                    print(f"  reached max iterations — L2 norm: {l2_norm:.2e}")

        return reconstructed_r, reconstructed_i, converged

    recon_spec = np.zeros_like(input_spec)
    reconstruct = _reconstruct_until_convergence if mode == 1 else _reconstruct_until_l2
    
    def _process_slice_3d(i):
        if verb:
            print(f"Doing IST slice {i}")
        
        slice_data = input_spec[i].copy()
        recon_real, recon_imag, converged = reconstruct(slice_data)
        recon_slice = retrieve_signal_ist_3d(recon_real, recon_imag)

        return i, recon_slice, converged


    
    results = Parallel(n_jobs = -1, return_as="generator")(delayed(_process_slice_3d)(i) for i in range(input_spec.shape[0]))
    for i, result, converged in results:
        if(ist_callback!=None):
            continue_reconstruction = ist_callback(converged)
            if(continue_reconstruction == False):
                return input_spec, converged_results
        recon_spec[i] = result
        if(converged==True):
            converged_results+=1


    return recon_spec, converged_results

def fid_from_absorption(S):
    """S: 2N complex spectrum, FFT order, last axis. Returns N-point complex FID."""
    N = S.shape[-1] // 2

    return fft.ifft(S, axis=-1)[..., :N].copy()


def get_thresh_signal(signal_ft: NDArray,
                        threshold: float):

    def soft_thresh(x, threshold):
        mag = np.abs(x)
        cutoff = threshold * np.max(mag)
        scale = np.maximum(0.0, mag - cutoff) / (mag + 1e-12)
        return x * scale

    thresh_real = soft_thresh(signal_ft, threshold)
    max_val = np.max(np.abs(thresh_real))
    #thresh_fid = fft.ifft(thresh_real)
    thresh_fid = fid_from_absorption(thresh_real)

    return thresh_real, thresh_fid, max_val


def ist_iteration_2d(nus_fid:NDArray, 
                  threshold:float, 
                  sampling_schedule:Union[list[int],NDArray])->tuple[NDArray,NDArray,np.floating,np.floating]:
    

    nus_fid = np.pad(nus_fid, pad_width = ([0,nus_fid.shape[0]])) 
    signal_ft = fft.fft(nus_fid)
    threshold_sig, thresh_fid, max_val = get_thresh_signal(signal_ft, threshold)
    leftover_sig = signal_ft-threshold_sig
    l2_norm = np.sqrt(np.vdot(leftover_sig, leftover_sig).real)

    # leftover_fid = fft.ifft(leftover_sig)
    leftover_fid = fid_from_absorption(leftover_sig)
    # leftover_fid = nus_fid - thresh_fid # commented this out it doesn't make sense
    leftover_fid = apply_sampling_schedule_nd(leftover_fid, sampling_schedule, (0,))

    return leftover_fid, thresh_fid, threshold_sig, l2_norm, leftover_max_val


def ist_2d(input_spec: NDArray,
           sampling_schedule: Union[list[int], NDArray],
           threshold: float = 0.9,
           terminate: float = 0.001,
           convergence_tol: float = 1e-6,
           max_iter: int = 4000,
           mode: int = 1,
           verb: bool = False,
           ist_callback = None) -> tuple[NDArray,int]:
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
    converged_results = 0

    def _reconstruct_until_convergence(nus_fid: NDArray) -> tuple[NDArray,bool]:
        reconstructed = np.zeros_like(nus_fid)
        prev_reconstructed = np.zeros_like(nus_fid)
        prev_norm = 0.0

        converged = False

        nus_fid_initial = copy.deepcopy(nus_fid)

        for iteration in range(1, max_iter + 1):
            nus_fid, threshold_signal, _, max_val = ist_iteration_2d(nus_fid, threshold, sampling_schedule)
            reconstructed += threshold_signal

            curr_norm = np.sqrt(np.vdot(reconstructed, reconstructed).real)
            relative_change = abs(curr_norm - prev_norm) / (curr_norm + 1e-10)
            prev_norm = curr_norm


            if(iteration==1):
                if((leftover_max_val+np.abs(np.max(threshold_ft)))/max_val < convergence_tol):
                    return nus_fid_initial, True
            
            if(leftover_max_val/max_val < convergence_tol):
                if verb:
                    print(f"  converged at iteration {iteration} — "
                        f"relative change: {relative_change:.2e}")
                converged = True
                break
            if iteration == max_iter:
                if verb:
                    print(f"  reached max iterations of {iteration+1} "
                        f"relative change: {relative_change:.2e}")
                break

        # return reconstructed + nus_fid, converged
        return reconstructed, converged

    def _reconstruct_until_l2(nus_fid: NDArray) -> tuple[NDArray,bool]:
        reconstructed = np.zeros_like(nus_fid)

        converged = False

        for iteration in range(1, max_iter + 1):
            nus_fid, threshold_signal, l2_norm, max_val = ist_iteration_2d(nus_fid, threshold, sampling_schedule)
            reconstructed += threshold_signal

            threshold*=0.3

            if l2_norm <= terminate:
                converged = True
                break
            if iteration == max_iter:
                if verb:
                    print(f"  reached max iterations — L2 norm: {l2_norm:.2e}")

        return reconstructed, converged
    
    def _process_slice(i):
        if verb:
            print(f"Doing IST slice {i}")
        reconstructed_slice, converged = reconstruct(input_spec[i].copy())
        return i, reconstructed_slice, converged


    recon_spec = np.zeros_like(input_spec)
    reconstruct = _reconstruct_until_convergence if mode == 1 else _reconstruct_until_l2
    # recon_buffer = np.zeros_like(input_spec)

    # for i in range(input_spec.shape[0]):
    #     print(f"IST slice {i + 1} / {input_spec.shape[0]}")
    #     recon_buffer[i] = reconstruct(input_spec[i].copy())
    
    
    results = Parallel(n_jobs = -1, return_as="generator")(delayed(_process_slice)(i) for i in range(input_spec.shape[0]))
    for i, result, converged in results:
        if(ist_callback!=None):
            continue_reconstruction = ist_callback(converged)
            if(continue_reconstruction == False):
                return input_spec, converged_results
        recon_spec[i] = result
        if(converged==True):
            converged_results+=1


    # recon_spec = fft.ifft(recon_buffer, axis = -1)

    return recon_spec, converged_results






def ist_2d_as_plane(input_spec: NDArray,
           sampling_schedule: Union[list[int], NDArray],
           max_time: float,
           r2: Optional[float] = None,
           threshold: float = 0.9,
           terminate: float = 0.001,
           convergence_tol: float = 1e-6,
           max_iter: int = 4000,
           verb: bool = False,
           ist_callback = None,
           max_val=1) -> NDArray:
    """
    IST reconstruction of 2D NUS data. IST is to be applied to the whole
    plane at once. All direct dimension points are replaced after each
    iteration.

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
    converged_results = 0

    if r2 is None:
        r2 = 1.2/(max_time)+1.0

    r2 = 0

    window_vals = np.arange(input_spec.shape[-1])/input_spec.shape[-1]
    window_dim2 = np.exp(-window_vals*r2*max_time)
    window_dim3 = np.exp(-window_vals*r2*max_time)


    def _reconstruct_until_convergence(nus_fid: NDArray) -> tuple[NDArray,NDArray]:
        reconstructed_r = None
        
        reconstructed_i = None
        prev_norm = 0.0

        converged = False


        for iteration in range(1, max_iter + 1):


            nus_fid, threshold_sig_real, threshold_sig_imag, _, leftover_max_val = ist_iteration_3d(nus_fid, threshold, sampling_schedule, window_dim2, window_dim3)
            if reconstructed_r is None:
                reconstructed_r = np.zeros_like(threshold_sig_real)
                reconstructed_i = np.zeros_like(threshold_sig_imag)
                
            reconstructed_r += threshold_sig_real 
            reconstructed_i += threshold_sig_imag

            curr_norm = np.sqrt(np.vdot(reconstructed_r, reconstructed_r).real)
            relative_change = abs(curr_norm - prev_norm) / (curr_norm + 1e-10)
            prev_norm = curr_norm


            print(leftover_max_val/max_val < 0.05)

            # if(leftover_max_val/max_val < convergence_tol):
            if(leftover_max_val/max_val < 0.05):
                if(verb):
                    print(f"  converged at iteration {iteration} — "
                        f"relative change: {relative_change:.2e}")
                converged = True
                break
            if iteration == max_iter:
                if(verb):
                    print(f"  reached max iterations of {iteration} — "
                        f"relative change: {relative_change:.2e}")


        leftover_real, leftover_imag = pack_signal_ist_3d(nus_fid)

        # return reconstructed_r, reconstructed_i, converged
        return reconstructed_r+leftover_real, reconstructed_i+leftover_imag, converged

    
    
    recon_spec = np.zeros_like(input_spec)
    reconstruct = _reconstruct_until_convergence
    
    def _process_slice_3d():
        
        recon_real, recon_imag, converged = reconstruct(input_spec)
        recon_data = retrieve_signal_ist_3d(recon_real, recon_imag)

        return recon_data, converged


    
    results = _process_slice_3d()
    result, converged = results
    if(ist_callback!=None):
        continue_reconstruction = ist_callback(converged)
        if(continue_reconstruction == False):
            return input_spec, converged_results
    recon_spec = result
    if(converged==True):
        converged_results+=1


    return recon_spec, converged_results