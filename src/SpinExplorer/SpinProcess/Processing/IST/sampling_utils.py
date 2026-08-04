import numpy as np
from numpy.typing import NDArray
from typing import Union
import itertools
from collections import defaultdict

def apply_sampling_schedule_to_signal(signal:NDArray,sampling_schedule:Union[list[int],NDArray])->NDArray:
    nus_sig = signal*1.0
    nus_sig[~np.isin(np.arange(len(nus_sig)), sampling_schedule)] = 0
    return nus_sig


def apply_sampling_schedule_to_2D_signal(spectrum:NDArray,sampling_schedule:Union[list[int],NDArray])->NDArray:
    # we have a 2D NMR spectrum. Direct dimension is fully sampled. Apply nus to indirect dimension
    # this is interleaved real imaginary and increases in number of points. If a point is not sampled both
    # the real and imaginary should be set to 0
    nus_sig = spectrum.copy()
    n_indirect = spectrum.shape[0] // 2

    # Convert indirect indices to row pairs
    # indirect point i → rows 2i (cos) and 2i+1 (sin)
    sampled_indirect = np.asarray(sampling_schedule)
    all_indirect = np.arange(n_indirect)
    unsampled = all_indirect[~np.isin(all_indirect, sampled_indirect)]

    # Zero both cos and sin rows for each unsampled indirect point
    nus_sig[unsampled * 2, :]     = 0  # cos rows
    nus_sig[unsampled * 2 + 1, :] = 0  # sin rows

    return nus_sig


def apply_sampling_schedule_to_2D_signal_and_collapse(spectrum: NDArray, sample_dict: dict, 
                                              sampling_schedule: Union[list[int], NDArray]) -> tuple[NDArray,dict]:
    """
    Prepare a 2D NMR dataset for IST reconstruction by:
    1. Extracting only the sampled rows from the interleaved spectrum
    2. Updating the spectral dictionary to reflect the reduced size
    3. Writing the collapsed dataset to an NMRPipe file for IST processing

    Parameters
    ----------
    spectrum          : full interleaved array, shape (n_indirect*2, n_direct)
    sample_dict       : NMRPipe header dictionary
    sampling_schedule : indices of sampled indirect points (0-based)
    outfile           : output NMRPipe filename
    """
    sampling_schedule = np.asarray(sampling_schedule)
    n_sampled = len(sampling_schedule)

    # Extract only sampled row pairs (cos + sin) from the interleaved spectrum
    # indirect point i → rows 2i (cos) and 2i+1 (sin)
    sampled_rows = np.empty((n_sampled * 2, spectrum.shape[1]), dtype=spectrum.dtype)
    sampled_rows[0::2] = spectrum[sampling_schedule * 2,     :]  # cos rows
    sampled_rows[1::2] = spectrum[sampling_schedule * 2 + 1, :]  # sin rows

    # Update header to reflect collapsed indirect dimension size
    ist_dict = sample_dict.copy()
    ist_dict['FDSPECNUM'] = float(n_sampled * 2)    # number of rows in collapsed file
    ist_dict['FDF1TDSIZE'] = float(n_sampled)        # indirect TD size
    ist_dict['FDF1APOD']   = float(n_sampled)        # apodisation size
    ist_dict['FDNUSDIM']   = 1.0                     # flag that indirect dim is non-uniform
    ist_dict['FDF1CENTER'] = float(n_sampled / 2 + 1)

    # output collapsed nus dataset and dictionary

    return sampled_rows, ist_dict

def apply_sampling_schedule_to_nd_signal_and_collapse(spectrum: NDArray, sample_dict: dict,
                                              sampling_schedule: Union[list[int], NDArray]) -> tuple[NDArray, dict]:
    """
    Prepare an nD NMR dataset for IST reconstruction by:
    1. Extracting only the sampled points from the interleaved spectrum
    2. Updating the spectral dictionary to reflect the reduced size

    Parameters
    ----------
    spectrum          : full interleaved array
                        2D: shape (n_indirect*2, n_direct)
                        3D: shape (n_indirect2*2, n_indirect1*2, n_direct)
                        4D: shape (n_indirect3*2, n_indirect2*2, indirect1*2, n_direct)
    sample_dict       : NMRPipe header dictionary
    sampling_schedule : 1D array of indices for 2D spectra, or
                        (n_sampled, n_indirect_dims) array of index tuples for nD

    Returns
    -------
    sampled_data : collapsed array containing only sampled points
    ist_dict     : updated NMRPipe header dictionary
    """
    sampling_schedule = np.asarray(sampling_schedule)
    ist_dict = sample_dict.copy()

    if sampling_schedule.ndim == 1:
        n_sampled = len(sampling_schedule)

        sampled_data = np.empty((n_sampled * 2, spectrum.shape[-1]), dtype=spectrum.dtype)
        sampled_data[0::2] = spectrum[sampling_schedule * 2,     :]
        sampled_data[1::2] = spectrum[sampling_schedule * 2 + 1, :]

        ist_dict['FDSPECNUM']  = float(n_sampled * 2)  # for 2D, FDSPECNUM = n_sampled * 2
        ist_dict['FDF1TDSIZE'] = float(n_sampled)
        ist_dict['FDF1APOD']   = float(n_sampled)
        ist_dict['FDF1CENTER'] = float(n_sampled / 2 + 1)
        ist_dict['FDNUSDIM']   = 1.0

    else:
        # nD spectrum — multiple indirect dimensions
        # sampling_schedule shape: (n_sampled, n_indirect_dims)
        n_sampled  = len(sampling_schedule)
        n_indirect_dims = sampling_schedule.shape[1]
        n_direct   = spectrum.shape[-1]

        # Each sampled point maps to 2^n_indirect_dims rows (cos/sin per indirect dim)
        # e.g. 3D: each (i2, i1) → 4 combinations: (cos2,cos1), (cos2,sin1), (sin2,cos1), (sin2,sin1)
        n_combos   = 2 ** n_indirect_dims
        sampled_data = np.empty((n_sampled * n_combos, n_direct), dtype=spectrum.dtype)

        for s_idx, indices in enumerate(sampling_schedule):
            # For each indirect dim, the interleaved row index is i*2 (cos) or i*2+1 (sin)
            # Generate all cos/sin combinations via itertools.product
            row_options = [(idx * 2, idx * 2 + 1) for idx in indices]

            for combo_idx, combo in enumerate(itertools.product(*row_options)):
                # combo is e.g. (cos_row_dim2, cos_row_dim1) — index into spectrum
                out_row = s_idx * n_combos + combo_idx
                sampled_data[out_row] = spectrum[combo + (slice(None),)]

        # Update header for each indirect dimension
        # Update header for each indirect dimension
        # sampling_schedule column 0 → faster indirect → FDF1 → spectrum axis 1
        # sampling_schedule column 1 → slower indirect → FDF3 → spectrum axis 0
        dim_map = {0: 'FDF1', 1: 'FDF3', 2: 'FDF4'}

        for dim_idx in range(n_indirect_dims):
            prefix = dim_map[dim_idx]
            unique_in_dim = len(np.unique(sampling_schedule[:, dim_idx]))
            ist_dict[f'{prefix}TDSIZE'] = float(unique_in_dim)
            ist_dict[f'{prefix}APOD']   = float(unique_in_dim)
            ist_dict[f'{prefix}CENTER'] = float(unique_in_dim / 2 + 1)

        # For a woven schedule sampled_data is flat: (n_sampled * n_combos, n_direct)
        # nmrglue needs to see this as a 2D array so treat it as a single
        # indirect dimension of size n_sampled * n_combos / 2
        ist_dict['FDSPECNUM']  = float(n_sampled * n_combos)
        ist_dict['FDDIMCOUNT'] = 2.0  # treat as 2D for NMRPipe IST processing
        ist_dict['FDNUSDIM']   = float(n_indirect_dims)

    return sampled_data, ist_dict




def inflate_spectra_2D_signal(spectrum: NDArray, sample_dict: dict,
                                   sampling_schedule: Union[list[int], NDArray],
                                   max_points: int) -> tuple[NDArray, dict]:
    """
    We have recorded n points of the indirect dimension according to a sampling schedule.
    We want to expand the indirect dimension up to max_points. The collected indirect points
    will match those in the sampling schedule. All other indirect points will be set to zero.

    Parameters
    ----------
    spectrum          : collapsed NUS array, shape (n_sampled*2, n_direct)
                        rows are interleaved cos/sin pairs
    sample_dict       : NMRPipe header dictionary
    sampling_schedule : indices of sampled indirect points (0-based), length n_sampled
    max_points        : full indirect dimension size to inflate to
    """
    sampling_schedule = np.asarray(sampling_schedule)
    n_sampled = len(sampling_schedule)
    n_direct  = spectrum.shape[1]

    if spectrum.shape[0] != n_sampled * 2:
        raise ValueError(f"spectrum has {spectrum.shape[0]} rows but sampling_schedule "
                         f"has {n_sampled} points — expected {n_sampled * 2} rows")

    # Allocate full-size zero array — shape (max_points*2, n_direct)
    inflated = np.zeros((max_points * 2, n_direct), dtype=spectrum.dtype)

    # Place sampled cos/sin row pairs at their correct indirect positions
    inflated[sampling_schedule * 2,     :] = spectrum[0::2]  # cos rows
    inflated[sampling_schedule * 2 + 1, :] = spectrum[1::2]  # sin rows

    # Update header to reflect inflated indirect dimension size
    inflated_dict = sample_dict.copy()
    inflated_dict['FDSPECNUM']  = float(max_points * 2)
    inflated_dict['FDF1TDSIZE'] = float(max_points)
    inflated_dict['FDF1APOD']   = float(max_points)
    inflated_dict['FDF1CENTER'] = float(max_points / 2 + 1)
    inflated_dict['FDNUSDIM']   = 1.0

    return inflated, inflated_dict


def apply_sampling_schedule_to_3d_signal_and_collapse(spectrum: NDArray, sample_dict: dict,
                                              sampling_schedule: Union[list[int], NDArray],
                                              acq_ord: int = 0) -> tuple[NDArray, dict]:
    """
    acq_ord = 0: this means we first do cos/sin modulation on the final dimension (t3) 
                 before the t2 dimension
    acq_ord = 1: This means we first do cos/sin modulation on the (t2) dimnesion before
                 the t3 dimension

    Prepare an nD NMR dataset for IST reconstruction by:
    1. Extracting only the sampled points from the interleaved spectrum
    2. Updating the spectral dictionary to reflect the reduced size

    Parameters
    ----------
    spectrum          : full interleaved array
                        3D: shape (n_indirect2*2, n_indirect1*2, n_direct)
    sample_dict       : NMRPipe header dictionary
    sampling_schedule : 1D array of indices for 2D spectra, or
                        (n_sampled, n_indirect_dims) array of index tuples for nD

    Returns
    -------
    sampled_data : collapsed array containing only sampled points
    ist_dict     : updated NMRPipe header dictionary
    """
    sampling_schedule = np.asarray(sampling_schedule)
    ist_dict = sample_dict.copy()


    # nD spectrum — multiple indirect dimensions
    # sampling_schedule shape: (n_sampled, n_indirect_dims)
    n_sampled  = len(sampling_schedule)
    n_indirect_dims = 2
    n_direct   = spectrum.shape[-1]

    n_combos = 4
    sampled_data = np.zeros((n_sampled * n_combos, n_direct), dtype=spectrum.dtype)
    idx = 0

    if acq_ord == 0:
        for sr0, sr1 in sampling_schedule:
            sampled_data[idx] = spectrum[2*sr0,   2*sr1,   :]; idx += 1
            sampled_data[idx] = spectrum[2*sr0+1, 2*sr1,   :]; idx += 1
            sampled_data[idx] = spectrum[2*sr0,   2*sr1+1, :]; idx += 1
            sampled_data[idx] = spectrum[2*sr0+1, 2*sr1+1, :]; idx += 1
    
    elif acq_ord == 1:
        for sr0, sr1 in sampling_schedule:
            sampled_data[idx] = spectrum[2*sr0,   2*sr1,   :]; idx += 1
            sampled_data[idx] = spectrum[2*sr0, 2*sr1+1,   :]; idx += 1
            sampled_data[idx] = spectrum[2*sr0+1,   2*sr1, :]; idx += 1
            sampled_data[idx] = spectrum[2*sr0+1, 2*sr1+1, :]; idx += 1
    

    dim_map = {0: 'FDF1', 1: 'FDF3', 2: 'FDF4'}

    for dim_idx in range(n_indirect_dims):
        prefix = dim_map[dim_idx]
        unique_in_dim = len(np.unique(sampling_schedule[:, dim_idx]))
        ist_dict[f'{prefix}TDSIZE'] = float(unique_in_dim)
        ist_dict[f'{prefix}APOD']   = float(unique_in_dim)
        ist_dict[f'{prefix}CENTER'] = float(unique_in_dim / 2 + 1)

    # For a woven schedule sampled_data is flat: (n_sampled * n_combos, n_direct)
    # nmrglue needs to see this as a 2D array so treat it as a single
    # indirect dimension of size n_sampled * n_combos / 2
    ist_dict['FDSPECNUM']  = float(n_sampled * n_combos)
    ist_dict['FDDIMCOUNT'] = 2.0  # treat as 2D for NMRPipe IST processing
    ist_dict['FDNUSDIM']   = float(n_indirect_dims)

    return sampled_data, ist_dict


def inflate_spectra_3D_signal(spectrum: NDArray, sample_dict: dict,
                                   sampling_schedule: Union[list[int], NDArray],
                                   max_points: Union[int, list[int]],
                                   acq_ord: int = 0) -> tuple[NDArray, dict]:
    inflated_dict     = sample_dict.copy()


    # nD spectrum — multiple indirect dimensions
    n_indirect_dims = 2 # this is for 3D spectra only

    if n_indirect_dims != sampling_schedule.shape[1]:
        raise ValueError(f"sampling schedule does not have correct shape for " 
                            f"3D spectra")

    if np.isscalar(max_points):
        max_points = [max_points] * n_indirect_dims
    max_points = list(max_points)

    if len(max_points) != n_indirect_dims:
        raise ValueError(f"max_points has {len(max_points)} entries but "
                            f"sampling_schedule has {n_indirect_dims} indirect dims")

    inflated = np.zeros(shape = (max_points[0]*2, max_points[1]*2, spectrum.shape[-1]), dtype = spectrum.dtype)
    
    idx = 0

    for sr0, sr1 in sampling_schedule:

        if acq_ord == 0:
            inflated[2*sr0,   2*sr1,   :] = spectrum[idx]; idx += 1
            inflated[2*sr0+1, 2*sr1,   :] = spectrum[idx]; idx += 1
            inflated[2*sr0,   2*sr1+1, :] = spectrum[idx]; idx += 1
            inflated[2*sr0+1, 2*sr1+1, :] = spectrum[idx]; idx += 1
        
        if acq_ord == 1: 
            inflated[2*sr0,   2*sr1,   :] = spectrum[idx]; idx += 1
            inflated[2*sr0+1, 2*sr1,   :] = spectrum[idx]; idx += 1
            inflated[2*sr0,   2*sr1+1, :] = spectrum[idx]; idx += 1
            inflated[2*sr0+1, 2*sr1+1, :] = spectrum[idx]; idx += 1


    # Update header
    dim_map = {0: 'FDF3', 1: 'FDF1', 2: 'FDF4'}
    for dim_idx, mp in enumerate(max_points):
        prefix = dim_map[dim_idx]
        inflated_dict[f'{prefix}TDSIZE'] = float(mp)
        inflated_dict[f'{prefix}APOD']   = float(mp)
        inflated_dict[f'{prefix}CENTER'] = float(mp / 2 + 1)

    inflated_dict['FDSPECNUM']  = float(max_points[1] * 2)  # slowest indirect dim * 2
    inflated_dict['FDDIMCOUNT'] = float(n_indirect_dims + 1)

    return inflated, inflated_dict



def inflate_spectra_nd_signal_ist(spectrum: NDArray, sample_dict: dict,
                                   sampling_schedule: Union[list[int], NDArray],
                                   max_points: Union[int, list[int]]) -> tuple[NDArray, dict]:
    sampling_schedule = np.asarray(sampling_schedule)

    inflated_dict     = sample_dict.copy()
    n_direct          = spectrum.shape[-1]

    if sampling_schedule.ndim == 1:
        # 2D spectrum — single indirect dimension
        n_sampled  = len(sampling_schedule)
        max_points = int(max_points)

        if spectrum.shape[0] != n_sampled * 2:
            raise ValueError(f"spectrum has {spectrum.shape[0]} rows but expected "
                             f"{n_sampled * 2} for {n_sampled} sampled points")

        inflated = np.zeros((max_points * 2, n_direct), dtype=spectrum.dtype)
        inflated[sampling_schedule * 2,     :] = spectrum[0::2]
        inflated[sampling_schedule * 2 + 1, :] = spectrum[1::2]

        inflated_dict['FDSPECNUM']  = float(max_points * 2)
        inflated_dict['FDF1TDSIZE'] = float(max_points)
        inflated_dict['FDF1APOD']   = float(max_points)
        inflated_dict['FDF1CENTER'] = float(max_points / 2 + 1)
        inflated_dict['FDNUSDIM']   = 1.0

    else:
        # nD spectrum — multiple indirect dimensions
        n_sampled       = len(sampling_schedule)
        n_indirect_dims = sampling_schedule.shape[1]
        n_combos        = 2 ** n_indirect_dims

        if np.isscalar(max_points):
            max_points = [max_points] * n_indirect_dims
        max_points = list(max_points)

        if len(max_points) != n_indirect_dims:
            raise ValueError(f"max_points has {len(max_points)} entries but "
                             f"sampling_schedule has {n_indirect_dims} indirect dims")

        # Handle flat 2D input from IST pipeline — shape (total_rows, n_direct)
        # where total_rows = product of all indirect interleaved sizes
        # We need to restore proper nD shape before inflating
        if spectrum.ndim == 2:
            total_max_rows = int(np.prod([mp * 2 for mp in max_points]))

            if spectrum.shape[0] == n_sampled * n_combos:
                # Still collapsed — inflate directly from flat representation
                #inflated_indirect_shape = tuple(mp * 2 for mp in reversed([num_points3, num_points2]))
                inflated_indirect_shape = tuple(mp * 2 for mp in max_points)

                inflated = np.zeros(inflated_indirect_shape + (n_direct,), dtype=spectrum.dtype)

                for s_idx, indices in enumerate(sampling_schedule):
                    row_options = [(int(idx) * 2, int(idx) * 2 + 1) for idx in indices]
                    for combo_idx, combo in enumerate(itertools.product(*row_options)):
                        src_row  = s_idx * n_combos + combo_idx
                        dest_idx = tuple(combo) + (slice(None),)
                        inflated[dest_idx] = spectrum[src_row]

            elif spectrum.shape[0] == total_max_rows:
                # Already full size but flat — just reshape to nD
                inflated_indirect_shape = tuple(mp * 2 for mp in reversed(max_points))
                inflated = spectrum.reshape(inflated_indirect_shape + (n_direct,))

            else:
                raise ValueError(
                    f"spectrum has {spectrum.shape[0]} rows — expected either "
                    f"{n_sampled * n_combos} (collapsed) or "
                    f"{total_max_rows} (full flat) rows"
                )

        else:
            # Already nD shape — just use as-is
            inflated = spectrum

        # Update header
        dim_map = {0: 'FDF1', 1: 'FDF3', 2: 'FDF4'}
        for dim_idx, mp in enumerate(max_points):
            prefix = dim_map[dim_idx]
            inflated_dict[f'{prefix}TDSIZE'] = float(mp)
            inflated_dict[f'{prefix}APOD']   = float(mp)
            inflated_dict[f'{prefix}CENTER'] = float(mp / 2 + 1)

        inflated_dict['FDSPECNUM']  = float(max_points[0] * 2)
        inflated_dict['FDDIMCOUNT'] = float(n_indirect_dims + 1)  # +1 for direct dim
        inflated_dict['FDNUSDIM']   = float(n_indirect_dims)

    return inflated, inflated_dict

def apply_sampling_schedule_nd(fid: NDArray,
                                sampling_schedule: Union[list[int], NDArray],
                                indirect_axes: tuple[int, ...],
                                acq_ord:int = 0) -> NDArray:
    """
    Zero out unsampled points along indirect axes.
    For 1D indirect: sampling_schedule is a list of indices.
    For 2D indirect: sampling_schedule is an array of (i1, i2) pairs.
    """
    result = fid.copy()
    sampling_schedule = np.asarray(sampling_schedule)
    n_indirect = len(indirect_axes)


    if n_indirect == 1:
        # 1D indirect — same as original apply_sampling_schedule_to_signal
        axis = indirect_axes[0]
        n_points = fid.shape[axis]
        all_indices = np.arange(n_points)
        unsampled = all_indices[~np.isin(all_indices, sampling_schedule)]
        idx = [slice(None)] * fid.ndim
        for u in unsampled:
            idx[axis] = u
            result[tuple(idx)] = 0

    elif n_indirect == 2:
        ax0, ax1 = indirect_axes
        n0 = fid.shape[ax0]
        n1 = fid.shape[ax1]

        max_i0 = sampling_schedule[:, 0].max()
        max_i1 = sampling_schedule[:, 1].max()

        if max_i0 >= n0 // 2 or max_i1 >= n1 // 2:
            print('swapping over sampling')
            sampling_schedule = sampling_schedule[:, ::-1]

        if acq_ord == 0:
            expanded = np.concatenate([
                sampling_schedule * 2,
                sampling_schedule * 2 + np.array([1, 0]),
                sampling_schedule * 2 + np.array([0, 1]),
                sampling_schedule * 2 + np.array([1, 1]),
            ], axis=0)

        elif acq_ord == 1:
            expanded = np.concatenate([
                sampling_schedule * 2,
                sampling_schedule * 2 + np.array([0, 1]),
                sampling_schedule * 2 + np.array([1, 0]),
                sampling_schedule * 2 + np.array([1, 1]),
            ], axis=0)

        mask = np.zeros((n0, n1), dtype=bool)
        mask[expanded[:, 0], expanded[:, 1]] = True

        result = np.moveaxis(result, (ax0, ax1), (0, 1))
        result[~mask] = 0
        result = np.moveaxis(result, (0, 1), (ax0, ax1))

    else:
        raise ValueError("Only 1D and 2D indirect dimensions currently supported")

    return result

def write_sched(sampling_schedule: Union[list, NDArray], outfile: str):
    """
    Write a sampling schedule to file.
    
    1D indirect: one integer per line  e.g. '42\n'
    nD indirect: space-separated indices per line e.g. '42 7\n'
    """
    sampling_schedule = np.asarray(sampling_schedule)
    
    with open(outfile, 'w') as f:
        if sampling_schedule.ndim == 1:
            # 1D — single index per line
            for val in sampling_schedule:
                f.write(f'{val}\n')
        else:
            # nD — space-separated indices per line
            for row in sampling_schedule:
                f.write(' '.join(str(v) for v in row) + '\n')

def read_sched(infile: str) -> NDArray:
    """
    Read a sampling schedule from file as written by write_sched.

    1D indirect: one integer per line  → returns shape (n_sampled,)
    nD indirect: space-separated indices per line → returns shape (n_sampled, n_dims)

    Parameters
    ----------
    infile : path to sampling schedule file
    """
    with open(infile, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    # Parse each line into a list of integers
    parsed = [list(map(int, line.split())) for line in lines]

    # Check dimensionality from first line
    n_dims = len(parsed[0])

    if n_dims == 1:
        return np.array([row[0] for row in parsed])  # shape (n_sampled,)
    else:
        return np.array(parsed)  # shape (n_sampled, n_dims)