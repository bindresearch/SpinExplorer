import numpy as np
from numpy.typing import NDArray
from typing import Union, Literal 
from dataclasses import dataclass, field
import matplotlib.pyplot as plt
import random
import math
import itertools
from scipy.stats import truncnorm
from scipy.signal import hilbert
import sys

import pyfftw
import pyfftw.interfaces.numpy_fft as fft


pyfftw.interfaces.cache.enable()
pyfftw.interfaces.cache.set_keepalive_time(30)  # keep plans cached for 30s

def generate_sampling_schedule_1d(num_points:int, sampling:float)->NDArray:
    """
    let's just do random sampling in the first instance
    """

    points_to_sample = int(num_points*sampling)-1
    sampled_points = random.sample(range(1, num_points), points_to_sample)
    sampled_points.append(0)

    return np.sort(sampled_points)

def poisson_func(lmbd:float)->int:
    """
    Generate a Poisson-distributed random integer using Knuth's algorithm.
    
    The Poisson distribution models the number of events occurring in a fixed
    interval, given an average rate (lambda). Knuth's method works by exploiting
    the relationship between the Poisson distribution and the exponential
    distribution: inter-arrival times of a Poisson process are exponentially
    distributed, so we multiply uniform random numbers until their product falls
    below exp(-lambda). The count of multiplications minus one is the sample.
    """
    L = math.exp(-lmbd)  # Threshold: stop when product of uniforms drops below this
    k = 0
    p = 1.0
    while True:
        u = random.random()  # Draw a uniform random number in [0, 1)
        p *= u               # Accumulate the product
        k += 1
        if p < L:            # Product has fallen below threshold — we have our sample
            break
    return k - 1             # Subtract 1 because we overshoot by one iteration

def generate_sampling_schedule_poisson(sampling_rate:float, num_points:int, sine_weighting: int)->list:
    """
    Generate a non-uniform (Poisson-gap) sampling schedule for NMR or similar
    spectroscopy applications.

    Non-uniform sampling (NUS) selects a sparse subset of p points from a full
    grid of z points. Rather than picking points uniformly at random, this
    algorithm uses Poisson-gap sampling: gaps between selected points follow a
    Poisson distribution whose rate varies sinusoidally across the grid. This
    produces a schedule that:
      - Avoids clustering (gaps are never zero, unlike pure random sampling)
      - Has a smoothly varying density (denser near the centre, sparser at edges)
      - Matches a target number of points p exactly via an iterative adjustment

    """
    p = int(sampling_rate*num_points)   # Target number of points to sample (e.g. 64)
    z = num_points    # Total grid size (e.g. 256)

    ld = z / p              # Average spacing between sampled points (decimation factor)
    adj = 2.0 * (ld - 1)   # Initial Poisson rate adjustment; (ld-1) because each
                            # step already advances by 1 before adding the Poisson gap

    v = [0] * z             # Storage for the selected grid indices

    # --- Outer loop: iteratively tune the Poisson rate until exactly p points
    # are selected. If we overshoot, increase adj to widen gaps; if we undershoot,
    # decrease adj to narrow them. Converges quickly (~1% adjustment per iteration).
    while True:
        i = 0  # Current position on the full grid
        n = 0  # Number of points selected so far

        # --- Inner loop: walk across the grid, selecting points with Poisson gaps
        while i < z:
            v[n] = i   # Record this grid point as selected
            i += 1     # Always advance by at least 1 (ensures no repeated indices)

            # Draw a Poisson-distributed gap. The rate varies sinusoidally:
            # sin(...) is low near the grid edges and peaks at the centre,
            # so gaps are smaller (denser sampling) near the centre and larger
            # (sparser sampling) near the edges — a common NUS strategy that
            # concentrates points where the signal decays most rapidly.
            
            if sine_weighting == 1:
                k = poisson_func(adj * math.sin((i + 0.5) / (z + 1) * np.pi/2.0))
            elif sine_weighting == 2:
                k = poisson_func(adj * (math.sin((i + 0.5) / (z + 1) * np.pi/2.0))**2.0)
            else:
                k = poisson_func(adj)

            i += k     # Skip forward by the Poisson gap
            n += 1     # One more point selected

        # --- Adjust the rate parameter and retry if point count is wrong
        if n > p:
            adj *= 1.02   # Too many points selected: widen gaps by 2%
        elif n < p:
            adj /= 1.02   # Too few points selected: narrow gaps by 2%
        else:
            break         # Exactly p points selected — schedule is complete

    return v[:p]


def generate_sampling_schedule_poisson_nd(sampling_rate: float,
                                           num_points: Union[list[int], int],
                                           sine_weighting: int,
                                           tolerance: float = 0.01) -> NDArray:
    if np.isscalar(num_points):
        num_points = [num_points]

    n_dims  = len(num_points)
    n_total = int(np.prod(num_points))
    p       = int(sampling_rate * n_total)

    if p < 1:
        raise ValueError(f"Sampling rate {sampling_rate} too low — "
                         f"would select 0 points from grid of {n_total}")

    z = n_total

    # Pre-compute sine weights for all positions at once
    positions = (np.arange(z) + 0.5) / (z + 1)
    if sine_weighting == 1:
        weights = np.sin(positions * math.pi / 2.0)
    elif sine_weighting == 2:
        weights = np.sin(positions * math.pi / 2.0) ** 2.0
    else:
        weights = np.ones(z)

    # Better initial adj estimate for low sampling rates
    ld  = z / p
    adj = max(2.0 * (ld - 1), (1.0 / sampling_rate) - 1.0)

    # Adaptive step — larger steps for low sampling rates
    adj_step = 1.0 + min(0.1, 5.0 * (1.0 - sampling_rate) / 100.0)

    def walk(adj: float) -> NDArray:
        """
        Vectorised walk: pre-generate all Poisson gaps at once and use
        cumsum to find selected positions — no Python loop needed.
        """
        # Draw more gaps than we'll need (p * safety factor)
        # If we run out we pad with large gaps to stay within z
        n_draw = min(int(p * 2.5), z)

        # Draw all gaps at once using vectorised Poisson
        # scipy.stats.poisson is faster than calling poisson_func in a loop
        from scipy.stats import poisson as scipy_poisson

        rates = adj * weights  # (z,) — rate at each position
        # Sample gaps using mean rates — approximate but fast
        # We use the rates at the first n_draw positions as proxy
        sample_rates = rates[:n_draw]
        # Replace zero rates with small epsilon to avoid degenerate gaps
        sample_rates = np.maximum(sample_rates, 1e-10)
        gaps = scipy_poisson.rvs(sample_rates) + 1  # +1 for mandatory advance

        # Cumulative sum gives selected positions
        positions_selected = np.cumsum(gaps) - gaps[0]  # start from 0
        positions_selected = positions_selected[positions_selected < z]

        return positions_selected

    # Iterative adjustment
    for _ in range(2000):
        selected = walk(adj)
        n = len(selected)

        if n == p:
            break
        elif n > p:
            adj *= adj_step
        else:
            adj /= adj_step
    
    # Trim or pad to exactly p if within tolerance
    selected = selected[:p]

    # Validate
    actual    = len(selected)
    deviation = abs(actual - p) / p
    if deviation > tolerance:
        raise ValueError(
            f"Could not achieve target sampling rate within tolerance. "
            f"Expected {p} points ({sampling_rate*100:.1f}% of {n_total}), "
            f"got {actual} "
            f"({deviation*100:.2f}% deviation, tolerance is {tolerance*100:.1f}%)"
        )

    if n_dims == 1:
        return selected.astype(int)

    # Unravel flat indices to nD coordinates
    nd_indices = np.array(np.unravel_index(selected.astype(int), num_points)).T

    return nd_indices  # shape (p, n_dims)

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

def apply_sampling_schedule_to_2D_signal_ist(spectrum: NDArray, sample_dict: dict, 
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

def apply_sampling_schedule_to_nd_signal_ist(spectrum: NDArray, sample_dict: dict,
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


def inflate_spectra_2D_signal_ist(spectrum: NDArray, sample_dict: dict,
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
                                indirect_axes: tuple[int, ...]) -> NDArray:
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
        if max_i0 >= n0 or max_i1 >= n1:
            print('swapping over sampling')
            sampling_schedule = sampling_schedule[:, ::-1]

        mask = np.zeros((n0, n1), dtype=bool)
        mask[sampling_schedule[:, 0], sampling_schedule[:, 1]] = True

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
    

def get_thresh_signal(signal_ft: NDArray, threshold: float)->NDArray:
    
    max = np.max(signal_ft)
    signal_cutoff = threshold*max
    thresh_sig = signal_ft - signal_cutoff
    thresh_sig[thresh_sig<0.0] = 0.0
    return thresh_sig


def ist_iteration_2d(nus_fid:NDArray, 
                  threshold:float, 
                  sampling_schedule:Union[list[int],NDArray])->tuple[NDArray,NDArray,np.floating]:
    # consider whether we want to put zero-filling into this
    
    signal_ft = np.fft.fft(nus_fid)
    threshold_sig = get_thresh_signal(signal_ft, threshold)
    leftover_sig = signal_ft-threshold_sig
    l2_norm = np.linalg.norm(leftover_sig)
    leftover_fid = np.fft.ifft(leftover_sig)
    leftover_fid = apply_sampling_schedule_nd(leftover_fid, sampling_schedule, (0,))
    
    return leftover_fid, threshold_sig, l2_norm

def ist_iteration_3d(data: NDArray,
                    threshold: float,
                    sampling_schedule: Union[list[int], NDArray]) -> tuple[NDArray, NDArray, np.floating]:

    # FT along faster indirect (axis=1): collapse cos/sin pairs → complex, FT, take real
    data = data[:, 0::2] + 1.0j * data[:, 1::2]        # (n3*2, n2)
    data = fft.fftshift(fft.fft(data, axis=1), axes=1)
    data = np.real(data)                                 # (n3*2, n2)

    # FT along slower indirect (axis=0): collapse cos/sin pairs → complex, FT, take real
    data = data[0::2, :] + 1.0j * data[1::2, :]         # (n3, n2)
    data = fft.fftshift(fft.fft(data, axis=0), axes=0)
    data = np.real(data)                                 # (n3, n2)

    thresh_sig = get_thresh_signal(data, threshold)
    leftover_sig = data - thresh_sig
    l2_norm = np.linalg.norm(leftover_sig)

    # Inverse: re-expand both axes back to interleaved form
    leftover_sig = retrieve_signal_ist_3d(leftover_sig)  # (n3*2, n2*2)
    leftover_sig = apply_sampling_schedule_nd(leftover_sig, sampling_schedule, (0, 1))
    return leftover_sig, thresh_sig, l2_norm

def retrieve_signal_ist_3d(leftover_sig):
    # leftover_sig shape: (n3*2, n2*2) — already interleaved

    # 1. Faster indirect (axis=1): collapse interleave → complex, ifft, re-expand
    leftover_complex1 = leftover_sig[:, 0::2] + 1j * leftover_sig[:, 1::2]
    n2 = leftover_complex1.shape[1]
    leftover_complex1 = fft.ifft(fft.ifftshift(leftover_complex1, axes=1), axis=1)
    new1 = np.zeros_like(leftover_sig, dtype=np.float32)
    new1[:, 0::2] = np.real(leftover_complex1)
    new1[:, 1::2] = np.imag(leftover_complex1)
    leftover_sig = new1

    # 2. Slower indirect (axis=0): collapse interleave → complex, ifft, re-expand
    leftover_complex0 = leftover_sig[0::2, :] + 1j * leftover_sig[1::2, :]
    n3 = leftover_complex0.shape[0]
    leftover_complex0 = fft.ifft(fft.ifftshift(leftover_complex0, axes=0), axis=0)
    new0 = np.zeros_like(leftover_sig, dtype=np.float32)
    new0[0::2, :] = np.real(leftover_complex0)
    new0[1::2, :] = np.imag(leftover_complex0)
    return new0

# def retrieve_signal_ist_3d(leftover_sig):
#     """Unscaled version — used inside IST loop to feed back into next iteration."""
#     leftover_sig = fft.ifft(fft.ifftshift(leftover_sig, axes=-1), axis=-1)
#     leftover_sig_new = np.zeros((leftover_sig.shape[0], leftover_sig.shape[1]*2))
#     leftover_sig_new[...,0::2] = np.real(leftover_sig[...,:])
#     leftover_sig_new[...,1::2] = np.imag(leftover_sig[...,:])
#     leftover_sig = leftover_sig_new

#     leftover_sig = np.transpose(leftover_sig)
#     leftover_sig = fft.ifft(fft.ifftshift(leftover_sig, axes=-1), axis=-1)
#     leftover_sig_new = np.zeros((leftover_sig.shape[0], leftover_sig.shape[1]*2))
#     leftover_sig_new[...,0::2] = np.real(leftover_sig[...,:]).astype(np.float32)
#     leftover_sig_new[...,1::2] = np.imag(leftover_sig[...,:]).astype(np.float32)
#     return leftover_sig_new


def ist_3d(input_spec: NDArray,
           sampling_schedule: Union[list[int], NDArray],
           threshold: float = 0.9,
           terminate: float = 1e-4,
           convergence_tol: float = 1e-10,
           max_iter: int = 400,
           mode: int = 1) -> NDArray:
    
    """
    IST reconstruction of 3D NUS data
    """

    def _reconstruct_until_convergence(nus_fid: NDArray) -> NDArray:
        reconstructed = None
        prev_reconstructed = None

        for iteration in range(1, max_iter + 1):
            nus_fid, threshold_signal, _ = ist_iteration_3d(nus_fid, threshold, sampling_schedule)
            if reconstructed is None:
                reconstructed = np.zeros_like(threshold_signal)
                prev_reconstructed = np.zeros_like(threshold_signal)
            reconstructed += threshold_signal

            change = np.linalg.norm(reconstructed - prev_reconstructed)
            relative_change = change / (np.linalg.norm(reconstructed) + 1e-10)
            prev_reconstructed = reconstructed.copy()

            if relative_change < convergence_tol:
                print(f"  converged at iteration {iteration} — "
                      f"relative change: {relative_change:.2e}")
                break
            if iteration == max_iter:
                print(f"  reached max iterations — "
                      f"relative change: {relative_change:.2e}")

        return reconstructed

    def _reconstruct_until_l2(nus_fid: NDArray) -> NDArray:
        reconstructed = np.zeros((nus_fid.shape[1]//2,nus_fid.shape[0]//2))
        l2_norm = terminate * 1000.0

        for iteration in range(1, max_iter + 1):
            nus_fid, threshold_signal, l2_norm = ist_iteration_3d(nus_fid, threshold, sampling_schedule)
            reconstructed += threshold_signal
            if l2_norm <= terminate:
                print(f"  converged at iteration {iteration} — "
                      f"l2 norm: {l2_norm:.2e}")
                break
            if iteration == max_iter:
                print(f"  reached max iterations — L2 norm: {l2_norm:.2e}")

        return reconstructed

    recon_spec = np.zeros_like(input_spec)
    reconstruct = _reconstruct_until_convergence if mode == 1 else _reconstruct_until_l2


    for i in range(input_spec.shape[0]):
        print(f"IST slice {i + 1} / {input_spec.shape[0]}")
        slice_data = input_spec[i].copy()
        reconstructed = reconstruct(slice_data)
        reconstructed_slice = retrieve_signal_ist_3d(reconstructed)
        recon_spec[i] = reconstructed_slice
    return recon_spec

def ist_2d(input_spec: NDArray,
           sampling_schedule: Union[list[int], NDArray],
           threshold: float = 0.9,
           terminate: float = 0.001,
           convergence_tol: float = 1e-4,
           max_iter: int = 400,
           mode: int = 1) -> NDArray:
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

    def _reconstruct_until_convergence(nus_fid: NDArray) -> NDArray:
        reconstructed = np.zeros_like(nus_fid)
        prev_reconstructed = np.zeros_like(nus_fid)

        for iteration in range(1, max_iter + 1):
            nus_fid, threshold_signal, _ = ist_iteration_2d(nus_fid, threshold, sampling_schedule)
            reconstructed += threshold_signal

            change = np.linalg.norm(reconstructed - prev_reconstructed)
            relative_change = change / (np.linalg.norm(reconstructed) + 1e-10)
            prev_reconstructed = reconstructed.copy()

            if relative_change < convergence_tol:
                print(f"  converged at iteration {iteration} — "
                      f"relative change: {relative_change:.2e}")
                break
            if iteration == max_iter:
                print(f"  reached max iterations — "
                      f"relative change: {relative_change:.2e}")

        return reconstructed

    def _reconstruct_until_l2(nus_fid: NDArray) -> NDArray:
        reconstructed = np.zeros_like(nus_fid)
        l2_norm = terminate * 1000.0

        for iteration in range(1, max_iter + 1):
            nus_fid, threshold_signal, l2_norm = ist_iteration_2d(nus_fid, threshold, sampling_schedule)
            reconstructed += threshold_signal

            if l2_norm <= terminate:
                break
            if iteration == max_iter:
                print(f"  reached max iterations — L2 norm: {l2_norm:.2e}")

        return reconstructed

    recon_spec = np.zeros_like(input_spec)
    reconstruct = _reconstruct_until_convergence if mode == 1 else _reconstruct_until_l2

    for i in range(input_spec.shape[0]):
        print(f"IST slice {i + 1} / {input_spec.shape[0]}")
        reconstructed_slice = reconstruct(input_spec[i].copy())
        recon_spec[i] = fft.ifft(reconstructed_slice)

    return recon_spec

def convert_space_separated_to_csv(infile: str, outfile: str) -> None:
    """
    Convert a space-separated file of numbers (1-4 per line) to CSV format.
    
    Parameters
    ----------
    infile  : path to input space-separated file
    outfile : path to output CSV file
    """
    with open(infile, 'r') as f_in, open(outfile, 'w') as f_out:
        for line in f_in:
            values = line.strip().split()
            if values:  # skip empty lines
                f_out.write(','.join(values) + '\n')
