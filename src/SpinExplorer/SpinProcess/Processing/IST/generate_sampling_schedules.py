import numpy as np
from numpy.typing import NDArray
import random
import math
from typing import Union

def generate_sampling_schedule_1d(num_points:int, sampling:float)->NDArray:
    """
    let's just do random sampling in the first instance
    """

    points_to_sample = int(num_points*sampling)-1
    sampled_points = random.sample(range(1, num_points), points_to_sample)
    sampled_points.append(0)

    return np.sort(sampled_points)


def generate_random_sampling_schedule_2d(
    num_points_indirect: tuple[int, int],
    sampling: float
) -> NDArray:
    """
    Random 2D NUS schedule with no repeated points, always including (0,0).
    
    Args:
        num_points_indirect: (n1, n2) — grid dimensions for the two indirect dims
        sampling: fractional sampling density (e.g. 0.25 for 25%)
    
    Returns:
        Sorted (M, 2) array of [t1, t2] index pairs
    """
    n1, n2 = num_points_indirect
    total_grid = n1 * n2
    points_to_sample = int(total_grid * sampling) - 1  # reserve one slot for (0,0)

    # All grid points except (0,0)
    all_points = [(i, j) for i in range(n1) for j in range(n2) if (i, j) != (0, 0)]
    sampled = random.sample(all_points, points_to_sample)
    sampled.append((0, 0))

    # Sort lexicographically: t1 first, then t2
    sampled.sort()
    return np.array(sampled)


def generate_random_sampling_schedule_nd(
    num_points_indirect: tuple[int, ...],
    sampling: float
) -> NDArray:
    """
    Random nD NUS schedule with no repeated points, always including the origin.

    Args:
        num_points_indirect: grid dimensions for each indirect dimension
        sampling: fractional sampling density (e.g. 0.25 for 25%)

    Returns:
        Sorted (M, n) array of index tuples, one row per sampled point
    """
    from itertools import product

    total_grid = math.prod(num_points_indirect)
    points_to_sample = int(total_grid * sampling) - 1  # reserve slot for origin

    origin = tuple(0 for _ in num_points_indirect)

    # All grid points except the origin
    all_points = [p for p in product(*[range(n) for n in num_points_indirect]) if p != origin]
    sampled = random.sample(all_points, points_to_sample)
    sampled.append(origin)

    sampled.sort()
    return np.array(sampled)



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



