"""
Combining the pairs of increments of a Rance-Kay (echo-antiecho) dimension.

A sensitivity enhanced dimension is collected as pairs of increments, an echo
and an anti-echo, which are added and subtracted to give the two halves of a
complex point. A spectrum can have more than one such dimension: a
triple-resonance 3D often has both of its indirect dimensions collected this
way, so each of them has to be combined in turn.

This is shared by the conversion in SpinConverter and by the conversion of the
command line tools, so that the two treat quadrature alike.
"""

import numpy as np


def shuffle_rance_kay_axis(data, axis, rotate_phase=True):
    """
    Combine the pairs of increments along one dimension of the data.

    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    axis : int
        Which axis of the data holds the pairs of increments.
    rotate_phase : bool
        True to remove the need for a 90 degree zero-order phase correction.

    Returns
    -------
    ndarray
        The data with that dimension combined.
    """
    values = np.swapaxes(data, 0, axis)
    shuffled = np.empty(values.shape, values.dtype)

    for i in range(0, values.shape[0] - 1, 2):
        first, second = values[i], values[i + 1]

        shuffled[i] = (
            1.0 * (first.real - second.real)
            + 1.0 * (first.imag - second.imag) * 1j
        )

        if rotate_phase is True:
            shuffled[i + 1] = (
                -1.0 * (first.imag + second.imag)
                + 1.0 * (first.real + second.real) * 1j
            )
        else:
            shuffled[i + 1] = (
                1.0 * (first.real + second.real)
                + 1.0 * (first.imag + second.imag) * 1j
            )

    if values.shape[0] % 2 != 0:
        # An odd number of increments leaves the last one without a partner, so
        # it is kept as it is rather than left as whatever the new array held
        shuffled[-1] = values[-1]

    return np.swapaxes(shuffled, 0, axis)


def shuffle_rance_kay(data, axes, rotate_phase=True):
    """
    Combine the pairs of increments of every Rance-Kay dimension of the data,
    one dimension after another. Data with no such dimension is given back
    unchanged.
    """
    for axis in axes:
        data = shuffle_rance_kay_axis(data, axis, rotate_phase=rotate_phase)

    return data
