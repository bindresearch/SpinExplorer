from dataclasses import dataclass, field
from enum import Enum
from typing import Optional
import numpy as np
import nmrglue as ng # type: ignore


class SolventSuppressionFilter(Enum):
    LOW_PASS = (0, "Low-pass")
    SPLINE = (1, "Spline")
    POLYNOMIAL = (2, "Polynomial")

    def __init__(self, number: int, display_name: str):
        self.number = number
        self.display_name = display_name


class LowpassChoices(Enum):
    BOXCAR = (0, "Boxcar")
    SINE = (1, "Sine")
    SINE_SQUARED = (2, "Sine Squared")  # Fixed typo

    def __init__(self, number: int, display_name: str):
        self.number = number
        self.display_name = display_name


class LPPredictedPoints(Enum):
    AFTER_FID = (0, "After FID")
    BEFORE_FID = (1, "Before FID")

    def __init__(self, number: int, display_name: str):
        self.number = number
        self.display_name = display_name


class LPPredictedCoefficients(Enum):  # Fixed typo
    FORWARD = (0, "Forward")
    BACKWARD = (1, "Backward")
    BOTH = (2, "Both")

    def __init__(self, number: int, display_name: str):
        self.number = number
        self.display_name = display_name


class LPChoices(Enum):
    NOTHING = (0, "None")
    LINEAR_PREDICTION = (1, "Linear Prediction")
    NUS_RECONSTRUCTION = (2, "NUS Reconstruction")

    def __init__(self, number: int, display_name: str):
        self.number = number
        self.display_name = display_name


class ApodTypes(Enum):
    NOTHING = (0, "None")
    EXPONENTIAL = (1, "Exponential")
    LORENTZ_TO_GAUSS = (2, "Lorentz to Gauss")
    SINEBELL = (3, "Sinebell")
    GAUSS_BROADENING = (4, "Gauss Broadening")
    TRAPAZOID = (5, "Trapazoid")  
    TRIANGLE = (6, "Triangle")  # Fixed typo

    def __init__(self, number: int, display_name: str):
        self.number = number
        self.display_name = display_name


class ZFTypes(Enum):
    DOUBLING_SPECTRUM_SIZE = (0, "Doubling spectrum size")
    ADDING_ZEROS = (1, "Adding Zeros")
    FINAL_DATA_SIZE = (2, "Final data size")

    def __init__(self, number: int, display_name: str):
        self.number = number
        self.display_name = display_name


class ZFAdditionalParams(Enum):
    DOUBLING_NUMBER = (0, "Doubling number")
    NUM_ZEROS = (1, "Number of zeros to add")
    FINAL_DATA_SIZE = (2, "Final data size")

    def __init__(self, number: int, display_name: str):
        self.number = number
        self.display_name = display_name


class FTOptions(Enum):
    STANDARD = (0, "Standard")
    AUTO = (1, "Auto (not recommended)")
    REAL = (2, "Real")
    INVERSE = (3, "Inverse")
    ALT = (4, "Sign alternation (alt)")
    NEG = (5, "Negative imaginaries (neg)")
    ALT_NEG = (6, "alt + neg")

    def __init__(self, number: int, display_name: str):
        self.number = number
        self.display_name = display_name


class BaselineOptions(Enum):
    LINEAR = (0, "Linear")
    POLYNOMIAL = (1, "Polynomial")

    def __init__(self, number: int, display_name: str):
        self.number = number
        self.display_name = display_name

# TODO: make separate utils area for functions like sol_general and sol_general_nd

def sol_general(data, filter, w=16, mode="same"):
    """
    Solvent filter with generic filter.

    Algorithm described in: Marion et al. JMR 1989 84 425-430

    Parameters
    ----------
    data : 1D or 2D ndarray
        Array of 1D or 2D NMR data.
    filter : ndarray
        Filter to convolve with data.  Not used in solvent filter functions
        which specific the filter, e.g. sol_boxcar.
    w : int, optional
        Filter length.  Not used here but is used in solent filter functions
        which specify the filter, e.g. sol_boxcar.
    mode : {'valid', 'same', 'full'}, optional
        Convolution mode, 'same' should be used.

    Returns
    -------
    ndata : 1D or 2D ndarray
        NMR data with solvent filter applied

    """
    import scipy # type: ignore

    A = filter.sum()
    if data.ndim == 2:
        filter = filter.reshape((1, -1))  # apply along axis=1
    elif data.ndim == 3:
        filter = filter.reshape((1, 1, -1))  # apply along axis=2
    return data - scipy.signal.convolve(data, filter, mode=mode) / A

def sol_general_nd(data, filter, axis=-1, mode="same"):
    """
    Generalized solvent suppression filter for N-D data.

    Applies solvent suppression along a specific axis using convolution.

    Parameters
    ----------
    data : ndarray
        N-D array of NMR data.
    filter : ndarray
        1D filter array to convolve with.
    axis : int, optional
        Axis along which to apply the filter (default: last axis).
    mode : {'valid', 'same', 'full'}, optional
        Convolution mode (usually 'same').

    Returns
    -------
    ndata : ndarray
        Filtered NMR data.
    """
    import numpy as np
    A = filter.sum()
    if A == 0:
        raise ValueError("Filter sum cannot be zero.")

    filtered_data = np.zeros(shape=data.shape)

    # Apply filter to each trace
    import scipy # type: ignore

    if len(data.shape) == 1:
        filtered_data = data - scipy.signal.convolve(data, filter, mode=mode) / A
    elif len(data.shape) == 2:
        for i, dat in enumerate(data):
            filtered_data[i] = (
                dat - scipy.signal.convolve(dat, filter, mode=mode) / A
            )
    elif len(data.shape) == 3:
        for j, dat in enumerate(data):
            for k, dat2 in enumerate(dat):
                filtered_data[j][k] = (
                    dat2 - scipy.signal.convolve(dat2, filter, mode=mode) / A
                )

    return filtered_data


def ext(dic, data, x1, xn, sw):
    """
    Extract a region. Adapted from nmrglue

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    x1 : int or 'default'
        Starting point of the X-axis extraction. 'default' will start the
        extraction at the first point.
    xn : int or 'default'
        Ending point of the X-axis extraction. 'default' will stop the
        extraction at the last point.
    sw : bool
        True to update the sweep width and ppm calibration parameters,
        recommended.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Extracted region of NMR data.

    """

    # store old sizes
    old_x = float(data.shape[-1])

    # slice find limits
    if x1 == "default":
        x_min = 0
    else:
        x_min = np.round(x1) - 1

    if xn == "default":
        x_max = data.shape[-1]
    else:
        x_max = np.round(xn)

    r_x = 1
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))

    # round size to be multiple of r_x when axis is cut
    if x1 != "default" or xn != "default":
        remain_x = (x_min - x_max) % r_x  # -len_x%r_x
        x_min = x_min - np.floor(remain_x / 2)
        x_max = x_max + remain_x - np.floor(remain_x / 2)

    if x_min < 0:
        x_max = x_max - x_min
        x_min = 0.0

    if x_max > data.shape[-1]:
        x_min = x_min - (x_max - data.shape[-1])
        x_max = data.shape[-1]

    no_of_dimensions = len(data.shape)
    if no_of_dimensions == 1:  # 1D Array
        data = data[int(x_min) : int(x_max)]
        dic["FDSIZE"] = x_max - x_min
        dic[fn + "SIZE"] = x_max - x_min
        dic[fn + "FTSIZE"] = x_max - x_min
        dic[fn + "TDSIZE"] = x_max - x_min
        dic[fn + "APODSIZE"] = x_max - x_min
    else:
        data = data[..., int(x_min) : int(x_max)]
        dic["FDSIZE"] = x_max - x_min
        dic[fn + "SIZE"] = x_max - x_min
        dic[fn + "FTSIZE"] = x_max - x_min
        dic[fn + "TDSIZE"] = x_max - x_min
        dic[fn + "APODSIZE"] = x_max - x_min

    # adjust sweep width and ppm calibration
    if sw:
        fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
        s = data.shape[-1]

        if dic[fn + "FTFLAG"] == 0:  # time domain
            dic[fn + "CENTER"] = float(int(s / 2.0 + 1))
            dic[fn + "APOD"] = s
            dic[fn + "TDSIZE"] = s
            dic = ng.pipe_proc.recalc_orig(dic, data, fn)
        else:  # freq domain
            dic[fn + "X1"] = x_min + 1
            dic[fn + "XN"] = x_max
            dic[fn + "APOD"] = np.floor(dic[fn + "APOD"] * s / old_x)
            dic[fn + "CENTER"] = dic[fn + "CENTER"] - x_min
            dic[fn + "SW"] = dic[fn + "SW"] * s / old_x
            dic = ng.pipe_proc.recalc_orig(dic, data, fn)

    dic = ng.pipe_proc.update_minmax(dic, data)
    return dic, data

def lp(
    dic,
    data,
    pred="default",
    x1="default",
    xn="default",
    ord=8,
    mode="f",
    append="after",
    bad_roots="auto",
    mirror=None,
    fix_mode="on",
    method="tls",
):
    """
    Linear Prediction

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    pred : int
        Number of points to predict, "default" chooses the vector size for
        forward prediction, 1 for backward prediction
    x1 : int or 'default'
        First point in 1D vector to use to extract LP filter. 'default' will
        use the first or last point depending on the mode.
    xn : int or 'default'
        Last point in 1D vector to use to extract LP filter. 'default' will use
        the first or last point depending on the mode.
    ord : int
        Prediction order, number of LP coefficients used in prediction.
    mode : {'f', 'b', 'fb'}
        Mode to generate LP filter, 'f' for forward,'b' for backward, 'fb' for
        forward-backward.
    append : {'before' or 'after'}
        Location to append predicted data, 'before' or 'after' the existing
        data.
    bad_roots {'incr', 'decr', None, 'auto'} :
        Type of roots which are will be marked as bad and stabilized. Choices
        are 'incr' for increasing roots, 'decr' for decreasing roots, or None
        for not root stabilization. The default 'auto' will set this parameter
        based upon the LP `mode` parameter: 'f' and 'fb' will results in an
        'incr' parameter. 'b' in 'decr'.
    mirror : {'90-180', '0-0', None}
        Mirror mode, option are '90-180' for a one point shifted mirror image,
        '0-0' for an exact mirror image, and None for no mirror imaging of the
        data.
    fix_mode : {'on', 'reflect'}
        Method used to stabilize bad roots, 'on' moves bad roots onto the unit
        circle, 'reflect' reflect bad roots across the unit circle.
    method : {'svd', 'qr', 'choleskey', 'tls'}
        Method to use to calculate the LP filter.

    Notes
    -----
    The results from this function do not match NMRPipe's LP function.  Also
    some additional parameter and different parameter in this function.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data with linear prediction applied.

    """
    # check parameter
    if mirror not in [None, "90-180", "0-0"]:
        raise ValueError("mirror must be None, '90-180' or '0-0'")

    # pred default values
    if pred == "default":
        if mode == "f" or mode == "fb":
            pred = data.shape[-1]  # double the number of points
        else:
            pred = 1  # predict 1 point before the data

    # remove first pred points if appending before data
    if append == "before":
        data = data[..., pred:]

    # create slice object
    if x1 == "default":
        x_min = 0
    elif mode == "before":
        x_min = x1 - pred - 1
    else:
        x_min = x1 - 1

    if xn == "default":
        x_max = data.shape[-1]
    else:
        x_max = xn - 1
    sl = slice(x_min, x_max)

    # mirror mode (remap to proc_lp names
    mirror = {None: None, "90-180": "180", "0-0": "0"}[mirror]

    # mode, append, bad_roots, fix_mode, and method are passed unchanged
    # use LP-TLS for best results
    data = lp2(
        data, pred, sl, ord, mode, append, bad_roots, fix_mode, mirror, method
    )

    # calculation for dictionary updates
    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
    s = data.shape[-1]
    s2 = s / 2.0 + 1

    # update the dictionary
    dic[fn + "CENTER"] = s2
    if dic["FD2DPHASE"] == 1 and fn != "FDF2":  # TPPI data
        dic[fn + "CENTER"] = np.round(s2 / 2.0 + 0.001)
    dic = ng.pipe_proc.recalc_orig(dic, data, fn)
    dic["FDSIZE"] = s
    dic[fn + "SIZE"] = s
    dic[fn + "APOD"] = s
    dic[fn + "TDSIZE"] = s

    dic = ng.pipe_proc.update_minmax(dic, data)
    return dic, data

def lp2(
        data,
        pred=1,
        slice=slice(None),
        order=8,
        mode="f",
        append="after",
        bad_roots="auto",
        fix_mode="on",
        mirror=None,
        method="svd",
    ):
        """
        Linear prediction extrapolation of 1D or 2D data.

        Parameters
        ----------
        data : ndarray
            1D or 2D NMR data with the last (-1) axis in the time domain.
        pred : int
            Number of points to predict along the last axis.
        slice : slice object, optional
            Slice object which selects the region along the last axis to use in LP
            equation.  The default (slice(None)) will use all points.
        order : int
            Prediction order, number of LP coefficients calculated.
        mode : {'f', 'b', 'fb' or 'bf'}
            Mode to generate LP filter. 'f' for forward,'b' for backward, fb for
            'forward-backward and 'bf' for backward-forward.
        append : {'before', 'after'}
            Location to append the data, either 'before' the current data, or
            'after' the existing data. This is independent of the `mode` parameter.
        bad_roots : {'incr', 'decr', None, 'auto'}
            Type of roots which to consider bad and to stabilize.  Option are those
            with increasing signals 'incr' or decreasing signals 'decr'.  None will
            perform no root stabilizing.  The default ('auto') will set the
            parameter based on the `mode` parameter.  'f' or 'fb' `mode` will
            results in a 'incr' `bad_roots` parameter, 'b' or 'bf` in 'decr'
        fix_mode : {'on', 'reflect'}
            Method used to stabilize bad roots, 'on' to move the roots onto the
            unit circle, 'reflect' to reflect bad roots across the unit circle.
            This parameter is ignored when `bad_roots` is None.
        mirror : {None, '0', '180'}
            Mode to form mirror image of data before processing.  None will
            process the data trace as provided (no mirror image). '0' or '180'
            forms a mirror image of the sliced trace to calculate the LP filter.
            '0' should be used with data with no delay, '180' with data
            with an initial half-point delay.
        method : {'svd', 'qr', 'choleskey', 'tls'}
            Method to use to calculate the LP filter. Choices are a SVD ('svd'), QR
            ('qr'), or Choleskey ('choleskey') decomposition, or Total Least
            Squares ('tls').

        Returns
        -------
        ndata : ndarray
            NMR data with `pred` number of points linear predicted and appended to
            the original data.

        Notes
        -----
        When given 2D data a series of 1D linear predictions are made to
        each row in the array, extending each by pred points. To perform a 2D
        linear prediction using a 2D prediction matrix use :py:func:`lp2d`.

        In forward-backward or backward-forward mode root stabilizing is done
        on both sets of signal roots as calculated in the first mode direction.
        After averaging the coefficient the roots are again stabilized.

        When the append parameter does not match the LP mode, for example
        if a backward linear prediction (mode='b') is used to predict points
        after the trace (append='after'), any root fixing is done before reversing
        the filter.

        """
        if data.ndim == 1:
            return ng.proc_lp.lp_1d(
                data,
                pred,
                slice,
                order,
                mode,
                append,
                bad_roots,
                fix_mode,
                mirror,
                method,
            )
        elif data.ndim == 2:
            # create empty array to hold output
            s = list(data.shape)
            s[-1] = s[-1] + pred
            new = np.empty(s, dtype=data.dtype)
            # vector-wise 1D LP
            for i, trace in enumerate(data):
                new[i] = ng.proc_lp.lp_1d(
                    trace,
                    pred,
                    slice,
                    order,
                    mode,
                    append,
                    bad_roots,
                    fix_mode,
                    mirror,
                    method,
                )
            return new
        else:
            # create empty array to hold output
            s = list(data.shape)
            s[-1] = s[-1] + pred
            new = np.empty(s, dtype=data.dtype)
            # vector-wise 1D LP
            for i, trace in enumerate(data):
                for j, trace2 in enumerate(trace):
                    new[i][j] = ng.proc_lp.lp_1d(
                        trace2,
                        pred,
                        slice,
                        order,
                        mode,
                        append,
                        bad_roots,
                        fix_mode,
                        mirror,
                        method,
                    )
            return new

def base(dic, data, nl=None, nw=0):
    """
    Linear baseline correction.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    nl : list
        List of baseline node points.
    nw : int
        Node width in points.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data with a linear baseline correction applied.

    """

    # change values in node list to start at 0
    nl = [i - 1 for i in nl]

    data = base2(data, nl, nw)
    dic = ng.pipe_proc.update_minmax(dic, data)
    return dic, data


def base2(data, nl, nw=0):
    """
    Linear (first-order) baseline correction based on node list.

    Parameters
    ----------
    data : 1D or 2D ndarray
        Array of 1D or 2D NMR data.
    nl : list
        List of baseline nodes.
    nw : float, optional
        Node half-width in points.

    Returns
    -------
    ndata : ndarray
        NMR data with first order baseline correction applied.  For 2D data
        baseline correction is applied for each trace along the last
        dimension.

    """
    if data.ndim == 1:
        data = data - ng.proc_bl.calc_bl_linear(data, nl, nw)
    elif data.ndim == 2:  # for 2D array loop over traces
        for i, vec in enumerate(data):
            data[i] = data[i] - ng.proc_bl.calc_bl_linear(vec, nl, nw)

    else:
        for i, vec in enumerate(data):
            for j, vec2 in enumerate(vec):
                data[i][j] = data[i][j] - ng.proc_bl.calc_bl_linear(vec2, nl, nw)

    return data



"""
    Obtained from nmrglue followed by customisation
    
    Copyright Notice and Statement for the nmrglue Project
    Copyright (c) 2010-2015 Jonathan J. Helmus
    All rights reserved.
    """

def remove_digital_filter(dic, data, truncate=True):
        """
        Remove the digital filter from Bruker data.

        Parameters
        ----------
        dic : dict
            Dictionary of Bruker parameters.
        data : ndarray
            Array of NMR data to remove digital filter from.
        truncate : bool, optional
            True to truncate the phase shift prior to removing the digital filter.
            This typically produces a better looking spectrum but may remove
            useful data.  False uses a non-truncated phase.
        post_proc : bool, optional
            True if the digital filter is to be removed post processing, i.e after
            fourier transformation. The corrected FID will not be returned, only a
            corrected spectrum in the frequency dimension will be returned

        Returns
        -------
        ndata : ndarray
            Array of NMR data with digital filter removed

        See Also
        ---------
        rm_dig_filter : Remove digital filter by specifying parameters.

        """
        if "acqus" not in dic:
            raise ValueError("dictionary does not contain acqus parameters")

        if "DECIM" not in dic["acqus"]:
            raise ValueError("dictionary does not contain DECIM parameter")
        decim = dic["acqus"]["DECIM"]

        if "DSPFVS" not in dic["acqus"]:
            raise ValueError("dictionary does not contain DSPFVS parameter")
        dspfvs = dic["acqus"]["DSPFVS"]

        if "GRPDLY" not in dic["acqus"]:
            grpdly = 0
        else:
            grpdly = dic["acqus"]["GRPDLY"]

        return rm_dig_filter(data, decim, dspfvs, grpdly, truncate)
    

"""
Obtained from nmrglue followed by customisation

Copyright Notice and Statement for the nmrglue Project
Copyright (c) 2010-2015 Jonathan J. Helmus
All rights reserved.
"""

def rm_dig_filter(data, decim, dspfvs, grpdly=0, truncate_grpdly=True):
        """
        Remove the digital filter from Bruker data.

        Parameters
        ----------
        data : ndarray
            Array of NMR data to remove digital filter from.
        decim : int
            Decimation rate (Bruker DECIM parameter).
        dspfvs : int
            Firmware version (Bruker DSPFVS parameter).
        grpdly : float, optional
            Group delay. (Bruker GRPDLY parameter). When non-zero decim and
            dspfvs are ignored.
        truncate_grpdly : bool, optional
            True to truncate the value of grpdly provided or determined from
            the decim and dspfvs parameters before removing the digital filter.
            This typically produces a better looking spectrum but may remove useful
            data.  False uses a non-truncated grpdly value.
        post_proc : bool, optional
            True if the digital filter is to be removed post processing, i.e after
            fourier transformation. The corrected time domain data will not be
            returned, only the corrected spectrum in the frequency dimension will
            be returned

        Returns
        -------
        ndata : ndarray
            Array of NMR data with digital filter removed.

        See Also
        --------
        remove_digital_filter : Remove digital filter using Bruker dictionary.

        """
        #    A first order phase correction equal to 2*PI*GRPDLY is applied to the
        #    data and the time-corrected FT data is returned

        # The frequency dimension will have the same number of points as the
        # original time domain data, but the time domain data will remain
        # uncorrected
        # -----------------------------------------------------------------------

        if grpdly > 0:  # use group delay value if provided (not 0 or -1)
            phase = grpdly

        # determine the phase correction
        else:
            if dspfvs >= 14:  # DSPFVS greater than 14 give no phase correction.
                phase = 0.0
            else:  # loop up the phase in the table
                bruker_dsp_table = {
                    10: {
                        2: 44.75,
                        3: 33.5,
                        4: 66.625,
                        6: 59.083333333333333,
                        8: 68.5625,
                        12: 60.375,
                        16: 69.53125,
                        24: 61.020833333333333,
                        32: 70.015625,
                        48: 61.34375,
                        64: 70.2578125,
                        96: 61.505208333333333,
                        128: 70.37890625,
                        192: 61.5859375,
                        256: 70.439453125,
                        384: 61.626302083333333,
                        512: 70.4697265625,
                        768: 61.646484375,
                        1024: 70.48486328125,
                        1536: 61.656575520833333,
                        2048: 70.492431640625,
                    },
                    11: {
                        2: 46.0,
                        3: 36.5,
                        4: 48.0,
                        6: 50.166666666666667,
                        8: 53.25,
                        12: 69.5,
                        16: 72.25,
                        24: 70.166666666666667,
                        32: 72.75,
                        48: 70.5,
                        64: 73.0,
                        96: 70.666666666666667,
                        128: 72.5,
                        192: 71.333333333333333,
                        256: 72.25,
                        384: 71.666666666666667,
                        512: 72.125,
                        768: 71.833333333333333,
                        1024: 72.0625,
                        1536: 71.916666666666667,
                        2048: 72.03125,
                    },
                    12: {
                        2: 46.0,
                        3: 36.5,
                        4: 48.0,
                        6: 50.166666666666667,
                        8: 53.25,
                        12: 69.5,
                        16: 71.625,
                        24: 70.166666666666667,
                        32: 72.125,
                        48: 70.5,
                        64: 72.375,
                        96: 70.666666666666667,
                        128: 72.5,
                        192: 71.333333333333333,
                        256: 72.25,
                        384: 71.666666666666667,
                        512: 72.125,
                        768: 71.833333333333333,
                        1024: 72.0625,
                        1536: 71.916666666666667,
                        2048: 72.03125,
                    },
                    13: {
                        2: 2.75,
                        3: 2.8333333333333333,
                        4: 2.875,
                        6: 2.9166666666666667,
                        8: 2.9375,
                        12: 2.9583333333333333,
                        16: 2.96875,
                        24: 2.9791666666666667,
                        32: 2.984375,
                        48: 2.9895833333333333,
                        64: 2.9921875,
                        96: 2.9947916666666667,
                    },
                }
                if dspfvs not in bruker_dsp_table:
                    raise ValueError("dspfvs not in lookup table")
                if decim not in bruker_dsp_table[dspfvs]:
                    raise ValueError("decim not in lookup table")
                phase = bruker_dsp_table[dspfvs][decim]

        if truncate_grpdly:  # truncate the phase
            phase = np.floor(phase)

        s = data.shape[-1]
        pdata = data * np.exp(-2.0j * np.pi * phase * np.arange(s) / s)
        pdata = pdata.astype(data.dtype)
        return pdata






@dataclass
class DimensionConfig:
    """Holds information on dimension processing parameters."""

    # Solvent suppression
    solvent_suppression_flag: bool
    solvent_suppression_filter_length: int
    solvent_suppression_choice: SolventSuppressionFilter
    solvent_suppression_low_pass: LowpassChoices

    # Linear prediction
    lp_flag: bool
    lp_choice: LPChoices
    lp_predicted_point: LPPredictedPoints
    lp_predicted_coefficient: LPPredictedCoefficients  # Fixed typo

    # NUS parameters
    nus_file: str
    nus_extension: int
    nus_cpu: int
    nus_iterations: int

    # Apodization
    apod_flag: bool
    apod_type: ApodTypes
    apod_correction: float  # First point correction

    # Zero filling
    zf_flag: bool
    zf_type: ZFTypes
    zf_additional_param: ZFAdditionalParams
    zf_additional_value: int
    zf_filling_round: bool

    # Fourier transform
    ft_flag: bool
    ft_option: FTOptions

    # Phasing
    ph_flag: bool
    ph_f1180_flag: bool
    ph_p0: float
    ph_p1: float
    ph_magnitude_mode: bool 

    # Extraction
    ex_flag: bool
    ex_start_ppm: float
    ex_end_ppm: float

    # Baseline correction
    bl_flag: bool
    bl_method: BaselineOptions
    bl_params_node_width: float
    bl_params_node_list: list
    bl_params_polynomial_order: int

    # Optional parameters for different apodization types
    # Exponential
    apod_lb: Optional[float] = None

    # Lorentz to Gauss
    apod_g1: Optional[float] = None
    apod_g2: Optional[float] = None
    apod_g3: Optional[float] = None

    # Sinebell
    apod_offset: Optional[float] = None
    apod_end: Optional[float] = None
    apod_power: Optional[int] = None

    # Gaussian broadening (shares apod_lb with exponential)
    apod_gb: Optional[float] = None

    # Trapezoid
    apod_ramp_up: Optional[int] = None
    apod_ramp_down: Optional[int] = None

    # Triangle
    apod_max_loc: Optional[float] = None

    def make_dimension_processing_dictionary(self, dim_num: int) -> dict:
        """
        Create a processing dictionary for this dimension.
        This can be appended to the main parameter dictionary
        and read in by SpinExplorer.
        """
        proc_dic: dict = {}
        proc_dic = self._write_solvent_suppression_dictionary(proc_dic)
        proc_dic = self._write_linear_prediction_dictionary(proc_dic)
        proc_dic = self._write_apodization_dictionary(proc_dic)
        proc_dic = self._write_zf_dictionary(proc_dic)
        proc_dic = self._write_ft_dictionary(proc_dic)
        proc_dic = self._write_phasing_dictionary(proc_dic, dim_num)
        proc_dic = self._write_extraction_parameters(proc_dic)
        proc_dic = self._write_baseline_parameters(proc_dic)
        return proc_dic

    def _write_solvent_suppression_dictionary(self, proc_dic: dict) -> dict:
        """Write solvent suppression part of processing dictionary."""
        proc_dic["Solvent Suppression"] = {
            "Solvent Suppression Flag": self.solvent_suppression_flag,
            "Filter Selection": [
                self.solvent_suppression_choice.number,
                self.solvent_suppression_choice.display_name
            ],
            "Lowpass Shape": [
                self.solvent_suppression_low_pass.number,
                self.solvent_suppression_low_pass.display_name
            ]
        }
        return proc_dic

    def _write_linear_prediction_dictionary(self, proc_dic: dict) -> dict:
        """Write linear prediction part of processing dictionary."""
        proc_dic[self.lp_choice.display_name] = {
            "Linear Prediction Flag": self.lp_flag
        }

        lp_dict = proc_dic[self.lp_choice.display_name]

        if self.lp_choice.number == 1:  # Linear Prediction
            lp_dict["Add predicted points"] = [
                self.lp_predicted_point.number,
                self.lp_predicted_point.display_name
            ]
            lp_dict["Predicted coefficients"] = [  
                self.lp_predicted_coefficient.number,
                self.lp_predicted_coefficient.display_name
            ]
        elif self.lp_choice.number == 2:  # NUS Reconstruction
            lp_dict.update({
                "NUS file": self.nus_file,
                "NUS extension": self.nus_extension,
                "NUS CPU number": self.nus_cpu,
                "NUS iterations": self.nus_iterations
            })

        return proc_dic

    def _write_apodization_dictionary(self, proc_dic: dict) -> dict:
        """Write apodization part of processing dictionary."""
        proc_dic["Apodization"] = {
            "Apodization flag": self.apod_flag,
            "Selection": self.apod_type.number,
            "Type": self.apod_type.display_name,
            "Parameters": {"First point correction": self.apod_correction}
        }

        params = proc_dic["Apodization"]["Parameters"]

        if self.apod_type.number == 1:  # Exponential
            params["Line broadening (Hz)"] = self.apod_lb
        elif self.apod_type.number == 2:  # Lorentz to Gauss
            params.update({
                "Inverse line broadening (Hz)": self.apod_g1,
                "Gaussian broadening (Hz)": self.apod_g2,
                "Gaussian shift": self.apod_g3
            })
        elif self.apod_type.number == 3:  # Sinebell
            params.update({
                "Offset (pi)": self.apod_offset,
                "End (pi)": self.apod_end,
                "Power": self.apod_power
            })
        elif self.apod_type.number == 4:  # Gaussian broadening
            params.update({
                "Line broadening (Hz)": self.apod_lb,
                "Gaussian broadening (Hz)": self.apod_gb
            })
        elif self.apod_type.number == 5:  # Trapezoid
            params.update({
                "Ramp up points": self.apod_ramp_up,
                "Ramp down points": self.apod_ramp_down
            })
        elif self.apod_type.number == 6:  # Triangle
            params["Location of maximum"] = self.apod_max_loc

        return proc_dic

    def _write_zf_dictionary(self, proc_dic: dict) -> dict:
        """Write zero filling part of processing dictionary."""
        proc_dic["Zero filling"] = {
            "Zero filling flag": self.zf_flag,
            "Selection": self.zf_type.number,
            "type": self.zf_type.display_name,  
            "Parameters": {
                self.zf_additional_param.display_name: self.zf_additional_value,
                "Round to nearest power of 2": self.zf_filling_round
            }
        }
        return proc_dic

    def _write_ft_dictionary(self, proc_dic: dict) -> dict:
        """Write Fourier transform part of processing dictionary."""
        proc_dic["Fourier transform"] = {
            "Fourier transform flag": self.ft_flag,
            "Fourier transform method selection": self.ft_option.number,
            "Fourier transform method type": self.ft_option.display_name
        }
        return proc_dic

    def _write_phasing_dictionary(self, proc_dic: dict, dimension: int) -> dict:
        """Write phasing part of processing dictionary."""
        proc_dic["Phasing"] = {
            "Phasing flag": self.ph_flag,
            "P0": self.ph_p0,
            "P1": self.ph_p1
        }

        if dimension == 0:
            proc_dic["Phasing"]["Magnitude mode"] = self.ph_magnitude_mode
        else:
            proc_dic["Phasing"]["f1180 flag"] = self.ph_f1180_flag

        return proc_dic

    def _write_extraction_parameters(self, proc_dic: dict) -> dict:
        """Write extraction part of processing dictionary."""
        proc_dic["Extraction"] = {
            "Extraction flag": self.ex_flag,
            "Start (ppm)": self.ex_start_ppm,
            "End (ppm)": self.ex_end_ppm
        }
        return proc_dic

    def _write_baseline_parameters(self, proc_dic: dict) -> dict:
        """Write baseline correction part of processing dictionary."""
        proc_dic["Baseline correction"] = {
            "Baseline correction flag": self.bl_flag,
            "Selection": self.bl_method.number,
            "Method": self.bl_method.display_name,
            "Parameters": {
                "Node width": self.bl_params_node_width,
                "Node list": self.bl_params_node_list,
                "Polynomial order": self.bl_params_polynomial_order
            }
        }
        return proc_dic
    
    def apply_processing(self, dic, data, dim, filter_removal):
        """
        
        Apply processing steps directly to dic and data.
        Returns the processed dic and data.
        """

        if self.solvent_suppression_flag and dim == 0:
            if self.solvent_suppression_choice == SolventSuppressionFilter.LOW_PASS:
                filter_size = self.solvent_suppression_filter_length
                if self.solvent_suppression_choice == LowpassChoices.BOXCAR:
                    from scipy.signal.windows import boxcar # type: ignore
                    filter = boxcar(filter_size)
                elif self.solvent_suppression_choice == LowpassChoices.SINE:
                    filter = np.cos(np.pi * np.linspace(-0.5, 0.5, filter_size))
                elif self.solvent_suppression_choice == LowpassChoices.SINE_SQUARED:
                    filter = np.cos(np.pi * np.linspace(-0.5, 0.5, filter_size))**2.0
            
                data = sol_general(data, filter, w=filter_size, mode = "same")

            else:
                print('other forms of solvent suppression not supported within nmrglue')

        if self.lp_flag:
            if self.lp_choice != LPChoices.LINEAR_PREDICTION:
                print('only Linear prediction is supported current not NUS reconstruction')
            
            # TODO: streamline this code so that these options can be taken care of within the enum class

            if self.lp_predicted_point == LPPredictedPoints.AFTER_FID:
                append = "after"
            else:
                append = "before"
            
            if self.lp_predicted_coefficient == LPPredictedCoefficients.FORWARD:
                mode = 'f'
            elif self.lp_predicted_coefficient == LPPredictedCoefficients.BACKWARD:
                mode = 'b'
            else:
                mode = 'fb'
            
            dic, data = lp(dic, data, pred = 'default', mode = mode, append = append)

        if self.apod_flag:
            if self.apod_type == ApodTypes.NOTHING:
                dic, data = ng.pipe_proc.em(dic, data, lb=0.0, c=self.apod_correction)
            if self.apod_type == ApodTypes.EXPONENTIAL:
                dic, data = ng.pipe_proc.em(dic, data, lb=self.apod_lb, c=self.apod_correction)
            elif self.apod_type == ApodTypes.LORENTZ_TO_GAUSS:
                dic, data = ng.pipe_proc.gm(dic, data, g1=self.apod_g1, g2=self.apod_g2, 
                                           g3=self.apod_g3, c=self.apod_correction)
            elif self.apod_type == ApodTypes.SINEBELL:
                dic, data = ng.pipe_proc.sp(dic, data, off = self.apod_offset, end = self.apod_end, 
                                            pow = self.apod_power, c = self.apod_correction)
            elif self.apod_type == ApodTypes.GAUSS_BROADENING:
                dic, data = ng.pipe_proc.gmb(dic, data, lb = self.apod_lb, gb=self.apod_gb,
                                              c = self.apod_correction)
            elif self.apod_type == ApodTypes.TRAPAZOID:
                dic, data = ng.pipe_proc.tp(dic, data, t1 = self.apod_ramp_up, t2 = self.apod_ramp_down,
                                            c = self.apod_correction)
            elif self.apod_type == ApodTypes.TRIANGLE:
                dic, data = ng.pipe_proc.tri(dic, data, loc = self.apod_max_loc,
                                             c = self.apod_correction)

        if self.zf_flag:
            if self.zf_type == ZFTypes.DOUBLING_SPECTRUM_SIZE:
                dic, data = ng.pipe_proc.zf(dic, data, zf=self.zf_additional_value, auto = self.zf_filling_round)
            elif self.zf_type == ZFTypes.ADDING_ZEROS:
                dic, data = ng.pipe_proc.zf(dic, data, pad=self.zf_additional_value, auto = self.zf_filling_round)
            elif self.zf_type == ZFTypes.FINAL_DATA_SIZE:
                dic, data = ng.pipe_proc.zf(dic, data, size=self.zf_additional_value, auto = self.zf_filling_round)
        

        if self.ft_flag:
            if self.ft_option == FTOptions.STANDARD:
                dic, data = ng.pipe_proc.ft(dic, data)
            elif self.ft_option == FTOptions.AUTO:
                dic, data = ng.pipe_proc.ft(dic, data, auto=True)
            elif self.ft_option == FTOptions.REAL:
                dic, data = ng.pipe_proc.ft(dic, data, real=True)
            elif self.ft_option == FTOptions.INVERSE:
                dic, data = ng.pipe_proc.ft(dic, data, inverse=True)
            elif self.ft_option == FTOptions.ALT:
                dic, data = ng.pipe_proc.ft(dic, data, alt=True)
            elif self.ft_option == FTOptions.NEG:
                dic, data = ng.pipe_proc.ft(dic, data, neg=True)
            elif self.ft_option == FTOptions.ALT_NEG:
                dic, data = ng.pipe_proc.ft(dic, data, alt = True, neg=True)
            
            if(dim==0 and filter_removal==True):
                dic_bruker, dat_bruker = ng.bruker.read("./")
                data = remove_digital_filter(dic_bruker, data, truncate=False)
        
        if self.ph_flag:
            dic, data = ng.pipe_proc.ps(dic, data, p0=self.ph_p0, p1=self.ph_p1)
        
        if dim == 0 and self.ph_magnitude_mode:
            dic,data = ng.pipe_proc.mc(dic,data)
        
        # Delete imaginaries
        dic, data = ng.pipe_proc.di(dic, data)

        if self.ex_flag:
            # Find the indexes of the ppm values selected
            # Get the ppm values from the data
            ppm_values = ng.pipe.make_uc(dic, data, dim=len(data.shape) - 1)
            ppm_values = ppm_values.ppm_scale()
            x_initial = np.abs(
                ppm_values - self.ex_start_ppm).argmin()
            x_final = np.abs(
                ppm_values - self.ex_end_ppm).argmin()
            if x_initial > x_final:
                x_initial, x_final = x_final, x_initial
            # Change x_initial and x_final so that the difference is an even number
            if (x_final - x_initial + 1) % 2 != 0:
                x_final += 1
            dic, data = ext(dic, data, x1=x_initial, xn=x_final, sw=True)

        if self.bl_flag:
            if self.bl_method == BaselineOptions.POLYNOMIAL:
                print('polynomial baselines not currently supported')
            node_list = (np.asarray(self.bl_params_node_list)/100.)*self.data.shaoe[-1]
            node_list = node_list.astype(int)
            node_list[node_list==0] = self.bl_params_node_width+1

            dic,data = base(dic,data,nl=node_list, nw = self.bl_params_node_width)

        return dic, data


    @classmethod
    def standard_proton(cls, **overrides) -> "DimensionConfig":
        """Create a standard 1H processing configuration."""
        defaults = {
            "solvent_suppression_flag": False,
            "solvent_suppression_filter_length": 0,
            "solvent_suppression_choice": SolventSuppressionFilter.LOW_PASS,
            "solvent_suppression_low_pass": LowpassChoices.BOXCAR,
            "lp_flag": False,
            "lp_choice": LPChoices.LINEAR_PREDICTION,
            "lp_predicted_point": LPPredictedPoints.AFTER_FID,
            "lp_predicted_coefficient": LPPredictedCoefficients.FORWARD,
            "nus_file": "",
            "nus_extension": 0,
            "nus_cpu": 0,
            "nus_iterations": 0,
            "apod_flag": True,
            "apod_type": ApodTypes.EXPONENTIAL,
            "apod_correction": 0.5,
            "apod_lb": 1.0,
            "zf_flag": True,
            "zf_type": ZFTypes.DOUBLING_SPECTRUM_SIZE,
            "zf_additional_param": ZFAdditionalParams.DOUBLING_NUMBER,
            "zf_additional_value": 1,
            "zf_filling_round": True,
            "ft_flag": True,
            "ft_option": FTOptions.STANDARD,
            "ph_flag": True,
            "ph_f1180_flag": False,
            "ph_p0": 0.0,
            "ph_p1": 0.0,
            "ph_magnitude_mode": False,
            "ex_flag": False,
            "ex_start_ppm": -1.0,
            "ex_end_ppm": 13.0,
            "bl_flag": False,
            "bl_method": BaselineOptions.LINEAR,
            "bl_params_node_width": 2,
            "bl_params_node_list": [0, 5, 95, 100],
            "bl_params_polynomial_order": 2,
        }
        defaults.update(overrides)
        return cls(**defaults)
    
    @classmethod
    def standard_nitrogen(cls, **overrides) -> "DimensionConfig":
        """Create a standard 15N processing configuration."""
        defaults = {
            "solvent_suppression_flag": False,
            "solvent_suppression_filter_length": 0,
            "solvent_suppression_choice": SolventSuppressionFilter.LOW_PASS,
            "solvent_suppression_low_pass": LowpassChoices.BOXCAR,
            "lp_flag": False,
            "lp_choice": LPChoices.LINEAR_PREDICTION,
            "lp_predicted_point": LPPredictedPoints.AFTER_FID,
            "lp_predicted_coefficient": LPPredictedCoefficients.FORWARD,
            "nus_file": "",
            "nus_extension": 0,
            "nus_cpu": 0,
            "nus_iterations": 0,
            "apod_flag": True,
            "apod_type": ApodTypes.EXPONENTIAL,
            "apod_correction": 0.5,
            "apod_lb": 1.0,
            "zf_flag": True,
            "zf_type": ZFTypes.DOUBLING_SPECTRUM_SIZE,
            "zf_additional_param": ZFAdditionalParams.DOUBLING_NUMBER,
            "zf_additional_value": 1,
            "zf_filling_round": True,
            "ft_flag": True,
            "ft_option": FTOptions.STANDARD,
            "ph_flag": True,
            "ph_f1180_flag": False,
            "ph_p0": 0.0,
            "ph_p1": 0.0,
            "ph_magnitude_mode": False,
            "ex_flag": False,
            "ex_start_ppm": 90.0,
            "ex_end_ppm": 140.0,
            "bl_flag": False,
            "bl_method": BaselineOptions.LINEAR,
            "bl_params_node_width": 2,
            "bl_params_node_list": [0, 5, 95, 100],
            "bl_params_polynomial_order": 2,
        }
        defaults.update(overrides)
        return cls(**defaults)
    
    @classmethod
    def nitrogen_alt(cls, **overrides) -> "DimensionConfig":
        """Create a 15N processing configuration with the sign alternation FTMODE."""
        defaults = {
            "solvent_suppression_flag": False,
            "solvent_suppression_filter_length": 0,
            "solvent_suppression_choice": SolventSuppressionFilter.LOW_PASS,
            "solvent_suppression_low_pass": LowpassChoices.BOXCAR,
            "lp_flag": False,
            "lp_choice": LPChoices.LINEAR_PREDICTION,
            "lp_predicted_point": LPPredictedPoints.AFTER_FID,
            "lp_predicted_coefficient": LPPredictedCoefficients.FORWARD,
            "nus_file": "",
            "nus_extension": 0,
            "nus_cpu": 0,
            "nus_iterations": 0,
            "apod_flag": True,
            "apod_type": ApodTypes.EXPONENTIAL,
            "apod_correction": 0.5,
            "apod_lb": 1.0,
            "zf_flag": True,
            "zf_type": ZFTypes.DOUBLING_SPECTRUM_SIZE,
            "zf_additional_param": ZFAdditionalParams.DOUBLING_NUMBER,
            "zf_additional_value": 1,
            "zf_filling_round": True,
            "ft_flag": True,
            "ft_option": FTOptions.ALT,
            "ph_flag": True,
            "ph_f1180_flag": False,
            "ph_p0": 0.0,
            "ph_p1": 0.0,
            "ph_magnitude_mode": False,
            "ex_flag": False,
            "ex_start_ppm": 90.0,
            "ex_end_ppm": 140.0,
            "bl_flag": False,
            "bl_method": BaselineOptions.LINEAR,
            "bl_params_node_width": 2,
            "bl_params_node_list": [0, 5, 95, 100],
            "bl_params_polynomial_order": 2,
        }
        defaults.update(overrides)
        return cls(**defaults)
    

    @classmethod
    def nus_nitrogen(cls, **overrides) -> "DimensionConfig":
        """Create a nus 15N processing configuration."""
        defaults = {
            "solvent_suppression_flag": False,
            "solvent_suppression_filter_length": 0,
            "solvent_suppression_choice": SolventSuppressionFilter.LOW_PASS,
            "solvent_suppression_low_pass": LowpassChoices.BOXCAR,
            "lp_flag": False,
            "lp_choice": LPChoices.NUS_RECONSTRUCTION,
            "lp_predicted_point": LPPredictedPoints.AFTER_FID,
            "lp_predicted_coefficient": LPPredictedCoefficients.FORWARD,
            "nus_file": "nuslist",
            "nus_extension": 0,
            "nus_cpu": 2,
            "nus_iterations": 50,
            "apod_flag": True,
            "apod_type": ApodTypes.EXPONENTIAL,
            "apod_correction": 0.5,
            "apod_lb": 1.0,
            "zf_flag": True,
            "zf_type": ZFTypes.DOUBLING_SPECTRUM_SIZE,
            "zf_additional_param": ZFAdditionalParams.DOUBLING_NUMBER,
            "zf_additional_value": 1,
            "zf_filling_round": True,
            "ft_flag": True,
            "ft_option": FTOptions.STANDARD,
            "ph_flag": True,
            "ph_f1180_flag": False,
            "ph_p0": 0.0,
            "ph_p1": 0.0,
            "ph_magnitude_mode": False,
            "ex_flag": False,
            "ex_start_ppm": 90.0,
            "ex_end_ppm": 140.0,
            "bl_flag": False,
            "bl_method": BaselineOptions.LINEAR,
            "bl_params_node_width": 2,
            "bl_params_node_list": [0, 5, 95, 100],
            "bl_params_polynomial_order": 2,
        }
        defaults.update(overrides)
        return cls(**defaults)
    

    @classmethod
    def standard_carbon(cls, **overrides) -> "DimensionConfig":
        """Create a standard 13C processing configuration."""
        defaults = {
            "solvent_suppression_flag": False,
            "solvent_suppression_filter_length": 0,
            "solvent_suppression_choice": SolventSuppressionFilter.LOW_PASS,
            "solvent_suppression_low_pass": LowpassChoices.BOXCAR,
            "lp_flag": False,
            "lp_choice": LPChoices.LINEAR_PREDICTION,
            "lp_predicted_point": LPPredictedPoints.AFTER_FID,
            "lp_predicted_coefficient": LPPredictedCoefficients.FORWARD,
            "nus_file": "",
            "nus_extension": 0,
            "nus_cpu": 0,
            "nus_iterations": 0,
            "apod_flag": True,
            "apod_type": ApodTypes.EXPONENTIAL,
            "apod_correction": 0.5,
            "apod_lb": 1.0,
            "zf_flag": True,
            "zf_type": ZFTypes.DOUBLING_SPECTRUM_SIZE,
            "zf_additional_param": ZFAdditionalParams.DOUBLING_NUMBER,
            "zf_additional_value": 1,
            "zf_filling_round": True,
            "ft_flag": True,
            "ft_option": FTOptions.STANDARD,
            "ph_flag": True,
            "ph_f1180_flag": False,
            "ph_p0": 0.0,
            "ph_p1": 0.0,
            "ph_magnitude_mode": False,
            "ex_flag": False,
            "ex_start_ppm": 90.0,
            "ex_end_ppm": 140.0,
            "bl_flag": False,
            "bl_method": BaselineOptions.LINEAR,
            "bl_params_node_width": 2,
            "bl_params_node_list": [0, 5, 95, 100],
            "bl_params_polynomial_order": 2,
        }
        defaults.update(overrides)
        return cls(**defaults)
    

    @classmethod
    def nus_carbon(cls, **overrides) -> "DimensionConfig":
        """Create a nus 13C processing configuration."""
        defaults = {
            "solvent_suppression_flag": False,
            "solvent_suppression_filter_length": 0,
            "solvent_suppression_choice": SolventSuppressionFilter.LOW_PASS,
            "solvent_suppression_low_pass": LowpassChoices.BOXCAR,
            "lp_flag": False,
            "lp_choice": LPChoices.NUS_RECONSTRUCTION,
            "lp_predicted_point": LPPredictedPoints.AFTER_FID,
            "lp_predicted_coefficient": LPPredictedCoefficients.FORWARD,
            "nus_file": "nuslist",
            "nus_extension": 0,
            "nus_cpu": 2,
            "nus_iterations": 50,
            "apod_flag": True,
            "apod_type": ApodTypes.EXPONENTIAL,
            "apod_correction": 0.5,
            "apod_lb": 1.0,
            "zf_flag": True,
            "zf_type": ZFTypes.DOUBLING_SPECTRUM_SIZE,
            "zf_additional_param": ZFAdditionalParams.DOUBLING_NUMBER,
            "zf_additional_value": 1,
            "zf_filling_round": True,
            "ft_flag": True,
            "ft_option": FTOptions.STANDARD,
            "ph_flag": True,
            "ph_f1180_flag": False,
            "ph_p0": 0.0,
            "ph_p1": 0.0,
            "ph_magnitude_mode": False,
            "ex_flag": False,
            "ex_start_ppm": 90.0,
            "ex_end_ppm": 140.0,
            "bl_flag": False,
            "bl_method": BaselineOptions.LINEAR,
            "bl_params_node_width": 2,
            "bl_params_node_list": [0, 5, 95, 100],
            "bl_params_polynomial_order": 2,
        }
        defaults.update(overrides)
        return cls(**defaults)


@dataclass
class ExperimentConfigStore:
    """
    Stores default parameters for one-step processing of NMR spectra with SpinExplorer.
    """
    dim_labels: list[str]
    dim_configs: list[DimensionConfig]
    fid_name: str 
    ft_name: str

    def __post_init__(self):
        """Validate that dim_labels and dim_configs have matching lengths."""
        if len(self.dim_labels) != len(self.dim_configs):
            raise ValueError(
                f"dim_labels length ({len(self.dim_labels)}) must match "
                f"dim_configs length ({len(self.dim_configs)})"
            )

    def make_processing_dictionary(self) -> dict:
        """Create the complete processing dictionary for all dimensions."""
        proc_dic = {"processing": {}}
        
        for i, label in enumerate(self.dim_labels):
            proc_dic["processing"][label] = (
                self.dim_configs[i].make_dimension_processing_dictionary(i)
            )
        
        return proc_dic
    
    def write_processing_dic(self, proc_dic)->None:
        import json

        with open('parameters.json','w') as outy:
            json.dump(proc_dic, outy, indent = 4)

    def set_pseudo_flag(self, flag: bool)->None:
        """
        Setting the pseudo flag to either true or false
        (required to get the right axis labels for nmrglue processing)
        """
        self.pseudo_flag = flag



    def process_data_pipe(self):
        """
        Process NMR data using nmrpipe (for data requiring
        NUS reconstruction)
        """

    
    def process_data(self, filter_removal=False, pseudo_flag=False):
        """
        Process NMR data directly using the configuration.
        
        Args:
            input_file: Path to input .fid file
            output_file: Optional path to save processed data
        
        Returns:
            Processed dic and data
        """

        # Load data
        dic, data = ng.pipe.read(self.fid_name)

        # Set the comment to nmrglue so that the fact nmrglue processing was used is noted in the processed spectrum header
        dic['FDCOMMENT'] = 'nmrglue'
        
        try:
            if(pseudo_flag==True):
                # Adding the fact that there is a pseudo axis to the FDCOMMENT
                dic['FDCOMMENT'] += '_pseudo'
        except:
            pass

        
        # Process each dimension
        for i, (label, config) in enumerate(zip(self.dim_labels, self.dim_configs)):
            print(f"Processing dimension {i} ({label})...")
            
            # Apply processing for this dimension
            dic, data = config.apply_processing(dic, data, i, filter_removal==False)

            # Transpose if not the last dimension
            if i < len(self.dim_configs) - 1:
                dic, data = ng.pipe_proc.tp(dic, data)
        
        # Save if output file specified
        if self.ft_name:
            ng.pipe.write(self.ft_name, dic, data, overwrite=True)
            print(f"Saved processed data to {self.ft_name}")
        
        return dic, data
    