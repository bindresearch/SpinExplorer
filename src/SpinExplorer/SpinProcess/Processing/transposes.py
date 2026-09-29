"""
Transposing 3D NMRPipe data.

A dimension is processed while it is the last axis of the data, so the data is
transposed between dimensions to bring each one into that place in turn. nmrglue
only transposes 2D data: its pipe_proc.tp works on the last two axes and its
ztp, which would exchange the first and last axes of a 3D, is not implemented.
These two functions do it for 3D data, updating the NMRPipe header along with
the data so that the two keep agreeing with one another.

They were moved here so that the processing in SpinProcess and the automatic
processing of the command line tools use the same code rather than each holding
their own version.
"""

import copy

import numpy as np
import nmrglue as ng  # type: ignore


def zero_transpose_3d(dic, data, nohyper=False):
    """
    Transpose NMRPipe-style data from (X, Y, Z) to (Z, Y, X),
    including correct updates to the NMRPipe dictionary.

    Parameters:
        dic (dict): NMRPipe dictionary
        data (ndarray): NMRPipe data, assumed shape (X, Y, Z)

    Returns:
        new_dic (dict): Transposed dictionary
        new_data (ndarray): Transposed data, shape (Z, Y, X)
    """
    # Transpose data from (X, Y, Z) to (Z, Y, X)
    new_data = np.transpose(data, axes=(2, 1, 0))

    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
    fn3 = "FDF" + str(int(dic["FDDIMORDER"][2]))  # F1, F2, etc

    # Create new dictionary
    new_dic = copy.deepcopy(dic)

    # swapping the FDDIMORDER1 and FDDIMORDER3 values
    order1 = dic["FDDIMORDER1"]
    order3 = dic["FDDIMORDER3"]
    new_dic["FDDIMORDER1"] = order3
    new_dic["FDDIMORDER3"] = order1

    new_dic["FDDIMORDER"] = [
        new_dic["FDDIMORDER1"],
        new_dic["FDDIMORDER2"],
        new_dic["FDDIMORDER3"],
        new_dic["FDDIMORDER4"],
    ]

    new_dic["FDSLICECOUNT"] = new_data.shape[-2]
    new_dic["FDSPECNUM"] = new_dic["FDSLICECOUNT"]
    new_dic["FDSIZE"] = new_data.shape[-1]

    if(nohyper==False):
        if dic[fn3 + "QUADFLAG"] != 1:
            # unpack complex as needed
            new_data = np.array(ng.proc_base.c2ri(new_data), dtype="complex64")
            if fn3 + "SIZE" in new_dic:
                # F1/F2 sizes are not stored in the nmrPipe header
                new_dic[fn3 + "SIZE"] = int(new_dic[fn3 + "SIZE"] / 2)

    return new_dic, new_data


"""
Obtained from nmrglue followed by customisation

Copyright Notice and Statement for the nmrglue Project
Copyright (c) 2010-2015 Jonathan J. Helmus
All rights reserved.
"""


def transpose_3d(dic, data, hyper=False, nohyper=False, auto=False, nohdr=False):
    """
    Exchange the last two axes of 3D data, so that the second dimension becomes
    the one which is processed.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters.
    data : ndarray
        Array of NMR data.
    hyper : bool
        True to perform hypercomplex transpose.
    nohyper : bool
        True to suppress hypercomplex transpose.
    auto : bool
        True to choose transpose mode automatically.
    nohdr : bool
        True to not update the transpose parameters in ndic.

    Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data which has been transposed.

    """
    # XXX test if works with TPPI
    if nohyper:
        hyper = False

    fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
    fn2 = "FDF" + str(int(dic["FDDIMORDER"][1]))  # F1, F2, etc

    if auto:
        if (dic[fn + "QUADFLAG"] != 1) and (dic[fn2 + "QUADFLAG"] != 1):
            hyper = True
        else:
            hyper = False

    if hyper:  # Hypercomplex transpose need type recast
        data = np.array(ng.proc_base.tp_hyper(data), dtype="complex64")
    else:
        data = np.transpose(data, axes=(0, 2, 1))
        if dic[fn2 + "QUADFLAG"] != 1 and nohyper is False:
            # unpack complex as needed
            data = np.array(ng.proc_base.c2ri(data), dtype="complex64")

    # update the dimensionality and order
    dic["FDSLICECOUNT"] = data.shape[-2]
    if (data.dtype == "float32") and (nohyper is True):
        # when nohyper is True and the new last dimension was complex
        # prior to transposing then FDSIZE is set as if the dimension was
        # converted to complex data, that is half the actual size.
        dic["FDSIZE"] = data.shape[-1] / 2
    else:
        dic["FDSIZE"] = data.shape[-1]

    dic["FDSPECNUM"] = dic["FDSLICECOUNT"]
    dic["FDDIMORDER1"], dic["FDDIMORDER2"] = (
        dic["FDDIMORDER2"],
        dic["FDDIMORDER1"],
    )
    dic["FDDIMORDER"] = [
        dic["FDDIMORDER1"],
        dic["FDDIMORDER2"],
        dic["FDDIMORDER3"],
        dic["FDDIMORDER4"],
    ]

    if dic["FD2DPHASE"] == 0:
        dic["FDF1QUADFLAG"], dic["FDF2QUADFLAG"] = (
            dic["FDF2QUADFLAG"],
            dic["FDF1QUADFLAG"],
        )

    if nohdr is not True:
        dic["FDTRANSPOSED"] = (dic["FDTRANSPOSED"] + 1) % 2

    dic = ng.pipe_proc.clean_minmax(dic)
    return dic, data
