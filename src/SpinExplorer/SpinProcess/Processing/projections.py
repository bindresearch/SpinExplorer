"""
Writing the projections of a 3D spectrum.

The bore view of a 3D shows a plane which is read from a projection file, so a
3D needs one for each of its three planes: a HNCO has H-N, H-CO and N-CO. They
are skyline projections, the largest and the smallest intensity along the axis
being collapsed added together, so that peaks of either sign are kept.

They are named after the two dimensions they hold, which is how the viewers find
them.

This is shared by the processing in SpinProcess and by the automatic processing,
so that the two write the same projections.
"""

import copy
import os
import shutil

import numpy as np
import nmrglue as ng  # type: ignore

from SpinExplorer.SpinProcess.Processing import transposes


def find_projection_name(dic) -> str:
    """
    The file name for a projection. Once the first axis of the data has
    been collapsed, the remaining axes of the plane are FDDIMORDER[1]
    (rows) and FDDIMORDER[0] (columns). The labels are written in that
    order so that they match the data when the projection is displayed.
    """
    rows = "FDF" + str(int(dic["FDDIMORDER"][1])) + "LABEL"
    columns = "FDF" + str(int(dic["FDDIMORDER"][0])) + "LABEL"

    return dic[rows] + "." + dic[columns] + ".dat"


def move_old_projections() -> None:
    """
    Put the projections which are already here into a folder of their own, so
    that the ones which are about to be written are the only ones found.
    """
    current_dir = os.getcwd()
    old_dir = os.path.join(current_dir, "OldProjections")

    # Create 'Old' directory if it doesn't exist
    os.makedirs(old_dir, exist_ok=True)

    # Loop through files in the current directory
    for filename in os.listdir(current_dir):
        if filename.endswith(".dat") and os.path.isfile(filename):
            source = os.path.join(current_dir, filename)
            destination = os.path.join(old_dir, filename)
            shutil.move(source, destination)


def write_projection(dic, data) -> str:
    """
    Collapse the first axis of the data and write the plane which is left as a
    projection, named after the two dimensions it holds.
    """
    plane = np.max(data, axis=0) + np.min(data, axis=0)

    plane_dic = copy.deepcopy(dic)

    fn = "FDF" + str(int(dic["FDDIMORDER"][2]))
    plane_dic["FDDIMCOUNT"] = 2
    for key in ["SIZE", "TDSIZE", "FTSIZE", "APOD", "APODSIZE", "SW", "CENTER"]:
        plane_dic[fn + key] = 0

    name = find_projection_name(plane_dic)
    ng.pipe.write(name, plane_dic, plane, overwrite=True)

    return name


def write_3d_projections(dic, data) -> list:
    """
    Write a projection of each of the three planes of a 3D spectrum, and say
    which files were written.
    """
    move_old_projections()

    written = [write_projection(dic, data)]

    # The other two planes are reached by turning the data so that a different
    # dimension is the one being collapsed
    first_dic, first_data = transposes.zero_transpose_3d(copy.deepcopy(dic), data)
    written.append(write_projection(first_dic, first_data))

    second_dic, second_data = transposes.transpose_3d(copy.deepcopy(dic), data)
    second_dic, second_data = transposes.zero_transpose_3d(second_dic, second_data)
    written.append(write_projection(second_dic, second_data))

    return written
