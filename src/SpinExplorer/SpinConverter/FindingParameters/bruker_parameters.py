#!/usr/bin/env python3

"""MIT License

Copyright (c) 2025 James Eaton, Andrew Baldwin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE."""


import wx
import nmrglue as ng
import os
import warnings

warnings.simplefilter("ignore", UserWarning)
from typing import List


class ParameterExtractorBruker:
    def __init__(self, nmrdata) -> None:
        """
        This class contains functions to search through the bruker acqus files
        to find relevant parameters needed for conversion to nmrPipe format.
        """
        self.nmrdata = nmrdata

        # Creating a hidden frame to be used as a parent for popout messages
        self.tempframe = wx.Frame(None, title="Temporary Parent", size=(1, 1))
        self.tempframe.Hide()  # Hide the frame since we don't need it to be visible

        self.acqus_file = open(self.nmrdata.parameter_file, "r")
        self.acqus_file_lines = self.acqus_file.readlines()
        self.acqus_file.close()
        self.udic = ng.bruker.guess_udic(self.nmrdata.nmr_dic, self.nmrdata.nmr_data)

    def find_size_bruker(self) -> None:
        """
        Search through Bruker acqus file to search for NMR spectrum sizes.
        The TD parameter denotes the direct dimension size and TD_INDIRECT
        for the indirect dimension sizes (at least as a first guess).
        Sometimes this varies depending on pulse sequence.
        """
        # Look for TD in the acqus file
        self.size_direct, self.size_direct_complex = self.find_direct_bruker()
        # Create an array of initial guess for the indirect dimension sizes
        self.size_indirect = self.find_indirect_bruker()
        # Checking that the direct dimension TD value agrees with the data size
        self.checking_direct_bruker()
        # Checking the indirect dimension
        self.checking_indirect_bruker()

    def find_direct_bruker(self) -> List[int]:
        """
        Reading through Bruker acqus file lines to find the TD entry for
        direct dimension size.
        """
        for i in range(len(self.acqus_file_lines)):
            if "##$TD=" in self.acqus_file_lines[i]:
                line = self.acqus_file_lines[i].split()
                size_direct = int(line[1])
                size_direct_complex = int(size_direct)
                break

        return size_direct, size_direct_complex

    def checking_direct_bruker(self) -> None:
        try:
            self.size_direct
            # Find the max value in nmr_data.shape
            if len(self.nmrdata.nmr_data.shape) == 1:
                self.size_1 = max(self.nmrdata.nmr_data.shape)
                self.size_2 = self.size_1
                for i in range(len(self.size_indirect)):
                    self.size_2 = self.size_2 / self.size_indirect[i]
                self.size_direct_complex = int(self.size_2 * 2)
            else:
                # Sometimes have the issue where stored complex data size is larger than TD, this ensures that the direct dimension
                # size is altered to the larger size to correctly be read
                self.size_1 = self.nmrdata.nmr_data.shape[-1]
                self.size_2 = self.size_1
                if self.size_1 * 2 > self.size_direct:
                    self.size_direct_complex = int(self.size_1 * 2)
                self.size_direct = self.size_1 * 2

        except:
            dlg = wx.MessageDialog(
                self.tempframe,
                "Error in finding data dimension sizes. Please check the acqus file and try again.",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            self.tempframe.Raise()
            self.tempframe.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            exit()

    def find_indirect_bruker(self) -> List:
        """
        Looking for initial guesses for the size of the Bruker indirect dimensions
        """
        size_indirect = []
        found_indirect_sizes = False
        for i in range(len(self.acqus_file_lines)):
            if "##$TD_INDIRECT" in self.acqus_file_lines[i]:
                found_indirect_sizes = True
                continue
            if found_indirect_sizes == True:
                line = self.acqus_file_lines[i].split()
                for j in range(len(line)):
                    if line[j] == "0":
                        pass
                    else:
                        size_indirect.append(int(line[j]))
                break

        return size_indirect

    def checking_indirect_bruker(self) -> None:
        """
        If have not found TD_INDIRECT values, check through acqu2s and acqu3s
        so see if there are any values present. If no other values can be
        found, the spectrum will be assumed to be 1D.
        """
        # Look to see if there is an acqu2s/acqu3s file
        if len(self.size_indirect) == 0:
            # Checking other spectrometer files if TD_INDIRECT is empty
            if "acqu2s" in os.listdir("./"):
                self.sizes_dim2 = []
                self.sizes_dim2_nus = []
                # Look for TD in the acqu2s file
                with open("acqu2s", "r") as file:
                    file_lines = file.readlines()
                    for i in range(len(file_lines)):
                        if "##$NUSTD=" in file_lines[i]:
                            line = file_lines[i].split("\n")[0].split()
                            self.sizes_dim2_nus = int(line[1])
                            break
                        if "##$TD=" in file_lines[i]:
                            line = file_lines[i].split()
                            self.sizes_dim2 = int(line[1])
                            break

                if self.sizes_dim2_nus != []:
                    self.size_indirect.append(self.sizes_dim2_nus)
                else:
                    self.size_indirect.append(self.sizes_dim2)

            if "acqu3s" in os.listdir("./"):
                self.sizes_dim3 = []
                self.sizes_dim3_nus = []
                # Look for TD in the acqu3s file
                with open("acqu3s", "r") as file:
                    file_lines = file.readlines()
                    for i in range(len(file_lines)):
                        if "##$NUSTD=" in file_lines[i]:
                            line = file_lines[i].split("\n")[0].split()
                            self.sizes_dim3_nus = int(line[1])
                            break
                        if "##$TD=" in file_lines[i]:
                            line = file_lines[i].split()
                            self.sizes_dim3 = int(line[1])
                            break

                if self.sizes_dim3_nus != []:
                    self.size_indirect.append(self.sizes_dim3_nus)
                else:
                    self.size_indirect.append(self.sizes_dim3)

        try:
            self.size_indirect
        except:
            dlg = wx.MessageDialog(
                self.tempframe,
                "Error: TD_INDIRECT not found in acqus file. Unable to determine size of data for indirect dimension. Unable to convert data to NMRPipe format. Please check the acqus file and try again.",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            self.tempframe.Raise()
            self.tempframe.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            exit()

        if self.size_indirect != [] and "acqu2s" not in os.listdir("./"):
            self.size_indirect = []

        # Remove values from size indirect if they are equal to 1
        sizes_new = []
        for size in self.size_indirect:
            if size != 1:
                sizes_new.append(size)

        self.size_indirect = sizes_new

        # Try to go through acqu2s and acqu3s and find the nucleus labels and corresponding TD values
        self.indirect_sizes_dict = {}
        if "acqu2s" in os.listdir("./"):
            # Look for TD in the acqu2s file
            nus = False
            with open("acqu2s", "r") as file:
                file_lines = file.readlines()

                for i in range(len(file_lines)):
                    if "##$NUSTD=" in file_lines[i] or "##$NusTD=" in file_lines[i]:
                        nus = True
                        line = file_lines[i].split("\n")[0].split()
                        self.sizes_dim2_nus = int(line[1])
                    if "##$TD=" in file_lines[i]:
                        line = file_lines[i].split()
                        self.sizes_dim2 = int(line[1])
                    if "##$NUC1=" in file_lines[i]:
                        line = file_lines[i].split("\n")[0]
                        nuc = line.split("<")[1].split(">")[0]

            try:
                if nus == False:
                    self.indirect_sizes_dict[nuc] = self.sizes_dim2
                else:
                    self.indirect_sizes_dict[nuc] = self.sizes_dim2_nus
            except:
                pass

        if "acqu3s" in os.listdir("./"):
            # Look for TD in the acqu2s file
            nus = False
            with open("acqu3s", "r") as file:
                file_lines = file.readlines()
                for i in range(len(file_lines)):
                    if "##$NUSTD=" in file_lines[i] or "##$NusTD=" in file_lines[i]:
                        nus = True
                        line = file_lines[i].split("\n")[0].split()
                        self.sizes_dim3_nus = int(line[1])
                    if "##$TD=" in file_lines[i]:
                        line = file_lines[i].split()
                        self.sizes_dim3 = int(line[1])
                    if "##$NUC1=" in file_lines[i]:
                        line = file_lines[i].split("\n")[0]
                        nuc = line.split("<")[1].split(">")[0]

            try:
                if nus == False:
                    if nuc not in self.indirect_sizes_dict.keys():
                        self.indirect_sizes_dict[nuc] = self.sizes_dim3
                    else:
                        self.indirect_sizes_dict[nuc + "_1"] = self.sizes_dim3
                else:
                    if nuc not in self.indirect_sizes_dict.keys():
                        self.indirect_sizes_dict[nuc] = self.sizes_dim3_nus
                    else:
                        self.indirect_sizes_dict[nuc + "_1"] = self.sizes_dim3_nus
            except:
                pass

    def find_sw_bruker(self) -> None:
        """
        Reading through the Bruker acqus file to find the sweep widths
        for each spectrum dimension.
        """
        for i in range(len(self.size_indirect) + 1):
            if i == 0:
                for j in range(len(self.acqus_file_lines)):
                    if "##$SW_h=" in self.acqus_file_lines[j]:
                        line = self.acqus_file_lines[j].split()
                        self.sw_direct = float(line[1])
                        break
            if i == 1:
                try:
                    self.sw_indirect = []
                    file = open("acqu2s", "r")
                    file_lines = file.readlines()
                    file.close()
                    for j in range(len(file_lines)):
                        if "##$SW_h=" in file_lines[j]:
                            line = file_lines[j].split()
                            self.sw_indirect.append(float(line[1]))
                            break
                except:
                    self.sw_indirect.append(0)
            if i == 2:
                try:
                    file = open("acqu3s", "r")
                    file_lines = file.readlines()
                    file.close()
                    for j in range(len(file_lines)):
                        if "##$SW_h=" in file_lines[j]:
                            line = file_lines[j].split()
                            self.sw_indirect.append(float(line[1]))
                            break
                except:
                    self.sw_indirect.append(0)
            if i == 3:
                try:
                    file = open("acqu4s", "r")
                    file_lines = file.readlines()
                    file.close()
                    for j in range(len(file_lines)):
                        if "##$SW_h=" in file_lines[j]:
                            line = file_lines[j].split()
                            self.sw_indirect.append(float(line[1]))
                            break
                except:
                    self.sw_indirect.append(0)

            if i > 3:
                dlg = wx.MessageDialog(
                    self.tempframe,
                    "Error: Only able to convert data with up to 4 indirect dimensions. Unable to convert data to NMRPipe format. Please check the acqus file and try again.",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                dlg.ShowModal()
                dlg.Destroy()
                exit()

    def find_nucleus_frequencies_bruker(self) -> None:
        """
        Searching through the Bruker acqus file for the Larmor frequency
        of each nucleus recorded.
        """
        self.nucleus_frequencies = []
        for i in range(len(self.size_indirect) + 1):
            if i == 0:
                for j in range(len(self.acqus_file_lines)):
                    if "##$SFO1=" in self.acqus_file_lines[j]:
                        line = self.acqus_file_lines[j].split()
                        self.nucleus_frequencies.append(float(line[1]))
                        break
            if i == 1:
                try:
                    file = open("acqu2s", "r")
                    file_lines = file.readlines()
                    file.close()
                    for j in range(len(file_lines)):
                        if "##$SFO1=" in file_lines[j]:
                            line = file_lines[j].split()
                            self.nucleus_frequencies.append(float(line[1]))
                            break
                except:
                    self.nucleus_frequencies.append(0)
            if i == 2:
                try:
                    file = open("acqu3s", "r")
                    file_lines = file.readlines()
                    file.close()
                    for j in range(len(file_lines)):
                        if "##$SFO1=" in file_lines[j]:
                            line = file_lines[j].split()
                            self.nucleus_frequencies.append(float(line[1]))
                            break
                except:
                    self.nucleus_frequencies.append(0)
            if i == 3:
                try:
                    file = open("acqu4s", "r")
                    file_lines = file.readlines()
                    file.close()
                    for j in range(len(file_lines)):
                        if "##$SFO1=" in file_lines[j]:
                            line = file_lines[j].split()
                            self.nucleus_frequencies.append(float(line[1]))
                            break
                except:
                    self.nucleus_frequencies.append(0)

    def find_labels_bruker(self) -> None:
        """
        Search through the Bruker acqus file to find the labels for each
        dimension.
        """
        self.labels = []
        for j in range(len(self.acqus_file_lines)):
            if (
                "##$NUC" in self.acqus_file_lines[j]
                and "##$NUCLEUS" not in self.acqus_file_lines[j]
                and "##$NUCLEI" not in self.acqus_file_lines[j]
            ):
                line = self.acqus_file_lines[j].split()
                self.labels.append(line[1].split("<")[1].split(">")[0])

        # If the labels are off, then change them to ID
        for i in range(len(self.labels)):
            if self.labels[i] == "off":
                self.labels[i] = "ID"

        # Work out what the correct order of labels is (i.e. if the nuclei match their frequencies etc)
        self.find_gamma_bruker()
        self.nuc1 = self.labels[0]
        # Find out the field strength using gamma and the carrier frequency
        field = self.nucleus_frequencies[0] / self.gamma[self.nuc1]
        # Find the correct order of labels
        self.labels_correct_order = []
        self.labels_correct_order.append(self.nuc1)
        for i in range(len(self.nucleus_frequencies)):
            if i == 0:
                continue
            else:
                if self.nucleus_frequencies[i] == 0:
                    self.labels_correct_order.append("ID")
                else:
                    for key in self.gamma:
                        if (
                            abs(
                                abs(field * self.gamma[key])
                                - self.nucleus_frequencies[i]
                            )
                            < 1
                        ):
                            if key in self.labels:
                                if key in self.labels_correct_order:
                                    self.labels_correct_order.append(key + "_1")
                                else:
                                    self.labels_correct_order.append(key)

    def find_gamma_bruker(self) -> None:
        """
        Producing a dictionary containing the gyromagnetic ratio of most NMR
        active nuclei.
        """
        self.gamma = {}
        self.gamma["1H"] = 267.5153151e6  # 267.5221877E6
        self.gamma["19F"] = 251.6628277e6
        self.gamma["13C"] = 67.262e6
        self.gamma["14N"] = 1.93297e7
        self.gamma["15N"] = -27.116e6
        self.gamma["31P"] = 108.282e6
        self.gamma["23Na"] = 70.882e6
        self.gamma["25Mg"] = -1.639e7
        self.gamma["39K"] = 1.2498e7
        self.gamma["41K"] = 0.686e7
        self.gamma["43Ca"] = -1.8025e7
        self.gamma["2H"] = 41.065e6
        self.gamma["7Li"] = 103.962e6
        self.gamma["17O"] = -36.264e6
        self.gamma["10B"] = 2.87471e7
        self.gamma["11B"] = 8.58406e7
        self.gamma["27Al"] = 6.97594e7
        self.gamma["29Si"] = -5.3146e7
        self.gamma["35Cl"] = 2.62401e7
        self.gamma["37Cl"] = 2.18428e7
        self.gamma["50V"] = 2.67164e7
        self.gamma["51V"] = 7.04578e7
        self.gamma["55Mn"] = 6.59777e7
        self.gamma["57Fe"] = 0.86399e7
        self.gamma["59Co"] = 6.3472e7
        self.gamma["63Cu"] = 7.0965e7
        self.gamma["65Cu"] = 7.6018e7
        self.gamma["67Zn"] = 16.767e6
        self.gamma["69Ga"] = 6.43685e7
        self.gamma["71Ga"] = 8.180163e7
        self.gamma["77Se"] = 5.115e7
        self.gamma["79Br"] = 6.70186e7
        self.gamma["81Br"] = 7.22421e7
        self.gamma["103Rh"] = -0.84579e7
        self.gamma["107Ag"] = -1.08718e7
        self.gamma["109Ag"] = -1.25001e7
        self.gamma["111Cd"] = -5.69259e7
        self.gamma["113Cd"] = -5.95504e7
        self.gamma["117Sn"] = -9.57865e7
        self.gamma["119Sn"] = -10.01926
        self.gamma["123Te"] = -7.04893e7
        self.gamma["125Te"] = -8.49722e7
        self.gamma["127I"] = 5.37937e7
        self.gamma["129Xe"] = -7.44069e7
        self.gamma["131Xe"] = 2.20564e7
        self.gamma["183W"] = 1.12070e7
        self.gamma["195Pt"] = 5.80466e7
        self.gamma["197Au"] = 0.4692e7
        self.gamma["199Hg"] = 4.81519e6
        self.gamma["201Hg"] = -1.77748e7
        self.gamma["203Tl"] = 15.43599e7
        self.gamma["205Tl"] = 15.58829e7
        self.gamma["207Pb"] = 5.64661e7

    def find_pseudo_axis_bruker(self) -> None:
        """
        Function to try to determine if a spectrum has a pseudo (non-complex)
        dimension. Read pulseprogram.precomp (or pulseprogram) file to see if the data
        has a pseudo axis. If pseudo, it should have a QF in the acqusition
        mode of that dimension. Note: this is not general, especially for
        custom written pulse sequences.
        """
        self.acqusition_modes = []
        self.pseudo_flag = 0
        if self.nmrdata.data_dimensions > 1 or self.size_indirect != []:
            self.pseudo_flag = 0
            try:
                self.pulseprogram_file = open("pulseprogram.precomp", "r")
                self.pulseprogram_file_lines = self.pulseprogram_file.readlines()
                self.pulseprogram_file.close()
            except:
                try:
                    self.pulseprogram_file = open("pulseprogram", "r")
                    self.pulseprogram_file_lines = self.pulseprogram_file.readlines()
                    self.pulseprogram_file.close()
                except:
                    dlg = wx.MessageDialog(
                        self.tempframe,
                        "Error: TD not found in acqus file. Unable to determine size of data for direct dimension. Unable to convert data to NMRPipe format. Please check the acqus file and try again.",
                        "Error",
                        wx.OK | wx.ICON_ERROR,
                    )
                    self.tempframe.Raise()
                    self.tempframe.SetFocus()
                    dlg.ShowModal()
                    dlg.Destroy()

            count = 0
            try:

                for i in range(len(self.pulseprogram_file_lines)):
                    if "AQ_mode" in self.pulseprogram_file_lines[i]:
                        line = self.pulseprogram_file_lines[i].split()
                        # Find the index of AQ_mode and delete everything before it
                        for j in range(len(line)):
                            if line[j] == "AQ_mode":
                                index = j
                        line = line[index + 1 :]
                        # Remove any terms containing brackets
                        new_line = []
                        for j in range(len(line)):
                            if "(" not in line[j]:
                                new_line.append(line[j])
                        break

                if new_line.count("QF") > 0:
                    # Then at least one of the dimensions is a pseudo axis
                    if new_line.count("QF") == 1:
                        # Then there is only one pseudo axis
                        self.pseudo_flag = 1
                    else:
                        # Then there are multiple pseudo axes
                        for i in range(len(new_line)):
                            if new_line[i] == "QF":
                                count += 1
                        self.pseudo_flag = count

                self.acqusition_modes = new_line

            except:
                try:
                    self.pseudo_flag = 0
                    # Read pulseprogram file to see if there is an FnMode option to input the acqusition mode
                    # This can be used as a guess for the acquisition modes
                    for i in range(len(self.pulseprogram_file_lines)):
                        if "FnMODE" in self.pulseprogram_file_lines[i]:
                            line = (
                                self.pulseprogram_file_lines[i]
                                .split("\n")[0]
                                .split("FnMODE:")[1]
                                .strip()
                            )
                            if len(line.split()) > 1:
                                self.acqusition_modes = line.split()
                            elif line.lower() == "echo-antiecho":
                                self.acqusition_modes = ["Echo-Antiecho"]
                            elif line.lower() == "complex":
                                self.acqusition_modes = ["Complex"]
                            elif line.lower() == "tppi":
                                self.acqusition_modes = ["TPPI"]
                            elif line.lower() == "qf":
                                self.acqusition_modes = ["QF"]

                except:
                    try:
                        self.pseudo_flag = 0
                        self.acqusition_modes = []
                        for line in self.pulseprogram_file_lines:
                            if (
                                "define list<frequency>" in line
                                or "define list<delay>" in line
                            ):
                                self.pseudo_flag = 1
                                self.acqusition_modes = ["QF"]

                                break

                    except:
                        self.pseudo_flag = 0
                        self.acqusition_modes = ["Complex"]

            while len(self.acqusition_modes) < len(self.size_indirect):
                for i in range(len(self.size_indirect) - len(self.acqusition_modes)):
                    self.acqusition_modes.append("Complex")

            if len(self.size_indirect) == 1:
                if self.pseudo_flag == 0:
                    if self.size_indirect[0] == 1:
                        self.pseudo_flag = 1
                        self.acqusition_modes[-1] = "QF"
            elif len(self.size_indirect) == 2:
                if self.pseudo_flag == 0:
                    if self.size_indirect[0] == 1:
                        self.pseudo_flag += 1
                        self.acqusition_modes[-2] = "QF"
                    if self.size_indirect[1] == 1:
                        self.pseudo_flag += 1
                        self.acqusition_modes[-1] = "QF"

            # If there is a pseudo axis, then the number of real points in that dimension is the same as the number of complex points
            if self.pseudo_flag == 1:
                sizes = []
                keys = []
                for key in self.udic:
                    if key == "ndim":
                        continue
                    else:
                        keys.append(key)
                        sizes.append(self.udic[key]["size"])

                # Find the index of the minimum size (this is default set to the pseudo axis)
                index = sizes.index(min(sizes))
                self.udic[keys[index]]["complex"] = False
                self.udic[keys[index]]["encoding"] = "Real"
                self.udic[keys[index]]["obs"] = 0
                self.udic[keys[index]]["sw"] = 0
                self.udic[keys[index]]["car"] = 0

            elif self.pseudo_flag == 2:
                sizes = []
                keys = []
                for key in self.udic:
                    if key == "ndim":
                        continue
                    else:
                        keys.append(key)
                        sizes.append(self.udic[key]["size"])

                # Find the index of the minimum size (this is default set to the pseudo axis)
                index = sizes.index(min(sizes))
                self.udic[keys[index]]["encoding"] = "Real"
                self.udic[keys[index]]["complex"] = False
                self.udic[keys[index]]["obs"] = 1
                self.udic[keys[index]]["sw"] = 1
                self.udic[keys[index]]["car"] = 1
                # Then find the next smallest size and set this to the pseudo axis
                sizes[index] = max(sizes)
                index = sizes.index(min(sizes))
                self.udic[keys[index]]["encoding"] = "Real"
                self.udic[keys[index]]["complex"] = False
                self.udic[keys[index]]["obs"] = 1
                self.udic[keys[index]]["sw"] = 1
                self.udic[keys[index]]["car"] = 1

    def find_temperature_bruker(self) -> None:
        """
        Find the temperature the spectrum was recorded at
        """
        self.temperature = 298.15  # Default temperature is 298.15K
        for i in range(len(self.acqus_file_lines)):
            if "##$TE=" in self.acqus_file_lines[i]:
                line = self.acqus_file_lines[i].split()
                self.temperature = float(line[1])
                break

    def calculate_carrier_frequency_bruker(self) -> None:
        """
        Calculate the carrier frequency for each dimension. If the direct
        dimension is proton, calculate the carrier based on water and include
        water referencing as an option.
        """

        if len(self.size_indirect) == 0:
            # return water chemical shift in range 0-100oC
            if (
                self.labels_correct_order[0] == "1H"
                or self.labels_correct_order[0] == "H1"
                or self.labels_correct_order[0] == "H"
            ):
                self.water_ppm = 7.83 - self.temperature / 96.9

                # Use O1/BF1 to calculate a second carrier frequency in case not centred on water
                for j in range(len(self.acqus_file_lines)):
                    if "##$O1=" in self.acqus_file_lines[j]:
                        line = self.acqus_file_lines[j].split()
                        self.O1 = float(line[1])
                        break
                for j in range(len(self.acqus_file_lines)):
                    if "##$BF1=" in self.acqus_file_lines[j]:
                        line = self.acqus_file_lines[j].split()
                        self.BF1 = float(line[1])
                        break
                self.carrier_frequency_1 = self.O1 / self.BF1

                self.references_proton = [self.water_ppm, self.carrier_frequency_1]
                self.references_proton_labels = ["H2O", "O1/BF1"]
            else:
                # Use O1/BF1 to calculate a second carrier frequency in case not centred on water
                for j in range(len(self.acqus_file_lines)):
                    if "##$O1=" in self.acqus_file_lines[j]:
                        line = self.acqus_file_lines[j].split()
                        self.O1 = float(line[1])
                        break
                for j in range(len(self.acqus_file_lines)):
                    if "##$BF1=" in self.acqus_file_lines[j]:
                        line = self.acqus_file_lines[j].split()
                        self.BF1 = float(line[1])
                        break
                self.carrier_frequency_1 = self.O1 / self.BF1

                self.references_proton = [self.carrier_frequency_1]
                self.references_proton_labels = ["O1/BF1"]

        else:
            if (
                self.labels_correct_order[0] == "1H"
                or self.labels_correct_order[0] == "H1"
                or self.labels_correct_order[0] == "H"
            ):
                # return water chemical shift in range 0-100oC
                self.water_ppm = 7.83 - self.temperature / 96.9

            # Use O1/BF1 to calculate a second carrier frequency in case not centred on water
            for j in range(len(self.acqus_file_lines)):
                if "##$O1=" in self.acqus_file_lines[j]:
                    line = self.acqus_file_lines[j].split()
                    self.O1 = float(line[1])
                    break
            for j in range(len(self.acqus_file_lines)):
                if "##$BF1=" in self.acqus_file_lines[j]:
                    line = self.acqus_file_lines[j].split()
                    self.BF1 = float(line[1])
                    break
            self.carrier_frequency_1 = self.O1 / self.BF1

            # Calculate carrier frequencies based on O2/BF2, O3/BF3

            for j in range(len(self.acqus_file_lines)):
                if "##$O2=" in self.acqus_file_lines[j]:
                    line = self.acqus_file_lines[j].split()
                    self.O2 = float(line[1])
                    break
            for j in range(len(self.acqus_file_lines)):
                if "##$O3=" in self.acqus_file_lines[j]:
                    line = self.acqus_file_lines[j].split()
                    self.O3 = float(line[1])
                    break
            for j in range(len(self.acqus_file_lines)):
                if "##$BF2=" in self.acqus_file_lines[j]:
                    line = self.acqus_file_lines[j].split()
                    self.BF2 = float(line[1])
                    break
            for j in range(len(self.acqus_file_lines)):
                if "##$BF3=" in self.acqus_file_lines[j]:
                    line = self.acqus_file_lines[j].split()
                    self.BF3 = float(line[1])
                    break
            # Calculate the carrier frequency
            self.carrier_frequency_2 = self.O2 / self.BF2
            self.carrier_frequency_3 = self.O3 / self.BF3

            self.ppms_referenced = []
            self.ppms_referenced_labels = []

            if (
                self.labels_correct_order[0] == "1H"
                or self.labels_correct_order[0] == "H1"
                or self.labels_correct_order[0] == "H"
            ):
                # Calculate referenced carrier frequencies based on water chemical shifts
                self.sfrq0 = self.nucleus_frequencies[0] / (1 + self.water_ppm * 1e-6)
                # Frequency for nucleus A = sfrq0 * gammaH/gammaA
                self.dfrq_13C = self.sfrq0 * 0.251449530
                self.dfrq_15N = self.sfrq0 * 0.101329118
                self.dfrq_P31 = self.sfrq0 * 0.4048064954
                self.dfrq_F19 = self.sfrq0 * 0.9412866605363297

                for i, label in enumerate(self.labels_correct_order):
                    # if(i+1>len(self.nucleus_frequencies)):
                    #     break
                    if label == "15N" or label == "N15":
                        self.ppms_referenced.append(
                            (
                                (self.nucleus_frequencies[i] - self.dfrq_15N)
                                / self.dfrq_15N
                            )
                            * 1e6
                        )
                        self.ppms_referenced_labels.append("15N (Referenced to H2O)")
                    if label == "13C" or label == "C13":
                        self.ppms_referenced.append(
                            (
                                (self.nucleus_frequencies[i] - self.dfrq_13C)
                                / self.dfrq_13C
                            )
                            * 1e6
                        )
                        self.ppms_referenced_labels.append("13C (Referenced to H2O)")
                    if label == "31P" or label == "P31":
                        self.ppms_referenced.append(
                            (
                                (self.nucleus_frequencies[i] - self.dfrq_P31)
                                / self.dfrq_P31
                            )
                            * 1e6
                        )
                        self.ppms_referenced_labels.append("31P (Referenced to H2O)")
                    if label == "19F" or label == "F19":
                        self.ppms_referenced.append(
                            (
                                (self.nucleus_frequencies[i] - self.dfrq_F19)
                                / self.dfrq_F19
                            )
                            * 1e6
                        )
                        self.ppms_referenced_labels.append("19F (Referenced to H2O)")

            if (
                self.labels_correct_order[0] == "1H"
                or self.labels_correct_order[0] == "H1"
                or self.labels_correct_order[0] == "H"
            ):
                self.references_proton = [
                    self.water_ppm,
                    self.carrier_frequency_1,
                    self.carrier_frequency_2,
                    self.carrier_frequency_3,
                ]
                self.references_proton_labels = ["H2O", "O1/BF1", "O2/BF2", "O3/BF3"]
            else:
                self.references_proton = [
                    self.carrier_frequency_1,
                    self.carrier_frequency_2,
                    self.carrier_frequency_3,
                ]
                self.references_proton_labels = ["O1/BF1", "O2/BF2", "O3/BF3"]

            self.references_other = []
            self.references_other_labels = []
            for i, ppm in enumerate(self.ppms_referenced):
                self.references_other.append(ppm)
                self.references_other_labels.append(self.ppms_referenced_labels[i])
            self.references_other.append(self.carrier_frequency_1)
            self.references_other.append(self.carrier_frequency_2)
            self.references_other.append(self.carrier_frequency_3)
            self.references_other_labels.append("O1/BF1")
            self.references_other_labels.append("O2/BF2")
            self.references_other_labels.append("O3/BF3")

    def find_bruker_digital_filter_parameters(self):
        # Search through the acqus file to find the parameters
        acqus_file_lines = open(self.nmrdata.parameter_file, "r").readlines()

        # Find the decim, dspfvs and grpdly parameters
        self.decim = 0
        self.dspfvs = 0
        self.grpdly = 0
        try:
            for i in range(len(acqus_file_lines)):
                if "##$DECIM=" in acqus_file_lines[i]:
                    line = acqus_file_lines[i].split()
                    self.decim = float(line[1])
                if "##$DSPFVS=" in acqus_file_lines[i]:
                    line = acqus_file_lines[i].split()
                    self.dspfvs = int(line[1])
                if "##$GRPDLY=" in acqus_file_lines[i]:
                    line = acqus_file_lines[i].split()
                    self.grpdly = float(line[1])
            self.include_digital_filter = True

        except:
            self.decim = 0
            self.dspfvs = 0
            self.grpdly = 0
            self.include_digital_filter = False

    def find_bruker_scaling_parameters(self):
        # Search through the acqus file to find the NS and NC parameters
        acqus_file_lines = open(self.nmrdata.parameter_file, "r").readlines()

        # Find the NS and NC parameters
        try:
            for i in range(len(acqus_file_lines)):
                if "##$NS=" in acqus_file_lines[i]:
                    line = acqus_file_lines[i].split()
                    self.NS = int(line[1])
                if "##$NC=" in acqus_file_lines[i]:
                    line = acqus_file_lines[i].split()
                    self.NC = int(line[1])
            self.include_scaling = True
        except:
            self.NS = 0
            self.NC = 0
            self.include_scaling = False

    def determine_byte_order(self):
        try:
            with open(self.nmrdata.parameter_file) as file:
                lines = file.readlines()
                for line in lines:
                    if ("##$BYTORDA") in line:
                        self.byte_order = line.split("\n")[0].split()[-1]
                        break
        except:
            self.byte_order = 0

    def determine_byte_size(self):
        try:
            with open(self.nmrdata.parameter_file) as file:
                lines = file.readlines()
                for line in lines:
                    if ("##$DTYPA") in line:
                        self.d_type = line.split("\n")[0].split()[-1]
                        break
        except:
            self.d_type = 0
