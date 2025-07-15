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
import numpy as np
import nmrglue as ng
import warnings

warnings.simplefilter("ignore", UserWarning)


class ParameterExtractorVarian:
    def __init__(self, converter) -> None:
        """
        This class contains functions to search through the varian
        procpar file to find relevant parameters needed for conversion
        to nmrPipe format.
        """
        self.converter = converter

        # Reading the Varian parameter file
        self.procpar_file = open("procpar", "r")
        self.procpar_file_lines = self.procpar_file.readlines()
        self.procpar_file.close()

        # Getting nmrglue to guess the intial parameter files
        self.udic = ng.varian.guess_udic(
            self.converter.nmr_dic, self.converter.nmr_data
        )
        self.find_size_varian()
        self.find_axes_pseudo_varian()
        self.find_sw_varian()
        self.find_nucleus_frequencies_varian()
        self.find_labels_varian()
        self.find_temperature_varian()
        self.calculate_carrier_frequency_varian()

    def find_size_varian(self):
        """
        Finding the np, ni, ni2 values to get spectrum dimensions
        """
        # Look for np in the procpar file
        include_np_value = False
        include_ni_value = False
        include_ni2_value = False
        include_ni3_value = False
        self.size_direct = 0
        self.size_indirect = []
        for i in range(len(self.procpar_file_lines)):
            line = self.procpar_file_lines[i].split()
            if include_np_value == True:
                self.size_direct = int(line[1])
                include_np_value = False
            elif line[0] == "np":
                include_np_value = True
            elif include_ni_value == True:
                self.size_indirect.append(int(line[1]))
                include_ni_value = False
            elif line[0] == "ni":
                include_ni_value = True
            elif include_ni2_value == True:
                self.size_indirect.append(int(line[1]))
                include_ni2_value = False
            elif line[0] == "ni2":
                include_ni2_value = True
            elif include_ni3_value == True:
                self.size_indirect.append(int(line[1]))
                include_ni3_value = False
            elif line[0] == "ni3":
                include_ni3_value = True

        if self.size_direct == 0:
            dlg = wx.MessageDialog(
                self.converter.tempframe,
                "Error: np not found in procpar file. Unable to determine size"
                " of data for direct dimension. Unable to convert data to NMRPipe"
                " format. Please check the procpar file and try again.",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            self.converter.tempframe.Raise()
            self.converter.tempframe.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            exit()

        return self.size_direct, self.size_indirect

    def find_sw_varian(self) -> None:
        """
        Looking for sw, sw1, sw2 in the procpar file
        """
        include_sw_value = False
        include_sw1_value = False
        include_sw2_value = False
        self.sw_direct = 0
        self.sw_indirect = {}
        for i in range(len(self.procpar_file_lines)):
            line = self.procpar_file_lines[i].split()
            if include_sw_value == True:
                self.sw_direct = float(line[1])
                include_sw_value = False
            elif line[0] == "sw":
                include_sw_value = True
            elif include_sw1_value == True:
                self.sw_indirect["sw1"] = float(line[1])
                include_sw1_value = False
            elif line[0] == "sw1":
                include_sw1_value = True
            elif include_sw2_value == True:
                self.sw_indirect["sw2"] = float(line[1])
                include_sw2_value = False
            elif line[0] == "sw2":
                include_sw2_value = True

        return self.sw_direct, self.sw_indirect

    def find_nucleus_frequencies_varian(self) -> None:
        """
        Looking for sfrq, dfrq2, dfrq3 in the procpar file
        """
        include_sfrq_value = False
        include_dfrq_value = False
        include_dfrq2_value = False
        include_dfrq3_value = False
        self.nucleus_frequency_direct = 0
        self.nucleus_frequencies_indirect = []
        self.nucleus_frequencies_indirect_order = []
        include_zero = False
        for i in range(len(self.procpar_file_lines)):
            line = self.procpar_file_lines[i].split()
            if include_sfrq_value == True:
                self.nucleus_frequency_direct = float(line[1])
                include_sfrq_value = False
            elif line[0] == "sfrq":
                include_sfrq_value = True
            elif include_dfrq_value == True:
                if line[1] != "0":
                    self.nucleus_frequencies_indirect.append(float(line[1]))
                    self.nucleus_frequencies_indirect_order.append("dfrq")
                else:
                    include_zero = True
                include_dfrq_value = False
            elif line[0] == "dfrq":
                include_dfrq_value = True
            elif include_dfrq2_value == True:
                if line[1] != "0":
                    self.nucleus_frequencies_indirect.append(float(line[1]))
                    self.nucleus_frequencies_indirect_order.append("dfrq2")
                else:
                    include_zero = True
                include_dfrq2_value = False
            elif line[0] == "dfrq2":
                include_dfrq2_value = True
            elif include_dfrq3_value == True:
                if line[1] != "0":
                    self.nucleus_frequencies_indirect.append(float(line[1]))
                    self.nucleus_frequencies_indirect_order.append("dfrq3")
                else:
                    include_zero = True
                include_dfrq3_value = False
            elif line[0] == "dfrq3":
                include_dfrq3_value = True

        if self.nucleus_frequency_direct == 0:
            dlg = wx.MessageDialog(
                self.converter.tempframe,
                "Error: sfrq not found in procpar file. Unable to determine"
                " nucleus frequency for direct dimension. Unable to convert "
                "data to NMRPipe format. Please check the procpar file and "
                "try again.",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            self.converter.tempframe.Raise()
            self.converter.tempframe.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            exit()

        if self.nucleus_frequencies_indirect == []:
            dlg = wx.MessageDialog(
                self.converter.tempframe,
                "Error: dfrq/dfrq2/dfrq3 not found in procpar file. Would you"
                " like to continue anyway?",
                "Error",
                wx.YES_NO | wx.ICON_ERROR,
            )
            self.converter.tempframe.Raise()
            self.converter.tempframe.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            if dlg.GetReturnCode() == wx.ID_NO:
                exit()
            else:
                self.nucleus_frequencies_indirect = [0.0, 0.0, 0.0]
                self.nucleus_frequencies_indirect_order = ["dfrq", "dfrq2", "dfrq3"]

        self.nucleus_frequencies_indirect.reverse()
        if include_zero == True:
            self.nucleus_frequencies_indirect.append(0.0)
        self.nucleus_frequencies = (
            [self.nucleus_frequency_direct] + self.nucleus_frequencies_indirect + [0.0]
        )

        return (
            self.nucleus_frequency_direct,
            self.nucleus_frequencies_indirect,
            self.nucleus_frequencies_indirect_order,
            self.nucleus_frequencies,
        )

    def find_labels_varian(self) -> None:
        """
        # Finding nucleus labels for each dimensions.
        # tn, dn, dn2, dn3 in the procpar file
        """

        include_tn_value = False
        include_dn_value = False
        include_dn2_value = False
        include_dn3_value = False
        self.label_direct = ""
        self.labels_indirect = []
        for i in range(len(self.procpar_file_lines)):
            line = self.procpar_file_lines[i].split()
            if include_tn_value == True:
                self.label_direct = line[1].split('"')[1]
                include_tn_value = False
            elif line[0] == "tn":
                include_tn_value = True
            elif include_dn_value == True:
                self.labels_indirect.append(line[1].split('"')[1])
                include_dn_value = False
            elif line[0] == "dn":
                include_dn_value = True
            elif include_dn2_value == True:
                self.labels_indirect.append(line[1].split('"')[1])
                include_dn2_value = False
            elif line[0] == "dn2":
                include_dn2_value = True
            elif include_dn3_value == True:
                self.labels_indirect.append(line[1].split('"')[1])
                include_dn3_value = False
            elif line[0] == "dn3":
                include_dn3_value = True

        if self.label_direct == "":
            dlg = wx.MessageDialog(
                self.converter.tempframe,
                "Error: Label for direct dimension (tn) not found in procpar file. Setting label as 1.",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            self.converter.tempframe.Raise()
            self.converter.tempframe.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            self.label_direct = "1"

        if self.labels_indirect == []:
            dlg = wx.MessageDialog(
                self.tempframe,
                "Error: Labels for indirect dimensions (dn/dn2/dn3) not found in procpar file. Setting labels as 2, 3, 4.",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            self.converter.tempframe.Raise()
            self.converter.tempframe.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            self.labels_indirect = ["2", "3", "4"]

        if self.converter.other_params == True:
            self.labels_correct_order = (
                [self.labels_indirect[0]]
                + [self.converter.arrayed_parameter]
                + self.labels_indirect[1:]
            )
        else:
            self.labels_correct_order = self.labels_indirect
            if self.reverse_acquisition_order == True:
                self.labels_correct_order.reverse()

    def find_temperature_varian(self) -> None:
        """
        Finding the temperature the spectrum was recorded at
        """
        self.temperature = 298.15  # Default temperature is 298.15K
        include_temp_value = False
        for i in range(len(self.procpar_file_lines)):
            line = self.procpar_file_lines[i].split()
            if include_temp_value == True:
                self.temperature = float(line[1]) + 273.15
                include_temp_value = False
            elif line[0] == "temp":
                include_temp_value = True

    def calculate_carrier_frequency_varian(self) -> None:
        """
        Calculate the carrier frequency for each dimension
        If the direct dimension is proton calculate based on
        water and include water referencing as an option.
        """

        if (
            self.label_direct == "1H"
            or self.label_direct == "H1"
            or self.label_direct == "H"
        ):
            # return water chemical shift in range 0-100oC
            self.water_ppm = 7.83 - self.temperature / 96.9

            self.references_proton = [self.water_ppm, 0]
            self.references_proton_labels = ["H2O", "Manual"]

            self.references_other = []
            self.references_other_labels = []

            # Determine frequency of other nuclei at 0ppm based on water ppm
            self.sfrq0 = self.nucleus_frequency_direct / (1 + self.water_ppm * 1e-6)
            self.dfrq_13C = self.sfrq0 * 0.251449530
            self.dfrq_15N = self.sfrq0 * 0.101329118
            self.dfrq_P31 = self.sfrq0 * 0.4048064954
            self.dfrq_F19 = self.sfrq0 * 0.941

            # Use dfrq2 and dfrq3 to calculate chemical shifts of the indirect dimensions

            # Determine order of nucleus labels based on dfrq2/dfrq3
            self.ordered_indirect_labels = []
            self.sw_ordered = []
            found_match = False
            for i, value in enumerate(self.nucleus_frequencies_indirect):
                if np.abs(value - self.dfrq_13C) / self.dfrq_13C * 100 < 1:
                    self.ordered_indirect_labels.append("C13")
                    self.references_other.append(
                        (value - self.dfrq_13C) / self.dfrq_13C * 1e6
                    )
                    self.references_other_labels.append("C13 (Referenced to H2O)")
                    found_match = True
                elif np.abs(value - self.dfrq_15N) / self.dfrq_15N * 100 < 1:
                    self.ordered_indirect_labels.append("N15")
                    self.references_other.append(
                        (value - self.dfrq_15N) / self.dfrq_15N * 1e6
                    )
                    self.references_other_labels.append("N15 (Referenced to H2O)")
                    found_match = True
                elif np.abs(value - self.dfrq_P31) / self.dfrq_P31 * 100 < 1:
                    self.ordered_indirect_labels.append("P31")
                    self.references_other.append(
                        (value - self.dfrq_P31) / self.dfrq_P31 * 1e6
                    )
                    self.references_other_labels.append("P31 (Referenced to H2O)")
                    found_match = True
                elif np.abs(value - self.dfrq_F19) / self.dfrq_F19 * 100 < 1:
                    self.ordered_indirect_labels.append("F19")
                    self.references_other.append(
                        (value - self.dfrq_F19) / self.dfrq_F19 * 1e6
                    )
                    self.references_other_labels.append("F19 (Referenced to H2O)")
                    found_match = True
                elif np.abs(value - self.sfrq0) / self.sfrq0 * 100 < 1:
                    self.ordered_indirect_labels.append("H1")
                    self.references_other.append(self.water_ppm)
                    self.references_other_labels.append("H2O")
                    found_match = True
                if found_match == True:
                    # self.sw_ordered.append(self.sw_indirect[i])
                    found_match = False

        else:
            self.references_proton = [0.0]
            self.references_proton_labels = ["Manual"]

            self.references_other = [0.0]
            self.references_other_labels = ["Manual"]

    def find_axes_pseudo_varian(self) -> None:
        """
        Search through procpar to find array values
        """
        self.phase = False
        self.phase2 = False
        self.other_params = False

        self.found_array = False
        for i, line in enumerate(self.procpar_file_lines):
            if self.found_array == True:
                array = line.split()[1]
                self.found_array = False
            if line.split()[0] == "array":
                self.found_array = True

        array = array.split('"')[1].split(",")
        if array == [""]:
            array = []

        self.reorder_array = False
        count = 0
        if "phase" in array:
            self.phase = True
            count += 1
        if "phase2" in array:
            self.phase2 = True
            count += 1
        if count != len(array):
            self.other_params = True

        # Delete phase and phase2 from array
        self.found_arrayed_parameter = False
        self.number_of_arrayed_parameters = 0
        self.phases = []
        if self.other_params == True:
            if self.phase == True:
                array.remove("phase")
                self.phases.append("phase")
            if self.phase2 == True:
                array.remove("phase2")
                self.phases.append("phase2")
            self.arrayed_parameter = array[0]

            # Search through procpar file to find the number of values of the arrayed parameter

            for i, line in enumerate(self.procpar_file_lines):
                if self.found_arrayed_parameter == True:
                    self.number_of_arrayed_parameters = int(line.split()[0])
                    self.found_arrayed_parameter = False
                if line.split()[0] == self.arrayed_parameter:
                    self.found_arrayed_parameter = True

        for i, line in enumerate(self.procpar_file_lines):
            if self.found_array == True:
                array = line.split()[1]
                self.found_array = False
            if line.split()[0] == "array":
                self.found_array = True

        array = array.split('"')[1].split(",")
        if array == [""]:
            array = []

        self.reorder_array = False
        count = 0
        if "phase" in array and "phase2" in array:
            self.phase = True
            count += 1
            if array[-1] == "phase" or array[-1] == "phase2" and len(array) > 2:
                self.reorder_array = True
        elif "phase2" in array:
            self.phase2 = True
            count += 1
            if array[-1] == "phase2" and len(array) > 1:
                self.reorder_array = True
        elif "phase" in array:
            self.phase = True
            count += 1
            if array[-1] == "phase" and len(array) > 1:
                self.reorder_array = True
        self.reverse_acquisition_order = False
        if "phase" in array and "phase2" in array:
            if array == ["phase", "phase2"]:
                self.reverse_acquisition_order = True
            else:
                self.reverse_acquisition_order = False

        if self.reverse_acquisition_order == False:
            self.size_indirect.reverse()

        if 0 in self.size_indirect:
            self.size_indirect.remove(0)
            self.size_indirect.append(0)

    def find_varian_scaling_parameters(self):
        found_nt = False
        self.include_scaling = False
        for i, line in enumerate(self.nmrdata.procpar_file_lines):
            if found_nt == True:
                self.NS = int(line.split()[1])
                self.include_scaling = True
                break
            elif line.split()[0] == "nt":
                self.NS = int(line.split()[1])
                found_nt = True
