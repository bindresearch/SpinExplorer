#!/usr/bin/env python3

"""MIT License

Copyright (c) 2025 James Eaton, Andrew Baldwin (University of Oxford)

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

import json
import pathlib
import wx
from typing import Union, Dict, Any


class Read_json:

    def __init__(self, notebook, nmrdata, dimension_tabs):
        """
        This class will read a parameter.json file and ensure that all of the
        parameters are valid (if not will go back to default parameters).
        It will then load all the correct parameters etc into the
        SpinProcess app.
        """

        self.notebook = notebook

        # Checking to see if parameters.json is found and if user wants to
        # read in these values
        if self.find_json() == True:
            parameter_dictionary = self.read()
            if parameter_dictionary == False:
                return
            else:
                parameter_inputs = InputParameters(
                    parameter_dictionary,
                    nmrdata,
                    notebook,
                    dimension_tabs,
                )

    def find_json(self) -> Union[bool, None]:
        """
        Code to determine if there is a parameters.json file in the current
        directory
        """
        if pathlib.Path("parameters.json").exists() == True:
            return self.ask_user()
        else:
            return None

    def ask_user(self) -> bool:
        """
        Code to ask the user if they want to load in parameters from the
        found parameters.json file.
        """
        # Asking the user if they want to overwrite the previous saved state
        dlg = wx.MessageDialog(
            self.notebook,
            "A set of saved parameters has been found (parameters.json). Would you like to read these values?",
            "Warning",
            wx.YES_NO | wx.ICON_QUESTION,
        )
        result = dlg.ShowModal()
        if result == wx.ID_NO:
            dlg.Destroy()
            return False
        else:
            dlg.Destroy()
            return True

    def read(self):
        """
        Reading the parameters.json file and giving an error to the user if this
        could not be read correctly.
        """
        try:
            with open("parameters.json", "r") as file:
                parameter_dictionary = json.load(file)["processing"]
            return parameter_dictionary
        except:
            # Unable to read the converter.json file effectively
            # Giving a warning to the user and returning
            message = "Unable to read the saved processing parameters in parameters.json. Continuing using default parameters."
            dlg = wx.MessageDialog(
                self.notebook, message, "Warning", wx.OK | wx.ICON_ERROR
            )
            dlg.ShowModal()
            dlg.Destroy()

            return False


class InputParameters:
    def __init__(self, parameter_dictionary, nmrdata, notebook, dimension_tabs):
        """
        This class will insert the saved processing parameters into the GUI
        and also update their stored variable values.

        1 - Check number of dimensions in parameters.json matches the
            number of dimensions in the graphical interface.

        2 - For each dimension, try to load the saved parameters, if that
            dimension fails, give an error to say the parameters could not
            be loaded for dimension XXX and that default values are being
            used.

        Each dimension loading is in a try/except look
        and if the loading doesn't work (if for example the
        json file has been tampered with) then this will
        tell the user which parts have not been loaded correctly
        and will instead use default values.
        """

        self.notebook = notebook
        self.nmrdata = nmrdata

        check = self.check_dimensions(parameter_dictionary, dimension_tabs)
        if check == False:
            return
        for i, dimension_tab in enumerate(dimension_tabs):
            # Loading in saved dimension values
            # try:
            label = (
                "Dimension {}".format(i) + " (" + str(self.nmrdata.axislabels[i]) + ")"
            )
            dictionary = parameter_dictionary[label]
            self.load_dimension(i, dimension_tab, dictionary)
            # except:
            #     self.error_dimension_loading(i)

    def check_dimensions(self, parameter_dictionary: Dict, dimension_tabs) -> bool:
        """
        Checking that the dimensionality of the processing parameters
        is the same as that of the number of graphical interface
        tabs.
        """
        if len(parameter_dictionary.keys()) != len(dimension_tabs):
            message = "The number of processing dimensions in parameters.json ({}) does not match the expected number of dimensions ({})".format(
                len(parameter_dictionary.keys()), len(dimension_tabs)
            )
            dlg = wx.MessageDialog(
                self.notebook,
                message,
                "Warning",
                wx.OK,
            )
            result = dlg.ShowModal()
            return False
        else:
            return True

    def load_dimension(self, dimension: int, dimension_tab, dictionary):
        """
        Loading in the parameters the the dimension.
        """
        if dimension == 0:
            # Load the solvent suppression values
            self.load_solvent_suppression(dimension, dimension_tab, dictionary)
        self.load_linear_prediction(dimension, dimension_tab, dictionary)
        self.load_apodization(dimension, dimension_tab, dictionary)
        self.load_zero_filling(dimension, dimension_tab, dictionary)
        self.load_fourier_transform(dimension, dimension_tab, dictionary)
        self.load_phasing(dimension, dimension_tab, dictionary)
        self.load_extraction(dimension, dimension_tab, dictionary)
        self.load_baseline_correction(dimension, dimension_tab, dictionary)

        self.notebook.Refresh()

    def load_solvent_suppression(self, dimension, dimension_tab, dictionary):
        """
        Reading in the saved parameters associated with solvent
        suppression.
        """
        suppression_flag = bool(
            dictionary["Solvent Suppression"]["Solvent Suppression Flag"]
        )
        dimension_tab.solvent_suppression.solvent_suppression_checkbox.SetValue(
            suppression_flag
        )
        dimension_tab.solvent_suppression.solvent_suppression_checkbox_value = (
            suppression_flag
        )

        filter_choice = int(dictionary["Solvent Suppression"]["Filter Selection"][0])
        dimension_tab.solvent_suppression.solvent_suppression_filter_selection = (
            filter_choice
        )

        lowpass_choice = int(dictionary["Solvent Suppression"]["Lowpass Shape"][0])
        dimension_tab.solvent_suppression.solvent_suppression_lowpass_shape_selection = (
            lowpass_choice
        )

    def load_linear_prediction(self, dimension, dimension_tab, dictionary):
        """
        Reading in the saved parameters associated with linear prediction.
        """
        if dimension == 0:
            linear_prediction_flag = bool(
                dictionary["Linear Prediction"]["Linear Prediction Flag"]
            )
            dimension_tab.linear_prediction.linear_prediction_checkbox.SetValue(
                linear_prediction_flag
            )
            dimension_tab.linear_prediction.linear_prediction_checkbox_value = (
                linear_prediction_flag
            )

            predicted_points_selection = int(
                dictionary["Linear Prediction"]["Add predicted points"][0]
            )
            dimension_tab.linear_prediction.linear_prediction_combobox.SetSelection(
                predicted_points_selection
            )

            dimension_tab.linear_prediction.linear_prediction_options_selection = (
                predicted_points_selection
            )

            predicted_coefficients_selection = int(
                dictionary["Linear Prediction"]["Predicted coefficients"][0]
            )
            dimension_tab.linear_prediction.linear_prediction_coefficients_combobox.SetSelection(
                predicted_coefficients_selection
            )
            dimension_tab.linear_prediction.linear_prediction_coefficients_selection = (
                predicted_coefficients_selection
            )

        else:
            key = list(dictionary["Linear Prediction"].keys())[0]
            if key == "None":
                dimension_tab.linear_prediction.linear_prediction_radio_box_indirect.SetSelection(
                    0
                )
                dimension_tab.linear_prediction.linear_prediction_radio_box_indirect_selection = (
                    0
                )

            elif key == "Linear Prediction":
                dimension_tab.linear_prediction.linear_prediction_radio_box_indirect.SetSelection(
                    1
                )
                dimension_tab.linear_prediction.linear_prediction_radio_box_indirect_selection = (
                    1
                )

                dimension_tab.linear_prediction.on_linear_prediction_radio_box_indirect(
                    wx.EVT_RADIOBOX
                )

                predicted_points_selection = int(
                    dictionary["Linear Prediction"][key]["Add predicted points"][0]
                )
                dimension_tab.linear_prediction.linear_prediction_radio_box_indirect.SetSelection(
                    predicted_points_selection
                )

                dimension_tab.linear_prediction.linear_prediction_indirect_options_selection = (
                    predicted_points_selection
                )

                predicted_coefficients_selection = int(
                    dictionary["Linear Prediction"][key]["Predicted coefficients"][0]
                )
                dimension_tab.linear_prediction.linear_prediction_coefficients_combobox_indirect.SetSelection(
                    predicted_coefficients_selection
                )

                dimension_tab.linear_prediction.linear_prediction_indirect_coefficients_selection = (
                    predicted_coefficients_selection
                )

                dimension_tab.linear_prediction.on_linear_prediction_combobox_coefficients_indirect(
                    wx.EVT_COMBOBOX
                )

            elif key == "NUS reconstruction":
                dimension_tab.linear_prediction.linear_prediction_radio_box_indirect.SetSelection(
                    2
                )
                dimension_tab.linear_prediction.linear_prediction_radio_box_indirect_selection = (
                    2
                )

                nusfile = dictionary["Linear Prediction"][key]["NUS file"]
                nus_extension = int(
                    dictionary["Linear Prediction"][key]["NUS extension"]
                )
                nus_cpu = int(dictionary["Linear Prediction"][key]["NUS CPU number"])
                nus_iterations = int(
                    dictionary["Linear Prediction"][key]["NUS iterations"]
                )

                dimension_tab.linear_prediction.nuslist_name_indirect = nusfile
                dimension_tab.linear_prediction.smile_nus_file_textcontrol_indirect.SetValue(
                    nusfile
                )
                dimension_tab.linear_prediction.smile_data_extension_number_indirect = (
                    nus_extension
                )
                dimension_tab.linear_prediction.smile_nus_extension_textcontrol_indirect.SetValue(
                    nus_extension
                )
                dimension_tab.linear_prediction.smile_data_extension_number_indirect = (
                    nus_cpu
                )
                dimension_tab.linear_prediction.smile_nus_cpu_textcontrol_indirect.SetValue(
                    nus_cpu
                )
                dimension_tab.linear_prediction.nus_iterations_indirect = nus_iterations
                dimension_tab.linear_prediction.smile_nus_iterations_textcontrol_indirect = (
                    nus_iterations
                )
                dimension_tab.linear_prediction.on_linear_prediction_radio_box_indirect(
                    wx.EVT_RADIOBOX
                )

    def load_apodization(self, dimension, dimension_tab, dictionary):
        """
        Reading in the saved parameters associated with apodization
        """

        apodization_flag = bool(dictionary["Apodization"]["Apodization flag"])
        apodization_value = int(dictionary["Apodization"]["Selection"])

        dimension_tab.apodization.apodization_checkbox.SetValue(apodization_flag)
        dimension_tab.apodization.apodization_checkbox_value = apodization_flag

        if apodization_value == 1:
            # exponential line broadening
            lb = dictionary["Apodization"]["Parameters"]["Line broadening (Hz)"]
            dimension_tab.apodization.exponential_line_broadening = float(lb)
        elif apodization_value == 2:
            # Lorentz to gauss
            g1 = dictionary["Apodization"]["Parameters"]["Inverse line broadening (Hz)"]
            g2 = dictionary["Apodization"]["Parameters"]["Gaussian broadening (Hz)"]
            g3 = dictionary["Apodization"]["Parameters"]["Gaussian shift"]

            dimension_tab.apodization.g1 = float(g1)
            dimension_tab.apodization.g2 = float(g2)
            dimension_tab.apodization.g3 = float(g3)

        elif apodization_value == 3:
            # Sinebell
            offset = dictionary["Apodization"]["Parameters"]["Offset (pi)"]
            end = dictionary["Apodization"]["Parameters"]["End (pi)"]
            power = dictionary["Apodization"]["Parameters"]["Power"]

            dimension_tab.apodization.offset = float(offset)
            dimension_tab.apodization.end = float(end)
            dimension_tab.apodization.power = float(power)

        elif apodization_value == 4:
            # gaussian broadening
            a = dictionary["Apodization"]["Parameters"]["Line broadening (Hz)"]
            b = dictionary["Apodization"]["Parameters"]["Gaussian broadening (Hz)"]

            dimension_tab.apodization.a = float(a)
            dimension_tab.apodization.b = float(b)

        elif apodization_value == 5:
            # trapezoid
            t1 = dictionary["Apodization"]["Parameters"]["Ramp up points"]
            t2 = dictionary["Apodization"]["Parameters"]["Ramp down points"]

            dimension_tab.apodization.t1 = int(t1)
            dimension_tab.apodization.t2 = int(t1)

        elif apodization_value == 6:
            # triangle
            loc = dictionary["Apodization"]["Parameters"]["Location of maximum"]
            dimension_tab.apodization.loc = float(loc)

        # Adding first point correction value
        c = dictionary["Apodization"]["Parameters"]["First point correction"]
        dimension_tab.apodization.apodization_first_point_textcontrol.SetValue(str(c))
        dimension_tab.apodization.apodization_first_point_scaling = float(c)

        dimension_tab.apodization.apodization_combobox.SetSelection(apodization_value)
        dimension_tab.apodization.on_apodization_combobox(wx.EVT_COMBOBOX)

    def load_zero_filling(self, dimension, dimension_tab, dictionary):
        """
        Reading in the saved parameters associated with zero filling.
        """

        zero_filling_flag = bool(dictionary["Zero filling"]["Zero filling flag"])
        zero_filling_selection = int(dictionary["Zero filling"]["Selection"])

        dimension_tab.zero_filling.zero_filling_combobox_selection = (
            zero_filling_selection
        )

        dimension_tab.zero_filling.zero_filling_checkbox_value = zero_filling_flag

        choices = ["Doubling number", "Zeros to add", "Final size"]
        choice = choices[zero_filling_selection]

        textbox_value = dictionary["Zero filling"]["Parameters"][choice]

        if zero_filling_selection == 0:
            dimension_tab.zero_filling.zero_filling_value_doubling_times = int(
                textbox_value
            )
        elif zero_filling_selection == 1:
            dimension_tab.zero_filling.zero_filling_value_zeros_to_add = int(
                textbox_value
            )
        elif zero_filling_selection == 2:
            dimension_tab.zero_filling.zero_filling_value_final_data_size = int(
                textbox_value
            )

        rounding = bool(
            dictionary["Zero filling"]["Parameters"]["Round to nearest power of 2"]
        )

        dimension_tab.zero_filling.zero_filling_checkbox.SetValue(zero_filling_flag)
        dimension_tab.zero_filling.zero_filling_combobox.SetSelection(
            zero_filling_selection
        )
        dimension_tab.zero_filling.zero_filling_textcontrol.SetValue(str(textbox_value))
        dimension_tab.zero_filling.zero_filling_round_checkbox.SetValue(rounding)

        dimension_tab.zero_filling.on_zero_filling_combobox(wx.EVT_COMBOBOX)

    def load_fourier_transform(self, dimension, dimension_tab, dictionary):
        """
        Reading in the saved parameters associated with fourier transform.
        """

        ft_flag = bool(dictionary["Fourier transform"]["Fourier transform flag"])
        ft_option = int(
            dictionary["Fourier transform"]["Fourier transform method selection"]
        )
        dimension_tab.fourier_transform.fourier_transform_checkbox.SetValue(ft_flag)
        dimension_tab.fourier_transform.ft_method_selection = ft_option

    def load_phasing(self, dimension, dimension_tab, dictionary):
        """
        Reading in the saved parameters associated with phase correction.
        """

        phasing_flag = bool(dictionary["Phasing"]["Phasing flag"])
        p0 = dictionary["Phasing"]["P0"]
        p1 = dictionary["Phasing"]["P1"]
        if dimension == 0:
            dimension_tab.phasing.phase_correction_checkbox.SetValue(phasing_flag)
            dimension_tab.phasing.phase_correction_checkbox_value = phasing_flag
            dimension_tab.phasing.p0_total = p0
            dimension_tab.phasing.phase_correction_p0_textcontrol.SetValue(str(p0))
            dimension_tab.phasing.p1_total = p1
            dimension_tab.phasing.phase_correction_p1_textcontrol.SetValue(str(p1))
            magnitude_mode = bool(dictionary["Phasing"]["Magnitude mode"])
            dimension_tab.phasing.magnitude_mode_checkbox.SetValue(magnitude_mode)
            dimension_tab.phasing.magnitude_mode_toggle = magnitude_mode

        else:
            dimension_tab.phasing.phase_correction_checkbox_indirect.SetValue(
                phasing_flag
            )
            dimension_tab.phasing.phase_correction_checkbox_value_indirect = (
                phasing_flag
            )
            dimension_tab.phasing.p0_total_indirect = p0
            dimension_tab.phasing.phase_correction_p0_textcontrol_indirect.SetValue(
                str(p0)
            )
            dimension_tab.phasing.p1_total_indirect = p1
            dimension_tab.phasing.phase_correction_p1_textcontrol_indirect.SetValue(
                str(p1)
            )
            f1180 = bool(dictionary["Phasing"]["f1180 flag"])
            dimension_tab.phasing.phase_correction_f1180_button_indirect.SetValue(f1180)

    def load_extraction(self, dimension, dimension_tab, dictionary):
        """
        Reading in the saved parameters associated with extraction
        """

        extraction_flag = bool(dictionary["Extraction"]["Extraction flag"])
        start = float(dictionary["Extraction"]["Start (ppm)"])
        end = float(dictionary["Extraction"]["End (ppm)"])

        dimension_tab.extraction.extraction_checkbox_value = extraction_flag
        dimension_tab.extraction.extraction_ppm_start = start
        dimension_tab.extraction.extraction_ppm_end = end

        dimension_tab.extraction.extraction_checkbox.SetValue(extraction_flag)
        dimension_tab.extraction.extraction_ppm_start_textcontrol.SetValue(str(start))
        dimension_tab.extraction.extraction_ppm_end_textcontrol.SetValue(str(end))

    def load_baseline_correction(self, dimension, dimension_tab, dictionary):
        """
        Reading in the saved parameters associated with baseline correction.
        """

        baseline_correction_flag = bool(
            dictionary["Baseline correction"]["Baseline correction flag"]
        )

        selection = int(dictionary["Baseline correction"]["Selection"])
        node_width = int(dictionary["Baseline correction"]["Parameters"]["Node width"])
        node_list = dictionary["Baseline correction"]["Parameters"]["Node list"]
        polynomial_order = int(
            dictionary["Baseline correction"]["Parameters"]["Polynomial order"]
        )

        dimension_tab.baseline_correction.baseline_correction_checkbox.SetValue(
            baseline_correction_flag
        )
        dimension_tab.baseline_correction.baseline_correction_radio_box.SetSelection(
            selection
        )
        dimension_tab.baseline_correction.node_width = node_width
        dimension_tab.baseline_correction.baseline_correction_nodes_textcontrol.SetValue(
            str(node_width)
        )
        dimension_tab.baseline_correction.node_list = node_list
        dimension_tab.baseline_correction.baseline_correction_node_list_textcontrol.SetValue(
            str(node_list)
        )
        dimension_tab.baseline_correction.polynomial_order = polynomial_order

        dimension_tab.baseline_correction.baseline_correction_polynomial_order_textcontrol.SetValue(
            str(polynomial_order)
        )

    def error_dimension_loading(self, dimension: int):
        """
        Outputting an error message if the dimension
        was not loaded correctly saying that some
        parameters may be default values.
        """
        message = "Parameters for dimension {} were not read correctly. Parameters may be default values".format(
            dimension
        )

        dlg = wx.MessageDialog(
            self.notebook,
            message,
            "Warning",
            wx.OK,
        )
        result = dlg.ShowModal()


def load_variables_from_nmrproc_com_1D(self):
    # Open processing_parameters.txt file and load the variables from it
    file = open("processing_parameters.txt", "r")
    lines = file.readlines()
    file.close()

    self.direct_solvent_suppression = False
    self.linear_prediction_checkbox_value = False
    self.apodization_checkbox_value = False
    self.zero_filling_checkbox_value = False
    self.fourier_transform_checkbox_value = False
    self.phase_correction_checkbox_value = False
    self.extraction_checkbox_value = False
    self.baseline_correction_checkbox_value = False

    include_line = False
    for line in lines:
        if "Dimension 1" in line:
            include_line = True
            continue
        if include_line == False:
            continue
        if include_line == True and "Dimension 2" in line:
            include_line = False
            break
        if include_line == True:
            line = line.split("\n")[0]
            if line.split(":")[0] == "Solvent Suppression":
                if line.split(": ")[1].strip() == "True" in line:
                    self.direct_solvent_suppression = True
                else:
                    self.direct_solvent_suppression = False
            elif line.split(":")[0] == "Filter Selection":
                self.solvent_suppression_filter_selection = int(line.split(": ")[1])
            elif line.split(":")[0] == "Lowpass Shape Selection":
                self.solvent_suppression_lowpass_shape_selection = int(
                    line.split(": ")[1]
                )
            elif line.split(":")[0] == "Linear Prediction":
                if line.split(": ")[1].strip() == "True":
                    self.linear_prediction_checkbox_value = True
                else:
                    self.linear_prediction_checkbox_value = False
            elif line.split(":")[0] == "Linear Prediction Options Selection":
                self.linear_prediction_options_selection = int(line.split(": ")[1])
            elif line.split(":")[0] == "Linear Prediction Coefficients Selection":
                self.linear_prediction_coefficients_selection = int(line.split(": ")[1])
            elif line.split(":")[0] == "Apodization":
                if line.split(": ")[1].strip() == "True" in line:
                    self.apodization_checkbox_value = True
                else:
                    self.apodization_checkbox_value = False
            elif line.split(":")[0] == "Apodization Combobox Selection":
                self.apodization_combobox_selection = int(line.split(": ")[1])
                self.apodization_combobox_selection_old = int(line.split(": ")[1])
            elif line.split(":")[0] == "Exponential Line Broadening":
                self.exponential_line_broadening = float(line.split(": ")[1])
            elif line.split(":")[0] == "Apodization First Point Scaling":
                self.apodization_first_point_scaling = float(line.split(": ")[1])
            elif line.split(":")[0] == "G1":
                self.g1 = float(line.split(": ")[1])
            elif line.split(":")[0] == "G2":
                self.g2 = float(line.split(": ")[1])
            elif line.split(":")[0] == "G3":
                self.g3 = float(line.split(": ")[1])
            elif line.split(":")[0] == "Offset":
                self.offset = float(line.split(": ")[1])
            elif line.split(":")[0] == "End":
                self.end = float(line.split(": ")[1])
            elif line.split(":")[0] == "Power":
                self.power = float(line.split(": ")[1])
            elif "A:" in line:
                self.a = float(line.split(": ")[1])
            elif "B:" in line:
                self.b = float(line.split(": ")[1])
            elif line.split(":")[0] == "T1":
                self.t1 = float(line.split(": ")[1])
            elif line.split(":")[0] == "T2":
                self.t2 = float(line.split(": ")[1])
            elif line.split(":")[0] == "Loc":
                self.loc = float(line.split(": ")[1])
            elif line.split(":")[0] == "Zero Filling":
                if "True" in line:
                    self.zero_filling_checkbox_value = True
                else:
                    self.zero_filling_checkbox_value = False
            elif line.split(":")[0] == "Zero Filling Combobox Selection":
                self.zero_filling_combobox_selection = int(line.split(": ")[1])
            elif line.split(":")[0] == "Zero Filling Value Doubling Times":
                self.zero_filling_value_doubling_times = int(line.split(": ")[1])
            elif line.split(":")[0] == "Zero Filling Value Zeros to Add":
                self.zero_filling_value_zeros_to_add = int(line.split(": ")[1])
            elif line.split(":")[0] == "Zero Filling Value Final Data Size":
                self.zero_filling_value_final_data_size = int(line.split(": ")[1])
            elif line.split(":")[0] == "Zero Filling Round Checkbox":
                if "True" in line:
                    self.zero_filling_round_checkbox_value = True
                else:
                    self.zero_filling_round_checkbox_value = False
            elif line.split(":")[0] == "Fourier Transform":
                if "True" in line:
                    self.fourier_transform_checkbox_value = True
                else:
                    self.fourier_transform_checkbox_value = False
            elif line.split(":")[0] == "Fourier Transform Method Selection":
                self.ft_method_selection = int(line.split(": ")[1])
            elif line.split(":")[0] == "Phase Correction":
                if "True" in line:
                    self.phase_correction_checkbox_value = True
                else:
                    self.phase_correction_checkbox_value = False
            elif line.split(":")[0] == "Phase Correction P0":
                self.p0_total = float(line.split(": ")[1])
            elif line.split(":")[0] == "Phase Correction P1":
                self.p1_total = float(line.split(": ")[1])
            elif line.split(":")[0] == "Magnitude Mode":
                if "True" in line:
                    self.magnitude_mode_toggle = True
                else:
                    self.magnitude_mode_toggle = False
            elif line.split(":")[0] == "Extraction":
                if "True" in line:
                    self.extraction_checkbox_value = True
                else:
                    self.extraction_checkbox_value = False
            elif line.split(":")[0] == "Extraction PPM Start":
                self.extraction_ppm_start = float(line.split(": ")[1])
            elif line.split(":")[0] == "Extraction PPM End":
                self.extraction_ppm_end = float(line.split(": ")[1])
            elif line.split(":")[0] == "Baseline Correction":
                if "True" in line:
                    self.baseline_correction_checkbox_value = True
                else:
                    self.baseline_correction_checkbox_value = False
            elif line.split(":")[0] == "Baseline Correction Radio Box Selection":
                self.baseline_correction_radio_box_selection = int(line.split(": ")[1])
            elif line.split(":")[0] == "Baseline Correction Nodes":
                self.baseline_correction_nodes = int(line.split(": ")[1])
            elif line.split(":")[0] == "Baseline Correction Node List":
                self.baseline_correction_node_list = line.split(": ")[1]
            elif line.split(":")[0] == "Baseline Correction Polynomial Order":
                self.baseline_correction_polynomial_order = int(line.split(": ")[1])
