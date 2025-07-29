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
from typing import Dict, Any, Union


class Save_json:

    def __init__(self, notebook, nmrdata, dimension_tabs):
        """
        This object will obtain the current parameters using the populate_dictionary
        classes and then save the resulting dictionary to a .json file
        """

        self.notebook = notebook

        continue_saving = self.check_saved_json()
        if continue_saving == False:
            return
        else:
            # Create dictionary
            parameter_dictionary_class = Populate_dictionary_global(
                nmrdata, dimension_tabs
            )
            parameter_dictionary = parameter_dictionary_class.parameter_dictionary
            if continue_saving == True:
                # checking to see if conversion parameters are saved
                # as well as processing parameters
                saved_conversion_params = self.check_saved_conversion_params()
                try:
                    parameter_dictionary["conversion"] = saved_conversion_params
                except:
                    self.error_conversion_parameters()

            self.write_json(parameter_dictionary)

    def write_json(self, parameter_dictionary):
        """
        Function to write the parameter dictionary to a .json file
        """
        filename = "parameters.json"
        with open(filename, "w") as file:
            json.dump(
                parameter_dictionary,
                file,
                indent=4,
            )

    def check_saved_json(self) -> Union[bool, None]:
        """
        Check to see if there are any saved converter.json files
        Ask the user if they want to overwrite the converter.json
        files
        """

        if pathlib.Path("parameters.json").exists() == True:
            # Asking the user if they want to overwrite the previous saved state
            dlg = wx.MessageDialog(
                self.notebook,
                "A previous set of saved parameters has been found (parameters.json). Would you like to overwrite this?",
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
        else:
            return None

    def check_saved_conversion_params(self) -> Union[Dict[str, Any], bool]:
        """
        Checking the processing.json file to see if there are connversion
        parameters saved as well as processing parameters.
        """

        try:
            with open("parameters.json", "r") as file:
                dictionary = json.load(file)
            if "conversion" in dictionary.keys():
                return dictionary["conversion"]

        except:
            # Unable to read the converter.json file effectively
            # Giving a warning to the user and returning
            message = (
                "Unable to read parameters.json. Saving without processing parameters."
            )
            dlg = wx.MessageDialog(
                self.notebook, message, "Warning", wx.OK | wx.ICON_ERROR
            )
            dlg.ShowModal()
            dlg.Destroy()

            return False

    def error_conversion_parameters(self):
        """
        Will be activated if conversion parameters are found but could
        not be read correctly. The conversion parameters will not be
        included in the parameters.json file.
        """
        message = "Conversion parameters were found in parameters.json but could not be read correctly. Continuing"
        dlg = wx.MessageDialog(self.notebook, message, "Warning", wx.OK | wx.ICON_ERROR)
        dlg.ShowModal()
        dlg.Destroy()


class Populate_dictionary_global:

    def __init__(self, nmrdata, dimension_tabs):
        """
        This class will find all the current parameters from the SpinConverter
        GUI which can then be saved.
        """

        self.nmrdata = nmrdata

        # Creating a dictionary to store the current parameters and adding
        # general information
        parameter_dictionary = self.create_dictionary()

        # For each dimension, add the processing parameters
        for i, dimension_tab in enumerate(dimension_tabs):
            parameter_dictionary = self.add_dimension_parameters(
                parameter_dictionary, i, dimension_tab
            )

        self.parameter_dictionary = parameter_dictionary

    def create_dictionary(self) -> Dict:
        """
        Create a dictionary to hold the relevant parameters for the
        processing parameters
        """
        dictionary = {}
        dictionary["processing"] = {}

        return dictionary

    def add_dimension_parameters(
        self, parameter_dictionary: Dict[str, Any], dimension: int, dimension_tab
    ) -> Dict[str, Any]:
        """
        This function will add the parameters associated with a dimension to the
        dictionary. The dictionary element will be named by the dimension number.
        """
        label = (
            "Dimension {}".format(dimension)
            + " ("
            + str(self.nmrdata.axislabels[dimension])
            + ")"
        )
        parameter_dictionary["processing"][label] = {}

        parameter_dictionary["processing"][label] = self.add_parameters(
            parameter_dictionary["processing"][label], dimension, dimension_tab
        )
        return parameter_dictionary

    def add_parameters(
        self, dimension_dictionary: Dict, dimension: int, dimension_tab
    ) -> Dict[str, Any]:
        """
        This function will add all the current processing parameters for the
        dimension to the dictionary.
        """
        if dimension == 0:
            dimension_dictionary = self.add_solvent_suppression_parameters(
                dimension_dictionary, dimension, dimension_tab
            )
        dimension_dictionary = self.add_linear_prediction_parameters(
            dimension_dictionary, dimension, dimension_tab
        )
        dimension_dictionary = self.add_apodization_parameters(
            dimension_dictionary, dimension, dimension_tab
        )
        dimension_dictionary = self.add_zero_filling_parameters(
            dimension_dictionary, dimension, dimension_tab
        )
        dimension_dictionary = self.add_fourier_transform_parameters(
            dimension_dictionary, dimension, dimension_tab
        )
        dimension_dictionary = self.add_phasing_parameters(
            dimension_dictionary, dimension, dimension_tab
        )
        dimension_dictionary = self.add_extraction_parameters(
            dimension_dictionary, dimension, dimension_tab
        )
        dimension_dictionary = self.add_baseline_correction_parameters(
            dimension_dictionary, dimension, dimension_tab
        )

        return dimension_dictionary

    def add_solvent_suppression_parameters(
        self, dimension_dictionary: Dict, dimension: int, dimension_tab
    ) -> Dict[str, Any]:
        """
        Adding the current solvent suppression parameters in SpinProcess
        to the dictionary.
        """

        dimension_dictionary["Solvent Suppression"] = {}
        dimension_dictionary["Solvent Suppression"][
            "Solvent Suppression Flag"
        ] = dimension_tab.solvent_suppression.solvent_suppression_checkbox.GetValue()

        filter_choices = ["Low-pass", "Spline", "Polynomial"]
        filter_index = (
            dimension_tab.solvent_suppression.solvent_suppression_filter_selection
        )
        dimension_dictionary["Solvent Suppression"]["Filter Selection"] = [
            filter_index,
            filter_choices[filter_index],
        ]

        lowpass_choices = ["Boxcar", "Sine", "Sine Squared"]
        lowpass_index = (
            dimension_tab.solvent_suppression.solvent_suppression_lowpass_shape_selection
        )
        dimension_dictionary["Solvent Suppression"]["Lowpass Shape"] = [
            lowpass_index,
            lowpass_choices[lowpass_index],
        ]

        return dimension_dictionary

    def add_linear_prediction_parameters(
        self, dimension_dictionary: Dict, dimension: int, dimension_tab
    ) -> Dict[str, Any]:
        """
        Adding the current linear prediction or NUS parameters in SpinProcess
        to the dictionary for dimension.
        """

        dimension_dictionary["Linear Prediction"] = {}
        if dimension == 0:
            dimension_dictionary["Linear Prediction"][
                "Linear Prediction Flag"
            ] = dimension_tab.linear_prediction.linear_prediction_checkbox.GetValue()
            # Add predicted points (0 is after FID, 1 is before FID)
            options = ["After FID", "Before FID"]
            value = (
                dimension_tab.linear_prediction.linear_prediction_combobox.GetSelection()
            )
            dimension_dictionary["Linear Prediction"]["Add predicted points"] = [
                value,
                options[value],
            ]
            # Predicted coefficients (0 is forward, 1 is backward, 2 is both)
            options = ["Forward", "Backward", "Both"]
            value = (
                dimension_tab.linear_prediction.linear_prediction_coefficients_combobox.GetSelection()
            )
            dimension_dictionary["Linear Prediction"]["Predicted coefficients"] = [
                value,
                options[value],
            ]
        else:
            # Linear prediction choices (0 is none, 1 is linear prediction, 2 is NUS reconstruction)
            value = (
                dimension_tab.linear_prediction.linear_prediction_radio_box_indirect.GetSelection()
            )
            choices = ["None", "Linear Prediction", "NUS Reconstruction"]
            dimension_dictionary["Linear Prediction"][choices[value]] = {}
            if value == 1:
                # Add predicted points (0 is after FID, 1 is before FID)
                options = ["After FID", "Before FID"]
                value1 = (
                    dimension_tab.linear_prediction.linear_prediction_combobox_indirect.GetSelection()
                )
                dimension_dictionary["Linear Prediction"][choices[value]][
                    "Add predicted points"
                ] = [value1, options[value1]]
                # Predicted coefficients (0 is forward, 1 is backward, 2 is both)
                options = ["Forward", "Backward", "Both"]
                value2 = (
                    dimension_tab.linear_prediction.linear_prediction_coefficients_combobox_indirect.GetSelection()
                )
                dimension_dictionary["Linear Prediction"][choices[value]][
                    "Predicted coefficients"
                ] = [value2, options[value2]]
            elif value == 2:
                nusfile = dimension_tab.linear_prediction.nuslist_name_indirect
                nus_extension = (
                    dimension_tab.linear_prediction.smile_data_extension_number_indirect
                )
                nus_cpu = dimension_tab.linear_prediction.number_of_nus_CPU_indirect
                nus_iterations = dimension_tab.linear_prediction.nus_iterations_indirect
                dimension_dictionary["Linear Prediction"][choices[value]][
                    "NUS file"
                ] = nusfile
                dimension_dictionary["Linear Prediction"][choices[value]][
                    "NUS extension"
                ] = nus_extension
                dimension_dictionary["Linear Prediction"][choices[value]][
                    "NUS CPU number"
                ] = nus_cpu
                dimension_dictionary["Linear Prediction"][choices[value]][
                    "NUS iterations"
                ] = nus_iterations

        return dimension_dictionary

    def add_apodization_parameters(
        self, dimension_dictionary: Dict[str, Any], dimension: int, dimension_tab
    ) -> Dict[str, Any]:
        """
        Adding the current apodization parameters in SpinProcess
        to the dictionary for dimension.
        """

        apodization_value = (
            dimension_tab.apodization.apodization_combobox.GetSelection()
        )
        apodization_type = dimension_tab.apodization.apodization_combobox.GetValue()
        checkbox_value = dimension_tab.apodization.apodization_checkbox_value

        dimension_dictionary["Apodization"] = {}
        dimension_dictionary["Apodization"]["Apodization flag"] = checkbox_value
        dimension_dictionary["Apodization"]["Selection"] = apodization_value
        dimension_dictionary["Apodization"]["Type"] = apodization_type
        dimension_dictionary["Apodization"]["Parameters"] = {}

        if apodization_value == 1:
            # exponential
            line_broadening = dimension_tab.apodization.exponential_line_broadening
            dimension_dictionary["Apodization"]["Parameters"][
                "Line broadening (Hz)"
            ] = line_broadening
        elif apodization_value == 2:
            # Lorentz to gauss
            g1 = dimension_tab.apodization.g1
            g2 = dimension_tab.apodization.g2
            g3 = dimension_tab.apodization.g3
            dimension_dictionary["Apodization"]["Parameters"][
                "Inverse line broadening (Hz)"
            ] = g1
            dimension_dictionary["Apodization"]["Parameters"][
                "Gaussian broadening (Hz)"
            ] = g2
            dimension_dictionary["Apodization"]["Parameters"]["Gaussian shift"] = g3
        elif apodization_value == 3:
            # Sinebell]
            offset = dimension_tab.apodization.offset
            end = dimension_tab.apodization.end
            power = dimension_tab.apodization.power
            dimension_dictionary["Apodization"]["Parameters"]["Offset (pi)"] = offset
            dimension_dictionary["Apodization"]["Parameters"]["End (pi)"] = end
            dimension_dictionary["Apodization"]["Parameters"]["Power"] = power
        elif apodization_value == 4:
            # Gaussian broadening
            line_broadening = dimension_tab.apodization.a
            gaussian_broadening = dimension_tab.apodization.b
            dimension_dictionary["Apodization"]["Parameters"][
                "Line broadening (Hz)"
            ] = line_broadening
            dimension_dictionary["Apodization"]["Parameters"][
                "Gaussian broadening (Hz)"
            ] = gaussian_broadening
        elif apodization_value == 5:
            # Trapezoid
            ramp_up = dimension_tab.apodization.t1
            ramp_down = dimension_tab.apodization.t2
            dimension_dictionary["Apodization"]["Parameters"][
                "Ramp up points"
            ] = ramp_up
            dimension_dictionary["Apodization"]["Parameters"][
                "Ramp down points"
            ] = ramp_down
        elif apodization_value == 6:
            loc = dimension_tab.apodization.loc
            dimension_dictionary["Apodization"]["Parameters"][
                "Location of maximum"
            ] = loc

        # Adding first point correction value
        c = dimension_tab.apodization.apodization_first_point_textcontrol.GetValue()
        dimension_dictionary["Apodization"]["Parameters"]["First point correction"] = c

        return dimension_dictionary

    def add_zero_filling_parameters(
        self, dimension_dictionary: Dict[str, Any], dimension: int, dimension_tab
    ) -> Dict[str, Any]:
        """
        Adding the current zero filling parameters in SpinProcess
        to the dictionary for dimension.
        """

        zero_filling_flag = dimension_tab.zero_filling.zero_filling_checkbox.GetValue()
        dimension_dictionary["Zero filling"] = {}
        dimension_dictionary["Zero filling"]["Zero filling flag"] = zero_filling_flag

        zero_filling_selection = (
            dimension_tab.zero_filling.zero_filling_combobox.GetSelection()
        )
        zero_filling_type = dimension_tab.zero_filling.zero_filling_combobox.GetValue()

        dimension_dictionary["Zero filling"]["Selection"] = zero_filling_selection
        dimension_dictionary["Zero filling"]["Type"] = zero_filling_type

        dimension_dictionary["Zero filling"]["Parameters"] = {}
        choices = ["Doubling number", "Zeros to add", "Final size"]
        value = dimension_tab.zero_filling.zero_filling_textcontrol.GetValue()
        dimension_dictionary["Zero filling"]["Parameters"][
            choices[zero_filling_selection]
        ] = value

        zero_filling_round = (
            dimension_tab.zero_filling.zero_filling_round_checkbox.GetValue()
        )
        dimension_dictionary["Zero filling"]["Parameters"][
            "Round to nearest power of 2"
        ] = zero_filling_round

        return dimension_dictionary

    def add_fourier_transform_parameters(
        self, dimension_dictionary: Dict[str, Any], dimension: int, dimension_tab
    ) -> Dict[str, Any]:
        """
        Adding the current fourier transform parameters in SpinProcess
        to the dictionary for dimension.
        """

        dimension_dictionary["Fourier transform"] = {}
        ft_flag = dimension_tab.fourier_transform.fourier_transform_checkbox.GetValue()
        dimension_dictionary["Fourier transform"]["Fourier transform flag"] = ft_flag

        ft_option = int(dimension_tab.fourier_transform.ft_method_selection)
        methods = [
            "Standard",
            "Auto (not recommended)",
            "Real",
            "Inverse",
            "Sign alternation (alt)",
            "Negate imaginaries (neg)",
            "alt + neg",
        ]
        dimension_dictionary["Fourier transform"][
            "Fourier transform method selection"
        ] = ft_option
        dimension_dictionary["Fourier transform"]["Fourier transform method type"] = (
            methods[ft_option]
        )

        return dimension_dictionary

    def add_phasing_parameters(
        self, dimension_dictionary: Dict[str, Any], dimension: int, dimension_tab
    ) -> Dict[str, Any]:
        """
        Adding the current phasing parameters in SpinProcess
        to the dictionary for dimension. If dimension is 0 (direct),
        then have an extra magnitude mode option. If dimension is
        greater than 0 (indirect), then there is an extra F1180 option.
        """

        dimension_dictionary["Phasing"] = {}
        if dimension == 0:
            phasing_flag = dimension_tab.phasing.phase_correction_checkbox.GetValue()
            p0 = dimension_tab.phasing.p0_total
            p1 = dimension_tab.phasing.p1_total
            magnitude_mode = dimension_tab.phasing.magnitude_mode_checkbox.GetValue()
            dimension_dictionary["Phasing"]["Magnitude mode"] = magnitude_mode
        else:
            phasing_flag = (
                dimension_tab.phasing.phase_correction_checkbox_indirect.GetValue()
            )
            p0 = dimension_tab.phasing.p0_total_indirect
            p1 = dimension_tab.phasing.p1_total_indirect
            f1180 = (
                dimension_tab.phasing.phase_correction_f1180_button_indirect.GetValue()
            )
            dimension_dictionary["Phasing"]["f1180 flag"] = f1180

        dimension_dictionary["Phasing"]["Phasing flag"] = phasing_flag
        dimension_dictionary["Phasing"]["P0"] = p0
        dimension_dictionary["Phasing"]["P1"] = p1

        return dimension_dictionary

    def add_extraction_parameters(
        self, dimension_dictionary: Dict[str, Any], dimension: int, dimension_tab
    ) -> Dict[str, Any]:
        """
        Adding the current extraction parameters in SpinProcess
        to the dictionary for dimension.
        """

        dimension_dictionary["Extraction"] = {}

        extraction_flag = dimension_tab.extraction.extraction_checkbox_value
        start = dimension_tab.extraction.extraction_ppm_start
        end = dimension_tab.extraction.extraction_ppm_end

        dimension_dictionary["Extraction"]["Extraction flag"] = extraction_flag
        dimension_dictionary["Extraction"]["Start (ppm)"] = start
        dimension_dictionary["Extraction"]["End (ppm)"] = end

        return dimension_dictionary

    def add_baseline_correction_parameters(
        self, dimension_dictionary: Dict[str, Any], dimension: int, dimension_tab
    ) -> Dict[str, Any]:
        """
        Adding the current baseline correction parameters in SpinProcess
        to the dictionary for dimension.
        """

        dimension_dictionary["Baseline correction"] = {}

        baseline_correction_flag = (
            dimension_tab.baseline_correction.baseline_correction_checkbox.GetValue()
        )
        selection = (
            dimension_tab.baseline_correction.baseline_correction_radio_box.GetSelection()
        )
        method = (
            dimension_tab.baseline_correction.baseline_correction_radio_box.GetStringSelection()
        )
        node_width = dimension_tab.baseline_correction.node_width
        node_list = dimension_tab.baseline_correction.node_list
        polynomial_order = dimension_tab.baseline_correction.polynomial_order

        dimension_dictionary["Baseline correction"][
            "Baseline correction flag"
        ] = baseline_correction_flag
        dimension_dictionary["Baseline correction"]["Selection"] = selection
        dimension_dictionary["Baseline correction"]["Method"] = method
        dimension_dictionary["Baseline correction"]["Parameters"] = {}
        dimension_dictionary["Baseline correction"]["Parameters"][
            "Node width"
        ] = node_width
        dimension_dictionary["Baseline correction"]["Parameters"][
            "Node list"
        ] = node_list
        dimension_dictionary["Baseline correction"]["Parameters"][
            "Polynomial order"
        ] = polynomial_order

        return dimension_dictionary
