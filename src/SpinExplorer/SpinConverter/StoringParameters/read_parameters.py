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

import json
import pathlib
import wx
from typing import Union, Dict, Any


class Read_json:
    def __init__(self, params, nmrdata, app):
        """
        This class will read a converter.json file and ensure that all of the
        parameters are valid (if not will go back to original read parameters).
        It will then load all the correct parameters etc into the
        SpinConverter app.
        """

        self.app = app

        # Checking to see if parameters.json is found and if user wants to
        # read in these values
        if self.find_json() == True:
            parameter_dictionary = self.read()
            if parameter_dictionary == False:
                return
            else:
                if self.check_data_mismatch(parameter_dictionary) == True:
                    # Data is mismatched, returning
                    return
                parameter_inputs = InputParameters(
                    parameter_dictionary, params, nmrdata, app
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
            self.app,
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
                parameter_dictionary = json.load(file)["conversion"]
            return parameter_dictionary
        except:
            # Unable to read the converter.json file effectively
            # Giving a warning to the user and returning
            message = "Unable to read the saved conversion parameters in parameters.json. Continuing using default parameters."
            dlg = wx.MessageDialog(self.app, message, "Warning", wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()

            return False

    def check_data_mismatch(self, parameter_dictionary: Dict[str, Any]) -> bool:
        """
        Checking for a data mismatch with the parameters.json file.
        i.e. Bruker data with Varian parameters.json file and vice
        versa.
        """
        if parameter_dictionary["general"]["spectrometer"] == "Bruker":
            if pathlib.Path("acqus").exists():
                # No data mismatch
                return False
            else:
                if pathlib.Path("procpar").exists():
                    self.data_mismatch_warning("Bruker", "Varian")
                else:
                    self.data_mismatch_warning("Bruker", "Unknown")
                return True

        elif parameter_dictionary["general"]["spectrometer"] == "Varian":
            if pathlib.Path("procpar").exists():
                # No data mismatch
                return False
            else:
                if pathlib.Path("acqus").exists():
                    self.data_mismatch_warning("Varian", "Bruker")
                else:
                    self.data_mismatch_warning("Varian", "Unknown")

                return True

        else:
            self.unknown_spectrometer_type_warning(
                parameter_dictionary["general"]["spectrometer"]
            )
            return True

    def data_mismatch_warning(self, spectrometer_type, data_type):
        """
        Giving an error message to the user that the spectrometer
        type in parameters.json and apparent data type do not match.
        """
        message = (
            "The spectrometer type ("
            + spectrometer_type
            + ") and apparent data type ("
            + data_type
            + ") do not match. Continuing with default parameters."
        )
        dlg = wx.MessageDialog(self.app, message, "Warning", wx.OK | wx.ICON_ERROR)
        dlg.ShowModal()
        dlg.Destroy()

    def unknown_spectrometer_type_warning(self, spectrometer_type):
        """
        Giving an error message to the user that the spectrometer
        type in parameters.json is unknown
        """
        message = (
            "The spectrometer type ("
            + spectrometer_type
            + ") is unknown. Continuing with default parameters."
        )
        dlg = wx.MessageDialog(self.app, message, "Warning", wx.OK | wx.ICON_ERROR)
        dlg.ShowModal()
        dlg.Destroy()


class InputParameters:
    def __init__(self, parameter_dictionary, params, nmrdata, app):
        """
        This class will insert all of the read parameters into the
        SpinConverter interface.
        """

        self.parameter_dictionary = parameter_dictionary
        self.params = params
        self.nmrdata = nmrdata
        self.app = app

        self.input_all_parameters()
        # try:
        #     self.input_all_parameters()
        # except:
        #     self.input_error()

    def input_all_parameters(self):
        """
        Function which runs all the functions to input the saved
        conversion parameters into the GUI.
        """

        # Input acquisition and referencing modes, then sizes etc
        self.input_acquisition_modes()
        self.input_nucleus_frequencies()
        self.input_temperature()
        self.input_carrier()
        self.input_complex_sizes()
        self.input_real_sizes()
        self.input_sweep_widths()
        self.input_labels()
        self.input_intensity_scaling()

        # Additional Bruker parameters
        if self.parameter_dictionary["general"]["spectrometer"] == "Bruker":
            self.input_digital_filter()
            self.input_other_options()

        self.input_nus()

    def input_error(self):
        """
        Giving an error message to the user that the conversion
        parameters were not correctly imported to the GUI.
        Continuing with default parameters.
        """
        message = "The conversion parameters were not correctly imported to the GUI. Continuing with default parameters."
        dlg = wx.MessageDialog(self.app, message, "Warning", wx.OK | wx.ICON_ERROR)
        dlg.ShowModal()
        dlg.Destroy()

    def input_acquisition_modes(self):
        """
        Inputting the saved acqusition modes into the GUI
        """
        acqusition_modes = self.parameter_dictionary["spectral parameters"][
            "acqusition modes"
        ]
        direct_mode_index = acqusition_modes["direct"]["index"]
        direct_mode = acqusition_modes["direct"]["mode"]

        self.app.format.acqusition_combo_boxes[0].SetSelection(int(direct_mode_index))
        self.app.format.acqusition_combo_boxes[0].SetValue(direct_mode)

        if acqusition_modes["indirect"] != {}:
            for i, mode in enumerate(acqusition_modes["indirect"]["mode"]):
                index = acqusition_modes["indirect"]["index"][i]

                self.app.format.acqusition_combo_boxes[i + 1].SetSelection(int(index))
                evt = wx.CommandEvent(
                    wx.EVT_COMBOBOX.typeId,
                    self.app.format.acqusition_combo_boxes[i + 1].GetId(),
                )
                evt.SetEventObject(self.app.format.acqusition_combo_boxes[i + 1])
                self.app.shared_format.on_acquisition_mode_change(evt)
                self.app.format.acqusition_combo_boxes[i + 1].SetValue(mode)

    def input_temperature(self):
        """
        Input the saved temperature value and update the corresponding
        carrier of water.
        """
        # Spawning the temperature change box and pressing the buttons silently
        # without user input
        self.app.shared_format.on_temperature_change_button(wx.EVT_BUTTON)
        self.app.shared_format.temperature_input.SetValue(
            str(self.parameter_dictionary["general"]["temperature"])
        )
        self.app.shared_format.on_save_and_exit_temperature_change(wx.EVT_BUTTON)

    def input_carrier(self):
        """
        Inputting the referencing mode selection for each dimension
        and then adding the saved carrier value for each dimension.
        """
        carrier_frequencies = self.parameter_dictionary["spectral parameters"][
            "carrier frequencies"
        ]["frequency"]
        carrier_frequency_combobox_selection = self.parameter_dictionary[
            "spectral parameters"
        ]["carrier frequencies"]["combobox index"]

        for i, frequency in enumerate(carrier_frequencies):
            self.app.format.carrier_combo_boxes[i].SetSelection(
                int(carrier_frequency_combobox_selection[i])
            )
            evt = wx.CommandEvent(
                wx.EVT_COMBOBOX.typeId,
                self.app.format.carrier_combo_boxes[i].GetId(),
            )
            evt.SetEventObject(self.app.format.carrier_combo_boxes[i])
            self.app.format.on_carrier_combo_box_change(evt)
            self.app.format.carrier_frequency_boxes[i].SetValue(str(frequency))

    def input_complex_sizes(self):
        """
        Inputting the saved complex sizes for each dimension.
        """
        sizes = self.parameter_dictionary["spectral parameters"]["sizes"]["complex"]
        for i, size in enumerate(sizes):
            self.app.format.N_complex_boxes[i].SetValue(str(size))

    def input_real_sizes(self):
        """
        Inputting the saved real sizes for each dimension.
        """
        sizes = self.parameter_dictionary["spectral parameters"]["sizes"]["real"]
        for i, size in enumerate(sizes):
            self.app.format.N_real_boxes[i].SetValue(str(size))

    def input_sweep_widths(self):
        """
        Inputting the saved sweep with values for each dimension.
        """
        sweep_widths = self.parameter_dictionary["spectral parameters"]["sweep widths"][
            "values"
        ]
        sweep_widths_selections = self.parameter_dictionary["spectral parameters"][
            "sweep widths"
        ]["indexes"]
        for i, sw in enumerate(sweep_widths):
            self.app.format.sweep_width_boxes[i].SetSelection(
                int(sweep_widths_selections[i])
            )
            self.app.format.sweep_width_boxes[i].SetValue(str(sw))

    def input_nucleus_frequencies(self):
        """
        Inputting the saved nucleus frequencies for each dimension.
        """
        nucleus_frequencies = self.parameter_dictionary["spectral parameters"][
            "nuclei frequencies"
        ]
        for i, frequency in enumerate(nucleus_frequencies):
            self.app.format.nuclei_frequency_boxes[i].SetValue(str(frequency))

    def input_labels(self):
        """
        Inputting the saved labels for each dimension
        """
        labels = self.parameter_dictionary["spectral parameters"]["labels"]
        for i, label in enumerate(labels):
            self.app.format.nucleus_type_boxes[i].SetValue(str(label))

    def input_intensity_scaling(self):
        """
        Inputting the saved intensity scaling values.
        """
        intensity_dict = self.parameter_dictionary["intensity scaling"]
        scale_by_ns_flag = bool(intensity_dict["Scale by number of scans (NS)"])
        if self.parameter_dictionary["general"]["spectrometer"] == "Bruker":
            scale_by_nc_flag = bool(
                intensity_dict["Scale by Bruker normalisation constant (NC)"]
            )
        scale_by_1000_flag = bool(intensity_dict["Scale by 1000"])
        scaling_number = intensity_dict["Scaling number"]

        self.app.shared_format.scaling_NS_checkbox.SetValue(scale_by_ns_flag)
        if self.parameter_dictionary["general"]["spectrometer"] == "Bruker":
            self.app.shared_format.scaling_NC.SetValue(scale_by_nc_flag)
        self.app.shared_format.scaling_by_number.SetValue(scale_by_1000_flag)
        self.app.shared_format.scaling_number.SetValue(scaling_number)

    def input_digital_filter(self):
        """
        Adding the saved digital filter option for saved Bruker
        data.
        """
        digital_filter_dict = self.parameter_dictionary["digital filter parameters"]
        digital_filter_removal = bool(digital_filter_dict["Remove Digital Filter"])
        self.app.format.digital_filter_checkbox.SetValue(digital_filter_removal)
        self.app.format.on_digital_filter_checkbox(wx.EVT_CHECKBOX)

        # Removal of digital filter before/after fourier transform
        removal = digital_filter_dict["Remove Before/After Fourier Transform"]
        if removal == "Before":
            removal_value = 0
        else:
            removal_value = 1

        decim = digital_filter_dict["Decimation Rate (decim)"]
        dspfvs = digital_filter_dict["DSP Firmware Version (dspfvs)"]
        grpdly = digital_filter_dict["Group Delay (grpdly)"]
        self.params.grpdly = float(grpdly)
        self.params.decim = float(decim)
        self.params.dspfvs = int(dspfvs)
        self.app.format.decim_value = float(decim)
        self.app.format.dspfvs_value = int(dspfvs)
        self.app.format.grpdly_value = float(grpdly)

        if digital_filter_removal == True:
            self.app.format.decim_textbox.SetValue(str(decim))
            self.app.format.dspfvs_textbox.SetValue(str(dspfvs))
            self.app.format.grpdly_textbox.SetValue(str(grpdly))

            self.app.format.digital_filter_radio_box.SetSelection(removal_value)

    def input_other_options(self):
        """
        Adding the other options for saved Bruker data
        """
        other_params_dict = self.parameter_dictionary["other parameters"]

        # Bad point removal threshold
        self.params.bad_point_threshold = float(
            other_params_dict["bad point threshold"]
        )
        self.bad_point_threshold_value = self.params.bad_point_threshold
        self.app.format.bad_point_threshold_textbox.SetValue(
            str(self.bad_point_threshold_value)
        )

        # Remove acquisition padding
        self.params.remove_acquisition_padding = bool(
            other_params_dict["remove acqusition padding"]
        )
        self.app.format.remove_acquisition_padding_checkbox.SetValue(
            self.params.remove_acquisition_padding
        )
        self.app.format.on_remove_acquisition_padding_checkbox(wx.EVT_CHECKBOX)

    def input_nus(self):
        """
        Adding the saved nus options if these are in the saved dictionary
        """
        nus_dic = self.parameter_dictionary["NUS information"]
        if nus_dic == "N/A":
            return
        else:
            self.app.shared_format.NUS_tickbox.SetValue(True)
            self.app.shared_format.nusfile = nus_dic["NUS file"]
            self.app.shared_format.nusfile_input.SetValue(str(nus_dic["NUS file"]))
            self.app.shared_format.nus_offset_box.SetValue(str(nus_dic["NUS offset"]))
            self.app.shared_format.NUS_sample_count_box.SetValue(
                str(nus_dic["NUS sample count"])
            )
            self.app.shared_format.reverse_NUS_tickbox.SetValue(
                bool(nus_dic["Reverse NUS schedule"])
            )
