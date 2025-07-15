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
from typing import Dict, Any, Union


class Save_json:

    def __init__(self, params, nmrdata, app):
        """
        This object will obtain the current parameters using the populate_dictionary
        classes and then save the resulting dictionary to a .json file
        """
        self.app = app
        continue_saving = self.check_saved_json()
        if continue_saving == False:
            return
        else:
            # Create dictionary
            parameter_dictionary_class = Populate_dictionary_global(
                params, nmrdata, app
            )
            parameter_dictionary = parameter_dictionary_class.parameter_dictionary
            if continue_saving == True:
                # checking to see if processing parameters are saved
                # as well as conversion parameters
                saved_processing_params = self.check_saved_processing_params()
                if type(saved_processing_params) == Dict[str, Any]:
                    parameter_dictionary["processing"] = saved_processing_params

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
                self.app,
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

    def check_saved_processing_params(self) -> Union[Dict[str, Any], bool]:
        """
        Checking the processing.json file to see if there are processing
        parameters saved as well as conversion parameters.
        """

        try:
            with open("parameters.json", "r") as file:
                dictionary = json.load(file)
            if "processing" in dictionary.keys():
                return dictionary["processing"]

        except:
            # Unable to read the converter.json file effectively
            # Giving a warning to the user and returning
            message = (
                "Unable to read parameters.json. Saving without processing parameters."
            )
            dlg = wx.MessageDialog(self.app, message, "Warning", wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()

            return False


class Populate_dictionary_global:

    def __init__(self, params, nmrdata, app):
        """
        This class will find all the current parameters from the SpinConverter
        GUI which can then be saved.
        """

        self.params = params
        self.nmrdata = nmrdata
        self.app = app

        # Creating a dictionary to store the current parameters and adding
        # general information
        parameter_dictionary = self.create_dictionary()
        parameter_dictionary["conversion"] = self.add_general(
            parameter_dictionary["conversion"]
        )

        # Adding spectral parameters
        parameter_dictionary["conversion"] = self.add_complex_sizes(
            parameter_dictionary["conversion"]
        )
        parameter_dictionary["conversion"] = self.add_real_sizes(
            parameter_dictionary["conversion"]
        )
        parameter_dictionary["conversion"] = self.add_acqusition_mode(
            parameter_dictionary["conversion"]
        )
        parameter_dictionary["conversion"] = self.add_sweep_widths(
            parameter_dictionary["conversion"]
        )
        parameter_dictionary["conversion"] = self.add_nuclei_frequency(
            parameter_dictionary["conversion"]
        )
        parameter_dictionary["conversion"] = self.add_labels(
            parameter_dictionary["conversion"]
        )
        parameter_dictionary["conversion"] = self.carrier_frequency(
            parameter_dictionary["conversion"]
        )

        # Adding intensity scaling options
        parameter_dictionary["conversion"] = self.scaling_information(
            parameter_dictionary["conversion"]
        )

        if parameter_dictionary["conversion"]["general"]["spectrometer"] == "Bruker":
            # Adding digital filter information and other options
            parameter_dictionary["conversion"] = self.digital_filter_information(
                parameter_dictionary["conversion"]
            )
            parameter_dictionary["conversion"] = self.other_options_information(
                parameter_dictionary["conversion"]
            )

            # If there is more than one dimension, save the non-uniform sampling (NUS) information
            parameter_dictionary["conversion"] = self.add_nus_information(
                parameter_dictionary["conversion"]
            )

        self.parameter_dictionary = parameter_dictionary

    def create_dictionary(self) -> Dict:
        """
        Create a dictionary to hold the relevant parameters for the
        conversion (note: for bruker data extra information such as
        digital filter parameters also need to be saved).
        """
        if self.nmrdata.spectrometer == "Bruker":
            dictionary = {
                "conversion": {
                    "general": {"spectrometer": "Bruker"},
                    "spectral parameters": {},
                    "intensity scaling": {},
                    "digital filter parameters": {},
                    "other parameters": {},
                }
            }
        else:
            dictionary = {
                "conversion": {
                    "general": {"spectrometer": "Varian"},
                    "spectral parameters": {},
                    "intensity scaling": {},
                }
            }
        return dictionary

    def add_general(self, dictionary: Dict[str, str]) -> Dict:
        """
        Adding general information such as spectrometer type
        and temperature and number of scans
        """
        dictionary["general"]["temperature"] = str(self.params.temperature)
        dictionary["general"]["number of scans (NS)"] = str(self.params.NS)
        return dictionary

    def add_complex_sizes(self, dictionary: Dict[str, Any]) -> Dict[str, Any]:
        """
        Adding the complex dimension sizes to the dictionary
        """
        sizes = []
        for size_box in self.app.format.N_complex_boxes:
            sizes.append(str(int(size_box.GetValue())))
        dictionary["spectral parameters"]["sizes"] = {}
        dictionary["spectral parameters"]["sizes"]["complex"] = sizes
        return dictionary

    def add_real_sizes(self, dictionary: Dict[str, Any]) -> Dict[str, Any]:
        """
        Adding the real dimension sizes to the dictionary
        """
        sizes = []
        for size_box in self.app.format.N_real_boxes:
            sizes.append(str(int(size_box.GetValue())))
        dictionary["spectral parameters"]["sizes"]["real"] = sizes
        return dictionary

    def add_acqusition_mode(self, dictionary: Dict[str, Any]) -> Dict[str, Any]:
        """
        Adding the acqusition mode for each dimension
        (both name and selection index)
        """
        indirect_modes = []
        indirect_mode_indexes = []
        dictionary["spectral parameters"]["acqusition modes"] = {}
        dictionary["spectral parameters"]["acqusition modes"]["direct"] = {}
        dictionary["spectral parameters"]["acqusition modes"]["indirect"] = {}
        for i, box in enumerate(self.app.format.acqusition_combo_boxes):
            if i == 0:
                box_selection = box.GetSelection()
                mode = self.app.format.acquisition_mode_options_direct[box_selection]
                dictionary["spectral parameters"]["acqusition modes"]["direct"][
                    "mode"
                ] = mode
                dictionary["spectral parameters"]["acqusition modes"]["direct"][
                    "index"
                ] = box_selection

            else:
                box_selection = box.GetSelection()
                mode = self.app.format.acquisition_mode_options_indirect[box_selection]
                indirect_modes.append(mode)
                indirect_mode_indexes.append(box_selection)

        if indirect_modes != []:
            dictionary["spectral parameters"]["acqusition modes"]["indirect"][
                "mode"
            ] = indirect_modes
            dictionary["spectral parameters"]["acqusition modes"]["indirect"][
                "index"
            ] = indirect_mode_indexes

        return dictionary

    def add_sweep_widths(self, dictionary: Dict[str, Any]) -> Dict[str, Any]:
        """
        Adding the sweep width from each dimension to the dictionary
        """
        sweep_widths = []
        sweep_width_indexes = []
        for box in self.app.format.sweep_width_boxes:
            sweep_width_selection = box.GetSelection()
            sweep_width_value = box.GetValue()
            sweep_widths.append(str(sweep_width_value))
            sweep_width_indexes.append(str(sweep_width_selection))

        dictionary["spectral parameters"]["sweep widths"] = {}
        dictionary["spectral parameters"]["sweep widths"]["values"] = sweep_widths
        dictionary["spectral parameters"]["sweep widths"][
            "indexes"
        ] = sweep_width_indexes

        return dictionary

    def add_nuclei_frequency(self, dictionary: Dict[str, Any]) -> Dict[str, Any]:
        """
        Adding the nucleus frequency from each dimension to the dictionary
        """
        nuclei_frequencies = []
        for box in self.app.format.nuclei_frequency_boxes:
            nuclei_frequencies.append(str(box.GetValue()))

        dictionary["spectral parameters"]["nuclei frequencies"] = nuclei_frequencies

        return dictionary

    def add_labels(self, dictionary: Dict) -> Dict:
        """
        Adding the labels from each dimension to the dictionary
        """
        labels = []
        for box in self.app.format.nucleus_type_boxes:
            labels.append(str(box.GetValue()))
        dictionary["spectral parameters"]["labels"] = labels
        return dictionary

    def carrier_frequency(self, dictionary: Dict[str, Any]) -> Dict[str, Any]:
        """
        Adding the carrier frequencies from each dimension to the dictionary
        """
        carrier_frequencies = []
        carrier_frequency_combobox = []
        carrier_frequency_combobox_selection = []

        for box in self.app.format.carrier_frequency_boxes:
            carrier_frequencies.append(str(box.GetValue()))

        for box in self.app.format.carrier_combo_boxes:
            carrier_frequency_combobox.append(str(box.GetValue()))
            carrier_frequency_combobox_selection.append(str(box.GetSelection()))

        dictionary["spectral parameters"]["carrier frequencies"] = {}
        dictionary["spectral parameters"]["carrier frequencies"][
            "frequency"
        ] = carrier_frequencies
        dictionary["spectral parameters"]["carrier frequencies"][
            "combobox"
        ] = carrier_frequency_combobox
        dictionary["spectral parameters"]["carrier frequencies"][
            "combobox index"
        ] = carrier_frequency_combobox_selection

        return dictionary

    def scaling_information(self, dictionary: Dict[str, Any]) -> Dict[str, Any]:
        """
        Adding the spectrum scaling information to the dictionary
        """
        scale_by_ns_flag = self.app.shared_format.scaling_NS_checkbox.GetValue()
        if dictionary["general"]["spectrometer"] == "Bruker":
            scale_by_nc_flag = self.app.shared_format.scaling_NC.GetValue()
        scale_by_1000_flag = self.app.shared_format.scaling_by_number.GetValue()
        scaling_number = self.app.shared_format.scaling_number.GetValue()

        dictionary["intensity scaling"][
            "Scale by number of scans (NS)"
        ] = scale_by_ns_flag
        dictionary["intensity scaling"][
            "Scale by Bruker normalisation constant (NC)"
        ] = scale_by_nc_flag
        dictionary["intensity scaling"]["Scale by 1000"] = scale_by_1000_flag
        dictionary["intensity scaling"]["Scaling number"] = str(scaling_number)

        return dictionary

    def digital_filter_information(self, dictionary: Dict[str, Any]) -> Dict[str, Any]:
        """
        Adding the digital filter information to the dictionary
        from the Bruker data
        """
        decim = self.params.decim
        dspfvs = self.params.dspfvs
        grpdly = self.params.grpdly

        # Remove digital filter checkbox
        remove_digital_filter = self.app.format.digital_filter_checkbox.GetValue()

        # Remove digital filter before or after fourier transform (0 is before, 1 is after)
        try:
            remove_after_ft = self.app.format.digital_filter_radio_box.GetSelection()
        except:
            remove_after_ft = 1
        if remove_after_ft == 0:
            remove = "Before"
        else:
            remove = "After"
        dictionary["digital filter parameters"][
            "Remove Digital Filter"
        ] = remove_digital_filter
        dictionary["digital filter parameters"][
            "Remove Before/After Fourier Transform"
        ] = remove
        dictionary["digital filter parameters"]["Decimation Rate (decim)"] = str(decim)
        dictionary["digital filter parameters"]["DSP Firmware Version (dspfvs)"] = str(
            dspfvs
        )
        dictionary["digital filter parameters"]["Group Delay (grpdly)"] = str(grpdly)

        return dictionary

    def other_options_information(self, dictionary: Dict[str, Any]) -> Dict[str, Any]:
        """
        Adding extra information to the dictionary for Bruker data
        """
        remove_acqusition_padding_flag = self.params.remove_acquisition_padding = True
        bad_point_threshold = str(
            self.app.format.bad_point_threshold_textbox.GetValue()
        )
        dictionary["other parameters"][
            "remove acqusition padding"
        ] = remove_acqusition_padding_flag
        dictionary["other parameters"]["bad point threshold"] = bad_point_threshold
        return dictionary

    def add_nus_information(self, dictionary: Dict[str, Any]) -> Dict[str, Any]:
        """
        Adding the current NUS information to the dictionary
        """
        if len(dictionary["spectral parameters"]["sizes"]["complex"]) > 1:
            if self.app.shared_format.NUS_tickbox.GetValue() == True:
                NUS_checkbox = True
                dictionary["NUS information"] = {}
                nuslist = self.app.shared_format.nusfile_input.GetValue()
                nus_sample_count = (
                    self.app.shared_format.NUS_sample_count_box.GetValue()
                )
                nus_offset = self.app.shared_format.nus_offset_box.GetValue()
                reverse_nus_schedule = (
                    self.app.shared_format.reverse_NUS_tickbox.GetValue()
                )
                dictionary["NUS information"]["Checkbox"] = NUS_checkbox
                dictionary["NUS information"]["NUS sample count"] = str(
                    nus_sample_count
                )
                dictionary["NUS information"]["NUS offset"] = str(nus_offset)
                dictionary["NUS information"]["NUS file"] = str(nuslist)
                dictionary["NUS information"][
                    "Reverse NUS schedule"
                ] = reverse_nus_schedule

                return dictionary

        dictionary["NUS information"] = "N/A"
        return dictionary
