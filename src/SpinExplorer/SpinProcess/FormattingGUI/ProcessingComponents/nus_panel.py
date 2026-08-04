#!/usr/bin/env python3

"""MIT License

Copyright (c) 2025 James Eaton, Andrew Baldwin (University of Oxford)
              2025, Bind Research

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
import os
import json


class NonUniformSampling:
    def __init__(self, app, nmr_data, parent, info_buttons):
        """
        This class produces all the graphical interface functionality
        relevent to non-uniform sampling (NUS) or linear prediction
        in the indirect dimensions.
        """

        self.app = app
        self.nmr_data = nmr_data
        self.info_buttons = info_buttons
        self.apodization_class = (
            0  # Placeholder for the class which will be added later
        )
        self.parent = parent

        self.set_initial_nus_variables()
        self.create_linear_prediction_sizer_indirect(parent)

    def set_initial_nus_variables(self):
        """
        Setting initial variables relevent to NUS and linear prediction
        in the indirect dimensions.
        """
        self.linear_prediction_radio_box_indirect_selection = 0
        self.linear_prediction_indirect_checkbox_value = False
        self.linear_prediction_indirect_options_selection = 0
        self.linear_prediction_indirect_coefficients_selection = 0
        self.linear_prediction_selection = 0

        # Check to see if the nuslist file exists in the current directory using os.path.isfile('nuslist')
        if os.path.isfile("nuslist"):
            self.nuslist_name_indirect = "nuslist"
        else:
            self.nuslist_name_indirect = ""

        self.number_of_nus_CPU_indirect = 1
        self.nus_iterations_indirect = 50
        self.smile_data_extension_number_indirect = (
            0  # int(self.nmr_data.number_of_points[1]*1.5)
        )
        self.ist_data_extension_number_indirect = 0
        self.ist_linear_prediction_only_flag = self.find_ist_linear_prediction_only_flag()
        self.ist_nus_iterations_indirect = 2000

    def find_ist_linear_prediction_only_flag(self):
        """
        Read through the parameters file and see if NUS reshuffling was performed
        during the conversion process. If it was, then return False, otherwise
        return True.
        """

        try:
            with open("parameters.json", "r") as file:
                parameter_dictionary = json.load(file)["conversion"]
                nus_dic = parameter_dictionary["NUS information"]
                if nus_dic == "N/A":
                    return True
                else:
                    return False
    
        except:
            return False

    def create_linear_prediction_sizer_indirect(self, parent):
        """
        Creating a sizer for the linear prediction options with a radio box to
        toggle between:
        - No linear prediction
        - Linear prediction
        - Non-uniform sampling
        """

        self.linear_prediction_sizer_indirect_label = wx.StaticBox(
            parent, -1, "Linear Prediction / NUS Reconstruction"
        )
        self.linear_prediction_sizer_indirect = wx.StaticBoxSizer(
            self.linear_prediction_sizer_indirect_label, wx.HORIZONTAL
        )
        self.linear_prediction_sizer_indirect.AddSpacer(10)

        # Have a radiobox for None, Linear Prediction and SMILE NUS Reconstruction
        self.linear_prediction_radio_box_indirect = wx.RadioBox(
            self.linear_prediction_sizer_indirect_label,
            -1,
            "",
            choices=["None", "Linear Prediction", "SMILE NUS Reconstruction", "SpinExplorer IST NUS Reconstruction"],
            style=wx.RA_SPECIFY_ROWS,
        )
        self.linear_prediction_radio_box_indirect.Bind(
            wx.EVT_RADIOBOX, self.on_linear_prediction_radio_box_indirect
        )
        self.linear_prediction_radio_box_indirect.SetSelection(
            self.linear_prediction_radio_box_indirect_selection
        )

        self.linear_prediction_sizer_indirect.Add(
            self.linear_prediction_radio_box_indirect, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.linear_prediction_sizer_indirect.AddSpacer(10)

        if self.linear_prediction_radio_box_indirect.GetSelection() == 1:
            # Have a combobox for linear prediction options
            self.linear_prediction_options_text = wx.StaticText(
                self.linear_prediction_sizer_indirect_label, -1, "Add Predicted Points:"
            )
            self.linear_prediction_sizer_indirect.Add(
                self.linear_prediction_options_text, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_indirect.AddSpacer(5)
            self.linear_prediction_options = ["After FID", "Before FID"]
            self.linear_prediction_combobox_indirect = wx.ComboBox(
                self.linear_prediction_sizer_indirect_label, -1, choices=self.linear_prediction_options, style=wx.CB_READONLY
            )
            self.linear_prediction_combobox_indirect.SetSelection(
                self.linear_prediction_indirect_options_selection
            )
            self.linear_prediction_combobox_indirect.Bind(
                wx.EVT_COMBOBOX, self.on_linear_prediction_combobox_indirect
            )
            self.linear_prediction_sizer_indirect.Add(
                self.linear_prediction_combobox_indirect, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_indirect.AddSpacer(10)
            # Have a combobox of predicted coefficient options
            self.linear_prediction_coefficients_text = wx.StaticText(
                self.linear_prediction_sizer_indirect_label, -1, "Predicted Coefficients:"
            )
            self.linear_prediction_sizer_indirect.Add(
                self.linear_prediction_coefficients_text, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_indirect.AddSpacer(5)
            self.linear_prediction_coefficients_options = [
                "Forward",
                "Backward",
                "Both",
            ]
            self.linear_prediction_coefficients_combobox_indirect = wx.ComboBox(
                self.linear_prediction_sizer_indirect_label,
                -1,
                choices=self.linear_prediction_coefficients_options,
                style=wx.CB_READONLY,
            )
            self.linear_prediction_coefficients_combobox_indirect.SetSelection(
                self.linear_prediction_indirect_coefficients_selection
            )
            self.linear_prediction_coefficients_combobox_indirect.Bind(
                wx.EVT_COMBOBOX,
                self.on_linear_prediction_combobox_coefficients_indirect,
            )
            self.linear_prediction_sizer_indirect.Add(
                self.linear_prediction_coefficients_combobox_indirect,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            self.linear_prediction_sizer_indirect.AddSpacer(10)
        elif self.linear_prediction_radio_box_indirect.GetSelection() == 2:
            # Have a set of options for SMILE NUS processing

            # NUS file
            self.smile_nus_file_text = wx.StaticText(self.linear_prediction_sizer_indirect_label, -1, "NUS File:")
            self.linear_prediction_sizer_indirect.Add(
                self.smile_nus_file_text, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_indirect.AddSpacer(5)

            self.smile_nus_file_textcontrol_indirect = wx.TextCtrl(
                self.linear_prediction_sizer_indirect_label, -1, self.nuslist_name_indirect, size=(100, 20)
            )
            self.smile_nus_file_textcontrol_indirect.Bind(
                wx.EVT_TEXT, self.on_smile_nus_file_textcontrol_indirect
            )

            self.linear_prediction_sizer_indirect.Add(
                self.smile_nus_file_textcontrol_indirect, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_indirect.AddSpacer(10)


            # Number of points to add to the data
            self.smile_nus_extension_text = wx.StaticText(self.linear_prediction_sizer_indirect_label, -1, "Data extension:")
            self.linear_prediction_sizer_indirect.Add(
                self.smile_nus_extension_text, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_indirect.AddSpacer(5)
            self.smile_nus_extension_textcontrol_indirect = wx.TextCtrl(
                self.linear_prediction_sizer_indirect_label,
                -1,
                str(self.smile_data_extension_number_indirect),
                size=(50, 20), style=wx.TE_PROCESS_ENTER
            )
            self.smile_nus_extension_textcontrol_indirect.Bind(
                wx.EVT_TEXT_ENTER, self.on_smile_nus_extension_textcontrol_indirect
            )
            self.linear_prediction_sizer_indirect.Add(
                self.smile_nus_extension_textcontrol_indirect,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            self.linear_prediction_sizer_indirect.AddSpacer(10)

            # Number of CPU's
            self.smile_nus_cpu_text = wx.StaticText(self.linear_prediction_sizer_indirect_label, -1, "Number of CPU's:")
            self.linear_prediction_sizer_indirect.Add(
                self.smile_nus_cpu_text, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_indirect.AddSpacer(5)
            self.smile_nus_cpu_textcontrol_indirect = wx.TextCtrl(
                self.linear_prediction_sizer_indirect_label, -1, str(self.number_of_nus_CPU_indirect), size=(30, 20), style=wx.TE_PROCESS_ENTER
            )
            self.smile_nus_cpu_textcontrol_indirect.Bind(
                wx.EVT_TEXT_ENTER, self.on_smile_nus_cpu_textcontrol_indirect
            )
            self.linear_prediction_sizer_indirect.Add(
                self.smile_nus_cpu_textcontrol_indirect, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_indirect.AddSpacer(10)

            # Number of iterations
            self.smile_nus_iterations_text = wx.StaticText(
                self.linear_prediction_sizer_indirect_label, -1, "Number of Iterations:"
            )
            self.linear_prediction_sizer_indirect.Add(
                self.smile_nus_iterations_text, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_indirect.AddSpacer(5)
            self.smile_nus_iterations_textcontrol_indirect = wx.TextCtrl(
                self.linear_prediction_sizer_indirect_label, -1, str(self.nus_iterations_indirect), size=(50, 20), style=wx.TE_PROCESS_ENTER
            )
            self.smile_nus_iterations_textcontrol_indirect.Bind(
                wx.EVT_TEXT_ENTER, self.on_smile_nus_iterations_textcontrol_indirect
            )
            self.linear_prediction_sizer_indirect.Add(
                self.smile_nus_iterations_textcontrol_indirect,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )

        elif self.linear_prediction_radio_box_indirect.GetSelection() == 3:
            # Have a set of options for IST NUS processing


            # NUS file
            self.ist_nus_file_text = wx.StaticText(self.linear_prediction_sizer_indirect_label, -1, "NUS File:")
            self.linear_prediction_sizer_indirect.Add(
                self.ist_nus_file_text, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_indirect.AddSpacer(5)

            self.ist_nus_file_textcontrol_indirect = wx.TextCtrl(
                self.linear_prediction_sizer_indirect_label, -1, self.nuslist_name_indirect, size=(100, 20), style=wx.TE_PROCESS_ENTER
            )
            self.ist_nus_file_textcontrol_indirect.Bind(
                wx.EVT_TEXT_ENTER, self.on_ist_nus_file_textcontrol_indirect
            )

            self.linear_prediction_sizer_indirect.Add(
                self.ist_nus_file_textcontrol_indirect, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_indirect.AddSpacer(10)


            # Number of points to add to the data
            self.ist_nus_extension_text = wx.StaticText(self.linear_prediction_sizer_indirect_label, -1, "Data extension:")
            self.linear_prediction_sizer_indirect.Add(
                self.ist_nus_extension_text, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_indirect.AddSpacer(5)
            self.ist_nus_extension_textcontrol_indirect = wx.TextCtrl(
                self.linear_prediction_sizer_indirect_label,
                -1,
                str(self.ist_data_extension_number_indirect),
                size=(50, 20), style=wx.TE_PROCESS_ENTER
            )
            self.ist_nus_extension_textcontrol_indirect.Bind(
                wx.EVT_TEXT_ENTER, self.on_ist_nus_extension_textcontrol_indirect
            )
            self.linear_prediction_sizer_indirect.Add(
                self.ist_nus_extension_textcontrol_indirect,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )

            self.linear_prediction_sizer_indirect.AddSpacer(10)

            # Number of iterations
            self.ist_nus_iterations_text = wx.StaticText(
                self.linear_prediction_sizer_indirect_label, -1, "Number of Iterations:"
            )
            self.linear_prediction_sizer_indirect.Add(
                self.ist_nus_iterations_text, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_indirect.AddSpacer(5)


            self.ist_nus_iterations_textcontrol_indirect = wx.TextCtrl(
                self.linear_prediction_sizer_indirect_label, -1, str(self.ist_nus_iterations_indirect), size=(50, 20), style=wx.TE_PROCESS_ENTER
            )
            self.ist_nus_iterations_textcontrol_indirect.Bind(
                wx.EVT_TEXT_ENTER, self.on_ist_nus_iterations_textcontrol_indirect
            )
            self.linear_prediction_sizer_indirect.Add(
                self.ist_nus_iterations_textcontrol_indirect,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )


            self.linear_prediction_sizer_indirect.AddSpacer(10)

            # Checkbox to determine if NUS reconstruction is to be applied or if IST is to be used
            # for only linear prediction
            self.ist_linear_prediction_only = wx.CheckBox(self.linear_prediction_sizer_indirect_label, -1, label='Data extension only')
            self.ist_linear_prediction_only.SetValue(self.ist_linear_prediction_only_flag)
            self.linear_prediction_sizer_indirect.Add(self.ist_linear_prediction_only, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_indirect.AddSpacer(5)
            self.ist_linear_prediction_only.Bind(wx.EVT_CHECKBOX, self.OnIST_LP_Only)

        # Have a button showing information on linear prediction
        self.linear_prediction_info = wx.Button(self.linear_prediction_sizer_indirect_label, -1, "\u24d8", size=(25, 32))
        self.linear_prediction_info.Bind(
            wx.EVT_BUTTON, self.info_buttons.on_linear_prediction_info_indirect
        )
        self.linear_prediction_sizer_indirect.AddSpacer(10)
        self.linear_prediction_sizer_indirect.Add(
            self.linear_prediction_info, 0, wx.ALIGN_CENTER_VERTICAL
        )
        parent.sizer_1.Add(self.linear_prediction_sizer_indirect)
        parent.sizer_1.AddSpacer(10)

    
    def on_linear_prediction_combobox_indirect(self, event):
        """
        When the linear prediction combobox is changed, update the
        stored value.
        """
        # Get the selection from the combobox and update the linear prediction options
        self.linear_prediction_indirect_options_selection = (
            self.linear_prediction_combobox_indirect.GetSelection()
        )

    def on_linear_prediction_combobox_coefficients_indirect(self, event):
        """
        Get the selection from the combobox and update the linear prediction options
        """
        self.linear_prediction_indirect_coefficients_selection = (
            self.linear_prediction_coefficients_combobox_indirect.GetSelection()
        )

    def on_smile_nus_file_textcontrol_indirect(self, event):
        """
        Get the value from the textcontrol
        """
        self.nuslist_name_indirect = self.smile_nus_file_textcontrol_indirect.GetValue()

        if(self.parent.parent.nmr_data.dim == 3 and self.parent.parent.nmr_data.pseudo_axis == False):
            # If IST linear prediction only is selected and there are more than 1 complex indirect dimensions, change the selection to the same for both indirect dimensions
            if(self.parent.parent.tabDim2!=self):
                self.parent.parent.tabDim2.linear_prediction.nuslist_name_indirect = self.nuslist_name_indirect
                self.parent.parent.tabDim2.smile_nus_file_textcontrol_indirect.SetValue(self.nuslist_name_indirect)
            if(self.parent.parent.tabDim3!=self):
                self.parent.parent.tabDim3.linear_prediction.nuslist_name_indirect = self.nuslist_name_indirect
                self.parent.parent.tabDim3.linear_prediction.smile_nus_file_textcontrol_indirect.SetValue(self.nuslist_name_indirect)


    def OnIST_LP_Only(self, event):
        """
        Update the parameter from the checkbox current value
        """
        self.ist_linear_prediction_only_flag = self.ist_linear_prediction_only.GetValue()

        if(self.parent.parent.nmr_data.dim == 3 and self.parent.parent.nmr_data.pseudo_axis == False):
            # If IST linear prediction only is selected and there are more than 1 complex indirect dimensions, change the selection to the same for both indirect dimensions
            if(self.parent.parent.tabDim2!=self):
                self.parent.parent.tabDim2.linear_prediction.ist_linear_prediction_only_flag = self.ist_linear_prediction_only_flag
                self.parent.parent.tabDim2.linear_prediction.ist_linear_prediction_only.SetValue(self.ist_linear_prediction_only_flag)
            if(self.parent.parent.tabDim3!=self):
                self.parent.parent.tabDim3.linear_prediction.ist_linear_prediction_only_flag = self.ist_linear_prediction_only_flag
                self.parent.parent.tabDim3.linear_prediction.ist_linear_prediction_only.SetValue(self.ist_linear_prediction_only_flag)


    
    def on_ist_nus_file_textcontrol_indirect(self, event):
        """
        Get the value from the textcontrol
        """
        self.nuslist_name_indirect = self.ist_nus_file_textcontrol_indirect.GetValue()

        if(self.parent.parent.nmr_data.dim == 3 and self.parent.parent.nmr_data.pseudo_axis == False):
            # If IST linear prediction only is selected and there are more than 1 complex indirect dimensions, change the selection to the same for both indirect dimensions
            if(self.parent.parent.tabDim2!=self):
                self.parent.parent.tabDim2.linear_prediction.nuslist_name_indirect = self.nuslist_name_indirect
                self.parent.parent.tabDim2.linear_prediction.ist_nus_file_textcontrol_indirect.SetValue(self.nuslist_name_indirect)
            if(self.parent.parent.tabDim3!=self):
                self.parent.parent.tabDim3.linear_prediction.nuslist_name_indirect = self.nuslist_name_indirect
                self.parent.parent.tabDim3.linear_prediction.ist_nus_file_textcontrol_indirect.SetValue(self.nuslist_name_indirect)


    def on_smile_nus_extension_textcontrol_indirect(self, event):
        """
        When changing the nus extension number, this function checks
        the parameter validity (must be an integer) and updates
        the stored value.
        """
        if self.smile_nus_extension_textcontrol_indirect.GetValue() != "":
            try:
                self.smile_data_extension_number_indirect = int(
                    self.smile_nus_extension_textcontrol_indirect.GetValue()
                )
            except:
                msg = wx.MessageDialog(
                    self.parent,
                    "The value entered for NUS data extension not a valid integer",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.smile_nus_extension_textcontrol_indirect.SetValue(
                    str(self.smile_data_extension_number_indirect)
                )
                return
        else:
            self.smile_data_extension_number_indirect = (
                self.smile_nus_extension_textcontrol_indirect.GetValue()
            )

    def on_ist_nus_extension_textcontrol_indirect(self, event):
        """
        When changing the nus extension number, this function checks
        the parameter validity (must be an integer) and updates
        the stored value.
        """
        if self.ist_nus_extension_textcontrol_indirect.GetValue() != "":
            try:
                self.ist_data_extension_number_indirect = int(
                    self.ist_nus_extension_textcontrol_indirect.GetValue()
                )
            except:
                msg = wx.MessageDialog(
                    self.parent,
                    "The value entered for NUS data extension not a valid integer",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.ist_nus_extension_textcontrol_indirect.SetValue(
                    str(self.ist_data_extension_number_indirect)
                )
                return
        else:
            self.ist_data_extension_number_indirect = (
                self.ist_nus_extension_textcontrol_indirect.GetValue()
            )

    def on_smile_nus_cpu_textcontrol_indirect(self, event):
        """
        When changing the nus CPU number, this function checks
        the parameter validity (must be an integer) and updates
        the stored value.
        """
        if self.smile_nus_cpu_textcontrol_indirect.GetValue() != "":
            try:
                self.number_of_nus_CPU_indirect = int(
                    self.smile_nus_cpu_textcontrol_indirect.GetValue()
                )

                if(self.parent.parent.nmr_data.dim == 3 and self.parent.parent.nmr_data.pseudo_axis == False):
                    # If SMILE NUS is selected and there are more than 1 complex indirect dimensions, change the number of CPU's to the same for both indirect dimensions
                    if(self.parent.parent.tabDim2!=self):
                        self.parent.parent.tabDim2.linear_prediction.number_of_nus_CPU_indirect = self.number_of_nus_CPU_indirect
                        self.parent.parent.tabDim2.linear_prediction.smile_nus_cpu_textcontrol_indirect.SetValue(str(self.number_of_nus_CPU_indirect))
                    if(self.parent.parent.tabDim3!=self):
                        self.parent.parent.tabDim3.linear_prediction.number_of_nus_CPU_indirect = self.number_of_nus_CPU_indirect
                        self.parent.parent.tabDim3.linear_prediction.smile_nus_cpu_textcontrol_indirect.SetValue(str(self.number_of_nus_CPU_indirect))
            except:
                msg = wx.MessageDialog(
                    self.parent,
                    "The value entered for number of CPU's is not a valid integer",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.smile_nus_cpu_textcontrol_indirect.SetValue(
                    str(self.number_of_nus_CPU_indirect)
                )
                return

    def on_smile_nus_iterations_textcontrol_indirect(self, event):
        """
        When changing the nus iteration number, this function checks
        the parameter validity (must be an integer) and updates
        the stored value.
        """
        if self.smile_nus_iterations_textcontrol_indirect.GetValue() != "":
            try:
                self.nus_iterations_indirect = int(
                    self.smile_nus_iterations_textcontrol_indirect.GetValue()
                )

                if(self.parent.parent.nmr_data.dim == 3 and self.parent.parent.nmr_data.pseudo_axis == False):
                    # If SMILE NUS is selected and there are more than 1 complex indirect dimensions, change the number of iterations to the same for both indirect dimensions
                    if(self.parent.parent.tabDim2!=self):
                        self.parent.parent.tabDim2.linear_prediction.nus_iterations_indirect = self.nus_iterations_indirect
                        self.parent.parent.tabDim2.linear_prediction.smile_nus_iterations_textcontrol_indirect.SetValue(str(self.nus_iterations_indirect))
                    if(self.parent.parent.tabDim3!=self):
                        self.parent.parent.tabDim3.linear_prediction.nus_iterations_indirect = self.nus_iterations_indirect
                        self.parent.parent.tabDim3.linear_prediction.smile_nus_iterations_textcontrol_indirect.SetValue(str(self.nus_iterations_indirect))
            except:
                msg = wx.MessageDialog(
                    self.parent,
                    "The value entered for number of iterations is not a valid integer",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.smile_nus_iterations_textcontrol_indirect.SetValue(
                    str(self.nus_iterations_indirect)
                )
                return
            


    def on_ist_nus_iterations_textcontrol_indirect(self, event):
        """
        When changing the nus iteration number, this function checks
        the parameter validity (must be an integer) and updates
        the stored value.
        """
        if self.ist_nus_iterations_textcontrol_indirect.GetValue() != "":
            try:
                self.ist_nus_iterations_indirect = int(
                    self.ist_nus_iterations_textcontrol_indirect.GetValue()
                )

                if(self.parent.parent.nmr_data.dim == 3 and self.parent.parent.nmr_data.pseudo_axis == False):
                    # If SMILE NUS is selected and there are more than 1 complex indirect dimensions, change the number of iterations to the same for both indirect dimensions
                    if(self.parent.parent.tabDim2!=self):
                        self.parent.parent.tabDim2.linear_prediction.ist_nus_iterations_indirect = self.ist_nus_iterations_indirect
                        self.parent.parent.tabDim2.linear_prediction.ist_nus_iterations_textcontrol_indirect.SetValue(str(self.ist_nus_iterations_indirect))
                    if(self.parent.parent.tabDim3!=self):
                        self.parent.parent.tabDim3.linear_prediction.nus_iterations_indirect = self.ist_nus_iterations_indirect
                        self.parent.parent.tabDim3.linear_prediction.ist_nus_iterations_textcontrol_indirect.SetValue(str(self.ist_nus_iterations_indirect))
            except:
                msg = wx.MessageDialog(
                    self.parent,
                    "The value entered for number of iterations is not a valid integer",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.ist_nus_iterations_textcontrol_indirect.SetValue(
                    str(self.nus_iterations_indirect)
                )
                return

    def on_linear_prediction_radio_box_indirect(self, event, match_dimensions=False):
        """
        Get the selection from the radio box and update the
        linear prediction options.
        """
        self.linear_prediction_radio_box_indirect_selection = (
            self.linear_prediction_radio_box_indirect.GetSelection()
        )

        if(match_dimensions == False):
            if(self.parent.parent.nmr_data.dim == 3 and self.parent.parent.nmr_data.pseudo_axis == False):
                if(self.linear_prediction_radio_box_indirect_selection == 2 or self.linear_prediction_radio_box_indirect_selection == 3):
                    # If SMILE or IST is selected and there are more than 1 complex indirect dimensions, change the selection to te same for both indirect dimensions
                    if(self.parent.parent.tabDim2!=self):
                        self.parent.parent.tabDim2.linear_prediction.linear_prediction_radio_box_indirect_selection = self.linear_prediction_radio_box_indirect_selection
                        self.parent.parent.tabDim2.linear_prediction.linear_prediction_radio_box_indirect.SetSelection(self.linear_prediction_radio_box_indirect_selection)
                        self.parent.parent.tabDim2.linear_prediction.on_linear_prediction_radio_box_indirect(wx.EVT_RADIOBOX, match_dimensions=True)
                    if(self.parent.parent.tabDim3!=self):
                        self.parent.parent.tabDim3.linear_prediction.linear_prediction_radio_box_indirect_selection = self.linear_prediction_radio_box_indirect_selection
                        self.parent.parent.tabDim3.linear_prediction.linear_prediction_radio_box_indirect.SetSelection(self.linear_prediction_radio_box_indirect_selection)
                        self.parent.parent.tabDim3.linear_prediction.on_linear_prediction_radio_box_indirect(wx.EVT_RADIOBOX, match_dimensions=True)

        # Remove all the old sizers and replot

        self.apodization_class.on_apodization_combobox(wx.EVT_COMBOBOX)
        # self.apodization_sizer.Clear(delete_windows=True)

        # Remove the linear prediction sizers
        self.linear_prediction_sizer_indirect.Clear(delete_windows=True)

        self.parent.sizer_1.Clear(delete_windows=True)

        self.parent.refresh_menu_bar()
        self.parent.Refresh()
        self.parent.Update()
        self.parent.Layout()

        

