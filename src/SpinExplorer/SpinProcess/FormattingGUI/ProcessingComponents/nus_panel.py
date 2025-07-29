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

    def create_linear_prediction_sizer_indirect(self, parent):
        """
        Creating a sizer for the linear prediction options with a radio box to
        toggle between:
        - No linear prediction
        - Linear prediction
        - Non-uniform sampling
        """

        self.linear_prediction_sizer_indirect_label = wx.StaticBox(
            parent, -1, "Linear Prediction/SMILE NUS Reconstruction"
        )
        self.linear_prediction_sizer_indirect = wx.StaticBoxSizer(
            self.linear_prediction_sizer_indirect_label, wx.HORIZONTAL
        )
        self.linear_prediction_sizer_indirect.AddSpacer(10)

        # Have a radiobox for None, Linear Prediction and SMILE NUS Reconstruction
        self.linear_prediction_radio_box_indirect = wx.RadioBox(
            parent,
            -1,
            "",
            choices=["None", "Linear Prediction", "SMILE NUS Reconstruction"],
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
                parent, -1, "Add Predicted Points:"
            )
            self.linear_prediction_sizer_indirect.Add(
                self.linear_prediction_options_text, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_indirect.AddSpacer(5)
            self.linear_prediction_options = ["After FID", "Before FID"]
            self.linear_prediction_combobox_indirect = wx.ComboBox(
                parent, -1, choices=self.linear_prediction_options, style=wx.CB_READONLY
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
                parent, -1, "Predicted Coefficients:"
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
                parent,
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
            self.smile_nus_file_text = wx.StaticText(parent, -1, "NUS File:")
            self.linear_prediction_sizer_indirect.Add(
                self.smile_nus_file_text, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_indirect.AddSpacer(5)

            self.smile_nus_file_textcontrol_indirect = wx.TextCtrl(
                parent, -1, self.nuslist_name_indirect, size=(100, 20)
            )
            self.smile_nus_file_textcontrol_indirect.Bind(
                wx.EVT_TEXT, self.on_smile_nus_file_textcontrol_indirect
            )

            self.linear_prediction_sizer_indirect.Add(
                self.smile_nus_file_textcontrol_indirect, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_indirect.AddSpacer(10)

            # # Zero order phase correction
            # self.smile_nus_p0_text = wx.StaticText(parent, -1, 'Zero Order Phase Correction (p0):')
            # self.linear_prediction_sizer_indirect.Add(self.smile_nus_p0_text, 0, wx.ALIGN_CENTER_VERTICAL)
            # self.linear_prediction_sizer_indirect.AddSpacer(5)
            # self.smile_nus_p0_textcontrol_indirect = wx.TextCtrl(parent, -1, str(self.p0_total_indirect), size=(50, 20))
            # self.smile_nus_p0_textcontrol_indirect.Bind(wx.EVT_TEXT, self.on_smile_nus_p0_textcontrol_indirect)
            # self.linear_prediction_sizer_indirect.Add(self.smile_nus_p0_textcontrol_indirect, 0, wx.ALIGN_CENTER_VERTICAL)
            # self.linear_prediction_sizer_indirect.AddSpacer(10)

            # # First order phase correction
            # self.smile_nus_p1_text = wx.StaticText(parent, -1, 'First Order Phase Correction (p1):')
            # self.linear_prediction_sizer_indirect.Add(self.smile_nus_p1_text, 0, wx.ALIGN_CENTER_VERTICAL)
            # self.linear_prediction_sizer_indirect.AddSpacer(5)
            # self.smile_nus_p1_textcontrol_indirect = wx.TextCtrl(parent, -1, str(self.p1_total_indirect), size=(50, 20))
            # self.smile_nus_p1_textcontrol_indirect.Bind(wx.EVT_TEXT, self.on_smile_nus_p1_textcontrol_indirect)
            # self.linear_prediction_sizer_indirect.Add(self.smile_nus_p1_textcontrol_indirect, 0, wx.ALIGN_CENTER_VERTICAL)
            # self.linear_prediction_sizer_indirect.AddSpacer(10)

            # Number of points to add to the data
            self.smile_nus_extension_text = wx.StaticText(parent, -1, "Data extension:")
            self.linear_prediction_sizer_indirect.Add(
                self.smile_nus_extension_text, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_indirect.AddSpacer(5)
            self.smile_nus_extension_textcontrol_indirect = wx.TextCtrl(
                parent,
                -1,
                str(self.smile_data_extension_number_indirect),
                size=(50, 20),
            )
            self.smile_nus_extension_textcontrol_indirect.Bind(
                wx.EVT_TEXT, self.on_smile_nus_extension_textcontrol_indirect
            )
            self.linear_prediction_sizer_indirect.Add(
                self.smile_nus_extension_textcontrol_indirect,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            self.linear_prediction_sizer_indirect.AddSpacer(10)

            # Number of CPU's
            self.smile_nus_cpu_text = wx.StaticText(parent, -1, "Number of CPU's:")
            self.linear_prediction_sizer_indirect.Add(
                self.smile_nus_cpu_text, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_indirect.AddSpacer(5)
            self.smile_nus_cpu_textcontrol_indirect = wx.TextCtrl(
                parent, -1, str(self.number_of_nus_CPU_indirect), size=(30, 20)
            )
            self.smile_nus_cpu_textcontrol_indirect.Bind(
                wx.EVT_TEXT, self.on_smile_nus_cpu_textcontrol_indirect
            )
            self.linear_prediction_sizer_indirect.Add(
                self.smile_nus_cpu_textcontrol_indirect, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_indirect.AddSpacer(10)

            # Number of iterations
            self.smile_nus_iterations_text = wx.StaticText(
                parent, -1, "Number of Iterations:"
            )
            self.linear_prediction_sizer_indirect.Add(
                self.smile_nus_iterations_text, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_indirect.AddSpacer(5)
            self.smile_nus_iterations_textcontrol_indirect = wx.TextCtrl(
                parent, -1, str(self.nus_iterations_indirect), size=(30, 20)
            )
            self.smile_nus_iterations_textcontrol_indirect.Bind(
                wx.EVT_TEXT, self.on_smile_nus_iterations_textcontrol_indirect
            )
            self.linear_prediction_sizer_indirect.Add(
                self.smile_nus_iterations_textcontrol_indirect,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )

        # Have a button showing information on linear prediction
        self.linear_prediction_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))
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

    # def on_smile_nus_p0_textcontrol_indirect(self, event):
    #     self.p0_total_indirect = self.smile_nus_p0_textcontrol_indirect.GetValue()
    #     self.phasing_from_smile = True
    #     # Update the phasing values in the phasing section too
    #     self.phase_correction_p0_textcontrol_indirect.SetValue(
    #         str(self.p0_total_indirect)
    #     )
    #     self.phasing_from_smile = False

    # def on_smile_nus_p1_textcontrol_indirect(self, event):
    #     self.p1_total_indirect = self.smile_nus_p1_textcontrol_indirect.GetValue()
    #     self.phasing_from_smile = True
    #     # Update the phasing values in the phasing section too
    #     self.phase_correction_p1_textcontrol_indirect.SetValue(
    #         str(self.p1_total_indirect)
    #     )
    #     self.phasing_from_smile = False

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

    def on_linear_prediction_radio_box_indirect(self, event):
        """
        Get the selection from the radio box and update the
        linear prediction options.
        """
        self.linear_prediction_radio_box_indirect_selection = (
            self.linear_prediction_radio_box_indirect.GetSelection()
        )

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
