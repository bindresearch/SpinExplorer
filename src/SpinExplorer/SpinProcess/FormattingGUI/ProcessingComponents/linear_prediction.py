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

import wx


class LinearPrediction:
    def __init__(self, app, nmr_data, parent, info_buttons):
        """
        This class contains all the functions related to linear prediction
        for the direct dimension.
        """

        self.app = app
        self.nmr_data = nmr_data
        self.info_buttons = info_buttons

        self.set_initial_linear_prediction_variables()
        self.create_linear_prediction_sizer(parent)

    def create_linear_prediction_sizer(self, parent):
        """
        Create a box for the linear prediction options
        """
        self.linear_prediction_box = wx.StaticBox(parent, -1, "Linear Prediction")
        self.linear_prediction_sizer = wx.StaticBoxSizer(
            self.linear_prediction_box, wx.HORIZONTAL
        )
        self.linear_prediction_checkbox = wx.CheckBox(
            parent, -1, "Apply linear prediction"
        )
        self.linear_prediction_checkbox.SetValue(self.linear_prediction_checkbox_value)
        self.linear_prediction_checkbox.Bind(
            wx.EVT_CHECKBOX, self.on_linear_prediction_checkbox
        )
        self.linear_prediction_sizer.Add(
            self.linear_prediction_checkbox, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.linear_prediction_sizer.AddSpacer(10)
        # Have a combobox for linear prediction options
        self.linear_prediction_options_text = wx.StaticText(
            parent, -1, "Add Predicted Points:"
        )
        self.linear_prediction_sizer.Add(
            self.linear_prediction_options_text, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.linear_prediction_sizer.AddSpacer(5)
        self.linear_prediction_options = ["After FID", "Before FID"]
        self.linear_prediction_combobox = wx.ComboBox(
            parent, -1, choices=self.linear_prediction_options, style=wx.CB_READONLY
        )
        self.linear_prediction_combobox.Bind(
            wx.EVT_COMBOBOX, self.on_linear_prediction_combobox_options
        )
        self.linear_prediction_combobox.SetSelection(
            self.linear_prediction_options_selection
        )
        self.linear_prediction_sizer.Add(
            self.linear_prediction_combobox, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.linear_prediction_sizer.AddSpacer(10)
        # Have a combobox of predicted coefficient options
        self.linear_prediction_coefficients_text = wx.StaticText(
            parent, -1, "Predicted Coefficients:"
        )
        self.linear_prediction_sizer.Add(
            self.linear_prediction_coefficients_text, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.linear_prediction_sizer.AddSpacer(5)
        self.linear_prediction_coefficients_options = ["Forward", "Backward", "Both"]
        self.linear_prediction_coefficients_combobox = wx.ComboBox(
            parent,
            -1,
            choices=self.linear_prediction_coefficients_options,
            style=wx.CB_READONLY,
        )
        self.linear_prediction_coefficients_combobox.Bind(
            wx.EVT_COMBOBOX, self.on_linear_prediction_coefficients_combobox
        )
        self.linear_prediction_coefficients_combobox.SetSelection(
            self.linear_prediction_coefficients_selection
        )
        self.linear_prediction_sizer.Add(
            self.linear_prediction_coefficients_combobox, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.linear_prediction_sizer.AddSpacer(10)

        # Have a button showing information on linear prediction
        self.linear_prediction_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))
        self.linear_prediction_info.Bind(
            wx.EVT_BUTTON, self.info_buttons.on_linear_prediction_info
        )
        self.linear_prediction_sizer.Add(
            self.linear_prediction_info, 0, wx.ALIGN_CENTER_VERTICAL
        )
        parent.sizer_1.Add(self.linear_prediction_sizer)
        parent.sizer_1.AddSpacer(10)

    def set_initial_linear_prediction_variables(self):
        """
        Setting the initial linear prediction processing parameters
        in the graphical interface to default values
        """
        self.linear_prediction_checkbox_value = False
        self.linear_prediction_options_selection = 0
        self.linear_prediction_coefficients_selection = 0

    def on_linear_prediction_checkbox(self, event):
        """
        Change the checkbox value parameter if the checkbox
        is clicked.
        """
        if self.linear_prediction_checkbox.GetValue() == True:
            self.linear_prediction_checkbox_value = True
        else:
            self.linear_prediction_checkbox_value = False

    def on_linear_prediction_combobox_options(self, event):
        """
        Change the current linear prediction option when the
        combobox is changed.
        """
        self.linear_prediction_options_selection = (
            self.linear_prediction_combobox.GetSelection()
        )

    def on_linear_prediction_coefficients_combobox(self, event):
        """
        Change the linear prediction coefficient selection when
        the combobox is changed.
        """
        self.linear_prediction_coefficients_selection = (
            self.linear_prediction_coefficients_combobox.GetSelection()
        )
