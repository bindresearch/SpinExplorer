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


class ZeroFilling:

    def __init__(self, app, nmr_data, parent, info_buttons, other_classes, dimension):
        """
        This class contains all the functions related to zero filling
        """

        self.app = app
        self.nmr_data = nmr_data
        self.info_buttons = info_buttons
        self.dimension_index = dimension
        if self.dimension_index == 0:
            self.solvent_suppression_class = other_classes[0]
            self.linear_prediction_class = other_classes[1]
            self.apodization_class = other_classes[2]
        else:
            self.linear_prediction_class = other_classes[0]
            self.apodization_class = other_classes[1]

        self.set_initial_zero_filling_variables()
        self.create_zero_filling_sizer(parent)

    def set_initial_zero_filling_variables(self):
        """
        Setting initial variables relevent to zero filling
        """
        self.zero_filling_checkbox_value = True
        self.zero_filling_combobox_selection = 0
        self.zero_filling_combobox_selection_old = 0
        self.zero_filling_value_doubling_times = 1
        self.zero_filling_value_zeros_to_add = 0
        self.zero_filling_value_final_data_size = (
            self.nmr_data.number_of_points[self.dimension_index] * 2
        )
        self.zero_filling_round_checkbox_value = True

    def create_zero_filling_sizer(self, parent):
        """
        Create a box for zero filling options
        """
        self.zero_filling_box = wx.StaticBox(parent, -1, "Zero Filling")
        self.zero_filling_sizer = wx.StaticBoxSizer(
            self.zero_filling_box, wx.HORIZONTAL
        )
        self.zero_filling_checkbox = wx.CheckBox(parent, -1, "Apply zero filling")
        self.zero_filling_checkbox.SetValue(self.zero_filling_checkbox_value)
        self.zero_filling_checkbox.Bind(wx.EVT_CHECKBOX, self.on_zero_filling_checkbox)
        self.zero_filling_sizer.Add(
            self.zero_filling_checkbox, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.zero_filling_sizer.AddSpacer(10)
        # Have a combobox for zero filling options
        self.zf_options_label = wx.StaticText(parent, -1, "Options:")
        self.zero_filling_sizer.Add(self.zf_options_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.zero_filling_sizer.AddSpacer(5)
        self.zero_filling_options = [
            "Doubling spectrum size",
            "Adding Zeros",
            "Final data size",
        ]
        self.zero_filling_combobox = wx.ComboBox(
            parent, -1, choices=self.zero_filling_options, style=wx.CB_READONLY
        )
        self.zero_filling_combobox.Bind(wx.EVT_COMBOBOX, self.on_zero_filling_combobox)
        self.zero_filling_combobox.SetSelection(self.zero_filling_combobox_selection)
        self.zero_filling_sizer.Add(
            self.zero_filling_combobox, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.zero_filling_sizer.AddSpacer(10)
        if self.zero_filling_combobox_selection == 0:
            # Have a textcontrol for the doubling number/number of zeros/final data size
            self.zf_value_label = wx.StaticText(parent, -1, "Doubling number:")
            self.zero_filling_sizer.Add(
                self.zf_value_label, 0, wx.ALIGN_CENTER_VERTICAL
            )

            self.zero_filling_textcontrol = wx.TextCtrl(
                parent, -1, str(self.zero_filling_value_doubling_times), size=(40, 20)
            )
            self.zero_filling_textcontrol.Bind(
                wx.EVT_TEXT, self.on_zero_filling_textcontrol_doubling_times
            )
            self.zero_filling_sizer.AddSpacer(5)
            self.zero_filling_sizer.Add(
                self.zero_filling_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.zero_filling_sizer.AddSpacer(20)
        elif self.zero_filling_combobox_selection == 1:
            # Have a textcontrol for the doubling number/number of zeros/final data size
            self.zf_value_label = wx.StaticText(parent, -1, "Number of zeros to add:")
            self.zero_filling_sizer.Add(
                self.zf_value_label, 0, wx.ALIGN_CENTER_VERTICAL
            )

            self.zero_filling_textcontrol = wx.TextCtrl(
                parent, -1, str(self.zero_filling_value_zeros_to_add), size=(40, 20)
            )
            self.zero_filling_textcontrol.Bind(
                wx.EVT_TEXT, self.on_zero_filling_textcontrol_zeros_to_add
            )
            self.zero_filling_sizer.AddSpacer(5)
            self.zero_filling_sizer.Add(
                self.zero_filling_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.zero_filling_sizer.AddSpacer(20)
        elif self.zero_filling_combobox_selection == 2:
            # Have a textcontrol for the doubling number/number of zeros/final data size
            self.zf_value_label = wx.StaticText(parent, -1, "Final data size:")
            self.zero_filling_sizer.Add(
                self.zf_value_label, 0, wx.ALIGN_CENTER_VERTICAL
            )

            self.zero_filling_textcontrol = wx.TextCtrl(
                parent, -1, str(self.zero_filling_value_final_data_size), size=(40, 20)
            )
            self.zero_filling_textcontrol.Bind(
                wx.EVT_TEXT, self.on_zero_filling_textcontrol_final_size
            )
            self.zero_filling_sizer.AddSpacer(5)
            self.zero_filling_sizer.Add(
                self.zero_filling_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.zero_filling_sizer.AddSpacer(20)

        # Have a checkbox for rounding to the nearest power of 2
        self.zero_filling_round_checkbox = wx.CheckBox(
            parent, -1, "Round to nearest power of 2"
        )
        self.zero_filling_round_checkbox.SetValue(True)
        self.zero_filling_sizer.Add(
            self.zero_filling_round_checkbox, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.zero_filling_sizer.AddSpacer(10)

        # Have a button showing information on zero filling
        self.zero_filling_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))
        self.zero_filling_info.Bind(wx.EVT_BUTTON, self.info_buttons.on_zero_fill_info)
        self.zero_filling_sizer.Add(self.zero_filling_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.zero_filling_sizer.AddSpacer(10)

        parent.sizer_1.Add(self.zero_filling_sizer)
        parent.sizer_1.AddSpacer(10)

    def on_zero_filling_checkbox(self, event):
        """
        When the zero filling checkbox is clicked, update
        the current stored value.
        """
        if self.zero_filling_checkbox.GetValue() == True:
            self.zero_filling_checkbox_value = True
        else:
            self.zero_filling_checkbox_value = False

    def on_zero_filling_textcontrol_doubling_times(self, event):
        """
        When a new value is typed into the doubling times box, check that
         the value is valid. If valid, update the stored parameter.
        """
        try:
            self.zero_filling_value_doubling_times = int(
                self.zero_filling_textcontrol.GetValue()
            )
        except:
            if self.zero_filling_textcontrol.GetValue() == "":
                self.zero_filling_value_doubling_times = 0
                self.zero_filling_textcontrol.SetValue("")
            else:
                # Give an popout error message saying that the value must be an integer
                self.zero_filling_value_doubling_times = 1
                self.zero_filling_textcontrol.SetValue(
                    str(self.zero_filling_value_doubling_times)
                )

                message = "The value for the zero filling doubling number must be an integer. Resetting value to 1."
                title = "Invalid value"
                style = wx.OK | wx.ICON_ERROR
                wx.MessageBox(message, title, style)

    def on_zero_filling_textcontrol_zeros_to_add(self, event):
        """
        When a new value is typed into the zeros to add box, check that
         the value is valid. If valid, update the stored parameter.
        """
        try:
            self.zero_filling_value_zeros_to_add = int(
                self.zero_filling_textcontrol.GetValue()
            )
        except:
            if self.zero_filling_textcontrol.GetValue() == "":
                self.zero_filling_value_zeros_to_add = 0
                self.zero_filling_textcontrol.SetValue("")
            else:
                # Give an popout error message saying that the value must be an integer
                self.zero_filling_value_zeros_to_add = 0
                self.zero_filling_textcontrol.SetValue(
                    str(self.zero_filling_value_zeros_to_add)
                )

                message = "The value for the zero filling (zeros to add) number must be an integer. Resetting value to 0."
                title = "Invalid value"
                style = wx.OK | wx.ICON_ERROR
                wx.MessageBox(message, title, style)

    def on_zero_filling_textcontrol_final_size(self, event):
        """
        When a new value is typed into the final size box, check that
         the value is valid. If valid, update the stored parameter.
        """
        try:
            self.zero_filling_value_final_data_size = int(
                self.zero_filling_textcontrol.GetValue()
            )
        except:
            if self.zero_filling_textcontrol.GetValue() == "":
                self.zero_filling_value_final_data_size = (
                    self.nmr_data.number_of_points[self.dimension_index] * 2
                )
                self.zero_filling_textcontrol.SetValue("")
            else:
                # Give an popout error message saying that the value must be an integer
                self.zero_filling_value_final_data_size = (
                    self.nmr_data.number_of_points[self.dimension_index] * 2
                )
                self.zero_filling_textcontrol.SetValue(
                    str(self.zero_filling_value_final_data_size)
                )

                message = """The value for the zero filling final data size must 
                be an integer. Resetting value to {}.""".format(
                    self.nmr_data.number_of_points[self.dimension_index] * 2
                )
                title = "Invalid value"
                style = wx.OK | wx.ICON_ERROR
                wx.MessageBox(message, title, style)

    def on_zero_filling_combobox(self, event):
        """
        When the zero fill combobox option is changed, need to clear
        all relevent sizers and then add them again (with the
        new version of the zero fill sizer)
        """

        self.zero_filling_combobox_selection = self.zero_filling_combobox.GetSelection()

        self.clear_zero_filling_sizer()

        # Remove the zf sizer
        self.zero_filling_sizer.Clear(delete_windows=True)
        # Within this apodization function is the necessary functionality
        # to refresh the sizes to the new values.
        self.apodization_class.on_apodization_combobox(wx.EVT_COMBOBOX)

        self.zero_filling_combobox_selection = self.zero_filling_combobox_selection_old

    def clear_zero_filling_sizer(self):
        """
        Function to clear the zero filling sizer
        """
        self.zero_filling_sizer.Clear()
        self.zero_filling_sizer.Detach(self.zero_filling_checkbox)
        self.zero_filling_checkbox.Destroy()
        self.zero_filling_sizer.Detach(self.zf_options_label)
        self.zf_options_label.Destroy()
        self.zero_filling_sizer.Detach(self.zero_filling_info)
        self.zero_filling_info.Destroy()
        self.zero_filling_sizer.Detach(self.zf_value_label)
        self.zf_value_label.Destroy()
        self.zero_filling_sizer.Detach(self.zero_filling_round_checkbox)
        self.zero_filling_round_checkbox.Destroy()
        self.zero_filling_sizer.Detach(self.zero_filling_textcontrol)
        self.zero_filling_textcontrol.Destroy()

        self.zero_filling_sizer.Detach(self.zero_filling_combobox)
        self.zero_filling_combobox.Hide()
