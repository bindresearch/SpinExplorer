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


class FourierTransform:

    def __init__(self, app, nmr_data, parent, info_buttons):
        """
        This class contains all the functions related to fourier transform
        parts of the graphical interface
        """

        self.app = app
        self.nmr_data = nmr_data
        self.info_buttons = info_buttons
        self.parent = parent

        self.set_initial_fourier_transform_variables()
        self.create_fourier_transform_sizer(parent)

    def set_initial_fourier_transform_variables(self):
        """
        Initialising variables relevent to the fourier transform
        processing interface
        """
        self.fourier_transform_checkbox_value = True
        self.ft_method_selection = 0  # Initially use the 'auto' method of FT as default

    def create_fourier_transform_sizer(self, parent):
        """
        Create a box for all the fourier transform options
        """
        self.fourier_transform_box = wx.StaticBox(parent, -1, "Fourier Transform")
        self.fourier_transform_sizer = wx.StaticBoxSizer(
            self.fourier_transform_box, wx.HORIZONTAL
        )
        self.fourier_transform_checkbox = wx.CheckBox(
            parent, -1, "Apply fourier transform"
        )
        self.fourier_transform_checkbox.SetValue(True)
        self.fourier_transform_sizer.Add(
            self.fourier_transform_checkbox, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.fourier_transform_sizer.AddSpacer(10)
        # Have a button for advanced options for fourier transform
        self.fourier_transform_advanced_options = wx.Button(
            parent, -1, "Advanced Options"
        )
        self.fourier_transform_advanced_options.Bind(
            wx.EVT_BUTTON, self.on_fourier_transform_advanced_options
        )
        self.fourier_transform_sizer.Add(
            self.fourier_transform_advanced_options, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.fourier_transform_sizer.AddSpacer(10)

        # Have a button showing information on fourier transform
        self.fourier_transform_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))
        self.fourier_transform_info.Bind(
            wx.EVT_BUTTON, self.info_buttons.on_fourier_transform_info
        )
        self.fourier_transform_sizer.Add(
            self.fourier_transform_info, 0, wx.ALIGN_CENTER_VERTICAL
        )
        parent.sizer_1.Add(self.fourier_transform_sizer)
        parent.sizer_1.AddSpacer(10)

    def on_fourier_transform_advanced_options(self, event):
        """
        Create a frame with a set of advanced options for the fourier transform
        implementation
        """
        self.fourier_transform_advanced_options_window = wx.Frame(
            self.parent, -1, "Fourier Transform Advanced Options", size=(700, 600)
        )

        self.fourier_transform_advanced_options_window_sizer = wx.BoxSizer(wx.VERTICAL)
        self.fourier_transform_advanced_options_window.SetSizer(
            self.fourier_transform_advanced_options_window_sizer
        )

        # Create a sizer for the fourier transform advanced options
        self.ft_label = wx.StaticBox(
            self.fourier_transform_advanced_options_window,
            -1,
            "Fourier Transform Method:",
        )
        self.fourier_transform_advanced_options_sizer = wx.StaticBoxSizer(
            self.ft_label, wx.VERTICAL
        )

        # Have a radiobox for auto, real, inverse, sign alternation
        self.fourier_transform_advanced_options_sizer.AddSpacer(10)
        self.fourier_transform_auto_real_inverse_sign_alternation_radio_box = (
            wx.RadioBox(
                self.fourier_transform_advanced_options_window,
                -1,
                choices=[
                    "Standard",
                    "Auto (not recommended)",
                    "Real",
                    "Inverse",
                    "Sign alternation (alt)",
                    "Negate imaginaries (neg)",
                    "alt + neg",
                ],
                style=wx.RA_SPECIFY_ROWS,
            )
        )
        self.fourier_transform_auto_real_inverse_sign_alternation_radio_box.SetSelection(
            self.ft_method_selection
        )
        self.fourier_transform_advanced_options_sizer.Add(
            self.fourier_transform_auto_real_inverse_sign_alternation_radio_box,
            0,
            wx.ALIGN_CENTER_HORIZONTAL,
        )
        self.fourier_transform_advanced_options_sizer.AddSpacer(10)

        self.ft_method_text = """        Standard: Perform a standard Fourier transform. \n\n
        Auto: The auto method will automatically select the best method for the fourier transform of the FID. \n\n
        Real: The Fourier Transform will be applied to the real part of the FID only. (For TPPI) \n\n
        Inverse: The inverse Fourier  Transform will be applied to the FID. \n\n
        Sign Alternation: The sign alternation method will be applied to the FID. (For States-TPPI) \n\n
        Negate imaginaries: Chnage sign of imaginaries before Fourier Transform\n\n"""

        self.ft_method_info = wx.StaticText(
            self.fourier_transform_advanced_options_window, -1, self.ft_method_text
        )
        self.fourier_transform_advanced_options_sizer.Add(
            self.ft_method_info, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.fourier_transform_advanced_options_sizer.AddSpacer(10)

        # Have a save and close button
        self.fourier_transform_advanced_options_save_button = wx.Button(
            self.fourier_transform_advanced_options_window, -1, "Save and Close"
        )
        self.fourier_transform_advanced_options_save_button.Bind(
            wx.EVT_BUTTON, self.on_fourier_transform_advanced_options_save
        )
        self.fourier_transform_advanced_options_sizer.Add(
            self.fourier_transform_advanced_options_save_button,
            0,
            wx.ALIGN_CENTER_HORIZONTAL,
        )

        self.fourier_transform_advanced_options_window_sizer.Add(
            self.fourier_transform_advanced_options_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL
        )

        self.fourier_transform_advanced_options_window.Show()

    def on_fourier_transform_advanced_options_save(self, event):
        """
        Saving the current selection and closing the window
        """
        self.ft_method_selection = (
            self.fourier_transform_auto_real_inverse_sign_alternation_radio_box.GetSelection()
        )
        self.fourier_transform_advanced_options_window.Close()
