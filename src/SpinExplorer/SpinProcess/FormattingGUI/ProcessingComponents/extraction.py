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


class Extraction:

    def __init__(self, app, nmr_data, parent, info_buttons):
        """
        This class contains all the functions related to data extraction
        parts of the graphical interface.
        """

        self.app = app
        self.nmr_data = nmr_data
        self.info_buttons = info_buttons

        self.set_initial_extraction_variables()
        self.create_extraction_sizer(parent)

    def set_initial_extraction_variables(self):
        """
        Initialising parameters necessary for the extraction
        section of the GUI
        """
        self.extraction_checkbox_value = False
        self.extraction_ppm_start = 0.0
        self.extraction_ppm_end = 0.0

    def create_extraction_sizer(self, parent):
        """
        A box for a user to choose the chemical shifts to extract
        between.
        """
        self.extraction_box = wx.StaticBox(parent, -1, "Extraction")
        self.extraction_sizer = wx.StaticBoxSizer(self.extraction_box, wx.HORIZONTAL)
        self.extraction_checkbox = wx.CheckBox(parent, -1, "Include data extraction")
        self.extraction_checkbox.Bind(wx.EVT_CHECKBOX, self.on_extraction_checkbox)
        self.extraction_checkbox.SetValue(self.extraction_checkbox_value)
        self.extraction_sizer.Add(self.extraction_checkbox, 0, wx.ALIGN_CENTER_VERTICAL)
        self.extraction_sizer.AddSpacer(10)
        # Have a textcontrol for the ppm start value
        self.extraction_ppm_start_label = wx.StaticText(
            parent, -1, "Start chemical shift (ppm):"
        )
        self.extraction_sizer.Add(
            self.extraction_ppm_start_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.extraction_ppm_start_textcontrol = wx.TextCtrl(
            parent, -1, str(self.extraction_ppm_start), size=(40, 20)
        )
        self.extraction_ppm_start_textcontrol.Bind(
            wx.EVT_TEXT, self.on_extraction_textcontrol
        )
        self.extraction_sizer.Add(
            self.extraction_ppm_start_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.extraction_sizer.AddSpacer(10)
        # Have a textcontrol for the ppm end value
        self.extraction_ppm_end_label = wx.StaticText(
            parent, -1, "End chemical shift (ppm):"
        )
        self.extraction_sizer.Add(
            self.extraction_ppm_end_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.extraction_ppm_end_textcontrol = wx.TextCtrl(
            parent, -1, str(self.extraction_ppm_end), size=(40, 20)
        )
        self.extraction_ppm_end_textcontrol.Bind(
            wx.EVT_TEXT, self.on_extraction_textcontrol
        )
        self.extraction_sizer.Add(
            self.extraction_ppm_end_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.extraction_sizer.AddSpacer(10)
        # Have a button showing information on extraction
        self.extraction_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))
        self.extraction_info.Bind(wx.EVT_BUTTON, self.info_buttons.on_extraction_info)
        self.extraction_sizer.Add(self.extraction_info, 0, wx.ALIGN_CENTER_VERTICAL)
        parent.sizer_1.Add(self.extraction_sizer)
        parent.sizer_1.AddSpacer(10)

    def on_extraction_checkbox(self, event):
        """
        When the extraction checkbox is clicked, update the stored value
        """
        self.extraction_checkbox_value = self.extraction_checkbox.GetValue()

    def on_extraction_textcontrol(self, event):
        """
        When a user updates the extraction start and end chemical shifts
        this function updates their stored values.
        """
        self.extraction_ppm_start = self.extraction_ppm_start_textcontrol.GetValue()
        self.extraction_ppm_end = self.extraction_ppm_end_textcontrol.GetValue()
