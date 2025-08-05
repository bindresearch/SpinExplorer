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


class BaselineCorrection:

    def __init__(self, app, nmr_data, parent, info_buttons):
        """
        This class contains all the functions related to baseline
        correction parts of the graphical interface.
        """

        self.app = app
        self.nmr_data = nmr_data
        self.info_buttons = info_buttons

        self.set_initial_baseline_correction_variables()
        self.create_baseline_correction_sizer(parent)

    def set_initial_baseline_correction_variables(self):
        """
        Initialising relevent variables for the baseline
        correction section of the graphical interface.
        """
        self.baseline_correction_radio_box_selection = 0
        self.node_width = "2"
        self.node_list = "0,5,95,100"
        self.polynomial_order = "4"

    def create_baseline_correction_sizer(self, parent):
        """
        Creating a box for baseline correction options.
        A user can select either linear or polynomial.
        There are also sections for the user to edit the
        node-width, node-list and polynomial order.
        """

        self.baseline_correction_box = wx.StaticBox(parent, -1, "Baseline Correction")
        self.baseline_correction_sizer = wx.StaticBoxSizer(
            self.baseline_correction_box, wx.HORIZONTAL
        )
        self.baseline_correction_checkbox = wx.CheckBox(
            parent, -1, "Apply baseline correction"
        )
        self.baseline_correction_checkbox.SetValue(False)
        self.baseline_correction_sizer.Add(
            self.baseline_correction_checkbox, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.baseline_correction_sizer.AddSpacer(10)
        # Have a radio box for linear or polynomial baseline correction
        self.baseline_correction_radio_box = wx.RadioBox(
            parent, -1, "Baseline Correction Method", choices=["Linear", "Polynomial"]
        )
        # Bind the radio box to a function that will update the baseline correction options
        self.baseline_correction_radio_box.Bind(
            wx.EVT_RADIOBOX, self.on_baseline_correction_radio_box
        )
        self.baseline_correction_radio_box.SetSelection(
            self.baseline_correction_radio_box_selection
        )
        self.baseline_correction_sizer.Add(
            self.baseline_correction_radio_box, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.baseline_correction_sizer.AddSpacer(10)

        # If linear baseline correction is selected, have a textcontrol for the node values to use
        self.baseline_correction_nodes_label = wx.StaticText(
            parent, -1, "Node width (pts):"
        )
        self.baseline_correction_sizer.Add(
            self.baseline_correction_nodes_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.baseline_correction_nodes_textcontrol = wx.TextCtrl(
            parent, -1, str(self.node_width), size=(30, 20)
        )
        self.baseline_correction_nodes_textcontrol.Bind(
            wx.EVT_TEXT, self.on_baseline_correction_textcontrol
        )
        self.baseline_correction_sizer.Add(
            self.baseline_correction_nodes_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.baseline_correction_sizer.AddSpacer(10)
        # Have a textcontrol for the node list (percentages)
        self.baseline_correction_node_list_label = wx.StaticText(
            parent, -1, "Node list (%):"
        )
        self.baseline_correction_sizer.Add(
            self.baseline_correction_node_list_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.baseline_correction_node_list_textcontrol = wx.TextCtrl(
            parent, -1, str(self.node_list), size=(100, 20)
        )
        self.baseline_correction_node_list_textcontrol.Bind(
            wx.EVT_TEXT, self.on_baseline_correction_textcontrol
        )
        self.baseline_correction_sizer.Add(
            self.baseline_correction_node_list_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.baseline_correction_sizer.AddSpacer(10)
        # If polynomial baseline correction is selected, have a textcontrol for the polynomial order

        self.baseline_correction_polynomial_order_label = wx.StaticText(
            parent, -1, "Polynomial order:"
        )
        self.baseline_correction_sizer.Add(
            self.baseline_correction_polynomial_order_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.baseline_correction_polynomial_order_textcontrol = wx.TextCtrl(
            parent, -1, str(self.polynomial_order), size=(30, 20)
        )
        self.baseline_correction_polynomial_order_textcontrol.Bind(
            wx.EVT_TEXT, self.on_baseline_correction_textcontrol
        )
        self.baseline_correction_sizer.Add(
            self.baseline_correction_polynomial_order_textcontrol,
            0,
            wx.ALIGN_CENTER_VERTICAL,
        )
        self.baseline_correction_sizer.AddSpacer(10)

        if self.baseline_correction_radio_box_selection == 0:
            self.baseline_correction_polynomial_order_label.Hide()
            self.baseline_correction_polynomial_order_textcontrol.Hide()

        # Have a button showing information on baseline correction
        self.baseline_correction_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))

        self.baseline_correction_info.Bind(
            wx.EVT_BUTTON, self.info_buttons.on_baseline_correction_info
        )
        self.baseline_correction_sizer.Add(
            self.baseline_correction_info, 0, wx.ALIGN_CENTER_VERTICAL
        )
        parent.sizer_1.Add(self.baseline_correction_sizer)
        parent.sizer_1.AddSpacer(10)

    def on_baseline_correction_textcontrol(self, event):
        """
        When a user updates the values in the boxes, the stored
        values are then updated.
        """
        self.node_width = self.baseline_correction_nodes_textcontrol.GetValue()
        self.node_list = self.baseline_correction_node_list_textcontrol.GetValue()
        self.polynomial_order = (
            self.baseline_correction_polynomial_order_textcontrol.GetValue()
        )

    def on_baseline_correction_radio_box(self, event):
        """
        When the radio box is changed between linear and polynomial then
        update the baseline correction box to only show the parameters
        relevent to that type of baseline correction.
        """
        # If the user selects linear or polynomial baseline correction, update the options
        self.baseline_correction_radio_box_selection = (
            self.baseline_correction_radio_box.GetSelection()
        )

        if self.baseline_correction_radio_box_selection == 0:
            # Remove the polynomial order textcontrol
            self.baseline_correction_sizer.Hide(
                self.baseline_correction_polynomial_order_label
            )
            self.baseline_correction_sizer.Hide(
                self.baseline_correction_polynomial_order_textcontrol
            )
            self.baseline_correction_sizer.Layout()
        elif self.baseline_correction_radio_box_selection == 1:
            # Add the polynomial order textcontrol
            self.baseline_correction_sizer.Show(
                self.baseline_correction_polynomial_order_label
            )
            self.baseline_correction_sizer.Show(
                self.baseline_correction_polynomial_order_textcontrol
            )
            self.baseline_correction_sizer.Layout()

        self.app.Refresh()
        self.app.Update()
        self.app.Layout()
