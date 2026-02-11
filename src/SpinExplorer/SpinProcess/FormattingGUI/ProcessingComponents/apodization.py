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

import matplotlib


import matplotlib.pyplot as plt
import numpy as np
import copy
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas


class Apodization:

    def __init__(self, app, nmr_data, parent, info_buttons, other_classes, dimension):
        """
        This class contains all the functions related to apodization
        or window functions
        """

        self.app = app
        self.nmr_data = nmr_data
        self.info_buttons = info_buttons
        if dimension == 0:
            self.solvent_suppression_class = other_classes[0]
            self.linear_prediction_class = other_classes[1]
        else:
            self.linear_prediction_class = other_classes[0]

        self.dimension_index = dimension
        self.parent = parent

        self.set_initial_apodization_variables()
        self.create_apodization_sizer(self.parent)

    def set_initial_apodization_variables(self):
        self.apodization_checkbox_value = True
        self.apodization_combobox_selection = 1
        self.apodization_combobox_selection_old = 1

        # Initial values for exponential apodization
        self.exponential_line_broadening = 0.5
        self.apodization_first_point_scaling = 0.5

        # Initial values for Lorentz to Gauss apodization
        self.g1 = 0.33
        self.g2 = 1
        self.g3 = 0.0

        # Initial values for Sinebell apodization
        self.offset = 0.5
        self.end = 0.98
        self.power = 1.0

        # Initial values for Gauss Broadening apodization
        self.a = -1.0 
        self.b = 0.2

        # Initial values for Trapezoid apodization
        self.t1 = int((self.nmr_data.number_of_points[self.dimension_index] / 2) / 4)
        self.t2 = int((self.nmr_data.number_of_points[self.dimension_index] / 2) / 4)

        # Initial values for Triangle apodization
        self.loc = 0.5

    def create_apodization_sizer(self, parent):
        """
        Create a box for apodization options
        """
        self.apodization_box = wx.StaticBox(parent, -1, "Apodization")
        self.apodization_sizer = wx.StaticBoxSizer(self.apodization_box, wx.HORIZONTAL)
        self.apodization_checkbox = wx.CheckBox(parent, -1, "Apply apodization")
        self.apodization_checkbox.SetValue(self.apodization_checkbox_value)
        self.apodization_checkbox.Bind(wx.EVT_CHECKBOX, self.on_apodization_checkbox)
        self.apodization_sizer.Add(
            self.apodization_checkbox, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.apodization_sizer.AddSpacer(10)
        # Have a combobox for apodization options
        self.apodization_options = [
            "None",
            "Exponential",
            "Lorentz to Gauss",
            "Sinebell",
            "Gauss Broadening",
            "Trapazoid",
            "Triangle",
        ]
        self.apodization_combobox = wx.ComboBox(
            parent, -1, choices=self.apodization_options, style=wx.CB_READONLY
        )
        self.apodization_combobox.SetSelection(self.apodization_combobox_selection)
        self.apodization_combobox.Bind(wx.EVT_COMBOBOX, self.on_apodization_combobox)
        self.apodization_sizer.Add(
            self.apodization_combobox, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.apodization_sizer.AddSpacer(10)
        if self.apodization_combobox_selection == 1:
            # Have a textcontrol for the line broadening
            self.apodization_line_broadening_label = wx.StaticText(
                parent, -1, "Line Broadening (Hz):"
            )
            self.apodization_sizer.Add(
                self.apodization_line_broadening_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_line_broadening_textcontrol = wx.TextCtrl(
                parent, -1, str(self.exponential_line_broadening), size=(30, 20)
            )
            self.apodization_line_broadening_textcontrol.Bind(
                wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol
            )
            self.apodization_sizer.Add(
                self.apodization_line_broadening_textcontrol,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(
                parent, -1, "First Point Scaling:"
            )
            self.apodization_sizer.Add(
                self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_first_point_textcontrol = wx.TextCtrl(
                parent, -1, str(self.apodization_first_point_scaling), size=(30, 20)
            )
            self.apodization_sizer.Add(
                self.apodization_first_point_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer.AddSpacer(10)
        elif self.apodization_combobox_selection == 2:
            # Have a textcontrol for the g1 value
            self.apodization_g1_label = wx.StaticText(
                parent, -1, "Inverse Lorentzian (Hz):"
            )
            self.apodization_sizer.Add(
                self.apodization_g1_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_g1_textcontrol = wx.TextCtrl(
                parent, -1, str(self.g1), size=(40, 20)
            )
            self.apodization_g1_textcontrol.Bind(
                wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol
            )
            self.apodization_sizer.Add(
                self.apodization_g1_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the g2 value
            self.apodization_g2_label = wx.StaticText(
                parent, -1, "Gaussian Broadening (Hz):"
            )
            self.apodization_sizer.Add(
                self.apodization_g2_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_g2_textcontrol = wx.TextCtrl(
                parent, -1, str(self.g2), size=(40, 20)
            )
            self.apodization_g2_textcontrol.Bind(
                wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol
            )
            self.apodization_sizer.Add(
                self.apodization_g2_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the g3 value
            self.apodization_g3_label = wx.StaticText(parent, -1, "Gaussian Shift:")
            self.apodization_sizer.Add(
                self.apodization_g3_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_g3_textcontrol = wx.TextCtrl(
                parent, -1, str(self.g3), size=(40, 20)
            )
            self.apodization_g3_textcontrol.Bind(
                wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol
            )
            self.apodization_sizer.Add(
                self.apodization_g3_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(
                parent, -1, "First Point Scaling:"
            )
            self.apodization_sizer.Add(
                self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_first_point_textcontrol = wx.TextCtrl(
                parent, -1, str(self.apodization_first_point_scaling), size=(30, 20)
            )
            self.apodization_sizer.Add(
                self.apodization_first_point_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer.AddSpacer(10)
        elif self.apodization_combobox_selection == 3:
            # Have a textcontrol for the offset value
            self.apodization_offset_label = wx.StaticText(
                parent, -1, "Offset (\u03c0):"
            )
            self.apodization_sizer.Add(
                self.apodization_offset_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_offset_textcontrol = wx.TextCtrl(
                parent, -1, str(self.offset), size=(40, 20)
            )
            self.apodization_offset_textcontrol.Bind(
                wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol
            )
            self.apodization_sizer.Add(
                self.apodization_offset_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the end value
            self.apodization_end_label = wx.StaticText(parent, -1, "End (\u03c0):")
            self.apodization_sizer.Add(
                self.apodization_end_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_end_textcontrol = wx.TextCtrl(
                parent, -1, str(self.end), size=(40, 20)
            )
            self.apodization_end_textcontrol.Bind(
                wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol
            )
            self.apodization_sizer.Add(
                self.apodization_end_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the power value
            self.apodization_power_label = wx.StaticText(parent, -1, "Power:")
            self.apodization_sizer.Add(
                self.apodization_power_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_power_textcontrol = wx.TextCtrl(
                parent, -1, str(self.power), size=(30, 20)
            )
            self.apodization_power_textcontrol.Bind(
                wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol
            )
            self.apodization_sizer.Add(
                self.apodization_power_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(
                parent, -1, "First Point Scaling:"
            )
            self.apodization_sizer.Add(
                self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_first_point_textcontrol = wx.TextCtrl(
                parent, -1, str(self.apodization_first_point_scaling), size=(30, 20)
            )
            self.apodization_sizer.Add(
                self.apodization_first_point_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer.AddSpacer(10)
        elif self.apodization_combobox_selection == 4:
            # Have a textcontrol for the a value
            self.apodization_a_label = wx.StaticText(
                parent, -1, "Line Broadening (Hz):"
            )
            self.apodization_sizer.Add(
                self.apodization_a_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_a_textcontrol = wx.TextCtrl(
                parent, -1, str(self.a), size=(40, 20)
            )
            self.apodization_a_textcontrol.Bind(
                wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol
            )
            self.apodization_sizer.Add(
                self.apodization_a_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the b value
            self.apodization_b_label = wx.StaticText(
                parent, -1, "GB factor (0.0-1.0):"
            )
            self.apodization_sizer.Add(
                self.apodization_b_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_b_textcontrol = wx.TextCtrl(
                parent, -1, str(self.b), size=(40, 20)
            )
            self.apodization_b_textcontrol.Bind(
                wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol
            )
            self.apodization_sizer.Add(
                self.apodization_b_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(
                parent, -1, "First Point Scaling:"
            )
            self.apodization_sizer.Add(
                self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_first_point_textcontrol = wx.TextCtrl(
                parent, -1, str(self.apodization_first_point_scaling), size=(30, 20)
            )
            self.apodization_sizer.Add(
                self.apodization_first_point_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer.AddSpacer(10)
        elif self.apodization_combobox_selection == 5:
            # Have a textcontrol for the t1 value
            self.apodization_t1_label = wx.StaticText(parent, -1, "Ramp up points:")
            self.apodization_sizer.Add(
                self.apodization_t1_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_t1_textcontrol = wx.TextCtrl(
                parent, -1, str(self.t1), size=(50, 20)
            )
            self.apodization_t1_textcontrol.Bind(
                wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol
            )
            self.apodization_sizer.Add(
                self.apodization_t1_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the t2 value
            self.apodization_t2_label = wx.StaticText(parent, -1, "Ramp down points:")
            self.apodization_sizer.Add(
                self.apodization_t2_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_t2_textcontrol = wx.TextCtrl(
                parent, -1, str(self.t2), size=(50, 20)
            )
            self.apodization_t2_textcontrol.Bind(
                wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol
            )
            self.apodization_sizer.Add(
                self.apodization_t2_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(
                parent, -1, "First Point Scaling:"
            )
            self.apodization_sizer.Add(
                self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_first_point_textcontrol = wx.TextCtrl(
                parent, -1, str(self.apodization_first_point_scaling), size=(30, 20)
            )
            self.apodization_sizer.Add(
                self.apodization_first_point_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer.AddSpacer(10)
        elif self.apodization_combobox_selection == 6:
            # Have a textcontrol for the loc value
            self.apodization_loc_label = wx.StaticText(
                parent, -1, "Location of maximum:"
            )
            self.apodization_sizer.Add(
                self.apodization_loc_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_loc_textcontrol = wx.TextCtrl(
                parent, -1, str(self.loc), size=(40, 20)
            )
            self.apodization_loc_textcontrol.Bind(
                wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol
            )
            self.apodization_sizer.Add(
                self.apodization_loc_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(
                parent, -1, "First Point Scaling:"
            )
            self.apodization_sizer.Add(
                self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_first_point_textcontrol = wx.TextCtrl(
                parent, -1, "0.5", size=(30, 20)
            )
            self.apodization_sizer.Add(
                self.apodization_first_point_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer.AddSpacer(10)

        # Have a button for information on currently selected apodization containing unicode i in a circle
        self.apodization_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))
        self.apodization_info.Bind(wx.EVT_BUTTON, self.info_buttons.on_apodization_info)
        self.apodization_sizer.Add(self.apodization_info, 0, wx.ALIGN_CENTER_VERTICAL)

        # Have a mini plot of the apodization function along with the FID first slice
        self.plot_window_function()

        self.parent.sizer_1.Add(self.apodization_sizer)
        self.parent.sizer_1.AddSpacer(10)

    def on_apodization_checkbox(self, event):
        """
        When the apodization checkbox is pressed, update
        the checkbox value
        """
        if self.apodization_checkbox.GetValue() == True:
            self.apodization_checkbox_value = True
        else:
            self.apodization_checkbox_value = False

    def on_apodization_textcontrol(self, event):
        """
        If the user presses enter after typing in
        a new apodization parameter value, update the plot
        """
        keycode = event.GetKeyCode()
        if keycode == wx.WXK_RETURN:
            self.update_window_function_plot()
        event.Skip()

    def on_apodization_combobox(self, event):
        """
        If the apodization combobox value is changed, update the GUI
        to represent the new apodization type.
        e.g. exponential line broadening needs the parameter lb (line broadening)
        but Gaussian broadening needs a different set of parameters
        """
        self.apodization_combobox_selection = self.apodization_combobox.GetSelection()

        # Destroy the combobox and textcontrols for the previous apodization function
        # self.apodization_sizer.Detach(self.apodization_combobox)
        # self.apodization_combobox.Destroy()

        if self.apodization_combobox_selection_old == 1:
            # Remove the previous textcontrols

            self.apodization_sizer.Detach(self.apodization_line_broadening_label)
            self.apodization_line_broadening_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_line_broadening_textcontrol)
            self.apodization_line_broadening_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_textcontrol)
            self.apodization_first_point_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)

        elif self.apodization_combobox_selection_old == 2:
            self.apodization_sizer.Detach(self.apodization_g1_label)
            self.apodization_g1_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_g1_textcontrol)
            self.apodization_g1_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_g2_label)
            self.apodization_g2_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_g2_textcontrol)
            self.apodization_g2_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_g3_label)
            self.apodization_g3_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_g3_textcontrol)
            self.apodization_g3_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_textcontrol)
            self.apodization_first_point_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)

        elif self.apodization_combobox_selection_old == 3:
            self.apodization_sizer.Detach(self.apodization_offset_label)
            self.apodization_offset_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_offset_textcontrol)
            self.apodization_offset_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_end_label)
            self.apodization_end_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_end_textcontrol)
            self.apodization_end_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_power_label)
            self.apodization_power_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_power_textcontrol)
            self.apodization_power_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_textcontrol)
            self.apodization_first_point_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)
        elif self.apodization_combobox_selection_old == 4:
            self.apodization_sizer.Detach(self.apodization_a_label)
            self.apodization_a_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_a_textcontrol)
            self.apodization_a_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_b_label)
            self.apodization_b_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_b_textcontrol)
            self.apodization_b_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_textcontrol)
            self.apodization_first_point_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)
        elif self.apodization_combobox_selection_old == 5:
            self.apodization_sizer.Detach(self.apodization_t1_label)
            self.apodization_t1_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_t1_textcontrol)
            self.apodization_t1_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_t2_label)
            self.apodization_t2_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_t2_textcontrol)
            self.apodization_t2_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_textcontrol)
            self.apodization_first_point_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)
        elif self.apodization_combobox_selection_old == 6:
            self.apodization_sizer.Detach(self.apodization_loc_label)
            self.apodization_loc_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_loc_textcontrol)
            self.apodization_loc_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_textcontrol)
            self.apodization_first_point_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)
        elif self.apodization_combobox_selection_old == 0:
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)

        self.apodization_sizer.Detach(self.apodization_checkbox)
        self.apodization_checkbox.Destroy()
        self.apodization_sizer.Detach(self.apodization_combobox)
        self.apodization_combobox.Hide()

        # Delete the current apodization sizer and then create a new one

        self.apodization_sizer.Detach(self.apodization_plot_sizer)
        self.apodization_plot_ax.clear()
        self.apodization_plot.clear()
        self.apodization_plot_sizer.Clear(True)

        self.parent.sizer_1.Remove(self.apodization_sizer)
        # self.apodization_sizer.Clear(delete_windows=True)

        # self.sizer_1.Remove(self.linear_prediction_sizer)

        if self.dimension_index == 0:
            # Remove the linear prediction sizers
            self.linear_prediction_class.linear_prediction_sizer.Clear(
                delete_windows=True
            )
            # Remove the solvent suppression sizers
            self.solvent_suppression_class.solvent_suppression_sizer.Clear(
                delete_windows=True
            )
        else:
            self.linear_prediction_class.linear_prediction_sizer_indirect.Clear(
                delete_windows=True
            )
            # self.sizer_1.Remove(self.solvent_suppression_sizer)

        self.parent.sizer_1.Clear(delete_windows=True)

        self.parent.refresh_menu_bar()
        self.app.Refresh()
        self.app.Update()
        self.app.Layout()

        self.apodization_combobox_selection_old = self.apodization_combobox_selection

    def OnPressWindow(self, event):
        # Create a matplotlib popout of the plot
        fig = plt.figure(facecolor="white")
        ax = fig.add_subplot(111)
        ax.set_facecolor("white")
        ax.spines["bottom"].set_color("k")
        ax.spines["top"].set_color("k")
        ax.spines["right"].set_color("k")
        ax.spines["left"].set_color("k")
        ax.tick_params(axis="x", colors="k")
        ax.tick_params(axis="y", colors="k")
        ax.yaxis.label.set_color("k")
        ax.xaxis.label.set_color("k")
        line1_x, line1_y = self.line1.get_data()
        ax.plot(line1_x, line1_y, color="#1f77b4")
        if self.dimension_index == 0:
            line2_x, line2_y = self.line2.get_data()
            ax.plot(line2_x, line2_y, color="k")
        ax.set_xlim(
            -(self.npoints)
            / self.nmr_data.spectral_width[self.dimension_index]
            / 20,
            (self.npoints)
            / self.nmr_data.spectral_width[self.dimension_index]
            * 21
            / 20,
        )
        ax.set_ylim(-1.5, 1.5)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Normalised Intensity (a.u.)")
        plt.show()

    def plot_window_function(self):
        self.apodization_plot_sizer = wx.BoxSizer(wx.VERTICAL)
        self.apodization_plot = Figure(figsize=(1, 0.5), facecolor="#e5e6e7")
        self.apodization_plot_ax = self.apodization_plot.add_subplot(111)
        self.apodization_plot_canvas = FigCanvas(self.parent, -1, self.apodization_plot)
        self.apodization_plot_canvas.mpl_connect(
            "button_press_event", self.OnPressWindow
        )

        self.apodization_plot_ax.set_xticks([])
        self.apodization_plot_ax.set_yticks([])

        # If the apodization function is None, make remove the axes of the plot
        if self.apodization_combobox_selection == 0:
            self.apodization_plot_ax.spines["top"].set_visible(False)
            self.apodization_plot_ax.spines["right"].set_visible(False)
            self.apodization_plot_ax.spines["bottom"].set_visible(False)
            self.apodization_plot_ax.spines["left"].set_visible(False)

        self.plot_window_function_input()

    def plot_window_function_input(self):

        # Is the digital filter removed before or after processing
        self.before_processing = False
        try:
            with open("fid.com", "r") as file:
                data = file.readlines()
                for line in data:
                    if "-AMX" in line.split():
                        self.before_processing = False
                        break
                    elif "-DMX" in line.split():
                        self.before_processing = True
                        break
        except:
            self.before_processing = True

        # Plot the NMR FID
        if self.nmr_data.fid_file == "fids":
            for i, (dic, plane) in enumerate(self.nmr_data.xiter):
                self.data = plane
                break
        else:
            self.data = self.nmr_data.data

        if self.dimension_index == 0:
            if self.nmr_data.dim == 1:
                self.data = self.data
            elif self.nmr_data.dim == 2 or self.nmr_data.fid_file == "fids":
                self.data = self.data[0]
            else:
                self.data = self.data[0][0]

        data = copy.deepcopy(self.data)

        

        if(self.dimension_index==0):
            self.npoints = self.nmr_data.number_of_points[self.dimension_index]
            x = np.linspace(
                0,
                (self.nmr_data.number_of_points[self.dimension_index])
                / self.nmr_data.spectral_width[self.dimension_index],
                int(self.nmr_data.number_of_points[self.dimension_index]),
            )
        else:
            self.npoints = int(self.nmr_data.number_of_points[self.dimension_index]/2)
            x = np.linspace(
                0,
                (self.nmr_data.number_of_points[self.dimension_index]/2)
                / self.nmr_data.spectral_width[self.dimension_index],
                int(self.nmr_data.number_of_points[self.dimension_index]/2),
            )

        x_data = np.linspace(
            0,
            (len(data)) / self.nmr_data.spectral_width[self.dimension_index],
            int(len(data)),
        )

        # if(self.before_processing == False):
        #     x=np.linspace(0, (self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0], int(self.nmr_data.number_of_points[0]/2))
        # else:
        #     x=np.linspace(0, (len(data))/self.nmr_data.spectral_width[0], int(len(data)))

        if self.dimension_index == 0:
            (self.line2,) = self.apodization_plot_ax.plot(
                x_data, data / max(data), color="k"
            )
        if self.apodization_combobox_selection == 1:
            # Exponential window function
            (self.line1,) = self.apodization_plot_ax.plot(
                x,
                np.exp(-(np.pi * x * self.exponential_line_broadening)),
                color="#1f77b4",
            )
            self.apodization_plot_ax.set_ylim(-1.5, 1.5)
            self.apodization_plot_ax.set_xlim(
                -(self.npoints)
                / self.nmr_data.spectral_width[self.dimension_index]
                / 20,
                (self.npoints)
                / self.nmr_data.spectral_width[self.dimension_index]
                * 21
                / 20,
            )

        elif self.apodization_combobox_selection == 2:
            # Lorentz to Gauss window function
            g1 = self.g1/self.nmr_data.spectral_width[self.dimension_index]
            g2 = self.g2/self.nmr_data.spectral_width[self.dimension_index]
            x1 = np.arange(self.npoints)
            e = (
                np.pi
                * x1
                * g1
            )
            g = (
                0.6
                * np.pi
                * g2
                * (
                    self.g3
                    * (
                        self.npoints
                        - 1
                    )
                    - x1
                )
            )
            func = np.exp(e - g * g)
            (self.line1,) = self.apodization_plot_ax.plot(x, func, color="#1f77b4")
            self.apodization_plot_ax.set_ylim(-1.5, 1.5)
            self.apodization_plot_ax.set_xlim(
                -(self.npoints)
                / self.nmr_data.spectral_width[self.dimension_index]
                / 20,
                (self.npoints)
                / self.nmr_data.spectral_width[self.dimension_index]
                * 21
                / 20,
            )
        elif self.apodization_combobox_selection == 3:
            # Sinebell window function
            x1 = np.arange(self.npoints)
            func = np.flip(
                np.sin(
                    (np.pi * self.offset + np.pi * (self.end - self.offset) * x1)
                    / (self.npoints-1
                    )
                )
                ** self.power
            )
            # func = (
            #     np.sin(
            #         (np.pi * self.offset + np.pi * (self.end - self.offset) * x)
            #         / (
            #             (
            #                 (
            #                     (
            #                         self.nmr_data.number_of_points[self.dimension_index]
            #                         / 2
            #                     )
            #                     / self.nmr_data.spectral_width[self.dimension_index]
            #                 )
            #             )
            #         )
            #     )
            #     ** self.power
            # )
            (self.line1,) = self.apodization_plot_ax.plot(x, func, color="#1f77b4")
            self.apodization_plot_ax.set_ylim(-1.5, 1.5)
            self.apodization_plot_ax.set_xlim(
                -(self.npoints)
                / self.nmr_data.spectral_width[self.dimension_index]
                / 20,
                (self.npoints)
                / self.nmr_data.spectral_width[self.dimension_index]
                * 21
                / 20,
            )
        elif self.apodization_combobox_selection == 4:
            # Gauss broadening window function
            x1 = np.arange(self.npoints)
            t = x1/self.nmr_data.spectral_width[self.dimension_index]
            aq = self.npoints/self.nmr_data.spectral_width[self.dimension_index]

            a = np.pi * self.a 
            b = -a / (2.0 * self.b * aq)
            x1 = np.arange(self.npoints)
            func = np.exp(-a * t - (b * (t**2)))
            (self.line1,) = self.apodization_plot_ax.plot(x, func, color="#1f77b4")
            self.apodization_plot_ax.set_ylim(-1.5, 1.5)
            self.apodization_plot_ax.set_xlim(
                -(self.npoints)
                / self.nmr_data.spectral_width[self.dimension_index]
                / 20,
                (self.npoints)
                / self.nmr_data.spectral_width[self.dimension_index]
                * 21
                / 20,
            )
        elif self.apodization_combobox_selection == 5:
            # Trapazoid window function
            func = np.concatenate(
                (
                    np.linspace(0, 1, int(self.t1)),
                    np.ones(
                        int(self.npoints)
                        - int(self.t1)
                        - int(self.t2)
                    ),
                    np.linspace(1, 0, int(self.t2)),
                )
            )
            (self.line1,) = self.apodization_plot_ax.plot(x, func, color="#1f77b4")
            self.apodization_plot_ax.set_ylim(-1.5, 1.5)
            self.apodization_plot_ax.set_xlim(
                -(self.npoints)
                / self.nmr_data.spectral_width[self.dimension_index]
                / 20,
                (self.npoints)
                / self.nmr_data.spectral_width[self.dimension_index]
                * 21
                / 20,
            )
        elif self.apodization_combobox_selection == 6:
            # Triangle window function
            func = np.concatenate(
                (
                    np.linspace(
                        0,
                        1,
                        int(
                            self.loc
                            * (self.npoints)
                        ),
                    ),
                    np.linspace(
                        1,
                        0,
                        int(
                            (1 - self.loc)
                            * (self.npoints)
                        ),
                    ),
                )
            )
            (self.line1,) = self.apodization_plot_ax.plot(x, func, color="#1f77b4")

            self.apodization_plot_ax.set_ylim(-1.5, 1.5)
            self.apodization_plot_ax.set_xlim(
                -(self.npoints)
                / self.nmr_data.spectral_width[self.dimension_index]
                / 20,
                (self.npoints)
                / self.nmr_data.spectral_width[self.dimension_index]
                * 21
                / 20,
            )

        self.apodization_plot_ax.set_xlim(
            -(self.npoints)
            / self.nmr_data.spectral_width[self.dimension_index]
            / 20,
            (self.npoints)
            / self.nmr_data.spectral_width[self.dimension_index]
            * 21
            / 20,
        )

        self.apodization_plot_sizer.Add(self.apodization_plot_canvas, 0, wx.EXPAND)

        self.apodization_sizer.Add(
            self.apodization_plot_sizer, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.apodization_sizer.AddSpacer(10)

    def update_window_function_plot(self):
        data = self.data
        if self.before_processing == False and self.dimension_index == 0:
            self.npoints = self.nmr_data.number_of_points[0]
            x = np.linspace(
                0,
                (self.nmr_data.number_of_points[0])
                / self.nmr_data.spectral_width[0],
                int(self.nmr_data.number_of_points[0]),
            )
        else:
            if(self.dimension_index==0):
                self.npoints = self.nmr_data.number_of_points[self.dimension_index]
                x = np.linspace(
                    0,
                    (self.nmr_data.number_of_points[self.dimension_index])
                    / self.nmr_data.spectral_width[self.dimension_index],
                    int(self.nmr_data.number_of_points[self.dimension_index]),
            )
            else:
                self.npoints = int(self.nmr_data.number_of_points[self.dimension_index]/2)
                x = np.linspace(
                    0,
                    (self.nmr_data.number_of_points[self.dimension_index]/2)
                    / self.nmr_data.spectral_width[self.dimension_index],
                    int(self.nmr_data.number_of_points[self.dimension_index]/2),
                )

            # if self.dimension_index == 0:
            #     x = np.linspace(
            #         0, (len(data)) / self.nmr_data.spectral_width[0], int(len(data))
            #     )
            # else:
            #     # x = np.linspace(
            #     #     0,
            #     #     (self.nmr_data.number_of_points[self.dimension_index])
            #     #     / self.nmr_data.spectral_width[self.dimension_index],
            #     #     self.nmr_data.number_of_points[self.dimension_index],
            #     # )
            #     x = np.linspace(
            #         0,
            #         (self.nmr_data.number_of_points[self.dimension_index] / 2)
            #         / self.nmr_data.spectral_width[self.dimension_index],
            #         int(self.nmr_data.number_of_points[self.dimension_index] / 2),
            #     )

        try:
            c = float(self.apodization_first_point_textcontrol.GetValue())
            self.apodization_first_point_scaling = c
        except:
            # Give a popout window saying that the values are not valid
            msg = wx.MessageDialog(
                self.app,
                "The value entered for apodization first point scaling is not valid (use 0.5 or 1.0)",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            self.apodization_first_point_textcontrol.SetValue(
                str(self.apodization_first_point_scaling)
            )
            return
        if c != 0.5 and c != 1.0:
            msg = wx.MessageDialog(
                self.app,
                "The value entered for apodization first point scaling is not valid (use 0.5 or 1.0)",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            self.apodization_first_point_textcontrol.SetValue(
                str(self.apodization_first_point_scaling)
            )
            return
        self.apodization_first_point_scaling = c
        if self.apodization_combobox_selection == 1:
            try:
                em = float(self.apodization_line_broadening_textcontrol.GetValue())
            except:
                # Give a popout window saying that the values are not valid
                msg = wx.MessageDialog(
                    self.app,
                    "The values entered are not valid",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_line_broadening_textcontrol.SetValue(
                    str(self.exponential_line_broadening)
                )
                return
            self.exponential_line_broadening = em
            self.line1.set_ydata(
                np.exp(-(np.pi * x * self.exponential_line_broadening))
            )
        elif self.apodization_combobox_selection == 2:
            try:
                g1 = float(self.apodization_g1_textcontrol.GetValue())
                g2 = float(self.apodization_g2_textcontrol.GetValue())
                g3 = float(self.apodization_g3_textcontrol.GetValue())
            except:
                # Give a popout window saying that the values are not valid
                msg = wx.MessageDialog(
                    self.app,
                    "The values entered are not valid",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_g1_textcontrol.SetValue(str(self.g1))
                self.apodization_g2_textcontrol.SetValue(str(self.g2))
                self.apodization_g3_textcontrol.SetValue(str(self.g3))
                return
            # Check to see if g3 is between 0 and 1
            if g3 < 0 or g3 > 1:
                msg = wx.MessageDialog(
                    self.app,
                    "Gaussian shift must be between 0 and 1",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_g3_textcontrol.SetValue(str(self.g3))
                return
            self.g1 = g1
            self.g2 = g2
            self.g3 = g3
            x1 = np.arange(self.npoints)


            g1 = g1/self.nmr_data.spectral_width[self.dimension_index]
            g2 = g2/self.nmr_data.spectral_width[self.dimension_index]

            e = (
                np.pi
                * x1
                * g1
            )
            g = (
                0.6
                * np.pi
                * g2
                * (
                    self.g3
                    * (
                        self.npoints
                        - 1
                    )
                    - x1
                )
            )



            func = np.exp(e - g * g)
            self.line1.set_ydata(func)

            self.apodization_plot_ax.set_xlim(
                0,
                (self.npoints)
                / self.nmr_data.spectral_width[self.dimension_index],
            )

        elif self.apodization_combobox_selection == 3:
            try:
                offset = float(self.apodization_offset_textcontrol.GetValue())
                end = float(self.apodization_end_textcontrol.GetValue())
                power = float(self.apodization_power_textcontrol.GetValue())
                power = int(power)
            except:
                msg = wx.MessageDialog(
                    self.app,
                    "The values entered are not valid",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_offset_textcontrol.SetValue(str(self.offset))
                self.apodization_end_textcontrol.SetValue(str(self.end))
                self.apodization_power_textcontrol.SetValue(str(self.power))
                return
            # Check that offset and end are between 0 and 1
            if offset < 0 or offset > 1:
                msg = wx.MessageDialog(
                    self.app,
                    "Offset values must be between 0 and 1",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_offset_textcontrol.SetValue(str(self.offset))
                return
            if end < 0 or end > 1:
                msg = wx.MessageDialog(
                    self.app,
                    "End values must be between 0 and 1",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_end_textcontrol.SetValue(str(self.end))
                return
            # Check that power is greater than 0
            if power < 0:
                msg = wx.MessageDialog(
                    self.app, "Power must be greater than 0", "Error", wx.OK | wx.ICON_ERROR
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_power_textcontrol.SetValue(str(self.power))
                return
            self.offset = offset
            self.end = end
            self.power = power
            x1 = np.arange(self.npoints)
            func = np.flip(
                np.sin(
                    (np.pi * self.offset + np.pi * (self.end - self.offset) * x1)
                    / (self.npoints-1
                    )
                )
                ** self.power
            )
            # func = (
            #     np.sin(
            #         (np.pi * self.offset + np.pi * (self.end - self.offset) * x)
            #         / (
            #             (
            #                 (
            #                     (
            #                         self.nmr_data.number_of_points[self.dimension_index]
            #                         / 2
            #                     )
            #                     / self.nmr_data.spectral_width[self.dimension_index]
            #                 )
            #             )
            #         )
            #     )
            #     ** self.power
            # )
            self.line1.set_ydata(func)
        elif self.apodization_combobox_selection == 4:
            try:
                a = float(self.apodization_a_textcontrol.GetValue())
                b = float(self.apodization_b_textcontrol.GetValue())
            except:
                msg = wx.MessageDialog(
                    self.app,
                    "The values entered are not valid",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_a_textcontrol.SetValue(str(self.a))
                self.apodization_b_textcontrol.SetValue(str(self.b))
                return
            self.a = a
            self.b = b

            x1 = np.arange(self.npoints)
            t = x1/self.nmr_data.spectral_width[self.dimension_index]
            aq = self.npoints/self.nmr_data.spectral_width[self.dimension_index]

            a = np.pi * a 
            b = -a / (2.0 * b * aq)
            x1 = np.arange(self.npoints)
            func = np.exp(-a * t - (b * (t**2)))
            self.line1.set_ydata(func)
        elif self.apodization_combobox_selection == 5:
            try:
                t1 = float(self.apodization_t1_textcontrol.GetValue())
                t2 = float(self.apodization_t2_textcontrol.GetValue())
            except:
                msg = wx.MessageDialog(
                    self.app,
                    "The values entered are not valid",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_t1_textcontrol.SetValue(str(self.t1))
                self.apodization_t2_textcontrol.SetValue(str(self.t2))
                return
            # Ensure that t1 and t2 are greater than 0
            if t1 < 0 or t2 < 0:
                msg = wx.MessageDialog(
                    self.app,
                    "Ramp up and ramp down points must be greater than 0",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_t1_textcontrol.SetValue(str(self.t1))
                self.apodization_t2_textcontrol.SetValue(str(self.t2))
                return
            # Ensure that t1 + t2 is less than the number of points
            if t1 + t2 > (self.npoints):
                message = (
                    "Ramp up and ramp down points must be less than the number of points ("
                    + str(self.npoints)
                    + ")"
                )
                msg = wx.MessageDialog(self.app, message, "Error", wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_t1_textcontrol.SetValue(str(self.t1))
                self.apodization_t2_textcontrol.SetValue(str(self.t2))
                return
            self.t1 = t1
            self.t2 = t2
            func = np.concatenate(
                (
                    np.linspace(0, 1, int(self.t1)),
                    np.ones(
                        int(self.npoints)
                        - int(self.t1)
                        - int(self.t2)
                    ),
                    np.linspace(1, 0, int(self.t2)),
                )
            )
            self.line1.set_ydata(func)
        elif self.apodization_combobox_selection == 6:
            try:
                loc = float(self.apodization_loc_textcontrol.GetValue())
            except:
                msg = wx.MessageDialog(
                    self.app,
                    "The values entered are not valid",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_loc_textcontrol.SetValue(str(self.loc))
                return
            # Ensure that loc is between 0 and 1
            if self.loc < 0 or self.loc > 1:
                msg = wx.MessageDialog(
                    self.app,
                    "Location of maximum must be between 0 and 1",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_loc_textcontrol.SetValue(str(self.loc))
                return
            self.loc = loc
            func = np.concatenate(
                (
                    np.linspace(
                        0,
                        1,
                        int(
                            self.loc
                            * (self.npoints)
                        ),
                    ),
                    np.linspace(
                        1,
                        0,
                        int(self.npoints)
                        - int(
                            self.loc
                            * int(
                                self.npoints
                            )
                        ),
                    ),
                )
            )
            self.line1.set_ydata(func)

        self.apodization_plot_canvas.draw()
