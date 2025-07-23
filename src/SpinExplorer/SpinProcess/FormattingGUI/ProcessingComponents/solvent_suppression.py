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
from typing import Dict


class SolventSuppression:
    def __init__(self, app, nmr_data, parent, info_buttons):
        """
        This class contains all the functions relevent for digital
        solvent suppression of the direct dimension
        """

        self.app = app
        self.nmr_data = nmr_data
        self.info_buttons = info_buttons
        self.parent = parent

        self.set_initial_solvent_suppression_variables()
        self.create_solvent_suppression_sizer(self.parent)

    def set_initial_solvent_suppression_variables(self):
        """
        Setting the initial solvent suppression processing parameters in the
        graphical interface to default values
        """
        if (
            self.nmr_data.axislabels[0] == "H1"
            or self.nmr_data.axislabels[0] == "1H"
            or self.nmr_data.axislabels[0] == "H"
        ):
            self.direct_solvent_suppression = True
        else:
            self.direct_solvent_suppression = False

        self.include_direct_linear_prediction = False

        self.solvent_suppression_filter_selection = 0
        self.solvent_suppression_lowpass_shape_selection = 0
        self.solvent_suppression_filter_length = 32
        self.solvent_suppression_polynomial_order = 2
        self.solvent_suppression_spline_noise = 1.0
        self.solvent_suppression_spline_smoothfactor = 1.1

    def load_solvent_suppression_variables(self, parameter_dictionary: Dict):
        """
        Function to load saved solvent suppression processing parameters
        into the graphical interface
        """
        pass

    def create_solvent_suppression_sizer(self, parent):
        # Create a box for solvent suppression options
        self.solvent_suppression_box = wx.StaticBox(parent, -1, "Solvent Suppression")
        self.solvent_suppression_sizer = wx.StaticBoxSizer(
            self.solvent_suppression_box, wx.HORIZONTAL
        )
        self.solvent_suppression_checkbox = wx.CheckBox(
            parent, -1, "Apply solvent suppression"
        )
        self.solvent_suppression_checkbox.SetValue(self.direct_solvent_suppression)
        self.solvent_suppression_sizer.Add(
            self.solvent_suppression_checkbox, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.solvent_suppression_sizer.AddSpacer(10)
        self.solvent_suppression_extra_options = wx.Button(
            parent, -1, "Advanced Options"
        )
        self.solvent_suppression_sizer.Add(
            self.solvent_suppression_extra_options, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.solvent_suppression_sizer.AddSpacer(10)
        self.solvent_suppression_extra_options.Bind(
            wx.EVT_BUTTON, self.solvent_suppression_extra_options_click
        )
        # Have a button showing information on solvent suppression
        self.solvent_suppression_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))
        self.solvent_suppression_info.Bind(
            wx.EVT_BUTTON, self.info_buttons.on_solvent_suppression_info
        )
        self.solvent_suppression_sizer.Add(
            self.solvent_suppression_info, 0, wx.ALIGN_CENTER_VERTICAL
        )
        parent.sizer_1.Add(self.solvent_suppression_sizer)
        parent.sizer_1.AddSpacer(10)

    def solvent_suppression_extra_options_click(self, event):
        """
        Creating a popup window with options for solvent suppression.
        Options include:
        - Low-pass filter
        - Spline filter (nmrPipe)
        - Polynomial filter (nmrPipe)

        For a low-pass filter, there are also the following options:
        Filter shape: Boxcar, Sine, Sine Squared
        Filter size: integer
        """

        self.solvent_suppression_extra_options_window = wx.Frame(
            self.parent, -1, "Solvent Suppression Advanced Options", size=(400, 300)
        )

        self.solvent_suppression_extra_options_window_sizer = wx.BoxSizer(wx.VERTICAL)
        self.solvent_suppression_extra_options_window.SetSizer(
            self.solvent_suppression_extra_options_window_sizer
        )

        # Create a sizer for the solvent suppression options
        self.solvent_suppression_extra_options_sizer = wx.BoxSizer(wx.VERTICAL)
        self.solvent_suppression_extra_options_window_sizer.Add(
            self.solvent_suppression_extra_options_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.solvent_suppression_extra_options_window_sizer.AddSpacer(10)

        # Have a radio button for low-pass filter, spline, or polynomial
        self.radio_box = wx.RadioBox(
            self.solvent_suppression_extra_options_window,
            -1,
            "Solvent Suppression Method",
            choices=["Low-pass filter", "Spline", "Polynomial"],
            style=wx.RA_SPECIFY_ROWS,
        )
        self.radio_box.Bind(wx.EVT_RADIOBOX, self.OnSolventSuppressionChoice)
        self.radio_box.SetSelection(self.solvent_suppression_filter_selection)
        self.solvent_suppression_extra_options_sizer.Add(
            self.radio_box, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.solvent_suppression_extra_options_sizer.AddSpacer(10)

        # Have a radiobox for the low-pass filter shape
        self.lowpass_shape_radio_box = wx.RadioBox(
            self.solvent_suppression_extra_options_window,
            -1,
            "Low-pass filter shape",
            choices=["Boxcar", "Sine", "Sine Squared"],
            style=wx.RA_SPECIFY_ROWS,
        )
        self.lowpass_shape_radio_box.Bind(wx.EVT_RADIOBOX, self.OnFilterShapeChoice)
        self.lowpass_shape_radio_box.SetSelection(
            self.solvent_suppression_lowpass_shape_selection
        )
        self.solvent_suppression_extra_options_sizer.Add(
            self.lowpass_shape_radio_box, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.solvent_suppression_extra_options_sizer.AddSpacer(10)

        # Have a textcontrol box for the filter size
        self.filter_size_label = wx.StaticBox(
            self.solvent_suppression_extra_options_window, -1, label="Filter Size:"
        )
        self.filter_size_sizer = wx.StaticBoxSizer(self.filter_size_label)
        self.filter_size_box = wx.TextCtrl(
            self.solvent_suppression_extra_options_window,
            -1,
            str(int(self.solvent_suppression_filter_length)),
        )
        self.filter_size_box.Bind(wx.EVT_KEY_DOWN, self.OnFilterSizeKey)
        self.filter_size_sizer.Add(self.filter_size_box)

        self.solvent_suppression_extra_options_sizer.Add(
            self.filter_size_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.solvent_suppression_extra_options_sizer.AddSpacer(10)

        self.solvent_suppression_extra_options_window.SetSizer(
            self.solvent_suppression_extra_options_window_sizer
        )
        self.solvent_suppression_extra_options_window.Show()

    def OnFilterSizeKey(self, event):
        """
        If the user presses enter, update solvent suppression
        filter length. Includes a check that the filter size
        is an integer.
        """
        keycode = event.GetKeyCode()
        if keycode == wx.WXK_RETURN:
            try:
                self.solvent_suppression_filter_length = int(
                    self.filter_size_box.GetValue()
                )
            except:
                # Providing a popout to say that the filter size could not be converted to an integer, please try again
                dlg = wx.MessageBox(
                    self,
                    "The filter size could not be converted to an integer, please try again.",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                dlg.ShowModal()
                dlg.Destroy()
        event.Skip()

    def OnSolventSuppressionChoice(self, event):
        """
        Change the type of digital solvent suppression filter
        """
        self.solvent_suppression_filter_selection = int(self.radio_box.GetSelection())

    def OnFilterShapeChoice(self, event):
        """
        Change the type of digital solvent suppression filter
        """
        self.solvent_suppression_lowpass_shape_selection = int(
            self.lowpass_shape_radio_box.GetSelection()
        )
