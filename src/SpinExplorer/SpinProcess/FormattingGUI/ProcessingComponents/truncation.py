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

from .dimension_size import update_dimension_size


class Truncation:

    def __init__(self, app, nmr_data, parent, info_buttons):
        """
        This class contains all the functions related to truncating the
        data, where only the first points of a dimension are processed
        and the points recorded after them are discarded.
        """

        self.app = app
        self.nmr_data = nmr_data
        self.parent = parent
        self.info_buttons = info_buttons

        self.set_initial_truncation_variables()
        self.create_truncation_sizer(parent)

    def set_initial_truncation_variables(self):
        """
        Initialising parameters necessary for the truncation section of
        the GUI. By default all of the recorded points are processed.
        """
        self.truncation_checkbox_value = False

        # None means all of the recorded points, so that the number shown
        # follows the size of the data until the user sets a value
        self.truncation_points = None

    def find_recorded_points(self) -> int:
        """
        The number of complex points recorded in this dimension.
        """
        dimension_size = getattr(self.parent, "dimension_size", None)
        if dimension_size == None:
            return 0

        return dimension_size.find_current_size()

    def find_truncation_points(self) -> int:
        """
        The number of complex points which will be processed.
        """
        if self.truncation_points == None:
            return self.find_recorded_points()

        return self.truncation_points

    def create_truncation_sizer(self, parent):
        """
        A box for a user to choose to process fewer points than were
        recorded in this dimension.
        """
        self.truncation_box = wx.StaticBox(parent, -1, "Truncation")
        self.truncation_sizer = wx.StaticBoxSizer(self.truncation_box, wx.HORIZONTAL)

        self.truncation_checkbox = wx.CheckBox(
            self.truncation_box, -1, "Process fewer points than recorded"
        )
        self.truncation_checkbox.SetValue(self.truncation_checkbox_value)
        self.truncation_checkbox.Bind(wx.EVT_CHECKBOX, self.on_truncation_checkbox)
        self.truncation_sizer.Add(
            self.truncation_checkbox, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.truncation_sizer.AddSpacer(10)

        self.truncation_points_label = wx.StaticText(
            self.truncation_box, -1, "Number of complex points:"
        )
        self.truncation_sizer.Add(
            self.truncation_points_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.truncation_sizer.AddSpacer(5)

        self.truncation_points_textcontrol = wx.TextCtrl(
            self.truncation_box, -1, str(self.find_truncation_points()), size=(60, 20)
        )
        self.truncation_points_textcontrol.Bind(
            wx.EVT_TEXT, self.on_truncation_textcontrol
        )
        self.truncation_sizer.Add(
            self.truncation_points_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.truncation_sizer.AddSpacer(10)

        self.truncation_recorded_label = wx.StaticText(
            self.truncation_box, -1, self.find_recorded_text()
        )
        self.truncation_sizer.Add(
            self.truncation_recorded_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.truncation_sizer.AddSpacer(10)

        # Have a button showing information on truncation
        self.truncation_info = wx.Button(self.truncation_box, -1, "\u24d8", size=(25, 32))
        self.truncation_info.Bind(wx.EVT_BUTTON, self.info_buttons.on_truncation_info)
        self.truncation_sizer.Add(self.truncation_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.truncation_sizer.AddSpacer(10)

        parent.sizer_1.Add(self.truncation_sizer)
        parent.sizer_1.AddSpacer(10)

    def find_recorded_text(self) -> str:
        """
        The text reminding the user how many points were recorded.
        """
        return "of {} recorded".format(self.find_recorded_points())

    def on_truncation_checkbox(self, event):
        """
        When the truncation checkbox is clicked, update the stored value
        """
        self.truncation_checkbox_value = self.truncation_checkbox.GetValue()

        update_dimension_size(self.parent)

    def on_truncation_textcontrol(self, event):
        """
        When a user changes the number of points to process, check that the
        value is a valid number of points and update the stored value.
        """
        value = self.truncation_points_textcontrol.GetValue()

        if value.strip() == "":
            # The box has been emptied, all of the recorded points are used
            # until a new value is typed
            self.truncation_points = None
            update_dimension_size(self.parent)
            return

        try:
            points = int(value)
        except ValueError:
            return

        if points < 1:
            return

        self.truncation_points = points

        update_dimension_size(self.parent)
