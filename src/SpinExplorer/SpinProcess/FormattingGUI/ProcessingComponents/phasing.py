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
import numpy as np
import nmrglue as ng
import os
from typing import Dict

# Importing package modules
from SpinExplorer.SpinProcess.FormattingGUI.ProcessingComponents.interactive_phasing import (
    InteractivePhasingFrame,
)
from SpinExplorer.SpinProcess.Processing.process_nmrglue import (
    ProcessNMRGlue,
)


class PhasingDirect:
    def __init__(self, app, nmr_data, parent, info_buttons, apodization_class):
        """
        This class contains all the functions related to phasing
        for the direct dimension. An additional interactive
        phasing module is included in interactive_phasing.py
        """

        self.app = app
        self.nmr_data = nmr_data
        self.info_buttons = info_buttons
        self.apodization_class = apodization_class
        self.parent = parent

        self.set_initial_phasing_variables()
        self.create_phase_correction_sizer(parent)

    def set_initial_phasing_variables(self):
        """
        Initialising relevent phasing variables
        """
        self.phase_correction_checkbox_value = True
        self.p0_total = 0.0
        self.p1_total = 0.0
        self.magnitude_mode_toggle = False

    def create_phase_correction_sizer(self, parent):
        """
        Create a box in the graphical interface for phase
        correction options
        """
        self.phase_correction_box = wx.StaticBox(parent, -1, "Phase Correction")
        self.phase_correction_sizer = wx.StaticBoxSizer(
            self.phase_correction_box, wx.HORIZONTAL
        )
        self.phase_correction_checkbox = wx.CheckBox(
            parent, -1, "Apply phase correction"
        )
        self.phase_correction_checkbox.SetValue(self.phase_correction_checkbox_value)
        self.phase_correction_checkbox.Bind(
            wx.EVT_CHECKBOX, self.on_phase_correction_checkbox
        )
        self.phase_correction_sizer.Add(
            self.phase_correction_checkbox, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.phase_correction_sizer.AddSpacer(10)
        # Have a textcontrol for p0 and p1 values
        self.phase_correction_p0_label = wx.StaticText(
            parent, -1, "Zero order correction (p0):"
        )
        self.phase_correction_sizer.Add(
            self.phase_correction_p0_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.phase_correction_p0_textcontrol = wx.TextCtrl(
            parent, -1, str(self.p0_total), size=(50, 20)
        )
        self.phase_correction_p0_textcontrol.Bind(
            wx.EVT_TEXT, self.on_phase_correction_textcontrol
        )
        self.phase_correction_sizer.Add(
            self.phase_correction_p0_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.phase_correction_sizer.AddSpacer(10)
        self.phase_correction_p1_label = wx.StaticText(
            parent, -1, "First order correction (p1):"
        )
        self.phase_correction_sizer.Add(
            self.phase_correction_p1_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.phase_correction_p1_textcontrol = wx.TextCtrl(
            parent, -1, str(self.p1_total), size=(50, 20)
        )
        self.phase_correction_p1_textcontrol.Bind(
            wx.EVT_TEXT, self.on_phase_correction_textcontrol
        )
        self.phase_correction_sizer.Add(
            self.phase_correction_p1_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.phase_correction_sizer.AddSpacer(10)

        # A button to toggle magnitude mode
        self.magnitude_mode_label = wx.StaticText(parent, -1, "Magnitude Mode:")
        self.phase_correction_sizer.Add(
            self.magnitude_mode_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.magnitude_mode_checkbox = wx.CheckBox(parent, -1)
        self.magnitude_mode_checkbox.SetValue(self.magnitude_mode_toggle)
        self.phase_correction_sizer.Add(
            self.magnitude_mode_checkbox, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.phase_correction_sizer.AddSpacer(10)

        # Have a button for automatic phase correction
        self.phase_correction_auto_button = wx.Button(
            parent, -1, "Interactive Phase Correction"
        )
        self.phase_correction_auto_button.Bind(
            wx.EVT_BUTTON, self.on_phase_correction_interactive
        )
        self.phase_correction_sizer.Add(
            self.phase_correction_auto_button, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.phase_correction_sizer.AddSpacer(10)

        # Have a button showing information on phase correction
        self.phase_correction_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))
        self.phase_correction_info.Bind(
            wx.EVT_BUTTON, self.info_buttons.on_phase_correction_info
        )
        self.phase_correction_sizer.Add(
            self.phase_correction_info, 0, wx.ALIGN_CENTER_VERTICAL
        )
        parent.sizer_1.Add(self.phase_correction_sizer)
        parent.sizer_1.AddSpacer(10)

    def on_phase_correction_checkbox(self, event):
        """
        If the phase correction checkbox is clicked, update
        the current stored checkbox value.
        """
        if self.phase_correction_checkbox.GetValue() == True:
            self.phase_correction_checkbox_value = True
        else:
            self.phase_correction_checkbox_value = False

    def on_phase_correction_textcontrol(self, event):
        """
        When the phase correction values are changed, update the stored values
        to represent the new values. If the first order phase correction (p1)
        is above 45 degrees, set the first point scaling to 1.0.
        """
        self.p0_total = self.phase_correction_p0_textcontrol.GetValue()
        self.p1_total = self.phase_correction_p1_textcontrol.GetValue()
        try:
            if np.abs(float(self.p1_total)) > 45:
                self.apodization_class.apodization_first_point_scaling = 1.0
                self.apodization_class.apodization_first_point_textcontrol.SetValue(
                    "1.0"
                )
            else:
                self.apodization_class.apodization_first_point_scaling = 0.5
                self.apodization_class.apodization_first_point_textcontrol.SetValue(
                    "0.5"
                )
        except:
            pass

    def find_nmr_data_for_phasing(self):
        # Check to see if nmrpipe fid file exists
        # # Check to see what the path of the original frame is
        # if self.parent.parent.original_frame != None:
        #     if self.parent.parent.original_frame.parent.path != "":
        #         path = self.parent.parent.original_frame.parent.path
        #         os.chdir(path)
        # if self.parent.parent.file_parser == True:
        #     os.chdir(self.parent.parent.path)

        if os.path.exists("./test.fid") == False and os.path.exists("fids") == False:
            # Error and exit
            self.error_message = wx.MessageDialog(
                self,
                """No NMRPipe FID file found in the current directory. 
                Please ensure that the file is named test.fid (or folder 
                fids for 3D data) and is in the current directory.""",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            self.error_message.ShowModal()
            return False

        dic, data = ng.pipe.read("test.fid")
        # Perform at fourier transform in the direct dimension
        dic, data = ng.pipe_proc.ft(dic, data, auto=True)
        dic_bruker, dat_bruker = ng.bruker.read("./")
        proc_glue = ProcessNMRGlue(
            self.parent.parent,
            self.parent.parent.tabs,
            self.nmr_data,
            interactive_phasing=True,
        )
        data = proc_glue.remove_digital_filter(dic_bruker, data)

        if len(self.parent.parent.tabs) == 2:
            dic, data = ng.pipe_proc.tp(dic, data)
            dic, data = ng.pipe_proc.ft(dic, data, auto=True)
            dic, data = ng.pipe_proc.tp(dic, data)

        self.nmr_d, self.nmr_spectrum = dic, data

        # # Check to see what the path of the original frame is
        # if self.parent.parent.original_frame != None:
        #     if self.parent.parent.original_frame.parent.cwd != "":
        #         cwd = self.parent.parent.original_frame.parent.cwd
        #         os.chdir(cwd)

        # if self.parent.parent.file_parser == True:
        #     os.chdir(self.parent.parent.cwd)

        # Get the ppm scale
        self.uc = ng.pipe.make_uc(self.nmr_d, self.nmr_spectrum, dim=-1)

        # Get the ppm scale
        self.ppm_scale = self.uc.ppm_scale()

    def on_phase_correction_interactive(self, event):

        v = self.find_nmr_data_for_phasing()
        if v == False:
            return

        # Make a new window with the interactive phase correction
        self.interactive_phase_correction_window = InteractivePhasingFrame(
            self, self.nmr_spectrum, self.ppm_scale, self.nmr_d
        )


class PhasingIndirect:
    def __init__(self, app, nmr_data, parent, info_buttons, apodization_class):
        """
        This class contains all the functions related to phasing
        for the indirect dimensions.
        """

        self.app = app
        self.nmr_data = nmr_data
        self.info_buttons = info_buttons
        self.apodization_class = apodization_class

        self.set_initial_phasing_variables_indirect()
        self.create_phase_correction_sizer_indirect(parent)

    def set_initial_phasing_variables_indirect(self):
        """
        Initialising variables needed for the phasing section of
        the graphical interface for the indirect dimensions.
        """
        self.phasing_indirect_checkbox_value = True
        self.p0_total_indirect = 0.0
        self.p1_total_indirect = 0.0
        self.p0_total_indirect_old = 0.0
        self.p1_total_indirect_old = 0.0
        self.phasing_from_smile = False
        self.f1180 = False

    def create_phase_correction_sizer_indirect(self, parent):
        """
        Create a box for phase correction options in the indirect dimensions.
        """
        self.phase_correction_box_indirect = wx.StaticBox(
            parent, -1, "Phase Correction"
        )
        self.phase_correction_sizer_indirect = wx.StaticBoxSizer(
            self.phase_correction_box_indirect, wx.HORIZONTAL
        )
        self.phase_correction_checkbox_indirect = wx.CheckBox(
            parent, -1, "Apply phase correction"
        )
        self.phase_correction_checkbox_indirect.Bind(
            wx.EVT_CHECKBOX, self.on_phase_correction_checkbox_indirect
        )
        self.phase_correction_checkbox_indirect.SetValue(
            self.phasing_indirect_checkbox_value
        )
        self.phase_correction_sizer_indirect.Add(
            self.phase_correction_checkbox_indirect, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.phase_correction_sizer_indirect.AddSpacer(10)
        # Have a textcontrol for p0 and p1 values
        self.phase_correction_p0_label = wx.StaticText(
            parent, -1, "Zero order correction (p0):"
        )
        self.phase_correction_sizer_indirect.Add(
            self.phase_correction_p0_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.phase_correction_p0_textcontrol_indirect = wx.TextCtrl(
            parent, -1, str(self.p0_total_indirect), size=(50, 20)
        )
        self.phase_correction_p0_textcontrol_indirect.Bind(
            wx.EVT_TEXT, self.on_phase_correction_p0_indirect
        )
        self.phase_correction_sizer_indirect.Add(
            self.phase_correction_p0_textcontrol_indirect, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.phase_correction_sizer_indirect.AddSpacer(10)
        self.phase_correction_p1_label = wx.StaticText(
            parent, -1, "First order correction (p1):"
        )
        self.phase_correction_sizer_indirect.Add(
            self.phase_correction_p1_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.phase_correction_p1_textcontrol_indirect = wx.TextCtrl(
            parent, -1, str(self.p1_total_indirect), size=(50, 20)
        )
        self.phase_correction_p1_textcontrol_indirect.Bind(
            wx.EVT_TEXT, self.on_phase_correction_p1_indirect
        )
        self.phase_correction_sizer_indirect.Add(
            self.phase_correction_p1_textcontrol_indirect, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.phase_correction_sizer_indirect.AddSpacer(10)

        # Have a checkbox for f1180
        self.phase_correction_f1180_button_indirect = wx.CheckBox(parent, -1, "F1180")
        self.phase_correction_f1180_button_indirect.Bind(
            wx.EVT_CHECKBOX, self.on_phase_correction_f1180
        )
        self.phase_correction_f1180_button_indirect.SetValue(self.f1180)
        self.phase_correction_sizer_indirect.Add(
            self.phase_correction_f1180_button_indirect, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.phase_correction_sizer_indirect.AddSpacer(10)

        # Have a button showing information on phase correction
        self.phase_correction_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))
        self.phase_correction_info.Bind(
            wx.EVT_BUTTON, self.info_buttons.on_phase_correction_info_indirect
        )
        self.phase_correction_sizer_indirect.Add(
            self.phase_correction_info, 0, wx.ALIGN_CENTER_VERTICAL
        )
        parent.sizer_1.Add(self.phase_correction_sizer_indirect)
        parent.sizer_1.AddSpacer(10)

    def on_phase_correction_checkbox_indirect(self, event):
        """
        When the phase correction checkbox is clicked, update
        the stored value.
        """
        self.phasing_indirect_checkbox_value = (
            self.phase_correction_checkbox_indirect.GetValue()
        )

    def on_phase_correction_p0_indirect(self, event):
        """
        When the p0 phase correction value is updated, update the
        stored value of p0
        """
        self.p0_total_indirect = (
            self.phase_correction_p0_textcontrol_indirect.GetValue()
        )

    def on_phase_correction_p1_indirect(self, event):
        """
        When the p1 phase correction value is updated, update the
        stored value of p1.
        If p1 is greater than 45 degrees, change the apodization first
        point scaling to 1.0
        """
        self.p1_total_indirect = (
            self.phase_correction_p1_textcontrol_indirect.GetValue()
        )
        try:
            if np.abs(float(self.p1_total_indirect)) > 45:
                self.apodization_class.apodization_first_point_scaling = 1.0
                self.apodization_class.apodization_first_point_textcontrol.SetValue(
                    str(self.apodization_first_point_scaling_indirect)
                )
            else:
                self.apodization_class.apodization_first_point_scaling = 0.5
                self.apodization_class.apodization_first_point_textcontrol_indirect.SetValue(
                    str(self.apodization_class.apodization_first_point_scaling)
                )
        except:
            pass

    def on_phase_correction_f1180(self, event):
        """
        When f1180 is checked, uodate the indirect phasing to
        p0=-90, p1=180.
        Also update the first point scaling factor to 1.0 in the
        apodization box.
        """
        if self.phase_correction_f1180_button_indirect.GetValue() == True:
            self.c_old = self.apodization_class.apodization_first_point_scaling
            self.p0_total_indirect_old = self.p0_total_indirect
            self.p1_total_indirect_old = self.p1_total_indirect
            # Apply -90 p0 and 180 p1 to the phase correction textcontrols
            self.p0_total_indirect = -90.0
            self.p1_total_indirect = 180.0
            self.phase_correction_p0_textcontrol_indirect.SetValue(
                str(self.p0_total_indirect)
            )
            self.phase_correction_p1_textcontrol_indirect.SetValue(
                str(self.p1_total_indirect)
            )
            # Disable the phase correction textcontrols
            self.phase_correction_p0_textcontrol_indirect.Disable()
            self.phase_correction_p1_textcontrol_indirect.Disable()
            self.apodization_class.apodization_first_point_scaling_indirect = 1.0
            self.apodization_class.apodization_first_point_textcontrol.SetValue(
                str(self.apodization_class.apodization_first_point_scaling)
            )
        else:
            self.p0_total_indirect = self.p0_total_indirect_old
            self.p1_total_indirect = self.p1_total_indirect_old
            self.phase_correction_p0_textcontrol_indirect.SetValue(
                str(self.p0_total_indirect)
            )
            self.phase_correction_p1_textcontrol_indirect.SetValue(
                str(self.p1_total_indirect)
            )
            self.phase_correction_p0_textcontrol_indirect.Enable()
            self.phase_correction_p1_textcontrol_indirect.Enable()
            try:
                self.apodization_class.apodization_first_point_scaling = self.c_old
            except:
                self.apodization_class.apodization_first_point_scaling = 0.5
            self.apodization_class.apodization_first_point_textcontrol.SetValue(
                str(self.apodization_class.apodization_first_point_scaling)
            )
