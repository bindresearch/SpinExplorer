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
import os
import pkgutil
import importlib
import inspect

# Importing relevant SpinProcess modules
import SpinExplorer.SpinProcess.FormattingGUI.ProcessingComponents as p


class DirectDimensionFrame(wx.Panel):

    def __init__(self, app, parent, info_buttons):
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 0.7 * self.monitorWidth
        self.height = 0.75 * self.monitorHeight
        self.parent = parent
        wx.Panel.__init__(self, parent, id=wx.ID_ANY, size=(self.width, self.height))
        # Create panel for processing dimension 1 of the data
        self.nmr_data = parent.nmr_data
        self.info_buttons = info_buttons
        self.app = app
        self.create_menu_bar()

    def load_variables(self):

        # See if NMR processing file (nmrproc.com) can be found, if it can try to load the variables from it
        if os.path.exists("processing_parameters.txt"):
            found_nmrproc_com = True
        else:
            found_nmrproc_com = False

        self.parent.load_variables = False
        if found_nmrproc_com == False:
            pass
        else:
            # Ask the user if they want to load the variables from the nmrproc.com file
            dlg = wx.MessageDialog(
                self,
                """A file containing NMR processing parameters has been found (processing_parameters.txt). 
                Do you want to load the variables from it?""",
                "Warning",
                wx.YES_NO | wx.ICON_WARNING,
            )
            self.Raise()
            self.SetFocus()
            result = dlg.ShowModal()
            if result == wx.ID_YES:
                try:
                    self.parent.load_variables = True
                    self.load_variables_from_nmrproc_com_1D()
                except:
                    pass

            else:
                pass

    def create_menu_bar(self):
        """
        Creating the main panel for direct dimension processing.
        """
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_sizer)

        self.sizer_1 = wx.BoxSizer(wx.VERTICAL)
        self.sizer_1.AddSpacer(10)

        # Create all the sizers
        self.solvent_suppression = p.SolventSuppression(
            self.app, self.nmr_data, self, self.info_buttons
        )
        self.linear_prediction = p.LinearPrediction(
            self.app, self.nmr_data, self, self.info_buttons
        )
        self.apodization = p.Apodization(
            self.app,
            self.nmr_data,
            self,
            self.info_buttons,
            [self.solvent_suppression, self.linear_prediction],
            0,
        )
        self.zero_filling = p.ZeroFilling(
            self.app,
            self.nmr_data,
            self,
            self.info_buttons,
            [self.solvent_suppression, self.linear_prediction, self.apodization],
            0,
        )
        self.fourier_transform = p.FourierTransform(
            self.app,
            self.nmr_data,
            self,
            self.info_buttons,
        )
        self.phasing = p.PhasingDirect(
            self.app, self.nmr_data, self, self.info_buttons, self.apodization
        )
        self.extraction = p.Extraction(
            self.app,
            self.nmr_data,
            self,
            self.info_buttons,
        )
        self.baseline_correction = p.BaselineCorrection(
            self.app,
            self.nmr_data,
            self,
            self.info_buttons,
        )

        self.main_sizer.Add(self.sizer_1, 0, wx.EXPAND)

        self.SetSizerAndFit(self.main_sizer)
        self.Layout()

        # Get the size of the main sizer and set the window size to 1.05 times the size of the main sizer
        self.width, self.height = self.main_sizer.GetSize()
        self.parent.parent.change_frame_size(
            int(self.width * 1.05), int(self.height * 1.25)
        )

    def refresh_menu_bar(self):
        """
        Refreshing the main panel for direct dimension processing.
        Will be activated when the comboboxes of linear prediction,
        apodization or zero filling are altered.
        """
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_sizer)
        self.sizer_1 = wx.BoxSizer(wx.VERTICAL)
        self.sizer_1.AddSpacer(10)

        self.solvent_suppression.create_solvent_suppression_sizer(self)
        self.linear_prediction.create_linear_prediction_sizer(self)
        self.apodization.create_apodization_sizer(self)
        self.zero_filling.create_zero_filling_sizer(self)
        self.fourier_transform.create_fourier_transform_sizer(self)
        self.phasing.create_phase_correction_sizer(self)
        self.extraction.create_extraction_sizer(self)
        self.baseline_correction.create_baseline_correction_sizer(self)

        self.main_sizer.Add(self.sizer_1, 0, wx.EXPAND)
        self.SetSizerAndFit(self.main_sizer)
        self.Layout()

        # Get the size of the main sizer and set the window size to 1.05 times the size of the main sizer
        self.width, self.height = self.main_sizer.GetSize()
        self.parent.parent.change_frame_size(
            int(self.width * 1.05), int(self.height * 1.25)
        )


class IndirectDimensionFrame(wx.Panel):

    def __init__(self, app, parent, info_buttons, direct_dimension_frame):
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 0.7 * self.monitorWidth
        self.height = 0.75 * self.monitorHeight
        self.parent = parent
        wx.Panel.__init__(self, parent, id=wx.ID_ANY, size=(self.width, self.height))
        # Create panel for processing dimension 1 of the data
        self.nmr_data = parent.nmr_data
        self.info_buttons = info_buttons
        self.app = app

        self.direct_dimension_frame = direct_dimension_frame

        self.create_menu_bar_indirect()

    def create_menu_bar_indirect(self):
        """
        Creating the main panel for indirect dimension processing.
        """
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_sizer)

        # Create a sizer for the processing options for the first dimension
        self.sizer_1 = wx.BoxSizer(wx.VERTICAL)
        self.sizer_1.AddSpacer(10)

        # Add all the processing modules

        self.linear_prediction = p.NonUniformSampling(
            self.app, self.nmr_data, self, self.info_buttons
        )
        self.apodization = p.Apodization(
            self.app,
            self.nmr_data,
            self,
            self.info_buttons,
            [self.linear_prediction],
            1,
        )
        self.linear_prediction.apodization_class = self.apodization

        self.zero_filling = p.ZeroFilling(
            self.app,
            self.nmr_data,
            self,
            self.info_buttons,
            [self.linear_prediction, self.apodization],
            1,
        )
        self.fourier_transform = p.FourierTransform(
            self.app,
            self.nmr_data,
            self,
            self.info_buttons,
        )
        self.phasing = p.PhasingIndirect(
            self.app, self.nmr_data, self, self.info_buttons, self.apodization
        )
        self.extraction = p.Extraction(
            self.app,
            self.nmr_data,
            self,
            self.info_buttons,
        )
        self.baseline_correction = p.BaselineCorrection(
            self.app,
            self.nmr_data,
            self,
            self.info_buttons,
        )

        self.main_sizer.Add(self.sizer_1, 0, wx.EXPAND)

        self.SetSizerAndFit(self.main_sizer)
        self.Layout()

        # Get the size of the main sizer and set the window size to 1.05 times the size of the main sizer
        self.width, self.height = self.main_sizer.GetSize()
        self.parent.parent.change_frame_size(
            int(self.width * 1.05), int(self.height * 1.25)
        )

    def refresh_menu_bar(self):
        """
        Refreshing the main panel for direct dimension processing.
        Will be activated when the comboboxes of linear prediction,
        apodization or zero filling are altered.
        """

        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_sizer)
        self.sizer_1 = wx.BoxSizer(wx.VERTICAL)
        self.sizer_1.AddSpacer(10)

        self.linear_prediction.create_linear_prediction_sizer_indirect(self)
        self.apodization.create_apodization_sizer(self)
        self.zero_filling.create_zero_filling_sizer(self)
        self.fourier_transform.create_fourier_transform_sizer(self)
        self.phasing.create_phase_correction_sizer_indirect(self)
        self.extraction.create_extraction_sizer(self)
        self.baseline_correction.create_baseline_correction_sizer(self)

        self.main_sizer.Add(self.sizer_1, 0, wx.EXPAND)
        self.SetSizerAndFit(self.main_sizer)
        self.Layout()

        # Get the size of the main sizer and set the window size to 1.05 times the size of the main sizer
        self.width, self.height = self.main_sizer.GetSize()
        self.parent.parent.change_frame_size(
            int(self.width * 1.05), int(self.height * 1.25)
        )
