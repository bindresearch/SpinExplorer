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

print("-------------------------------------------------------------")
print("                        SpinConverter                        ")
print("-------------------------------------------------------------")
print("                (version 1.4) 20th June 2025                 ")
print(" (c) 2025 James Eaton, Andrew Baldwin (University of Oxford) ")
print("                        MIT License                          ")
print("-------------------------------------------------------------")
print("            Converting NMR data to nmrPipe format            ")
print("-------------------------------------------------------------")
print(" Documentation at:")
print(" https://github.com/james-eaton-1/SpinExplorer")
print("-------------------------------------------------------------")
print("")


# Import relevant external modules
import sys
import wx
import numpy as np
import nmrglue as ng
import subprocess
import os
import darkdetect
import warnings

# Importing internal classes
from SpinExplorer.SpinConverter.FindingParameters.parameters import FindingParameters
from SpinExplorer.SpinConverter.FormattingGUI.bruker_formatting import (
    FormatParametersBruker,
)
from SpinExplorer.SpinConverter.FormattingGUI.varian_formatting import (
    FormatParametersVarian,
)
from SpinExplorer.SpinConverter.FormattingGUI.shared_formatting import (
    SharedFormatting,
)
from SpinExplorer.SpinConverter.StoringParameters.save_parameters import Save_json
from SpinExplorer.SpinConverter.StoringParameters.read_parameters import Read_json
from SpinExplorer.SpinConverter.Conversion.convert_pipe import Convert_pipe
from SpinExplorer.SpinConverter.Conversion.convert_nmrglue import Convert_nmrglue


warnings.simplefilter("ignore", UserWarning)


# Check to see if using mac, linux or windows
if sys.platform == "darwin":
    platform = "mac"
elif sys.platform == "linux":
    platform = "linux"
else:
    platform = "windows"

# See if the nmrPipe command works, if not set the platform to windows
if platform == "mac" or platform == "linux":
    try:
        p = subprocess.Popen(
            "nmrPipe", stdout=subprocess.DEVNULL, stderr=subprocess.PIPE
        )
        out, err = p.communicate()
        if "NMRPipe System Version" in str(err):
            platform = platform
        else:
            platform = "windows"

    except:
        platform = "windows"


# James Eaton, 10/06/2025, University of Oxford
# This program is designed to allow the user to convert NMR data from Bruker/Varian into NMRPipe format so it can be viewed using
# SpinView.py. It is designed to be used with the SpinProcess.py program used to process the converted nmrPipe FID to produce
# an NMR spectrum. These spectra can then be viewed using SpinView.py, a GUI for viewing NMR data.


class MyApp(wx.Frame):
    def __init__(self):
        """
        This class creates the GUI showing the found parameters with scope
        for changing the parameters.
        """

        # Get the monitor size and set the window size to 85% of the monitor size
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 0.6 * self.monitorWidth
        self.height = 0.6 * self.monitorHeight
        self.app_frame = wx.Frame.__init__(
            self,
            None,
            wx.ID_ANY,
            "SpinConverter",
            wx.DefaultPosition,
            size=(int(self.width), int(self.height)),
        )

        self.file_parser = False

        # Initialise the NMR data and parameter class
        self.nmrdata = FindingParameters()

        # Creating a canvas and formatting the app on it
        self.create_canvas()
        self.format_app()

        # Reading previously saved parameters if present
        read = Read_json(self.nmrdata.params, self.nmrdata, self)

        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self.Show()
        self.Centre()

    def OnClose(self, event):
        """
        Ensuring the application is closed after pressing close
        """
        self.Destroy()
        sys.exit()

    def create_canvas(self) -> None:
        """
        Creating a canvas for the application
        """
        # Create the main sizer
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)

        if darkdetect.isDark() == True and platform != "windows":
            self.SetBackgroundColour((53, 53, 53, 255))
        else:
            self.SetBackgroundColour("White")

    def create_sizers(self) -> None:
        self.parameters_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.menu_bar = wx.BoxSizer(wx.VERTICAL)
        self.extra_sizers = wx.BoxSizer(wx.VERTICAL)
        self.parameters_sizer.Add(self.menu_bar)
        self.extra_boxes_total = wx.BoxSizer(wx.VERTICAL)
        self.extra_boxes_0 = wx.BoxSizer(wx.HORIZONTAL)
        self.extra_boxes_0_total = wx.BoxSizer(wx.HORIZONTAL)
        self.extra_boxes = wx.BoxSizer(wx.HORIZONTAL)
        self.extra_boxes_total_1 = wx.BoxSizer(wx.HORIZONTAL)
        self.extra_boxes.AddSpacer(10)
        self.extra_boxes_0.AddSpacer(20)
        self.extra_boxes_0_total.Add(self.extra_boxes_0, 0, wx.CENTER)
        self.extra_boxes_total.Add(self.extra_boxes_0_total, 0, wx.CENTER)
        self.extra_boxes_total.AddSpacer(20)
        self.extra_boxes_total_1.Add(self.extra_boxes, 0, wx.CENTER)
        self.extra_boxes_total.Add(self.extra_boxes_total_1, 0, wx.CENTER)

    def format_app(self) -> None:
        """
        Adding the TextControl sizers and buttons to the app and
        populating them with parameters.
        """

        self.create_sizers()

        if self.nmrdata.spectrometer == "Bruker":
            self.format = FormatParametersBruker(
                self, self.nmrdata.params, self.nmrdata
            )
            self.shared_format = SharedFormatting(
                self, self.nmrdata.params, self.nmrdata
            )
            self.format.input_sizes_bruker()
            self.format.input_acquisition_modes_bruker()
            self.format.input_sweep_widths_bruker()
            self.format.get_nuclei_frequency_bruker()
            self.format.get_nuclei_labels_bruker()
            self.format.get_carrier_frequencies_bruker()
            if len(self.format.N_complex_boxes) > 1:
                self.shared_format.acquisition_2D_mode_combo_box()
            self.shared_format.create_temperature_box()
            self.format.create_bruker_digital_filter_box()
            self.shared_format.create_conversion_box()
            self.format.create_other_options_box()
            self.shared_format.create_intensity_scaling_box()
            if self.nmrdata.params.size_indirect != []:
                self.shared_format.find_nus_file()
                self.shared_format.input_NUS_list_box()
            else:
                self.shared_format.include_NUS = False
        elif self.nmrdata.spectrometer == "Varian":
            self.format = FormatParametersVarian(self)
            self.shared_format = SharedFormatting(
                self, self.nmrdata.params, self.nmrdata
            )
            self.format.input_sizes_varian()
            self.format.input_acquisition_modes_varian()
            self.format.input_sweep_widths_varian()
            self.format.get_nuclei_frequency_varian()
            self.format.get_nuclei_labels_varian()
            self.format.get_carrier_frequencies_varian()
            if len(self.shared_format.N_complex_boxes) > 1:
                self.shared_format.acquisition_2D_mode_combo_box()
            self.shared_format.create_temperature_box()
            self.shared_format.create_conversion_box()
            self.shared_format.create_intensity_scaling_box()
            if self.nmrdata.phase != False or self.nmrdata.phase2 != False:
                self.shared_format.find_nus_file()
                self.shared_format.input_NUS_list_box()
            else:
                self.shared_format.include_NUS = False

        self.main_sizer.Add(self.parameters_sizer, 0, wx.CENTER)
        self.main_sizer.Add(self.extra_sizers, 0, wx.CENTER)

        self.SetSizerAndFit(self.main_sizer)

        # Get the width and height of the main_sizer
        self.width, self.height = self.main_sizer.GetSize()
        self.SetSize((int(self.width * 1.25), int(self.height * 1.25)))
        self.Centre()

    def on_save_parameters(self, event) -> None:
        """
        Saving the current SpinConverter parameters to parameters.json
        """
        save = Save_json(self.nmrdata.params, self.nmrdata, self)

    def on_convert_pipe(self, event) -> None:
        """
        Checking to see that nmrPipe is installed and then performing
        nmrPipe conversion.
        """
        if platform == "windows":
            # Outputting a message saying that nmrPipe conversion is not possible on windows
            dlg = wx.MessageDialog(
                self,
                "It seems like nmrPipe is not installed. Please use the nmrglue convert button instead.",
                "Warning",
                wx.OK | wx.ICON_WARNING,
            )
            self.Raise()
            self.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            return
        else:
            pipe_conversion = Convert_pipe(self, self.nmrdata.params, self.nmrdata)

    def on_convert_glue(self, event) -> None:
        glue_conversion = Convert_nmrglue(self, self.nmrdata.params, self.nmrdata)


def main():
    app = wx.App()
    frame = MyApp()
    app.MainLoop()


if __name__ == "__main__":
    main()
