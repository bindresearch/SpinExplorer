#!/usr/bin/env python3


"""MIT License

Copyright (c) 2025 James Eaton, Andrew Baldwin

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
import nmrglue as ng
import os
import warnings
from SpinExplorer.SpinConverter.FindingParameters.varian_parameters import (
    ParameterExtractorVarian,
)
from SpinExplorer.SpinConverter.FindingParameters.bruker_parameters import (
    ParameterExtractorBruker,
)

warnings.simplefilter("ignore", UserWarning)


class FindingParameters:

    def __init__(self) -> None:
        """
        Initialising the converter class. The class reads the spectrometer
        parameter files to obtain relevant spectra for nmrpipe conversion
        scripts.
        """

        # Creating a hidden frame to be used as a parent for popout messages
        self.tempframe = wx.Frame(None, title="Temporary Parent", size=(1, 1))
        self.tempframe.Hide()  # Hide the frame since we don't need it to be visible

        # Get the NMR data and parameters
        self.find_nmr_files()
        self.read_nmr_data()
        self.find_parameters()

    def find_nmr_files(self) -> None:
        """
        Finding the relevant NMR parameter and data files in the current
        directory.
        """
        self.parameter_file = ""
        if "acqus" in os.listdir("."):
            self.spectrometer = "Bruker"
            self.parameter_file = "acqus"
        elif "procpar" in os.listdir("."):
            self.spectrometer = "Varian"
            self.parameter_file = "procpar"
        else:
            # Try to find acqu file
            if "acqu" in os.listdir("."):
                self.spectrometer = "Bruker"
                self.parameter_file = "acqu"
                # Inform the user that an acqu parameter file but not an acqus parameter file has been found, would you like to continue?
                dlg = wx.MessageDialog(
                    self.tempframe,
                    "An acqu parameter file but not an acqus parameter file has been found, would you like to continue?",
                    "Continue",
                    wx.YES_NO | wx.ICON_ERROR,
                )
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                if dlg.ShowModal() == wx.ID_YES:
                    dlg.Destroy()
                else:
                    dlg.Destroy()
                    exit()
            else:
                # Give a popout error message saying that there are no bruker NMR files in the current directory, but found an acqus file
                dlg = wx.MessageDialog(
                    self.tempframe,
                    "No Bruker (acqus) or Varian (procpar) NMR files found in the current directory. Please check the current directory and try again.",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                dlg.ShowModal()
                dlg.Destroy()
                exit()
        if self.spectrometer == "Bruker":
            self.files = []
            for file in os.listdir("."):
                if file.endswith("ser") or file == "fid":
                    self.files.append(file)
            if len(self.files) == 0:
                # Give a popout error message saying that there are no bruker NMR files in the current directory, but found an acqus file
                dlg = wx.MessageDialog(
                    self.tempframe,
                    "No Bruker NMR files found in the current directory."
                    "Please check the current directory and try again.",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                dlg.ShowModal()
                dlg.Destroy()
                exit()

            elif len(self.files) == 2:
                # Give a popout error message saying that there are two NMR files in the current directory
                dlg = wx.MessageDialog(
                    self.tempframe,
                    "An acqus file was found and two NMR files (fid and ser)"
                    "found in the current directory. Please check the current"
                    "directory and try again.",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                dlg.ShowModal()
                dlg.Destroy()
                exit()
        elif self.spectrometer == "Varian":
            self.files = []
            for file in os.listdir("."):
                if file == "fid" or file == "origfid":
                    self.files.append(file)
            if len(self.files) == 0:
                # Give a popout error message saying that there are no bruker NMR files in the current directory, but found a procpar file
                dlg = wx.MessageDialog(
                    self.tempframe,
                    "No Varian NMR files found in the current directory, but a procpar file was found. Please press no to exit, or press yes to manually select the raw FID file.",
                    "Error",
                    wx.YES_NO | wx.ICON_ERROR,
                )
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                if dlg.ShowModal() == wx.ID_YES:
                    dlg = wx.FileDialog(
                        self.tempframe,
                        "Select the raw FID file",
                        wildcard="",
                        style=wx.FD_OPEN,
                    )
                    dlg.SetDirectory(os.getcwd())
                    if dlg.ShowModal() == wx.ID_OK:
                        self.files = [dlg.GetPath()]
                    else:
                        dlg.Destroy()
                        exit()
                else:
                    dlg.Destroy()
                    exit()

    def read_nmr_data(self) -> None:
        """
        Reading the NMR data file to obtain data dimensions
        """
        if self.spectrometer == "Bruker":
            # if pdata directory exists, if it is empty then change its name to pdata_original
            if "pdata" in os.listdir("."):
                if os.listdir("pdata") == []:
                    os.rename("pdata", "pdata_original")
            try:
                self.nmr_dic, self.nmr_data = ng.bruker.read(
                    dir="./", bin_file=self.files[0]
                )
            except:
                dlg = wx.MessageDialog(
                    self.tempframe,
                    "Error: Unable to read the Bruker NMR data. Please check the current directory and try again.",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                dlg.ShowModal()
                dlg.Destroy()
                exit()
        if self.spectrometer == "Varian":
            try:
                self.nmr_dic, self.nmr_data = ng.varian.read(
                    dir="./", fid_file=self.files[0]
                )
            except:
                dlg = wx.MessageDialog(
                    self.tempframe,
                    "Error: Unable to read the Varian NMR data. Please check the current directory and try again.",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                dlg.ShowModal()
                dlg.Destroy()
                exit()

        self.data_dimensions = len(self.nmr_data.shape)

    def find_parameters(self) -> None:
        """
        Finding the relevant parameters from the spectrometer parameter files
        """
        if self.spectrometer == "Bruker":
            # Initialising the class to extract relevant Bruker parameters
            self.params = ParameterExtractorBruker(self)
            # Extracting spectrum dimensions
            self.params.find_size_bruker()
            # Determine if any axis is a pseudo (non-complex) axis
            self.params.find_pseudo_axis_bruker()
            # Finding the spectrum sweep widths
            self.params.find_sw_bruker()
            # Finding the nucleus spectrometer frequencies
            self.params.find_nucleus_frequencies_bruker()
            # Finding the dimension labels
            self.params.find_labels_bruker()
            # Finding the temperature the experiment was performed at
            self.params.find_temperature_bruker()
            # Finding the carrier frequency of each dimension
            self.params.calculate_carrier_frequency_bruker()
            # Finding the digital filter parameters for bruker spectra
            self.params.find_bruker_digital_filter_parameters()
            # Find bruker scaling parameters
            self.params.find_bruker_scaling_parameters()
            # Finding bruker byte order and byte size
            self.params.determine_byte_order()
            self.params.determine_byte_size()

        else:
            # Spectrometer is Varian
            # Initialising the class to extract relevant varian parameters
            self.params = ParameterExtractorVarian(self)
            # Extracting spectrum dimensions
            self.params.find_size_varian()
            # Determine if any axis is a pseudo (non-complex) axis
            self.params.find_axes_pseudo_varian()
            # Finding the spectrum sweep widths
            self.params.find_sw_varian()
            # Finding the nucleus spectrometer frequencies
            self.params.find_nucleus_frequencies_varian()
            # Finding the dimension labels
            self.params.find_labels_varian()
            # Finding the temperature the experiment was performed at
            self.params.find_temperature_varian()
            # Finding the carrier frequency of each dimension
            self.params.calculate_carrier_frequency_varian()
            # Find varian scaling parameters
            self.params.find_varian_scaling_parameters()
