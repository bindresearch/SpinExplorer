#!/usr/bin/env python3

"""
MIT License

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
SOFTWARE.
"""

import wx
import os
import nmrglue as ng
from typing import Union
import json


class ReadFID:
    def __init__(self, app):
        """
        This class is for reading any fid file using nmrglue.
        If multiple nmrpipe format fid's are present then
        a pop out will ask the user which ones they would like
        to process.
        """

        self.app = app

        self.find_fid()
        self.read_fid()
        self.get_dimensions()
        self.find_pseudo_axes()
        self.guessing_FT_modes()
        self.find_sweep_widths()

    def find_fid(self) -> None:
        """
        Searching through the current directories to find converter fids
        (with extension .fid)
        """
        # Find nmrPipe fid files in the current directory
        fid_files = [file for file in os.listdir() if file.endswith(".fid")]
        if len(fid_files) == 0:
            dlg = wx.MessageDialog(
                self.app,
                "No nmrPipe file found in current directory",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            self.app.Raise()
            self.app.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            exit()
        elif len(fid_files) > 1:
            res = ChooseFile(fid_files, self)
            res.ShowModal()
            res.Destroy()
        else:
            self.fid_file = fid_files[0]

    def read_fid(self) -> None:
        """
        Reading the nmrPipe FID file using nmrglue
        """
        self.dic, self.data = ng.pipe.read(self.fid_file)

    def guessing_FT_modes(self) -> None:
        """
        This function will guess the Fourier transform modes
        for the indirect dimensions depending on the acquisition
        modes in parameters.json.

        By default the selection will be 0 (normal fourier transform)
        If acquisition mode:
        = TPPI - selection is 1 (real fourier transform)
        = States-TPPI - selection is 4 (sign alternation Fourier transform)

        Direct dimension will always be 0
        """

        self.ft_options = [0]
        if len(self.data.shape) > 1:
            try:
                with open("parameters.json", "r") as file:
                    acqusition_modes = json.load(file)["conversion"][
                        "spectral parameters"
                    ]["acqusition modes"]["indirect"]["mode"]
                    size = len(self.data.shape) - 1
                    for i in range(size):
                        mode = acqusition_modes[i]
                        if self.pseudo_axis == True:
                            if i == self.index and len(self.data.shape) == 3:
                                continue
                        if mode == "TPPI":
                            self.ft_options.append(1)
                        elif mode == "States-TPPI":
                            self.ft_options.append(4)
                        else:
                            self.ft_options.append(0)
            except:
                for i in range(len(self.data.shape) - 1):
                    self.ft_options.append(0)

    def get_dimensions(self):
        """
        Get the data size from the data and read the nmrglue dictionary
        to get axislabels and the size of each nmr dimension.
        """

        self.dim = len(self.data.shape)
        self.number_of_points = list(self.data.shape[::-1])

        self.axislabels = []
        if self.dim == 1:
            # If 1D take FDF1LABEL
            self.axislabels.append(self.dic["FDF2LABEL"])
        elif self.dim == 2:
            # If 2D take FDF2LABEL as direct and FDF1LABEL as indirect
            self.axislabels.append(self.dic["FDF2LABEL"])
            self.axislabels.append(self.dic["FDF1LABEL"])
        else:
            # If 3D take FDF3LABEL as direct, FDF2LABEL as indirect1 and FDF3LABEL as indirect2
            self.axislabels.append(self.dic["FDF2LABEL"])
            self.axislabels.append(self.dic["FDF1LABEL"])
            self.axislabels.append(self.dic["FDF3LABEL"])

    def find_sweep_widths(self) -> list[Union[int, float]]:
        """
        Finding the spectral sweep widths from the nmrglue
        dictionary. This is needed to properly show the
        window functions/apodization aligned with the FIDs.
        """
        self.spectral_width = []
        if self.dim == 1:
            # If 1D take FDF1SW
            self.spectral_width.append(self.dic["FDF2SW"])
        elif self.dim == 2:
            # If 2D take FDF2SW as direct and FDF1SW as indirect
            self.spectral_width.append(self.dic["FDF2SW"])
            self.spectral_width.append(self.dic["FDF1SW"])
        else:
            # If 3D take FDF3SW as direct, FDF2SW as indirect1 and FDF3SW as indirect2
            self.spectral_width.append(self.dic["FDF2SW"])
            self.spectral_width.append(self.dic["FDF1SW"])
            self.spectral_width.append(self.dic["FDF3SW"])

    def find_pseudo_axes(self):
        """
        Read through the nmrglue dictionary to find if there
        are any pseudo axes
        """
        self.pseudo_axis = False
        if self.dim == 1:
            # If 1D, check FDF1QUADFLAG (0 is complex, 1 is real)
            if int(self.dic["FDF2QUADFLAG"]) == 1:
                self.pseudo_axis = True
        elif self.dim == 2:
            # Check FDF1QUADFLAG and FDF2QUADFLAG (0 is complex, 1 is real)
            if int(self.dic["FDF1QUADFLAG"]) == 1:
                self.pseudo_axis = True
                self.index = 0
        else:
            # Check FDF1QUADFLAG/FDF2QUADFLAG/FDF3QUADFLAG (0 is complex, 1 is real)
            if int(self.dic["FDF2QUADFLAG"]) == 1:
                self.pseudo_axis = True
                self.index = 2
            elif int(self.dic["FDF1QUADFLAG"]) == 1:
                self.pseudo_axis = True
                self.index = 1
            elif int(self.dic["FDF3QUADFLAG"]) == 1:
                self.pseudo_axis = True
                self.index = 0


class ChooseFile(wx.Dialog):
    def __init__(self, spectrum_file, parent):
        wx.Dialog.__init__(
            self,
            None,
            wx.ID_ANY,
            "Select FID Data",
            wx.DefaultPosition,
            size=(300, 200),
        )
        self.spectrum_file = spectrum_file
        self.parent = parent
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.main_sizer.AddSpacer(10)
        self.message = wx.StaticText(
            self,
            label="Multiple FID files in current directory. Please select an an option to continue.\n",
        )
        self.main_sizer.Add(self.message, 0, wx.ALL, 5)
        self.file_combobox = wx.ComboBox(
            self, choices=spectrum_file, style=wx.CB_READONLY
        )
        self.main_sizer.Add(
            self.file_combobox, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 5
        )
        self.ok_button = wx.Button(self, label="OK")
        self.ok_button.Bind(wx.EVT_BUTTON, self.OnOK)
        self.main_sizer.Add(self.ok_button, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 5)
        self.SetSizer(self.main_sizer)
        self.Centre()

    def OnOK(self, event):
        file_selection = self.file_combobox.GetSelection()
        self.parent.fid_file = self.spectrum_file[file_selection]
        self.Close()
