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
SOFTWARE."""


print("-------------------------------------------------------------")
print("                          SpinView                           ")
print("-------------------------------------------------------------")
print("                (version 1.2) 20th June 2025                 ")
print(" (c) 2025 James Eaton, Andrew Baldwin (University of Oxford) ")
print("                        MIT License                          ")
print("-------------------------------------------------------------")
print("              Viewing and analysing NMR spectra              ")
print("-------------------------------------------------------------")
print(" Documentation at:")
print(" https://github.com/james-eaton-1/SpinExplorer")
print("-------------------------------------------------------------")
print("")


import sys
import wx

# Import relevant modules
import numpy as np
import matplotlib
from scipy.optimize import leastsq

matplotlib.use("wxAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import (
    NavigationToolbar2WxAgg as NavigationToolbar,
)
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import nmrglue as ng
import subprocess
import os
from scipy.interpolate import make_interp_spline

matplotlib.rcParams["font.sans-serif"] = "Arial"
matplotlib.rcParams["font.family"] = "sans-serif"


# Find out the version of operating system being used (Mac, Linux, Windows)
if sys.platform == "linux":
    platform = "linux"
    height = 30
elif sys.platform == "darwin":
    platform = "mac"
    height = 16
else:
    platform = "windows"
    height = 30

platform = "windows"

# James Eaton, 10/06/2025, University of Oxford
# This code will read in 1D, 2D and 3D NMRPipe (as well as topspin processed) data files and plot them
# It will also allow spectra to be phased, moved, multiplied, overlaid, etc.
# In addition, it allows the analysis of pseudo2D diffusion and relaxation (R1,R2) data


# This class reads in the NMRPipe data
class GetData:
    def __init__(self, file=""):

        # Create a hidden frame to be used as a parent for popout messages
        self.tempframe = wx.Frame(None, title="Temporary Parent", size=(1, 1))
        self.tempframe.Hide()  # Hide the frame since we don't need it to be visible

        self.file = file
        self.path = os.getcwd()
        if self.file == "":
            self.get_filename()
        self.read_data()
        self.dim = self.get_dimensions()
        if(self.file != "."):
            # NMRPipe data
            self.get_axislabels_nmrglue()
        else:
            # Bruker Topspin processed data
            self.generic_labels_bruker()

    # Get the filename of the NMRPipe data file
    def get_filename(self):
        self.found_file = False
        current_directory = os.getcwd()
        files = os.listdir(current_directory)
        spectrum_file = []
        self.brukerdata = False
        for file in files:
            if file.endswith(".ft"):
                spectrum_file.append(file)
            if file.endswith(".ft2"):
                spectrum_file.append(file)
            if file.endswith(".ft3"):
                spectrum_file.append(file)
            if file in [
                "1r",
                "1i",
                "2rr",
                "2ri",
                "3rrr",
                "3rri",
                "3rir",
                "3rii",
                "3irr",
                "3iri",
                "3iir",
                "3iii",
            ]:
                # Topspin processed Bruker data is present
                spectrum_file.append(".")
                break

        if len(spectrum_file) == 0:
            dlg = wx.MessageDialog(
                self.tempframe,
                "No NMRPipe or Bruker data files in current directory. Please convert and process data first.",
                "Error",
                wx.OK | wx.ICON_INFORMATION,
            )
            self.tempframe.Raise()
            self.tempframe.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            exit()
        if len(spectrum_file) == 1:
            self.file = spectrum_file[0]
        if len(spectrum_file) > 1:
            res = ChooseFile(spectrum_file, self)
            res.Raise()
            res.SetFocus()
            res.ShowModal()
            res.Destroy()

    # Read in the NMRPipe data file
    def read_data(self):
        self.found_file = False
        try:
            if self.file != ".":
                self.dic, self.data = ng.pipe.read(self.file)
            else:
                self.dic, self.data = ng.bruker.read_pdata(self.file)
            if len(self.data) == 0:
                # Give a popout saying the NMRPipe file has not been read properly. Retry processing
                dlg = wx.MessageDialog(
                    self.tempframe,
                    "Data file was read but data array is empty. Ensure raw data is downloaded to local device.",
                    "Error",
                    wx.OK | wx.ICON_INFORMATION,
                )
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                dlg.ShowModal()
                dlg.Destroy()
                self.found_file = True
                exit()

        except:
            if self.found_file == False:
                # Give a popout saying the NMRPipe file has not been read properly. Retry processing
                dlg = wx.MessageDialog(
                    self.tempframe,
                    "NMRPipe file not read properly. Please retry processing the file then try again.",
                    "Error",
                    wx.OK | wx.ICON_INFORMATION,
                )
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                dlg.ShowModal()
                dlg.Destroy()
            exit()

    # Work out NMR spectrum dimensions in order to get the plotting correct (need contour plot for 2D/3D but not for 1D)
    def get_dimensions(self):
        if (
            type(self.data[0]) == np.float32
            or type(self.data[0]) == np.float64
            or type(self.data[0]) == np.complex64
        ):
            return 1
        if len(self.data.shape) == 2:
            pseudo = False
            for val in self.data.shape:
                if val == 1:
                    pseudo = True
            if pseudo == True:
                self.data = self.data[0]
                return 1
            else:
                return 2
        if len(self.data.shape) == 3:
            pseudo = False
            for val in self.data.shape:
                if val == 1:
                    pseudo = True
            if pseudo == True:
                self.data_new = []
                for i, val2 in enumerate(self.data):
                    if self.data.shape[i] != 1:
                        self.data_new.append(val2)
                self.data = self.data_new
                return 3
            else:
                return 3

    # def read_labels_file(self):
    #     file = open("labels.txt", "r")
    #     label = file.readlines()
    #     for i, line in enumerate(label):
    #         if i == 0:
    #             line = line.split("\n")[0].split(",")
    #             self.axislabels = line
    #     file.close()

    # def read_spectrum_header(self):
    #     self.axislabels = []

    #     if(self.dim == 1):
    #         # If 1D take FDF1LABEL
    #         self.axislabels.append(self.dic["FDF1LABEL"])
    #     elif(self.dim == 2):
    #         # If 2D take FDF2LABEL as direct and FDF1LABEL as indirect
    #         self.axislabels.append(self.dic['FDF1LABEL'])
    #         self.axislabels.append(self.dic['FDF2LABEL'])
    #     else:
    #         # If 3D take FDF3LABEL as direct, FDF2LABEL as indirect1 and FDF3LABEL as indirect3
    #         self.axislabels.append(self.dic['FDF1LABEL'])
    #         self.axislabels.append(self.dic['FDF2LABEL'])
    #         self.axislabels.append(self.dic['FDF3LABEL'])

        

    #     # # Try the command showhdr to get the axis labels
    #     # command = "showhdr " + self.file
    #     # # output = subprocess.check_output(command)
    #     # p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    #     # (output, err) = p.communicate()
    #     # ## Wait for command to terminate. Get return returncode ##
    #     # p_status = p.wait()
    #     # output = output.decode()
    #     # output = output.split("\n")
    #     # if output != [""]:
    #     #     for line in output:
    #     #         if "NAME" in line:
    #     #             line = line.split()[1:]
    #     #             self.axislabels = line
    #     #             if self.dim == 3:
    #     #                 self.axislabels.reverse()

    # def read_fid_com_file(self):
    #     # Open fid.com file and get the axis labels
    #     file = open("fid.com", "r")
    #     fid_com = file.readlines()
    #     file.close()
    #     for line in fid_com:
    #         if "LAB" in line:
    #             line = line.split("\n")[0].split()
    #             # deleting the first element of the list which is the word 'LAB'
    #             del line[0]
    #             # deleting the last element of the list which is the '\' character
    #             del line[-1]
    #             if self.dim == 1:
    #                 self.labels.append(line[0])
    #             elif self.dim == 2:
    #                 self.labels.append(line[0])
    #                 self.labels.append(line[2])
    #             elif self.dim == 3:
    #                 self.labels.append(line[0])
    #                 self.labels.append(line[2])
    #                 self.labels.append(line[4])
    #     self.axislabels = self.labels
    #     if self.dim == 3:
    #         self.axislabels.reverse()

    def get_axislabels_nmrglue(self):
        """
        Reading the nmrglue dictionary (self.dic) to obtain the correct axis
        labels associated with the data.
        """

        self.axislabels = []

        if(self.dim == 1):
            # If 1D take FDF1LABEL
            self.axislabels.append(self.dic["FDF1LABEL"])
        elif(self.dim == 2):
            # If 2D take FDF2LABEL as direct and FDF1LABEL as indirect
            self.axislabels.append(self.dic['FDF1LABEL'])
            self.axislabels.append(self.dic['FDF2LABEL'])
        else:
            # If 3D take FDF3LABEL as direct, FDF2LABEL as indirect1 and FDF3LABEL as indirect3
            self.axislabels.append(self.dic['FDF1LABEL'])
            self.axislabels.append(self.dic['FDF2LABEL'])
            self.axislabels.append(self.dic['FDF3LABEL'])



    def generic_labels_bruker(self):
        """
        Input generic dim1, dim2, dim3 axis labels for Topspin
        processed data. This is temporary and should be updated
        to include correct labels from the Bruker dictionary.
        """
        self.labels = []

        if self.dim == 1:
            self.labels = ["dim1"]
        if self.dim == 2:
            self.labels = ["dim1", "dim2"]
        if self.dim == 3:
            self.labels = ["dim1", "dim2", "dim3"]
        self.axislabels = self.labels

    # def get_axislabels(self):
    #     self.labels = []
    #     # Try to find the labels in fid.com file, then try to match the dimension size of the label in the fid.com file to the dimension size of the data
    #     # This is to ensure that the correct labels are used for the correct dimension
    #     if platform == "linux" or platform == "mac":
    #         try:
    #             # If the user has already opened and customised the labels they will be in the labels.txt file
    #             self.read_labels_file()
    #         except:
    #             self.read_spectrum_header()
    #         else:
    #             try:
    #                 # Open fid.com file and get the axis labels
    #                 self.read_fid_com_file()
    #             except:
    #                 self.axislabels = []

    #         if self.axislabels == [] or len(self.axislabels) != self.dim:
    #             # Uable to find axis labels automatically so set as dim1 and dim2 etc
    #             if self.dim == 1:
    #                 self.labels = ["dim1"]
    #             if self.dim == 2:
    #                 self.labels = ["dim1", "dim2"]
    #             if self.dim == 3:
    #                 self.labels = ["dim1", "dim2", "dim3"]
    #             self.axislabels = self.labels

    #     else:
    #         try:
    #             self.read_labels_file()
    #         except:
    #             try:
    #                 self.read_fid_com_file()
    #                 if len(self.axislabels) != self.dim:
    #                     if self.dim == 1:
    #                         self.labels = ["dim1"]
    #                     if self.dim == 2:
    #                         self.labels = ["dim1", "dim2"]
    #                     if self.dim == 3:
    #                         self.labels = ["dim1", "dim2", "dim3"]

    #             except:
    #                 # Uable to find axis labels automatically so set as dim1 and dim2 etc
    #                 if self.dim == 1:
    #                     self.labels = ["dim1"]
    #                 if self.dim == 2:
    #                     self.labels = ["dim1", "dim2"]
    #                 if self.dim == 3:
    #                     self.labels = ["dim1", "dim2", "dim3"]

    #                 self.axislabels = self.labels


class ChooseFile(wx.Dialog):
    def __init__(self, spectrum_file, parent, session_choice=False):
        if session_choice == False:
            name = "Select NMRPipe Data File"
        else:
            name = "Select Session File"
        dialog = wx.Dialog.__init__(
            self,
            None,
            wx.ID_ANY,
            name,
            wx.DefaultPosition,
            size=(300, 200),
            style=wx.DEFAULT_DIALOG_STYLE,
        )
        self.spectrum_file = spectrum_file
        self.parent = parent
        self.session_choice = session_choice
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.main_sizer.AddSpacer(10)
        if self.session_choice == False:
            self.message = wx.StaticText(
                self,
                label="Multiple NMRPipe data files in current directory. Please select an NMRPipe file to show.\n",
            )
        else:
            self.message = wx.StaticText(
                self,
                label="Multiple session files in current directory. Please select a session file to load.\n",
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
        self.Show()

    def OnOK(self, event):
        file_selection = self.file_combobox.GetSelection()
        self.parent.file = self.spectrum_file[file_selection]
        self.parent.session_file = self.spectrum_file[file_selection]
        self.Close()
        if self.session_choice == False:
            self.parent.read_data()
            self.parent.dim = self.parent.get_dimensions()
            self.parent.get_axislabels()


# This class creates sliders which can contain floating point values
# Source: https://stackoverflow.com/questions/4709087/wxslider-with-floating-point-values
class FloatSlider(wx.Slider):

    def __init__(
        self,
        parent,
        id,
        value,
        minval,
        maxval,
        res,
        size=wx.DefaultSize,
        style=wx.SL_HORIZONTAL,
        name="floatslider",
    ):
        self._value = value
        self._min = minval
        self._max = maxval
        self._res = res
        ival, imin, imax = [round(v / res) for v in (value, minval, maxval)]
        self._islider = super(FloatSlider, self)
        self._islider.__init__(
            parent, id, ival, imin, imax, size=size, style=style, name=name
        )
        self.Bind(wx.EVT_SCROLL, self._OnScroll)

    def _OnScroll(self, event):
        ival = self._islider.GetValue()
        imin = self._islider.GetMin()
        imax = self._islider.GetMax()
        if ival == imin:
            self._value = self._min
        elif ival == imax:
            self._value = self._max
        else:
            self._value = ival * self._res
        event.Skip()

    def GetValue(self):
        return self._value

    def GetMin(self):
        return self._min

    def GetMax(self):
        return self._max

    def GetRes(self):
        return self._res

    def SetValue(self, value):
        self._islider.SetValue(round(value / self._res))
        self._value = value

    def SetMin(self, minval):
        self._islider.SetMin(round(minval / self._res))
        self._min = minval

    def SetMax(self, maxval):
        self._islider.SetMax(round(maxval / self._res))
        self._max = maxval

    def SetRes(self, res):
        self._islider.SetRange(round(self._min / res), round(self._max / res))
        self._islider.SetValue(round(self._value / res))
        self._res = res

    def SetRange(self, minval, maxval):
        self._islider.SetRange(round(minval / self._res), round(maxval / self._res))
        self._min = minval
        self._max = maxval


# This class creates the GUI main frame


class MyApp(wx.Frame):
    def __init__(self):
        # Get the monitor size and set the window size to 85% of the monitor size
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        self.width = int(1.0 * sizes[0][0])
        self.height = int(0.875 * sizes[0][1])
        self.reprocess = False

        self.app_frame = wx.Frame.__init__(
            self,
            None,
            wx.ID_ANY,
            "SpinView",
            wx.DefaultPosition,
            size=(int(self.width), int(self.height)),
        )
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.main_sizer.AddSpacer(10)
        self.display_index_current = wx.Display.GetFromWindow(self)

        # Variables needed to set the correct path so code can be used with unidecFile parser
        self.path = ""
        self.cwd = ""
        self.file_parser = False

        # Find if there are any sessions saved in the current directory
        self.find_sessions()
        if self.session_file != "":
            ReadSession(self, self.session_file)
        else:
            self.nmrdata = GetData()
            if self.nmrdata.dim == 1:
                self.viewer = OneDViewer(parent=self, nmrdata=self.nmrdata)
                self.main_sizer.Add(self.viewer, 1, wx.EXPAND)
            elif self.nmrdata.dim == 2:
                self.viewer = TwoDViewer(parent=self, nmrdata=self.nmrdata)
                self.main_sizer.Add(self.viewer, 1, wx.EXPAND)
            elif self.nmrdata.dim == 3:
                self.viewer = ThreeDViewer(parent=self, nmrdata=self.nmrdata)
                self.main_sizer.Add(self.viewer, 1, wx.EXPAND)

        self.SetSizer(self.main_sizer)

        # Make negative contour lines solid
        matplotlib.rc("contour", negative_linestyle="solid")

        self.Show()
        self.Centre()

        # Bind method to check/resize the window when the frame is moved
        self.Bind(wx.EVT_MOVE, self.OnMoveFrame)

        # Bind method to resize the window when the frame is resized
        self.Bind(wx.EVT_SIZE, self.OnSizeFrame)

        self.Bind(wx.EVT_CLOSE, self.OnClose)

    def OnClose(self, event):
        # Save the session file if the user wants to save it
        if self.reprocess == True:
            return
        else:
            try:
                if self.nmrdata.dim == 1 or self.nmrdata.dim == 2:
                    dlg = wx.MessageDialog(
                        None,
                        "Do you want to save the session file?",
                        "Save Session File",
                        wx.YES_NO | wx.ICON_INFORMATION,
                    )
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    if result == wx.ID_YES:
                        dlg.Destroy()
                        try:
                            if self.nmrdata.dim == 1:
                                self.viewer.OnSaveSessionButton(wx.EVT_BUTTON)
                            else:
                                self.viewer.OnSaveSessionButton2D(wx.EVT_BUTTON)
                        except:
                            dlg2 = wx.MessageDialog(
                                None,
                                "Session file not saved properly.",
                                "Save Session File",
                                wx.OK | wx.ICON_INFORMATION,
                            )
                            self.Raise()
                            self.SetFocus()
                            result2 = dlg2.ShowModal()
                            dlg2.Destroy()

                    dlg.Destroy()
                self.Destroy()
                sys.exit()
            except:
                self.Destroy()
                sys.exit()

    def OnMoveFrame(self, event):
        # Get the new default display if the frame is moved
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        display_index = wx.Display.GetFromWindow(self)
        if display_index != self.display_index_current:
            self.width = int(1.0 * sizes[display_index][0])
            self.height = int(0.875 * sizes[display_index][1])
            self.SetSize((self.width, self.height))
            self.viewer.canvas.SetSize(
                (
                    self.width * 0.0104,
                    (self.height - self.viewer.bottom_sizer.GetMinSize()[1] - 100)
                    * 0.0104,
                )
            )
            self.viewer.fig.set_size_inches(
                self.width * 0.0104,
                (self.height - self.viewer.bottom_sizer.GetMinSize()[1] - 100) * 0.0104,
            )
            self.viewer.UpdateFrame()
            self.display_index_current = display_index

        event.Skip()

    def OnSizeFrame(self, event):
        # Get the new frame size
        self.width, self.height = self.GetSize()
        self.SetSize((self.width, self.height))
        self.viewer.canvas.SetSize(
            (
                self.width * 0.0104,
                (self.height - self.viewer.bottom_sizer.GetMinSize()[1] - 100) * 0.0104,
            )
        )
        self.viewer.fig.set_size_inches(
            self.width * 0.0104,
            (self.height - self.viewer.bottom_sizer.GetMinSize()[1] - 100) * 0.0104,
        )
        self.viewer.UpdateFrame()
        event.Skip()

    def find_sessions(self):
        self.sessions = []
        if self.path != "":
            os.chdir(self.path)
            files = os.listdir(self.path)
        else:
            files = os.listdir()
        for file in files:
            if file.endswith(".session"):
                self.sessions.append(file)
        if self.path != "":
            os.chdir(self.cwd)

        # If there are no found sessions then session flag needs to be set to False
        if len(self.sessions) == 0:
            self.session_file = ""
        elif len(self.sessions) == 1:
            # Ask the user if they want to load the session file
            dlg = wx.MessageDialog(
                None,
                "Session file found ({}). Do you want to load the session file?".format(
                    self.sessions[0]
                ),
                "Session File Found",
                wx.YES_NO | wx.ICON_INFORMATION,
            )
            self.Raise()
            self.SetFocus()
            result = dlg.ShowModal()
            if result == wx.ID_YES:
                self.session_file = self.sessions[0]
            else:
                self.session_file = ""
            dlg.Destroy()
        else:
            # Asking the user if they wish to load a session file
            dlg = wx.MessageDialog(
                None,
                "Multiple session files found. Do you want to load a session file?",
                "Session Files Found",
                wx.YES_NO | wx.ICON_INFORMATION,
            )
            self.Raise()
            self.SetFocus()
            result = dlg.ShowModal()
            if result == wx.ID_YES:
                dlg.Destroy()
                # Asking the user to select the session file they want to load
                res = ChooseFile(self.sessions, self, session_choice=True)
                res.ShowModal()
                res.Destroy()
            else:
                self.session_file = ""
                dlg.Destroy()

    # Initialising global app variables variables variables
    def set_variables(self):
        # Colours for 1D lines
        self.colours = [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
        ]
        self.colour_value = self.colours[0]

        # Initial 1D slice colour for 2D/3D spectra is set to navy
        self.colour_slice = "navy"

        # List of cmap colours for when overlaying multiple spectra
        self.cmap = "#e41a1c"
        self.cmap_neg = "#377eb8"
        self.twoD_colours = [
            "#e41a1c",
            "#377eb8",
            "#4daf4a",
            "#984ea3",
            "#ff7f00",
            "#ff33eb",
        ]
        self.twoD_label_colours = self.twoD_colours

        self.twoD_slices_horizontal = []
        self.twoD_slices_vertical = []

        # Range of the sliders to for moving spectra left/right/up/down
        self.reference_range_values = [
            "0.01",
            "0.1",
            "0.5",
            "1.0",
            "5.0",
            "10.0",
            "50.0",
        ]
        self.reference_range = float(self.reference_range_values[0])
        self.reference_rangeX = float(self.reference_range_values[0])
        self.reference_rangeY = float(self.reference_range_values[0])

        # Range of the sliders to for moving spectra up/down in 1D spectra
        self.vertical_range_values = [
            "0.01",
            "0.1",
            "0.5",
            "1.0",
            "10.0",
            "50.0",
            "100.0",
            "1000.0",
            "10000",
        ]

        # Range of the sliders to for multiplying 1D spectra
        self.multiply_range_values = [
            "1.01",
            "1.1",
            "1.5",
            "2",
            "5",
            "10",
            "50",
            "100",
            "1000",
            "10000",
            "100000",
            "1000000",
            "10000000",
            "100000000",
            "1000000000",
        ]

        # Initial x,y movements for referencing are set to zero
        self.x_movement = 0
        self.y_movement = 0

        # Multiplot mode is initially set to off
        self.multiplot_mode = False

        # Dictionary to store the values of the sliders for each spectrum in multiplot mode
        self.values_dictionary = {}

        # Initial multiply factor is 1
        self.multiply_factor = 1

        # 1D slice color of 2D spectra is initially set to navy
        self.slice_colour = "navy"

        # Initial colour/reference/vertical index from list of colours is set to 0
        self.index = 0
        self.ref_index = 0
        self.vertical_index = 0

        # List to hold the multiple 2D spectra in multiplot mode
        self.twoD_spectra = []

        self.linewidth = 1.0
        self.linewidth1D = 1.5

        self.x_difference = 0
        self.y_difference = 0

        # Initially set the transpose flag to False
        self.transpose = False
        self.transposed2D = False

        # Default options for pivot point for P1 phasing
        self.pivot_x_default = 0
        self.pivot_x = self.pivot_x_default

        self.pivot_y_default = 0
        self.pivot_y = self.pivot_y_default

        self.slice_mode = None

        # Suppressing complex warning from numpy - prevents the complex warning from being printed to terminal when phasing
        import warnings

        # warnings.simplefilter("ignore", np.ComplexWarning)  # For old numpy versions
        warnings.simplefilter(
            "ignore", np.exceptions.ComplexWarning
        )  # For new numpy versions

    def UpdateFrame(self):
        self.canvas.draw()
        self.canvas.Refresh()
        self.canvas.Update()
        self.panel.Refresh()
        self.panel.Update()

    def OnSysColourChanged(self, event):

        pass


# Frame for One-Dimensional NMR Spectra
class OneDViewer(wx.Panel):
    def __init__(self, parent, nmrdata, uc0=None):
        # Getting the monitor size and set the window size to 85% of the monitor size
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        # Getting the current display index
        self.display_index = wx.Display.GetFromWindow(parent)
        self.width = int(1.0 * sizes[self.display_index][0])
        self.height = int(0.875 * sizes[self.display_index][1])
        self.parent = parent
        self.uc0_initial = uc0
        self.stack = False
        self.uc0 = uc0
        wx.Panel.__init__(self, parent, id=wx.ID_ANY, size=(self.width, self.height))
        self.nmrdata = nmrdata
        self.set_initial_variables_1D()
        self.create_button_panel_1D()
        self.create_hidden_button_panel_1D()
        self.create_canvas_1D()
        self.add_to_main_sizer1D()
        self.draw_figure_1D()

        # if(self.multiplot_mode == True):
        #     for i in range(len(self.values_dictionary)):
        #         self.plot_combobox.SetSelection(i)
        #         self.OnSelectPlot()
        # else:
        #     self.colour_chooser.SetSelection(self.values_dictionary[self.active_plot_index]['color index'])
        #     self.linewidth_slider.SetValue(self.values_dictionary[self.active_plot_index]['linewidth'])
        #     self.reference_range_chooser.SetSelection(self.values_dictionary[self.active_plot_index]['move left/right range index'])
        #     self.OnReferenceCombo(wx.EVT_SCROLL)
        #     self.reference_slider.SetValue(self.values_dictionary[self.active_plot_index]['move left/right'])
        #     self.vertical_range_chooser.SetSelection(self.values_dictionary[self.active_plot_index]['move up/down range index'])
        #     self.OnVerticalCombo(wx.EVT_SCROLL)
        #     self.vertical_slider.SetValue(self.values_dictionary[self.active_plot_index]['move up/down'])
        #     self.multiply_range_chooser.SetSelection(int(self.values_dictionary[self.active_plot_index]['multiply range index']))
        #     self.OnMultiplyCombo(wx.EVT_SCROLL)
        #     self.multiply_slider.SetValue(self.values_dictionary[self.active_plot_index]['multiply value'])
        #     self.P0_slider.SetValue(self.values_dictionary[self.active_plot_index]['p0 Coarse'])
        #     self.P1_slider.SetValue(self.values_dictionary[self.active_plot_index]['p1 Coarse'])
        #     self.P0_slider_fine.SetValue(self.values_dictionary[self.active_plot_index]['p0 Fine'])
        #     self.P1_slider_fine.SetValue(self.values_dictionary[self.active_plot_index]['p1 Fine'])

        #     # Update the plot to reflect the previously saved values for the active plot
        #     self.OnColourChoice1D(wx.EVT_SCROLL)
        #     self.OnLinewidthScroll1D(wx.EVT_SCROLL)
        #     self.OnReferenceScroll1D(wx.EVT_SCROLL)
        #     self.OnVerticalScroll1D(wx.EVT_SCROLL)
        #     self.OnMultiplyScroll1D(wx.EVT_SCROLL)
        #     self.OnSliderScroll1D(wx.EVT_SCROLL)

    def add_to_main_sizer1D(self):
        # Creating the main sizer
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.main_sizer.Add(self.canvas, 10, wx.EXPAND)
        self.main_sizer.Add(self.toolbar, 0, wx.EXPAND)
        # Adding all sizers to the main sizer
        self.main_sizer.Add(self.bottom_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.show_bottom_sizer = True
        self.main_sizer.Add(self.show_button_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.main_sizer.Hide(self.show_button_sizer)
        self.SetSizer(self.main_sizer)

    def create_canvas_1D(self):
        # Creating the figure and canvas to draw on
        self.panel = wx.Panel(self)
        self.fig = Figure(
            figsize=(
                self.width * 0.0104,
                (self.height - self.bottom_sizer.GetMinSize()[1] - 100) * 0.0104,
            )
        )
        self.canvas = FigCanvas(self, -1, self.fig)
        self.toolbar = NavigationToolbar(self.canvas)

    # Initialising variables for the 1D frame
    def set_initial_variables_1D(self):

        # Colours for 1D lines
        self.colours = [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
            "black",
            "navy",
            "tan",
            "lightcoral",
            "maroon",
            "lightgreen",
            "deeppink",
            "fuchsia",
        ]
        self.colour_value = self.colours[0]

        # Range of the sliders to for moving spectra left/right/up/down
        self.reference_range_values = [
            "0.01",
            "0.1",
            "0.5",
            "1.0",
            "5.0",
            "10.0",
            "50.0",
            "100",
            "200",
            "500",
            "1000",
            "2000",
            "5000",
            "10000",
        ]
        self.reference_range = float(self.reference_range_values[0])
        self.reference_rangeX = float(self.reference_range_values[0])
        self.reference_rangeY = float(self.reference_range_values[0])

        # Range of the sliders to for moving spectra up/down in 1D spectra
        self.vertical_range_values = [
            "0.01",
            "0.1",
            "0.5",
            "1.0",
            "10.0",
            "50.0",
            "100.0",
            "1000.0",
            "10000",
            "100000",
            "1000000",
        ]

        # Range of the sliders to for multiplying 1D spectra
        self.multiply_range_values = [
            "1.01",
            "1.1",
            "1.5",
            "2",
            "5",
            "10",
            "50",
            "100",
            "1000",
            "10000",
            "100000",
            "1000000",
            "10000000",
            "100000000",
            "1000000000",
            "10000000000",
            "100000000000",
            "1000000000000",
        ]

        # Initial x,y movements for referencing are set to zero
        self.x_movement = 0
        self.y_movement = 0

        # Multiplot mode is initially set to off
        self.multiplot_mode = False

        # Dictionary to store the values of the sliders for each spectrum in multiplot mode
        self.values_dictionary = {}

        # Initial multiply factor is 1
        self.multiply_factor = 1

        # Initial colour/reference/vertical index from list of colours is set to 0
        self.index = 0
        self.ref_index = 0
        self.vertical_index = 0

        self.linewidth = 1.0
        self.linewidth1D = 1.5

        self.x_difference = 0
        self.y_difference = 0

        # Default options for pivot point for P1 phasing
        self.pivot_x_default = 0
        self.pivot_x = self.pivot_x_default

        self.pivot_y_default = 0
        self.pivot_y = self.pivot_y_default

        self.total_P0 = 0.0
        self.total_P1 = 0.0

        # Initially have no baseline spline
        self.data_spline = [0]

        # # # Suppress complex warning from numpy
        import warnings

        # warnings.simplefilter("ignore", np.ComplexWarning)  # For old numpy versions
        warnings.simplefilter(
            "ignore", np.exceptions.ComplexWarning
        )  # For new numpy versions

    def UpdateFrame(self):
        self.canvas.draw()
        self.canvas.Refresh()
        self.canvas.Update()
        self.panel.Refresh()
        self.panel.Update()

    def create_button_panel_1D(self):
        # Creating a button to choose between plots in 1D spectra
        self.select_plot_label = wx.StaticBox(self, -1, "Select Plot:")
        self.select_plot_sizer = wx.StaticBoxSizer(self.select_plot_label, wx.VERTICAL)
        self.plot_combobox = wx.ComboBox(
            self, choices=["Main Plot"], style=wx.CB_READONLY
        )
        self.plot_combobox.Bind(wx.EVT_COMBOBOX, self.OnSelectPlot)
        self.select_plot_sizer.Add(self.plot_combobox, 0, wx.ALL, 5)
        # Checkbox where can select all plots to be edited at the same time
        self.select_all_checkbox = wx.CheckBox(self, label="Select All")
        self.select_plot_sizer.Add(
            self.select_all_checkbox, 0, wx.ALIGN_CENTER_HORIZONTAL, 5
        )

        # Creating the phasing 1D sizer
        self.phasing_label = wx.StaticBox(self, -1, "Phasing:")
        self.phasing_sizer = wx.StaticBoxSizer(self.phasing_label, wx.VERTICAL)
        self.P0_label = wx.StaticText(self, label="P0 (Coarse):", size=(70, height))
        self.P1_label = wx.StaticText(self, label="P1 (Coarse):", size=(70, height))
        self.P0_slider = FloatSlider(
            self,
            id=-1,
            value=0,
            minval=-180,
            maxval=180,
            res=0.1,
            size=(int(self.parent.width / 5), height),
        )
        self.P1_slider = FloatSlider(
            self,
            id=-1,
            value=0,
            minval=-180,
            maxval=180,
            res=0.1,
            size=(int(self.parent.width / 5), height),
        )
        self.P0_slider.Bind(wx.EVT_SLIDER, self.OnSliderScroll1D)
        self.P1_slider.Bind(wx.EVT_SLIDER, self.OnSliderScroll1D)
        self.P0_label_fine = wx.StaticText(self, label="P0 (Fine):", size=(70, height))
        self.P1_label_fine = wx.StaticText(self, label="P1 (Fine):", size=(70, height))
        self.P0_slider_fine = FloatSlider(
            self,
            id=-1,
            value=0,
            minval=-10,
            maxval=10,
            res=0.01,
            size=(int(self.parent.width / 5), height),
        )
        self.P1_slider_fine = FloatSlider(
            self,
            id=-1,
            value=0,
            minval=-10,
            maxval=10,
            res=0.01,
            size=(int(self.parent.width / 5), height),
        )
        self.P0_slider_fine.Bind(wx.EVT_SLIDER, self.OnSliderScroll1D)
        self.P1_slider_fine.Bind(wx.EVT_SLIDER, self.OnSliderScroll1D)
        self.P0_total = wx.StaticText(self, label="P0 (Total):", size=(70, height))
        self.P1_total = wx.StaticText(self, label="P1 (Total):", size=(70, height))
        self.P0_total_value = wx.StaticText(self, label="0.00", size=(70, height))
        self.P1_total_value = wx.StaticText(self, label="0.00", size=(70, height))

        # Adding a button to set the pivot point for phasing
        self.pivot_button = wx.Button(self, label="Set Pivot Point")
        self.pivot_button.Bind(wx.EVT_BUTTON, self.OnPivotButton)
        self.pivot_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.pivot_sizer.Add(self.pivot_button)

        # Adding a button to remove the pivot point
        self.remove_pivot_button = wx.Button(self, label="Remove Pivot Point")
        self.remove_pivot_button.Bind(wx.EVT_BUTTON, self.OnRemovePivotButton)
        self.pivot_sizer.AddSpacer(10)
        self.pivot_sizer.Add(self.remove_pivot_button)

        self.p0_sizer_labels = wx.BoxSizer(wx.VERTICAL)
        self.p0_sizer_labels.Add(self.P0_label)
        self.p0_sizer_labels.AddSpacer(10)
        self.p0_sizer_labels.Add(self.P0_label_fine)
        self.p0_sizer_labels.AddSpacer(10)
        self.p0_sizer_labels.Add(self.P0_total)

        self.p0_sizer_sliders = wx.BoxSizer(wx.VERTICAL)
        self.p0_sizer_sliders.Add(self.P0_slider, wx.ALIGN_CENTER_HORIZONTAL, 0)
        self.p0_sizer_sliders.AddSpacer(10)
        self.p0_sizer_sliders.Add(self.P0_slider_fine, wx.ALIGN_CENTER_HORIZONTAL, 0)
        self.p0_sizer_sliders.AddSpacer(10)
        self.p0_sizer_sliders.Add(self.P0_total_value, wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.p1_sizer_labels = wx.BoxSizer(wx.VERTICAL)
        self.p1_sizer_labels.Add(self.P1_label)
        self.p1_sizer_labels.AddSpacer(10)
        self.p1_sizer_labels.Add(self.P1_label_fine)
        self.p1_sizer_labels.AddSpacer(10)
        self.p1_sizer_labels.Add(self.P1_total)

        self.p1_sizer_sliders = wx.BoxSizer(wx.VERTICAL)
        self.p1_sizer_sliders.Add(self.P1_slider, wx.ALIGN_CENTER_HORIZONTAL, 0)
        self.p1_sizer_sliders.AddSpacer(10)
        self.p1_sizer_sliders.Add(self.P1_slider_fine, wx.ALIGN_CENTER_HORIZONTAL, 0)
        self.p1_sizer_sliders.AddSpacer(10)
        self.p1_sizer_sliders.Add(self.P1_total_value, wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.phasing_sizer1 = wx.BoxSizer(wx.HORIZONTAL)
        self.phasing_sizer1.AddSpacer(5)
        self.phasing_sizer1.Add(self.p0_sizer_labels, wx.ALIGN_TOP)
        self.phasing_sizer1.AddSpacer(10)
        self.phasing_sizer1.Add(self.p0_sizer_sliders, wx.ALIGN_TOP)
        self.phasing_sizer1.AddSpacer(10)
        self.phasing_sizer1.Add(self.p1_sizer_labels, wx.ALIGN_TOP)
        self.phasing_sizer1.AddSpacer(10)
        self.phasing_sizer1.Add(self.p1_sizer_sliders, wx.ALIGN_TOP)
        self.phasing_sizer1.AddSpacer(5)

        self.phasing_sizer.Add(self.phasing_sizer1)
        self.phasing_sizer.AddSpacer(10)
        self.phasing_sizer.Add(self.pivot_sizer, wx.ALIGN_CENTER, 1)

        # Creating a sizer for changing the y axis limits in the spectrum
        self.contour_label = wx.StaticBox(self, -1, "Y Axis Zoom (%):")
        self.contour_sizer = wx.StaticBoxSizer(self.contour_label, wx.VERTICAL)
        width = int(self.parent.width / 4.5)
        self.intensity_slider = FloatSlider(
            self, id=-1, value=0, minval=-1, maxval=10, res=0.01, size=(width, height)
        )
        self.intensity_slider.Bind(wx.EVT_SLIDER, self.OnIntensityScroll1D)
        self.contour_sizer.AddSpacer(5)
        self.contour_sizer.Add(self.intensity_slider)

        # Creating a slider for referencing a 1D spectrum (move spectrum left/right)
        total_zoom_width = self.contour_sizer.GetMinSize()[0] - 15
        if total_zoom_width < 150:
            width = int(total_zoom_width * 0.4)
            slider_width = int(total_zoom_width * 0.6)
        else:
            width = 55
            slider_width = total_zoom_width - 15 - 55
        self.reference_label = wx.StaticBox(self, -1, "Move \u2190/\u2192 (ppm):")
        self.reference_total = wx.StaticBoxSizer(self.reference_label, wx.VERTICAL)
        self.reference_sizer_full = wx.BoxSizer(wx.HORIZONTAL)
        self.reference_sizer = wx.BoxSizer(wx.VERTICAL)
        self.reference_sizer2 = wx.BoxSizer(wx.VERTICAL)
        self.reference_slider = FloatSlider(
            self,
            id=-1,
            value=0,
            minval=-self.reference_range,
            maxval=self.reference_range,
            res=2 * self.reference_range / 1000,
            size=(slider_width, height),
        )
        self.reference_slider.Bind(wx.EVT_SLIDER, self.OnReferenceScroll1D)
        self.reference_sizer.Add(self.reference_slider, wx.ALIGN_CENTER_HORIZONTAL, 5)
        self.reference_range_chooser = wx.ComboBox(
            self,
            value=self.reference_range_values[0],
            choices=self.reference_range_values,
            size=(width, height),
        )
        self.reference_range_chooser.Bind(wx.EVT_COMBOBOX, self.OnReferenceCombo)
        self.reference_range_chooser.SetSelection(0)
        self.reference_sizer2.Add(
            self.reference_range_chooser, wx.ALIGN_CENTER_HORIZONTAL, 5
        )
        self.reference_value_label = wx.StaticText(self, label="0.0")
        self.reference_sizer.AddSpacer(5)
        self.reference_sizer.Add(
            self.reference_value_label, wx.ALIGN_CENTER_HORIZONTAL, 5
        )
        self.reference_range_text = wx.StaticText(self, label="Range")
        self.reference_sizer2.AddSpacer(5)
        self.reference_sizer2.Add(
            self.reference_range_text, wx.ALIGN_CENTER_HORIZONTAL, 5
        )
        self.reference_sizer_full.Add(self.reference_sizer)
        self.reference_sizer_full.AddSpacer(5)
        self.reference_sizer_full.Add(self.reference_sizer2)

        self.reference_total.Add(self.reference_sizer_full)

        if self.reference_total.GetMinSize()[0] < self.contour_sizer.GetMinSize()[0]:
            self.reference_sizer_full.AddSpacer(
                self.contour_sizer.GetMinSize()[0]
                - self.reference_total.GetMinSize()[0]
            )

        # Create a sizer to move the data vertically
        self.vertical_range = int(max(self.nmrdata.data))
        self.vertical_label = wx.StaticBox(self, -1, "Move \u2191/\u2193 (%):")
        self.vertical_sizer = wx.StaticBoxSizer(self.vertical_label, wx.HORIZONTAL)
        self.vertical_sizer2 = wx.BoxSizer(wx.VERTICAL)
        self.vertical_slider = FloatSlider(
            self,
            id=-1,
            value=0,
            minval=-self.vertical_range * float(self.vertical_range_values[0]) / 100,
            maxval=self.vertical_range * float(self.vertical_range_values[0]) / 100,
            res=self.vertical_range * float(self.vertical_range_values[0]) / 10000,
            size=(slider_width, height),
        )
        self.vertical_slider.Bind(wx.EVT_SLIDER, self.OnVerticalScroll1D)
        self.vertical_sizer.Add(self.vertical_slider)
        self.vertical_sizer.AddSpacer(5)
        self.vertical_range_chooser = wx.ComboBox(
            self,
            value=self.vertical_range_values[0],
            choices=self.vertical_range_values,
            size=(width, height),
        )
        self.vertical_range_chooser.Bind(wx.EVT_COMBOBOX, self.OnVerticalCombo)
        self.vertical_range_chooser.SetSelection(0)
        self.vertical_sizer2.Add(
            self.vertical_range_chooser, wx.ALIGN_CENTER_HORIZONTAL, 5
        )
        self.vertical_range_label = wx.StaticText(self, label="Range")
        self.vertical_sizer2.AddSpacer(5)
        self.vertical_sizer2.Add(
            self.vertical_range_label, wx.ALIGN_CENTER_HORIZONTAL, 5
        )

        self.vertical_sizer.Add(self.vertical_sizer2)

        # Creating a combobox to change the colour of the 1D spectrum
        self.colour_label = wx.StaticBox(self, -1, "1D Line Colour")
        self.colour_sizer = wx.StaticBoxSizer(self.colour_label, wx.VERTICAL)
        self.options = [
            "Blue",
            "Orange",
            "Green",
            "Red",
            "Purple",
            "Brown",
            "Pink",
            "Gray",
            "Lime",
            "Turquoise         ",
            "Black",
            "Navy",
            "Tan",
            "Light Coral",
            "Maroon",
            "Light Green",
            "Deep Pink",
            "Fuchsia",
        ]
        self.colour_chooser = wx.ComboBox(
            self, value=self.options[0], choices=self.options, size=(100, height)
        )
        self.colour_chooser.Bind(wx.EVT_COMBOBOX, self.OnColourChoice1D)
        self.colour_chooser.SetSelection(0)
        spacer = 15
        self.colour_sizer.AddSpacer(spacer)
        self.colour_sizer.Add(self.colour_chooser)
        self.colour_sizer.AddSpacer(spacer)

        # Creating a slider to change the linewidth of the 1D spectrum
        self.linewidth_label = wx.StaticBox(self, -1, "1D Line Width")
        self.linewidth_sizer = wx.StaticBoxSizer(self.linewidth_label, wx.VERTICAL)
        self.linewidth_slider = FloatSlider(
            self, id=-1, value=0.5, minval=0.1, maxval=2, res=0.1, size=(100, height)
        )
        self.linewidth_slider.Bind(wx.EVT_SLIDER, self.OnLinewidthScroll1D)
        spacer = 15
        self.linewidth_sizer.AddSpacer(spacer)
        self.linewidth_sizer.Add(self.linewidth_slider)
        self.linewidth_sizer.AddSpacer(spacer)

        # Creating a slider to multiply of the 1D spectrum, with a combobox to choose the range of the slider
        total_phasing_width = self.phasing_sizer.GetMinSize()[0]
        leftover_width = (
            total_phasing_width
            - self.select_plot_sizer.GetMinSize()[0]
            - self.colour_sizer.GetMinSize()[0]
            - self.linewidth_sizer.GetMinSize()[0]
            - 10
            - 3 * int(self.parent.width / 100)
            - 20
        )
        if leftover_width < 200:
            range_width = int(leftover_width * 0.4)
            slider_width = int(leftover_width * 0.6)
        else:
            leftover_width = leftover_width - 100
            slider_width = int(leftover_width)
            range_width = 100
        self.multiply_label = wx.StaticBox(self, -1, "Multiplication Factor:")
        self.multiply_total = wx.StaticBoxSizer(self.multiply_label, wx.VERTICAL)
        self.multiply_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.multiply_value = 1
        self.multiply_slider = FloatSlider(
            self,
            id=-1,
            value=1,
            minval=0.1,
            maxval=float(self.multiply_range_values[0]),
            res=float(self.multiply_range_values[0]) / 1000,
            size=(slider_width, height),
        )
        self.multiply_slider.Bind(wx.EVT_SLIDER, self.OnMultiplyScroll1D)
        self.multiply_sizer_column1 = wx.BoxSizer(wx.VERTICAL)
        self.multiply_sizer_column2 = wx.BoxSizer(wx.VERTICAL)
        self.multiply_sizer_column1.Add(
            self.multiply_slider, wx.ALIGN_CENTER_HORIZONTAL, border=5
        )
        self.multiply_sizer.AddSpacer(5)
        self.multiply_range_chooser = wx.ComboBox(
            self,
            value=self.multiply_range_values[0],
            choices=self.multiply_range_values,
            size=(range_width, height),
        )
        self.multiply_range_chooser.Bind(wx.EVT_COMBOBOX, self.OnMultiplyCombo)
        self.multiply_range_chooser.SetSelection(0)
        self.multiply_sizer_column2.Add(self.multiply_range_chooser, border=5)
        self.multiply_label_value = wx.StaticText(self, label="1.000")
        self.multiply_combobox_label = wx.StaticText(self, label="Range")
        self.multiply_sizer_column1.AddSpacer(5)
        self.multiply_sizer_column1.Add(
            self.multiply_label_value, wx.ALIGN_CENTER_HORIZONTAL, 5
        )
        self.multiply_sizer_column2.AddSpacer(5)
        self.multiply_sizer_column2.Add(
            self.multiply_combobox_label, wx.ALIGN_CENTER_HORIZONTAL, 5
        )
        self.multiply_sizer.Add(self.multiply_sizer_column1)
        self.multiply_sizer.AddSpacer(5)
        self.multiply_sizer.Add(self.multiply_sizer_column2)
        self.multiply_total.AddSpacer(5)
        self.multiply_total.Add(self.multiply_sizer)
        current_height = self.multiply_total.GetMinSize()[1]
        linewidth_height = self.linewidth_sizer.GetMinSize()[1]
        if current_height < linewidth_height:
            self.multiply_total.AddSpacer(linewidth_height - current_height)

        # Making button to find the maximum intensity of the 1D spectrum
        self.max_button = wx.Button(self, label="Calculate Intensity", size=(130, 30))
        self.max_button.Bind(wx.EVT_BUTTON, self.OnMaxButton)

        self.baseline = wx.Button(self, label="Baseline", size=(130, 30))
        self.baseline.Bind(wx.EVT_BUTTON, self.OnBaseline)

        # Making button to subtract one spectrum from another
        self.subtract_button = wx.Button(self, label="Subtract Spectra", size=(130, 30))
        self.subtract_button.Bind(wx.EVT_BUTTON, self.OnSubtractButton)

        # Button to reset the parameters
        self.reset_button = wx.Button(self, label="Reset Parameters", size=(130, 30))
        self.reset_button.Bind(wx.EVT_BUTTON, self.OnResetButton1D)

        # Button to reprocess a spectrum
        self.reprocess_button = wx.Button(self, label="Re-process", size=(130, 30))
        self.reprocess_button.Bind(wx.EVT_BUTTON, self.OnReprocessButton1D)

        # Button to save a spectrum as a new nmrpipe .ft file
        self.save_button = wx.Button(self, label="Save Spectrum", size=(130, 30))
        self.save_button.Bind(wx.EVT_BUTTON, self.OnSaveButton)

        # Button to save the current session
        self.save_session_button = wx.Button(self, label="Save Session", size=(130, 30))
        self.save_session_button.Bind(wx.EVT_BUTTON, self.OnSaveSessionButton)

        # Button to hide the options for viewing
        self.hide_button = wx.Button(self, label="Hide Options", size=(130, 30))
        self.hide_button.Bind(wx.EVT_BUTTON, self.OnHideButton)

        self.button_sizers = wx.BoxSizer(wx.VERTICAL)
        self.button_sizers.Add(self.max_button)
        self.button_sizers.AddSpacer(5)
        self.button_sizers.Add(self.baseline)
        self.button_sizers.AddSpacer(5)
        self.button_sizers.Add(self.reset_button)
        self.button_sizers.AddSpacer(5)
        self.button_sizers.Add(self.subtract_button)
        self.button_sizers.AddSpacer(5)
        self.button_sizers.Add(self.reprocess_button)
        self.button_sizers.AddSpacer(5)
        self.button_sizers.Add(self.save_button)
        self.button_sizers.AddSpacer(5)
        self.button_sizers.Add(self.save_session_button)
        self.button_sizers.AddSpacer(5)
        self.button_sizers.Add(self.hide_button)
        self.button_sizers.AddSpacer(5)

        # Put all sizers together
        self.intensity_reference_sizer = wx.BoxSizer(wx.VERTICAL)
        self.intensity_reference_sizer.Add(self.contour_sizer)
        if platform == "linux":
            spacer = 15
        else:
            spacer = 10
        self.intensity_reference_sizer.AddSpacer(spacer)
        self.intensity_reference_sizer.Add(self.reference_total)
        self.intensity_reference_sizer.AddSpacer(spacer)
        self.intensity_reference_sizer.Add(self.vertical_sizer)
        self.intensity_reference_sizer.AddSpacer(spacer)

        # Create a sizer for the left side of the panel and add the select plot and phasing sizers to it
        self.left_sizer = wx.BoxSizer(wx.VERTICAL)
        self.top_left_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.top_left_sizer.Add(self.select_plot_sizer)
        spacer1 = int(self.parent.width / 100)
        self.top_left_sizer.AddSpacer(spacer1)
        self.top_left_sizer.Add(self.colour_sizer)
        self.top_left_sizer.AddSpacer(spacer1)
        self.top_left_sizer.Add(self.linewidth_sizer)
        self.top_left_sizer.AddSpacer(spacer1)
        self.top_left_sizer.Add(self.multiply_total)

        self.left_sizer.Add(self.top_left_sizer)
        self.left_sizer.AddSpacer(5)
        self.left_sizer.Add(self.phasing_sizer)
        self.left_sizer.AddSpacer(5)
        self.bottom_right_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.bottom_right_sizer.Add(self.intensity_reference_sizer)
        self.bottom_right_sizer.AddSpacer(5)
        self.bottom_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.bottom_sizer.Add(self.left_sizer)
        self.bottom_sizer.AddSpacer(10)
        self.bottom_sizer.Add(self.bottom_right_sizer)
        self.bottom_sizer.AddSpacer(5)
        self.bottom_sizer.Add(self.button_sizers)

    def create_hidden_button_panel_1D(self):
        # Creating a button to show the options
        # All other buttons/sliders are hidden when in hidden mode
        self.show_button = wx.Button(self, label="Show Options")
        self.show_button.Bind(wx.EVT_BUTTON, self.OnHideButton)
        self.show_button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.show_button_sizer.Add(self.show_button, wx.ALIGN_CENTER, 5)
        self.show_button_sizer.AddSpacer(5)

    def OnHideButton(self, event):
        if self.show_bottom_sizer == True:
            # Hide the panel
            self.main_sizer.Hide(self.bottom_sizer)
            self.main_sizer.Show(self.show_button_sizer)
            self.UpdateFrame()
            self.Layout()
            self.show_bottom_sizer = False
        else:
            # Show the panel
            self.main_sizer.Show(self.bottom_sizer)
            self.main_sizer.Hide(self.show_button_sizer)
            self.show_bottom_sizer = True
            self.UpdateFrame()
            self.Layout()

    def OnBaseline(self, event):
        # Checking to see if in multiplot mode (baselining mode only allowed when viewing a single plot)
        if self.multiplot_mode == True:
            message = "Currently in multiplot mode - manual baselining is not available currently in multiplot mode."
            dlg = wx.MessageDialog(self, message, "Warning", wx.OK)
            dlg.ShowModal()
            dlg.Destroy()
            return

        # Opening up a message to ask the user to select on nodes for the spline baseline
        message = "Manual baselining: click on points in the spectrum to be used as nodes for the baselining. Then press the key b to calculate the spline and subtract the baseline. Press the key c to cancel and clear the baseline."
        dlg = wx.MessageDialog(self, message, "Manual Baseline", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()

        self.points = []

        self.points_plots = []

        cid = self.fig.canvas.mpl_connect("button_press_event", self.onclickbaseline)
        key_press_connect = self.fig.canvas.mpl_connect(
            "key_press_event", self.on_key_baseline
        )

    def onclickbaseline(self, event):
        if event.xdata != None and event.ydata != None:
            self.points.append([event.xdata, event.ydata])
            self.points_plots.append(self.ax.plot(event.xdata, event.ydata, "ro"))
            self.fig.canvas.draw()

    def on_key_baseline(self, event):
        if event.key == "b":
            self.points = np.array(self.points)
            self.points = self.points[self.points[:, 0].argsort()]
            # Interpolate the points
            x = self.points[:, 0]
            y = self.points[:, 1]
            # Find the indexes of max(x) and min(x) in the data
            max_index = np.abs(self.ppms - max(x)).argmin()
            min_index = np.abs(self.ppms - min(x)).argmin()

            xnew = np.linspace(
                max(x), min(x), num=np.abs(int(max_index - min_index)), endpoint=True
            )
            spl = make_interp_spline(x, y, k=3)
            self.data_spline = spl(xnew)
            (self.spline_plot,) = self.ax.plot(
                xnew, self.data_spline, color="tab:orange"
            )
            self.fig.canvas.draw()

            # To any points not within spline region, make sure to add these as zero values to self.data_spline
            before_zeros = np.zeros(max_index)
            after_zeros = np.zeros((len(self.data) - min_index))
            self.data_spline = np.concatenate(
                (before_zeros, self.data_spline, after_zeros)
            )

            self.phase1D()
            # Then disable baseline mode
            self.fig.canvas.mpl_disconnect("button_press_event")
            self.fig.canvas.mpl_disconnect("key_press_event")

        if event.key == "c":
            # Clear the baseline and disable the button/key presses
            # Then disable baseline mode
            self.data_spline = [0]
            self.points = []

            try:
                for plot in self.points_plots:
                    plot[0].remove()
            except:
                pass
            try:
                self.spline_plot.remove()
            except:
                pass

            self.phase1D()

    def OnLoadPeakList(self, event):
        # Opening up a file window asking the user to select the 1D peak list - must be in the format of 1st column = peak_name, 2nd column = peak_position
        dlg = wx.FileDialog(self, "Select the peak list", wildcard="", style=wx.FD_OPEN)
        dlg.SetDirectory(os.getcwd())
        if dlg.ShowModal() == wx.ID_OK:
            self.peaklist_file = dlg.GetPath()
        else:
            dlg.Destroy()
            return

        self.ReadPeakList()

    def ReadPeakList(self):
        try:
            file = open(self.peaklist_file)
            lines = file.readlines()
            file.close()
        except:
            message = "Unable to open and read peak list. Please ensure the peaklist selected is correct."
            dlg = wx.MessageDialog(self, message, "Warning", wx.OK)
            dlg.ShowModal()
            dlg.Destroy()
            return

        self.peak_names = []
        self.peak_locations = []

        # Need to determine if the peaklist is in NMR-STAR format or something else

    def OnSelectPlot(self, event):
        # Saving the updated values for the previous plot for colour, linewidth, referencing, vertical scroll, and phasing
        if self.multiplot_mode == True:
            self.values_dictionary[self.active_plot_index][
                "color index"
            ] = self.colour_chooser.GetSelection()
            self.values_dictionary[self.active_plot_index][
                "linewidth"
            ] = self.linewidth_slider.GetValue()
            self.values_dictionary[self.active_plot_index][
                "move left/right range index"
            ] = self.reference_range_chooser.GetSelection()
            self.values_dictionary[self.active_plot_index][
                "move left/right"
            ] = self.reference_slider.GetValue()
            self.values_dictionary[self.active_plot_index][
                "move up/down range index"
            ] = self.vertical_range_chooser.GetSelection()
            self.values_dictionary[self.active_plot_index][
                "move up/down"
            ] = self.vertical_slider.GetValue()
            self.values_dictionary[self.active_plot_index][
                "multiply value"
            ] = self.multiply_slider.GetValue()
            self.values_dictionary[self.active_plot_index][
                "multiply range index"
            ] = self.multiply_range_chooser.GetSelection()
            self.values_dictionary[self.active_plot_index][
                "p0 Coarse"
            ] = self.P0_slider.GetValue()
            self.values_dictionary[self.active_plot_index][
                "p1 Coarse"
            ] = self.P1_slider.GetValue()
            self.values_dictionary[self.active_plot_index][
                "p0 Fine"
            ] = self.P0_slider_fine.GetValue()
            self.values_dictionary[self.active_plot_index][
                "p1 Fine"
            ] = self.P1_slider_fine.GetValue()

        # Function to change the active plot when a user selects a new plot from the combobox
        self.multiplot_mode = True
        self.active_plot_index = self.plot_combobox.GetSelection()

        # Updating the values in the GUI to reflect the previously saved values for the active plot
        self.colour_chooser.SetSelection(
            self.values_dictionary[self.active_plot_index]["color index"]
        )
        self.linewidth_slider.SetValue(
            self.values_dictionary[self.active_plot_index]["linewidth"]
        )
        self.reference_range_chooser.SetSelection(
            self.values_dictionary[self.active_plot_index][
                "move left/right range index"
            ]
        )
        self.OnReferenceCombo(wx.EVT_SCROLL)
        self.reference_slider.SetValue(
            self.values_dictionary[self.active_plot_index]["move left/right"]
        )
        self.vertical_range_chooser.SetSelection(
            self.values_dictionary[self.active_plot_index]["move up/down range index"]
        )
        self.OnVerticalCombo(wx.EVT_SCROLL)
        self.vertical_slider.SetValue(
            self.values_dictionary[self.active_plot_index]["move up/down"]
        )
        self.multiply_range_chooser.SetSelection(
            int(self.values_dictionary[self.active_plot_index]["multiply range index"])
        )
        self.OnMultiplyCombo(wx.EVT_SCROLL)
        self.multiply_slider.SetValue(
            self.values_dictionary[self.active_plot_index]["multiply value"]
        )
        self.P0_slider.SetValue(
            self.values_dictionary[self.active_plot_index]["p0 Coarse"]
        )
        self.P1_slider.SetValue(
            self.values_dictionary[self.active_plot_index]["p1 Coarse"]
        )
        self.P0_slider_fine.SetValue(
            self.values_dictionary[self.active_plot_index]["p0 Fine"]
        )
        self.P1_slider_fine.SetValue(
            self.values_dictionary[self.active_plot_index]["p1 Fine"]
        )

        # Updating the plot to reflect the previously saved values for the active plot
        self.OnColourChoice1D(wx.EVT_SCROLL)
        self.OnLinewidthScroll1D(wx.EVT_SCROLL)
        self.OnReferenceScroll1D(wx.EVT_SCROLL)
        self.OnVerticalScroll1D(wx.EVT_SCROLL)
        self.OnMultiplyScroll1D(wx.EVT_SCROLL)
        self.OnSliderScroll1D(wx.EVT_SCROLL)

    def OnSaveSessionButton(self, event):
        # Function to save the current session
        # File menu popout to ask the user which directory to save the session in
        dlg = wx.FileDialog(
            self,
            "Save Session",
            wildcard="Session files (*.session)|*.session",
            style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT,
        )
        dlg.SetDirectory(os.getcwd())
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.save_session(path)
            dlg.Destroy()
        else:
            return

    def save_session(self, path):
        # Function to save the current session
        # Saving the current session as a .session file

        if self.uc0_initial == None:  # If not coming from stack mode, save normally
            with open(path, "w") as f:
                if self.stack == False:
                    f.write("1D\n")
                else:
                    f.write("1D stack\n")
                if self.multiplot_mode == False:
                    f.write("MultiplotMode:False\n")
                    # Check if the file exists
                    if platform == "windows":
                        if (
                            os.path.exists(
                                str(self.parent.nmrdata.path)
                                + "\\"
                                + str(self.parent.nmrdata.file)
                            )
                            == True
                        ):
                            f.write(
                                "file_path:"
                                + str(self.parent.nmrdata.path)
                                + "\\"
                                + str(self.parent.nmrdata.file)
                                + "\n"
                            )
                        else:
                            f.write("file_path:" + str(self.parent.nmrdata.file) + "\n")
                    else:
                        if (
                            os.path.exists(
                                str(self.parent.nmrdata.path)
                                + "/"
                                + str(self.parent.nmrdata.file)
                            )
                            == True
                        ):
                            f.write(
                                "file_path:"
                                + str(self.parent.nmrdata.path)
                                + "/"
                                + str(self.parent.nmrdata.file)
                                + "\n"
                            )
                        else:
                            f.write("file_path:" + str(self.parent.nmrdata.file) + "\n")
                    f.write("p0_coarse:" + str(self.P0_slider.GetValue()) + "\n")
                    f.write("p1_coarse:" + str(self.P1_slider.GetValue()) + "\n")
                    f.write("p0_fine:" + str(self.P0_slider_fine.GetValue()) + "\n")
                    f.write("p1_fine:" + str(self.P1_slider_fine.GetValue()) + "\n")
                    f.write("colour:" + str(self.colour_chooser.GetSelection()) + "\n")
                    f.write("linewidth:" + str(self.linewidth_slider.GetValue()) + "\n")
                    f.write(
                        "reference_range:"
                        + str(self.reference_range_chooser.GetSelection())
                        + "\n"
                    )
                    f.write(
                        "reference_value:"
                        + str(self.reference_slider.GetValue())
                        + "\n"
                    )
                    f.write(
                        "vertical_range:"
                        + str(self.vertical_range_chooser.GetSelection())
                        + "\n"
                    )
                    f.write(
                        "vertical_value:" + str(self.vertical_slider.GetValue()) + "\n"
                    )
                    f.write(
                        "multiply_range:"
                        + str(self.multiply_range_chooser.GetSelection())
                        + "\n"
                    )
                    f.write(
                        "multiply_value:" + str(self.multiply_slider.GetValue()) + "\n"
                    )
                    f.write("pivot_point:" + str(self.pivot_line.get_xdata()[0]) + "\n")
                    f.write("pivot_x:" + str(self.pivot_x) + "\n")
                    f.write(
                        "pivot_visible:" + str(self.pivot_line.get_visible()) + "\n"
                    )

                else:
                    f.write("MultiplotMode:True\n")
                    for i in range(len(self.values_dictionary)):
                        f.write(
                            "file_path:" + str(self.values_dictionary[i]["path"]) + "\n"
                        )
                        f.write(
                            "title:" + str(self.values_dictionary[i]["title"]) + "\n"
                        )
                        f.write(
                            "p0_coarse:"
                            + str(self.values_dictionary[i]["p0 Coarse"])
                            + "\n"
                        )
                        f.write(
                            "p1_coarse:"
                            + str(self.values_dictionary[i]["p1 Coarse"])
                            + "\n"
                        )
                        f.write(
                            "p0_fine:"
                            + str(self.values_dictionary[i]["p0 Fine"])
                            + "\n"
                        )
                        f.write(
                            "p1_fine:"
                            + str(self.values_dictionary[i]["p1 Fine"])
                            + "\n"
                        )
                        f.write(
                            "colour:"
                            + str(self.values_dictionary[i]["color index"])
                            + "\n"
                        )
                        f.write(
                            "linewidth:"
                            + str(self.values_dictionary[i]["linewidth"])
                            + "\n"
                        )
                        f.write(
                            "reference_range:"
                            + str(
                                self.values_dictionary[i]["move left/right range index"]
                            )
                            + "\n"
                        )
                        f.write(
                            "reference_value:"
                            + str(self.values_dictionary[i]["move left/right"])
                            + "\n"
                        )
                        f.write(
                            "vertical_range:"
                            + str(self.values_dictionary[i]["move up/down range index"])
                            + "\n"
                        )
                        f.write(
                            "vertical_value:"
                            + str(self.values_dictionary[i]["move up/down"])
                            + "\n"
                        )
                        f.write(
                            "multiply_range:"
                            + str(self.values_dictionary[i]["multiply range index"])
                            + "\n"
                        )
                        f.write(
                            "multiply_value:"
                            + str(self.values_dictionary[i]["multiply value"])
                            + "\n"
                        )
                        f.write(
                            "pivot_point:" + str(self.pivot_line.get_xdata()[0]) + "\n"
                        )
                        f.write("pivot_x:" + str(self.pivot_x) + "\n")
                        f.write(
                            "pivot_visible:" + str(self.pivot_line.get_visible()) + "\n"
                        )

        else:  # If coming from stack mode, save the session as a stack session
            # Asking the user if they want to save the session as a stack session
            if platform == "windows":
                directory = "\\stack_session"
            else:
                directory = "/stack_session"
            try:
                if self.parent.parent.parent.path != "":
                    path1 = self.parent.parent.parent.path + directory
                else:
                    path1 = os.getcwd() + directory
            except:
                path1 = os.getcwd() + directory
            if os.path.exists(path1) == False:
                os.mkdir(path1)
            f = open(path, "w")
            f.write("1D stack\n")
            f.write("MultiplotMode:True\n")
            for i in range(len(self.values_dictionary)):
                data = np.array(
                    self.values_dictionary[i]["original_data"]
                    * self.values_dictionary[i]["multiply value"]
                    + self.values_dictionary[i]["move up/down"]
                    * np.ones(len(self.values_dictionary[i]["original_data"]))
                )
                data = data.astype(np.float32)

                dic = self.values_dictionary[i]["dictionary"]
                obs = dic["FDF2OBS"]
                sw = dic["FDF2SW"]
                car = dic["FDF2CAR"]
                size = dic["FDF2TDSIZE"]
                label = dic["FDF2LABEL"]
                orig = dic["FDF2ORIG"]
                center = dic["FDF2CENTER"]
                udic = {
                    "ndim": 1,
                    0: {
                        "size": size,
                        "complex": False,
                        "encoding": "int",
                        "sw": sw,
                        "obs": obs,
                        "car": car,
                        "label": label,
                        "time": False,
                        "freq": True,
                    },
                }
                dic = ng.pipe.create_dic(udic)
                dic["FDF2LABEL"] = label
                dic["FDF2OBS"] = obs
                dic["FDF2SW"] = sw
                dic["FDF2CAR"] = car
                dic["FDF2SIZE"] = size
                dic["FDF2ORIG"] = orig
                dic["FDF2CENTER"] = center
                ng.pipe.write(
                    path1 + "/" + self.values_dictionary[i]["title"] + ".ft",
                    dic,
                    data,
                    overwrite=True,
                )

                f.write(
                    "file_path:"
                    + path1
                    + "/"
                    + str(self.values_dictionary[i]["title"])
                    + ".ft"
                    + "\n"
                )
                f.write("title:" + str(self.values_dictionary[i]["title"]) + "\n")
                f.write(
                    "p0_coarse:" + str(self.values_dictionary[i]["p0 Coarse"]) + "\n"
                )
                f.write(
                    "p1_coarse:" + str(self.values_dictionary[i]["p1 Coarse"]) + "\n"
                )
                f.write("p0_fine:" + str(self.values_dictionary[i]["p0 Fine"]) + "\n")
                f.write("p1_fine:" + str(self.values_dictionary[i]["p1 Fine"]) + "\n")
                f.write(
                    "colour:" + str(self.values_dictionary[i]["color index"]) + "\n"
                )
                f.write(
                    "linewidth:" + str(self.values_dictionary[i]["linewidth"]) + "\n"
                )
                f.write(
                    "reference_range:"
                    + str(self.values_dictionary[i]["move left/right range index"])
                    + "\n"
                )
                f.write(
                    "reference_value:"
                    + str(self.values_dictionary[i]["move left/right"])
                    + "\n"
                )
                f.write(
                    "vertical_range:"
                    + str(self.values_dictionary[i]["move up/down range index"])
                    + "\n"
                )
                f.write(
                    "vertical_value:"
                    + str(self.values_dictionary[i]["move up/down"])
                    + "\n"
                )
                f.write(
                    "multiply_range:"
                    + str(self.values_dictionary[i]["multiply range index"])
                    + "\n"
                )
                f.write(
                    "multiply_value:"
                    + str(self.values_dictionary[i]["multiply value"])
                    + "\n"
                )
                f.write("pivot_point:" + str(self.pivot_line.get_xdata()[0]) + "\n")
                f.write("pivot_x:" + str(self.pivot_x) + "\n")
                f.write("pivot_visible:" + str(self.pivot_line.get_visible()) + "\n")
            f.close()
            return

    def OnResetButton1D(self, event):
        if self.multiplot_mode == True:
            # Check to see if select all is checked
            if self.select_all_checkbox.GetValue() == True:
                message = "Select all checkbox is set to True. This action will reset all the parameters for all the plots. Do you want to continue?"
                dlg = wx.MessageDialog(
                    self, message, "Reset Parameters", wx.YES_NO | wx.ICON_QUESTION
                )
                result = dlg.ShowModal()
                if result == wx.ID_YES:
                    for i in range(len(self.values_dictionary)):
                        self.values_dictionary[i]["color index"] = i
                        self.values_dictionary[i]["linewidth"] = 0.5
                        self.values_dictionary[i]["move left/right range index"] = 0
                        self.values_dictionary[i]["move left/right"] = 0
                        self.values_dictionary[i]["move up/down range index"] = 0
                        self.values_dictionary[i]["move up/down"] = 0
                        self.values_dictionary[i]["multiply value"] = 1
                        self.values_dictionary[i]["multiply range index"] = 0
                        self.values_dictionary[i]["p0 Coarse"] = 0
                        self.values_dictionary[i]["p1 Coarse"] = 0
                        self.values_dictionary[i]["p0 Fine"] = 0
                        self.values_dictionary[i]["p1 Fine"] = 0
                        if i == 0:
                            self.line1.set_color(self.colours[0])
                        else:
                            self.extra_plots[i - 1][0].set_color(self.colours[i])

                    self.colour_chooser.SetSelection(self.active_plot_index)
                    self.linewidth_slider.SetValue(0.5)
                    self.reference_range_chooser.SetSelection(0)
                    self.OnReferenceCombo(wx.EVT_SCROLL)
                    self.reference_slider.SetValue(0)
                    self.vertical_range_chooser.SetSelection(0)
                    self.OnVerticalCombo(wx.EVT_SCROLL)
                    self.vertical_slider.SetValue(0)
                    self.multiply_range_chooser.SetSelection(0)
                    self.OnMultiplyCombo(wx.EVT_SCROLL)
                    self.multiply_slider.SetValue(1)
                    self.P0_slider.SetValue(0)
                    self.P1_slider.SetValue(0)
                    self.P0_slider_fine.SetValue(0)
                    self.P1_slider_fine.SetValue(0)

                    self.select_all_checkbox.SetValue(False)
                    self.OnColourChoice1D(wx.EVT_SCROLL)
                    self.select_all_checkbox.SetValue(True)

                    self.OnLinewidthScroll1D(wx.EVT_SCROLL)
                    self.OnReferenceScroll1D(wx.EVT_SCROLL)
                    self.OnVerticalScroll1D(wx.EVT_SCROLL)
                    self.OnMultiplyScroll1D(wx.EVT_SCROLL)
                    self.OnSliderScroll1D(wx.EVT_SCROLL)
                else:
                    return

            else:
                message = "This action will reset all the parameters for the current selected plot. Do you want to continue?"
                dlg = wx.MessageDialog(
                    self, message, "Reset Parameters", wx.YES_NO | wx.ICON_QUESTION
                )
                result = dlg.ShowModal()
                if result == wx.ID_YES:
                    self.values_dictionary[self.active_plot_index][
                        "color index"
                    ] = self.active_plot_index
                    self.values_dictionary[self.active_plot_index]["linewidth"] = 0.5
                    self.values_dictionary[self.active_plot_index][
                        "move left/right range index"
                    ] = 0
                    self.values_dictionary[self.active_plot_index][
                        "move left/right"
                    ] = 0
                    self.values_dictionary[self.active_plot_index][
                        "move up/down range index"
                    ] = 0
                    self.values_dictionary[self.active_plot_index]["move up/down"] = 0
                    self.values_dictionary[self.active_plot_index]["multiply value"] = 1
                    self.values_dictionary[self.active_plot_index][
                        "multiply range index"
                    ] = 0
                    self.values_dictionary[self.active_plot_index]["p0 Coarse"] = 0
                    self.values_dictionary[self.active_plot_index]["p1 Coarse"] = 0
                    self.values_dictionary[self.active_plot_index]["p0 Fine"] = 0
                    self.values_dictionary[self.active_plot_index]["p1 Fine"] = 0
                    self.colour_chooser.SetSelection(0)
                    self.linewidth_slider.SetValue(0.5)
                    self.reference_range_chooser.SetSelection(0)
                    self.OnReferenceCombo(wx.EVT_SCROLL)
                    self.reference_slider.SetValue(0)
                    self.vertical_range_chooser.SetSelection(0)
                    self.OnVerticalCombo(wx.EVT_SCROLL)
                    self.vertical_slider.SetValue(0)
                    self.multiply_range_chooser.SetSelection(0)
                    self.OnMultiplyCombo(wx.EVT_SCROLL)
                    self.multiply_slider.SetValue(1)
                    self.P0_slider.SetValue(0)
                    self.P1_slider.SetValue(0)
                    self.P0_slider_fine.SetValue(0)
                    self.P1_slider_fine.SetValue(0)
                    self.OnColourChoice1D(wx.EVT_SCROLL)
                    self.OnLinewidthScroll1D(wx.EVT_SCROLL)
                    self.OnReferenceScroll1D(wx.EVT_SCROLL)
                    self.OnVerticalScroll1D(wx.EVT_SCROLL)
                    self.OnMultiplyScroll1D(wx.EVT_SCROLL)
                    self.OnSliderScroll1D(wx.EVT_SCROLL)
                else:
                    return
        else:
            message = "This action will reset all the parameters for the plot. Do you want to continue?"
            dlg = wx.MessageDialog(
                self, message, "Reset Parameters", wx.YES_NO | wx.ICON_QUESTION
            )
            result = dlg.ShowModal()
            if result == wx.ID_YES:
                self.colour_chooser.SetSelection(0)
                self.linewidth_slider.SetValue(0.5)
                self.reference_range_chooser.SetSelection(0)
                self.OnReferenceCombo(wx.EVT_SCROLL)
                self.reference_slider.SetValue(0)
                self.vertical_range_chooser.SetSelection(0)
                self.OnVerticalCombo(wx.EVT_SCROLL)
                self.vertical_slider.SetValue(0)
                self.multiply_range_chooser.SetSelection(0)
                self.OnMultiplyCombo(wx.EVT_SCROLL)
                self.multiply_slider.SetValue(1)
                self.P0_slider.SetValue(0)
                self.P1_slider.SetValue(0)
                self.P0_slider_fine.SetValue(0)
                self.P1_slider_fine.SetValue(0)
                self.OnColourChoice1D(wx.EVT_SCROLL)
                self.OnLinewidthScroll1D(wx.EVT_SCROLL)
                self.OnReferenceScroll1D(wx.EVT_SCROLL)
                self.OnVerticalScroll1D(wx.EVT_SCROLL)
                self.OnMultiplyScroll1D(wx.EVT_SCROLL)
                self.OnSliderScroll1D(wx.EVT_SCROLL)
            else:
                return

    def OnReprocessButton1D(self, event):
        if self.multiplot_mode == False:
            # Opening an instance of SpinProcess
            if self.parent.path != "":
                os.chdir(self.parent.path)
            try:
                from SpinExplorer.SpinProcess.SpinProcess import SpinProcess
            except:
                # Output saying that SpinProcess is not available
                message = (
                    "Cannot find SpinProcess module - reprocessing is not possible"
                )
                dlg = wx.MessageDialog(
                    self, message, "Reprocess", wx.OK | wx.ICON_INFORMATION
                )
                dlg.ShowModal()
                dlg.Destroy()
                return

            reprocessing_frame = SpinProcess(
                self, path=self.parent.path, cwd=self.parent.cwd, reprocess=True
            )
            if self.parent.cwd != "":
                os.chdir(self.parent.cwd)
        else:
            # Checking to see if data has been originated from stack mode
            if self.uc0_initial != None or self.stack == True:
                # Popout saying that 1D reprocessing is not possible when 1D's are generated from stacking 2D data. Please reprocess the original 2D data and stack again
                message = "1D reprocessing is not possible when 1D's are generated from stacking 2D data. Please reprocess the original 2D data and stack again."
                dlg = wx.MessageDialog(
                    self, message, "Reprocess", wx.OK | wx.ICON_INFORMATION
                )
                dlg.ShowModal()
                dlg.Destroy()
                return

            # Give a popout saying that this will allow reprocessing of the currently selected plot (Do you want to continue)
            message = "This action will allow reprocessing of the currently selected plot. Do you want to continue?"
            dlg = wx.MessageDialog(
                self, message, "Reprocess", wx.YES_NO | wx.ICON_QUESTION
            )
            result = dlg.ShowModal()
            if result == wx.ID_YES:
                # Open an instance of SpinProcess
                path = self.values_dictionary[self.active_plot_index]["path"].split(
                    "/"
                )[0:-1]
                path = "/".join(path)
                try:
                    os.chdir(path)
                except:
                    # Give an error saying that the path was not found
                    dlg = wx.MessageDialog(
                        self,
                        "Path not found for plot ({}). Unable to reprocess this data.".format(
                            self.values_dictionary[self.active_plot_index]["title"]
                        ),
                        "Error",
                        wx.OK | wx.ICON_ERROR,
                    )
                    dlg.ShowModal()
                    dlg.Destroy()
                    return
                try:
                    from SpinExplorer.SpinProcess.SpinProcess import SpinProcess
                except:
                    # Output saying that SpinProcess is not available
                    message = (
                        "Cannot find SpinProcess module - reprocessing is not possible"
                    )
                    dlg = wx.MessageDialog(
                        self, message, "Reprocess", wx.OK | wx.ICON_INFORMATION
                    )
                    dlg.ShowModal()
                    dlg.Destroy()
                    return

                reprocessing_frame = SpinProcess(self, path=path, cwd=self.parent.cwd)
                reprocessing_frame.reprocess = True
                if self.parent.cwd != "":
                    os.chdir(self.parent.cwd)
            else:
                return

    def OnMaxButton(self, event):

        # Asking the user to select a region of the spectrum where they want to find the intensity
        dlg = wx.MessageBox(
            "Click and drag to select a region of the spectrum to find the maximum intensity/integral.",
            "Max Intensity",
            wx.OK | wx.ICON_INFORMATION,
        )

        self.intensity_region.set_visible(True)

        self.UpdateFrame()

        self.press = False
        self.move = False
        self.noise_select_press = self.canvas.mpl_connect(
            "button_press_event", self.OnPress
        )
        self.noise_select_release = self.canvas.mpl_connect(
            "button_release_event", self.OnReleaseNoise
        )
        self.noise_select_move = self.canvas.mpl_connect(
            "motion_notify_event", self.OnMove
        )

    def OnPress(self, event):
        if event.inaxes == self.ax:
            self.press = True
            self.x0 = event.xdata

    def OnMove(self, event):
        if event.inaxes == self.ax:
            self.move_intensity(event)

    def move_intensity(self, event):
        if self.press:
            self.move = True
            self.x1 = event.xdata
            if self.x1 > self.x0:
                xmax = self.x1
                xmin = self.x0
            else:
                xmax = self.x0
                xmin = self.x1

            self.intensity_region.set_x(xmin)
            self.intensity_region.set_width(xmax - xmin)
            # self.intensity_region.set_xy(xy=[[xmin,0],[xmin,1],[xmax,1],[xmax,0]]) # no longer works in recent matplotlib versions
            self.intensity_region.set_visible(True)

            self.UpdateFrame()

    def release_intensity(self, event):
        if self.press:
            self.x2 = event.xdata
            if self.x2 > self.x0:
                xmax = self.x2
                xmin = self.x0
            else:
                xmax = self.x0
                xmin = self.x2
            self.intensity_x_initial = xmin
            self.intensity_x_final = xmax
            self.intensity_region.set_x(xmin)
            self.intensity_region.set_width(xmax - xmin)
            # self.intensity_region.set(xy=[[xmin,0],[xmin,1],[xmax,1],[xmax,0]]) # no longer works in recent matplotlib versions

            self.UpdateFrame()
        self.press = False
        self.move = False
        self.canvas.mpl_disconnect(self.noise_select_press)
        self.canvas.mpl_disconnect(self.noise_select_move)
        self.canvas.mpl_disconnect(self.noise_select_release)

        # If in multiplot mode, output the values for the active plot
        if self.multiplot_mode == True:
            # Check to see if select all is checked
            if self.select_all_checkbox.GetValue() == True:
                message = "Maximum intensities in selected region for each slice:\n"
                for i in range(len(self.values_dictionary)):

                    # Find the index of the ppms in the region selected
                    self.intensity_index_initial = np.abs(
                        self.values_dictionary[i]["original_ppms"]
                        + self.values_dictionary[i]["move left/right"]
                        - self.intensity_x_final
                    ).argmin()
                    self.intensity_index_final = np.abs(
                        self.values_dictionary[i]["original_ppms"]
                        + self.values_dictionary[i]["move left/right"]
                        - self.intensity_x_initial
                    ).argmin()
                    try:
                        max_intensity = max(
                            np.real(
                                self.values_dictionary[i]["original_data"][
                                    self.intensity_index_initial : self.intensity_index_final
                                ]
                                * self.values_dictionary[i]["multiply value"]
                            )
                            + self.values_dictionary[i]["move up/down"]
                            * np.ones(
                                len(
                                    self.values_dictionary[i]["original_data"][
                                        self.intensity_index_initial : self.intensity_index_final
                                    ]
                                )
                            )
                        )
                        # Find the ppm of the maximum intensity
                        max_intensity_index = np.argmax(
                            np.real(
                                self.values_dictionary[i]["original_data"][
                                    self.intensity_index_initial : self.intensity_index_final
                                ]
                                * self.values_dictionary[i]["multiply value"]
                            )
                            + self.values_dictionary[i]["move up/down"]
                            * np.ones(
                                len(
                                    self.values_dictionary[i]["original_data"][
                                        self.intensity_index_initial : self.intensity_index_final
                                    ]
                                )
                            )
                        )
                        max_intensity_ppm = (
                            self.values_dictionary[i]["original_ppms"][
                                self.intensity_index_initial : self.intensity_index_final
                            ][max_intensity_index]
                            + self.values_dictionary[i]["move left/right"]
                        )
                        mean_intensity = np.mean(
                            np.real(
                                self.values_dictionary[i]["original_data"][
                                    self.intensity_index_initial : self.intensity_index_final
                                ]
                                * self.values_dictionary[i]["multiply value"]
                            )
                            + self.values_dictionary[i]["move up/down"]
                            * np.ones(
                                len(
                                    self.values_dictionary[i]["original_data"][
                                        self.intensity_index_initial : self.intensity_index_final
                                    ]
                                )
                            )
                        )
                        integral = np.sum(
                            np.real(
                                self.values_dictionary[i]["original_data"][
                                    self.intensity_index_initial : self.intensity_index_final
                                ]
                                * self.values_dictionary[i]["multiply value"]
                            )
                            + +self.values_dictionary[i]["move up/down"]
                            * np.ones(
                                len(
                                    self.values_dictionary[i]["original_data"][
                                        self.intensity_index_initial : self.intensity_index_final
                                    ]
                                )
                            )
                        )
                        stdev = np.std(
                            np.real(
                                self.values_dictionary[i]["original_data"][
                                    self.intensity_index_initial : self.intensity_index_final
                                ]
                                * self.values_dictionary[i]["multiply value"]
                            )
                            + +self.values_dictionary[i]["move up/down"]
                            * np.ones(
                                len(
                                    self.values_dictionary[i]["original_data"][
                                        self.intensity_index_initial : self.intensity_index_final
                                    ]
                                )
                            )
                        )
                        message += self.values_dictionary[i][
                            "title"
                        ] + ":\nPPM of Max Intensity: {:.3f}\nMax Intensity: {:E} \nMean Intensity {:E} \nIntegral: {:E} \nStandard deviation of intensity: {:E} \n\n".format(
                            max_intensity_ppm,
                            max_intensity,
                            mean_intensity,
                            integral,
                            stdev,
                        )
                        min_ppm = (
                            min(
                                self.values_dictionary[i]["original_ppms"][
                                    self.intensity_index_initial : self.intensity_index_final
                                ]
                            )
                            + self.values_dictionary[0]["move left/right"]
                        )
                        max_ppm = (
                            max(
                                self.values_dictionary[i]["original_ppms"][
                                    self.intensity_index_initial : self.intensity_index_final
                                ]
                            )
                            + self.values_dictionary[0]["move left/right"]
                        )
                    except:
                        message += (
                            self.values_dictionary[i]["title"]
                            + ":\nNo NMR data found in selected chemical shift range \n"
                        )

                message += "Selected PPM Range:\n{:.3f}-{:.3f}".format(min_ppm, max_ppm)
            else:
                # Find only the max intensity of the current active plot
                self.intensity_index_initial = np.abs(
                    self.values_dictionary[self.active_plot_index]["original_ppms"]
                    + self.values_dictionary[self.active_plot_index]["move left/right"]
                    - self.intensity_x_final
                ).argmin()
                self.intensity_index_final = np.abs(
                    self.values_dictionary[self.active_plot_index]["original_ppms"]
                    + self.values_dictionary[self.active_plot_index]["move left/right"]
                    - self.intensity_x_initial
                ).argmin()
                try:
                    max_intensity = max(
                        np.real(
                            self.values_dictionary[self.active_plot_index][
                                "original_data"
                            ][self.intensity_index_initial : self.intensity_index_final]
                            * self.values_dictionary[self.active_plot_index][
                                "multiply value"
                            ]
                        )
                        + self.values_dictionary[self.active_plot_index]["move up/down"]
                        * np.ones(
                            len(
                                self.values_dictionary[self.active_plot_index][
                                    "original_data"
                                ][
                                    self.intensity_index_initial : self.intensity_index_final
                                ]
                            )
                        )
                    )
                    max_intensity_index = np.argmax(
                        np.real(
                            self.values_dictionary[self.active_plot_index][
                                "original_data"
                            ][self.intensity_index_initial : self.intensity_index_final]
                            * self.values_dictionary[self.active_plot_index][
                                "multiply value"
                            ]
                        )
                        + self.values_dictionary[self.active_plot_index]["move up/down"]
                        * np.ones(
                            len(
                                self.values_dictionary[self.active_plot_index][
                                    "original_data"
                                ][
                                    self.intensity_index_initial : self.intensity_index_final
                                ]
                            )
                        )
                    )
                    max_intensity_ppm = (
                        self.values_dictionary[self.active_plot_index]["original_ppms"][
                            self.intensity_index_initial : self.intensity_index_final
                        ][max_intensity_index]
                        + self.values_dictionary[self.active_plot_index][
                            "move left/right"
                        ]
                    )
                    mean_intensity = np.mean(
                        np.real(
                            self.values_dictionary[self.active_plot_index][
                                "original_data"
                            ][self.intensity_index_initial : self.intensity_index_final]
                            * self.values_dictionary[self.active_plot_index][
                                "multiply value"
                            ]
                        )
                        + self.values_dictionary[self.active_plot_index]["move up/down"]
                        * np.ones(
                            len(
                                self.values_dictionary[self.active_plot_index][
                                    "original_data"
                                ][
                                    self.intensity_index_initial : self.intensity_index_final
                                ]
                            )
                        )
                    )
                    integral = np.sum(
                        np.real(
                            self.values_dictionary[self.active_plot_index][
                                "original_data"
                            ][self.intensity_index_initial : self.intensity_index_final]
                            * self.values_dictionary[self.active_plot_index][
                                "multiply value"
                            ]
                        )
                        + self.values_dictionary[self.active_plot_index]["move up/down"]
                        * np.ones(
                            len(
                                self.values_dictionary[self.active_plot_index][
                                    "original_data"
                                ][
                                    self.intensity_index_initial : self.intensity_index_final
                                ]
                            )
                        )
                    )
                    stdev = np.std(
                        np.real(
                            self.values_dictionary[self.active_plot_index][
                                "original_data"
                            ][self.intensity_index_initial : self.intensity_index_final]
                            * self.values_dictionary[self.active_plot_index][
                                "multiply value"
                            ]
                        )
                        + self.values_dictionary[self.active_plot_index]["move up/down"]
                        * np.ones(
                            len(
                                self.values_dictionary[self.active_plot_index][
                                    "original_data"
                                ][
                                    self.intensity_index_initial : self.intensity_index_final
                                ]
                            )
                        )
                    )
                    min_ppm = (
                        min(
                            self.values_dictionary[self.active_plot_index][
                                "original_ppms"
                            ][self.intensity_index_initial : self.intensity_index_final]
                        )
                        + self.values_dictionary[self.active_plot_index][
                            "move left/right"
                        ]
                    )
                    max_ppm = (
                        max(
                            self.values_dictionary[self.active_plot_index][
                                "original_ppms"
                            ][self.intensity_index_initial : self.intensity_index_final]
                        )
                        + self.values_dictionary[self.active_plot_index][
                            "move left/right"
                        ]
                    )
                    message = "PPM of Max Intensity:\n{:.3f}\nMaximum Intensity:\n{:E}\nMean Intensity:\n{:E}\nIntegral:\n{:E}\nStandard deviation:\n{:E}\nSelected PPM Range:\n{:.3f}-{:.3f} ".format(
                        max_intensity_ppm,
                        max_intensity,
                        mean_intensity,
                        integral,
                        stdev,
                        min_ppm,
                        max_ppm,
                    )
                except:
                    message = "No NMR data found in selected chemical shift range for current selected plot"

            # Message box showing the max intensity
            wx.MessageBox(message, "Max Intensity", wx.OK | wx.ICON_INFORMATION)
        else:
            # Finding the index of the ppms in the region selected
            self.intensity_index_initial = np.abs(
                self.ppms - self.intensity_x_final
            ).argmin()
            self.intensity_index_final = np.abs(
                self.ppms - self.intensity_x_initial
            ).argmin()
            # Message box showing the max intensity
            try:
                max_intensity = max(
                    np.real(
                        self.data[
                            self.intensity_index_initial : self.intensity_index_final
                        ]
                        * self.multiply_value
                    )
                    + self.vertical_slider.GetValue()
                    * np.ones(
                        len(
                            self.data[
                                self.intensity_index_initial : self.intensity_index_final
                            ]
                        )
                    )
                )
                max_intensity_index = np.argmax(
                    np.real(
                        self.data[
                            self.intensity_index_initial : self.intensity_index_final
                        ]
                        * self.multiply_value
                    )
                    + self.vertical_slider.GetValue()
                    * np.ones(
                        len(
                            self.data[
                                self.intensity_index_initial : self.intensity_index_final
                            ]
                        )
                    )
                )
                max_intensity_ppm = self.ppms[
                    self.intensity_index_initial : self.intensity_index_final
                ][max_intensity_index]
                mean_intensity = np.mean(
                    np.real(
                        self.data[
                            self.intensity_index_initial : self.intensity_index_final
                        ]
                        * self.multiply_value
                    )
                    + self.vertical_slider.GetValue()
                    * np.ones(
                        len(
                            self.data[
                                self.intensity_index_initial : self.intensity_index_final
                            ]
                        )
                    )
                )
                integral = np.sum(
                    np.real(
                        self.data[
                            self.intensity_index_initial : self.intensity_index_final
                        ]
                        * self.multiply_value
                    )
                    + self.vertical_slider.GetValue()
                    * np.ones(
                        len(
                            self.data[
                                self.intensity_index_initial : self.intensity_index_final
                            ]
                        )
                    )
                )
                stdev = np.std(
                    np.real(
                        self.data[
                            self.intensity_index_initial : self.intensity_index_final
                        ]
                        * self.multiply_value
                    )
                    + self.vertical_slider.GetValue()
                    * np.ones(
                        len(
                            self.data[
                                self.intensity_index_initial : self.intensity_index_final
                            ]
                        )
                    )
                )
                min_ppm = min(
                    self.ppms[self.intensity_index_initial : self.intensity_index_final]
                )
                max_ppm = max(
                    self.ppms[self.intensity_index_initial : self.intensity_index_final]
                )

                # If on Varian, try to find the frequency in Hz too
                try:
                    dic_v, data_v = ng.varian.read("./")
                    # getting ppm values for the offsets used
                    self.tof = float(dic_v["procpar"]["tof"]["values"][0])
                    # getting the sfrq
                    self.sfrq = float(dic_v["procpar"]["sfrq"]["values"][0])
                    # From the fid.com file finding the carrier
                    file = open("fid.com", "r")
                    fid_com = file.readlines()
                    for line in fid_com:
                        if "CAR" in line:
                            line = line.split("\n")[0].split()
                            del line[0]
                            # deleting the last element of the list which is the '\' character
                            del line[-1]
                            self.carrier = float(line[0])

                    def find_Hz(ppm):
                        Hz = (ppm - self.carrier) * self.sfrq + self.tof
                        return Hz

                    max_intensity_Hz = find_Hz(max_intensity_ppm)
                    min_Hz = find_Hz(min_ppm)
                    max_Hz = find_Hz(max_ppm)
                    wx.MessageBox(
                        "Location of Max Intensity (ppm/Hz):\n{:.3f}/{:.3f}\nMaximum Intensity:\n{:E}\nMean Intensity:{:E}\nIntegral:\n{:E}\nStandard deviation:\n{:E}\nSelected PPM Range:\n{:.3f}-{:.3f}\nDifference(Hz)\n{:.3f}".format(
                            max_intensity_ppm,
                            max_intensity_Hz,
                            max_intensity,
                            mean_intensity,
                            integral,
                            stdev,
                            min_ppm,
                            max_ppm,
                            np.abs(max_Hz - min_Hz),
                        ),
                        "Max Intensity",
                        wx.OK | wx.ICON_INFORMATION,
                    )

                except:
                    wx.MessageBox(
                        "PPM of Max Intensity:\n{:.3f}\nMaximum Intensity:\n{:E}\nMean Intensity:{:E}\nIntegral:\n{:E}\nStandard deviation:\n{:E}\nSelected PPM Range:\n{:.3f}-{:.3f}".format(
                            max_intensity_ppm,
                            max_intensity,
                            mean_intensity,
                            integral,
                            stdev,
                            min_ppm,
                            max_ppm,
                        ),
                        "Max Intensity",
                        wx.OK | wx.ICON_INFORMATION,
                    )
            except:
                wx.MessageBox(
                    "No NMR data found in selected chemical shift range",
                    "Max Intensity",
                    wx.OK | wx.ICON_INFORMATION,
                )

    def OnReleaseNoise(self, event):
        if event.inaxes == self.ax:
            self.release_intensity(event)

    def OnPivotButton(self, event):
        # Getting the user to select a pivot point for phasing by clicking on the spectrum
        wx.MessageBox(
            "Click on the spectrum to set the location of the pivot point for P1 phasing.",
            "Pivot Point",
            wx.OK | wx.ICON_INFORMATION,
        )
        self.pivot_press = self.canvas.mpl_connect(
            "button_press_event", self.OnPivotClick
        )

    def OnPivotClick(self, event):
        # Function to get the x value of the pivot point for phasing
        self.pivot_x = event.xdata
        self.pivot_line.set_xdata([self.pivot_x])

        # Finding the index of the point closest to the pivot point
        self.pivot_index = np.abs(self.ppm_original - self.pivot_x).argmin()
        self.pivot_x = self.pivot_index
        self.canvas.mpl_disconnect(self.pivot_press)
        self.pivot_line.set_visible(True)
        self.OnSliderScroll1D(wx.EVT_SCROLL)

    def OnPivotLoad(self, pivot_x):
        # Function to load the pivot point from a saved session
        self.pivot_x = pivot_x
        self.pivot_line.set_xdata([self.pivot_x])

    def OnRemovePivotButton(self, event):
        if self.pivot_line.get_visible() != True:
            # Give a message saying there is no pivot point to remove
            wx.MessageBox(
                "There is no pivot point to remove.",
                "Remove Pivot Point",
                wx.OK | wx.ICON_INFORMATION,
            )
        else:
            # Function to remove the pivot point for phasing
            self.pivot_x = self.pivot_x_default
            self.pivot_line.set_visible(False)
            self.OnSliderScroll1D(wx.EVT_SCROLL)

    def OnColourChoice1D(self, event):
        # Function to change the colour of the 1D spectrum when a user selects a new colour from the combobox
        self.index = self.colour_chooser.GetSelection()
        self.colour_value = self.colours[self.index]
        if self.multiplot_mode == True:
            self.values_dictionary[self.active_plot_index]["color index"] = self.index
            if self.active_plot_index == 0:
                self.line1.set_color(self.colour_value)
            else:
                self.extra_plots[int(self.active_plot_index) - 1][0].set_color(
                    self.colour_value
                )
            self.ax.legend()
        else:
            self.line1.set_color(self.colour_value)
        self.UpdateFrame()

    def OnLinewidthScroll1D(self, event):
        # Function to change the linewidth of the 1D spectrum
        linewidth_value = float(self.linewidth_slider.GetValue())
        if self.multiplot_mode == True:
            if self.select_all_checkbox.GetValue() == False:
                self.values_dictionary[self.active_plot_index][
                    "linewidth"
                ] = linewidth_value
                if self.active_plot_index == 0:
                    self.line1.set_linewidth(linewidth_value)
                else:
                    self.extra_plots[int(self.active_plot_index) - 1][0].set_linewidth(
                        linewidth_value
                    )
            else:
                self.values_dictionary[0]["linewidth"] = linewidth_value
                self.line1.set_linewidth(linewidth_value)
                for i in range(len(self.extra_plots)):
                    self.values_dictionary[i + 1]["linewidth"] = linewidth_value
                    self.extra_plots[i][0].set_linewidth(linewidth_value)

            self.ax.legend()
        else:
            self.line1.set_linewidth(linewidth_value)
        self.UpdateFrame()

    def OnMultiplyScroll1D(self, event):
        # Function to multiply the 1D spectrum by a constant value
        self.multiply_value = float(self.multiply_slider.GetValue())
        self.multiply_label_value.SetLabel("{:.3f}".format(self.multiply_value))
        if self.multiplot_mode == True:
            if self.select_all_checkbox.GetValue() == False:
                self.values_dictionary[self.active_plot_index][
                    "multiply value"
                ] = self.multiply_value
            else:
                for i in range(len(self.values_dictionary)):
                    self.values_dictionary[i]["multiply value"] = self.multiply_value

        if self.multiplot_mode == True:
            if self.select_all_checkbox.GetValue() == False:
                self.y_data = self.values_dictionary[self.active_plot_index][
                    "original_data"
                ]
                if self.active_plot_index == 0:
                    self.line1.set_ydata(
                        (
                            self.y_data
                            + self.vertical_slider.GetValue()
                            * np.ones(len(self.y_data))
                        )
                        * self.multiply_value
                    )
                else:
                    self.extra_plots[int(self.active_plot_index) - 1][0].set_ydata(
                        (
                            self.y_data
                            + self.vertical_slider.GetValue()
                            * np.ones(len(self.y_data))
                        )
                        * self.multiply_value
                    )
            else:
                self.line1.set_ydata(
                    (
                        self.data
                        + self.values_dictionary[0]["move up/down"]
                        * np.ones(len(self.data))
                    )
                    * self.multiply_value
                )
                for i in range(len(self.extra_plots)):
                    self.extra_plots[i][0].set_ydata(
                        (
                            self.data
                            + self.values_dictionary[i + 1]["move up/down"]
                            * np.ones(len(self.data))
                        )
                        * self.multiply_value
                    )
        else:
            self.line1.set_ydata(
                (self.data + self.vertical_slider.GetValue() * np.ones(len(self.data)))
                * self.multiply_value
            )
        self.OnSliderScroll1D(wx.EVT_SCROLL)
        self.UpdateFrame()

    def OnMultiplyCombo(self, event):
        self.multiply_index = int(self.multiply_range_chooser.GetSelection())
        if self.multiplot_mode == True:
            if self.select_all_checkbox.GetValue() == False:
                self.values_dictionary[self.active_plot_index][
                    "multiply range index"
                ] = self.multiply_index
            else:
                for i in range(len(self.values_dictionary)):
                    self.values_dictionary[i][
                        "multiply range index"
                    ] = self.multiply_index
        self.multiply_range = float(self.multiply_range_values[self.multiply_index])
        self.multiply_slider.SetRange(0.1, self.multiply_range)
        self.multiply_slider.SetRes(self.multiply_range / 1000)
        self.multiply_slider.Bind(wx.EVT_SLIDER, self.OnMultiplyScroll1D)

    def OnReferenceScroll1D(self, event):
        # Function to move the spectrum left/right in the ppm scale when the slider position is changed
        reference_value = float(self.reference_slider.GetValue())
        self.reference_value_label.SetLabel("{:.4f}".format(reference_value))

        if self.multiplot_mode == False:
            self.ppms = (
                self.ppm_original + np.ones(len(self.ppm_original)) * reference_value
            )
            self.line1.set_xdata(self.ppms)
        else:
            if self.select_all_checkbox.GetValue() == False:
                self.values_dictionary[self.active_plot_index][
                    "move left/right"
                ] = reference_value
                if self.active_plot_index == 0:
                    self.ppm_original = self.values_dictionary[self.active_plot_index][
                        "original_ppms"
                    ]
                    self.ppms = (
                        self.ppm_original
                        + np.ones(len(self.ppm_original)) * reference_value
                    )
                    self.line1.set_xdata(self.ppms)
                else:
                    self.ppm_original = self.values_dictionary[self.active_plot_index][
                        "original_ppms"
                    ]
                    self.ppms = (
                        self.ppm_original
                        + np.ones(len(self.ppm_original)) * reference_value
                    )

                    self.extra_plots[int(self.active_plot_index) - 1][0].set_xdata(
                        self.ppms
                    )
            else:
                for i in range(len(self.values_dictionary)):
                    self.values_dictionary[i]["move left/right"] = reference_value
                    self.ppm_original = self.values_dictionary[i]["original_ppms"]
                    self.ppms = (
                        self.ppm_original
                        + np.ones(len(self.ppm_original)) * reference_value
                    )
                    if i == 0:
                        self.line1.set_xdata(self.ppms)
                    else:
                        self.extra_plots[i - 1][0].set_xdata(self.ppms)
        self.UpdateFrame()

    def OnVerticalScroll1D(self, event):
        if self.multiplot_mode == True:
            if self.select_all_checkbox.GetValue() == False:
                self.values_dictionary[self.active_plot_index]["move up/down"] = float(
                    self.vertical_slider.GetValue()
                )
            else:
                for i in range(len(self.values_dictionary)):
                    self.values_dictionary[i]["move up/down"] = float(
                        self.vertical_slider.GetValue()
                    )

        self.OnSliderScroll1D(wx.EVT_SCROLL)
        self.UpdateFrame()

    def OnReferenceCombo(self, event):
        # Function to change the slider limits for the move left/right slider
        self.ref_index = int(self.reference_range_chooser.GetSelection())
        if self.multiplot_mode == True:
            if self.select_all_checkbox.GetValue() == False:
                self.values_dictionary[self.active_plot_index][
                    "move left/right range index"
                ] = self.ref_index
            else:
                for i in range(len(self.values_dictionary)):
                    self.values_dictionary[i][
                        "move left/right range index"
                    ] = self.ref_index
        self.reference_range = float(self.reference_range_values[self.ref_index])
        self.reference_slider.SetRange(-self.reference_range, self.reference_range)
        self.reference_slider.SetRes(self.reference_range / 1000)
        self.reference_slider.Bind(wx.EVT_SLIDER, self.OnReferenceScroll1D)

    def OnVerticalCombo(self, event):
        # Function to change the slider limits for the vertical shift slider
        self.vertical_index = int(self.vertical_range_chooser.GetSelection())
        if self.multiplot_mode == True:
            if self.select_all_checkbox.GetValue() == False:
                self.values_dictionary[self.active_plot_index][
                    "move up/down range index"
                ] = self.vertical_index
            else:
                for i in range(len(self.values_dictionary)):
                    self.values_dictionary[i][
                        "move up/down range index"
                    ] = self.vertical_index
        self.vertical_percentage = float(
            self.vertical_range_values[self.vertical_index]
        )
        self.vertical_slider.SetRange(
            -self.vertical_range * self.vertical_percentage / 100,
            self.vertical_range * self.vertical_percentage / 100,
        )
        self.vertical_slider.SetRes(
            self.vertical_range * self.vertical_percentage / 10000
        )
        self.vertical_slider.Bind(wx.EVT_SLIDER, self.OnVerticalScroll1D)

    def OnSubtractButton(self, event):
        # If not in multiplot mode, create pop up window saying that there is only one spectrum loaded so subtraction is not possible
        if self.multiplot_mode == False:
            msg = "Only one spectrum loaded. Subtraction not possible"
            dlg = wx.MessageDialog(None, msg, "Error", wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
            return
        else:
            titles = []
            for i in range(len(self.values_dictionary)):
                titles.append(self.values_dictionary[i]["title"])
            self.subraction_selection_frame = wx.Frame(
                self, title="Subtraction Selection", size=(400, 200)
            )
            self.subtraction_sizer_total = wx.BoxSizer(wx.HORIZONTAL)
            self.subtraction_sizer_main = wx.BoxSizer(wx.VERTICAL)
            self.subtraction_sizer1 = wx.BoxSizer(wx.HORIZONTAL)
            self.subtraction_sizer2 = wx.BoxSizer(wx.HORIZONTAL)
            self.subtraction_sizer3 = wx.BoxSizer(wx.HORIZONTAL)
            self.spectrum_subraction_text = wx.StaticText(
                self.subraction_selection_frame,
                label="Input spectra to be subtracted (Spectrum 1 - Spectrum 2):",
            )
            self.spectrum1_label = wx.StaticText(
                self.subraction_selection_frame, label="Spectrum 1:"
            )
            self.spectrum2_label = wx.StaticText(
                self.subraction_selection_frame, label="Spectrum 2:"
            )
            self.subracted_spectrum_label = wx.StaticText(
                self.subraction_selection_frame, label="Name of subtracted spectrum:"
            )
            self.spectrum1_combobox = wx.ComboBox(
                self.subraction_selection_frame, choices=titles, style=wx.CB_READONLY
            )
            self.spectrum1_combobox.SetSelection(1)
            self.spectrum2_combobox = wx.ComboBox(
                self.subraction_selection_frame, choices=titles, style=wx.CB_READONLY
            )
            self.spectrum2_combobox.SetSelection(0)
            self.subracted_spectrum_text = wx.TextCtrl(
                self.subraction_selection_frame,
                size=(100, self.spectrum1_combobox.GetSize().GetHeight()),
            )

            self.subtraction_sizer1.Add(self.spectrum1_label)
            self.subtraction_sizer1.AddSpacer(5)
            self.subtraction_sizer1.Add(self.spectrum1_combobox)
            self.subtraction_sizer2.Add(self.spectrum2_label)
            self.subtraction_sizer2.AddSpacer(5)
            self.subtraction_sizer2.Add(self.spectrum2_combobox)
            self.subtraction_sizer3.Add(self.subracted_spectrum_label)
            self.subtraction_sizer3.AddSpacer(5)
            self.subtraction_sizer3.Add(self.subracted_spectrum_text)

            self.subtraction_sizer_main.Add(
                self.spectrum_subraction_text, 0, wx.ALIGN_CENTER_HORIZONTAL
            )
            self.subtraction_sizer_main.AddSpacer(5)
            self.subtraction_sizer_main.Add(
                self.subtraction_sizer1, 0, wx.ALIGN_CENTER_HORIZONTAL
            )
            self.subtraction_sizer_main.AddSpacer(5)
            self.subtraction_sizer_main.Add(
                self.subtraction_sizer2, 0, wx.ALIGN_CENTER_HORIZONTAL
            )
            self.subtraction_sizer_main.AddSpacer(5)
            self.subtraction_sizer_main.Add(
                self.subtraction_sizer3, 0, wx.ALIGN_CENTER_HORIZONTAL
            )
            self.subtraction_sizer_main.AddSpacer(10)

            self.subtract_button = wx.Button(
                self.subraction_selection_frame, label="Subtract"
            )
            self.subtract_button.Bind(wx.EVT_BUTTON, self.OnSubtractSpectra)
            self.subtraction_sizer_main.Add(
                self.subtract_button, 0, wx.ALIGN_CENTER_HORIZONTAL
            )

            self.subtraction_sizer_total.Add(
                self.subtraction_sizer_main, 5, wx.ALIGN_CENTER
            )

            self.subraction_selection_frame.SetSizer(self.subtraction_sizer_main)
            self.subraction_selection_frame.Show()

    def OnSubtractSpectra(self, event):
        # Get the index of the spectra to be subtracted
        spectrum1_index = self.spectrum1_combobox.GetSelection()
        spectrum2_index = self.spectrum2_combobox.GetSelection()
        # Get the current state of the data to be subtracted (including any baseline subtraction, multiplication and movement in the ppm scale)
        spectrum1_data = self.values_dictionary[spectrum1_index]["original_data"]
        spectrum2_data = self.values_dictionary[spectrum2_index]["original_data"]
        spectrum1_ppms = self.values_dictionary[spectrum1_index]["original_ppms"]
        spectrum2_ppms = self.values_dictionary[spectrum2_index]["original_ppms"]
        spectrum1_vertical = self.values_dictionary[spectrum1_index]["move up/down"]
        spectrum2_vertical = self.values_dictionary[spectrum2_index]["move up/down"]
        spectrum1_reference = self.values_dictionary[spectrum1_index]["move left/right"]
        spectrum2_reference = self.values_dictionary[spectrum2_index]["move left/right"]
        spectrum1_multiply = self.values_dictionary[spectrum1_index]["multiply value"]
        spectrum2_multiply = self.values_dictionary[spectrum2_index]["multiply value"]

        modified_spectrum1_data = (
            spectrum1_data * spectrum1_multiply
            + np.ones(len(spectrum1_data)) * spectrum1_vertical
        )
        modified_spectrum2_data = (
            spectrum2_data * spectrum2_multiply
            + np.ones(len(spectrum2_data)) * spectrum2_vertical
        )
        modified_spectrum1_ppms = (
            spectrum1_ppms + np.ones(len(spectrum1_ppms)) * spectrum1_reference
        )
        modified_spectrum2_ppms = (
            spectrum2_ppms + np.ones(len(spectrum2_ppms)) * spectrum2_reference
        )

        # Find the overlapping ppms between the two spectra
        min_ppms = max(modified_spectrum1_ppms[-1], modified_spectrum2_ppms[-1])
        max_ppms = min(modified_spectrum1_ppms[0], modified_spectrum2_ppms[0])

        # Get the index of all modified_spectrum1_ppms and modified_spectrum2_ppms that are within the overlapping range
        spectrum1_index_initial = np.abs(modified_spectrum1_ppms - max_ppms).argmin()
        spectrum1_index_final = np.abs(modified_spectrum1_ppms - min_ppms).argmin()
        spectrum2_index_initial = np.abs(modified_spectrum2_ppms - max_ppms).argmin()
        spectrum2_index_final = np.abs(modified_spectrum2_ppms - min_ppms).argmin()

        # Get the data for the common ppm values
        common_ppms = modified_spectrum1_ppms[
            spectrum1_index_initial:spectrum1_index_final
        ]
        common_ppms_2 = modified_spectrum2_ppms[
            spectrum2_index_initial:spectrum2_index_final
        ]
        common_spectrum1_data = modified_spectrum1_data[
            spectrum1_index_initial:spectrum1_index_final
        ]
        common_spectrum2_data = modified_spectrum2_data[
            spectrum2_index_initial:spectrum2_index_final
        ]
        self.subtracted_ppms = common_ppms

        if len(common_ppms) == 0:
            # Give a message box saying that the spectra do not overlap
            message = "The spectra do not overlap. Subtraction not possible"
            dlg = wx.MessageDialog(self, message, "Error", wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
            return

        # Check to see if the length of common_spectrum1_data and common_spectrum2_data are the same
        if len(common_spectrum1_data) != len(common_spectrum2_data):
            # Give a message box saying that the spectra are not the same length. Ask the user if they want to interpolate the data to the same length
            message = "The spectra are not the same length. Do you want to interpolate the data to the same length?"
            dlg = wx.MessageDialog(
                self, message, "Interpolate Data", wx.YES_NO | wx.ICON_QUESTION
            )
            result = dlg.ShowModal()
            if result == wx.ID_YES:
                # Interpolate the data to the same length
                if len(common_spectrum1_data) > len(common_spectrum2_data):
                    common_spectrum2_data = np.interp(
                        np.flip(common_ppms),
                        np.flip(common_ppms_2),
                        np.flip(common_spectrum2_data),
                    )
                    common_spectrum2_data = np.flip(common_spectrum2_data)
                    self.subtracted_ppms = common_ppms
                else:
                    common_spectrum1_data = np.interp(
                        np.flip(common_ppms_2),
                        np.flip(common_ppms),
                        np.flip(common_spectrum1_data),
                    )
                    common_spectrum1_data = np.flip(common_spectrum1_data)
                    self.subtracted_ppms = common_ppms_2
            else:
                return

        # subtract the two spectra
        subtracted_data = common_spectrum1_data - common_spectrum2_data

        # Save the subtracted data to the same directory as the original spectra
        obs = float(self.nmrdata.dic["FDF2OBS"])
        # Set the carrier to the middle of the ppm range
        car = (max(self.subtracted_ppms) + min(self.subtracted_ppms)) / 2
        # Set the size to the length of the subtracted data
        size = len(subtracted_data)
        # Set the label
        label = self.nmrdata.dic["FDF2LABEL"]
        # Set the sweep width to the difference between the maximum and minimum ppm values
        sw = (max(self.subtracted_ppms) - min(self.subtracted_ppms)) * obs
        udic = {
            "ndim": 1,
            0: {
                "size": size,
                "complex": False,
                "encoding": "int",
                "sw": sw,
                "obs": obs,
                "car": car,
                "label": label,
                "time": False,
                "freq": True,
            },
        }
        dic = ng.pipe.create_dic(udic)
        dic["FDF2OBS"] = obs
        dic["FDF2CAR"] = car
        dic["FDF2SIZE"] = size
        dic["FDF2LABEL"] = label
        dic["FDF2SW"] = sw
        orig = min(self.subtracted_ppms) * obs
        center = ((max(self.subtracted_ppms) + min(self.subtracted_ppms)) / 2) * obs
        dic["FDF2ORIG"] = orig
        dic["FDF2CENTER"] = center

        cwd = os.getcwd()
        try:
            if self.parent.path != "":
                os.chdir(self.parent.path)
        except:
            # Change to the path of the original spectra
            path = self.values_dictionary[0]["path"].split("/")[0:-1]
            path = "/".join(path)
            os.chdir(path)
        subtracted_data_32 = subtracted_data.astype(np.float32)
        ng.pipe.write(
            self.subracted_spectrum_text.GetValue() + ".ft",
            dic,
            subtracted_data_32,
            overwrite=True,
        )

        try:
            if self.parent.cwd != "":
                os.chdir(self.parent.cwd)
        except:
            os.chdir(cwd)

        new_spectrum_name = self.subracted_spectrum_text.GetValue()
        # Add the new spectrum to the dictionary of spectra
        path = self.values_dictionary[0]["path"].split("/")[0:-1]
        path = "/".join(path)
        self.values_dictionary[len(self.values_dictionary.keys())] = {
            "title": new_spectrum_name,
            "original_data": subtracted_data,
            "original_ppms": modified_spectrum2_ppms,
            "move up/down": 0,
            "move left/right": 0,
            "multiply value": 1,
            "color index": len(self.extra_plots) + 1,
            "move up/down range index": 0,
            "move left/right range index": 0,
            "multiply range index": 0,
            "linewidth": 0.5,
            "p0 Coarse": 0,
            "p0 Fine": 0,
            "p1 Coarse": 0,
            "p1 Fine": 0,
            "path": path + "/" + new_spectrum_name + ".ft",
        }

        # Plot the new spectrum
        self.extra_plots.append(
            self.ax.plot(
                self.subtracted_ppms,
                subtracted_data,
                color=self.files.color_list[
                    len(self.extra_plots) + 1 - len(self.files.color_list)
                ],
                label=new_spectrum_name,
                linewidth=0.5,
            )
        )

        # # Input the subtracted spectrum into the values_dictionary
        # self.values_dictionary[len(self.values_dictionary.keys())]['color index'] = len(self.extra_plots)
        # self.values_dictionary[len(self.values_dictionary.keys())]['linewidth'] = 0.5
        # self.values_dictionary[len(self.values_dictionary.keys())]['p0 Coarse'] = 0
        # self.values_dictionary[len(self.values_dictionary.keys())]['p0 Fine'] = 0
        # self.values_dictionary[len(self.values_dictionary.keys())]['p1 Coarse'] = 0
        # self.values_dictionary[len(self.values_dictionary.keys())]['p1 Fine'] = 0
        # self.values_dictionary[len(self.values_dictionary.keys())]['move up/down range index'] = 0
        # self.values_dictionary[len(self.values_dictionary.keys())]['move left/right range index'] = 0
        # self.values_dictionary[len(self.values_dictionary.keys())]['multiply value index'] = 0
        # self.values_dictionary[len(self.values_dictionary.keys())]['original_data'] = subtracted_data
        # self.values_dictionary[len(self.values_dictionary.keys())]['original_ppms'] = self.subtracted_ppms
        # self.values_dictionary[len(self.values_dictionary.keys())]['move up/down'] = 0
        # self.values_dictionary[len(self.values_dictionary.keys())]['move left/right'] = 0
        # self.values_dictionary[len(self.values_dictionary.keys())]['multiply value'] = 1
        # self.values_dictionary[len(self.values_dictionary.keys())]['title'] = new_spectrum_name
        # self.values_dictionary[len(self.values_dictionary.keys())]['path'] = file_path

        # Add labels of the extra plots to the select plot box
        self.plot_labels = self.plot_combobox.GetItems()
        self.plot_labels.append(new_spectrum_name)
        self.plot_combobox.Clear()
        self.plot_combobox.AppendItems(self.plot_labels)
        self.plot_combobox.SetSelection(0)
        self.ax.legend()
        self.UpdateFrame()

        self.subraction_selection_frame.Destroy()

        # Save the new spectrum to the same directory as the original spectra

    def OnSaveButton(self, event):
        if self.multiplot_mode == False:
            # Have a popout asking the user what the want the spectrum to be saved as
            msg = "Input a name for the spectrum to be saved as"
            dlg = wx.TextEntryDialog(None, msg)
            res = dlg.ShowModal()
            spectrum_name = dlg.GetValue()
            dlg.Destroy()
            # Get the current state of the data to be saved
            data = self.line1.get_ydata()
            data_float32 = data.astype(np.float32)
            dic_new = self.nmrdata.dic
            car_old = self.nmrdata.dic["FDF2CAR"]
            car = float(car_old) + float(self.reference_slider.GetValue())
            obs = self.nmrdata.dic["FDF2OBS"]
            orig = self.nmrdata.dic["FDF2ORIG"]
            center = self.nmrdata.dic["FDF2CENTER"]

            orig = float(orig) + float(self.reference_slider.GetValue()) * float(obs)
            center = float(center) + float(self.reference_slider.GetValue()) * float(
                obs
            )
            dic_new["FDF2CAR"] = car
            dic_new["FDF2ORIG"] = orig
            dic_new["FDF2CENTER"] = center

            if self.parent.path != "":
                os.chdir(self.parent.path)
            ng.pipe.write(spectrum_name + ".ft", dic_new, data_float32, overwrite=True)
            if self.parent.cwd != "":
                os.chdir(self.parent.cwd)

        else:
            # Have a popout asking the user to pick from a combobox which spectrum they want to save
            titles = []
            for i in range(len(self.values_dictionary.keys())):
                titles.append(self.values_dictionary[i]["title"])
            self.save_frame = wx.Frame(self, title="Save Spectrum", size=(400, 200))
            self.save_sizer_main = wx.BoxSizer(wx.VERTICAL)
            self.save_sizer1 = wx.BoxSizer(wx.HORIZONTAL)
            self.save_sizer2 = wx.BoxSizer(wx.HORIZONTAL)
            self.save_text = wx.StaticText(
                self.save_frame, label="Select spectrum to be saved:"
            )
            self.save_combobox = wx.ComboBox(
                self.save_frame, choices=titles, style=wx.CB_READONLY
            )
            self.save_combobox.SetSelection(0)
            self.save_textcontrol_text = wx.StaticText(
                self.save_frame, label="Input name for spectrum to be saved as:"
            )
            self.save_textcontrol = wx.TextCtrl(
                self.save_frame, size=(100, self.save_combobox.GetSize().GetHeight())
            )
            self.save_button = wx.Button(self.save_frame, label="Save")
            self.save_button.Bind(wx.EVT_BUTTON, self.OnSaveSpectrum)
            self.save_sizer1.Add(self.save_text)
            self.save_sizer1.AddSpacer(5)
            self.save_sizer1.Add(self.save_combobox)
            self.save_sizer2.Add(self.save_textcontrol_text)
            self.save_sizer2.AddSpacer(5)
            self.save_sizer2.Add(self.save_textcontrol)
            self.save_sizer_main.AddSpacer(10)
            self.save_sizer_main.Add(self.save_sizer1, 0, wx.ALIGN_CENTER_HORIZONTAL)
            self.save_sizer_main.AddSpacer(10)
            self.save_sizer_main.Add(self.save_sizer2, 0, wx.ALIGN_CENTER_HORIZONTAL)
            self.save_sizer_main.AddSpacer(10)
            self.save_sizer_main.Add(self.save_button, 0, wx.ALIGN_CENTER_HORIZONTAL)
            self.save_frame.SetSizer(self.save_sizer_main)
            self.save_frame.Show()

    def OnSaveSpectrum(self, event):
        # Get the index of the spectrum to be saved
        spectrum_index = self.save_combobox.GetSelection()
        # Get the name of the spectrum to be saved as
        spectrum_name = self.save_textcontrol.GetValue()
        # Get the current state of the data to be saved
        if spectrum_index == 0:
            data = self.line1.get_ydata()
        else:
            data = self.extra_plots[spectrum_index - 1][0].get_ydata()
        data_float32 = data.astype(np.float32)

        # See if data has come from stackmode or not
        dic = self.nmrdata.dic

        if float(len(data_float32)) == dic["FDF1TDSIZE"]:
            obs = float(self.nmrdata.dic["FDF1OBS"])
            # Set the carrier to the middle of the ppm range
            car = self.nmrdata.dic["FDF1CAR"]
            # Set the size to the length of the subtracted data
            size = len(data_float32)
            # Set the label
            label = self.nmrdata.dic["FDF1LABEL"]
            # Set the sweep width to the difference between the maximum and minimum ppm values
            sw = self.nmrdata.dic["FDF1SW"]
            orig = self.nmrdata.dic["FDF1ORIG"]
            center = self.nmrdata.dic["FDF1CENTER"]
            udic = {
                "ndim": 1,
                0: {
                    "size": size,
                    "complex": False,
                    "encoding": "int",
                    "sw": sw,
                    "obs": obs,
                    "car": car,
                    "label": label,
                    "time": False,
                    "freq": True,
                },
            }
            dic = ng.pipe.create_dic(udic)
            dic["FDF2OBS"] = obs
            dic["FDF2CAR"] = car
            dic["FDF2SIZE"] = size
            dic["FDF2LABEL"] = label
            dic["FDF2SW"] = sw
            dic["FDF2ORIG"] = orig
            dic["FDF2CENTER"] = center

        elif float(len(data_float32)) == dic["FDF2TDSIZE"]:

            obs = float(self.nmrdata.dic["FDF2OBS"])
            # Set the carrier to the middle of the ppm range
            car = self.nmrdata.dic["FDF2CAR"]
            # Set the size to the length of the subtracted data
            size = len(data_float32)
            # Set the label
            label = self.nmrdata.dic["FDF2LABEL"]
            # Set the sweep width to the difference between the maximum and minimum ppm values
            sw = self.nmrdata.dic["FDF2SW"]
            orig = self.nmrdata.dic["FDF2ORIG"]
            center = self.nmrdata.dic["FDF2CENTER"]
            udic = {
                "ndim": 1,
                0: {
                    "size": size,
                    "complex": False,
                    "encoding": "int",
                    "sw": sw,
                    "obs": obs,
                    "car": car,
                    "label": label,
                    "time": False,
                    "freq": True,
                },
            }
            dic = ng.pipe.create_dic(udic)
            dic["FDF2OBS"] = obs
            dic["FDF2CAR"] = car
            dic["FDF2SIZE"] = size
            dic["FDF2LABEL"] = label
            dic["FDF2SW"] = sw
            dic["FDF2ORIG"] = orig
            dic["FDF2CENTER"] = center

        ng.pipe.write(spectrum_name + ".ft", dic, data_float32, overwrite=True)

        self.save_frame.Destroy()

    def draw_figure_1D(self):
        # Function to plot the 1D spectrum
        self.ax = self.fig.add_subplot(111)
        # Get ppm values for x axis
        if self.uc0 == None:
            if self.nmrdata.file != ".":
                self.uc0 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data)
            else:
                udic = ng.bruker.guess_udic(self.nmrdata.dic, self.nmrdata.data)
                self.uc0 = ng.fileiobase.uc_from_udic(udic)

        self.ppm_original = self.uc0.ppm_scale()
        self.ppms = self.ppm_original
        self.data = self.nmrdata.data
        (self.line1,) = self.ax.plot(self.ppms, self.data, linewidth=0.5)
        self.ax.set_xlabel(self.nmrdata.axislabels[0])
        self.ax.set_ylabel("Intensity")
        self.ax.set_xlim(max(self.ppms), min(self.ppms))
        self.line1.set_color(self.colour_value)

        # Create a pivot line and set to invisible
        self.pivot_line = self.ax.axvline(
            self.pivot_x_default, color="black", linestyle="--"
        )
        self.pivot_line.set_visible(False)

        # Create an intensity region and set to invisible
        self.intensity_region = self.ax.axvspan(
            min(self.ppm_original), min(self.ppm_original), alpha=0.2, color="gray"
        )
        self.intensity_region.set_visible(False)

        self.UpdateFrame()

        self.active_plot = self.line1

        self.files = FileDrop(self.canvas, self.ax, self)
        self.canvas.SetDropTarget(self.files)

    def OnSliderScroll1D(self, event):
        # Get all the slider values for P0 and P1 (coarse and fine), put the combined coarse and fine values on the screen
        self.total_P0 = self.P0_slider.GetValue() + self.P0_slider_fine.GetValue()
        self.total_P1 = self.P1_slider.GetValue() + self.P1_slider_fine.GetValue()
        self.P0_total_value.SetLabel("{:.2f}".format(self.total_P0))
        self.P1_total_value.SetLabel("{:.2f}".format(self.total_P1))
        self.phase1D()

    def phase1D(self):
        # Function to phase the data using the combined course/fine phasing values and plot
        self.multiply_value = float(self.multiply_slider.GetValue())
        if self.multiplot_mode == False:
            imaginary_data = ng.process.proc_base.ht(
                self.nmrdata.data, self.nmrdata.data.shape[0]
            )
            self.data = imaginary_data * np.exp(
                1j
                * (
                    self.total_P0 * np.pi / 180
                    + self.total_P1
                    * (np.pi / 180)
                    * (
                        np.arange(
                            -self.pivot_x, -self.pivot_x + self.nmrdata.data.shape[0]
                        )
                        / self.nmrdata.data.shape[0]
                    )
                )
            ) + np.ones(len(self.nmrdata.data)) * float(self.vertical_slider.GetValue())
            if len(self.data_spline) > 1:
                try:
                    self.data = self.data - self.data_spline
                except:
                    print("Baseline subtraction unsuccessful - continuing")
                    pass
            self.line1.set_ydata(
                self.data * self.multiply_value
                + np.ones(len(self.data)) * float(self.vertical_slider.GetValue())
            )
        else:
            if self.select_all_checkbox.GetValue() == False:
                self.values_dictionary[self.active_plot_index][
                    "p0 Coarse"
                ] = self.P0_slider.GetValue()
                self.values_dictionary[self.active_plot_index][
                    "p0 Fine"
                ] = self.P0_slider_fine.GetValue()
                self.values_dictionary[self.active_plot_index][
                    "p1 Coarse"
                ] = self.P1_slider.GetValue()
                self.values_dictionary[self.active_plot_index][
                    "p1 Fine"
                ] = self.P1_slider_fine.GetValue()
                if self.active_plot_index == 0:
                    original_data = self.values_dictionary[self.active_plot_index][
                        "original_data"
                    ]
                    imaginary_data = ng.process.proc_base.ht(
                        original_data, original_data.shape[0]
                    )
                    self.data = imaginary_data * np.exp(
                        1j
                        * (
                            self.total_P0 * np.pi / 180
                            + self.total_P1
                            * (np.pi / 180)
                            * (
                                np.arange(
                                    -self.pivot_x,
                                    -self.pivot_x + original_data.shape[0],
                                )
                                / original_data.shape[0]
                            )
                        )
                    )
                    self.line1.set_ydata(
                        self.data
                        * self.values_dictionary[self.active_plot_index][
                            "multiply value"
                        ]
                        + np.ones(len(self.data))
                        * float(
                            self.values_dictionary[self.active_plot_index][
                                "move up/down"
                            ]
                        )
                    )
                else:
                    original_data = self.values_dictionary[self.active_plot_index][
                        "original_data"
                    ]
                    imaginary_data = ng.process.proc_base.ht(
                        original_data, original_data.shape[0]
                    )
                    self.data = imaginary_data * np.exp(
                        1j
                        * (
                            self.total_P0 * np.pi / 180
                            + self.total_P1
                            * (np.pi / 180)
                            * (
                                np.arange(
                                    -self.pivot_x,
                                    -self.pivot_x + original_data.shape[0],
                                )
                                / original_data.shape[0]
                            )
                        )
                    )
                    self.extra_plots[self.active_plot_index - 1][0].set_ydata(
                        self.data
                        * self.values_dictionary[self.active_plot_index][
                            "multiply value"
                        ]
                        + np.ones(len(self.data))
                        * self.values_dictionary[self.active_plot_index]["move up/down"]
                    )
            else:
                for i in range(len(self.values_dictionary)):
                    original_data = self.values_dictionary[i]["original_data"]
                    self.values_dictionary[i]["p0 Coarse"] = self.P0_slider.GetValue()
                    self.values_dictionary[i][
                        "p0 Fine"
                    ] = self.P0_slider_fine.GetValue()
                    self.values_dictionary[i]["p1 Coarse"] = self.P1_slider.GetValue()
                    self.values_dictionary[i][
                        "p1 Fine"
                    ] = self.P1_slider_fine.GetValue()
                    imaginary_data = ng.process.proc_base.ht(
                        original_data, original_data.shape[0]
                    )
                    self.data = imaginary_data * np.exp(
                        1j
                        * (
                            self.total_P0 * np.pi / 180
                            + self.total_P1
                            * (np.pi / 180)
                            * (
                                np.arange(
                                    -self.pivot_x,
                                    -self.pivot_x + original_data.shape[0],
                                )
                                / original_data.shape[0]
                            )
                        )
                    )
                    if i == 0:
                        self.line1.set_ydata(
                            self.data * self.values_dictionary[i]["multiply value"]
                            + np.ones(len(self.data))
                            * self.values_dictionary[i]["move up/down"]
                        )
                    else:
                        self.extra_plots[i - 1][0].set_ydata(
                            self.data * self.values_dictionary[i]["multiply value"]
                            + np.ones(len(self.data))
                            * self.values_dictionary[i]["move up/down"]
                        )
        self.UpdateFrame()

    def OnIntensityScroll1D(self, event):
        # Function to change the y axis limits
        intensity_percent = 10 ** float(self.intensity_slider.GetValue())

        if self.nmrdata.dim == 1:
            self.ax.set_ylim(
                -(np.max(self.nmrdata.data) / 8) / (intensity_percent / 100),
                np.max(self.nmrdata.data) / (intensity_percent / 100),
            )
            self.UpdateFrame()


# A class to create a panel for viewing 2D NMR spectra
class TwoDViewer(wx.Panel):
    def __init__(self, parent, nmrdata, threeDprojection=False):
        # Get the monitor size and set the window size to 85% of the monitor size
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        self.display_index = wx.Display.GetFromWindow(parent)
        self.width = int(1.0 * sizes[self.display_index][0])
        self.height = int(0.875 * sizes[self.display_index][1])
        self.parent = parent
        self.threeDprojection = threeDprojection
        wx.Panel.__init__(self, parent, id=wx.ID_ANY, size=(self.width, self.height))
        self.nmrdata = nmrdata
        self.set_initial_variables_2D()
        self.create_button_panel_2D()
        self.create_hidden_button_panel_2D()
        self.create_canvas_2D()
        self.add_to_main_sizer_2D()
        self.draw_figure_2D()

    def add_to_main_sizer_2D(self):
        # Create the main sizer
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.main_sizer.Add(self.canvas, 10, wx.EXPAND)
        self.main_sizer.Add(self.toolbar, 0, wx.EXPAND)
        self.main_sizer.Add(self.bottom_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.main_sizer.Add(self.show_button_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.main_sizer.Hide(self.show_button_sizer)
        self.SetSizer(self.main_sizer)

    def create_canvas_2D(self):
        # Create the figure and canvas to draw on
        self.panel = wx.Panel(self)
        self.fig = Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        self.toolbar = NavigationToolbar(self.canvas)

    def create_hidden_button_panel_2D(self):
        # Create a button to show the options
        self.show_button = wx.Button(self, label="Show Options")
        self.show_button.Bind(wx.EVT_BUTTON, self.OnHideButton)
        self.show_button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.show_button_sizer.Add(self.show_button, wx.ALIGN_CENTER, 5)
        self.show_button_sizer.AddSpacer(5)

    def set_initial_variables_2D(self):
        # Colours for 1D lines
        self.colours = [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
        ]
        self.colour_value = self.colours[0]

        # Initial 1D slice colour for 2D/3D spectra is set to navy
        self.colour_slice = "navy"

        # List of cmap colours for when overlaying multiple spectra
        self.cmap = "#e41a1c"
        self.cmap_neg = "#377eb8"
        self.twoD_colours = [
            "#e41a1c",
            "#377eb8",
            "#4daf4a",
            "#984ea3",
            "#ff7f00",
            "#ff33eb",
            "tan",
            "lightcoral",
            "maroon",
            "lightgreen",
            "deeppink",
            "fuchsia",
        ]
        self.twoD_label_colours = self.twoD_colours

        self.twoD_slices_horizontal = []
        self.twoD_slices_vertical = []

        # Range of the sliders to for moving spectra left/right/up/down
        self.reference_range_values = [
            "0.01",
            "0.1",
            "0.5",
            "1.0",
            "5.0",
            "10.0",
            "50.0",
        ]
        self.reference_range = float(self.reference_range_values[0])
        self.reference_rangeX = float(self.reference_range_values[0])
        self.reference_rangeY = float(self.reference_range_values[0])

        # Range of the sliders to for moving spectra up/down in 1D spectra
        self.vertical_range_values = [
            "0.01",
            "0.1",
            "0.5",
            "1.0",
            "10.0",
            "50.0",
            "100.0",
            "1000.0",
            "10000",
        ]

        # Range of the sliders to for multiplying 1D spectra
        self.multiply_range_values = [
            "1.01",
            "1.1",
            "1.5",
            "2",
            "5",
            "10",
            "50",
            "100",
            "1000",
            "10000",
            "100000",
            "1000000",
            "10000000",
            "100000000",
            "1000000000",
        ]

        # Initial x,y movements for referencing are set to zero
        self.x_movement = 0
        self.y_movement = 0

        # Multiplot mode is initially set to off
        self.multiplot_mode = False

        # Dictionary to store the values of the sliders for each spectrum in multiplot mode
        self.values_dictionary = {}

        # Initial multiply factor is 1
        self.multiply_factor = 1

        # 1D slice color of 2D spectra is initially set to green
        self.slice_colour = "navy"

        # Initial colour/reference/vertical index from list of colours is set to 0
        self.index = 0
        self.ref_index = 0
        self.vertical_index = 0

        # List to hold the multiple 2D spectra in multiplot mode
        self.twoD_spectra = []

        self.linewidth = 1.0
        self.linewidth1D = 1.5

        self.x_difference = 0
        self.y_difference = 0

        # Initially set the transpose flag to False
        self.transpose = False
        self.transposed2D = False

        # Default options for pivot point for P1 phasing
        self.pivot_x_default = 0
        self.pivot_x = self.pivot_x_default

        self.pivot_y_default = 0
        self.pivot_y = self.pivot_y_default

        self.slice_mode = None

        self.do_not_update = False

        self.show_bottom_sizer = True

        # Suppress complex warning from numpy
        import warnings

        # warnings.simplefilter("ignore", np.ComplexWarning)  # For old numpy versions
        warnings.simplefilter(
            "ignore", np.exceptions.ComplexWarning
        )  # For new numpy versions

    def UpdateFrame(self):
        if self.do_not_update == False:
            # If the do_not_update flag is not True, update the frame
            self.canvas.draw()
            self.canvas.Refresh()
            self.canvas.Update()
            self.panel.Refresh()
            self.panel.Update()

    def create_button_panel_2D(self):

        # Create a sizer to choose a plot when in multiplot mode
        self.select_plot_label = wx.StaticBox(self, -1, "Select Plot:")
        self.select_plot_sizer = wx.StaticBoxSizer(self.select_plot_label, wx.VERTICAL)
        self.select_plot_sizer.AddSpacer(5)
        # Create a checkbox to select all plots
        self.select_all_checkbox = wx.CheckBox(self, label="Select All")
        self.select_all_checkbox.SetValue(False)

        self.plot_combobox = wx.ComboBox(
            self, choices=["Main Plot"], style=wx.CB_READONLY
        )
        self.plot_combobox.Bind(wx.EVT_COMBOBOX, self.OnSelectPlot2D)
        self.select_plot_sizer.Add(self.plot_combobox, 1, wx.ALIGN_CENTER_HORIZONTAL, 5)
        self.select_plot_sizer.AddSpacer(5)
        self.select_plot_sizer.Add(
            self.select_all_checkbox, 1, wx.ALIGN_CENTER_HORIZONTAL, 5
        )
        self.select_plot_sizer.AddSpacer(5)

        # Create a button to change the labels of the x and y axes (Don't include this button if in 3D mode)
        width = 100
        if self.threeDprojection == False:
            self.label_button = wx.Button(
                self, label="Change Labels", size=(width, height)
            )
            self.label_button.Bind(wx.EVT_BUTTON, self.OnLabelButton)

            self.save_session_button = wx.Button(
                self, label="Save Session", size=(width, height)
            )
            self.save_session_button.Bind(wx.EVT_BUTTON, self.OnSaveSessionButton2D)

        self.reset_button = wx.Button(self, label="Reset", size=(width, height))
        self.reset_button.Bind(wx.EVT_BUTTON, self.OnResetButton2D)

        # Create a button to transpose the given NMR spectrum
        self.transpose_button = wx.Button(self, label="Transpose", size=(width, height))
        self.transpose_button.Bind(wx.EVT_BUTTON, self.OnTransposeButton)

        # Create a button to stack the slices of the given NMR spectrum
        self.stack_button = wx.Button(self, label="Stack Slices", size=(width, height))
        self.stack_button.Bind(wx.EVT_BUTTON, self.OnStackButton)

        # Create a button for Re-Processing
        self.reprocess_button = wx.Button(
            self, label="Re-Process", size=(width, height)
        )
        self.reprocess_button.Bind(wx.EVT_BUTTON, self.OnReprocessButton)

        # Create a button to fit the diffusion data of the given NMR spectrum
        self.fit_diffusion_button = wx.Button(
            self, label="Fit Diffusion", size=(width, height)
        )
        self.fit_diffusion_button.Bind(wx.EVT_BUTTON, self.OnFitDiffusionButton)

        # Create a button to fit the relaxation data of the given NMR spectrum
        self.fit_relax_button = wx.Button(
            self, label="Fit Relaxation", size=(width, height)
        )
        self.fit_relax_button.Bind(wx.EVT_BUTTON, self.OnFitRelaxButton)

        # Create a button which will open a CESTView panel to analyse pseudo2D CEST data
        self.CEST_button = wx.Button(self, label="CEST Analysis", size=(width, height))
        self.CEST_button.Bind(wx.EVT_BUTTON, self.OnCESTButton)

        # Create a button which will make the correct files in order to perform uSTA analysis
        self.uSTA_button = wx.Button(self, label="uSTA", size=(width, height))
        self.uSTA_button.Bind(wx.EVT_BUTTON, self.OnuSTAButton)

        # Create a button to toggle the main sizer between shown and hidden
        self.toggle_button = wx.Button(self, label="Hide Options", size=(width, height))
        self.toggle_button.Bind(wx.EVT_BUTTON, self.OnHideButton)

        # Add the buttons to a sizer
        self.general_options_sizer = wx.BoxSizer(wx.HORIZONTAL)
        if self.threeDprojection == False:
            self.general_options_sizer.Add(self.label_button)
            self.general_options_sizer.AddSpacer(5)
            self.general_options_sizer.Add(self.reset_button)
            self.general_options_sizer.AddSpacer(5)
            self.general_options_sizer.Add(self.save_session_button)
            self.general_options_sizer.AddSpacer(5)
            self.general_options_sizer.Add(self.reprocess_button)

            # Add the diffusion and relaxation fit sizers to their own sizer

            self.fit_sizer = wx.BoxSizer(wx.HORIZONTAL)
            self.fit_sizer.Add(self.transpose_button)
            self.fit_sizer.AddSpacer(5)
            self.fit_sizer.Add(self.stack_button)
            self.fit_sizer.AddSpacer(5)
            self.fit_sizer.Add(self.fit_diffusion_button)
            self.fit_sizer.AddSpacer(5)
            self.fit_sizer.Add(self.fit_relax_button)

            self.hide_sizer = wx.BoxSizer(wx.HORIZONTAL)
            self.hide_sizer.Add(self.toggle_button)
            self.hide_sizer.AddSpacer(5)
            self.hide_sizer.Add(self.CEST_button)
            self.hide_sizer.AddSpacer(5)
            self.hide_sizer.Add(self.uSTA_button)

        else:
            self.general_options_sizer.Add(
                self.transpose_button, 1, wx.EXPAND | wx.ALL, 5
            )
            self.general_options_sizer.AddSpacer(15)
            self.general_options_sizer.Add(self.stack_button, 1, wx.EXPAND | wx.ALL, 5)
            self.general_options_sizer.AddSpacer(15)
            self.general_options_sizer.Add(
                self.fit_diffusion_button, 1, wx.EXPAND | wx.ALL, 5
            )
            self.general_options_sizer.AddSpacer(15)
            self.general_options_sizer.Add(
                self.fit_relax_button, 1, wx.EXPAND | wx.ALL, 5
            )

            self.hide_sizer = wx.BoxSizer(wx.HORIZONTAL)
            self.hide_sizer.Add(self.toggle_button)
            self.hide_sizer.AddSpacer(5)
            self.hide_sizer.Add(self.reset_button)

        # Create a sizer to phase the data
        self.phasing_label = wx.StaticBox(self, -1, "Phasing:")
        self.phasing_sizer = wx.StaticBoxSizer(self.phasing_label, wx.HORIZONTAL)
        self.P0_label = wx.StaticText(self, label="P0 (Coarse):")
        self.P1_label = wx.StaticText(self, label="P1 (Coarse):")
        self.P0_slider = FloatSlider(
            self,
            id=-1,
            value=0,
            minval=-180,
            maxval=180,
            res=0.1,
            size=(int(self.width / 6.5), height),
        )
        self.P1_slider = FloatSlider(
            self,
            id=-1,
            value=0,
            minval=-180,
            maxval=180,
            res=0.1,
            size=(int(self.width / 6.5), height),
        )
        self.P0_slider.Bind(wx.EVT_SLIDER, self.OnSliderScroll2D)
        self.P1_slider.Bind(wx.EVT_SLIDER, self.OnSliderScroll2D)
        self.P0_label_fine = wx.StaticText(self, label="P0 (Fine):")
        self.P1_label_fine = wx.StaticText(self, label="P1 (Fine):")
        self.P0_slider_fine = FloatSlider(
            self,
            id=-1,
            value=0,
            minval=-10,
            maxval=10,
            res=0.01,
            size=(int(self.width / 6.5), height),
        )
        self.P1_slider_fine = FloatSlider(
            self,
            id=-1,
            value=0,
            minval=-10,
            maxval=10,
            res=0.01,
            size=(int(self.width / 6.5), height),
        )
        self.P0_slider_fine.Bind(wx.EVT_SLIDER, self.OnSliderScroll2D)
        self.P1_slider_fine.Bind(wx.EVT_SLIDER, self.OnSliderScroll2D)
        self.P0_total = wx.StaticText(self, label="P0 (Total):")
        self.P1_total = wx.StaticText(self, label="P1 (Total):")
        self.P0_total_value = wx.StaticText(self, label="0")
        self.P1_total_value = wx.StaticText(self, label="0")

        self.P0_label_sizer = wx.BoxSizer(wx.VERTICAL)
        self.P0_label_sizer.Add(self.P0_label)
        self.P0_label_sizer.AddSpacer(10)
        self.P0_label_sizer.Add(self.P0_label_fine)
        self.P0_label_sizer.AddSpacer(10)
        self.P0_label_sizer.Add(self.P0_total)

        self.P0_slider_sizer = wx.BoxSizer(wx.VERTICAL)
        self.P0_slider_sizer.Add(self.P0_slider, wx.ALIGN_CENTER_HORIZONTAL, 0)
        self.P0_slider_sizer.AddSpacer(10)
        self.P0_slider_sizer.Add(self.P0_slider_fine, wx.ALIGN_CENTER_HORIZONTAL, 0)
        self.P0_slider_sizer.AddSpacer(10)
        self.P0_slider_sizer.Add(self.P0_total_value, wx.ALIGN_CENTER_HORIZONTAL, 5)

        self.P1_label_sizer = wx.BoxSizer(wx.VERTICAL)
        self.P1_label_sizer.Add(self.P1_label)
        self.P1_label_sizer.AddSpacer(10)
        self.P1_label_sizer.Add(self.P1_label_fine)
        self.P1_label_sizer.AddSpacer(10)
        self.P1_label_sizer.Add(self.P1_total)

        self.P1_slider_sizer = wx.BoxSizer(wx.VERTICAL)
        self.P1_slider_sizer.Add(self.P1_slider, wx.ALIGN_CENTER_HORIZONTAL, 0)
        self.P1_slider_sizer.AddSpacer(10)
        self.P1_slider_sizer.Add(self.P1_slider_fine, wx.ALIGN_CENTER_HORIZONTAL, 0)
        self.P1_slider_sizer.AddSpacer(10)
        self.P1_slider_sizer.Add(self.P1_total_value, wx.ALIGN_CENTER_HORIZONTAL, 5)
        self.P1_slider_sizer.AddSpacer(10)

        # Add a button to set the pivot point for phasing
        self.pivot_button = wx.Button(self, label="Set Pivot Point")
        self.pivot_button.Bind(wx.EVT_BUTTON, self.OnPivotButton2D)
        self.pivot_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.pivot_sizer.Add(self.pivot_button)

        # Add a button to remove the pivot point
        self.remove_pivot_button = wx.Button(self, label="Remove Pivot Point")
        self.remove_pivot_button.Bind(wx.EVT_BUTTON, self.OnRemovePivotButton2D)
        self.pivot_sizer.AddSpacer(20)
        self.pivot_sizer.Add(self.remove_pivot_button)

        self.P1_slider_sizer.AddSpacer(5)
        self.P1_slider_sizer.Add(self.pivot_sizer, wx.ALIGN_CENTER_HORIZONTAL, 1)

        self.phasing_sizer.Add(self.P0_label_sizer, wx.ALIGN_TOP)
        self.phasing_sizer.AddSpacer(10)
        self.phasing_sizer.Add(self.P0_slider_sizer, wx.ALIGN_TOP)
        self.phasing_sizer.AddSpacer(10)
        self.phasing_sizer.Add(self.P1_label_sizer, wx.ALIGN_TOP)
        self.phasing_sizer.AddSpacer(10)
        self.phasing_sizer.Add(self.P1_slider_sizer, wx.ALIGN_TOP)

        # Create a sizer for changing the contour levels of the spectrum
        self.contour_label = wx.StaticBox(self, -1, "Contour Start = max(data)/x")
        self.contour_sizer = wx.StaticBoxSizer(self.contour_label, wx.VERTICAL)
        self.csizer = wx.BoxSizer(wx.HORIZONTAL)
        self.contour2_label = wx.StaticText(self, label="x:")
        self.contour_slider = FloatSlider(
            self, id=-1, value=1, minval=0, maxval=3, res=0.01, size=(200, height)
        )
        self.contour_slider.Bind(wx.EVT_SLIDER, self.OnMinContour2D)
        self.csizer.Add(self.contour2_label)
        self.csizer.AddSpacer(5)
        self.csizer.Add(self.contour_slider)
        self.contour_sizer.AddSpacer(5)
        self.contour_sizer.Add(self.csizer)
        self.contour_sizer.AddSpacer(5)
        self.contour_value_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.contour_value_sizer.AddSpacer(100)
        self.contour_value_label = wx.StaticText(self, label="10")
        self.contour_value_sizer.Add(self.contour_value_label)
        self.contour_sizer.Add(self.contour_value_sizer)

        # Create a sizer for changing the y axis limits for a 1D slice
        self.intensity_label = wx.StaticBox(self, -1, "1D Y Axis Zoom (%):")
        self.intensity_sizer = wx.StaticBoxSizer(self.intensity_label, wx.VERTICAL)
        self.intensity_slider = FloatSlider(
            self, id=-1, value=0, minval=-1, maxval=10, res=0.01, size=(250, height)
        )
        self.intensity_slider.Bind(wx.EVT_SLIDER, self.OnIntensityScroll2D)
        self.intensity_sizer.AddSpacer(5)
        self.intensity_sizer.Add(self.intensity_slider)
        self.intensity_sizer.AddSpacer(5)

        # Create a sizer for multiplying the 2D data by a constant, this is useful when overlaying different datasets with different intensities
        self.multiply_label = wx.StaticBox(
            self, -1, "Multiply 2D Data by " + "n\u207f" + ":"
        )
        self.multiply_sizer = wx.StaticBoxSizer(self.multiply_label, wx.VERTICAL)
        self.multiply_inner_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.multiply_ranges = [
            "1.01",
            "1.1",
            "1.5",
            "2",
            "5",
            "10",
            "50",
            "100",
            "1000",
            "10000",
            "100000",
            "1000000",
            "10000000",
            "100000000",
            "1000000000",
        ]
        self.multiply_slider = FloatSlider(
            self,
            id=-1,
            value=1.0,
            minval=0,
            maxval=float(self.multiply_ranges[0]),
            res=0.01,
            size=(230, height),
        )
        self.multiply_slider.Bind(wx.EVT_SLIDER, self.OnMultiplyScroll2D)
        self.multiply_inner_sizer.AddSpacer(5)
        self.multiply_inner_sizer.Add(self.multiply_slider)
        self.multiply_sizer.AddSpacer(5)
        self.multiply_sizer.Add(self.multiply_inner_sizer)
        self.multiply_sizer.AddSpacer(5)
        self.multiply_value_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.multiply_value_n_label = wx.StaticText(self, label="n: ")
        self.multiply_value_label = wx.StaticText(self, label="1.00")
        self.multiply_value_sizer.Add(self.multiply_value_n_label)
        self.multiply_value_sizer.AddSpacer(5)
        self.multiply_value_sizer.Add(self.multiply_value_label)

        self.multiply_value_range_label = wx.StaticText(self, label="Range:")
        self.multiply_value_sizer.AddSpacer(30)
        self.multiply_value_sizer.Add(self.multiply_value_range_label)

        # Make a combobox to select the multiply range
        self.multiply_range_chooser2d = wx.ComboBox(
            self, value=self.multiply_ranges[0], choices=self.multiply_ranges
        )
        self.multiply_range_chooser2d.Bind(wx.EVT_COMBOBOX, self.OnMultiplyCombo2D)
        self.multiply_value_sizer.AddSpacer(5)
        self.multiply_value_sizer.Add(self.multiply_range_chooser2d)

        self.multiply_sizer.Add(
            self.multiply_value_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL
        )

        # Add a slider to change the number of contour levels
        self.contour_levels_label = wx.StaticBox(self, -1, "Contour Levels:")
        self.contour_levels = wx.StaticBoxSizer(self.contour_levels_label, wx.VERTICAL)
        self.contour_levels.AddSpacer(5)
        self.contour_levels_slider = FloatSlider(
            self, id=-1, value=20, minval=1, maxval=30, res=1, size=(215, height)
        )
        self.contour_levels_slider.Bind(wx.EVT_SLIDER, self.OnContourLevels)
        self.contour_levels.Add(self.contour_levels_slider)
        self.contour_levels.AddSpacer(5)

        # Add sliders to move the 2D plots left/right/up/down with a combobox to choose the scale of the slider
        self.move_label = wx.StaticBox(self, -1, "Move 2D Plot:")
        self.move_sizer_total = wx.StaticBoxSizer(self.move_label, wx.HORIZONTAL)
        self.move_sizer = wx.BoxSizer(wx.VERTICAL)
        self.move_x = wx.BoxSizer(wx.HORIZONTAL)
        self.move_y = wx.BoxSizer(wx.HORIZONTAL)
        self.move_ranges = wx.BoxSizer(wx.VERTICAL)
        self.move_x.Add(wx.StaticText(self, label="X:"))
        self.move_y.Add(wx.StaticText(self, label="Y:"))
        self.move_x.AddSpacer(5)
        self.move_y.AddSpacer(5)

        self.move_x_slider = FloatSlider(
            self,
            id=-1,
            value=0,
            minval=-self.reference_rangeX,
            maxval=self.reference_rangeX,
            res=self.reference_rangeX / 1000,
            size=(int(self.width / 3.5), height),
        )
        self.move_y_slider = FloatSlider(
            self,
            id=-1,
            value=0,
            minval=-self.reference_rangeY,
            maxval=self.reference_rangeY,
            res=self.reference_rangeY / 1000,
            size=(int(self.width / 3.5), height),
        )
        self.move_x_slider.Bind(wx.EVT_SLIDER, self.OnMoveX)
        self.move_y_slider.Bind(wx.EVT_SLIDER, self.OnMoveY)
        self.reference_range_chooserX = wx.ComboBox(
            self,
            value=self.reference_range_values[0],
            choices=self.reference_range_values,
        )
        self.reference_range_chooserX.Bind(wx.EVT_COMBOBOX, self.OnReferenceComboX)
        self.reference_range_chooserY = wx.ComboBox(
            self,
            value=self.reference_range_values[0],
            choices=self.reference_range_values,
        )
        self.reference_range_chooserY.Bind(wx.EVT_COMBOBOX, self.OnReferenceComboY)

        self.move_x.Add(self.move_x_slider)
        self.move_x.AddSpacer(5)
        self.move_ranges.Add(
            self.reference_range_chooserX, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.move_y.Add(self.move_y_slider)
        self.move_y.AddSpacer(5)
        self.move_ranges.AddSpacer(10)
        self.move_ranges.Add(
            self.reference_range_chooserY, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.move_ranges.AddSpacer(5)
        self.move_ranges.Add(
            wx.StaticText(self, label="Range (ppm)"), 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.move_sizer.Add(self.move_x)
        self.move_sizer.AddSpacer(10)
        self.move_sizer.Add(self.move_y)
        self.move_sizer.AddSpacer(5)
        self.move_values = wx.BoxSizer(wx.HORIZONTAL)

        self.move_values.Add(wx.StaticText(self, label="X Movement (ppm):"))
        self.move_values.AddSpacer(10)
        self.move_x_value_label = wx.StaticText(self, label="0.0")
        self.move_values.Add(self.move_x_value_label)
        self.move_values.AddSpacer(50)
        self.move_values.Add(wx.StaticText(self, label="Y Movement (ppm):"))
        self.move_values.AddSpacer(10)
        self.move_y_value_label = wx.StaticText(self, label="0.0")
        self.move_values.Add(self.move_y_value_label)
        self.move_sizer.Add(self.move_values, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.move_sizer_total.Add(self.move_sizer, 0, wx.ALIGN_CENTER_VERTICAL)
        self.move_sizer_total.AddSpacer(5)
        self.move_sizer_total.Add(self.move_ranges, 0, wx.ALIGN_CENTER_VERTICAL)
        leftover_space = (
            self.phasing_sizer.GetSize()[0]
            - self.select_plot_sizer.GetSize()[0]
            - self.move_sizer_total.GetSize()[0]
        )
        self.move_sizer_total.AddSpacer(leftover_space)

        # Create a slider to adjust contour linewiths
        self.contour_width_label = wx.StaticBox(self, -1, "Contour Linewidth:")
        self.contour_width = wx.StaticBoxSizer(self.contour_width_label, wx.VERTICAL)
        self.contour_width.AddSpacer(5)
        self.contour_width_slider = FloatSlider(
            self, id=-1, value=1, minval=0.1, maxval=2, res=0.1, size=(215, height)
        )
        self.contour_width_slider.Bind(wx.EVT_SLIDER, self.OnContourWidth)
        self.contour_width.Add(self.contour_width_slider)

        # Create a slider to adjust 1D slice linewidths
        self.linewidth_label = wx.StaticBox(self, -1, "1D Slice Linewidth:")
        self.line_width = wx.StaticBoxSizer(self.linewidth_label, wx.VERTICAL)
        self.line_width.AddSpacer(5)
        self.line_width_slider = FloatSlider(
            self, id=-1, value=1, minval=0.1, maxval=2, res=0.1, size=(250, height)
        )
        self.line_width_slider.Bind(wx.EVT_SLIDER, self.On2DLinewidth)
        self.line_width.Add(self.line_width_slider)

        # Put all the sizers together

        # Sizer to hold all options/sliders/buttons for 2D contour scaling, movement left/right/up/down etc
        self.twoD_sizer = wx.BoxSizer(wx.VERTICAL)
        self.twoD_sizer.Add(self.contour_sizer)
        self.twoD_sizer.AddSpacer(10)
        self.twoD_sizer.Add(self.contour_width)
        self.twoD_sizer.AddSpacer(10)
        self.twoD_sizer.Add(self.contour_levels)
        # Sizer to hold options/sliders/buttons for 1D line plot slices
        self.oneD_line_sizer = wx.BoxSizer(wx.VERTICAL)
        self.oneD_line_sizer.Add(self.intensity_sizer)
        self.oneD_line_sizer.AddSpacer(10)
        self.oneD_line_sizer.Add(self.line_width)
        self.oneD_line_sizer.AddSpacer(10)
        self.oneD_line_sizer.Add(self.multiply_sizer)

        # Sizer to hold all sliders/buttons for spectrum selection and moving spectra left/right/up/down etc
        self.top_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.top_left_sizer = wx.BoxSizer(wx.VERTICAL)
        self.top_left_sizer.Add(self.select_plot_sizer)
        self.top_left_sizer.AddSpacer(10)
        self.top_sizer.Add(self.top_left_sizer)
        self.top_sizer.AddSpacer(20)
        self.top_sizer.Add(self.move_sizer_total)

        # Sizer to hold all the phasing options
        self.bottom_left_sizer = wx.BoxSizer(wx.VERTICAL)
        self.bottom_left_sizer.Add(self.top_sizer)
        self.bottom_left_sizer.AddSpacer(10)
        self.bottom_left_sizer.Add(self.phasing_sizer)
        self.bottom_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.bottom_sizer.Add(self.bottom_left_sizer)
        self.bottom_sizer.AddSpacer(20)

        self.bottom_right_sizer = wx.BoxSizer(wx.VERTICAL)
        self.right_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.right_sizer.Add(self.twoD_sizer, 0, wx.EXPAND)
        self.right_sizer.AddSpacer(10)
        self.right_sizer.Add(self.oneD_line_sizer, 0, wx.EXPAND)
        self.bottom_right_sizer.Add(self.right_sizer)
        self.bottom_right_sizer.AddSpacer(10)
        self.buttons_sizer = wx.BoxSizer(wx.VERTICAL)
        self.buttons_sizer.Add(
            self.general_options_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.buttons_sizer.AddSpacer(10)
        if self.threeDprojection == False:
            self.buttons_sizer.Add(self.fit_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.buttons_sizer.AddSpacer(10)
        self.buttons_sizer.Add(self.hide_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.bottom_right_sizer.Add(self.buttons_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)

        self.bottom_right_sizer.AddSpacer(20)

        self.bottom_sizer.Add(self.bottom_right_sizer)

        self.slice_mode = None

    def OnHideButton(self, event):
        if self.show_bottom_sizer == True:
            self.main_sizer.Hide(self.bottom_sizer)
            self.main_sizer.Show(self.show_button_sizer)
            self.UpdateFrame()
            self.Layout()
            self.show_bottom_sizer = False
        else:
            self.main_sizer.Show(self.bottom_sizer)
            self.main_sizer.Hide(self.show_button_sizer)
            self.show_bottom_sizer = True
            self.UpdateFrame()
            self.Layout()

    def OnReprocessButton(self, event):
        # Open an instance of SpinProcess
        if self.parent.path != "":
            os.chdir(self.parent.path)
        from SpinExplorer.SpinProcess import SpinProcess

        reprocessing_frame = SpinProcess(self)
        reprocessing_frame.reprocess = True
        if self.parent.cwd != "":
            os.chdir(self.parent.cwd)

    def OnSaveSessionButton2D(self, event):
        # Function to save the current session
        # Give a file menu popout to ask the user which directory to save the session in
        dlg = wx.FileDialog(
            self,
            "Save Session",
            wildcard="Session files (*.session)|*.session",
            style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT,
        )
        dlg.SetDirectory(os.getcwd())
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.save_session2D(path)
            dlg.Destroy()
        else:
            return

    def save_session2D(self, path):
        # Function to save the current session
        # Save the current session to a file
        with open(path, "w") as f:
            f.write("2D\n")
            if self.multiplot_mode == False:
                f.write("MultiplotMode:False\n")
                f.write("Transposed2D:" + str(self.transposed2D) + "\n")
                if platform == "windows":
                    if (
                        os.path.exists(
                            str(self.parent.nmrdata.path)
                            + "\\"
                            + str(self.parent.nmrdata.file)
                        )
                        == True
                    ):
                        f.write(
                            "file_path:"
                            + str(self.parent.nmrdata.path)
                            + "\\"
                            + str(self.parent.nmrdata.file)
                            + "\n"
                        )
                    else:
                        f.write("file_path:" + str(self.parent.nmrdata.file) + "\n")
                else:
                    if (
                        os.path.exists(
                            str(self.parent.nmrdata.path)
                            + "/"
                            + str(self.parent.nmrdata.file)
                        )
                        == True
                    ):
                        f.write(
                            "file_path:"
                            + str(self.parent.nmrdata.path)
                            + "/"
                            + str(self.parent.nmrdata.file)
                            + "\n"
                        )
                    else:
                        f.write("file_path:" + str(self.parent.nmrdata.file) + "\n")

                f.write("p0 Coarse:" + str(self.P0_slider.GetValue()) + "\n")
                f.write("p1 Coarse:" + str(self.P1_slider.GetValue()) + "\n")
                f.write("p0 Fine:" + str(self.P0_slider_fine.GetValue()) + "\n")
                f.write("p1 Fine:" + str(self.P1_slider_fine.GetValue()) + "\n")
                f.write("move x:" + str(self.move_x_slider.GetValue()) + "\n")
                f.write("move y:" + str(self.move_y_slider.GetValue()) + "\n")
                f.write(
                    "move x range index:"
                    + str(self.reference_range_chooserX.GetSelection())
                    + "\n"
                )
                f.write(
                    "move y range index:"
                    + str(self.reference_range_chooserY.GetSelection())
                    + "\n"
                )
                f.write(
                    "contour linewidth:"
                    + str(self.contour_width_slider.GetValue())
                    + "\n"
                )
                f.write(
                    "multiply factor:" + str(self.multiply_slider.GetValue()) + "\n"
                )
                f.write(
                    "contour levels:"
                    + str(self.contour_levels_slider.GetValue())
                    + "\n"
                )
                f.write("transposed:False\n")
            else:
                f.write("MultiplotMode:True\n")
                f.write("Transposed2D:" + str(self.transposed2D) + "\n")
                for i in range(len(self.values_dictionary)):
                    f.write("file_path:" + self.values_dictionary[i]["path"] + "\n")
                    f.write("title:" + self.values_dictionary[i]["title"] + "\n")
                    f.write(
                        "p0 Coarse:"
                        + str(self.values_dictionary[i]["p0 Coarse"])
                        + "\n"
                    )
                    f.write(
                        "p1 Coarse:"
                        + str(self.values_dictionary[i]["p1 Coarse"])
                        + "\n"
                    )
                    f.write(
                        "p0 Fine:" + str(self.values_dictionary[i]["p0 Fine"]) + "\n"
                    )
                    f.write(
                        "p1 Fine:" + str(self.values_dictionary[i]["p1 Fine"]) + "\n"
                    )
                    f.write("move x:" + str(self.values_dictionary[i]["move x"]) + "\n")
                    f.write("move y:" + str(self.values_dictionary[i]["move y"]) + "\n")
                    f.write(
                        "move x range index:"
                        + str(self.values_dictionary[i]["move x range index"])
                        + "\n"
                    )
                    f.write(
                        "move y range index:"
                        + str(self.values_dictionary[i]["move y range index"])
                        + "\n"
                    )
                    f.write(
                        "contour linewidth:"
                        + str(self.values_dictionary[i]["contour linewidth"])
                        + "\n"
                    )
                    f.write(
                        "multiply factor:"
                        + str(self.values_dictionary[i]["multiply factor"])
                        + "\n"
                    )
                    f.write(
                        "contour levels:"
                        + str(self.values_dictionary[i]["contour levels"])
                        + "\n"
                    )
                    try:
                        f.write(
                            "transposed:"
                            + str(self.values_dictionary[i]["transposed"])
                            + "\n"
                        )
                    except:
                        f.write("transposed:False\n")
            f.close()

    def OnResetButton2D(self, event):
        if self.multiplot_mode == False:
            # Get the user to confirm if they want to reset plot
            dlg = wx.MessageDialog(
                self,
                "This will reset all the parameters to their default values. Do you want to continue?",
                "Reset Plot",
                wx.YES_NO | wx.ICON_QUESTION,
            )
            result = dlg.ShowModal()
            if result == wx.ID_YES:
                # Reset the plot
                self.P0_slider.SetValue(0)
                self.P1_slider.SetValue(0)
                self.P0_slider_fine.SetValue(0)
                self.P1_slider_fine.SetValue(0)
                self.P0_total_value.SetLabel("0.00")
                self.P1_total_value.SetLabel("0.00")
                self.contour_width_slider.SetValue(1)
                self.contour_slider.SetValue(1)
                self.contour_value_label.SetLabel("10")
                self.contour_levels_slider.SetValue(20)
                self.move_x_slider.SetValue(0)
                self.move_y_slider.SetValue(0)
                self.move_x_value_label.SetLabel("0.00")
                self.move_y_value_label.SetLabel("0.00")
                self.multiply_slider.SetValue(0)
                self.multiply_value_label.SetLabel("0")
                self.line_width_slider.SetValue(1)
                # if(self.transposed2D==True):
                #     self.OnTransposeButton(event)
                # self.transposed2D = False
                self.OnMoveX(event)
                self.OnMoveY(event)
                self.OnMultiplyScroll2D(event)
                self.OnMinContour2D(event)
                self.OnSliderScroll2D(event)
                self.OnIntensityScroll2D(event)
                self.OnContourWidth(event)
                self.OnContourLevels(event)
                self.On2DLinewidth(event)
                self.UpdateFrame()
        else:
            if self.select_all_checkbox.IsChecked() == False:
                # Get the user to confirm if they want to reset plot
                dlg = wx.MessageDialog(
                    self,
                    "This will reset all the parameters to their default values for the selected plot. Do you want to continue?",
                    "Reset Plot",
                    wx.YES_NO | wx.ICON_QUESTION,
                )
                result = dlg.ShowModal()
                if result == wx.ID_YES:
                    # Reset the plot
                    self.P0_slider.SetValue(0)
                    self.P1_slider.SetValue(0)
                    self.P0_slider_fine.SetValue(0)
                    self.P1_slider_fine.SetValue(0)
                    self.P0_total_value.SetLabel("0.00")
                    self.P1_total_value.SetLabel("0.00")
                    self.contour_width_slider.SetValue(1)
                    self.contour_slider.SetValue(1)
                    self.contour_value_label.SetLabel("10")
                    self.contour_levels_slider.SetValue(20)
                    self.move_x_slider.SetValue(0)
                    self.move_y_slider.SetValue(0)
                    self.move_x_value_label.SetLabel("0.00")
                    self.move_y_value_label.SetLabel("0.00")
                    self.multiply_slider.SetValue(0)
                    self.multiply_value_label.SetLabel("0")
                    self.line_width_slider.SetValue(1)
                    self.values_dictionary[self.active_plot_index]["p0 Coarse"] = 0
                    self.values_dictionary[self.active_plot_index]["p1 Coarse"] = 0
                    self.values_dictionary[self.active_plot_index]["p0 Fine"] = 0
                    self.values_dictionary[self.active_plot_index]["p1 Fine"] = 0
                    self.values_dictionary[self.active_plot_index][
                        "contour linewidth"
                    ] = 1
                    self.values_dictionary[self.active_plot_index][
                        "contour levels"
                    ] = 20
                    self.values_dictionary[self.active_plot_index]["move x"] = 0
                    self.values_dictionary[self.active_plot_index]["move y"] = 0
                    self.values_dictionary[self.active_plot_index][
                        "multiply factor"
                    ] = 0
                    self.values_dictionary[self.active_plot_index]["linewidth 1D"] = 1
                    self.values_dictionary[self.active_plot_index][
                        "move x range index"
                    ] = 0
                    self.values_dictionary[self.active_plot_index][
                        "move y range index"
                    ] = 0

                    # if(self.transposed2D==True):
                    #     self.OnTransposeButton(event)
                    # self.transposed2D = False
                    self.OnMoveX(event)
                    self.OnMoveY(event)
                    self.OnMinContour2D(event)
                    self.OnSliderScroll2D(event)
                    self.OnIntensityScroll2D(event)
                    self.OnContourWidth(event)
                    self.OnContourLevels(event)
                    self.OnMultiplyScroll2D(event)
                    self.On2DLinewidth(event)
                    titles = []
                    for i in range(len(self.values_dictionary.keys())):
                        titles.append(self.values_dictionary[i]["title"])

                    self.ax.legend(self.files.custom_lines, titles)
                    self.UpdateFrame()
            else:
                # Get the user to confirm if they want to reset plot
                dlg = wx.MessageDialog(
                    self,
                    "This will reset all the parameters to their default values for all plots. Do you want to continue?",
                    "Reset Plot",
                    wx.YES_NO | wx.ICON_QUESTION,
                )
                result = dlg.ShowModal()
                if result == wx.ID_YES:
                    # Reset the plot
                    self.P0_slider.SetValue(0)
                    self.P1_slider.SetValue(0)
                    self.P0_slider_fine.SetValue(0)
                    self.P1_slider_fine.SetValue(0)
                    self.P0_total_value.SetLabel("0.00")
                    self.P1_total_value.SetLabel("0.00")
                    self.contour_width_slider.SetValue(1)
                    self.contour_slider.SetValue(1)
                    self.contour_value_label.SetLabel("10")
                    self.contour_levels_slider.SetValue(20)
                    self.move_x_slider.SetValue(0)
                    self.move_y_slider.SetValue(0)
                    self.move_x_value_label.SetLabel("0.00")
                    self.move_y_value_label.SetLabel("0.00")
                    self.multiply_slider.SetValue(0)
                    self.multiply_value_label.SetLabel("0")
                    self.line_width_slider.SetValue(1)
                    for key in self.values_dictionary:
                        self.values_dictionary[key]["p0 Coarse"] = 0
                        self.values_dictionary[key]["p1 Coarse"] = 0
                        self.values_dictionary[key]["p0 Fine"] = 0
                        self.values_dictionary[key]["p1 Fine"] = 0
                        self.values_dictionary[key]["contour linewidth"] = 1
                        self.values_dictionary[key]["contour levels"] = 20
                        self.values_dictionary[key]["move x"] = 0
                        self.values_dictionary[key]["move y"] = 0
                        self.values_dictionary[key]["multiply factor"] = 0
                        self.values_dictionary[key]["linewidth 1D"] = 1
                        self.values_dictionary[key]["move x range index"] = 0
                        self.values_dictionary[key]["move y range index"] = 0

                    # if(self.transposed2D==True):
                    #     self.OnTransposeButton(event)
                    # self.transposed2D = False
                    self.OnMoveX(event)
                    self.OnMoveY(event)
                    self.OnMinContour2D(event)
                    self.OnSliderScroll2D(event)
                    self.OnIntensityScroll2D(event)
                    self.OnContourWidth(event)
                    self.OnContourLevels(event)
                    self.OnMultiplyScroll2D(event)
                    self.On2DLinewidth(event)
                    titles = []
                    for i in range(len(self.values_dictionary.keys())):
                        titles.append(self.values_dictionary[i]["title"])

                    self.ax.legend(self.files.custom_lines, titles)
                    self.UpdateFrame()

    def OnLabelButton(self, event):
        # Get the current labels of the x and y axes
        if self.transposed2D == False:
            x_label = self.nmrdata.axislabels[1]
            y_label = self.nmrdata.axislabels[0]
        else:
            x_label = self.nmrdata.axislabels[0]
            y_label = self.nmrdata.axislabels[1]

        # Get the ppm values for the x and y axes
        x_ppms = self.ppms_0
        y_ppms = self.ppms_1

        # Create a window to allow the user to see the current labels and ppm values and change the labels accordingly
        self.dlg = wx.Dialog(self, title="Change Labels")
        self.dlg.SetSize(500, 200)

        # Create a sizer to hold the labels and ppm values
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.AddSpacer(10)

        self.label_change_label = wx.StaticBox(self.dlg, -1, "Input desired labels:")
        self.total_label_change_sizer = wx.StaticBoxSizer(
            self.label_change_label, wx.VERTICAL
        )
        self.total_label_change_sizer.AddSpacer(10)

        # Create a sizer to hold the x axis labels and ppm values
        x_sizer = wx.BoxSizer(wx.HORIZONTAL)
        x_sizer.AddSpacer(10)
        x_sizer.Add(wx.StaticText(self.dlg, label="X Axis Label:"))
        x_sizer.AddSpacer(5)
        self.xlabel_box = wx.TextCtrl(self.dlg, value=x_label, size=(100, 20))
        x_sizer.Add(self.xlabel_box)
        x_sizer.AddSpacer(10)
        x_ppm_limits = "{:.2f}".format(min(x_ppms)) + "-{:.2f}".format(max(x_ppms))
        x_sizer.Add(wx.StaticText(self.dlg, label="X Axis Limits (ppm):"))
        x_sizer.AddSpacer(5)
        x_sizer.Add(wx.StaticText(self.dlg, label=x_ppm_limits))
        self.total_label_change_sizer.Add(x_sizer)
        self.total_label_change_sizer.AddSpacer(10)

        # Create a sizer to hold the y axis labels and ppm values
        y_sizer = wx.BoxSizer(wx.HORIZONTAL)
        y_sizer.AddSpacer(10)
        y_sizer.Add(wx.StaticText(self.dlg, label="Y Axis Label:"))
        y_sizer.AddSpacer(5)
        self.ylabel_box = wx.TextCtrl(self.dlg, value=y_label, size=(100, 20))
        y_sizer.Add(self.ylabel_box)
        y_sizer.AddSpacer(10)
        y_ppm_limits = "{:.2f}".format(min(y_ppms)) + "-{:.2f}".format(max(y_ppms))
        y_sizer.Add(wx.StaticText(self.dlg, label="Y Axis Limits (ppm):"))
        y_sizer.AddSpacer(5)
        y_sizer.Add(wx.StaticText(self.dlg, label=y_ppm_limits))
        self.total_label_change_sizer.Add(y_sizer)
        self.total_label_change_sizer.AddSpacer(10)

        # Add a save and close button to the sizer
        save_button = wx.Button(self.dlg, label="Save and Close")
        save_button.Bind(wx.EVT_BUTTON, self.OnSaveLabels)
        self.total_label_change_sizer.Add(save_button, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.total_label_change_sizer.AddSpacer(10)

        sizer.Add(self.total_label_change_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)

        # Show the sizer in the dialog
        self.dlg.SetSizer(sizer)

        # Show the dialog
        self.dlg.ShowModal()

    def OnSaveLabels(self, event):
        # Get the new labels for the x and y axes
        x_label = self.xlabel_box.GetValue()
        y_label = self.ylabel_box.GetValue()

        # Update the labels in the plot
        if self.transposed2D == False:
            self.nmrdata.axislabels[1] = x_label
            self.nmrdata.axislabels[0] = y_label
        else:
            self.nmrdata.axislabels[0] = x_label
            self.nmrdata.axislabels[1] = y_label

        self.ax.set_xlabel(self.nmrdata.axislabels[1])
        self.ax.set_ylabel(self.nmrdata.axislabels[0])

        # Save the labels to a labels.txt file
        if self.parent.path != "":
            os.chdir(self.parent.path)
        with open("labels.txt", "w") as f:
            # write labels as label1, label2
            f.write(self.nmrdata.axislabels[0] + "," + self.nmrdata.axislabels[1])
        if self.parent.cwd != "":
            os.chdir(self.parent.cwd)

        self.dlg.Destroy()

        self.UpdateFrame()

    def OnuSTAButton(self, event):
        # Producing an output to the user telling them this will produce a .data file required for uSTA analysis. Full uSTA implementation is in development.
        dlg = wx.MessageDialog(
            self,
            "This will produce a .data file required for uSTA analysis. Full uSTA implementation is in development. Do you want to continue?",
            "uSTA",
            wx.YES_NO | wx.ICON_WARNING,
        )
        self.Raise()
        self.SetFocus()
        result = dlg.ShowModal()
        if result == wx.ID_NO:
            dlg.Destroy()
            return
        else:
            dlg.Destroy()
            # Determining if data has only two spectra (on resonance/off resonance)
            if len(self.nmrdata.data) > 2:
                # Try transposing the data
                usta_data = self.nmrdata.data.T
                if len(usta_data) > 2:
                    dlg = wx.MessageDialog(
                        self,
                        "There are more than 2 spectra in the pseudo2D data. Expected 2 spectra, one for on resonance and one for off resonance.",
                        "Error",
                        wx.OK | wx.ICON_WARNING,
                    )
                    self.Raise()
                    self.SetFocus()
                    dlg.ShowModal()
                    dlg.Destroy()
                    return

            # Try to read the acqus file to get D20 and PL10 value

            self.mixing_time = 0
            self.power_level = 0
            try:
                dvals_next_line = False
                plvals_next_line = False
                with open("acqus", "r") as file:
                    file_lines = file.readlines()
                    for line in file_lines:
                        if dvals_next_line == True:
                            self.mixing_time = line.split()[20]
                            dvals_next_line = False
                        if plvals_next_line == True:
                            self.power_level = line.split()[10]
                            plvals_next_line = False
                        if "##$D=" in line:
                            dvals_next_line = True
                        if "##$PL=" in line:
                            plvals_next_line = True
            except:
                pass

            # Getting the user to input/confirm the uSTA spectral parameters

            uSTA_input = uSTA_Dialog(title="Input uSTA parameters", parent=self)

    def draw_figure_2D(self):
        self.ax = self.fig.add_subplot(111)
        self.axes1D = self.ax.twinx()
        self.axes1D_2 = self.ax.twiny()

        self.pivot_line = self.axes1D_2.axvline(
            self.pivot_x_default, color="black", linestyle="--"
        )
        self.pivot_line.set_visible(False)

        self.pivot_line_y = self.axes1D_2.axhline(
            self.pivot_y_default, color="black", linestyle="--"
        )
        self.pivot_line_y.set_visible(False)

        self.key_press_connect = self.fig.canvas.mpl_connect(
            "key_press_event", self.on_key_2d
        )
        self.click_press_connect = self.fig.canvas.mpl_connect(
            "button_press_event", self.on_click_2d
        )

        contour_start = np.max(self.nmrdata.data) / 10  # contour level start value
        self.contour_num = 20  # number of contour levels
        self.contour_factor = 1.20  # scaling factor between contour levels
        # calculate contour levels
        self.cl = contour_start * self.contour_factor ** np.arange(self.contour_num)
        self.cl_neg = -contour_start * self.contour_factor ** np.flip(
            np.arange(self.contour_num)
        )

        # Get ppm values for x and y axis
        if self.nmrdata.file != ".":
            self.uc0 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=0)
            self.uc1 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=1)
        else:
            udic = ng.bruker.guess_udic(self.nmrdata.dic, self.nmrdata.data)
            self.uc0 = ng.fileiobase.uc_from_udic(udic, dim=0)
            self.uc1 = ng.fileiobase.uc_from_udic(udic, dim=1)
        self.ppms_0 = self.uc0.ppm_scale()
        self.ppms_1 = self.uc1.ppm_scale()
        self.new_x_ppms = self.ppms_0
        self.new_y_ppms = self.ppms_1
        self.X, self.Y = np.meshgrid(self.ppms_1, self.ppms_0)
        self.contour1 = self.ax.contour(
            self.Y,
            self.X,
            self.nmrdata.data * self.multiply_factor,
            self.cl,
            colors=self.cmap,
            linewidths=self.linewidth,
        )
        self.contour1_neg = self.ax.contour(
            self.Y,
            self.X,
            self.nmrdata.data * self.multiply_factor,
            self.cl_neg,
            colors=self.cmap_neg,
            linewidths=self.linewidth,
        )
        self.ax.set_xlabel(self.nmrdata.axislabels[1])
        self.ax.set_ylabel(self.nmrdata.axislabels[0])
        self.ax.set_xlim(max(self.ppms_0), min(self.ppms_0))
        self.ax.set_ylim(max(self.ppms_1), min(self.ppms_1))
        (self.line1,) = self.axes1D.plot(
            self.ppms_0,
            self.nmrdata.data[:, 1] * self.multiply_factor,
            color=self.slice_colour,
            linewidth=self.linewidth1D,
        )
        self.axes1D.set_yticks([])
        self.line2 = self.ax.axhline(self.ppms_1[1], color="k")
        intensity_percent = 10 ** (float(self.intensity_slider.GetValue()))
        self.axes1D.set_ylim(
            -(np.max(self.nmrdata.data) / 8) / (intensity_percent / 100),
            np.max(self.nmrdata.data) / (intensity_percent / 100),
        )
        self.line1.set_visible(False)
        self.line2.set_visible(False)
        (self.line3,) = self.axes1D_2.plot(
            self.nmrdata.data[1, :] * self.multiply_factor,
            self.ppms_1,
            color=self.slice_colour,
            linewidth=self.linewidth1D,
        )
        self.line4 = self.ax.axvline(self.ppms_0[1], color="k")
        self.axes1D_2.set_xlim(
            -(np.max(self.nmrdata.data) / 8) / (intensity_percent / 100),
            np.max(self.nmrdata.data) / (intensity_percent / 100),
        )
        self.axes1D_2.set_xticks([])
        self.line3.set_visible(False)
        self.line4.set_visible(False)

        self.files = FileDrop(self.canvas, self.ax, self)
        self.canvas.SetDropTarget(self.files)

        self.UpdateFrame()

    def OnPivotButton2D(self, event):
        # If the user has not selected a horizontal or vertical slice, give a message box to tell them to do so
        if self.slice_mode == None:
            wx.MessageBox(
                "This mode requires that a user has already selected their desired horizontal or vertical slice for P1 phasing. Please select slice and repeat.",
                "Pivot Point",
                wx.OK | wx.ICON_INFORMATION,
            )
            return

        # Get the user to select a pivot point for phasing by clicking on the spectrum
        # Give a message box to tell the user to click on the spectrum where they want the pivot point to be
        wx.MessageBox(
            "Click on the spectrum to set the location of the pivot point for P1 phasing.",
            "Pivot Point",
            wx.OK | wx.ICON_INFORMATION,
        )

        # Deactivate the click_press and key_press events
        self.fig.canvas.mpl_disconnect(self.click_press_connect)
        self.fig.canvas.mpl_disconnect(self.key_press_connect)

        # Allow the key press to select the pivot point
        self.pivot_press = self.fig.canvas.mpl_connect(
            "button_press_event", self.OnPivotClick2D
        )

    def OnPivotClick2D(self, event):
        self.x1, self.y1 = self.ax.transData.inverted().transform((event.x, event.y))
        # Check to see if currently in x or y slice mode
        if self.slice_mode == "x":
            # Function to get the x value of the pivot point for phasing
            self.pivot_x = event.xdata
            self.pivot_line.set_xdata([self.pivot_x])

            # Find the index of the point closest to the pivot point
            self.pivot_index = np.abs(self.new_x_ppms - self.x1).argmin()
            self.pivot_x = self.pivot_index
            self.pivot_line.set_visible(True)
        elif self.slice_mode == "y":
            # Function to get the y value of the pivot point for phasing
            self.pivot_y = event.ydata
            self.pivot_line_y.set_ydata([self.pivot_y])

            # Find the index of the point closest to the pivot point
            self.pivot_index_y = np.abs(self.new_y_ppms - self.y1).argmin()
            self.pivot_y = self.pivot_index_y
            self.pivot_line_y.set_visible(True)

        self.fig.canvas.mpl_disconnect(self.pivot_press)
        # Reactivate the click_press and key_press events
        self.key_press_connect = self.fig.canvas.mpl_connect(
            "key_press_event", self.on_key_2d
        )
        self.click_press_connect = self.fig.canvas.mpl_connect(
            "button_press_event", self.on_click_2d
        )

        self.UpdateFrame()

    def OnRemovePivotButton2D(self, event):

        if self.slice_mode == "x":
            if self.pivot_line.get_visible() != True:
                # There is no pivot point to remove
                wx.MessageBox(
                    "There is no pivot point to remove.",
                    "Remove Pivot Point",
                    wx.OK | wx.ICON_INFORMATION,
                )
            else:
                # Function to remove the pivot point for phasing
                self.pivot_x = self.pivot_x_default
                self.pivot_line.set_visible(False)
                self.OnSliderScroll2D(wx.EVT_SCROLL)
                self.key_press_connect = self.fig.canvas.mpl_connect(
                    "key_press_event", self.on_key_2d
                )
                self.click_press_connect = self.fig.canvas.mpl_connect(
                    "button_press_event", self.on_click_2d
                )
        elif self.slice_mode == "y":
            if self.pivot_line_y.get_visible() != True:
                # There is no pivot point to remove
                wx.MessageBox(
                    "There is no pivot point to remove.",
                    "Remove Pivot Point",
                    wx.OK | wx.ICON_INFORMATION,
                )
            else:
                # Function to remove the pivot point for phasing
                self.pivot_y = self.pivot_y_default
                self.pivot_line_y.set_visible(False)
                self.OnSliderScroll2D(wx.EVT_SCROLL)
                self.key_press_connect = self.fig.canvas.mpl_connect(
                    "key_press_event", self.on_key_2d
                )
                self.click_press_connect = self.fig.canvas.mpl_connect(
                    "button_press_event", self.on_click_2d
                )
        else:
            if self.pivot_line.get_visible() == True:
                self.pivot_x = self.pivot_x_default
                self.pivot_line.set_visible(False)
            if self.pivot_line_y.get_visible() == True:
                self.pivot_y = self.pivot_y_default
                self.pivot_line_y.set_visible(False)
            self.UpdateFrame()

    def OnSelectPlot2D(self, event):
        # Save the updated values for the previous plot for colour, linewidth, referencing, vertical scroll, and phasing
        if self.multiplot_mode == True:
            if self.reference_range_chooserX.GetSelection() < 0:
                self.values_dictionary[self.active_plot_index]["move x range index"] = 0
            else:
                self.values_dictionary[self.active_plot_index][
                    "move x range index"
                ] = self.reference_range_chooserX.GetSelection()
            self.values_dictionary[self.active_plot_index][
                "move x"
            ] = self.move_x_slider.GetValue()
            if self.reference_range_chooserY.GetSelection() < 0:
                self.values_dictionary[self.active_plot_index]["move y range index"] = 0
            else:
                self.values_dictionary[self.active_plot_index][
                    "move y range index"
                ] = self.reference_range_chooserY.GetSelection()
            self.values_dictionary[self.active_plot_index][
                "move y"
            ] = self.move_y_slider.GetValue()
            self.values_dictionary[self.active_plot_index][
                "p0 Coarse"
            ] = self.P0_slider.GetValue()
            self.values_dictionary[self.active_plot_index][
                "p1 Coarse"
            ] = self.P1_slider.GetValue()
            self.values_dictionary[self.active_plot_index][
                "p0 Fine"
            ] = self.P0_slider_fine.GetValue()
            self.values_dictionary[self.active_plot_index][
                "p1 Fine"
            ] = self.P1_slider_fine.GetValue()
            self.values_dictionary[self.active_plot_index][
                "contour linewidth"
            ] = self.contour_width_slider.GetValue()
            self.values_dictionary[self.active_plot_index][
                "linewidth 1D"
            ] = self.line_width_slider.GetValue()

        # Function to change the active plot when a user selects a new plot from the combobox
        self.multiplot_mode = True
        self.active_plot_index = self.plot_combobox.GetSelection()

        # Update the values in the GUI to reflect the previously saved values for the active plot
        self.reference_range_chooserX.SetSelection(
            self.values_dictionary[self.active_plot_index]["move x range index"]
        )
        self.OnReferenceComboX(
            wx.EVT_SCROLL
        )  # Using scroll as a random event to trigger the function
        self.move_x_slider.SetValue(
            self.values_dictionary[self.active_plot_index]["move x"]
        )
        self.reference_range_chooserY.SetSelection(
            self.values_dictionary[self.active_plot_index]["move y range index"]
        )
        self.OnReferenceComboY(wx.EVT_SCROLL)
        self.move_y_slider.SetValue(
            self.values_dictionary[self.active_plot_index]["move y"]
        )
        self.P0_slider.SetValue(
            self.values_dictionary[self.active_plot_index]["p0 Coarse"]
        )
        self.P1_slider.SetValue(
            self.values_dictionary[self.active_plot_index]["p1 Coarse"]
        )
        self.P0_slider_fine.SetValue(
            self.values_dictionary[self.active_plot_index]["p0 Fine"]
        )
        self.P1_slider_fine.SetValue(
            self.values_dictionary[self.active_plot_index]["p1 Fine"]
        )
        self.contour_width_slider.SetValue(
            self.values_dictionary[self.active_plot_index]["contour linewidth"]
        )
        self.line_width_slider.SetValue(
            self.values_dictionary[self.active_plot_index]["linewidth 1D"]
        )
        self.multiply_slider.SetValue(
            np.log10(self.values_dictionary[self.active_plot_index]["multiply factor"])
        )

        # Update the plot to reflect the previously saved values for the active plot
        self.OnMoveX(wx.EVT_SCROLL)
        self.OnMoveY(wx.EVT_SCROLL)
        self.OnSliderScroll2D(wx.EVT_SCROLL)
        self.OnMinContour2D(wx.EVT_SCROLL)

    def OnTransposeButton(self, event):
        if self.transposed2D == False:
            self.transposed2D = True
        else:
            self.transposed2D = False
        if self.multiplot_mode == False:
            xlim_old, ylim_old = self.ax.get_xlim(), self.ax.get_ylim()
            self.X_old, self.Y_old = self.X, self.Y
            self.new_x_ppms_old = self.new_x_ppms
            self.new_y_ppms_old = self.new_y_ppms
            self.new_x_ppms = self.new_y_ppms_old
            self.new_y_ppms = self.new_x_ppms_old
            self.X, self.Y = np.meshgrid(self.new_y_ppms, self.new_x_ppms)
            self.nmr_data_old = self.nmrdata.data
            self.nmrdata.data = self.nmr_data_old.T
            self.ax.clear()
            self.contour1 = self.ax.contour(
                self.Y,
                self.X,
                self.nmrdata.data * self.multiply_factor,
                self.cl,
                colors=self.cmap,
                linewidths=self.linewidth,
            )
            self.contour1_neg = self.ax.contour(
                self.Y,
                self.X,
                self.nmrdata.data * self.multiply_factor,
                self.cl_neg,
                colors=self.cmap_neg,
                linewidths=self.linewidth,
            )
            self.ax.set_xlim([max(self.new_x_ppms), min(self.new_x_ppms)])
            self.ax.set_ylim([max(self.new_y_ppms), min(self.new_y_ppms)])
            self.axislabels_old = self.nmrdata.axislabels[0], self.nmrdata.axislabels[1]
            self.nmrdata.axislabels[1] = self.axislabels_old[0]
            self.nmrdata.axislabels[0] = self.axislabels_old[1]

            uc0, uc1 = self.uc0, self.uc1

            self.uc0 = uc1
            self.uc1 = uc0

            self.ax.set_xlabel(self.nmrdata.axislabels[1])
            self.ax.set_ylabel(self.nmrdata.axislabels[0])

            # Swap the move x and move y sliders and comboboxes
            # Get the x and y movement selections and slider values
            move_x_range_index = self.reference_range_chooserX.GetSelection()
            move_x_value = self.move_x_slider.GetValue()
            move_y_range_index = self.reference_range_chooserY.GetSelection()
            move_y_value = self.move_y_slider.GetValue()

            # Set the x and y movement selections and slider values
            self.reference_range_chooserX.SetSelection(move_y_range_index)
            self.OnReferenceComboX(wx.EVT_SCROLL)
            self.move_x_slider.SetValue(move_y_value)
            self.reference_range_chooserY.SetSelection(move_x_range_index)
            self.OnReferenceComboY(wx.EVT_SCROLL)
            self.move_y_slider.SetValue(move_x_value)

            self.OnMinContour2D(wx.EVT_SCROLL)
            self.toolbar.update()

        else:

            # Add in the ability to transpose the data in multiplot mode
            self.ax.clear()
            self.twoD_spectra = []
            self.twoD_slices_horizontal = []
            self.twoD_slices_vertical = []
            for i in range(len(self.values_dictionary.keys())):
                self.values_dictionary[i]["new_x_ppms_old"] = self.values_dictionary[i][
                    "new_x_ppms"
                ]
                self.values_dictionary[i]["new_y_ppms_old"] = self.values_dictionary[i][
                    "new_y_ppms"
                ]
                self.values_dictionary[i]["new_x_ppms"] = self.values_dictionary[i][
                    "new_y_ppms_old"
                ]
                self.values_dictionary[i]["new_y_ppms"] = self.values_dictionary[i][
                    "new_x_ppms_old"
                ]
                self.X, self.Y = np.meshgrid(
                    self.values_dictionary[i]["new_y_ppms"],
                    self.values_dictionary[i]["new_x_ppms"],
                )
                self.values_dictionary[i]["z_data_old"] = self.values_dictionary[i][
                    "z_data"
                ]
                try:
                    self.values_dictionary[i]["z_data"] = self.values_dictionary[i][
                        "z_data_old"
                    ].T
                    self.twoD_spectra.append(
                        self.ax.contour(
                            self.Y,
                            self.X,
                            self.values_dictionary[i]["z_data"]
                            * self.values_dictionary[i]["multiply factor"],
                            self.cl,
                            colors=self.cmap,
                            linewidths=self.linewidth,
                        )
                    )
                except:
                    self.values_dictionary[i]["z_data"] = self.values_dictionary[i][
                        "z_data_old"
                    ]
                    self.twoD_spectra.append(
                        self.ax.contour(
                            self.Y,
                            self.X,
                            self.values_dictionary[i]["z_data"]
                            * self.values_dictionary[i]["multiply factor"],
                            self.cl,
                            colors=self.cmap,
                            linewidths=self.linewidth,
                        )
                    )

                self.twoD_slices_horizontal.append(
                    self.axes1D.plot(
                        self.values_dictionary[i]["new_x_ppms"],
                        self.values_dictionary[i]["z_data"][:, 1]
                        * self.values_dictionary[i]["multiply factor"],
                        color=self.twoD_label_colours[i],
                        linewidth=self.values_dictionary[i]["linewidth 1D"],
                    )
                )
                self.twoD_slices_vertical.append(
                    self.axes1D_2.plot(
                        self.values_dictionary[i]["new_y_ppms"],
                        self.values_dictionary[i]["z_data"][1, :]
                        * self.values_dictionary[i]["multiply factor"],
                        color=self.twoD_label_colours[i],
                        linewidth=self.values_dictionary[i]["linewidth 1D"],
                    )
                )

            self.line_h = self.ax.axhline(
                y=self.values_dictionary[i]["new_x_ppms"][1], color="black", lw=1.5
            )
            self.line_v = self.ax.axvline(
                x=self.values_dictionary[i]["new_y_ppms"][1], color="black", lw=1.5
            )
            self.line_h.set_visible(False)
            self.line_v.set_visible(False)

            for i in range(len(self.twoD_slices_horizontal)):
                self.twoD_slices_horizontal[i][0].set_visible(False)
                self.twoD_slices_vertical[i][0].set_visible(False)

            self.ax.set_xlim(
                [
                    max(self.values_dictionary[0]["new_x_ppms"]),
                    min(self.values_dictionary[0]["new_x_ppms"]),
                ]
            )
            self.ax.set_ylim(
                [
                    max(self.values_dictionary[0]["new_y_ppms"]),
                    min(self.values_dictionary[0]["new_y_ppms"]),
                ]
            )
            self.axislabels_old = self.nmrdata.axislabels[0], self.nmrdata.axislabels[1]
            self.nmrdata.axislabels[1] = self.axislabels_old[0]
            self.nmrdata.axislabels[0] = self.axislabels_old[1]
            self.ax.set_xlabel(self.nmrdata.axislabels[1])
            self.ax.set_ylabel(self.nmrdata.axislabels[0])

            # Update all the x and y movement selections and slider values
            for i in range(len(self.values_dictionary.keys())):
                self.move_x_range_index = self.values_dictionary[i][
                    "move x range index"
                ]
                self.move_x_value = self.values_dictionary[i]["move x"]
                self.move_y_range_index = self.values_dictionary[i][
                    "move y range index"
                ]
                self.move_y_value = self.values_dictionary[i]["move y"]
                self.values_dictionary[i][
                    "move x range index"
                ] = self.move_y_range_index
                self.values_dictionary[i]["move x"] = self.move_y_value
                self.values_dictionary[i][
                    "move y range index"
                ] = self.move_x_range_index
                self.values_dictionary[i]["move y"] = self.move_x_value

            self.reference_range_chooserX.SetSelection(
                self.values_dictionary[self.active_plot_index]["move y range index"]
            )
            self.OnReferenceComboX(wx.EVT_SCROLL)
            self.move_x_slider.SetValue(
                self.values_dictionary[self.active_plot_index]["move y"]
            )
            self.reference_range_chooserY.SetSelection(
                self.values_dictionary[self.active_plot_index]["move x range index"]
            )
            self.OnReferenceComboY(wx.EVT_SCROLL)
            self.move_y_slider.SetValue(
                self.values_dictionary[self.active_plot_index]["move x"]
            )

            self.OnMinContour2D(wx.EVT_SCROLL)
            self.toolbar.update()
            titles = []
            for i in range(len(self.values_dictionary.keys())):
                titles.append(self.values_dictionary[i]["title"])

            self.ax.legend(self.files.custom_lines, titles)
            self.UpdateFrame()

    def OnStackButton(self, event):

        if self.multiplot_mode == False:
            # If the number of slices is greater than 30, pop up a window to ask the user if they want to continue
            if len(self.nmrdata.data.T) > 30:
                self.continue_window = wx.MessageDialog(
                    self,
                    "There are "
                    + str(len(self.nmrdata.data.T))
                    + " slices along y axis. Stacking may take a long time. Consider transposing the spectrum and trying again. Do you want to continue?",
                    "Warning",
                    wx.YES_NO | wx.ICON_WARNING,
                )
                if self.continue_window.ShowModal() == wx.ID_NO:
                    self.continue_window.Destroy()
                    return
                else:
                    self.continue_window.Destroy()
            if self.transposed2D == False:
                self.stacks = Stack2D(title="Stacked Slices", parent=self)
            else:
                self.stacks = Stack2D(title="Stacked Slices", parent=self)
        else:
            # Pop up a window to say that this feature is not available in multiplot mode
            self.error_window = wx.MessageDialog(
                self,
                "Stacking is not available in multiplot mode",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            self.error_window.ShowModal()
            self.error_window.Destroy()

    def OnFitDiffusionButton(self, event):
        if self.multiplot_mode == False:
            # If the number of slices is greater than 30, pop up a window to ask the user if they want to continue
            if len(self.nmrdata.data.T) > 30:
                self.continue_window = wx.MessageDialog(
                    self,
                    "There are "
                    + str(len(self.nmrdata.data.T))
                    + " slices along y axis. Consider transposing the spectrum and trying again. Do you want to continue?",
                    "Warning",
                    wx.YES_NO | wx.ICON_WARNING,
                )
                if self.continue_window.ShowModal() == wx.ID_NO:
                    self.continue_window.Destroy()
                    return
                else:
                    self.continue_window.Destroy()
            self.diffusion = DiffusionFit(title="Diffusion Fit", parent=self)
        else:
            # Pop up a window to say that this feature is not available in multiplot mode
            self.error_window = wx.MessageDialog(
                self,
                "Diffusion fitting is not available in multiplot mode",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            self.error_window.ShowModal()
            self.error_window.Destroy()

    def OnCESTButton(self, event):
        if self.multiplot_mode == False:
            # If the number of slices is greater than 30, pop up a window to ask the user if they want to continue
            if len(self.nmrdata.data.T) > 100:
                self.continue_window = wx.MessageDialog(
                    self,
                    "There are "
                    + str(len(self.nmrdata.data.T))
                    + " slices along y axis. Consider transposing the spectrum and trying again. Do you want to continue?",
                    "Warning",
                    wx.YES_NO | wx.ICON_WARNING,
                )
                if self.continue_window.ShowModal() == wx.ID_NO:
                    self.continue_window.Destroy()
                    return
                else:
                    self.continue_window.Destroy()
            # See whether the user has selected a vertical slice
            if self.slice_mode == "y":
                # Give a popout saying this feature will provide the CEST profile for the selected vertical slice
                self.continue_window = wx.MessageDialog(
                    self,
                    "This feature creates the CEST profile for the currently selected vertical slice. Continue?",
                    "Message",
                    wx.YES_NO,
                )
                if self.continue_window.ShowModal() == wx.ID_NO:
                    self.continue_window.Destroy()
                    return
                else:
                    self.continue_window.Destroy()
            else:
                # Give a popout asking the user to select a vertical slice before continuing
                self.continue_window = wx.MessageDialog(
                    self,
                    "Please select a vertical slice by pressing v and clicking on desired location. Then press CEST Analysis again.",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                self.continue_window.ShowModal()
                self.continue_window.Destroy()
                return

            # Ask the user if the CEST data is arrayed as 'On-Resonance, Off-Resonance' or 'Off-Resonance, On-Resonance'
            self.CESTArrayOrder = ""
            self.CESTArray_order_selection = CESTOrder_Dialog(
                title="CEST Array Order", parent=self
            )

        else:
            # Pop up a window to say that this feature is not available in multiplot mode
            self.error_window = wx.MessageDialog(
                self,
                "CEST data plotting is not available in multiplot mode",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            self.error_window.ShowModal()
            self.error_window.Destroy()

    def continue_deletion(self):
        self.CESTArray_order_selection.Destroy()
        if self.CESTArrayOrder == "":
            # User made no selection so return
            self.return_window = wx.MessageDialog(
                self,
                "No selection made. Returning to main window.",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            self.return_window.ShowModal()
            self.return_window.Destroy()
            return

        self.CEST = CESTFrame(
            title="CEST", parent=self, CESTArrayOrder=self.CESTArrayOrder
        )

    def OnFitRelaxButton(self, event):
        if self.multiplot_mode == False:
            # If the number of slices is greater than 30, pop up a window to ask the user if they want to continue
            if len(self.nmrdata.data.T) > 30:
                self.continue_window = wx.MessageDialog(
                    self,
                    "There are "
                    + str(len(self.nmrdata.data.T))
                    + " slices along y axis. Consider transposing the spectrum and trying again. Do you want to continue?",
                    "Warning",
                    wx.YES_NO | wx.ICON_WARNING,
                )
                if self.continue_window.ShowModal() == wx.ID_NO:
                    self.continue_window.Destroy()
                    return
                else:
                    self.continue_window.Destroy()
            self.RelaxFit = RelaxFit(title="Relaxation Fit", parent=self)
        else:
            # Pop up a window to say that this feature is not available in multiplot mode
            self.error_window = wx.MessageDialog(
                self,
                "Relaxation fitting is not available in multiplot mode",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            self.error_window.ShowModal()
            self.error_window.Destroy()

    def OnMinContour2D(self, event):
        # Function to update the contour levels when the user changes the number of contour levels
        self.x_val = 10 ** float(self.contour_slider.GetValue())
        intensity_percent = 10 ** (float(self.intensity_slider.GetValue()))

        if self.multiplot_mode == False:
            # update contour levels
            self.contour_start = np.max(np.abs(self.nmrdata.data)) / self.x_val
            self.cl = self.contour_start * self.contour_factor ** np.arange(
                self.contour_num
            )
            self.cl_neg = -self.contour_start * self.contour_factor ** np.flip(
                np.arange(self.contour_num)
            )
            xlim, ylim = self.ax.get_xlim(), self.ax.get_ylim()
            self.ax.clear()
            self.contour1 = self.ax.contour(
                self.Y,
                self.X,
                self.nmrdata.data * self.multiply_factor,
                self.cl,
                colors=self.cmap,
                linewidths=self.linewidth,
            )
            self.contour1_neg = self.ax.contour(
                self.Y,
                self.X,
                self.nmrdata.data * self.multiply_factor,
                self.cl_neg,
                colors=self.cmap_neg,
                linewidths=self.linewidth,
            )

            if self.line1.get_visible() == True:
                self.line2 = self.ax.axhline(self.y1, color="k")

            if self.line3.get_visible() == True:
                self.line4 = self.ax.axvline(self.x1, color="k")

            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)
            self.ax.set_xlabel(self.nmrdata.axislabels[1])
            self.ax.set_ylabel(self.nmrdata.axislabels[0])
        else:
            self.contour_start = np.max(np.abs(self.nmrdata.data)) / self.x_val
            xlim, ylim = self.ax.get_xlim(), self.ax.get_ylim()
            xlabel, ylabel = self.ax.get_xlabel(), self.ax.get_ylabel()
            self.ax.clear()
            for i in range(len(self.values_dictionary.keys())):
                self.cl = self.contour_start * self.contour_factor ** np.arange(
                    self.values_dictionary[i]["contour levels"]
                )
                multiply_factor = self.values_dictionary[i]["multiply factor"]
                x, y = np.meshgrid(
                    self.values_dictionary[i]["new_y_ppms"],
                    self.values_dictionary[i]["new_x_ppms"],
                )
                self.ax.contour(
                    y,
                    x,
                    self.values_dictionary[i]["z_data"] * multiply_factor,
                    self.cl,
                    colors=self.twoD_colours[i],
                    linewidths=self.values_dictionary[i]["contour linewidth"],
                )
                self.ax.legend(self.files.custom_lines, self.files.custom_labels)

            if self.twoD_slices_horizontal[0][0].get_visible() == True:
                # for i in range(len(self.twoD_slices_horizontal)):
                self.line_h = self.ax.axhline(self.y1, color="k")
                self.axes1D.set_ylim(
                    -(np.max(self.nmrdata.data) / 8) / (intensity_percent / 100),
                    np.max(self.nmrdata.data) / (intensity_percent / 100),
                )

            if self.twoD_slices_vertical[0][0].get_visible() == True:
                # for i in range(len(self.twoD_slices_horizontal)):
                self.line_v = self.ax.axhline(self.x1, color="k")
                self.axes1D_2.set_ylim(
                    -(np.max(self.nmrdata.data) / 8) / (intensity_percent / 100),
                    np.max(self.nmrdata.data) / (intensity_percent / 100),
                )

            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)
            self.ax.set_xlabel(xlabel)
            self.ax.set_ylabel(ylabel)

        self.contour_value_label.SetLabel(
            str(int(10 ** float(self.contour_slider.GetValue())))
        )
        self.OnSliderScroll2D(wx.EVT_SCROLL)
        self.OnIntensityScroll2D(wx.EVT_SCROLL)

        self.UpdateFrame()

    def OnContourLevels(self, event):
        # update number of contour levels
        self.contour_num = round(float(self.contour_levels_slider.GetValue()))
        if self.multiplot_mode == True:
            if self.select_all_checkbox.IsChecked() == False:
                self.values_dictionary[self.active_plot_index][
                    "contour levels"
                ] = self.contour_num
            else:
                for i in range(len(self.twoD_slices_horizontal)):
                    self.values_dictionary[i]["contour levels"] = self.contour_num
        self.OnMinContour2D(event)

    def OnMultiplyScroll2D(self, event):
        self.multiply_factor = float(self.multiply_slider.GetValue())
        self.multiply_value_label.SetLabel(
            "{:.2f}".format(float(self.multiply_slider.GetValue()))
        )
        if self.multiplot_mode == True:
            if self.select_all_checkbox.IsChecked() == False:
                self.values_dictionary[self.active_plot_index][
                    "multiply factor"
                ] = self.multiply_factor
            else:
                for i in range(len(self.twoD_slices_horizontal)):
                    self.values_dictionary[i]["multiply factor"] = self.multiply_factor
        self.OnMinContour2D(event)

    def OnMultiplyCombo2D(self, event):
        self.multiply_factor = float(self.multiply_slider.GetValue())
        self.multiply_range = float(self.multiply_range_chooser2d.GetValue())
        self.multiply_slider.SetMax(self.multiply_range)
        self.multiply_slider.SetValue(self.multiply_factor)

        if self.multiplot_mode == True:
            if self.select_all_checkbox.IsChecked() == False:
                self.values_dictionary[self.active_plot_index][
                    "multiply factor"
                ] = self.multiply_factor
            else:
                for i in range(len(self.twoD_slices_horizontal)):
                    self.values_dictionary[i]["multiply factor"] = self.multiply_factor
        self.OnMinContour2D(event)

    def OnContourWidth(self, event):
        # update contour linewidth
        self.linewidth = float(self.contour_width_slider.GetValue())
        if self.multiplot_mode == True:
            if self.select_all_checkbox.IsChecked() == False:
                self.values_dictionary[self.active_plot_index][
                    "contour linewidth"
                ] = self.linewidth
            else:
                for i in range(len(self.twoD_slices_horizontal)):
                    self.values_dictionary[i]["contour linewidth"] = self.linewidth

        self.OnMinContour2D(event)

    def On2DLinewidth(self, event):
        # update contour linewidth
        self.linewidth1D = float(self.line_width_slider.GetValue())
        if self.multiplot_mode == True:
            if self.select_all_checkbox.IsChecked() == False:
                self.values_dictionary[self.active_plot_index][
                    "linewidth 1D"
                ] = self.linewidth1D
                self.twoD_slices_horizontal[self.active_plot_index][0].set_linewidth(
                    self.linewidth1D
                )
                self.twoD_slices_vertical[self.active_plot_index][0].set_linewidth(
                    self.linewidth1D
                )
            else:
                for i in range(len(self.twoD_slices_horizontal)):
                    self.values_dictionary[i]["linewidth 1D"] = self.linewidth1D
                    self.twoD_slices_horizontal[i][0].set_linewidth(self.linewidth1D)
                    self.twoD_slices_vertical[i][0].set_linewidth(self.linewidth1D)
        else:
            self.line1.set_linewidth(self.linewidth1D)
            self.line3.set_linewidth(self.linewidth1D)
        self.UpdateFrame()

    def OnMoveX(self, event):
        # update x-axis
        self.x_movement = float(self.move_x_slider.GetValue())
        self.move_x_value_label.SetLabel("{:.4f}".format(self.x_movement))
        if self.multiplot_mode == False:
            if self.transposed2D == False:
                self.new_x_ppms = (
                    self.ppms_0 + np.ones(len(self.ppms_0)) * self.x_movement
                )
            else:
                self.new_x_ppms = (
                    self.ppms_1 + np.ones(len(self.ppms_1)) * self.x_movement
                )
            self.X, self.Y = np.meshgrid(self.new_y_ppms, self.new_x_ppms)
            self.OnMinContour2D(wx.EVT_SCROLL)
            self.UpdateFrame()

        else:
            if self.transposed2D == False:
                if self.select_all_checkbox.IsChecked() == False:
                    self.values_dictionary[self.active_plot_index]["new_x_ppms"] = (
                        self.values_dictionary[self.active_plot_index][
                            "original_x_ppms"
                        ]
                        + np.ones(
                            len(
                                self.values_dictionary[self.active_plot_index][
                                    "original_x_ppms"
                                ]
                            )
                        )
                        * self.x_movement
                    )
                    self.values_dictionary[self.active_plot_index][
                        "move x"
                    ] = self.x_movement
                else:
                    for i in range(len(self.twoD_slices_horizontal)):
                        self.values_dictionary[i]["new_x_ppms"] = (
                            self.values_dictionary[i]["original_x_ppms"]
                            + np.ones(len(self.values_dictionary[i]["original_x_ppms"]))
                            * self.x_movement
                        )
                        self.values_dictionary[i]["move x"] = self.x_movement
            else:
                if self.select_all_checkbox.IsChecked() == False:
                    self.values_dictionary[self.active_plot_index]["new_x_ppms"] = (
                        self.values_dictionary[self.active_plot_index][
                            "original_y_ppms"
                        ]
                        + np.ones(
                            len(
                                self.values_dictionary[self.active_plot_index][
                                    "original_y_ppms"
                                ]
                            )
                        )
                        * self.x_movement
                    )
                    self.values_dictionary[self.active_plot_index][
                        "move x"
                    ] = self.x_movement
                else:
                    for i in range(len(self.twoD_slices_horizontal)):
                        self.values_dictionary[i]["new_x_ppms"] = (
                            self.values_dictionary[i]["original_y_ppms"]
                            + np.ones(len(self.values_dictionary[i]["original_y_ppms"]))
                            * self.x_movement
                        )
                        self.values_dictionary[i]["move x"] = self.x_movement
            self.OnMinContour2D(wx.EVT_SCROLL)
            self.UpdateFrame()

    def OnMoveY(self, event):
        # update y-axis
        self.y_movement = float(self.move_y_slider.GetValue())
        self.move_y_value_label.SetLabel("{:.4f}".format(self.y_movement))
        if self.multiplot_mode == False:
            if self.transposed2D == False:
                self.new_y_ppms = (
                    self.ppms_1 + np.ones(len(self.ppms_1)) * self.y_movement
                )
            else:
                self.new_y_ppms = (
                    self.ppms_0 + np.ones(len(self.ppms_0)) * self.y_movement
                )
            self.X, self.Y = np.meshgrid(self.new_y_ppms, self.new_x_ppms)
            self.OnMinContour2D(wx.EVT_SCROLL)
            self.UpdateFrame()

        else:
            if self.transposed2D == False:
                if self.select_all_checkbox.IsChecked() == False:
                    self.values_dictionary[self.active_plot_index]["new_y_ppms"] = (
                        self.values_dictionary[self.active_plot_index][
                            "original_y_ppms"
                        ]
                        + np.ones(
                            len(
                                self.values_dictionary[self.active_plot_index][
                                    "original_y_ppms"
                                ]
                            )
                        )
                        * self.y_movement
                    )
                    self.values_dictionary[self.active_plot_index][
                        "move y"
                    ] = self.y_movement
                else:
                    for i in range(len(self.twoD_slices_horizontal)):
                        self.values_dictionary[i]["new_y_ppms"] = (
                            self.values_dictionary[i]["original_y_ppms"]
                            + np.ones(len(self.values_dictionary[i]["original_y_ppms"]))
                            * self.y_movement
                        )
                        self.values_dictionary[i]["move y"] = self.y_movement
            else:
                if self.select_all_checkbox.IsChecked() == False:
                    self.values_dictionary[self.active_plot_index]["new_y_ppms"] = (
                        self.values_dictionary[self.active_plot_index][
                            "original_x_ppms"
                        ]
                        + np.ones(
                            len(
                                self.values_dictionary[self.active_plot_index][
                                    "original_x_ppms"
                                ]
                            )
                        )
                        * self.y_movement
                    )
                    self.values_dictionary[self.active_plot_index][
                        "move y"
                    ] = self.y_movement
                else:
                    for i in range(len(self.twoD_slices_horizontal)):
                        self.values_dictionary[i]["new_y_ppms"] = (
                            self.values_dictionary[i]["original_x_ppms"]
                            + np.ones(len(self.values_dictionary[i]["original_x_ppms"]))
                            * self.y_movement
                        )
                        self.values_dictionary[i]["move y"] = self.y_movement
            self.OnMinContour2D(wx.EVT_SCROLL)
            self.UpdateFrame()

    def OnReferenceComboX(self, event):
        # Change the range for the move-x slider
        index = int(self.reference_range_chooserX.GetSelection())
        self.reference_rangeX = float(self.reference_range_values[index])
        if self.multiplot_mode == True:
            if self.select_all_checkbox.IsChecked() == False:
                self.values_dictionary[self.active_plot_index][
                    "move x range index"
                ] = index
            else:
                for i in range(len(self.twoD_slices_horizontal)):
                    self.values_dictionary[i]["move x range index"] = index
        self.move_x_slider.SetMin(-self.reference_rangeX)
        self.move_x_slider.SetMax(self.reference_rangeX)
        self.move_x_slider.SetRes(self.reference_rangeX / 1000)
        self.move_x_slider.Bind(wx.EVT_SLIDER, self.OnMoveX)

    def OnReferenceComboY(self, event):
        # Change the range for the move-y slider
        index = int(self.reference_range_chooserY.GetSelection())
        self.reference_rangeY = float(self.reference_range_values[index])
        if self.multiplot_mode == True:
            if self.select_all_checkbox.IsChecked() == False:
                self.values_dictionary[self.active_plot_index][
                    "move y range index"
                ] = index
            else:
                for i in range(len(self.twoD_slices_horizontal)):
                    self.values_dictionary[i]["move y range index"] = index
        self.move_y_slider.SetMin(-self.reference_rangeY)
        self.move_y_slider.SetMax(self.reference_rangeY)
        self.move_y_slider.SetRes(self.reference_rangeY / 1000)
        self.move_y_slider.Bind(wx.EVT_SLIDER, self.OnMoveY)

    def on_key_2d(self, event):
        # key press event for 2D plot (Plot horizontal and vertical slices)
        if self.multiplot_mode == False:
            if event.key == "h":
                self.axes1D.set_ylim(
                    -np.max(self.nmrdata.data / 8), np.max(self.nmrdata.data)
                )
                # plot a horizontal slice of the data
                if self.line1.get_visible() == True:
                    self.slice_mode = None
                    self.line1.set_visible(False)
                    self.line2.set_visible(False)
                    self.UpdateFrame()

                else:
                    if self.line3.get_visible() == True:
                        self.slice_mode = None
                        self.line3.set_visible(False)
                        self.line4.set_visible(False)
                        self.UpdateFrame()

                    else:
                        self.slice_mode = "x"
                        (self.line1,) = self.axes1D.plot(
                            self.new_x_ppms,
                            self.nmrdata.data[
                                :, self.uc1(str(self.new_y_ppms[1]) + "ppm")
                            ]
                            * self.multiply_factor,
                            color=self.slice_colour,
                        )
                        self.line2 = self.ax.axhline(self.new_y_ppms[1], color="k")
                        self.UpdateFrame()

            if event.key == "v":
                self.axes1D_2.set_xlim(
                    -np.max(self.nmrdata.data / 8), np.max(self.nmrdata.data)
                )
                if self.line3.get_visible() == True:
                    self.slice_mode = None
                    self.line3.set_visible(False)
                    self.line4.set_visible(False)
                    self.UpdateFrame()
                else:
                    if self.line1.get_visible() == True:
                        self.slice_mode = None
                        self.line1.set_visible(False)
                        self.line2.set_visible(False)
                        self.UpdateFrame()
                    else:
                        self.line3.set_visible = True
                        self.line4.set_visible = True
                        (self.line3,) = self.axes1D_2.plot(
                            self.nmrdata.data[
                                self.uc0(str(self.new_x_ppms[1]) + "ppm"), :
                            ]
                            * self.multiply_factor,
                            self.new_y_ppms,
                            color=self.slice_colour,
                        )
                        self.line4 = self.ax.axvline(self.new_x_ppms[1], color="k")
                        self.slice_mode = "y"
                        self.UpdateFrame()

        else:
            if event.key == "h":
                self.axes1D.set_ylim(
                    -np.max(self.nmrdata.data / 8), np.max(self.nmrdata.data)
                )
                # plot a horizontal slice of the data
                if self.twoD_slices_horizontal[0][0].get_visible() == True:
                    for i in range(len(self.twoD_slices_horizontal)):
                        self.twoD_slices_horizontal[i][0].set_visible(False)
                    self.line_h.set_visible(False)
                    self.UpdateFrame()
                else:
                    if self.twoD_slices_vertical[0][0].get_visible() == True:
                        for i in range(len(self.twoD_slices_vertical)):
                            self.twoD_slices_vertical[i][0].set_visible(False)
                        self.line_v.set_visible(False)
                        self.UpdateFrame()
                    else:
                        for i in range(len(self.twoD_slices_horizontal)):
                            multiply_factor = self.values_dictionary[i][
                                "multiply factor"
                            ]
                            try:
                                self.twoD_slices_horizontal[i] = self.axes1D.plot(
                                    self.values_dictionary[i]["new_x_ppms"],
                                    self.values_dictionary[i]["z_data"][
                                        :,
                                        self.values_dictionary[i]["uc1"](
                                            str(self.new_y_ppms[1]) + "ppm"
                                        ),
                                    ]
                                    * multiply_factor,
                                    color=self.twoD_label_colours[i],
                                    linewidth=self.values_dictionary[i]["linewidth 1D"],
                                )
                            except:
                                self.twoD_slices_horizontal[i] = self.axes1D.plot(
                                    self.values_dictionary[i]["new_x_ppms"],
                                    self.values_dictionary[i]["z_data"][
                                        :,
                                        self.values_dictionary[i]["uc0"](
                                            str(self.new_y_ppms[1]) + "ppm"
                                        ),
                                    ]
                                    * multiply_factor,
                                    color=self.twoD_label_colours[i],
                                    linewidth=self.values_dictionary[i]["linewidth 1D"],
                                )
                        self.line_h = self.ax.axhline(self.new_y_ppms[1], color="k")
                        self.UpdateFrame()

            if event.key == "v":
                self.axes1D_2.set_xlim(
                    -np.max(self.nmrdata.data / 8), np.max(self.nmrdata.data)
                )
                if self.twoD_slices_vertical[0][0].get_visible() == True:
                    for i in range(len(self.twoD_slices_vertical)):
                        self.twoD_slices_vertical[i][0].set_visible(False)
                    self.line_v.set_visible(False)
                    self.UpdateFrame()
                else:
                    if self.twoD_slices_horizontal[0][0].get_visible() == True:
                        for i in range(len(self.twoD_slices_horizontal)):
                            self.twoD_slices_horizontal[i][0].set_visible(False)
                        self.line_h.set_visible(False)
                        self.UpdateFrame()
                    else:
                        for i in range(len(self.twoD_slices_vertical)):
                            multiply_factor = self.values_dictionary[i][
                                "multiply factor"
                            ]
                            try:
                                self.twoD_slices_vertical[i] = self.axes1D_2.plot(
                                    self.values_dictionary[i]["z_data"][
                                        self.values_dictionary[i]["uc0"](
                                            str(self.new_x_ppms[1]) + "ppm"
                                        ),
                                        :,
                                    ]
                                    * multiply_factor,
                                    self.values_dictionary[i]["new_y_ppms"],
                                    color=self.twoD_label_colours[i],
                                    linewidth=self.values_dictionary[i]["linewidth 1D"],
                                )
                            except:
                                self.twoD_slices_vertical[i] = self.axes1D_2.plot(
                                    self.values_dictionary[i]["z_data"][
                                        self.values_dictionary[i]["uc1"](
                                            str(self.new_x_ppms[1]) + "ppm"
                                        ),
                                        :,
                                    ]
                                    * multiply_factor,
                                    self.values_dictionary[i]["new_y_ppms"],
                                    color=self.twoD_label_colours[i],
                                    linewidth=self.values_dictionary[i]["linewidth 1D"],
                                )
                        self.line_v = self.ax.axvline(self.new_x_ppms[1], color="k")
                        self.UpdateFrame()

    def on_click_2d(self, event):

        # mouse click event for 2D plot (Plot horizontal and vertical slices for given mouse position on-click)

        self.x1, self.y1 = self.ax.transData.inverted().transform((event.x, event.y))

        if self.x1 != None and self.y1 != None:

            if self.multiplot_mode == False:

                if self.line1.get_visible() == True:
                    self.line1.set_ydata(
                        self.nmrdata.data[
                            :, self.uc1(str(self.y1 - self.y_movement) + "ppm")
                        ]
                        * self.multiply_factor
                    )
                    self.line2.set_ydata([self.y1])
                    self.line1.set_xdata(self.ppms_0 + self.x_movement)
                    self.OnSliderScroll2D(wx.EVT_SCROLL)
                    self.UpdateFrame()
                if self.line3.get_visible() == True:
                    self.line3.set_xdata(
                        self.nmrdata.data[
                            self.uc0(str(self.x1 - self.x_movement) + "ppm"), :
                        ]
                        * self.multiply_factor
                    )
                    self.line4.set_xdata([self.x1])
                    self.line3.set_ydata(self.ppms_1 + self.y_movement)
                    self.OnSliderScroll2D(wx.EVT_SCROLL)
                    self.UpdateFrame()

            else:
                if self.twoD_slices_horizontal[0][0].get_visible() == True:
                    for i in range(len(self.twoD_slices_horizontal)):
                        multiply_factor = self.values_dictionary[i]["multiply factor"]
                        self.y_difference = self.values_dictionary[i]["move y"]
                        try:
                            if self.transposed2D == False:
                                self.twoD_slices_horizontal[i][0].set_ydata(
                                    self.values_dictionary[i]["z_data"][
                                        :,
                                        self.values_dictionary[i]["uc1"](
                                            str(self.y1 - self.y_difference) + "ppm"
                                        ),
                                    ]
                                    * multiply_factor
                                )
                                self.twoD_slices_horizontal[i][0].set_xdata(
                                    self.values_dictionary[i]["new_x_ppms"]
                                )
                            else:
                                self.twoD_slices_horizontal[i][0].set_ydata(
                                    self.values_dictionary[i]["z_data"][
                                        :,
                                        self.values_dictionary[i]["uc0"](
                                            str(self.y1 - self.y_difference) + "ppm"
                                        ),
                                    ]
                                    * multiply_factor
                                )
                                self.twoD_slices_horizontal[i][0].set_xdata(
                                    self.values_dictionary[i]["new_x_ppms"]
                                )
                        except:
                            self.twoD_slices_vertical[i][0].set_xdata(
                                0
                                * np.ones(
                                    len(
                                        self.values_dictionary[i]["z_data"][:, 0]
                                        * multiply_factor
                                    )
                                )
                            )
                            self.twoD_slices_vertical[i][0].set_ydata(
                                0
                                * np.ones(len(self.values_dictionary[i]["new_x_ppms"]))
                            )
                    self.line_h.set_ydata([self.y1])
                    self.OnSliderScroll2D(wx.EVT_SCROLL)
                    self.UpdateFrame()
                if self.twoD_slices_vertical[0][0].get_visible() == True:
                    for i in range(len(self.twoD_slices_vertical)):
                        multiply_factor = self.values_dictionary[i]["multiply factor"]
                        self.x_difference = self.values_dictionary[i]["move x"]
                        try:
                            if self.transposed2D == False:
                                self.twoD_slices_vertical[i][0].set_xdata(
                                    self.values_dictionary[i]["z_data"][
                                        self.values_dictionary[i]["uc0"](
                                            str(self.x1 - self.x_difference) + "ppm"
                                        ),
                                        :,
                                    ]
                                    * multiply_factor
                                )
                                self.twoD_slices_vertical[i][0].set_ydata(
                                    self.values_dictionary[i]["new_y_ppms"]
                                )
                            else:
                                self.twoD_slices_vertical[i][0].set_xdata(
                                    self.values_dictionary[i]["z_data"][
                                        self.values_dictionary[i]["uc1"](
                                            str(self.x1 - self.x_difference) + "ppm"
                                        ),
                                        :,
                                    ]
                                    * multiply_factor
                                )
                                self.twoD_slices_vertical[i][0].set_ydata(
                                    self.values_dictionary[i]["new_y_ppms"]
                                )
                        except:
                            self.twoD_slices_vertical[i][0].set_xdata(
                                0
                                * np.ones(
                                    len(
                                        self.values_dictionary[i]["z_data"][0, :]
                                        * multiply_factor
                                    )
                                )
                            )
                            self.twoD_slices_vertical[i][0].set_ydata(
                                0
                                * np.ones(len(self.values_dictionary[i]["new_y_ppms"]))
                            )
                    self.line_v.set_xdata([self.x1])
                    self.OnSliderScroll2D(wx.EVT_SCROLL)
                    self.UpdateFrame()

    def OnSliderScroll2D(self, event):
        # Get all the slider values for P0 and P1 (coarse and fine), put the combined coarse and fine values on the screen
        self.total_P0 = self.P0_slider.GetValue() + self.P0_slider_fine.GetValue()
        self.total_P1 = self.P1_slider.GetValue() + self.P1_slider_fine.GetValue()
        self.P0_total_value.SetLabel("{:.2f}".format(self.total_P0))
        self.P1_total_value.SetLabel("{:.2f}".format(self.total_P1))
        self.phase2D()

    def phase2D(self):
        # Phase the 2D data with the combined coarse/fine phasing values and plot the result
        if self.multiplot_mode == False:
            try:
                if self.line1.get_visible() == True:
                    data = (
                        self.nmrdata.data[
                            :, self.uc1(str(self.y1 - self.y_movement) + "ppm")
                        ]
                        * self.multiply_factor
                    )
                    complex_data = ng.process.proc_base.ht(
                        data, self.nmrdata.data.shape[0]
                    )
                    phased_data = complex_data * np.exp(
                        1j
                        * (
                            self.total_P0 * np.pi / 180
                            + self.total_P1
                            * (np.pi / 180)
                            * (
                                np.arange(
                                    -self.pivot_x,
                                    -self.pivot_x + self.nmrdata.data.shape[0],
                                )
                                / self.nmrdata.data.shape[0]
                            )
                        )
                    )
                    # phased_data = ng.process.proc_base.ps(complex_data, p0=self.total_P0, p1=self.total_P1)
                    self.line1.set_ydata(phased_data)
                    self.line1.set_xdata(self.new_x_ppms)
                    self.line1.set_linewidth(self.linewidth1D)
                    # self.line2 = self.ax.axhline(self.y1, color='k')
                    self.UpdateFrame()
                if self.line3.get_visible() == True:
                    data = (
                        self.nmrdata.data[
                            self.uc0(str(self.x1 - self.x_movement) + "ppm"), :
                        ]
                        * self.multiply_factor
                    )
                    complex_data = ng.process.proc_base.ht(
                        data, self.nmrdata.data.shape[1]
                    )
                    phased_data = complex_data * np.exp(
                        1j
                        * (
                            self.total_P0 * np.pi / 180
                            + self.total_P1
                            * (np.pi / 180)
                            * (
                                np.arange(
                                    -self.pivot_y,
                                    -self.pivot_y + self.nmrdata.data.shape[1],
                                )
                                / self.nmrdata.data.shape[1]
                            )
                        )
                    )
                    # phased_data = ng.process.proc_base.ps(complex_data, p0=self.total_P0, p1=self.total_P1)
                    self.line3.set_xdata(phased_data)
                    self.line3.set_ydata(self.new_y_ppms)
                    self.line3.set_linewidth(self.linewidth1D)
                    # self.line4 = self.ax.axvline(self.x1, color='k')
                    self.UpdateFrame()
            except:
                self.OnTransposeButton(wx.EVT_BUTTON)
                self.OnSliderScroll2D(wx.EVT_SCROLL)
                # Give a pop-up window to say that transposing is not supported whilst horizontal or vertical slices are plotted
                self.error_window = wx.MessageDialog(
                    self,
                    "Transposing is not supported whilst horizontal or vertical slices are plotted.",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                self.error_window.ShowModal()
                self.error_window.Destroy()

        else:
            if self.twoD_slices_horizontal[0][0].get_visible() == True:
                if self.select_all_checkbox.IsChecked() == False:
                    multiply_factor = self.values_dictionary[self.active_plot_index][
                        "multiply factor"
                    ]
                    self.values_dictionary[self.active_plot_index][
                        "p0 Coarse"
                    ] = self.P0_slider.GetValue()
                    self.values_dictionary[self.active_plot_index][
                        "p0 Fine"
                    ] = self.P0_slider_fine.GetValue()
                    self.values_dictionary[self.active_plot_index][
                        "p1 Coarse"
                    ] = self.P1_slider.GetValue()
                    self.values_dictionary[self.active_plot_index][
                        "p1 Fine"
                    ] = self.P1_slider_fine.GetValue()
                    if self.transposed2D == False:
                        data = (
                            self.values_dictionary[self.active_plot_index]["z_data"][
                                :,
                                self.values_dictionary[self.active_plot_index]["uc1"](
                                    str(self.y1 - self.y_difference) + "ppm"
                                ),
                            ]
                            * multiply_factor
                        )
                    else:
                        data = (
                            self.values_dictionary[self.active_plot_index]["z_data"][
                                :,
                                self.values_dictionary[self.active_plot_index]["uc0"](
                                    str(self.y1 - self.y_difference) + "ppm"
                                ),
                            ]
                            * multiply_factor
                        )

                    complex_data = ng.process.proc_base.ht(
                        data,
                        self.values_dictionary[self.active_plot_index]["z_data"].shape[
                            0
                        ],
                    )
                    phased_data = ng.process.proc_base.ps(
                        complex_data, p0=self.total_P0, p1=self.total_P1
                    )
                    self.twoD_slices_horizontal[self.active_plot_index][0].set_ydata(
                        phased_data
                    )
                    self.twoD_slices_horizontal[self.active_plot_index][0].set_xdata(
                        self.values_dictionary[self.active_plot_index]["new_x_ppms"]
                    )
                    self.twoD_slices_horizontal[self.active_plot_index][
                        0
                    ].set_linewidth(
                        self.values_dictionary[self.active_plot_index]["linewidth 1D"]
                    )
                else:
                    for i in range(len(self.twoD_slices_horizontal)):
                        multiply_factor = self.values_dictionary[i]["multiply factor"]
                        self.values_dictionary[i][
                            "p0 Coarse"
                        ] = self.P0_slider.GetValue()
                        self.values_dictionary[i][
                            "p0 Fine"
                        ] = self.P0_slider_fine.GetValue()
                        self.values_dictionary[i][
                            "p1 Coarse"
                        ] = self.P1_slider.GetValue()
                        self.values_dictionary[i][
                            "p1 Fine"
                        ] = self.P1_slider_fine.GetValue()
                        if self.transposed2D == False:
                            data = (
                                self.values_dictionary[i]["z_data"][
                                    :,
                                    self.values_dictionary[i]["uc1"](
                                        str(self.y1 - self.y_difference) + "ppm"
                                    ),
                                ]
                                * multiply_factor
                            )
                        else:
                            data = (
                                self.values_dictionary[i]["z_data"][
                                    :,
                                    self.values_dictionary[i]["uc0"](
                                        str(self.y1 - self.y_difference) + "ppm"
                                    ),
                                ]
                                * multiply_factor
                            )
                        complex_data = ng.process.proc_base.ht(
                            data, self.values_dictionary[i]["z_data"].shape[0]
                        )
                        phased_data = ng.process.proc_base.ps(
                            complex_data, p0=self.total_P0, p1=self.total_P1
                        )
                        self.twoD_slices_horizontal[i][0].set_ydata(phased_data)
                        self.twoD_slices_horizontal[i][0].set_xdata(
                            self.values_dictionary[i]["new_x_ppms"]
                        )
                        self.twoD_slices_horizontal[i][0].set_linewidth(
                            self.values_dictionary[i]["linewidth 1D"]
                        )

                self.UpdateFrame()
            if self.twoD_slices_vertical[0][0].get_visible() == True:
                if self.select_all_checkbox.IsChecked() == False:
                    multiply_factor = self.values_dictionary[self.active_plot_index][
                        "multiply factor"
                    ]
                    self.x_difference = self.values_dictionary[self.active_plot_index][
                        "move x"
                    ]
                    if self.transposed2D == False:
                        data = (
                            self.values_dictionary[self.active_plot_index]["z_data"][
                                self.values_dictionary[self.active_plot_index]["uc0"](
                                    str(self.x1 - self.x_difference) + "ppm"
                                ),
                                :,
                            ]
                            * multiply_factor
                        )
                    else:
                        data = (
                            self.values_dictionary[self.active_plot_index]["z_data"][
                                self.values_dictionary[self.active_plot_index]["uc1"](
                                    str(self.x1 - self.x_difference) + "ppm"
                                ),
                                :,
                            ]
                            * multiply_factor
                        )
                    complex_data = ng.process.proc_base.ht(
                        data,
                        self.values_dictionary[self.active_plot_index]["z_data"].shape[
                            1
                        ],
                    )
                    phased_data = ng.process.proc_base.ps(
                        complex_data, p0=self.total_P0, p1=self.total_P1
                    )
                    self.twoD_slices_vertical[self.active_plot_index][0].set_xdata(
                        phased_data
                    )
                    self.twoD_slices_vertical[self.active_plot_index][0].set_ydata(
                        self.values_dictionary[self.active_plot_index]["new_y_ppms"]
                    )
                    self.twoD_slices_vertical[self.active_plot_index][0].set_linewidth(
                        self.values_dictionary[self.active_plot_index]["linewidth 1D"]
                    )
                else:
                    for i in range(len(self.twoD_slices_vertical)):
                        multiply_factor = self.values_dictionary[i]["multiply factor"]
                        self.values_dictionary[i][
                            "p0 Coarse"
                        ] = self.P0_slider.GetValue()
                        self.values_dictionary[i][
                            "p0 Fine"
                        ] = self.P0_slider_fine.GetValue()
                        self.values_dictionary[i][
                            "p1 Coarse"
                        ] = self.P1_slider.GetValue()
                        self.values_dictionary[i][
                            "p1 Fine"
                        ] = self.P1_slider_fine.GetValue()
                        self.x_difference = self.values_dictionary[i]["move x"]
                        if self.transposed2D == False:
                            data = (
                                self.values_dictionary[i]["z_data"][
                                    self.values_dictionary[i]["uc0"](
                                        str(self.x1 - self.x_difference) + "ppm"
                                    ),
                                    :,
                                ]
                                * multiply_factor
                            )
                        else:
                            data = (
                                self.values_dictionary[i]["z_data"][
                                    self.values_dictionary[i]["uc1"](
                                        str(self.x1 - self.x_difference) + "ppm"
                                    ),
                                    :,
                                ]
                                * multiply_factor
                            )
                        complex_data = ng.process.proc_base.ht(
                            data, self.values_dictionary[i]["z_data"].shape[1]
                        )
                        phased_data = ng.process.proc_base.ps(
                            complex_data, p0=self.total_P0, p1=self.total_P1
                        )
                        self.twoD_slices_vertical[i][0].set_xdata(phased_data)
                        self.twoD_slices_vertical[i][0].set_ydata(
                            self.values_dictionary[i]["new_y_ppms"]
                        )
                        self.twoD_slices_vertical[i][0].set_linewidth(
                            self.values_dictionary[i]["linewidth 1D"]
                        )
                self.UpdateFrame()

    def OnIntensityScroll2D(self, event):

        # Change the y-axis limits of the 1D slices in the 2D plot
        intensity_percent = 10 ** (float(self.intensity_slider.GetValue()))

        if self.multiplot_mode == False:
            if self.line1.get_visible() == True:
                self.axes1D.set_ylim(
                    -(np.max(self.nmrdata.data) / 8) / (intensity_percent / 100),
                    np.max(self.nmrdata.data) / (intensity_percent / 100),
                )
                self.UpdateFrame()
            if self.line3.get_visible() == True:
                self.axes1D_2.set_xlim(
                    -(np.max(self.nmrdata.data) / 8) / (intensity_percent / 100),
                    np.max(self.nmrdata.data) / (intensity_percent / 100),
                )
                self.UpdateFrame()
        else:
            if self.twoD_slices_horizontal[0][0].get_visible() == True:
                self.axes1D.set_ylim(
                    -(np.max(self.nmrdata.data) / 8) / (intensity_percent / 100),
                    np.max(self.nmrdata.data) / (intensity_percent / 100),
                )
                self.UpdateFrame()
            if self.twoD_slices_vertical[0][0].get_visible() == True:
                self.axes1D_2.set_xlim(
                    -(np.max(self.nmrdata.data) / 8) / (intensity_percent / 100),
                    np.max(self.nmrdata.data) / (intensity_percent / 100),
                )
                self.UpdateFrame()


# The class for viewing 3D spectra
class ThreeDViewer(wx.Panel):
    def __init__(self, parent, nmrdata):
        # Get the monitor size and set the window size to 85% of the monitor size
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        self.display_index = wx.Display.GetFromWindow(parent)
        self.width = int(1.0 * sizes[self.display_index][0])
        self.height = int(0.875 * sizes[self.display_index][1])
        self.parent = parent
        wx.Panel.__init__(self, parent, id=wx.ID_ANY, size=(self.width, self.height))
        self.nmrdata = nmrdata

        self.set_initial_variables_3D()
        self.create_button_panel_3D()
        self.create_hidden_button_panel_3D()
        self.create_canvas_3D()
        self.add_to_main_sizer_3D()
        self.draw_figure_3D()

    def add_to_main_sizer_3D(self):
        # Create the main sizer
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.main_sizer.Add(self.canvas, 10, wx.EXPAND)
        self.main_sizer.Add(self.toolbar, 0, wx.EXPAND)
        self.main_sizer.Add(self.bottom_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.main_sizer.Add(self.show_button_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.main_sizer.Hide(self.show_button_sizer)
        self.SetSizer(self.main_sizer)

    def create_hidden_button_panel_3D(self):
        # Create a button to show the options
        self.show_button = wx.Button(self, label="Show Options")
        self.show_button.Bind(wx.EVT_BUTTON, self.OnHideButton)
        self.show_button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.show_button_sizer.Add(self.show_button, wx.ALIGN_CENTER, 5)
        self.show_button_sizer.AddSpacer(5)

    def create_canvas_3D(self):
        # Create the figure and canvas to draw on
        self.panel = wx.Panel(self)
        self.fig = Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        self.toolbar = NavigationToolbar(self.canvas)

    def set_initial_variables_3D(self):
        # Colours for 1D lines
        self.colours = [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
        ]
        self.colour_value = self.colours[0]

        # Initial 1D slice colour for 2D/3D spectra is set to navy
        self.colour_slice = "navy"

        # List of cmap colours for when overlaying multiple spectra
        self.cmap = "#e41a1c"
        self.cmap_neg = "#377eb8"
        self.twoD_colours = [
            "#e41a1c",
            "#377eb8",
            "#4daf4a",
            "#984ea3",
            "#ff7f00",
            "#ff33eb",
        ]
        self.twoD_label_colours = self.twoD_colours

        self.twoD_slices_horizontal = []
        self.twoD_slices_vertical = []

        # Range of the sliders to for moving spectra left/right/up/down
        self.reference_range_values = [
            "0.01",
            "0.1",
            "0.5",
            "1.0",
            "5.0",
            "10.0",
            "50.0",
        ]
        self.reference_range = float(self.reference_range_values[0])
        self.reference_rangeX = float(self.reference_range_values[0])
        self.reference_rangeY = float(self.reference_range_values[0])

        # Range of the sliders to for moving spectra up/down in 1D spectra
        self.vertical_range_values = [
            "0.01",
            "0.1",
            "0.5",
            "1.0",
            "10.0",
            "50.0",
            "100.0",
            "1000.0",
            "10000",
        ]

        # Range of the sliders to for multiplying 1D spectra
        self.multiply_range_values = [
            "1.01",
            "1.1",
            "1.5",
            "2",
            "5",
            "10",
            "50",
            "100",
            "1000",
            "10000",
            "100000",
            "1000000",
            "10000000",
            "100000000",
            "1000000000",
        ]

        # Initial x,y movements for referencing are set to zero
        self.x_movement = 0
        self.y_movement = 0

        # Multiplot mode is initially set to off
        self.multiplot_mode = False

        # Dictionary to store the values of the sliders for each spectrum in multiplot mode
        self.values_dictionary = {}

        # Initial multiply factor is 1
        self.multiply_factor = 1

        # 1D slice color of 2D spectra is initially set to navy
        self.slice_colour = "navy"

        # Initial colour/reference/vertical index from list of colours is set to 0
        self.index = 0
        self.ref_index = 0
        self.vertical_index = 0

        # List to hold the multiple 2D spectra in multiplot mode
        self.twoD_spectra = []

        self.linewidth = 1.0
        self.linewidth1D = 1.5

        self.x_difference = 0
        self.y_difference = 0

        # Initially set the transpose flag to False
        self.transpose = False
        self.transposed2D = False
        self.transposed3D = [2, 3, 1]

        # Default options for pivot point for P1 phasing
        self.pivot_x_default = 0
        self.pivot_x = self.pivot_x_default

        self.pivot_y_default = 0
        self.pivot_y = self.pivot_y_default

        self.slice_mode = None

        self.show_bottom_sizer = True

        # Suppress complex warning from numpy
        import warnings

        # warnings.simplefilter("ignore", np.ComplexWarning)  # For old numpy versions
        warnings.simplefilter(
            "ignore", np.exceptions.ComplexWarning
        )  # For new numpy versions

    def UpdateFrame(self):
        self.canvas.draw()
        self.canvas.Refresh()
        self.canvas.Update()
        self.panel.Refresh()
        self.panel.Update()

    def create_button_panel_3D(self):
        self.phasing_label = wx.StaticBox(self, -1, "Phasing:")
        self.phasing_sizer = wx.StaticBoxSizer(self.phasing_label, wx.VERTICAL)

        self.P0_label = wx.StaticText(self, label="P0 (Coarse):")
        self.P1_label = wx.StaticText(self, label="P1 (Coarse):")
        self.P0_slider = FloatSlider(
            self, id=-1, value=0, minval=-180, maxval=180, res=0.1, size=(257, height)
        )
        self.P1_slider = FloatSlider(
            self, id=-1, value=0, minval=-180, maxval=180, res=0.1, size=(257, height)
        )
        self.P0_slider.Bind(wx.EVT_SLIDER, self.OnSliderScroll3D)
        self.P1_slider.Bind(wx.EVT_SLIDER, self.OnSliderScroll3D)

        self.P0_label_fine = wx.StaticText(self, label="P0 (Fine):     ")
        self.P1_label_fine = wx.StaticText(self, label="P1 (Fine):     ")
        self.P0_slider_fine = FloatSlider(
            self, id=-1, value=0, minval=-10, maxval=10, res=0.01, size=(257, height)
        )
        self.P1_slider_fine = FloatSlider(
            self, id=-1, value=0, minval=-10, maxval=10, res=0.01, size=(257, height)
        )
        self.P0_slider_fine.Bind(wx.EVT_SLIDER, self.OnSliderScroll3D)
        self.P1_slider_fine.Bind(wx.EVT_SLIDER, self.OnSliderScroll3D)

        self.sizer_coarse = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer_coarse.Add(self.P0_label)
        self.sizer_coarse.AddSpacer(5)
        self.sizer_coarse.Add(self.P0_slider)
        self.sizer_coarse.AddSpacer(20)
        self.sizer_coarse.Add(self.P1_label)
        self.sizer_coarse.AddSpacer(5)
        self.sizer_coarse.Add(self.P1_slider)

        self.sizer_fine = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer_fine.Add(self.P0_label_fine)
        self.sizer_fine.AddSpacer(5)
        self.sizer_fine.Add(self.P0_slider_fine)
        self.sizer_fine.AddSpacer(20)
        self.sizer_fine.Add(self.P1_label_fine)
        self.sizer_fine.AddSpacer(5)
        self.sizer_fine.Add(self.P1_slider_fine)

        self.phasing_combined = wx.BoxSizer(wx.HORIZONTAL)
        self.P0_total = wx.StaticText(self, label="P0 (Total):")
        self.P1_total = wx.StaticText(self, label="P1 (Total):")
        self.P0_total_value = wx.StaticText(self, label="0")
        self.P1_total_value = wx.StaticText(self, label="0")
        self.phasing_combined.Add(self.P0_total)
        self.phasing_combined.AddSpacer(135)
        self.phasing_combined.Add(self.P0_total_value)
        self.phasing_combined.AddSpacer(152)
        self.phasing_combined.Add(self.P1_total)
        self.phasing_combined.AddSpacer(135)
        self.phasing_combined.Add(self.P1_total_value)

        self.phasing_sizer.AddSpacer(5)
        self.phasing_sizer.Add(self.sizer_coarse)
        self.phasing_sizer.AddSpacer(10)
        self.phasing_sizer.Add(self.sizer_fine)
        self.phasing_sizer.AddSpacer(10)
        self.phasing_sizer.Add(self.phasing_combined)

        # Add sliders to move the 2D plots left/right/up/down with a combobox to choose the scale of the slider
        self.move_label = wx.StaticBox(
            self,
            -1,
            "Move 2D Plot:                                                                              Range(ppm):",
        )
        self.move_sizer = wx.StaticBoxSizer(self.move_label, wx.VERTICAL)
        self.move_x = wx.BoxSizer(wx.HORIZONTAL)
        self.move_y = wx.BoxSizer(wx.HORIZONTAL)
        self.move_x.Add(wx.StaticText(self, label="X:"))
        self.move_y.Add(wx.StaticText(self, label="Y:"))
        self.move_x.AddSpacer(5)
        self.move_y.AddSpacer(5)
        self.move_x_slider = FloatSlider(
            self,
            id=-1,
            value=0,
            minval=-self.reference_rangeX,
            maxval=self.reference_rangeX,
            res=self.reference_rangeX / 1000,
            size=(300, height),
        )
        self.move_y_slider = FloatSlider(
            self,
            id=-1,
            value=0,
            minval=-self.reference_rangeY,
            maxval=self.reference_rangeY,
            res=self.reference_rangeY / 1000,
            size=(300, height),
        )
        self.move_x_slider.Bind(wx.EVT_SLIDER, self.OnMoveX_3D)
        self.move_y_slider.Bind(wx.EVT_SLIDER, self.OnMoveY_3D)
        self.reference_range_chooserX = wx.ComboBox(
            self,
            value=self.reference_range_values[0],
            choices=self.reference_range_values,
        )
        self.reference_range_chooserX.Bind(wx.EVT_COMBOBOX, self.OnReferenceComboX_3D)
        self.reference_range_chooserY = wx.ComboBox(
            self,
            value=self.reference_range_values[0],
            choices=self.reference_range_values,
        )
        self.reference_range_chooserY.Bind(wx.EVT_COMBOBOX, self.OnReferenceComboY_3D)
        self.move_x.Add(self.move_x_slider)
        self.move_x.AddSpacer(5)
        self.move_x.Add(self.reference_range_chooserX)
        self.move_y.Add(self.move_y_slider)
        self.move_y.AddSpacer(5)
        self.move_y.Add(self.reference_range_chooserY)
        self.move_sizer.Add(self.move_x)
        self.move_sizer.AddSpacer(5)
        self.move_sizer.Add(self.move_y)
        self.move_sizer.AddSpacer(5)
        self.move_val_box = wx.BoxSizer(wx.HORIZONTAL)
        self.move_val_box.AddSpacer(20)
        self.move_val_box.Add(wx.StaticText(self, label="Move X (ppm):"))
        self.move_val_box.AddSpacer(5)
        self.move_val_x = wx.StaticText(self, label="0.00")
        self.move_val_box.Add(self.move_val_x)
        self.move_val_box.AddSpacer(35)
        self.move_val_box.Add(wx.StaticText(self, label="Move Y (ppm):"))
        self.move_val_box.AddSpacer(5)
        self.move_val_y = wx.StaticText(self, label="0.00")
        self.move_val_box.Add(self.move_val_y)
        self.move_sizer.Add(self.move_val_box)

        self.bottom_left_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.bottom_left_sizer.Add(self.move_sizer)

        # Create a slider for changing the linewidth of the contour lines
        self.linewidth_label = wx.StaticBox(self, -1, "Contour Line Width:")
        self.linewidth_sizer = wx.StaticBoxSizer(self.linewidth_label, wx.VERTICAL)
        self.linewidth_slider = FloatSlider(
            self, id=-1, value=0.5, minval=0.1, maxval=2, res=0.1, size=(265, height)
        )
        self.linewidth_slider.Bind(wx.EVT_SLIDER, self.OnLinewidthScroll3D)
        self.linewidth_sizer.AddSpacer(5)
        self.linewidth_sizer.Add(self.linewidth_slider)
        self.linewidth_sizer.AddSpacer(5)
        self.contour_linewidth = 0.5

        # Create a sizer to slide through the z axis levels
        self.z_label = wx.StaticBox(
            self, -1, "Z Value (" + str(self.nmrdata.axislabels[0]) + "):"
        )
        z_values = ng.pipe.make_uc(
            self.nmrdata.dic, self.nmrdata.data, dim=0
        ).ppm_scale()
        self.z_slider = FloatSlider(
            self,
            id=-1,
            value=0,
            minval=0,
            maxval=len(z_values) - 1,
            res=1,
            size=(265, height),
        )
        self.z_slider.Bind(wx.EVT_SLIDER, self.OnZScroll3D)
        self.z_sizer = wx.StaticBoxSizer(self.z_label, wx.VERTICAL)
        self.z_sizer.AddSpacer(15)
        self.z_sizer.Add(self.z_slider)
        self.z_sizer.AddSpacer(15)
        self.z_val_box = wx.BoxSizer(wx.HORIZONTAL)
        self.z_val_box.AddSpacer(132)
        self.z_val = wx.StaticText(self, label="0")
        self.z_val_box.Add(self.z_val)
        self.z_sizer.Add(self.z_val_box)
        self.z_sizer.AddSpacer(4)

        self.bottom_left_sizer.AddSpacer(10)
        self.bottom_left_sizer.Add(self.z_sizer)

        self.left_sizer = wx.BoxSizer(wx.VERTICAL)
        self.left_sizer.Add(self.phasing_sizer)
        self.left_sizer.AddSpacer(10)
        self.left_sizer.Add(self.bottom_left_sizer)
        self.left_sizer.AddSpacer(10)

        # Create a sizer for changing the contour levels of the peaks
        self.contour_label = wx.StaticBox(self, -1, "Contour Start = max(data)/x")
        self.contour_sizer = wx.StaticBoxSizer(self.contour_label, wx.VERTICAL)
        self.csizer = wx.BoxSizer(wx.HORIZONTAL)
        self.contour2_label = wx.StaticText(self, label="x:")
        self.contour_slider = FloatSlider(
            self, id=-1, value=1, minval=0, maxval=3, res=0.01, size=(250, height)
        )
        self.contour_slider.Bind(wx.EVT_SLIDER, self.OnMinContour3D)
        self.csizer.Add(self.contour2_label)
        self.csizer.AddSpacer(5)
        self.csizer.Add(self.contour_slider)
        self.contour_sizer.Add(self.csizer)
        self.contour_sizer.AddSpacer(5)
        self.contour_val_box = wx.BoxSizer(wx.HORIZONTAL)
        self.contour_val_box.AddSpacer(125)
        self.contour_val = wx.StaticText(self, label="10")
        self.contour_val_box.Add(self.contour_val)
        self.contour_sizer.Add(self.contour_val_box)

        # Create a sizer for changing the y axis limits of a selected 1D slice in the 2D plot
        self.intensity_label = wx.StaticBox(self, -1, "Intensity Scaling 1D (%):")
        self.intensity_sizer = wx.StaticBoxSizer(self.intensity_label, wx.VERTICAL)
        self.intensity_slider = FloatSlider(
            self, id=-1, value=2, minval=0, maxval=6, res=0.01, size=(265, height)
        )
        self.intensity_slider.Bind(wx.EVT_SLIDER, self.OnIntensityScroll3D)
        self.intensity_sizer.AddSpacer(5)
        self.intensity_sizer.Add(self.intensity_slider)
        self.intensity_sizer.AddSpacer(5)

        # Create a button called Projections, which when clicking it will open up a new window with the projections of the 3D plot
        self.projection_button = wx.Button(self, label="Projections", size=(120, 30))
        self.projection_button.Bind(wx.EVT_BUTTON, self.OnProjectionButton)

        # Create a button called 3D Plot which when clicking it will open up a new window with a full 3D plot
        self.plot3D_button = wx.Button(self, label="3D Plot", size=(120, 30))
        self.plot3D_button.Bind(wx.EVT_BUTTON, self.OnPlot3DButton)

        # This button will create a waterfall plot along the pseudo axis of the currently highlighted slice contour plot (either horizontal or vertical)
        self.waterfall_button = wx.Button(self, label="Waterfall Plot", size=(120, 30))
        self.waterfall_button.Bind(wx.EVT_BUTTON, self.OnWaterfallButton)

        # Create a combobox with the different possible data orientations (e.g. X: H, Y: C13, Z: N15)) which the user can select from
        self.orientation_label = wx.StaticBox(self, -1, "Data Orientation:")
        self.orientation_sizer = wx.StaticBoxSizer(self.orientation_label, wx.VERTICAL)
        labels = self.nmrdata.axislabels
        # Set the initial label to 1,2,0
        options = ["(" + labels[1] + "," + labels[2] + ")," + labels[0]]
        options.append("(" + labels[2] + "," + labels[1] + ")," + labels[0])
        options.append("(" + labels[1] + "," + labels[0] + ")," + labels[2])
        options.append("(" + labels[0] + "," + labels[1] + ")," + labels[2])
        options.append("(" + labels[2] + "," + labels[0] + ")," + labels[1])
        options.append("(" + labels[0] + "," + labels[2] + ")," + labels[1])

        self.orientation_chooser = wx.ComboBox(self, value=options[0], choices=options)
        self.orientation_chooser.Bind(wx.EVT_COMBOBOX, self.OnOrientationCombo)
        self.orientation_chooser.SetSelection(0)
        self.orientation_sizer.Add(self.orientation_chooser)

        # Create a button for changing the labels of the axes
        self.label_button = wx.Button(self, label="Change Labels", size=(120, 30))
        self.label_button.Bind(wx.EVT_BUTTON, self.OnLabelButton3D)

        # Create a button for re-processing
        self.reprocess_button = wx.Button(self, label="Re-Process", size=(120, 30))
        self.reprocess_button.Bind(wx.EVT_BUTTON, self.OnReprocessButton)

        # Create a button to show bore
        self.show_bore_button = wx.Button(self, label="Show Bore", size=(120, 30))
        self.show_bore_button.Bind(wx.EVT_BUTTON, self.OnShowBoreButton)

        # Create a button to show/hide options
        self.show_hide_button = wx.Button(self, label="Hide Options", size=(120, 30))
        self.show_hide_button.Bind(wx.EVT_BUTTON, self.OnHideButton)

        # Put all sizers together
        self.bottom_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.bottom_sizer.Add(self.left_sizer)
        self.bottom_sizer.AddSpacer(5)
        rightbox = wx.BoxSizer(wx.VERTICAL)
        rightbox.Add(self.contour_sizer)
        rightbox.AddSpacer(10)
        rightbox.Add(self.linewidth_sizer)
        rightbox.AddSpacer(10)
        rightbox.Add(self.intensity_sizer)
        right_right_sizer = wx.BoxSizer(wx.VERTICAL)
        right_right_sizer.AddSpacer(5)
        self.button_sizer = wx.BoxSizer(wx.VERTICAL)
        self.button_sizer.Add(self.orientation_sizer)
        self.button_sizer.AddSpacer(5)
        self.button_sizer.Add(self.projection_button)
        self.button_sizer.AddSpacer(5)
        self.button_sizer.Add(self.plot3D_button)
        self.button_sizer.AddSpacer(5)
        self.button_sizer.Add(self.waterfall_button)
        self.button_sizer.AddSpacer(5)
        self.button_sizer.Add(self.label_button)
        self.button_sizer.AddSpacer(5)
        self.button_sizer.Add(self.reprocess_button)
        self.button_sizer.AddSpacer(5)
        self.button_sizer.Add(self.show_bore_button)
        self.button_sizer.AddSpacer(5)
        self.button_sizer.Add(self.show_hide_button)
        right_right_sizer.Add(self.button_sizer)
        self.bottom_sizer.Add(rightbox)
        self.bottom_sizer.AddSpacer(5)
        self.bottom_sizer.Add(right_right_sizer)

    def OnHideButton(self, event):
        if self.show_bottom_sizer == True:
            self.main_sizer.Hide(self.bottom_sizer)
            self.main_sizer.Show(self.show_button_sizer)
            self.UpdateFrame()
            self.Layout()
            self.show_bottom_sizer = False
        else:
            self.main_sizer.Show(self.bottom_sizer)
            self.main_sizer.Hide(self.show_button_sizer)
            self.show_bottom_sizer = True
            self.UpdateFrame()
            self.Layout()

    def OnReprocessButton(self, event):
        # Open an instance of SpinProcess
        if self.parent.path != "":
            os.chdir(self.parent.path)
        from SpinExplorer.SpinProcess import SpinProcess

        reprocessing_frame = SpinProcess(self)
        reprocessing_frame.reprocess = True
        if self.parent.cwd != "":
            os.chdir(self.parent.cwd)

    def OnShowBoreButton(self, event):
        # Open a SpinBore frame

        # Find out which projection is currently selected
        if self.orientation_chooser.GetSelection() == 0:
            # projection is x_name.y_name.dat
            projection = (
                self.nmrdata.axislabels[1] + "." + self.nmrdata.axislabels[2] + ".dat"
            )
        elif self.orientation_chooser.GetSelection() == 1:
            # projection is y_name.x_name.dat
            projection = (
                self.nmrdata.axislabels[2] + "." + self.nmrdata.axislabels[1] + ".dat"
            )
        elif self.orientation_chooser.GetSelection() == 2:
            # projection is x_name.z_name.dat
            projection = (
                self.nmrdata.axislabels[1] + "." + self.nmrdata.axislabels[0] + ".dat"
            )
        elif self.orientation_chooser.GetSelection() == 3:
            # projection is z_name.x_name.dat
            projection = (
                self.nmrdata.axislabels[0] + "." + self.nmrdata.axislabels[1] + ".dat"
            )
        elif self.orientation_chooser.GetSelection() == 4:
            # projection is z_name.y_name.dat
            projection = (
                self.nmrdata.axislabels[0] + "." + self.nmrdata.axislabels[2] + ".dat"
            )
        elif self.orientation_chooser.GetSelection() == 5:
            # projection is y_name.z_name.dat
            projection = (
                self.nmrdata.axislabels[2] + "." + self.nmrdata.axislabels[0] + ".dat"
            )

        # Check to see if the projection file exists
        if os.path.exists(projection) == False:
            # Swap the axis labels
            name = projection.split(".dat")[0].split(".")
            projection = name[1] + "." + name[0] + ".dat"

        # Check to see if the projection file exists
        if os.path.exists(projection) == False:
            # Give a warning that the projection file does not exist
            dlg = wx.MessageDialog(
                self,
                "The projection file does not exist. ",
                "Warning",
                wx.OK | wx.ICON_WARNING,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return

        frame = SpinBore(title="SpinBore", projection=projection, parent=self)
        frame.Show()

    def OnOrientationCombo(self, event):
        self.nmrdata.data = self.data_original
        if self.orientation_chooser.GetSelection() == 0:
            self.z_label.SetLabel("Z Value (" + str(self.nmrdata.axislabels[0]) + "):")
            self.ax.clear()
            self.axes1D.clear()
            self.axes1D_2.clear()
            self.fig.clear()
            self.draw_figure_3D()
            self.ax.set_xlabel(self.nmrdata.axislabels[1])
            self.ax.set_ylabel(self.nmrdata.axislabels[2])
        elif self.orientation_chooser.GetSelection() == 1:
            self.z_label.SetLabel("Z Value (" + str(self.nmrdata.axislabels[0]) + "):")
            self.ax.clear()
            # Get ppm values for x and y axis
            self.uc0 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=2)
            self.uc1 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=1)
            self.uc2 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=0)
            self.ppms_0 = self.uc0.ppm_scale()
            self.ppms_1 = self.uc1.ppm_scale()
            self.ppms_2 = self.uc2.ppm_scale()

            # Transpose the data to the right format
            self.nmrdata.data = np.transpose(self.data_original, (0, 2, 1))

            # Find the plane of the 3D data that has the highest total intensity
            self.total_intensity = []
            for i in range(len(self.nmrdata.data)):
                self.total_intensity.append(np.sum(np.abs(self.nmrdata.data[i])))

            self.max_intensity_index = np.argmax(self.total_intensity)
            # Set the z slider to the index of the plane with the highest total intensity
            self.z_slider.SetMax(len(self.ppms_2) - 1)
            self.z_slider.SetValue(self.max_intensity_index)

            # Replot the data
            self.replot_3D()

            self.ax.set_xlabel(self.nmrdata.axislabels[2])
            self.ax.set_ylabel(self.nmrdata.axislabels[1])
            self.UpdateFrame()

        elif self.orientation_chooser.GetSelection() == 2:
            self.z_label.SetLabel("Z Value (" + str(self.nmrdata.axislabels[2]) + "):")
            self.ax.clear()
            # Get ppm values for x and y axis
            self.uc0 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=1)
            self.uc1 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=0)
            self.uc2 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=2)
            self.ppms_0 = self.uc0.ppm_scale()
            self.ppms_1 = self.uc1.ppm_scale()
            self.ppms_2 = self.uc2.ppm_scale()

            # Transpose the data to the right format
            self.nmrdata.data = np.transpose(self.data_original, (2, 1, 0))
            # Find the plane of the 3D data that has the highest total intensity
            self.total_intensity = []
            for i in range(len(self.nmrdata.data)):
                self.total_intensity.append(np.sum(np.abs(self.nmrdata.data[i])))

            self.max_intensity_index = np.argmax(self.total_intensity)
            # Set the z slider to the index of the plane with the highest total intensity
            self.z_slider.SetMax(len(self.ppms_2) - 1)
            self.z_slider.SetValue(self.max_intensity_index)

            self.replot_3D()
            self.ax.set_xlabel(self.nmrdata.axislabels[1])
            self.ax.set_ylabel(self.nmrdata.axislabels[0])
            self.UpdateFrame()

        elif self.orientation_chooser.GetSelection() == 3:
            self.z_label.SetLabel("Z Value (" + str(self.nmrdata.axislabels[2]) + "):")
            self.ax.clear()
            # Get ppm values for x and y axis
            self.uc0 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=0)
            self.uc1 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=1)
            self.uc2 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=2)
            self.ppms_0 = self.uc0.ppm_scale()
            self.ppms_1 = self.uc1.ppm_scale()
            self.ppms_2 = self.uc2.ppm_scale()

            # Transpose the data to the right format
            self.nmrdata.data = np.transpose(self.data_original, (2, 0, 1))
            # Find the plane of the 3D data that has the highest total intensity
            self.total_intensity = []
            for i in range(len(self.nmrdata.data)):
                self.total_intensity.append(np.sum(np.abs(self.nmrdata.data[i])))

            self.max_intensity_index = np.argmax(self.total_intensity)
            # Set the z slider to the index of the plane with the highest total intensity
            self.z_slider.SetMax(len(self.ppms_2) - 1)
            self.z_slider.SetValue(self.max_intensity_index)

            self.replot_3D()
            self.ax.set_xlabel(self.nmrdata.axislabels[0])
            self.ax.set_ylabel(self.nmrdata.axislabels[1])
            self.UpdateFrame()

        elif self.orientation_chooser.GetSelection() == 4:
            self.z_label.SetLabel("Z Value (" + str(self.nmrdata.axislabels[1]) + "):")
            self.ax.clear()
            # Get ppm values for x and y axis
            self.uc0 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=2)
            self.uc1 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=0)
            self.uc2 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=1)
            self.ppms_0 = self.uc0.ppm_scale()
            self.ppms_1 = self.uc1.ppm_scale()
            self.ppms_2 = self.uc2.ppm_scale()

            # Transpose the data to the right format
            self.nmrdata.data = np.transpose(self.data_original, (1, 2, 0))
            # Find the plane of the 3D data that has the highest total intensity
            self.total_intensity = []
            for i in range(len(self.nmrdata.data)):
                self.total_intensity.append(np.sum(np.abs(self.nmrdata.data[i])))

            self.max_intensity_index = np.argmax(self.total_intensity)
            # Set the z slider to the index of the plane with the highest total intensity
            self.z_slider.SetMax(len(self.ppms_2) - 1)
            self.z_slider.SetValue(self.max_intensity_index)

            self.replot_3D()
            self.ax.set_xlabel(self.nmrdata.axislabels[2])
            self.ax.set_ylabel(self.nmrdata.axislabels[0])
            self.UpdateFrame()
        elif self.orientation_chooser.GetSelection() == 5:
            self.z_label.SetLabel("Z Value (" + str(self.nmrdata.axislabels[1]) + "):")
            self.ax.clear()
            # Get ppm values for x and y axis
            self.uc0 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=0)
            self.uc1 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=2)
            self.uc2 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=1)
            self.ppms_0 = self.uc0.ppm_scale()
            self.ppms_1 = self.uc1.ppm_scale()
            self.ppms_2 = self.uc2.ppm_scale()

            # Transpose the data to the right format
            self.nmrdata.data = np.transpose(self.data_original, (1, 0, 2))
            # Find the plane of the 3D data that has the highest total intensity
            self.total_intensity = []
            for i in range(len(self.nmrdata.data)):
                self.total_intensity.append(np.sum(np.abs(self.nmrdata.data[i])))

            self.max_intensity_index = np.argmax(self.total_intensity)
            # Set the z slider to the index of the plane with the highest total intensity
            self.z_slider.SetMax(len(self.ppms_2) - 1)
            self.z_slider.SetValue(self.max_intensity_index)

            self.replot_3D()
            self.ax.set_xlabel(self.nmrdata.axislabels[0])
            self.ax.set_ylabel(self.nmrdata.axislabels[2])
            self.UpdateFrame()

    def replot_3D(self):
        self.new_x_ppms = self.ppms_0
        self.new_y_ppms = self.ppms_1
        self.X, self.Y = np.meshgrid(self.ppms_1, self.ppms_0)
        self.ax.contour(
            self.Y,
            self.X,
            self.nmrdata.data[self.max_intensity_index],
            self.cl,
            colors=self.cmap,
            linewidths=self.contour_linewidth,
        )
        self.ax.contour(
            self.Y,
            self.X,
            self.nmrdata.data[self.max_intensity_index],
            self.cl_neg,
            colors=self.cmap_neg,
            linewidths=self.contour_linewidth,
        )
        self.ax.set_xlim(max(self.ppms_0), min(self.ppms_0))
        self.ax.set_ylim(max(self.ppms_1), min(self.ppms_1))
        (self.line1,) = self.axes1D.plot(
            self.ppms_0,
            self.nmrdata.data[self.max_intensity_index][:, 1],
            color=self.slice_colour,
        )
        self.line2 = self.ax.axhline(self.ppms_1[1], color="k")
        self.axes1D.set_ylim(
            -np.max(self.nmrdata.data[self.max_intensity_index] / 10),
            np.max(self.nmrdata.data[self.max_intensity_index]),
        )
        self.axes1D.set_yticks([])
        self.axes1D.set_xticks([])
        self.line1.set_visible(False)
        self.line2.set_visible(False)
        (self.line3,) = self.axes1D_2.plot(
            self.nmrdata.data[self.max_intensity_index][1, :],
            self.ppms_1,
            color=self.slice_colour,
        )
        self.line4 = self.ax.axvline(self.ppms_0[1], color="k")
        self.axes1D_2.set_xlim(
            -np.max(self.nmrdata.data[self.max_intensity_index] / 10),
            np.max(self.nmrdata.data[self.max_intensity_index]),
        )
        self.axes1D_2.set_xticks([])
        self.axes1D_2.set_yticks([])
        self.line3.set_visible(False)
        self.line4.set_visible(False)

        self.UpdateFrame()
        self.OnZScroll3D(None)

    def OnLabelButton3D(self, event):
        # Get the current labels of the x and y axes
        x_label = self.ax.get_xlabel()
        y_label = self.ax.get_ylabel()
        # Z label is the element of self.nmrdata.axislabels that is not x or y
        z_label = self.nmrdata.axislabels[0]
        if z_label == x_label or z_label == y_label:
            z_label = self.nmrdata.axislabels[1]
            if z_label == x_label or z_label == y_label:
                z_label = self.nmrdata.axislabels[2]

        # Get the ppm values for the x and y axes
        x_ppms = self.ppms_0
        y_ppms = self.ppms_1
        z_ppms = self.ppms_2

        # Create a window to allow the user to see the current labels and ppm values and change the labels accordingly
        self.dlg = wx.Dialog(self, title="Change Labels")
        self.dlg.SetSize(500, 200)

        # Create a sizer to hold the labels and ppm values
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.AddSpacer(10)

        self.label_change_label = wx.StaticBox(self.dlg, -1, "Input desired labels:")
        self.total_label_change_sizer = wx.StaticBoxSizer(
            self.label_change_label, wx.VERTICAL
        )
        self.total_label_change_sizer.AddSpacer(10)

        # Create a sizer to hold the x axis labels and ppm values
        x_sizer = wx.BoxSizer(wx.HORIZONTAL)
        x_sizer.AddSpacer(10)
        x_sizer.Add(wx.StaticText(self.dlg, label="X Axis Label:"))
        x_sizer.AddSpacer(5)
        self.xlabel_box = wx.TextCtrl(self.dlg, value=x_label, size=(100, 20))
        x_sizer.Add(self.xlabel_box)
        x_sizer.AddSpacer(10)
        x_ppm_limits = "{:.2f}".format(min(x_ppms)) + "-{:.2f}".format(max(x_ppms))
        x_sizer.Add(wx.StaticText(self.dlg, label="X Axis Limits (ppm):"))
        x_sizer.AddSpacer(5)
        x_sizer.Add(wx.StaticText(self.dlg, label=x_ppm_limits))
        self.total_label_change_sizer.Add(x_sizer)
        self.total_label_change_sizer.AddSpacer(10)

        # Create a sizer to hold the y axis labels and ppm values
        y_sizer = wx.BoxSizer(wx.HORIZONTAL)
        y_sizer.AddSpacer(10)
        y_sizer.Add(wx.StaticText(self.dlg, label="Y Axis Label:"))
        y_sizer.AddSpacer(5)
        self.ylabel_box = wx.TextCtrl(self.dlg, value=y_label, size=(100, 20))
        y_sizer.Add(self.ylabel_box)
        y_sizer.AddSpacer(10)
        y_ppm_limits = "{:.2f}".format(min(y_ppms)) + "-{:.2f}".format(max(y_ppms))
        y_sizer.Add(wx.StaticText(self.dlg, label="Y Axis Limits (ppm):"))
        y_sizer.AddSpacer(5)
        y_sizer.Add(wx.StaticText(self.dlg, label=y_ppm_limits))
        self.total_label_change_sizer.Add(y_sizer)
        self.total_label_change_sizer.AddSpacer(10)

        # Create a sizer to hold the y axis labels and ppm values
        z_sizer = wx.BoxSizer(wx.HORIZONTAL)
        z_sizer.AddSpacer(10)
        z_sizer.Add(wx.StaticText(self.dlg, label="Z Axis Label:"))
        z_sizer.AddSpacer(5)
        self.zlabel_box = wx.TextCtrl(self.dlg, value=z_label, size=(100, 20))
        z_sizer.Add(self.zlabel_box)
        z_sizer.AddSpacer(10)
        z_ppm_limits = "{:.2f}".format(min(z_ppms)) + "-{:.2f}".format(max(z_ppms))
        z_sizer.Add(wx.StaticText(self.dlg, label="Z Axis Limits (ppm):"))
        z_sizer.AddSpacer(5)
        z_sizer.Add(wx.StaticText(self.dlg, label=z_ppm_limits))
        self.total_label_change_sizer.Add(z_sizer)
        self.total_label_change_sizer.AddSpacer(10)

        # Add a save and close button to the sizer
        save_button = wx.Button(self.dlg, label="Save")
        save_button.Bind(wx.EVT_BUTTON, self.OnSaveLabels3D)
        self.total_label_change_sizer.Add(save_button, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.total_label_change_sizer.AddSpacer(10)

        sizer.Add(self.total_label_change_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)

        # Show the sizer in the dialog
        self.dlg.SetSizer(sizer)

        # Show the dialog
        self.dlg.ShowModal()
        self.dlg.Destroy()

    def OnSaveLabels3D(self, event):
        # Get the new labels for the x and y axes
        x_label = self.xlabel_box.GetValue()
        y_label = self.ylabel_box.GetValue()
        z_label = self.zlabel_box.GetValue()

        orientation_chooser_selection = self.orientation_chooser.GetSelection()

        # Update the labels in the plot, z_slider label and data orientation options
        if self.orientation_chooser.GetSelection() == 0:
            self.nmrdata.axislabels = [z_label, x_label, y_label]
            self.ax.set_xlabel(self.nmrdata.axislabels[1])
            self.ax.set_ylabel(self.nmrdata.axislabels[2])
            self.z_label.SetLabel("Z Value (" + str(self.nmrdata.axislabels[0]) + "):")
        elif self.orientation_chooser.GetSelection() == 1:
            self.nmrdata.axislabels = [z_label, y_label, x_label]
            self.ax.set_xlabel(self.nmrdata.axislabels[2])
            self.ax.set_ylabel(self.nmrdata.axislabels[1])
            self.z_label.SetLabel("Z Value (" + str(self.nmrdata.axislabels[0]) + "):")
        elif self.orientation_chooser.GetSelection() == 2:
            self.nmrdata.axislabels = [y_label, x_label, z_label]
            self.ax.set_xlabel(self.nmrdata.axislabels[1])
            self.ax.set_ylabel(self.nmrdata.axislabels[0])
            self.z_label.SetLabel("Z Value (" + str(self.nmrdata.axislabels[2]) + "):")
        elif self.orientation_chooser.GetSelection() == 3:
            self.nmrdata.axislabels = [x_label, y_label, z_label]
            self.ax.set_xlabel(self.nmrdata.axislabels[0])
            self.ax.set_ylabel(self.nmrdata.axislabels[1])
            self.z_label.SetLabel("Z Value (" + str(self.nmrdata.axislabels[2]) + "):")
        elif self.orientation_chooser.GetSelection() == 4:
            self.nmrdata.axislabels = [y_label, z_label, x_label]
            self.ax.set_xlabel(self.nmrdata.axislabels[2])
            self.ax.set_ylabel(self.nmrdata.axislabels[0])
            self.z_label.SetLabel("Z Value (" + str(self.nmrdata.axislabels[1]) + "):")
        elif self.orientation_chooser.GetSelection() == 5:
            self.nmrdata.axislabels = [x_label, z_label, y_label]
            self.ax.set_xlabel(self.nmrdata.axislabels[0])
            self.ax.set_ylabel(self.nmrdata.axislabels[2])
            self.z_label.SetLabel("Z Value (" + str(self.nmrdata.axislabels[1]) + "):")

        # Save the labels to a labels.txt file
        if self.parent.path != "":
            os.chdir(self.parent.path)
        with open("labels.txt", "w") as f:
            # write labels as label1, label2
            f.write(
                self.nmrdata.axislabels[0]
                + ","
                + self.nmrdata.axislabels[1]
                + ","
                + self.nmrdata.axislabels[2]
            )
        if self.parent.cwd != "":
            os.chdir(self.parent.cwd)

        self.orientation_chooser.Clear()
        self.orientation_chooser.Append(
            [
                "("
                + self.nmrdata.axislabels[1]
                + ","
                + self.nmrdata.axislabels[2]
                + "),"
                + self.nmrdata.axislabels[0],
                "("
                + self.nmrdata.axislabels[2]
                + ","
                + self.nmrdata.axislabels[1]
                + "),"
                + self.nmrdata.axislabels[0],
                "("
                + self.nmrdata.axislabels[1]
                + ","
                + self.nmrdata.axislabels[0]
                + "),"
                + self.nmrdata.axislabels[2],
                "("
                + self.nmrdata.axislabels[0]
                + ","
                + self.nmrdata.axislabels[1]
                + "),"
                + self.nmrdata.axislabels[2],
                "("
                + self.nmrdata.axislabels[2]
                + ","
                + self.nmrdata.axislabels[0]
                + "),"
                + self.nmrdata.axislabels[1],
                "("
                + self.nmrdata.axislabels[0]
                + ","
                + self.nmrdata.axislabels[2]
                + "),"
                + self.nmrdata.axislabels[1],
            ]
        )
        self.orientation_chooser.SetSelection(orientation_chooser_selection)
        self.OnOrientationCombo(None)

    def draw_figure_3D(self):
        self.ax = self.fig.add_subplot(111)
        self.axes1D = self.ax.twinx()
        self.axes1D_2 = self.ax.twiny()

        self.fig.canvas.mpl_connect("key_press_event", self.on_key_3d)
        self.fig.canvas.mpl_connect("button_press_event", self.on_click_3d)

        # plot parameters
        contour_start = (
            np.max(np.abs(self.nmrdata.data)) / 10
        )  # contour level start value
        self.contour_num = 20  # number of contour levels
        self.contour_factor = 1.2  # scaling factor between contour levels
        # calculate contour levels
        self.cl = contour_start * self.contour_factor ** np.arange(self.contour_num)
        self.cl_neg = -contour_start * self.contour_factor ** np.flip(
            np.arange(self.contour_num)
        )

        # Find the plane of the 3D data that has the highest total intensity
        self.total_intensity = []
        for i in range(len(self.nmrdata.data)):
            self.total_intensity.append(np.sum(np.abs(self.nmrdata.data[i])))

        self.max_intensity_index = np.argmax(self.total_intensity)
        # Set the z slider to the index of the plane with the highest total intensity
        self.z_slider.SetValue(self.max_intensity_index)

        # Get ppm values for x and y axis
        self.data_original = self.nmrdata.data

        if self.nmrdata.file != ".":
            self.uc0 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=1)
            self.uc1 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=2)
            self.uc2 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=0)
        else:
            udic = ng.bruker.guess_udic(self.nmrdata.dic, self.nmrdata.data)
            self.uc0 = ng.fileiobase.uc_from_udic(udic, dim=1)
            self.uc1 = ng.fileiobase.uc_from_udic(udic, dim=2)
            self.uc2 = ng.fileiobase.uc_from_udic(udic, dim=0)

        self.ppms_0 = self.uc0.ppm_scale()
        self.ppms_1 = self.uc1.ppm_scale()
        self.ppms_2 = self.uc2.ppm_scale()
        self.new_x_ppms = self.ppms_0
        self.new_y_ppms = self.ppms_1
        self.X, self.Y = np.meshgrid(self.ppms_1, self.ppms_0)
        self.ax.contour(
            self.Y,
            self.X,
            self.nmrdata.data[self.max_intensity_index],
            self.cl,
            colors=self.cmap,
            linewidths=self.contour_linewidth,
        )
        self.ax.contour(
            self.Y,
            self.X,
            self.nmrdata.data[self.max_intensity_index],
            self.cl_neg,
            colors=self.cmap_neg,
            linewidths=self.contour_linewidth,
        )
        self.ax.set_xlabel(self.nmrdata.axislabels[1])
        self.ax.set_ylabel(self.nmrdata.axislabels[2])
        self.ax.set_xlim(max(self.ppms_0), min(self.ppms_0))
        self.ax.set_ylim(max(self.ppms_1), min(self.ppms_1))

        (self.line1,) = self.axes1D.plot(
            self.ppms_0,
            self.nmrdata.data[self.max_intensity_index][:, 1],
            color=self.slice_colour,
        )
        self.line2 = self.ax.axhline(self.ppms_1[1], color="k")
        self.axes1D.set_ylim(
            -np.max(self.nmrdata.data[self.max_intensity_index] / 10),
            np.max(self.nmrdata.data[self.max_intensity_index]),
        )
        self.axes1D.set_yticks([])
        self.axes1D.set_xticks([])
        self.line1.set_visible(False)
        self.line2.set_visible(False)
        (self.line3,) = self.axes1D_2.plot(
            self.nmrdata.data[self.max_intensity_index][1, :],
            self.ppms_1,
            color=self.slice_colour,
        )
        self.line4 = self.ax.axvline(self.ppms_0[1], color="k")
        self.axes1D_2.set_xlim(
            -np.max(self.nmrdata.data[self.max_intensity_index] / 10),
            np.max(self.nmrdata.data[self.max_intensity_index]),
        )
        self.axes1D_2.set_xticks([])
        self.axes1D_2.set_yticks([])
        self.line3.set_visible(False)
        self.line4.set_visible(False)

        self.z_slider.SetMax(len(self.ppms_2) - 1)

        self.UpdateFrame()
        self.OnZScroll3D(None)

    def OnProjectionButton(self, event):
        # Make projection window
        self.projection_window = ProjectionFrame(parent=self, title="Projections")

    def OnPlot3DButton(self, event):
        # Make a 3D plot window
        self.threeD_warning = wx.MessageDialog(
            self,
            "3D plotting can take a while, do you want to continue?",
            "Warning",
            wx.YES_NO | wx.ICON_WARNING,
        )
        self.threeD_warning.ShowModal()
        if self.threeD_warning.ShowModal() == wx.ID_YES:
            self.plot3D_window = Plot3DFrame(parent=self, title="3D Plot")
            self.threeD_warning.Destroy()
        else:
            self.threeD_warning.Destroy()
            return

    def OnWaterfallButton(self, event):
        # See if the user has selected a slice
        if self.line1.get_visible() == False and self.line3.get_visible() == False:
            # Give a warning that the user needs to select a slice
            dlg = wx.MessageDialog(
                self,
                "Please select a slice to produce a waterfall plot.",
                "Warning",
                wx.OK | wx.ICON_WARNING,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return

        # Give a text message to say this will produce a waterfall plot along the pseudo axis of the slice highlighted in the contour plot
        self.waterfall_warning = wx.MessageDialog(
            self,
            "This will produce a waterfall plot along the pseudo axis of the slice highlighted in the contour plot. Do you want to continue?",
            "Warning",
            wx.YES_NO | wx.ICON_WARNING,
        )
        self.waterfall_warning.ShowModal()
        if self.waterfall_warning.ShowModal() == wx.ID_YES:
            if self.line1.get_visible() == True:
                visible = "line1"
            else:
                visible = "line3"
            self.waterfall_window = WaterfallFrame(
                parent=self, title="Waterfall Plot", visible=visible
            )
            self.waterfall_warning.Destroy()
        else:
            self.waterfall_warning.Destroy()
            return

    def OnLinewidthScroll3D(self, event):
        self.contour_linewidth = float(self.linewidth_slider.GetValue())
        self.OnMinContour3D(event)

    def OnMoveX_3D(self, event):
        # update x-axis
        z_index = int(self.z_slider.GetValue())
        self.x_movement = float(self.move_x_slider.GetValue())
        self.move_val_x.SetLabel(str(round(self.x_movement, 4)))
        self.new_x_ppms = self.ppms_0 + np.ones(len(self.ppms_0)) * self.x_movement
        self.X, self.Y = np.meshgrid(self.new_y_ppms, self.new_x_ppms)
        xlim, ylim = self.ax.get_xlim(), self.ax.get_ylim()
        self.ax.clear()
        self.ax.contour(
            self.Y,
            self.X,
            self.nmrdata.data[z_index],
            self.cl,
            colors=self.cmap,
            linewidths=self.contour_linewidth,
        )
        self.ax.contour(
            self.Y,
            self.X,
            self.nmrdata.data[z_index],
            self.cl_neg,
            colors=self.cmap_neg,
            linewidths=self.contour_linewidth,
        )
        self.ax.set_xlabel(self.nmrdata.axislabels[1])
        self.ax.set_ylabel(self.nmrdata.axislabels[2])
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        if self.line1.get_visible() == True:
            self.line1.set_ydata(
                self.nmrdata.data[z_index][:, self.uc1(str(self.y1) + "ppm")]
            )
            self.line1.set_xdata(self.new_x_ppms)
            self.line2 = self.ax.axhline(self.y1 + self.y_movement, color="k")
            self.axes1D.set_ylim(
                -np.max(self.nmrdata.data[z_index] / 10),
                np.max(self.nmrdata.data[z_index]),
            )
        if self.line3.get_visible() == True:
            self.line3.set_xdata(
                self.nmrdata.data[z_index][self.uc0(str(self.x1) + "ppm"), :]
            )
            self.line3.set_ydata(self.new_y_ppms)
            self.line4 = self.ax.axvline(self.x1 + self.x_movement, color="k")
            self.axes1D_2.set_xlim(
                -np.max(self.nmrdata.data[z_index] / 10),
                np.max(self.nmrdata.data[z_index]),
            )
        self.UpdateFrame()

    def OnMoveY_3D(self, event):
        # update y-axis
        z_index = int(self.z_slider.GetValue())
        self.y_movement = float(self.move_y_slider.GetValue())
        self.move_val_y.SetLabel(str(round(self.y_movement, 4)))
        self.new_y_ppms = self.ppms_1 + np.ones(len(self.ppms_1)) * self.y_movement
        self.X, self.Y = np.meshgrid(self.new_y_ppms, self.new_x_ppms)
        xlim, ylim = self.ax.get_xlim(), self.ax.get_ylim()
        self.ax.clear()
        self.ax.contour(
            self.Y,
            self.X,
            self.nmrdata.data[z_index],
            self.cl,
            colors=self.cmap,
            linestyles="solid",
            linewidths=self.contour_linewidth,
        )
        self.ax.contour(
            self.Y,
            self.X,
            self.nmrdata.data[z_index],
            self.cl_neg,
            colors=self.cmap_neg,
            linestyles="solid",
            linewidths=self.contour_linewidth,
        )
        self.ax.set_xlabel(self.nmrdata.axislabels[1])
        self.ax.set_ylabel(self.nmrdata.axislabels[2])
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        if self.line1.get_visible() == True:
            self.line1.set_ydata(
                self.nmrdata.data[z_index][:, self.uc1(str(self.y1) + "ppm")]
            )
            self.line1.set_xdata(self.new_x_ppms)
            self.line2 = self.ax.axhline(self.y1, color="k")
            self.axes1D.set_ylim(
                -np.max(self.nmrdata.data[z_index] / 10),
                np.max(self.nmrdata.data[z_index]),
            )
        if self.line3.get_visible() == True:
            self.line3.set_xdata(
                self.nmrdata.data[z_index][self.uc0(str(self.x1) + "ppm"), :]
            )
            self.line3.set_ydata(self.new_y_ppms)
            self.line4 = self.ax.axvline(self.x1 + self.x_movement, color="k")
            self.axes1D_2.set_xlim(
                -np.max(self.nmrdata.data[z_index] / 10),
                np.max(self.nmrdata.data[z_index]),
            )
        self.UpdateFrame()

    def OnReferenceComboX_3D(self, event):
        # Change the range for the move-x slider
        index = int(self.reference_range_chooserX.GetSelection())
        self.reference_rangeX = float(self.reference_range_values[index])
        self.move_x_slider.SetMin(-self.reference_rangeX)
        self.move_x_slider.SetMax(self.reference_rangeX)
        self.move_x_slider.SetRes(self.reference_rangeX / 1000)
        self.move_x_slider.Bind(wx.EVT_SLIDER, self.OnMoveX_3D)

    def OnReferenceComboY_3D(self, event):
        # Change the range for the move-y slider
        index = int(self.reference_range_chooserY.GetSelection())
        self.reference_rangeY = float(self.reference_range_values[index])
        self.move_y_slider.SetMin(-self.reference_rangeY)
        self.move_y_slider.SetMax(self.reference_rangeY)
        self.move_y_slider.SetRes(self.reference_rangeY / 1000)
        self.move_y_slider.Bind(wx.EVT_SLIDER, self.OnMoveY_3D)

    def on_key_3d(self, event):
        # Plot horizontal/vertical slices of the data
        if event.key == "h":
            z_index = int(self.z_slider.GetValue())
            self.axes1D.set_ylim(
                -np.max(self.nmrdata.data[z_index] / 8),
                np.max(self.nmrdata.data[z_index]),
            )
            # plot a horizontal slice of the data
            if self.line1.get_visible() == True:
                self.line1.set_visible(False)
                self.line2.set_visible(False)
                self.UpdateFrame()
            else:
                if self.line3.get_visible() == True:
                    self.line3.set_visible(False)
                    self.line4.set_visible(False)
                    self.UpdateFrame()
                else:
                    (self.line1,) = self.axes1D.plot(
                        self.ppms_0,
                        self.nmrdata.data[z_index][
                            :, self.uc1(str(self.ppms_1[1]) + "ppm")
                        ],
                        color=self.slice_colour,
                    )
                    self.line2 = self.ax.axhline(self.ppms_1[1], color="k")
                    self.UpdateFrame()

        if event.key == "v":
            z_index = int(self.z_slider.GetValue())
            self.axes1D_2.set_xlim(
                -np.max(self.nmrdata.data[z_index] / 8),
                np.max(self.nmrdata.data[z_index]),
            )
            if self.line3.get_visible() == True:
                self.line3.set_visible(False)
                self.line4.set_visible(False)
                self.UpdateFrame()
            else:
                if self.line1.get_visible() == True:
                    self.line1.set_visible(False)
                    self.line2.set_visible(False)
                    self.UpdateFrame()
                else:
                    self.line3.set_visible = True
                    self.line4.set_visible = True
                    (self.line3,) = self.axes1D_2.plot(
                        self.nmrdata.data[z_index][
                            self.uc0(str(self.ppms_0[1]) + "ppm"), :
                        ],
                        self.ppms_1,
                        color=self.slice_colour,
                    )
                    self.line4 = self.ax.axvline(self.ppms_0[1], color="k")
                    self.UpdateFrame()

    def on_click_3d(self, event):
        # Get the x and y values of the click and plot the horizontal/vertical slices at that point
        z_index = int(self.z_slider.GetValue())
        self.x1, self.y1 = self.ax.transData.inverted().transform((event.x, event.y))
        if self.x1 != None and self.y1 != None:
            if self.line1.get_visible() == True:
                self.line1.set_ydata(
                    self.nmrdata.data[z_index][
                        :, self.uc1(str(self.y1 - self.y_movement) + "ppm")
                    ]
                )
                self.line2.set_ydata([self.y1])
                self.line1.set_xdata(self.new_x_ppms)
                self.OnSliderScroll3D(None)
            if self.line3.get_visible() == True:
                self.line3.set_xdata(
                    self.nmrdata.data[z_index][
                        self.uc0(str(self.x1 - self.x_movement) + "ppm"), :
                    ]
                )
                self.line4.set_xdata([self.x1])
                self.line3.set_ydata(self.new_y_ppms)
                self.OnSliderScroll3D(None)

    def OnSliderScroll3D(self, event):
        # Get all the slider values for P0 and P1 (coarse and fine), put the combined coarse and fine values on the screen
        self.total_P0 = self.P0_slider.GetValue() + self.P0_slider_fine.GetValue()
        self.total_P1 = self.P1_slider.GetValue() + self.P1_slider_fine.GetValue()
        self.P0_total_value.SetLabel("{:.2f}".format(self.total_P0))
        self.P1_total_value.SetLabel("{:.2f}".format(self.total_P1))
        self.phase3D()

    def phase3D(self):
        z_index = int(self.z_slider.GetValue())
        if self.line1.get_visible() == True:
            data = self.nmrdata.data[z_index][:, self.uc1(str(self.y1) + "ppm")]
            complex_data = ng.process.proc_base.ht(data, self.nmrdata.data.shape[1])
            self.phased_data = ng.process.proc_base.ps(
                complex_data, p0=self.total_P0, p1=self.total_P1
            )
            self.line1.set_ydata(self.phased_data)
        if self.line3.get_visible() == True:
            data = self.nmrdata.data[z_index][self.uc0(str(self.x1) + "ppm"), :]
            complex_data = ng.process.proc_base.ht(data, self.nmrdata.data.shape[2])
            self.phased_data2 = ng.process.proc_base.ps(
                complex_data, p0=self.total_P0, p1=self.total_P1
            )
            self.line3.set_xdata(self.phased_data2)
        self.UpdateFrame()

    def OnMinContour3D(self, event):
        # Get the new contour limits and redraw the plot
        z_index = int(self.z_slider.GetValue())
        contour_val = 10 ** float(self.contour_slider.GetValue())
        self.contour_val.SetLabel(str(int(contour_val)))
        self.contour_start = (
            np.max(np.abs(self.nmrdata.data[int(z_index)])) / contour_val
        )
        self.cl = self.contour_start * self.contour_factor ** np.arange(
            self.contour_num
        )
        self.cl_neg = -self.contour_start * self.contour_factor ** np.flip(
            np.arange(self.contour_num)
        )
        xlim, ylim = self.ax.get_xlim(), self.ax.get_ylim()
        xlabel = self.ax.get_xlabel()
        ylabel = self.ax.get_ylabel()
        self.ax.clear()
        self.ax.contour(
            self.Y,
            self.X,
            self.nmrdata.data[int(z_index)],
            self.cl,
            colors=self.cmap,
            linewidths=self.contour_linewidth,
            linestyles="solid",
        )
        self.ax.contour(
            self.Y,
            self.X,
            self.nmrdata.data[int(z_index)],
            self.cl_neg,
            colors=self.cmap_neg,
            linewidths=self.contour_linewidth,
            linestyles="solid",
        )
        if self.line1.get_visible() == True:
            self.line1.set_ydata(
                self.nmrdata.data[int(z_index)][:, self.uc1(str(self.y1) + "ppm")]
            )
            self.line2 = self.ax.axhline(self.y1 + self.y_movement, color="k")
            self.axes1D.set_ylim(
                -np.max(self.nmrdata.data[int(z_index)] / 10),
                np.max(self.nmrdata.data[int(z_index)]),
            )
        if self.line3.get_visible() == True:
            self.line3.set_xdata(
                self.nmrdata.data[int(z_index)][self.uc0(str(self.x1) + "ppm"), :]
            )
            self.line4 = self.ax.axvline(self.x1 + self.x_movement, color="k")
            self.axes1D_2.set_xlim(
                -np.max(self.nmrdata.data[int(z_index)] / 10),
                np.max(self.nmrdata.data[int(z_index)]),
            )
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        self.UpdateFrame()

    def OnZScroll3D(self, event):
        # Get the new z value and redraw the plot
        z_index = int(self.z_slider.GetValue())
        self.z_val.SetLabel(
            "Index: "
            + str(z_index)
            + " , "
            + "{:.2f}".format(self.ppms_2[z_index - 1])
            + "ppm"
        )
        xlim, ylim = self.ax.get_xlim(), self.ax.get_ylim()
        xlim_1, ylim_1 = self.axes1D.get_xlim(), self.axes1D.get_ylim()
        xlim_2, ylim_2 = self.axes1D_2.get_xlim(), self.axes1D_2.get_ylim()
        xlabel = self.ax.get_xlabel()
        ylabel = self.ax.get_ylabel()
        self.ax.clear()
        self.ax.contour(
            self.Y,
            self.X,
            self.nmrdata.data[int(z_index)],
            self.cl,
            colors=self.cmap,
            linewidths=self.contour_linewidth,
            linestyles="solid",
        )
        self.ax.contour(
            self.Y,
            self.X,
            self.nmrdata.data[int(z_index)],
            self.cl_neg,
            colors=self.cmap_neg,
            linewidths=self.contour_linewidth,
            linestyles="solid",
        )
        if self.line1.get_visible() == True:
            self.line1.set_ydata(
                self.nmrdata.data[int(z_index)][:, self.uc1(str(self.y1) + "ppm")]
            )
            self.line1.set_xdata(self.new_x_ppms)
            self.line2 = self.ax.axhline(self.y1 + self.y_movement, color="k")
            self.axes1D.set_ylim(
                -np.max(self.nmrdata.data[int(z_index)] / 10),
                np.max(self.nmrdata.data[int(z_index)]),
            )
        if self.line3.get_visible() == True:
            self.line3.set_xdata(
                self.nmrdata.data[int(z_index)][self.uc0(str(self.x1) + "ppm"), :]
            )
            self.line3.set_ydata(self.new_y_ppms)
            self.line4 = self.ax.axvline(self.x1 + self.x_movement, color="k")
            self.axes1D_2.set_xlim(
                -np.max(self.nmrdata.data[int(z_index)] / 10),
                np.max(self.nmrdata.data[int(z_index)]),
            )

        # self.axes1D.set_ylim(-np.max(self.nmrdata.data[int(z_index)]/10), np.max(self.nmrdata.data[int(z_index)]))
        # self.axes1D_2.set_xlim(-np.max(self.nmrdata.data[int(z_index)]/10), np.max(self.nmrdata.data[int(z_index)]))
        self.axes1D.set_ylim(ylim_1)
        self.axes1D_2.set_xlim(xlim_2)
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        self.OnSliderScroll3D(event)

    def OnIntensityScroll3D(self, event):
        # Get the new y axis limits of the 1D slice and redraw the plot
        intensity_percent = 10 ** float(self.intensity_slider.GetValue())
        z_index = int(self.z_slider.GetValue())
        if self.line1.get_visible() == True:
            self.axes1D.set_ylim(
                -(np.max(self.nmrdata.data[z_index]) / 8) / (intensity_percent / 100),
                np.max(self.nmrdata.data[z_index]) / (intensity_percent / 100),
            )
            self.UpdateFrame()
        if self.line3.get_visible() == True:
            self.axes1D_2.set_xlim(
                -(np.max(self.nmrdata.data[z_index]) / 8) / (intensity_percent / 100),
                np.max(self.nmrdata.data[z_index]) / (intensity_percent / 100),
            )
            self.UpdateFrame()


# A class which will overlay pseudo2D stacks on a OneDPlot
class StackOverlay:
    def __init__(self, parent, nmrdata, axis):
        self.axis = axis
        self.ucs = []
        self.data = []
        self.parent = parent
        self.parent.extra_plots = (
            []
        )  # This is a list of extra plot objects that are added to the canvas

        if self.parent.nmrdata.dim > 1:
            self.axes1D_H = self.parent.axes1D
            self.axes1D_V = self.parent.axes1D_2
        self.color_list = [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
        ]
        self.custom_lines = []
        self.custom_labels = []
        self.parent.active_plot_index = 0

    def plot_initial(self):
        # Input current values of the sliders for the first plot into the dictionary
        self.parent.values_dictionary[0] = {}
        self.parent.values_dictionary[0]["title"] = "1"
        self.parent.values_dictionary[0][
            "linewidth"
        ] = self.parent.line1.get_linewidth()
        self.parent.values_dictionary[0]["color index"] = self.parent.index
        self.parent.values_dictionary[0]["original_ppms"] = self.parent.ppm_original
        self.parent.values_dictionary[0]["original_data"] = self.parent.nmrdata.data
        self.parent.values_dictionary[0]["move left/right"] = float(
            self.parent.reference_slider.GetValue()
        )
        self.parent.values_dictionary[0][
            "move left/right range index"
        ] = self.parent.ref_index
        self.parent.values_dictionary[0]["move up/down"] = float(
            self.parent.vertical_slider.GetValue()
        )
        self.parent.values_dictionary[0][
            "move up/down range index"
        ] = self.parent.vertical_index
        self.parent.values_dictionary[0]["multiply value"] = float(
            self.parent.multiply_slider.GetValue()
        )
        self.parent.values_dictionary[0]["multiply range index"] = int(
            self.parent.multiply_range_chooser.GetSelection()
        )
        self.parent.values_dictionary[0]["p0 Coarse"] = float(
            self.parent.P0_slider.GetValue()
        )
        self.parent.values_dictionary[0]["p0 Fine"] = float(
            self.parent.P0_slider_fine.GetValue()
        )
        self.parent.values_dictionary[0]["p1 Coarse"] = float(
            self.parent.P1_slider.GetValue()
        )
        self.parent.values_dictionary[0]["p1 Fine"] = float(
            self.parent.P1_slider_fine.GetValue
        )
        self.parent.line1.set_label("1")
        self.linewidth = self.parent.line1.get_linewidth()
        self.choices = []
        self.choices.append("1")

    def plot_extra(self, data, title, index):
        # Input current values of the sliders for the first plot into the dictionary
        self.parent.values_dictionary[index] = {}
        self.parent.values_dictionary[index]["title"] = title
        self.parent.values_dictionary[index][
            "linewidth"
        ] = self.parent.line1.get_linewidth()
        self.parent.values_dictionary[index]["color index"] = self.parent.index
        self.parent.values_dictionary[index]["original_ppms"] = self.parent.ppm_original
        self.parent.values_dictionary[index]["original_data"] = self.parent.nmrdata.data
        self.parent.values_dictionary[index]["move left/right"] = float(
            self.parent.reference_slider.GetValue()
        )
        self.parent.values_dictionary[index][
            "move left/right range index"
        ] = self.parent.ref_index
        self.parent.values_dictionary[index]["move up/down"] = float(
            self.parent.vertical_slider.GetValue()
        )
        self.parent.values_dictionary[index][
            "move up/down range index"
        ] = self.parent.vertical_index
        self.parent.values_dictionary[index]["multiply value"] = float(
            self.parent.multiply_slider.GetValue()
        )
        self.parent.values_dictionary[index]["multiply value"] = float(
            self.parent.multiply_slider.GetValue()
        )
        self.parent.values_dictionary[index]["multiply range index"] = int(
            self.parent.multiply_range_chooser.GetSelection()
        )
        self.parent.values_dictionary[index]["p0 Coarse"] = float(
            self.parent.P0_slider.GetValue()
        )
        self.parent.values_dictionary[index]["p0 Fine"] = float(
            self.parent.P0_slider_fine.GetValue()
        )
        self.parent.values_dictionary[index]["p1 Coarse"] = float(
            self.parent.P1_slider.GetValue()
        )
        self.parent.values_dictionary[index]["p1 Fine"] = float(
            self.parent.P1_slider_fine.GetValue
        )
        self.parent.line1.set_label(title)
        self.linewidth = 1.5
        self.choices.append(title)


# This class is used to drop files onto the canvas
class FileDrop(wx.FileDropTarget):

    def __init__(self, canvas, axis, parent):

        wx.FileDropTarget.__init__(self)
        self.canvas = canvas
        self.axis = axis
        self.ucs = []
        self.data = []
        self.first_drop = True
        self.parent = parent
        self.parent.extra_plots = (
            []
        )  # This is a list of extra plot objects that are added to the canvas
        self.stackmode = False
        self.transposed_stack = False
        self.nmrdata_original = []

        # Create a hidden frame to be used as a parent for popout messages
        self.tempframe = wx.Frame(None, title="Temporary Parent", size=(1, 1))
        self.tempframe.Hide()  # Hide the frame since we don't need it to be visible

        if self.parent.nmrdata.dim > 1:
            self.axes1D_H = self.parent.axes1D
            self.axes1D_V = self.parent.axes1D_2
        self.color_list = [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "black",
            "navy",
            "tan",
            "lightcoral",
            "maroon",
            "lightgreen",
            "deeppink",
            "fuchsia",
        ]
        self.custom_lines = []
        self.custom_labels = []
        self.parent.active_plot_index = 0

    def OnDropFiles(self, x, y, filenames):

        for name in filenames:
            bruker = False
            if os.path.isdir(name):
                files = os.listdir(name)
                self.brukerdata = False
                for file in files:
                    if file in [
                        "1r",
                        "1i",
                        "2rr",
                        "2ri",
                        "3rrr",
                        "3rri",
                        "3rir",
                        "3rii",
                        "3irr",
                        "3iri",
                        "3iir",
                        "3iii",
                    ]:
                        bruker = True
                        break
            if ".dat" in name or ".ft" in name or bruker == True:
                if self.stackmode == False:
                    if bruker == False:
                        dic, data = ng.pipe.read(name)
                    else:
                        dic, data = ng.bruker.read_pdata(name)
                    if len(data.shape) == 1:
                        if self.parent.nmrdata.dim != 1:
                            msg = "Cannot drop 1D data onto a 2D/3D plot"
                            dlg = wx.MessageDialog(
                                None, msg, "Error", wx.OK | wx.ICON_ERROR
                            )
                            dlg.ShowModal()
                            dlg.Destroy()
                            return False
                        if self.first_drop:
                            msg = "Entering multiple plot mode: Please enter title of the first dataset"
                            dlg = wx.TextEntryDialog(None, msg)
                            res = dlg.ShowModal()
                            if res == wx.ID_CANCEL:
                                return False

                            # Input current values of the sliders for the first plot into the dictionary
                            self.parent.values_dictionary[0] = {}
                            self.parent.values_dictionary[0]["title"] = dlg.GetValue()
                            self.parent.values_dictionary[0][
                                "linewidth"
                            ] = self.parent.line1.get_linewidth()
                            self.parent.values_dictionary[0][
                                "color index"
                            ] = self.parent.index
                            self.parent.values_dictionary[0][
                                "original_ppms"
                            ] = self.parent.ppm_original
                            self.parent.values_dictionary[0][
                                "original_data"
                            ] = self.parent.nmrdata.data
                            self.parent.values_dictionary[0][
                                "dictionary"
                            ] = self.parent.nmrdata.dic
                            self.parent.values_dictionary[0]["move left/right"] = float(
                                self.parent.reference_slider.GetValue()
                            )
                            self.parent.values_dictionary[0][
                                "move left/right range index"
                            ] = self.parent.ref_index
                            self.parent.values_dictionary[0]["move up/down"] = float(
                                self.parent.vertical_slider.GetValue()
                            )
                            self.parent.values_dictionary[0][
                                "move up/down range index"
                            ] = self.parent.vertical_index
                            self.parent.values_dictionary[0]["multiply value"] = float(
                                self.parent.multiply_slider.GetValue()
                            )
                            self.parent.values_dictionary[0]["multiply range index"] = (
                                int(self.parent.multiply_range_chooser.GetSelection())
                            )
                            self.parent.values_dictionary[0]["p0 Coarse"] = float(
                                self.parent.P0_slider.GetValue()
                            )
                            self.parent.values_dictionary[0]["p0 Fine"] = float(
                                self.parent.P0_slider_fine.GetValue()
                            )
                            self.parent.values_dictionary[0]["p1 Coarse"] = float(
                                self.parent.P1_slider.GetValue()
                            )
                            self.parent.values_dictionary[0]["p1 Fine"] = float(
                                self.parent.P1_slider_fine.GetValue()
                            )
                            self.parent.values_dictionary[0]["dictionary"] = dic
                            try:
                                if self.parent.parent.parent.path != "":
                                    path = self.parent.parent.parent.path
                                else:
                                    path = os.getcwd()
                            except:
                                path = os.getcwd()
                            if platform == "windows":
                                self.parent.values_dictionary[0]["path"] = (
                                    path + "\\" + self.parent.nmrdata.file
                                )
                            else:
                                self.parent.values_dictionary[0]["path"] = (
                                    path + "/" + self.parent.nmrdata.file
                                )

                            self.parent.line1.set_label(dlg.GetValue())
                            self.linewidth = self.parent.line1.get_linewidth()
                            self.choices = []
                            self.choices.append(dlg.GetValue())
                            self.first_drop = False

                        if bruker == False:
                            uc0 = ng.pipe.make_uc(dic, data, dim=0)
                        else:
                            udic = ng.bruker.guess_udic(dic, data)
                            uc0 = ng.fileiobase.uc_from_udic(udic)
                        self.data.append(data)
                        x0, x1 = uc0.ppm_limits()
                        uc0.ppms_scale = np.linspace(x0, x1, int(uc0._size))
                        msg = "Please enter title of this data!"
                        dlg = wx.TextEntryDialog(self.tempframe, msg)
                        self.tempframe.Raise()
                        self.tempframe.SetFocus()
                        res = dlg.ShowModal()
                        if res == wx.ID_CANCEL:
                            self.canvas.draw_idle()
                            return False

                        # Add default values for the new plot to the values dictionary
                        self.parent.values_dictionary[
                            len(self.parent.extra_plots) + 1
                        ] = {}
                        self.parent.values_dictionary[len(self.parent.extra_plots) + 1][
                            "title"
                        ] = dlg.GetValue()
                        self.parent.values_dictionary[len(self.parent.extra_plots) + 1][
                            "linewidth"
                        ] = 0.5
                        self.parent.values_dictionary[len(self.parent.extra_plots) + 1][
                            "color index"
                        ] = (len(self.parent.extra_plots) + 1)
                        self.parent.values_dictionary[len(self.parent.extra_plots) + 1][
                            "original_ppms"
                        ] = uc0.ppms_scale
                        self.parent.values_dictionary[len(self.parent.extra_plots) + 1][
                            "original_data"
                        ] = data
                        self.parent.values_dictionary[len(self.parent.extra_plots) + 1][
                            "dictionary"
                        ] = dic
                        self.parent.values_dictionary[len(self.parent.extra_plots) + 1][
                            "move left/right"
                        ] = 0
                        self.parent.values_dictionary[len(self.parent.extra_plots) + 1][
                            "move left/right range index"
                        ] = self.parent.ref_index
                        self.parent.values_dictionary[len(self.parent.extra_plots) + 1][
                            "move up/down"
                        ] = 0
                        self.parent.values_dictionary[len(self.parent.extra_plots) + 1][
                            "move up/down range index"
                        ] = self.parent.vertical_index
                        self.parent.values_dictionary[len(self.parent.extra_plots) + 1][
                            "multiply value"
                        ] = 1
                        self.parent.values_dictionary[len(self.parent.extra_plots) + 1][
                            "multiply range index"
                        ] = 0
                        self.parent.values_dictionary[len(self.parent.extra_plots) + 1][
                            "p0 Coarse"
                        ] = 0
                        self.parent.values_dictionary[len(self.parent.extra_plots) + 1][
                            "p0 Fine"
                        ] = 0
                        self.parent.values_dictionary[len(self.parent.extra_plots) + 1][
                            "p1 Coarse"
                        ] = 0
                        self.parent.values_dictionary[len(self.parent.extra_plots) + 1][
                            "p1 Fine"
                        ] = 0
                        self.parent.values_dictionary[len(self.parent.extra_plots) + 1][
                            "path"
                        ] = name
                        self.parent.values_dictionary[len(self.parent.extra_plots) + 1][
                            "dictionary"
                        ] = dic

                        # Add labels of the extra plots to the select plot box
                        self.choices.append(dlg.GetValue())
                        self.parent.plot_combobox.Clear()
                        self.parent.plot_combobox.AppendItems(self.choices)
                        self.parent.plot_combobox.SetSelection(0)
                        xlim, ylim = self.axis.get_xlim(), self.axis.get_ylim()
                        if len(self.parent.extra_plots) + 1 < len(self.color_list):
                            self.parent.extra_plots.append(
                                self.axis.plot(
                                    uc0.ppms_scale,
                                    data,
                                    color=self.color_list[
                                        len(self.parent.extra_plots) + 1
                                    ],
                                    label=dlg.GetValue(),
                                    linewidth=0.5,
                                )
                            )
                        else:
                            self.parent.values_dictionary[
                                len(self.parent.extra_plots) + 1
                            ]["color index"] = (
                                len(self.parent.extra_plots) + 1 - len(self.color_list)
                            )
                            self.parent.extra_plots.append(
                                self.axis.plot(
                                    uc0.ppms_scale,
                                    data,
                                    color=self.color_list[
                                        len(self.parent.extra_plots)
                                        + 1
                                        - len(self.color_list)
                                    ],
                                    label=dlg.GetValue(),
                                    linewidth=0.5,
                                )
                            )

                        self.axis.legend()
                        self.axis.set_xlim(xlim)
                        self.axis.set_ylim(ylim)

                        self.parent.OnSelectPlot(wx.EVT_COMBOBOX)

                        self.canvas.draw()

                    elif len(data.shape) == 2:
                        if self.parent.nmrdata.dim != 2:
                            msg = "Please enter a 2D dataset"
                            dlg = wx.MessageDialog(
                                self.tempframe, msg, "Error", wx.OK | wx.ICON_ERROR
                            )
                            self.tempframe.Raise()
                            self.tempframe.SetFocus()
                            dlg.ShowModal()
                            dlg.Destroy()
                            return False
                        if len(self.parent.twoD_spectra) == 12:
                            msg = "Maximum number of overlayed 2D plots reached (12)"
                            dlg = wx.MessageDialog(
                                self.tempframe, msg, "Error", wx.OK | wx.ICON_ERROR
                            )
                            self.tempframe.Raise()
                            self.tempframe.SetFocus()
                            dlg.ShowModal()
                            dlg.Destroy()
                            return False

                        if self.parent.transposed2D == True:
                            self.transposed = True
                            self.parent.do_not_update = True
                            self.parent.OnTransposeButton(wx.EVT_BUTTON)
                        else:
                            self.transposed = False

                        msg = "Please enter title of the dropped data!"
                        dlg_1 = wx.TextEntryDialog(self.tempframe, msg)
                        self.tempframe.Raise()
                        self.tempframe.SetFocus()
                        res = dlg_1.ShowModal()

                        if res == wx.ID_CANCEL:
                            return False

                        try:
                            self.parent.parent.parent.Raise()
                            dlg_1.Raise()
                        except:
                            pass

                        if self.first_drop:
                            msg = "Please enter title of the first dataset"
                            dlg = wx.TextEntryDialog(self.tempframe, msg)
                            self.tempframe.Raise()
                            self.tempframe.SetFocus()
                            res = dlg.ShowModal()
                            if res == wx.ID_CANCEL:
                                self.first_drop = True
                                return False
                            self.parent.values_dictionary[0] = {}
                            self.parent.values_dictionary[0]["title"] = dlg.GetValue()
                            self.parent.values_dictionary[0][
                                "move x"
                            ] = self.parent.move_x_slider.GetValue()
                            self.parent.values_dictionary[0][
                                "move y"
                            ] = self.parent.move_y_slider.GetValue()
                            if self.parent.reference_range_chooserX.GetSelection() < 0:
                                self.parent.values_dictionary[0][
                                    "move x range index"
                                ] = 0
                            else:
                                self.parent.values_dictionary[0][
                                    "move x range index"
                                ] = self.parent.reference_range_chooserX.GetSelection()
                            if self.parent.reference_range_chooserY.GetSelection() < 0:
                                self.parent.values_dictionary[0][
                                    "move y range index"
                                ] = 0
                            else:
                                self.parent.values_dictionary[0][
                                    "move y range index"
                                ] = self.parent.reference_range_chooserY.GetSelection()
                            self.parent.values_dictionary[0][
                                "p0 Coarse"
                            ] = self.parent.P0_slider.GetValue()
                            self.parent.values_dictionary[0][
                                "p0 Fine"
                            ] = self.parent.P0_slider_fine.GetValue()
                            self.parent.values_dictionary[0][
                                "p1 Coarse"
                            ] = self.parent.P1_slider.GetValue()
                            self.parent.values_dictionary[0][
                                "p1 Fine"
                            ] = self.parent.P1_slider_fine.GetValue()
                            self.parent.values_dictionary[0][
                                "original_x_ppms"
                            ] = self.parent.ppms_0
                            self.parent.values_dictionary[0][
                                "original_y_ppms"
                            ] = self.parent.ppms_1
                            self.parent.values_dictionary[0]["new_x_ppms"] = (
                                self.parent.ppms_0
                                + np.ones(len(self.parent.ppms_0))
                                * self.parent.move_x_slider.GetValue()
                            )
                            self.parent.values_dictionary[0]["new_y_ppms"] = (
                                self.parent.ppms_1
                                + np.ones(len(self.parent.ppms_1))
                                * self.parent.move_y_slider.GetValue()
                            )
                            self.parent.values_dictionary[0][
                                "z_data"
                            ] = self.parent.nmrdata.data
                            self.parent.values_dictionary[0][
                                "contour linewidth"
                            ] = self.parent.contour_width_slider.GetValue()
                            self.parent.values_dictionary[0][
                                "linewidth 1D"
                            ] = self.parent.line_width_slider.GetValue()
                            self.parent.values_dictionary[0]["uc0"] = self.parent.uc0
                            self.parent.values_dictionary[0]["uc1"] = self.parent.uc1
                            self.parent.values_dictionary[0][
                                "multiply factor"
                            ] = self.parent.multiply_factor
                            self.parent.values_dictionary[0][
                                "contour levels"
                            ] = self.parent.contour_levels_slider.GetValue()
                            self.parent.values_dictionary[0]["transposed"] = False
                            try:
                                if self.parent.parent.parent.path != "":
                                    path = self.parent.parent.parent.path
                                else:
                                    path = os.getcwd()
                            except:
                                path = os.getcwd()
                            if platform == "windows":
                                self.parent.values_dictionary[0]["path"] = (
                                    path + "\\" + self.parent.nmrdata.file
                                )
                            else:
                                self.parent.values_dictionary[0]["path"] = (
                                    path + "/" + self.parent.nmrdata.file
                                )

                            # Turn on multiplot mode
                            self.parent.multiplot_mode = True

                            # Create labels for the 2D contour plots

                            self.custom_lines.append(
                                Line2D(
                                    [0],
                                    [0],
                                    color=self.parent.twoD_label_colours[0],
                                    lw=1.5,
                                )
                            )
                            self.custom_labels.append(dlg.GetValue())

                            self.first_drop = False

                        try:
                            self.parent.parent.parent.Raise()
                        except:
                            pass

                        if bruker == False:
                            uc0 = ng.pipe.make_uc(dic, data, dim=0)
                            uc1 = ng.pipe.make_uc(dic, data, dim=1)
                        else:
                            udic = ng.bruker.guess_udic(dic, data)
                            uc0 = ng.fileiobase.uc_from_udic(udic, dim=0)
                            uc1 = ng.fileiobase.uc_from_udic(udic, dim=1)
                        ppm0 = uc0.ppm_scale()
                        ppm1 = uc1.ppm_scale()
                        x, y = np.meshgrid(ppm1, ppm0)

                        if len(self.parent.twoD_spectra) == 0:
                            index = 1
                        else:
                            index = len(self.parent.twoD_spectra)
                        self.parent.values_dictionary[index] = {}
                        self.parent.values_dictionary[index]["title"] = dlg_1.GetValue()
                        self.parent.values_dictionary[index]["move x"] = 0
                        self.parent.values_dictionary[index]["move y"] = 0
                        self.parent.values_dictionary[index]["move x range index"] = 0
                        self.parent.values_dictionary[index]["move y range index"] = 0
                        self.parent.values_dictionary[index]["p0 Coarse"] = 0
                        self.parent.values_dictionary[index]["p0 Fine"] = 0
                        self.parent.values_dictionary[index]["p1 Coarse"] = 0
                        self.parent.values_dictionary[index]["p1 Fine"] = 0
                        self.parent.values_dictionary[index]["original_x_ppms"] = ppm0
                        self.parent.values_dictionary[index]["original_y_ppms"] = ppm1
                        self.parent.values_dictionary[index]["new_x_ppms"] = ppm0
                        self.parent.values_dictionary[index]["new_y_ppms"] = ppm1
                        self.parent.values_dictionary[index]["z_data"] = data
                        self.parent.values_dictionary[index]["contour linewidth"] = 1.0
                        self.parent.values_dictionary[index]["linewidth 1D"] = 1.0
                        self.parent.values_dictionary[index]["uc0"] = uc0
                        self.parent.values_dictionary[index]["uc1"] = uc1
                        self.parent.values_dictionary[index]["multiply factor"] = 1.0
                        self.parent.values_dictionary[index]["contour levels"] = 20
                        self.parent.values_dictionary[index]["path"] = name

                        # Work out the difference in max intensities between the first and the added spectra
                        max_intensity = np.max(data)
                        max_intensity_0 = np.max(self.parent.nmrdata.data)
                        if max_intensity_0 > max_intensity:
                            max_intensity_diff = max_intensity / max_intensity_0

                            if max_intensity_diff < 0.1:
                                # The max intensity of the added spectrum is less than 10% of the max intensity of the first spectrum
                                dlg = wx.MessageDialog(
                                    None,
                                    "The maximum intensity of the new spectrum is less than 10% of the maximum intensity of the first spectrum. Do you want to scale the new spectrum to the first spectrum?",
                                    "Warning",
                                    wx.YES_NO | wx.ICON_WARNING,
                                )
                                res = dlg.ShowModal()
                                if res == wx.ID_YES:
                                    self.parent.values_dictionary[index][
                                        "multiply factor"
                                    ] = (1 / max_intensity_diff)
                                else:
                                    self.parent.values_dictionary[index][
                                        "multiply factor"
                                    ] = 1
                            else:
                                self.parent.values_dictionary[index][
                                    "multiply factor"
                                ] = 1

                        else:
                            max_intensity_diff = max_intensity_0 / max_intensity
                            if max_intensity_diff < 0.1:
                                # The max intensity of the first spectrum is less than 10% of the max intensity of the added spectrum
                                dlg = wx.MessageDialog(
                                    None,
                                    "The maximum intensity of the first spectrum is less than 10% of the maximum intensity of the new spectrum. Do you want to scale the first spectrum to the new spectrum?",
                                    "Warning",
                                    wx.YES_NO | wx.ICON_WARNING,
                                )
                                res = dlg.ShowModal()
                                if res == wx.ID_YES:
                                    self.parent.values_dictionary[0][
                                        "multiply factor"
                                    ] = max_intensity_diff
                                else:
                                    self.parent.values_dictionary[0][
                                        "multiply factor"
                                    ] = 1
                            else:
                                self.parent.values_dictionary[0]["multiply factor"] = 1

                        if len(self.parent.twoD_spectra) == 0:
                            self.custom_lines.append(
                                Line2D(
                                    [0],
                                    [0],
                                    color=self.parent.twoD_label_colours[
                                        len(self.parent.twoD_spectra) + 1
                                    ],
                                    lw=1.5,
                                )
                            )
                        else:
                            self.custom_lines.append(
                                Line2D(
                                    [0],
                                    [0],
                                    color=self.parent.twoD_label_colours[
                                        len(self.parent.twoD_spectra)
                                    ],
                                    lw=1.5,
                                )
                            )
                        self.custom_labels.append(dlg_1.GetValue())

                        xlim, ylim = self.axis.get_xlim(), self.axis.get_ylim()
                        xlabel = self.axis.get_xlabel()
                        ylabel = self.axis.get_ylabel()
                        self.axis.clear()
                        self.parent.axes1D.clear()
                        self.parent.axes1D_2.clear()
                        self.parent.axes1D.set_yticks([])
                        self.parent.axes1D_2.set_xticks([])

                        self.parent.twoD_spectra = []
                        self.parent.twoD_slices_horizontal = []
                        self.parent.twoD_slices_vertical = []
                        length = len(self.parent.values_dictionary.keys())
                        for i in range(len(self.parent.values_dictionary)):
                            multiply_factor = self.parent.values_dictionary[i][
                                "multiply factor"
                            ]

                            # If transpose is false, then the x-axis is the first axis and the y-axis is the second axis
                            if self.parent.transposed2D == False:
                                self.parent.values_dictionary[i]["new_x_ppms_old"] = (
                                    self.parent.values_dictionary[i]["new_x_ppms"]
                                )
                                self.parent.values_dictionary[i]["new_y_ppms_old"] = (
                                    self.parent.values_dictionary[i]["new_y_ppms"]
                                )
                                if (
                                    np.abs(
                                        (
                                            self.parent.values_dictionary[i][
                                                "new_x_ppms"
                                            ][0]
                                            - self.parent.values_dictionary[0][
                                                "new_x_ppms"
                                            ][0]
                                        )
                                        / self.parent.values_dictionary[0][
                                            "new_x_ppms"
                                        ][0]
                                    )
                                    > 0.2
                                ):
                                    # More than 20% difference in the x-axis (consider transposing new spectra)
                                    # Give a popout asking if the user wants to still add the new spectrum
                                    msg = "The x-axis of the new spectrum is significantly different from the x-axis of the first spectrum. Do you want to transpose the new spectrum?"
                                    dlg = wx.MessageDialog(
                                        None,
                                        msg,
                                        "Warning",
                                        wx.YES_NO | wx.ICON_WARNING,
                                    )
                                    res = dlg.ShowModal()
                                    if res == wx.ID_YES:
                                        transpose = True
                                    else:
                                        transpose = False

                                    self.parent.values_dictionary[i][
                                        "transposed"
                                    ] = transpose

                                    if transpose == True:
                                        self.parent.values_dictionary[i][
                                            "new_x_ppms_old"
                                        ] = self.parent.values_dictionary[i][
                                            "new_x_ppms"
                                        ]
                                        self.parent.values_dictionary[i][
                                            "new_y_ppms_old"
                                        ] = self.parent.values_dictionary[i][
                                            "new_y_ppms"
                                        ]
                                        self.parent.values_dictionary[i][
                                            "new_x_ppms"
                                        ] = self.parent.values_dictionary[i][
                                            "new_y_ppms_old"
                                        ]
                                        self.parent.values_dictionary[i][
                                            "new_y_ppms"
                                        ] = self.parent.values_dictionary[i][
                                            "new_x_ppms_old"
                                        ]
                                        self.parent.values_dictionary[i][
                                            "original_x_ppms"
                                        ] = self.parent.values_dictionary[i][
                                            "new_x_ppms"
                                        ]
                                        self.parent.values_dictionary[i][
                                            "original_y_ppms"
                                        ] = self.parent.values_dictionary[i][
                                            "original_y_ppms"
                                        ]
                                        uc0 = self.parent.values_dictionary[i]["uc1"]
                                        uc1 = self.parent.values_dictionary[i]["uc0"]
                                        self.parent.values_dictionary[i]["uc0"] = uc0
                                        self.parent.values_dictionary[i]["uc1"] = uc1
                                        self.parent.values_dictionary[i]["z_data"] = (
                                            self.parent.values_dictionary[i]["z_data"].T
                                        )

                                x, y = np.meshgrid(
                                    self.parent.values_dictionary[i]["new_y_ppms"],
                                    self.parent.values_dictionary[i]["new_x_ppms"],
                                )
                                self.parent.twoD_spectra.append(
                                    self.axis.contour(
                                        y,
                                        x,
                                        self.parent.values_dictionary[i]["z_data"]
                                        * multiply_factor,
                                        colors=self.parent.twoD_colours[i],
                                        levels=self.parent.cl,
                                        linewidths=self.parent.values_dictionary[i][
                                            "contour linewidth"
                                        ],
                                    )
                                )
                            else:
                                self.parent.values_dictionary[i]["new_x_ppms_old"] = (
                                    self.parent.values_dictionary[i]["new_x_ppms"]
                                )
                                self.parent.values_dictionary[i]["new_y_ppms_old"] = (
                                    self.parent.values_dictionary[i]["new_y_ppms"]
                                )
                                self.parent.values_dictionary[i]["new_x_ppms"] = (
                                    self.parent.values_dictionary[i]["new_y_ppms_old"]
                                )
                                self.parent.values_dictionary[i]["new_y_ppms"] = (
                                    self.parent.values_dictionary[i]["new_x_ppms_old"]
                                )
                                x, y = np.meshgrid(
                                    self.parent.values_dictionary[i]["new_y_ppms"],
                                    self.parent.values_dictionary[i]["new_x_ppms"],
                                )
                                if i > len(self.parent.values_dictionary.keys()) - 1:
                                    self.parent.twoD_spectra.append(
                                        self.axis.contour(
                                            x,
                                            y,
                                            self.parent.values_dictionary[i]["z_data"].T
                                            * multiply_factor,
                                            colors=self.parent.twoD_colours[i],
                                            levels=self.parent.cl,
                                            linewidths=self.parent.values_dictionary[i][
                                                "contour linewidth"
                                            ],
                                        )
                                    )
                                else:
                                    try:
                                        self.parent.twoD_spectra.append(
                                            self.axis.contour(
                                                x,
                                                y,
                                                self.parent.values_dictionary[i][
                                                    "z_data"
                                                ]
                                                * multiply_factor,
                                                colors=self.parent.twoD_colours[i],
                                                levels=self.parent.cl,
                                                linewidths=self.parent.values_dictionary[
                                                    i
                                                ][
                                                    "contour linewidth"
                                                ],
                                            )
                                        )
                                    except:
                                        self.parent.twoD_spectra.append(
                                            self.axis.contour(
                                                x,
                                                y,
                                                self.parent.values_dictionary[i][
                                                    "z_data"
                                                ].T
                                                * multiply_factor,
                                                colors=self.parent.twoD_colours[i],
                                                levels=self.parent.cl,
                                                linewidths=self.parent.values_dictionary[
                                                    i
                                                ][
                                                    "contour linewidth"
                                                ],
                                            )
                                        )

                            if self.parent.transposed2D == False:
                                self.parent.twoD_slices_horizontal.append(
                                    self.parent.axes1D.plot(
                                        self.parent.values_dictionary[i]["new_x_ppms"],
                                        self.parent.values_dictionary[i]["z_data"][:, 1]
                                        * multiply_factor,
                                        color=self.parent.twoD_label_colours[i],
                                        linewidth=self.parent.values_dictionary[i][
                                            "linewidth 1D"
                                        ],
                                    )
                                )
                                self.parent.twoD_slices_vertical.append(
                                    self.parent.axes1D_2.plot(
                                        self.parent.values_dictionary[i]["new_y_ppms"],
                                        self.parent.values_dictionary[i]["z_data"][1, :]
                                        * multiply_factor,
                                        color=self.parent.twoD_label_colours[i],
                                        linewidth=self.parent.values_dictionary[i][
                                            "linewidth 1D"
                                        ],
                                    )
                                )

                            else:
                                if i == 0:
                                    self.parent.twoD_slices_horizontal.append(
                                        self.parent.axes1D.plot(
                                            self.parent.values_dictionary[i][
                                                "new_x_ppms"
                                            ],
                                            self.parent.values_dictionary[i][
                                                "z_data"
                                            ].T[1, :]
                                            * multiply_factor,
                                            color=self.parent.twoD_label_colours[i],
                                            linewidth=self.parent.values_dictionary[i][
                                                "linewidth 1D"
                                            ],
                                        )
                                    )
                                    self.parent.twoD_slices_vertical.append(
                                        self.parent.axes1D_2.plot(
                                            self.parent.values_dictionary[i][
                                                "new_y_ppms"
                                            ],
                                            self.parent.values_dictionary[i][
                                                "z_data"
                                            ].T[:, 1]
                                            * multiply_factor,
                                            color=self.parent.twoD_label_colours[i],
                                            linewidth=self.parent.values_dictionary[i][
                                                "linewidth 1D"
                                            ],
                                        )
                                    )
                                else:
                                    self.parent.twoD_slices_horizontal.append(
                                        self.parent.axes1D.plot(
                                            self.parent.values_dictionary[i][
                                                "new_x_ppms"
                                            ],
                                            self.parent.values_dictionary[i][
                                                "z_data"
                                            ].T[:, 1]
                                            * multiply_factor,
                                            color=self.parent.twoD_label_colours[i],
                                            linewidth=self.parent.values_dictionary[i][
                                                "linewidth 1D"
                                            ],
                                        )
                                    )
                                    self.parent.twoD_slices_vertical.append(
                                        self.parent.axes1D_2.plot(
                                            self.parent.values_dictionary[i][
                                                "new_y_ppms"
                                            ],
                                            self.parent.values_dictionary[i][
                                                "z_data"
                                            ].T[1, :]
                                            * multiply_factor,
                                            color=self.parent.twoD_label_colours[i],
                                            linewidth=self.parent.values_dictionary[i][
                                                "linewidth 1D"
                                            ],
                                        )
                                    )
                        self.parent.line_h = self.axis.axhline(
                            y=self.parent.values_dictionary[i]["new_x_ppms"][1],
                            color="black",
                            lw=1.5,
                        )
                        self.parent.line_v = self.axis.axvline(
                            x=self.parent.values_dictionary[i]["new_y_ppms"][1],
                            color="black",
                            lw=1.5,
                        )
                        self.parent.line_h.set_visible(False)
                        self.parent.line_v.set_visible(False)

                        # Set all vertical and horizontal slices to invisible initially

                        for i in range(len(self.parent.twoD_slices_horizontal)):
                            self.parent.twoD_slices_horizontal[i][0].set_visible(False)
                            self.parent.twoD_slices_vertical[i][0].set_visible(False)

                        self.axis.legend(self.custom_lines, self.custom_labels)
                        if self.parent.transposed2D == False:
                            self.axis.set_xlim(xlim)
                            self.axis.set_ylim(ylim)
                            self.axis.set_xlabel(xlabel)
                            self.axis.set_ylabel(ylabel)
                        else:
                            self.axis.set_xlim(ylim)
                            self.axis.set_ylim(xlim)
                            self.axis.set_xlabel(ylabel)
                            self.axis.set_ylabel(xlabel)

                        # Add labels of the extra plots to the select plot box
                        self.parent.plot_combobox.Clear()
                        self.parent.plot_combobox.AppendItems(self.custom_labels)
                        self.parent.plot_combobox.SetSelection(0)

                        self.parent.UpdateFrame()

                        if self.transposed == True:
                            self.parent.do_not_update = False
                            self.parent.OnTransposeButton(wx.EVT_BUTTON)

                    else:
                        msg = "This is not 1D or 2D data - currently more dimensions are not supported..."
                        dlg = wx.MessageDialog(self.tempframe, msg)
                        self.tempframe.Raise()
                        self.tempframe.SetFocus()
                        dlg.ShowModal()
                        return False
                # If in stack mode overlay the 1D spectra without asking the user for a title
                elif self.stackmode == True:
                    # Input current values of the sliders for the first plot into the dictionary
                    dic, data_original = ng.pipe.read(name)
                    # data_original = data_original.T
                    if self.transposed_stack == True:
                        uc0 = ng.pipe.make_uc(dic, data_original, dim=1)
                    else:
                        uc0 = ng.pipe.make_uc(dic, data_original, dim=0)
                        data_original = data_original.T
                    while len(data_original) > len(self.color_list):
                        self.color_list = self.color_list * 2
                    x0, x1 = uc0.ppm_limits()
                    uc0.ppms_scale = np.linspace(x0, x1, int(uc0._size))
                    uc0_ppms = uc0.ppm_scale()
                    data = []
                    data.append(data_original[0])
                    self.stackfirstpoint()
                    self.parent.multiplot_mode = True
                    for i in range(len(data_original)):
                        if i == 0:
                            continue
                        else:
                            self.data.append(data_original[i])
                            # Add default values for the new plot to the values dictionary
                            self.parent.values_dictionary[
                                len(self.parent.extra_plots) + 1
                            ] = {}
                            self.parent.values_dictionary[
                                len(self.parent.extra_plots) + 1
                            ]["title"] = str(i + 1)
                            self.parent.values_dictionary[
                                len(self.parent.extra_plots) + 1
                            ]["linewidth"] = self.linewidth
                            self.parent.values_dictionary[
                                len(self.parent.extra_plots) + 1
                            ]["color index"] = (len(self.parent.extra_plots) + 1)
                            self.parent.values_dictionary[
                                len(self.parent.extra_plots) + 1
                            ]["original_ppms"] = uc0.ppm_scale()
                            self.parent.values_dictionary[
                                len(self.parent.extra_plots) + 1
                            ]["original_data"] = data_original[i]
                            self.parent.values_dictionary[
                                len(self.parent.extra_plots) + 1
                            ]["move left/right"] = 0
                            self.parent.values_dictionary[
                                len(self.parent.extra_plots) + 1
                            ]["move left/right range index"] = self.parent.ref_index
                            self.parent.values_dictionary[
                                len(self.parent.extra_plots) + 1
                            ]["move up/down"] = 0
                            self.parent.values_dictionary[
                                len(self.parent.extra_plots) + 1
                            ]["move up/down range index"] = self.parent.vertical_index
                            self.parent.values_dictionary[
                                len(self.parent.extra_plots) + 1
                            ]["multiply value"] = 1
                            self.parent.values_dictionary[
                                len(self.parent.extra_plots) + 1
                            ]["multiply range index"] = 0
                            self.parent.values_dictionary[
                                len(self.parent.extra_plots) + 1
                            ]["p0 Coarse"] = 0
                            self.parent.values_dictionary[
                                len(self.parent.extra_plots) + 1
                            ]["p0 Fine"] = 0
                            self.parent.values_dictionary[
                                len(self.parent.extra_plots) + 1
                            ]["p1 Coarse"] = 0
                            self.parent.values_dictionary[
                                len(self.parent.extra_plots) + 1
                            ]["p1 Fine"] = 0
                            self.parent.values_dictionary[
                                len(self.parent.extra_plots) + 1
                            ]["path"] = name
                            self.parent.values_dictionary[
                                len(self.parent.extra_plots) + 1
                            ]["dictionary"] = dic
                            # Add labels of the extra plots to the select plot box
                            self.choices.append(str(i + 1))
                            self.parent.plot_combobox.Clear()
                            self.parent.plot_combobox.AppendItems(self.choices)
                            self.parent.plot_combobox.SetSelection(0)
                            if len(self.parent.extra_plots) + 1 < len(self.color_list):
                                self.parent.extra_plots.append(
                                    self.axis.plot(
                                        uc0_ppms,
                                        data_original[i],
                                        color=self.color_list[
                                            len(self.parent.extra_plots) + 1
                                        ],
                                        label=str(i + 1),
                                        linewidth=self.linewidth,
                                    )
                                )
                            else:
                                self.parent.values_dictionary[
                                    len(self.extra_plots) + 1
                                ]["color index"] = (
                                    len(self.parent.extra_plots)
                                    + 1
                                    - len(self.color_list)
                                )
                                self.parent.extra_plots.append(
                                    self.axis.plot(
                                        uc0_ppms,
                                        data_original[i],
                                        color=self.color_list[
                                            len(self.parent.extra_plots)
                                            + 1
                                            - len(self.color_list)
                                        ],
                                        label=str(i + 1),
                                        linewidth=self.linewidth,
                                    )
                                )

                    self.axis.legend()
                    self.canvas.draw()
            else:

                msg = "Can only deal with nmrPipe *.ft* files!"
                dlg = wx.MessageDialog(self.tempframe, msg)
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                dlg.ShowModal()

                return False

        return True

    def stackfirstpoint(self):
        self.parent.values_dictionary[0] = {}
        self.parent.values_dictionary[0]["title"] = "1"
        self.parent.values_dictionary[0][
            "linewidth"
        ] = self.parent.line1.get_linewidth()
        self.parent.values_dictionary[0]["color index"] = self.parent.index
        self.parent.values_dictionary[0]["original_ppms"] = self.parent.ppm_original
        self.parent.values_dictionary[0]["original_data"] = self.parent.nmrdata.data
        self.parent.values_dictionary[0]["move left/right"] = float(
            self.parent.reference_slider.GetValue()
        )
        self.parent.values_dictionary[0][
            "move left/right range index"
        ] = self.parent.ref_index
        self.parent.values_dictionary[0]["move up/down"] = float(
            self.parent.vertical_slider.GetValue()
        )
        self.parent.values_dictionary[0][
            "move up/down range index"
        ] = self.parent.vertical_index
        self.parent.values_dictionary[0]["multiply value"] = float(
            self.parent.multiply_slider.GetValue()
        )
        self.parent.values_dictionary[0]["multiply range index"] = int(
            self.parent.multiply_range_chooser.GetSelection()
        )
        self.parent.values_dictionary[0]["p0 Coarse"] = float(
            self.parent.P0_slider.GetValue()
        )
        self.parent.values_dictionary[0]["p0 Fine"] = float(
            self.parent.P0_slider_fine.GetValue()
        )
        self.parent.values_dictionary[0]["p1 Coarse"] = float(
            self.parent.P1_slider.GetValue()
        )
        self.parent.values_dictionary[0]["p1 Fine"] = float(
            self.parent.P1_slider_fine.GetValue()
        )
        self.parent.values_dictionary[0]["dictionary"] = self.parent.nmrdata.dic
        self.parent.line1.set_label("1")
        self.linewidth = self.parent.line1.get_linewidth()
        self.choices = []
        self.choices.append("1")
        self.first_drop = False

        try:
            if self.parent.parent.parent.path != "":
                path = self.parent.parent.parent.path
            else:
                path = os.getcwd()
        except:
            path = os.getcwd()
        if platform == "windows":
            self.parent.values_dictionary[0]["path"] = (
                path + "\\" + self.parent.nmrdata.file
            )
        else:
            self.parent.values_dictionary[0]["path"] = (
                path + "/" + self.parent.nmrdata.file
            )


class Projection3DNotebook(wx.Notebook):
    def __init__(self, parent):
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 0.7 * self.monitorWidth
        self.height = 0.75 * self.monitorHeight
        self.parent = parent
        wx.Notebook.__init__(
            self,
            parent,
            id=wx.ID_ANY,
            style=wx.BK_DEFAULT,
            size=(self.width, self.height),
        )

        if self.parent.parent.parent.path != "":
            os.chdir(self.parent.parent.parent.path)
        # Search for the projections in the current directory (.dat files)
        self.projection_files = []
        for file in os.listdir():
            if file.endswith(".dat"):
                self.projection_files.append(file)
        for file in self.projection_files:
            if "prof" in file:
                self.projection_files.remove(file)

        self.nmrdata = []
        for file in self.projection_files:
            self.nmrdata.append(ReadProjection(file))

        if self.parent.parent.parent.cwd != "":
            os.chdir(self.parent.parent.parent.cwd)

        self.projection_selection_index = 0
        self.projection_selection_index_old = 0
        self.projection_parameters = {}
        self.projection_parameters["tab1"] = {}
        self.projection_parameters["tab2"] = {}
        self.projection_parameters["tab3"] = {}

        self.toolbar_selection = 0
        self.projection_panel1 = TwoDViewer(
            self, self.nmrdata[0], threeDprojection=True
        )
        self.toolbar_selection = 1
        self.projection_panel2 = TwoDViewer(
            self, self.nmrdata[1], threeDprojection=True
        )
        self.toolbar_selection = 2
        self.projection_panel3 = TwoDViewer(
            self, self.nmrdata[2], threeDprojection=True
        )

        self.AddPage(self.projection_panel1, self.nmrdata[0].filename.split(".dat")[0])
        self.AddPage(self.projection_panel2, self.nmrdata[1].filename.split(".dat")[0])
        self.AddPage(self.projection_panel3, self.nmrdata[2].filename.split(".dat")[0])

        self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.OnPageChanged)

    def OnPageChanged(self, event):
        self.projection_selection_index = self.GetSelection()


class ReadProjection:
    def __init__(self, filename):
        self.filename = filename
        self.file = filename

        self.read_data()
        self.dim = self.get_dimensions()
        self.get_axislabels()

    # Read in the NMRPipe data file
    def read_data(self):
        self.dic, self.data = ng.pipe.read(self.filename)

    # Work out NMR spectrum dimensions in order to get the plotting correct (need contour plot for 2D/3D but not for 1D)
    def get_dimensions(self):
        if type(self.data[0]) == np.float32:
            return 1
        if len(self.data.shape) == 2:
            return 2
        if len(self.data.shape) == 3:
            return 3

    def get_axislabels(self):
        self.axislabels = []
        file_split = self.filename.split(".dat")[0].split(".")
        for i in range(len(file_split)):
            self.axislabels.append(file_split[i])


class ProjectionFrame(wx.Frame):
    def __init__(self, title, parent=None):
        self.main_frame = parent
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()

        width = int(1.0 * self.monitorWidth)
        height = int(0.85 * self.monitorHeight)
        self.parent = parent
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.display_index_current = wx.Display.GetFromWindow(self)
        self.notebook = Projection3DNotebook(self)
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.main_sizer.AddSpacer(10)
        self.main_sizer.Add(self.notebook, 1, wx.EXPAND)

        self.SetSizerAndFit(self.main_sizer)
        self.Show()
        self.Centre()

        # Bind method to check/resize the window when the frame is moved
        self.Bind(wx.EVT_MOVE, self.OnMoveFrame)

        # Bind method to resize the window when the frame is resized
        self.Bind(wx.EVT_SIZE, self.OnSizeFrame)

    def OnMoveFrame(self, event):
        # Get the new default display if the frame is moved
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        display_index = wx.Display.GetFromWindow(self)
        if display_index != self.display_index_current:
            self.display_index_current = display_index
            self.width = int(1.0 * sizes[display_index][0])
            self.height = int(0.875 * sizes[display_index][1])
            self.SetSize((self.width, self.height))
            self.notebook.projection_panel1.canvas.SetSize(
                (
                    self.width * 0.0104,
                    (
                        self.height
                        - self.notebook.projection_panel1.bottom_sizer.GetMinSize()[1]
                        - 100
                    )
                    * 0.0104,
                )
            )
            self.notebook.projection_panel1.fig.set_size_inches(
                self.width * 0.0104,
                (
                    self.height
                    - self.notebook.projection_panel1.bottom_sizer.GetMinSize()[1]
                    - 100
                )
                * 0.0104,
            )
            self.notebook.projection_panel2.canvas.SetSize(
                (
                    self.width * 0.0104,
                    (
                        self.height
                        - self.notebook.projection_panel2.bottom_sizer.GetMinSize()[1]
                        - 100
                    )
                    * 0.0104,
                )
            )
            self.notebook.projection_panel2.fig.set_size_inches(
                self.width * 0.0104,
                (
                    self.height
                    - self.notebook.projection_panel2.bottom_sizer.GetMinSize()[1]
                    - 100
                )
                * 0.0104,
            )
            self.notebook.projection_panel3.canvas.SetSize(
                (
                    self.width * 0.0104,
                    (
                        self.height
                        - self.notebook.projection_panel3.bottom_sizer.GetMinSize()[1]
                        - 100
                    )
                    * 0.0104,
                )
            )
            self.notebook.projection_panel3.fig.set_size_inches(
                self.width * 0.0104,
                (
                    self.height
                    - self.notebook.projection_panel3.bottom_sizer.GetMinSize()[1]
                    - 100
                )
                * 0.0104,
            )
            self.UpdateProjectionFrame()
        event.Skip()

    def OnSizeFrame(self, event):
        # Get the new frame size
        self.width, self.height = self.GetSize()
        self.SetSize((self.width, self.height))
        self.notebook.projection_panel1.canvas.SetSize(
            (
                self.width * 0.0104,
                (
                    self.height
                    - self.notebook.projection_panel1.bottom_sizer.GetMinSize()[1]
                    - 100
                )
                * 0.0104,
            )
        )
        self.notebook.projection_panel1.fig.set_size_inches(
            self.width * 0.0104,
            (
                self.height
                - self.notebook.projection_panel1.bottom_sizer.GetMinSize()[1]
                - 100
            )
            * 0.0104,
        )
        self.notebook.projection_panel2.canvas.SetSize(
            (
                self.width * 0.0104,
                (
                    self.height
                    - self.notebook.projection_panel2.bottom_sizer.GetMinSize()[1]
                    - 100
                )
                * 0.0104,
            )
        )
        self.notebook.projection_panel2.fig.set_size_inches(
            self.width * 0.0104,
            (
                self.height
                - self.notebook.projection_panel2.bottom_sizer.GetMinSize()[1]
                - 100
            )
            * 0.0104,
        )
        self.notebook.projection_panel3.canvas.SetSize(
            (
                self.width * 0.0104,
                (
                    self.height
                    - self.notebook.projection_panel3.bottom_sizer.GetMinSize()[1]
                    - 100
                )
                * 0.0104,
            )
        )
        self.notebook.projection_panel3.fig.set_size_inches(
            self.width * 0.0104,
            (
                self.height
                - self.notebook.projection_panel3.bottom_sizer.GetMinSize()[1]
                - 100
            )
            * 0.0104,
        )
        self.UpdateProjectionFrame()
        event.Skip()

    def UpdateProjectionFrame(self):
        self.notebook.projection_panel1.UpdateFrame()
        self.notebook.projection_panel2.UpdateFrame()
        self.notebook.projection_panel3.UpdateFrame()


class WaterfallFrame(wx.Frame):
    def __init__(self, title, parent=None, visible="line1"):
        self.main_frame = parent
        self.visible = visible
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        width = int(1.0 * self.monitorWidth)
        height = int(0.85 * self.monitorHeight)
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.display_index_current = wx.Display.GetFromWindow(self)
        self.panel_waterfall = wx.Panel(self, -1)
        self.main_waterfall_sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_waterfall_sizer)

        self.fig_waterfall = Figure()
        self.canvas_waterfall = FigCanvas(self, -1, self.fig_waterfall)
        self.main_waterfall_sizer.Add(self.canvas_waterfall, 10, flag=wx.GROW)
        self.toolbar_waterfall = NavigationToolbar(self.canvas_waterfall)
        self.main_waterfall_sizer.Add(self.toolbar_waterfall, 0, wx.EXPAND)
        self.sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.create_waterfall_plot_sizer()

        self.titlecolor = "black"

        self.plot_waterfall()
        self.Show()

        # Bind method to check/resize the window when the frame is moved
        self.Bind(wx.EVT_MOVE, self.OnMoveFrame)

        # Bind method to resize the window when the frame is resized
        self.Bind(wx.EVT_SIZE, self.OnSizeFrame)

    def OnMoveFrame(self, event):
        # Get the new default display if the frame is moved
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        display_index = wx.Display.GetFromWindow(self)
        if display_index != self.display_index_current:
            self.display_index_current = display_index
            self.width = int(1.0 * sizes[display_index][0])
            self.height = int(0.875 * sizes[display_index][1])
            self.SetSize((self.width, self.height))
            self.canvas_waterfall.SetSize(
                (
                    self.width * 0.0104,
                    (self.height - self.sizer.GetMinSize()[1] - 100) * 0.0104,
                )
            )
            self.fig_waterfall.set_size_inches(
                self.width * 0.0104,
                (self.height - self.sizer.GetMinSize()[1] - 100) * 0.0104,
            )
            self.UpdateWaterfallFrame()
        event.Skip()

    def OnSizeFrame(self, event):
        # Get the new frame size
        self.width, self.height = self.GetSize()
        self.SetSize((self.width, self.height))
        self.canvas_waterfall.SetSize(
            (
                self.width * 0.0104,
                (self.height - self.sizer.GetMinSize()[1] - 100) * 0.0104,
            )
        )
        self.fig_waterfall.set_size_inches(
            self.width * 0.0104,
            (self.height - self.sizer.GetMinSize()[1] - 100) * 0.0104,
        )
        self.UpdateWaterfallFrame()
        event.Skip()

    def create_waterfall_plot_sizer(self):
        # Have a slider for adjusting the y axis range
        self.y_range_label = wx.StaticBox(self, -1, "Y-axis zoom")
        self.y_range_sizer = wx.StaticBoxSizer(self.y_range_label, wx.VERTICAL)
        self.y_range_slider = FloatSlider(
            self, id=-1, value=1, minval=0, maxval=3, res=0.01, style=wx.SL_HORIZONTAL
        )
        self.y_range_slider.Bind(wx.EVT_SLIDER, self.OnYRangeSlider)
        self.y_range_sizer.Add(self.y_range_slider)
        self.sizer.Add(self.y_range_sizer)
        self.main_waterfall_sizer.Add(self.sizer)
        # # Make a slider to change the contour levels
        # self.contour_label = wx.StaticBox(self, -1, "Contour levels")
        # self.contour_sizer = wx.StaticBoxSizer(self.contour_label, wx.VERTICAL)
        # self.contour_slider = FloatSlider(self, id=-1, value=1, minval=0, maxval=3, res=0.01,style=wx.SL_HORIZONTAL)
        # self.contour_slider.Bind(wx.EVT_SLIDER, self.OnContourSlider)
        # self.contour_sizer.Add(self.contour_slider)
        # self.sizer.Add(self.contour_sizer)
        # self.main_waterfall_sizer.Add(self.sizer)

    def plot_waterfall(self):
        self.ax = self.fig_waterfall.add_subplot(111)

        # Get all the slices along the pseudo3D for the currently selected slice
        if self.visible == "line1":
            vals = []
            for i in range(len(self.main_frame.nmrdata.data)):
                vals.append(
                    self.main_frame.nmrdata.data[i][
                        :, self.main_frame.uc1(str(self.main_frame.y1) + "ppm")
                    ]
                )
            for i in range(len(vals)):
                self.ax.plot(
                    self.main_frame.line1.get_xdata(), vals[i], label=str(i + 1)
                )
            self.ax.set_xlabel("ppm")
            self.ax.set_ylabel("Intensity")
            self.ax.legend()
        elif self.visible == "line3":
            vals = []
            for i in range(len(self.main_frame.nmrdata.data)):
                vals.append(
                    self.main_frame.nmrdata.data[i][
                        self.main_frame.uc0(str(self.main_frame.x1) + "ppm"), :
                    ]
                )
            for i in range(len(vals)):
                self.ax.plot(
                    self.main_frame.line3.get_ydata(), vals[i], label=str(i + 1)
                )
            self.ax.set_xlabel("ppm")
            self.ax.set_ylabel("Intensity")
            self.ax.legend()

        self.UpdateWaterfallFrame()

    def OnYRangeSlider(self, event):
        pass

    def UpdateWaterfallFrame(self):
        self.canvas_waterfall.draw()
        self.canvas_waterfall.Refresh()
        self.canvas_waterfall.Update()
        self.panel_waterfall.Refresh()
        self.panel_waterfall.Update()


class Plot3DFrame(wx.Frame):
    def __init__(self, title, parent=None):
        self.main_frame = parent
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        width = int(1.0 * self.monitorWidth)
        height = int(0.85 * self.monitorHeight)
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.display_index_current = wx.Display.GetFromWindow(self)
        self.panel_3d = wx.Panel(self, -1)
        self.main_3d_sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_3d_sizer)

        self.fig_3d = Figure()
        self.canvas_3d = FigCanvas(self, -1, self.fig_3d)
        self.main_3d_sizer.Add(self.canvas_3d, 10, flag=wx.GROW)
        self.toolbar_3d = NavigationToolbar(self.canvas_3d)
        self.main_3d_sizer.Add(self.toolbar_3d, 0, wx.EXPAND)
        self.sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.create_3D_plot_sizer()

        self.titlecolor = "black"

        self.plot3d()
        self.Show()

        # Bind method to check/resize the window when the frame is moved
        self.Bind(wx.EVT_MOVE, self.OnMoveFrame)

        # Bind method to resize the window when the frame is resized
        self.Bind(wx.EVT_SIZE, self.OnSizeFrame)

    def OnMoveFrame(self, event):
        # Get the new default display if the frame is moved
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        display_index = wx.Display.GetFromWindow(self)
        if display_index != self.display_index_current:
            self.display_index_current = display_index
            self.width = int(1.0 * sizes[display_index][0])
            self.height = int(0.875 * sizes[display_index][1])
            self.SetSize((self.width, self.height))
            self.canvas_3d.SetSize(
                (
                    self.width * 0.0104,
                    (self.height - self.sizer.GetMinSize()[1] - 100) * 0.0104,
                )
            )
            self.fig_3d.set_size_inches(
                self.width * 0.0104,
                (self.height - self.sizer.GetMinSize()[1] - 100) * 0.0104,
            )
            self.Update3DFrame()
        event.Skip()

    def OnSizeFrame(self, event):
        # Get the new frame size
        self.width, self.height = self.GetSize()
        self.SetSize((self.width, self.height))
        self.canvas_3d.SetSize(
            (
                self.width * 0.0104,
                (self.height - self.sizer.GetMinSize()[1] - 100) * 0.0104,
            )
        )
        self.fig_3d.set_size_inches(
            self.width * 0.0104,
            (self.height - self.sizer.GetMinSize()[1] - 100) * 0.0104,
        )
        self.Update3DFrame()
        event.Skip()

    def create_3D_plot_sizer(self):

        # Make a slider to change the contour levels
        self.contour_label = wx.StaticBox(self, -1, "Contour levels")
        self.contour_sizer = wx.StaticBoxSizer(self.contour_label, wx.VERTICAL)
        self.contour_slider = FloatSlider(
            self, id=-1, value=1, minval=0, maxval=3, res=0.01, style=wx.SL_HORIZONTAL
        )
        self.contour_slider.Bind(wx.EVT_SLIDER, self.OnContourSlider)
        self.contour_sizer.Add(self.contour_slider)
        self.sizer.Add(self.contour_sizer)
        self.main_3d_sizer.Add(self.sizer)

    def Update3DFrame(self):
        self.canvas_3d.draw()
        self.canvas_3d.Refresh()
        self.canvas_3d.Update()
        self.panel_3d.Refresh()
        self.panel_3d.Update()

    def plot3d(self):

        self.ax = self.fig_3d.add_subplot(111, projection="3d")

        contour_start = (
            np.max(self.main_frame.nmrdata.data) / 10
        )  # contour level start value
        self.contour_num = 20  # number of contour levels
        self.contour_factor = 1.2  # scaling factor between contour levels
        # calculate contour levels
        self.cl = contour_start * self.contour_factor ** np.arange(self.contour_num)
        self.cl_neg = -contour_start * self.contour_factor ** np.flip(
            np.arange(self.contour_num)
        )

        for i in range(self.main_frame.nmrdata.data.shape[0]):
            x = self.main_frame.ppms_1
            y = self.main_frame.ppms_0
            z = self.main_frame.nmrdata.data[i]
            x, y = np.meshgrid(x, y)
            self.ax.contour(
                x,
                y,
                z,
                zdir="z",
                offset=self.main_frame.ppms_2[i],
                levels=self.cl,
                colors="tab:orange",
                linewidths=0.5,
            )

        self.ax.set_zlim3d(
            np.min(self.main_frame.ppms_2), np.max(self.main_frame.ppms_2)
        )
        self.ax.set_xlim3d(
            np.min(self.main_frame.ppms_1), np.max(self.main_frame.ppms_1)
        )
        self.ax.set_ylim3d(
            np.min(self.main_frame.ppms_0), np.max(self.main_frame.ppms_0)
        )

        # Customize the plot by adding labels and adjusting the viewing angle
        self.ax.set_xlabel(self.main_frame.nmrdata.axislabels[2])
        self.ax.set_ylabel(self.main_frame.nmrdata.axislabels[1])
        self.ax.set_zlabel(self.main_frame.nmrdata.axislabels[0])

        self.ax.view_init(azim=200)

        self.ax.xaxis.pane.set_visible(False)
        self.ax.yaxis.pane.set_visible(False)
        self.ax.zaxis.pane.set_visible(False)

        # Remove grid lines
        self.ax.grid(False)

    def OnContourSlider(self, event):
        contour_val = 10 ** float(self.contour_slider.GetValue())
        contour_start = (
            np.max(self.main_frame.nmrdata.data) / contour_val
        )  # contour level start value
        self.contour_num = 20  # number of contour levels
        self.contour_factor = 1.2  # scaling factor between contour levels
        # calculate contour levels
        self.cl = contour_start * self.contour_factor ** np.arange(self.contour_num)
        self.cl_neg = -contour_start * self.contour_factor ** np.flip(
            np.arange(self.contour_num)
        )

        xlim = self.ax.get_xlim3d()
        ylim = self.ax.get_ylim3d()
        zlim = self.ax.get_zlim3d()
        xlabel = self.ax.get_xlabel()
        ylabel = self.ax.get_ylabel()
        zlabel = self.ax.get_zlabel()
        view = self.ax.azim

        self.ax.clear()

        for i in range(self.main_frame.nmrdata.data.shape[0]):
            x = self.main_frame.ppms_1
            y = self.main_frame.ppms_0
            z = self.main_frame.nmrdata.data[i]
            x, y = np.meshgrid(x, y)
            self.ax.contour(
                x,
                y,
                z,
                zdir="z",
                offset=self.main_frame.ppms_2[i],
                levels=self.cl,
                colors="red",
                linewidths=0.5,
            )

        self.ax.set_zlim3d(zlim[0], zlim[1])
        self.ax.set_xlim3d(xlim[0], xlim[1])
        self.ax.set_ylim3d(ylim[0], ylim[1])
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        self.ax.set_zlabel(zlabel)
        self.ax.view_init(azim=view)

        self.Update3DFrame()


class Stack2D(wx.Frame):
    def __init__(self, title, parent):
        self.main_frame = parent
        self.width = wx.GetDisplaySize()[0]
        try:
            if self.main_frame.parent.file_parser == True:
                os.chdir(self.main_frame.parent.path)
        except:
            pass
        if self.main_frame.parent.nmrdata.dim == 2:
            nmr_data_0 = GetData(file=self.main_frame.parent.nmrdata.file)
        else:
            projection_files = self.main_frame.parent.projection_files
            # Get current projection file
            file = projection_files[self.main_frame.parent.projection_selection_index]

            nmr_data_0 = ReadProjection(filename=file)
        self.nmr_data_old = nmr_data_0.data
        nmr_data_0.dim = 1

        if parent.transposed2D == True:
            nmr_data_0.data = nmr_data_0.data[0]
            nmr_data_0.axislabels = nmr_data_0.axislabels[0]
        else:
            nmr_data_0.data = nmr_data_0.data.T[0]
            nmr_data_0.axislabels = nmr_data_0.axislabels[1]
        # Get the monitor size and set the window size to 85% of the monitor size
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        self.display_index = wx.Display.GetFromWindow(parent)
        self.display_index_current = self.display_index
        self.width = 1.0 * sizes[self.display_index][0]
        self.height = 0.875 * sizes[self.display_index][1]
        wx.Frame.__init__(
            self, parent=parent, title=title, size=(int(self.width), int(self.height))
        )
        self.panel_stack = wx.Panel(self, -1)
        self.main_stack_sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_stack_sizer)

        try:
            dic, dat = ng.pipe.read(nmr_data_0.file)
        except:
            dic, dat = ng.pipe.read(nmr_data_0.filename)
        if parent.transposed2D == True:
            uc0 = ng.pipe.make_uc(dic, dat, dim=1)
        else:
            uc0 = ng.pipe.make_uc(dic, dat, dim=0)

        self.viewer_oneD = OneDViewer(parent=self, nmrdata=nmr_data_0, uc0=uc0)
        self.main_stack_sizer.Add(self.viewer_oneD, 1, wx.EXPAND)

        self.SetSizer(self.main_stack_sizer)

        # Make negative contour lines solid
        matplotlib.rc("contour", negative_linestyle="solid")

        self.viewer_oneD.files.stackmode = True
        if parent.transposed2D == True:
            self.viewer_oneD.files.transposed_stack = True
        self.viewer_oneD.files.nmrdata_original = parent.nmrdata
        try:
            self.viewer_oneD.files.OnDropFiles(0, 0, [nmr_data_0.file])
        except:
            self.viewer_oneD.files.OnDropFiles(0, 0, [nmr_data_0.filename])
        self.viewer_oneD.files.stackmode = False

        self.Show()
        self.Centre()

        try:
            if self.main_frame.parent.file_parser == True:
                os.chdir(self.main_frame.parent.cwd)
        except:
            pass

        # Bind method to check/resize the window when the frame is moved
        self.Bind(wx.EVT_MOVE, self.OnMoveFrame)

        # Bind method to resize the window when the frame is resized
        self.Bind(wx.EVT_SIZE, self.OnSizeFrame)

    def OnMoveFrame(self, event):
        # Get the new default display if the frame is moved
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        display_index = wx.Display.GetFromWindow(self)
        if display_index != self.display_index_current:
            self.display_index_current = display_index
            self.width = int(1.0 * sizes[display_index][0])
            self.height = int(0.875 * sizes[display_index][1])
            self.SetSize((self.width, self.height))
            self.viewer_oneD.canvas.SetSize(
                (
                    self.width * 0.0104,
                    (self.height - self.viewer_oneD.bottom_sizer.GetMinSize()[1] - 100)
                    * 0.0104,
                )
            )
            self.viewer_oneD.fig.set_size_inches(
                self.width * 0.0104,
                (self.height - self.viewer_oneD.bottom_sizer.GetMinSize()[1] - 100)
                * 0.0104,
            )
            self.viewer_oneD.UpdateFrame()
        event.Skip()

    def OnSizeFrame(self, event):
        # Get the new frame size
        self.width, self.height = self.GetSize()
        self.SetSize((self.width, self.height))
        self.viewer_oneD.canvas.SetSize(
            (
                self.width * 0.0104,
                (self.height - self.viewer_oneD.bottom_sizer.GetMinSize()[1] - 100)
                * 0.0104,
            )
        )
        self.viewer_oneD.fig.set_size_inches(
            self.width * 0.0104,
            (self.height - self.viewer_oneD.bottom_sizer.GetMinSize()[1] - 100)
            * 0.0104,
        )
        self.viewer_oneD.UpdateFrame()
        event.Skip()


class SpinBore(wx.Frame):
    def __init__(self, title, projection, parent=None):
        self.main_frame = parent
        # Get the monitor size and set the window size to 85% of the monitor size
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        self.display_index = wx.Display.GetFromWindow(parent)
        self.display_index_current = self.display_index
        self.width = int(1.0 * sizes[self.display_index][0])
        self.height = int(0.875 * sizes[self.display_index][1])
        wx.Frame.__init__(
            self, parent=parent, title=title, size=(self.width, self.height)
        )
        self.panel_bore = wx.Panel(self, -1)
        self.main_bore_sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_bore_sizer)

        self.fig_bore = Figure()
        self.fig_bore.tight_layout()
        self.canvas_bore = FigCanvas(self, -1, self.fig_bore)
        self.main_bore_sizer.Add(self.canvas_bore, 10, flag=wx.GROW)
        self.toolbar_bore = NavigationToolbar(self.canvas_bore)
        self.main_bore_sizer.Add(self.toolbar_bore, 0, wx.EXPAND)

        # Read the projection file
        self.nmrdata = ReadProjection(projection)

        self.sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.titlecolor = "black"

        self.make_bore_sizer()
        self.plot_bore_data()
        self.Show()
        self.Centre()

        # Bind method to check/resize the window when the frame is moved
        self.Bind(wx.EVT_MOVE, self.OnMoveFrame)

        # Bind method to resize the window when the frame is resized
        self.Bind(wx.EVT_SIZE, self.OnSizeFrame)

    def OnMoveFrame(self, event):
        # Get the new default display if the frame is moved
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        display_index = wx.Display.GetFromWindow(self)
        if display_index != self.display_index_current:
            self.display_index_current = display_index
            self.width = int(1.0 * sizes[display_index][0])
            self.height = int(0.875 * sizes[display_index][1])
            self.SetSize((self.width, self.height))
            self.canvas_bore.SetSize(
                (
                    self.width * 0.0104,
                    (self.height - self.bore_sizer.GetMinSize()[1] - 100) * 0.0104,
                )
            )
            self.fig_bore.set_size_inches(
                self.width * 0.0104,
                (self.height - self.bore_sizer.GetMinSize()[1] - 100) * 0.0104,
            )
            self.UpdateBoreFrame()
        event.Skip()

    def OnSizeFrame(self, event):
        # Get the new frame size
        self.width, self.height = self.GetSize()
        self.SetSize((self.width, self.height))
        self.canvas_bore.SetSize(
            (
                self.width * 0.0104,
                (self.height - self.bore_sizer.GetMinSize()[1] - 100) * 0.0104,
            )
        )
        self.fig_bore.set_size_inches(
            self.width * 0.0104,
            (self.height - self.bore_sizer.GetMinSize()[1] - 100) * 0.0104,
        )
        self.UpdateBoreFrame()
        event.Skip()

    def UpdateBoreFrame(self):
        # Updates the plots in the frame
        self.canvas_bore.draw()
        self.canvas_bore.Refresh()
        self.canvas_bore.Update()
        self.panel_bore.Refresh()
        self.panel_bore.Update()

    def make_bore_sizer(self):

        # Make sizer related to the 2D plot
        self.bore_sizer_2D_label = wx.StaticBox(self, -1, "2D Plot")
        self.bore_sizer_2D = wx.StaticBoxSizer(self.bore_sizer_2D_label, wx.HORIZONTAL)

        # Make a slider to change the contour levels
        self.bore_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.bore_contour_label = wx.StaticBox(self, -1, "Contour Max")
        self.bore_contour_sizer = wx.StaticBoxSizer(
            self.bore_contour_label, wx.VERTICAL
        )
        self.bore_slider = FloatSlider(
            self, id=-1, value=1, minval=0, maxval=3, res=0.01, style=wx.SL_HORIZONTAL
        )
        self.bore_slider.Bind(wx.EVT_SLIDER, self.OnBoreSlider)
        self.bore_contour_sizer.Add(self.bore_slider)
        self.bore_sizer_2D.Add(self.bore_contour_sizer)

        # Button to transpose the 2D data
        self.bore_transpose_button = wx.Button(self, -1, "Transpose")
        self.bore_transpose_button.Bind(wx.EVT_BUTTON, self.OnTransposeButtonBore)
        self.bore_sizer_2D.AddSpacer(10)
        self.bore_sizer_2D.Add(self.bore_transpose_button, 0, wx.ALIGN_CENTER_VERTICAL)

        # Sizer containing all 1D bore related items
        self.bore_sizer_1D_label = wx.StaticBox(self, -1, "1D Plot")
        self.bore_sizer_1D = wx.StaticBoxSizer(self.bore_sizer_1D_label, wx.HORIZONTAL)

        # Slider to change the scaling of the bore intensity
        self.bore_intensity_label = wx.StaticBox(self, -1, "Intensity")
        self.bore_intensity_sizer = wx.StaticBoxSizer(
            self.bore_intensity_label, wx.VERTICAL
        )
        self.bore_intensity_slider = FloatSlider(
            self, id=-1, value=1, minval=-1, maxval=10, res=0.01, style=wx.SL_HORIZONTAL
        )
        self.bore_intensity_slider.Bind(wx.EVT_SLIDER, self.OnIntensitySlider)
        self.bore_intensity_sizer.Add(self.bore_intensity_slider)
        self.bore_sizer_1D.Add(self.bore_intensity_sizer)

        self.bore_overlay_sizer_label = wx.StaticBox(self, -1, "Overlay")
        self.bore_overlay_sizer = wx.StaticBoxSizer(
            self.bore_overlay_sizer_label, wx.HORIZONTAL
        )

        # Toggle amino acid projections
        self.bore_toggle_button = wx.CheckBox(self, -1, "Show Amino Acid Predictions")
        self.bore_toggle_button.Bind(wx.EVT_CHECKBOX, self.OnToggleAminoAcid)

        # Have a combo box with 1H, 13C, 15N
        self.bore_combo_box = wx.ComboBox(
            self, -1, choices=["1H", "13C"], style=wx.CB_READONLY
        )
        self.bore_combo_box.Bind(wx.EVT_COMBOBOX, self.OnNucleusSelection)
        self.bore_overlay_sizer.Add(
            self.bore_toggle_button, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.bore_overlay_sizer.AddSpacer(10)
        self.bore_overlay_sizer.Add(self.bore_combo_box, 0, wx.ALIGN_CENTER_VERTICAL)

        # Have a combobox for free/protein
        self.bore_free_protein_combo_box = wx.ComboBox(
            self, -1, choices=["Free", "Protein"], style=wx.CB_READONLY
        )
        self.bore_free_protein_combo_box.Bind(
            wx.EVT_COMBOBOX, self.OnFreeProteinSelection
        )
        self.bore_overlay_sizer.AddSpacer(10)
        self.bore_overlay_sizer.Add(
            self.bore_free_protein_combo_box, 0, wx.ALIGN_CENTER_VERTICAL
        )

        # Combobox for amino acid selection
        self.bore_amino_acid_combo_box = wx.ComboBox(
            self,
            -1,
            choices=[
                "Alanine (A)",
                "Arginine (R)",
                "Asparagine (N)",
                "Aspartic Acid (D)",
                "Cysteine (C)",
                "Glutamic Acid (E)",
                "Glutamine (Q)",
                "Glycine (G)",
                "Histidine (H)",
                "Isoleucine (I)",
                "Leucine (L)",
                "Lysine (K)",
                "Methionine (M)",
                "Phenylalanine (F)",
                "Proline (P)",
                "Serine (S)",
                "Threonine (T)",
                "Tryptophan (W)",
                "Tyrosine (Y)",
                "Valine (V)",
            ],
            style=wx.CB_READONLY,
        )
        self.bore_amino_acid_combo_box.Bind(wx.EVT_COMBOBOX, self.OnAminoAcidSelection)
        self.bore_overlay_sizer.AddSpacer(10)
        self.bore_overlay_sizer.Add(
            self.bore_amino_acid_combo_box, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.include_overlay = False
        self.free_protein = "Free"

        self.bore_sizer_1D.AddSpacer(10)
        self.bore_sizer_1D.Add(self.bore_overlay_sizer)

        # Make a sizer for the strip plot
        self.bore_sizer_stripplot_label = wx.StaticBox(self, -1, "Strip Plot")
        self.bore_sizer_strip = wx.StaticBoxSizer(
            self.bore_sizer_stripplot_label, wx.HORIZONTAL
        )

        # Make a slider to change the contour levels
        self.bore_strip_contour_label = wx.StaticBox(self, -1, "Contour Max")
        self.bore_strip_contour_sizer = wx.StaticBoxSizer(
            self.bore_strip_contour_label, wx.VERTICAL
        )
        self.bore_strip_slider = FloatSlider(
            self, id=-1, value=1, minval=0, maxval=3, res=0.01, style=wx.SL_HORIZONTAL
        )
        self.bore_strip_slider.Bind(wx.EVT_SLIDER, self.OnBoreSliderStripPlot)
        self.bore_strip_contour_sizer.Add(self.bore_strip_slider)
        self.bore_sizer_strip.Add(self.bore_strip_contour_sizer)

        # Have a textcontrol for the strip width (default to 1ppm)
        self.bore_strip_width_label = wx.StaticBox(self, -1, "Strip Width (ppm)")
        self.bore_strip_width_sizer = wx.StaticBoxSizer(
            self.bore_strip_width_label, wx.VERTICAL
        )
        self.strip_width = 1.00
        self.bore_strip_width_text = wx.TextCtrl(
            self, -1, str(self.strip_width), style=wx.TE_PROCESS_ENTER
        )
        self.bore_strip_width_text.Bind(wx.EVT_TEXT_ENTER, self.OnStripWidthEnter)
        self.bore_strip_width_sizer.Add(self.bore_strip_width_text)
        self.bore_sizer_strip.AddSpacer(10)
        self.bore_sizer_strip.Add(self.bore_strip_width_sizer)

        # Button to swap the x axis for the strip plot
        self.bore_stripswap_button = wx.Button(self, -1, "Swap Axes")
        self.bore_stripswap_button.Bind(wx.EVT_BUTTON, self.OnStripSwapButtonBore)
        self.bore_sizer_strip.AddSpacer(10)
        self.bore_sizer_strip.Add(
            self.bore_stripswap_button, 0, wx.ALIGN_CENTER_VERTICAL
        )

        self.bore_sizer.Add(self.bore_sizer_2D)
        self.bore_sizer.AddSpacer(10)
        self.bore_sizer.Add(self.bore_sizer_1D)
        self.bore_sizer_row2 = wx.BoxSizer(wx.HORIZONTAL)
        self.bore_sizer_row2.AddSpacer(10)
        self.bore_sizer_row2.Add(self.bore_sizer_strip)

        self.main_bore_sizer.Add(self.bore_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.main_bore_sizer.Add(self.bore_sizer_row2, 0, wx.ALIGN_CENTER_HORIZONTAL)

    def OnStripSwapButtonBore(self, event):
        pass

    def OnStripWidthEnter(self, event):
        pass

    def plot_bore_data(self):
        # Make a figure containing 2 plots, one large 2D contour plot and a vertical smaller plot showing the bore down a selected 2D coordinate
        self.ax_bore, self.ax_bore_2, self.ax_bore_3 = self.fig_bore.subplots(
            1, 3, gridspec_kw={"width_ratios": [3, 1, 1]}
        )

        self.click_press_connect = self.fig_bore.canvas.mpl_connect(
            "button_press_event", self.on_click_bore
        )

        self.cmap = "#e41a1c"
        self.cmap_neg = "#377eb8"
        self.transposed2D = False

        contour_start = np.max(self.nmrdata.data) / 10  # contour level start value
        self.contour_num = 20  # number of contour levels
        self.contour_factor = 1.20  # scaling factor between contour levels
        # calculate contour levels
        self.cl = contour_start * self.contour_factor ** np.arange(self.contour_num)
        self.cl_neg = -contour_start * self.contour_factor ** np.flip(
            np.arange(self.contour_num)
        )

        # Get ppm values for x and y axis
        self.uc0 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=0)
        self.uc1 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=1)
        self.ppms_0 = self.uc0.ppm_scale()
        self.ppms_1 = self.uc1.ppm_scale()
        self.new_x_ppms = self.ppms_0
        self.new_y_ppms = self.ppms_1
        self.X, self.Y = np.meshgrid(self.ppms_1, self.ppms_0)
        self.contour1 = self.ax_bore.contour(
            self.Y, self.X, self.nmrdata.data, self.cl, colors=self.cmap, linewidths=0.5
        )
        self.contour1_neg = self.ax_bore.contour(
            self.Y,
            self.X,
            self.nmrdata.data,
            self.cl_neg,
            colors=self.cmap_neg,
            linewidths=0.5,
        )
        self.ax_bore.set_xlabel(self.nmrdata.axislabels[1])
        self.ax_bore.set_ylabel(self.nmrdata.axislabels[0])
        self.ax_bore.set_xlim(max(self.ppms_0), min(self.ppms_0))
        self.ax_bore.set_ylim(max(self.ppms_1), min(self.ppms_1))

        self.bore_initial = self.ppms_0[0], self.ppms_1[0]
        bore_initial_index = 0, 0
        # For each value in the bore data, find the intensity of the bore position
        self.bore_data = []
        for i in range(len(self.main_frame.ppms_2)):
            self.bore_data.append(
                self.main_frame.nmrdata.data[i][bore_initial_index[1]][
                    bore_initial_index[0]
                ]
            )

        # Plot the bore data
        self.ax_bore_2.plot(
            self.bore_data, self.main_frame.ppms_2, color="red", linewidth=0.5
        )
        self.ax_bore_2.set_ylim(
            max(self.main_frame.ppms_2), min(self.main_frame.ppms_2)
        )

        # Find the label of the 3rd dimension
        labels = self.main_frame.nmrdata.axislabels
        for i, label in enumerate(labels):
            if (
                label != self.nmrdata.axislabels[0]
                and label != self.nmrdata.axislabels[1]
            ):
                self.ax_bore_2.set_ylabel(label)

        (self.cross,) = self.ax_bore.plot(
            self.bore_initial[0], self.bore_initial[1], marker="X", color="k"
        )

        self.line = self.ax_bore_2.axhline(
            y=self.bore_initial[1], color="black", linewidth=0.5
        )
        self.line2 = self.ax_bore_2.axhline(
            y=self.bore_initial[0], color="black", linewidth=0.5
        )

        self.ax_bore_2.set_title("1D Bore")

        # Plot the strip plot contour plot
        contour_start_strip = (
            np.max(self.nmrdata.data) / 10
        )  # contour level start value
        self.contour_num_strip = 20  # number of contour levels
        self.contour_factor_strip = 1.20  # scaling factor between contour levels
        # calculate contour levels
        self.cl_strip = contour_start_strip * self.contour_factor_strip ** np.arange(
            self.contour_num_strip
        )
        self.cl_neg_strip = -contour_start_strip * self.contour_factor_strip ** np.flip(
            np.arange(self.contour_num_strip)
        )

        # # Get the bore data for the strip plot
        # self.bore_data_strip = []
        # # Get the x axis ppm values for the strip plot (bore_initial[0]+/-strip_width)
        # self.bore_initial0_indexes = np.where(np.abs(self.ppms_0 - self.bore_initial[0]) < self.strip_width)
        # self.bore_initial1_indexes = np.where(np.abs(self.ppms_1 - self.bore_initial[1]) < self.strip_width)

        self.bore_data_strip1 = []
        self.bore_data_strip2 = []
        for i in range(len(self.main_frame.ppms_2)):
            # Get the contour data for the strip plot
            self.bore_data_strip1.append(
                self.main_frame.nmrdata.data[i][bore_initial_index[1]]
            )
            # self.bore_data_strip2.append(self.main_frame.nmrdata.data[i][bore_initial_index[1]][self.bore_initial1_indexes[0][0]:self.bore_initial1_indexes[0][-1]])

        self.bore_data_strip1 = np.array(self.bore_data_strip1)

        # Get the ppm values for the strip plot
        self.ppms_2 = self.main_frame.ppms_2

        self.Xstrip, self.Ystrip = np.meshgrid(self.ppms_0, self.ppms_2)

        self.ax_bore_3.contour(
            self.Xstrip,
            self.Ystrip,
            self.bore_data_strip1,
            self.cl_strip,
            colors=self.cmap,
            linewidths=0.5,
        )
        self.ax_bore_3.contour(
            self.Xstrip,
            self.Ystrip,
            self.bore_data_strip1,
            self.cl_neg_strip,
            colors=self.cmap_neg,
            linewidths=0.5,
        )
        self.line3 = self.ax_bore_3.axvline(
            x=self.bore_initial[0], color="black", linewidth=0.5
        )
        self.ax_bore_3.set_xlabel(self.nmrdata.axislabels[1])
        self.ax_bore_3.set_xlim(max(self.ppms_0), min(self.ppms_0))
        self.ax_bore_3.set_ylim(max(self.ppms_2), min(self.ppms_2))

        self.ax_bore_3.set_title("Strip Plot")

        self.UpdateBoreFrame()

    def OnIntensitySlider(self, event):
        # Function to change the y axis limits
        intensity_percent = 10 ** float(self.bore_intensity_slider.GetValue())
        self.ax_bore_2.set_xlim(
            -(np.max(self.nmrdata.data) / 8) / (intensity_percent / 100),
            np.max(self.nmrdata.data) / (intensity_percent / 100),
        )
        self.UpdateBoreFrame()

    def on_click_bore(self, event):
        intensity_percent = 10 ** float(self.bore_intensity_slider.GetValue())
        if self.ax_bore_2.get_title() == "":
            title = ""
        else:
            title = self.ax_bore_2.get_title()
        if event.inaxes == self.ax_bore:
            # print(event.xdata, event.ydata)
            self.cross.set_xdata([event.xdata])
            self.cross.set_ydata([event.ydata])

            # Change the bore slice shown on the plot on the right
            if len(self.new_x_ppms) != len(self.main_frame.ppms_0):
                self.bore_initial = event.xdata, event.ydata
                self.bore_initial_index = np.argmin(
                    np.abs(self.main_frame.ppms_1 - self.bore_initial[0])
                ), np.argmin(np.abs(self.main_frame.ppms_0 - self.bore_initial[1]))
                self.bore_data = []
                for i in range(len(self.main_frame.ppms_2)):
                    self.bore_data.append(
                        self.main_frame.nmrdata.data[i][self.bore_initial_index[1]][
                            self.bore_initial_index[0]
                        ]
                    )
                self.bore_data = np.array(self.bore_data)
                ylabel = self.ax_bore_2.get_ylabel()
                self.ax_bore_2.clear()
                self.ax_bore_2.set_title(title)
                self.ax_bore_2.plot(
                    self.bore_data, self.main_frame.ppms_2, color="red", linewidth=0.5
                )
                self.ax_bore_2.set_ylim(
                    max(self.main_frame.ppms_2), min(self.main_frame.ppms_2)
                )
                self.ax_bore_2.set_xlim(
                    -(np.max(self.nmrdata.data) / 8) / (intensity_percent / 100),
                    np.max(self.nmrdata.data) / (intensity_percent / 100),
                )
                self.ax_bore_2.set_ylabel(ylabel)

                self.line1 = self.ax_bore_2.axhline(
                    y=event.xdata, color="black", linewidth=0.5
                )
                self.line2 = self.ax_bore_2.axhline(
                    y=event.ydata, color="black", linewidth=0.5
                )

                self.bore_data_strip1 = []
                for i in range(len(self.main_frame.ppms_2)):
                    # Get the contour data for the strip plot
                    self.bore_data_strip1.append(
                        self.main_frame.nmrdata.data[i][self.bore_initial_index[1]]
                    )

                self.bore_data_strip1 = np.array(self.bore_data_strip1)

                # Get the ppm values for the strip plot
                self.ppms_2 = self.main_frame.ppms_2

                self.Xstrip, self.Ystrip = np.meshgrid(self.ppms_0, self.ppms_2)

                title = self.ax_bore_3.get_title()
                xlim3, ylim3 = self.ax_bore_3.get_xlim(), self.ax_bore_3.get_ylim()
                self.ax_bore_3.clear()

                self.ax_bore_3.contour(
                    self.Xstrip,
                    self.Ystrip,
                    self.bore_data_strip1,
                    self.cl_strip,
                    colors=self.cmap,
                    linewidths=0.5,
                )
                self.ax_bore_3.contour(
                    self.Xstrip,
                    self.Ystrip,
                    self.bore_data_strip1,
                    self.cl_neg_strip,
                    colors=self.cmap_neg,
                    linewidths=0.5,
                )

                self.line3 = self.ax_bore_3.axvline(
                    x=self.bore_initial[0], color="black", linewidth=0.5
                )
                self.ax_bore_3.set_xlabel(self.nmrdata.axislabels[1])
                self.ax_bore_3.set_xlim(xlim3)
                self.ax_bore_3.set_title(title)
                self.ax_bore_3.set_ylim(ylim3)

            else:
                self.bore_initial = event.xdata, event.ydata
                self.bore_initial_index = np.argmin(
                    np.abs(self.main_frame.ppms_0 - self.bore_initial[0])
                ), np.argmin(np.abs(self.main_frame.ppms_1 - self.bore_initial[1]))
                self.bore_data = []
                for i in range(len(self.main_frame.ppms_2)):
                    self.bore_data.append(
                        self.main_frame.nmrdata.data[i][self.bore_initial_index[0]][
                            self.bore_initial_index[1]
                        ]
                    )
                self.bore_data = np.array(self.bore_data)
                ylabel = self.ax_bore_2.get_ylabel()
                self.ax_bore_2.clear()
                self.ax_bore_2.set_title(title)
                self.ax_bore_2.plot(
                    self.bore_data, self.main_frame.ppms_2, color="red", linewidth=0.5
                )
                self.ax_bore_2.set_ylim(
                    max(self.main_frame.ppms_2), min(self.main_frame.ppms_2)
                )
                self.ax_bore_2.set_xlim(
                    -(np.max(self.nmrdata.data) / 8) / (intensity_percent / 100),
                    np.max(self.nmrdata.data) / (intensity_percent / 100),
                )
                self.ax_bore_2.set_ylabel(ylabel)
                self.line1 = self.ax_bore_2.axhline(
                    y=event.xdata, color="black", linewidth=0.5
                )
                self.line2 = self.ax_bore_2.axhline(
                    y=event.ydata, color="black", linewidth=0.5
                )

                self.bore_data_strip1 = []
                for i in range(len(self.main_frame.ppms_2)):
                    # Get the contour data for the strip plot
                    self.bore_data_strip1.append(
                        self.main_frame.nmrdata.data[i][self.bore_initial_index[0]]
                    )

                self.bore_data_strip1 = np.array(self.bore_data_strip1)

                # Get the ppm values for the strip plot
                self.ppms_2 = self.main_frame.ppms_2

                self.Xstrip, self.Ystrip = np.meshgrid(self.ppms_0, self.ppms_2)

                title = self.ax_bore_3.get_title()
                xlim3, ylim3 = self.ax_bore_3.get_xlim(), self.ax_bore_3.get_ylim()
                self.ax_bore_3.clear()

                self.ax_bore_3.contour(
                    self.Xstrip,
                    self.Ystrip,
                    self.bore_data_strip1,
                    self.cl_strip,
                    colors=self.cmap,
                    linewidths=0.5,
                )
                self.ax_bore_3.contour(
                    self.Xstrip,
                    self.Ystrip,
                    self.bore_data_strip1,
                    self.cl_neg_strip,
                    colors=self.cmap_neg,
                    linewidths=0.5,
                )

                self.line3 = self.ax_bore_3.axvline(
                    x=self.bore_initial[1], color="black", linewidth=0.5
                )
                self.ax_bore_3.set_xlabel(self.nmrdata.axislabels[0])
                self.ax_bore_3.set_xlim(xlim3)
                self.ax_bore_3.set_title(title)
                self.ax_bore_3.set_ylim(ylim3)

            self.OverlayBore()

    def OverlayBore(self):
        if self.include_overlay == True:
            intensity_percent = 10 ** float(self.bore_intensity_slider.GetValue())
            # Plot current amino acid selection
            data = self.bmrb_free[self.amino_acid][self.nucleus][0]
            labels = self.bmrb_free[self.amino_acid][self.nucleus][1]

            xdata = self.bore_data
            ydata = self.main_frame.ppms_2
            line1_data = self.bore_initial[0]
            line2_data = self.bore_initial[1]
            ylabel = self.ax_bore_2.get_ylabel()
            self.ax_bore_2.clear()
            self.ax_bore_2.set_title(self.amino_acid)
            self.ax_bore_2.plot(
                self.bore_data, self.main_frame.ppms_2, color="red", linewidth=0.5
            )
            self.ax_bore_2.set_ylim(
                max(self.main_frame.ppms_2), min(self.main_frame.ppms_2)
            )
            self.ax_bore_2.set_xlim(
                -(np.max(self.nmrdata.data) / 8) / (intensity_percent / 100),
                np.max(self.nmrdata.data) / (intensity_percent / 100),
            )
            self.ax_bore_2.set_ylabel(ylabel)
            self.line1 = self.ax_bore_2.axhline(
                y=line1_data, color="black", linewidth=0.5
            )
            self.line2 = self.ax_bore_2.axhline(
                y=line2_data, color="black", linewidth=0.5
            )
            # Plot a horizontal line for each peak with a label
            for i in range(len(data)):
                self.ax_bore_2.axhline(y=data[i], color="black", linewidth=0.5)
                self.ax_bore_2.text(
                    np.max(self.bore_data),
                    data[i],
                    labels[i],
                    color="black",
                    fontsize=6,
                )

        self.UpdateBoreFrame()

    def OnBoreSliderStripPlot(self, event):
        self.x_val3 = 10 ** float(self.bore_strip_slider.GetValue())
        # update contour levels for strip plot
        self.contour_start_strip = np.max(np.abs(self.nmrdata.data)) / self.x_val3
        self.cl_strip = (
            self.contour_start_strip
            * self.contour_factor_strip ** np.arange(self.contour_num_strip)
        )
        self.cl_neg_strip = (
            -self.contour_start_strip
            * self.contour_factor_strip ** np.flip(np.arange(self.contour_num_strip))
        )
        xlim3, ylim3 = self.ax_bore_3.get_xlim(), self.ax_bore_3.get_ylim()
        xlabel = self.ax_bore_3.get_xlabel()
        title = self.ax_bore_3.get_title()
        self.ax_bore_3.clear()
        self.contour1 = self.ax_bore_3.contour(
            self.Xstrip,
            self.Ystrip,
            self.bore_data_strip1,
            self.cl_strip,
            colors=self.cmap,
            linewidths=0.5,
        )
        self.contour1_neg = self.ax_bore_3.contour(
            self.Xstrip,
            self.Ystrip,
            self.bore_data_strip1,
            self.cl_neg_strip,
            colors=self.cmap_neg,
            linewidths=0.5,
        )
        self.line3 = self.ax_bore_3.axvline(
            x=self.bore_initial[1], color="black", linewidth=0.5
        )
        self.ax_bore_3.set_xlim(xlim3)
        self.ax_bore_3.set_ylim(ylim3)
        self.ax_bore_3.set_xlabel(self.nmrdata.axislabels[1])
        self.ax_bore_3.set_title(title)
        self.UpdateBoreFrame()

    def OnBoreSlider(self, event):
        # Function to update the contour levels when the user changes the number of contour levels
        self.x_val = 10 ** float(self.bore_slider.GetValue())

        # update contour levels
        self.contour_start = np.max(np.abs(self.nmrdata.data)) / self.x_val
        self.cl = self.contour_start * self.contour_factor ** np.arange(
            self.contour_num
        )
        self.cl_neg = -self.contour_start * self.contour_factor ** np.flip(
            np.arange(self.contour_num)
        )

        xlim, ylim = self.ax_bore.get_xlim(), self.ax_bore.get_ylim()
        cross_data = self.cross.get_data()
        self.ax_bore.clear()
        self.contour1 = self.ax_bore.contour(
            self.Y, self.X, self.nmrdata.data, self.cl, colors=self.cmap, linewidths=0.5
        )
        self.contour1_neg = self.ax_bore.contour(
            self.Y,
            self.X,
            self.nmrdata.data,
            self.cl_neg,
            colors=self.cmap_neg,
            linewidths=0.5,
        )
        (self.cross,) = self.ax_bore.plot(
            cross_data[0], cross_data[1], marker="X", color="k"
        )

        self.ax_bore.set_xlim(xlim)
        self.ax_bore.set_ylim(ylim)
        self.ax_bore.set_xlabel(self.nmrdata.axislabels[1])
        self.ax_bore.set_ylabel(self.nmrdata.axislabels[0])

        self.UpdateBoreFrame()

    def OnTransposeButtonBore(self, event):
        if self.transposed2D == False:
            self.transposed2D = True
        else:
            self.transposed2D = False
        xlim_old, ylim_old = self.ax_bore.get_xlim(), self.ax_bore.get_ylim()
        self.X_old, self.Y_old = self.X, self.Y
        self.new_x_ppms_old = self.new_x_ppms
        self.new_y_ppms_old = self.new_y_ppms
        self.new_x_ppms = self.new_y_ppms_old
        self.new_y_ppms = self.new_x_ppms_old
        self.X, self.Y = np.meshgrid(self.new_y_ppms, self.new_x_ppms)
        self.nmr_data_old = self.nmrdata.data
        self.nmrdata.data = self.nmr_data_old.T
        cross_data = self.cross.get_data()
        self.ax_bore.clear()
        self.contour1 = self.ax_bore.contour(
            self.Y, self.X, self.nmrdata.data, self.cl, colors=self.cmap, linewidths=0.5
        )
        self.contour1_neg = self.ax_bore.contour(
            self.Y,
            self.X,
            self.nmrdata.data,
            self.cl_neg,
            colors=self.cmap_neg,
            linewidths=0.5,
        )
        self.ax_bore.set_xlim([max(self.new_x_ppms), min(self.new_x_ppms)])
        self.ax_bore.set_ylim([max(self.new_y_ppms), min(self.new_y_ppms)])
        self.axislabels_old = self.nmrdata.axislabels[0], self.nmrdata.axislabels[1]
        self.nmrdata.axislabels[1] = self.axislabels_old[0]
        self.nmrdata.axislabels[0] = self.axislabels_old[1]
        uc0, uc1 = self.uc0, self.uc1
        self.uc0 = uc1
        self.uc1 = uc0
        self.ax_bore.set_xlabel(self.nmrdata.axislabels[1])
        self.ax_bore.set_ylabel(self.nmrdata.axislabels[0])
        (self.cross,) = self.ax_bore.plot(
            cross_data[1], cross_data[0], marker="X", color="k"
        )
        self.OnBoreSlider(wx.EVT_SCROLL)
        self.toolbar_bore.update()

    def OnToggleAminoAcid(self, event):
        # Read in the amino acid/BMRB protein statistics
        if self.bore_toggle_button.GetValue() == True:
            self.include_overlay = True
            self.read_bmrb()
        else:
            self.include_overlay = False

        self.OnAminoAcidSelection(event)

    def OnAminoAcidSelection(self, event):
        # Get the amino acid selection
        amino_acid = self.bore_amino_acid_combo_box.GetValue()
        if amino_acid == "Alanine (A)":
            amino_acid = "ALA"
        elif amino_acid == "Arginine (R)":
            amino_acid = "ARG"
        elif amino_acid == "Asparagine (N)":
            amino_acid = "ASN"
        elif amino_acid == "Aspartic Acid (D)":
            amino_acid = "ASP"
        elif amino_acid == "Cysteine (C)":
            amino_acid = "CYS"
        elif amino_acid == "Glutamic Acid (E)":
            amino_acid = "GLU"
        elif amino_acid == "Glutamine (Q)":
            amino_acid = "GLN"
        elif amino_acid == "Glycine (G)":
            amino_acid = "GLY"
        elif amino_acid == "Histidine (H)":
            amino_acid = "HIS"
        elif amino_acid == "Isoleucine (I)":
            amino_acid = "ILE"
        elif amino_acid == "Leucine (L)":
            amino_acid = "LEU"
        elif amino_acid == "Lysine (K)":
            amino_acid = "LYS"
        elif amino_acid == "Methionine (M)":
            amino_acid = "MET"
        elif amino_acid == "Phenylalanine (F)":
            amino_acid = "PHE"
        elif amino_acid == "Proline (P)":
            amino_acid = "PRO"
        elif amino_acid == "Serine (S)":
            amino_acid = "SER"
        elif amino_acid == "Threonine (T)":
            amino_acid = "THR"
        elif amino_acid == "Tryptophan (W)":
            amino_acid = "TRP"
        elif amino_acid == "Tyrosine (Y)":
            amino_acid = "TYR"
        elif amino_acid == "Valine (V)":
            amino_acid = "VAL"

        self.amino_acid = amino_acid
        # Get the nucleus selection
        self.OnNucleusSelection(event)

    def OnNucleusSelection(self, event):
        if self.bore_combo_box.GetSelection() == 0:
            self.nucleus = "H"
        else:
            self.nucleus = "C"
        self.OnFreeProteinSelection(event)

    def OnFreeProteinSelection(self, event):
        if self.bore_free_protein_combo_box.GetSelection == 0:
            self.free_protein = "Free"
        else:
            self.free_protein = "Protein"

        self.OverlayBore()

    def read_bmrb(self):
        directory = __file__.split("SpinView.py")[0]
        file1 = directory + "bmrb_free.txt"  # Free amino acid chemical shift data
        file2 = directory + "bmrb.txt"  # Protein amino acid chemical shift data

        self.read_free()
        # else:
        #     # Give a warning saying that the 'bmrb_free.txt' file is missing
        #     message = "The file 'bmrb_free.txt' is missing from the SpinView directory. Unable to show free amino acid overlays."
        #     dlg = wx.MessageDialog(None, message, 'File Missing', wx.OK | wx.ICON_WARNING)
        #     dlg.ShowModal()
        #     dlg.Destroy()
        #     self.include_overlay = False
        #     return

    def read_free(self):
        self.bmrb_free = {}
        data = """ALA CO 178.56
ALA CA 53.2
ALA CB 18.83
ALA HA 3.771
ALA HB 1.471

ARG CO 177.238
ARG CA 57.002
ARG CB 30.281
ARG CG 26.577
ARG CD 43.201
ARG CZ 159.504
ARG HA 3.764
ARG HB 1.909
ARG HG 1.679
ARG HD 3.236

ASP CO 177
ASP CA 54.91
ASP CB 180.3
ASP HA 3.889
ASP HB1 2.786
ASP HB2 2.703

ASN CO 177.173
ASN CA 53.990
ASN CB 37.278
ASN CG 176.206
ASN HA 3.991
ASN HB1 2.940
ASN HB2 2.840

CYS CO 175.486
CYS CA 58.680
CYS CB 27.647
CYS HA 3.952
CYS HB 3.044

GLU CO 177.360
GLU CA 57.357
GLU CB 29.728
GLU CG 36.20
GLU CD 184.088
GLU HA 3.747
GLU HB 2.078
GLU HG 2.339

GLN CO 176.83
GLN CA 56.83
GLN CB 28.95
GLN CG 33.52
GLN CD 180.37
GLN HA 3.764
GLN HB 2.13
GLN HG 2.447

GLY CO 175.225
GLY CA 44.133
GLY HA 3.545

HIS CO 176.642
HIS CA 57.435
HIS CB 30.696
HIS CG 134.45
HIS HA 3.98
HIS HB1 3.234
HIS HB2 3.131

ILE CO 176.972
ILE CA 62.249
ILE CB 38.614
ILE CG1 17.411
ILE CG2 27.174
ILE CD 13.834
ILE HA 3.657
ILE HB 1.969
ILE HG1 0.998
ILE HG2a 1.458
ILE HG2b 1.249
ILE HD 0.927

LEU CO 178.382
LEU CA 56.112
LEU CB 42.526
LEU CG 26.871
LEU CD1 24.751
LEU CD2 23.589
LEU HA 3.719
LEU HB/G 1.701
LEU HD 0.949

LYS CO 177.484
LYS CA 57.190
LYS CB 32.628
LYS CG 24.145
LYS CD 29.137
LYS CE 41.754
LYS HA 3.754
LYS HB 1.895
LYS HG 1.465
LYS HD 1.716
LYS HE 3.012

MET CO 177.093
MET CA 56.584
MET CB 32.395
MET CG 31.513
MET CD 16.636
MET HA 3.850
MET HB1 2.183
MET HB2 2.122
MET HG 2.629
MET HD 2.122

PHE CO 176.774
PHE CA 58.744
PHE CB 39.095
PHE HA 3.975
PHE HB1 3.271
PHE HB2 3.110

PRO CO 177.483
PRO CA 63.922
PRO CB 31.728
PRO HA 4.119
PRO HB1 2.337
PRO HB2 2.022
PRO HG 2.022
PRO HD 3.366

SER CO 175.227
SER CA 59.096
SER CB 62.906
SER HA 3.828
SER HB 3.952

THR CO 175.689
THR CA 63.172
THR CB 68.679
THR CG 22.179
THR HA 3.573
THR HB 4.241
THR HG 1.318

TRP CO 177.332
TRP CA 57.764
TRP CB 29.152
TRP HA 4.036
TRP HB1 3.471
TRP HB2 3.292

TYR CO 176.964
TYR CA 58.838
TYR CB 38.277
TYR HA 3.936
TYR HB1 3.200
TYR HB2 3.055

VAL CO 177.086
VAL CA 63.083
VAL CB 31.834
VAL CG1 20.696
VAL CG2 19.368
VAL HA 3.599
VAL HB 2.266
VAL HG1 1.034
VAL HG2 0.981"""

        for line in data.splitlines():
            line = line.split("\n")[0].split()
            if len(line) != 3:
                continue
            if line[0] not in self.bmrb_free.keys():
                self.bmrb_free[line[0]] = {}
            if "H" not in self.bmrb_free[line[0]].keys():
                self.bmrb_free[line[0]]["H"] = []
            if "C" not in self.bmrb_free[line[0]].keys():
                self.bmrb_free[line[0]]["C"] = []
            if line[1][0] == "H":
                self.bmrb_free[line[0]]["H"].append([line[1], float(line[2])])
            elif line[1][0] == "C":
                self.bmrb_free[line[0]]["C"].append([line[1], float(line[2])])

        for key in self.bmrb_free.keys():
            H_values = []
            H_labels = []
            C_values = []
            C_labels = []
            for i in range(len(self.bmrb_free[key]["H"])):
                H_values.append(self.bmrb_free[key]["H"][i][1])
                H_labels.append(self.bmrb_free[key]["H"][i][0])
            for i in range(len(self.bmrb_free[key]["C"])):
                C_values.append(self.bmrb_free[key]["C"][i][1])
                C_labels.append(self.bmrb_free[key]["C"][i][0])

            self.bmrb_free[key]["H"] = [H_values, H_labels]
            self.bmrb_free[key]["C"] = [C_values, C_labels]


class uSTA_Dialog(wx.Dialog):
    def __init__(self, title, parent):
        self.main_frame = parent
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        width = 450
        height = 150
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.panel_uSTAparams = wx.Panel(self, -1)
        self.main_uSTAparams = wx.BoxSizer(wx.HORIZONTAL)
        self.SetSizer(self.main_uSTAparams)

        self.make_uSTAparams_sizer()
        self.Show()

    def make_uSTAparams_sizer(self):
        # Make a sizer to hold the text box and button
        self.uSTAparams_sizer = wx.BoxSizer(wx.VERTICAL)
        self.uSTAparams_sizer.AddSpacer(30)

        self.uSTAparams_sizer2 = wx.BoxSizer(wx.HORIZONTAL)
        self.uSTAparams_sizer2.AddSpacer(5)

        self.mixing_time_label = wx.StaticText(self, -1, "Mixing time (s):")
        self.uSTAparams_sizer2.Add(self.mixing_time_label)

        self.uSTA_mixing_time = wx.TextCtrl(
            self, -1, value=self.main_frame.mixing_time, size=(100, 20)
        )
        self.uSTAparams_sizer2.AddSpacer(10)
        self.uSTAparams_sizer2.Add(self.uSTA_mixing_time)
        self.uSTAparams_sizer2.AddSpacer(10)

        self.power_level_label = wx.StaticText(self, -1, "Power level:")

        self.uSTAparams_sizer2.Add(self.power_level_label)
        self.uSTA_power_level = wx.TextCtrl(
            self, -1, value=self.main_frame.power_level, size=(100, 20)
        )
        self.uSTAparams_sizer2.AddSpacer(10)
        self.uSTAparams_sizer2.Add(self.uSTA_power_level)
        self.uSTAparams_sizer2.AddSpacer(10)

        self.uSTAparams_sizer.Add(self.uSTAparams_sizer2)

        self.uSTAparams_sizer.AddSpacer(30)

        # Have a button to confirm the selection
        self.confirm_button = wx.Button(self, -1, "Confirm")
        self.confirm_button.Bind(wx.EVT_BUTTON, self.OnConfirm)
        self.uSTAparams_sizer.Add(self.confirm_button, 0, wx.ALIGN_CENTER | wx.ALL, 10)
        self.uSTAparams_sizer.AddSpacer(10)

        self.main_uSTAparams.AddSpacer(10)
        self.main_uSTAparams.Add(self.uSTAparams_sizer)
        self.main_uSTAparams.AddSpacer(10)

    def OnConfirm(self, event):

        # Getting the ppm values
        ppm_values = self.main_frame.ppms_0
        if len(ppm_values) == 2:
            ppm_values = self.main_frame.ppms_1

        # Getting the current directory name
        current_dir_name = os.path.basename(os.getcwd())
        # Getting the intensities of on and on resonance spectra
        usta_data = self.main_frame.nmrdata.data
        if len(usta_data) > 2:
            usta_data = self.main_frame.nmrdata.data.T

        on_data = usta_data[0]
        off_data = usta_data[1]

        data_file_name = current_dir_name + ".data"
        data_file = open(data_file_name, "w")
        for i in range(len(ppm_values)):
            data_file.write(
                current_dir_name
                + "\t"
                + self.uSTA_mixing_time.GetValue()
                + "\t"
                + self.uSTA_power_level.GetValue()
                + "\t"
                + str(ppm_values[i])
                + "\t"
                + str(on_data[i])
                + "\t"
                + str(off_data[i])
                + "\n"
            )
        data_file.close()

        self.Destroy()


class CESTOrder_Dialog(wx.Dialog):
    def __init__(self, title, parent):
        self.main_frame = parent
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        width = 300
        height = 150
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.panel_CESTOrder = wx.Panel(self, -1)
        self.main_CESTOrder = wx.BoxSizer(wx.HORIZONTAL)
        self.SetSizer(self.main_CESTOrder)

        self.make_CESTOrder_sizer()
        self.Show()

    def make_CESTOrder_sizer(self):
        # Make a sizer to hold the text box and button
        self.CESTOrder_sizer = wx.BoxSizer(wx.VERTICAL)
        self.CESTOrder_sizer.AddSpacer(5)

        self.CESTOrder_label = wx.StaticText(
            self, -1, "Interleaved order:", style=wx.ALIGN_LEFT
        )

        choices = ["OnResonance, OffResonance", "OffResonance, OnResonance"]
        self.CESTOrder_radiobox = wx.RadioBox(
            self, -1, choices=choices, style=wx.CB_READONLY | wx.RA_SPECIFY_ROWS
        )
        self.CESTOrder_sizer.AddSpacer(10)
        self.CESTOrder_sizer.Add(self.CESTOrder_radiobox, wx.ALIGN_LEFT)
        self.CESTOrder_sizer.AddSpacer(10)
        # Have a button to confirm the selection
        self.confirm_button = wx.Button(self, -1, "Confirm")
        self.confirm_button.Bind(wx.EVT_BUTTON, self.OnConfirm)
        self.CESTOrder_sizer.Add(self.confirm_button, wx.ALIGN_CENTER)
        self.CESTOrder_sizer.AddSpacer(10)

        self.main_CESTOrder.AddSpacer(10)
        self.main_CESTOrder.Add(self.CESTOrder_sizer, wx.ALIGN_CENTER)
        self.main_CESTOrder.AddSpacer(10)

    def OnConfirm(self, event):
        self.main_frame.CESTArrayOrder = self.CESTOrder_radiobox.GetSelection()
        self.main_frame.continue_deletion()
        self.Destroy()


class CESTFrame(wx.Frame):
    def __init__(self, title, parent=None, CESTArrayOrder=0):
        self.main_frame = parent
        # Get the monitor size and set the window size to 85% of the monitor size
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        self.display_index = wx.Display.GetFromWindow(parent)
        self.display_index_current = self.display_index
        self.width = int(1.0 * sizes[self.display_index][0])
        self.height = int(0.875 * sizes[self.display_index][1])
        wx.Frame.__init__(
            self, parent=parent, title=title, size=(self.width, self.height)
        )
        self.panel_CEST = wx.Panel(self, -1)
        self.main_CEST_sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_CEST_sizer)

        self.CESTArrayOrder = CESTArrayOrder

        self.fig_CEST = Figure()
        self.fig_CEST.tight_layout()
        self.canvas_CEST = FigCanvas(self, -1, self.fig_CEST)
        self.main_CEST_sizer.Add(self.canvas_CEST, 1, flag=wx.EXPAND | wx.ALL)
        self.toolbar_CEST = NavigationToolbar(self.canvas_CEST)
        self.main_CEST_sizer.Add(self.toolbar_CEST, 0, wx.EXPAND)

        self.make_CEST_sizer()

        self.titlecolor = "black"

        self.find_CEST_frequencies()
        self.organise_CEST_data()

        self.plot_CEST_data()
        self.Show()

        # Bind method to check/resize the window when the frame is moved
        self.Bind(wx.EVT_MOVE, self.OnMoveFrame)

        # Bind method to resize the window when the frame is resized
        self.Bind(wx.EVT_SIZE, self.OnSizeFrame)

    def find_CEST_frequencies(self):
        # Try to find procpar file in current directory, if can find it, work out frequencies based on tof_sel
        # If cannot find it will just have to the frequency indexes of 0 to n
        try:
            dic_v, data_v = ng.varian.read("./")

            # find the CEST offset values used
            offsets = dic_v["procpar"]["tof_sel"]["values"]
            self.offsets_Hz = []
            for i in range(len(offsets)):
                self.offsets_Hz.append(float(offsets[i]))

            # get ppm values for the offsets used
            self.tof = float(dic_v["procpar"]["tof"]["values"][0])

            # get the sfrq
            self.sfrq = float(dic_v["procpar"]["sfrq"]["values"][0])

            # open fid.com to find the carrier frequency
            with open("fid.com", "r") as file:
                lines = file.readlines()
                for line in lines:
                    if "-xCAR" in line:
                        self.carrier = float(line.split()[1])

            self.offsets_ppm = []
            for tof_sel in self.offsets_Hz:
                offset_ppm = (
                    tof_sel - self.tof
                ) / self.sfrq + self.carrier  # get actual ppm values
                self.offsets_ppm.append(offset_ppm)

        except:
            self.offsets_ppm = np.arange(0, len(self.main_frame.ppms_0), 1)

    def organise_CEST_data(self):
        # Get the CEST data from the main frame
        self.CEST_data = self.main_frame.nmrdata.data.T
        self.cest_on_data = []
        self.cest_off_data = []
        for i in range(len(self.CEST_data)):
            if i % 2 == 0:
                if self.CESTArrayOrder == 0:
                    self.cest_on_data.append(self.CEST_data[i])
                else:
                    self.cest_off_data.append(self.CEST_data[i])
            else:
                if self.CESTArrayOrder == 0:
                    self.cest_off_data.append(self.CEST_data[i])
                else:
                    self.cest_on_data.append(self.CEST_data[i])

        # Find the selected 1H chemical shift range in the main frame
        self.selected_shift = self.main_frame.line4.get_xdata()[0]
        self.selected_shift_index = np.argmin(
            np.abs(self.main_frame.ppms_1 - self.selected_shift)
        )

        self.non_normalized_cest_data = []
        self.normalized_cest_data = []
        for i in range(len(self.cest_on_data)):
            self.non_normalized_cest_data.append(
                self.cest_on_data[i][self.selected_shift_index]
            )
            self.normalized_cest_data.append(
                self.cest_on_data[i][self.selected_shift_index]
                / self.cest_off_data[i][self.selected_shift_index]
            )

    def OnMoveFrame(self, event):
        # Get the new default display if the frame is moved
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        display_index = wx.Display.GetFromWindow(self)
        if display_index != self.display_index_current:
            self.display_index_current = display_index
            self.width = int(1.0 * sizes[display_index][0])
            self.height = int(0.875 * sizes[display_index][1])
            self.SetSize((self.width, self.height))
            self.canvas_CEST.SetSize(
                (
                    self.width * 0.0104,
                    (self.height - self.CEST_sizer.GetMinSize()[1] - 100) * 0.0104,
                )
            )
            self.fig_CEST.set_size_inches(
                self.width * 0.0104,
                (self.height - self.CEST_sizer.GetMinSize()[1] - 100) * 0.0104,
            )
            self.UpdateCESTFrame()
        event.Skip()

    def OnSizeFrame(self, event):
        # Get the new frame size
        self.width, self.height = self.GetSize()
        self.SetSize((self.width, self.height))
        self.canvas_CEST.SetSize(
            (
                self.width * 0.0104,
                (self.height - self.CEST_sizer.GetMinSize()[1] - 100) * 0.0104,
            )
        )
        self.fig_CEST.set_size_inches(
            self.width * 0.0104,
            (self.height - self.CEST_sizer.GetMinSize()[1] - 100) * 0.0104,
        )
        self.UpdateCESTFrame()
        event.Skip()

    def UpdateCESTFrame(self):
        self.canvas_CEST.draw()
        self.canvas_CEST.Refresh()
        self.canvas_CEST.Update()
        self.panel_CEST.Refresh()
        self.panel_CEST.Update()

    def make_CEST_sizer(self):
        self.CEST_sizer = wx.BoxSizer(wx.HORIZONTAL)

        # Create a button that opens a file for a user to input the delay times
        self.CEST_ppm_range_label = wx.StaticBox(
            self, -1, "3D plot chemical shift range"
        )
        self.CEST_ppm_sizer_total = wx.StaticBoxSizer(
            self.CEST_ppm_range_label, wx.VERTICAL
        )
        self.CEST_ppm_sizer = wx.BoxSizer(wx.HORIZONTAL)

        # ppm min value
        self.CEST_ppm_min_label = wx.StaticText(self, -1, "Min ppm:")
        self.min_val = min(self.main_frame.ppms_1)
        self.CEST_ppm_min_text = wx.TextCtrl(
            self, -1, str(self.min_val), style=wx.TE_PROCESS_ENTER
        )
        self.CEST_ppm_min_text.Bind(wx.EVT_TEXT_ENTER, self.OnCEST_ppm_change)
        self.CEST_ppm_sizer.Add(self.CEST_ppm_min_label, wx.ALIGN_CENTER)
        self.CEST_ppm_sizer.AddSpacer(5)
        self.CEST_ppm_sizer.Add(self.CEST_ppm_min_text)

        # ppm max value
        self.CEST_ppm_max_label = wx.StaticText(self, -1, "Max ppm:")
        self.max_val = max(self.main_frame.ppms_1)
        self.CEST_ppm_max_text = wx.TextCtrl(
            self, -1, str(self.max_val), style=wx.TE_PROCESS_ENTER
        )
        self.CEST_ppm_max_text.Bind(wx.EVT_TEXT_ENTER, self.OnCEST_ppm_change)
        self.CEST_ppm_sizer.AddSpacer(10)
        self.CEST_ppm_sizer.Add(self.CEST_ppm_max_label)
        self.CEST_ppm_sizer.AddSpacer(5)
        self.CEST_ppm_sizer.Add(self.CEST_ppm_max_text)

        self.CEST_ppm_sizer_total.Add(self.CEST_ppm_sizer)
        self.CEST_sizer.Add(self.CEST_ppm_sizer_total, wx.ALIGN_CENTER_HORIZONTAL)

        # Make a button which will save the normalised CEST data as a spectrum
        self.save_CEST_button = wx.Button(self, -1, "Save normalised CEST data")
        self.save_CEST_button.Bind(wx.EVT_BUTTON, self.OnSaveCESTData)
        self.CEST_sizer.AddSpacer(10)
        self.CEST_sizer.Add(self.save_CEST_button, wx.ALIGN_CENTER_HORIZONTAL)

        self.main_CEST_sizer.AddSpacer(10)
        self.main_CEST_sizer.Add(self.CEST_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.main_CEST_sizer.AddSpacer(10)

    def OnSaveCESTData(self, event):
        # The ppms are self.offsets_ppm
        # The data is self.normalized_cest_data
        # Save the data as an nmrPipe spectrum using nmrglue

        try:
            # Ask the user what the file name should be
            file_name = ""
            message = "Input the file name to save the CEST data as"
            dlg = wx.TextEntryDialog(None, message, "Save CEST data", "CEST_data.ft")
            if dlg.ShowModal() == wx.ID_OK:
                file_name = dlg.GetValue()
                dlg.Destroy()

            if file_name == "":
                # Give an error message to say that the user must input a file name
                message = "Error: No file name given"
                dlg = wx.MessageDialog(None, message, "Error", wx.OK | wx.ICON_ERROR)
                dlg.ShowModal()
                dlg.Destroy()
                return

            # Create the nmrPipe file
            data = np.flip(np.array(self.normalized_cest_data) * 100)
            data = data.astype(np.float32)

            obs = self.sfrq
            sw = max(self.offsets_Hz) - min(self.offsets_Hz)
            car = self.carrier
            size = len(data)
            label = "ppm"
            orig = car * obs - sw * (size / 2 - 1) / size
            center = self.main_frame.nmrdata.dic["FDF2CENTER"]
            udic = {
                "ndim": 1,
                0: {
                    "size": size,
                    "complex": False,
                    "encoding": "int",
                    "sw": sw,
                    "obs": obs,
                    "car": car,
                    "label": label,
                    "time": False,
                    "freq": True,
                },
            }

            dic = ng.pipe.create_dic(udic)

            dic["FDF2LABEL"] = label
            dic["FDF2OBS"] = obs
            dic["FDF2SW"] = sw
            dic["FDF2CAR"] = car
            dic["FDF2SIZE"] = size
            dic["FDF2ORIG"] = orig
            dic["FDF2CENTER"] = center
            ng.pipe.write(file_name, dic, data, overwrite=True)

            message = "CEST data saved as 'CEST_data.ft'"
            dlg = wx.MessageDialog(
                None, message, "Save successful", wx.OK | wx.ICON_INFORMATION
            )
            dlg.ShowModal()
            dlg.Destroy()

        except:
            message = "Error saving CEST data"
            dlg = wx.MessageDialog(None, message, "Error", wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()

    def plot_CEST_data(self):

        self.ax = self.fig_CEST.add_subplot(121, projection="3d")

        for i, y in enumerate(self.offsets_ppm):
            self.ax.plot(
                self.main_frame.ppms_1,
                np.full_like(self.main_frame.ppms_1, y),
                self.cest_on_data[i],
                color="tab:grey",
                linewidth=1.5,
                alpha=0.5,
            )
            self.ax.scatter(
                self.main_frame.ppms_1[self.selected_shift_index],
                y,
                self.non_normalized_cest_data[i],
                color="red",
                s=10,
            )

        self.ax.set_title("CEST data")
        self.ax.set_xlabel(self.main_frame.nmrdata.axislabels[1])
        self.ax.set_ylabel(self.main_frame.nmrdata.axislabels[0])

        self.ax.set_xlim([max(self.main_frame.ppms_1), min(self.main_frame.ppms_1)])
        self.ax.set_ylim([max(self.offsets_ppm), min(self.offsets_ppm)])
        self.ax.view_init(elev=10.0, azim=-45)

        self.ax2 = self.fig_CEST.add_subplot(122)
        self.ax2.plot(
            self.offsets_ppm, self.normalized_cest_data, color="tab:red", linewidth=1.5
        )
        self.ax2.set_title("Normalized CEST data")
        self.ax2.set_xlabel("CEST offset (ppm)")
        self.ax2.set_ylabel("Normalized CEST data")
        self.ax2.set_xlim([max(self.offsets_ppm), min(self.offsets_ppm)])
        self.ax2.set_ylim([-0.1, 1.1])

        self.UpdateCESTFrame()

    def OnCEST_ppm_change(self, event):
        min_val = float(self.CEST_ppm_min_text.GetValue())
        max_val = float(self.CEST_ppm_max_text.GetValue())

        # Find the indexes of the min and max values in ppms_1
        min_index = np.argmin(np.abs(self.main_frame.ppms_1 - min_val))
        max_index = np.argmin(np.abs(self.main_frame.ppms_1 - max_val))

        if min_index > max_index:
            min_index, max_index = max_index, min_index

        # Get the data between the min and max values
        self.cest_on_data_new = []
        self.cest_off_data_new = []
        for i in range(len(self.cest_on_data)):
            self.cest_on_data_new.append(
                self.cest_on_data[i].tolist()[min_index:max_index]
            )
            self.cest_off_data_new.append(self.cest_off_data[i][min_index:max_index])

        self.ax.clear()
        for i, y in enumerate(self.offsets_ppm):
            self.ax.plot(
                self.main_frame.ppms_1[min_index:max_index],
                np.full_like(self.main_frame.ppms_1[min_index:max_index], y),
                self.cest_on_data_new[i],
                color="tab:grey",
                linewidth=1.5,
                alpha=0.5,
            )
            if min_index <= self.selected_shift_index <= max_index:
                self.ax.scatter(
                    self.main_frame.ppms_1[self.selected_shift_index],
                    y,
                    self.non_normalized_cest_data[i],
                    color="red",
                    s=10,
                )

        self.ax.set_title("CEST data")
        self.ax.set_xlabel(self.main_frame.nmrdata.axislabels[1])
        self.ax.set_ylabel(self.main_frame.nmrdata.axislabels[0])

        self.ax.set_xlim(
            [
                max(self.main_frame.ppms_1[min_index:max_index]),
                min(self.main_frame.ppms_1[min_index:max_index]),
            ]
        )
        self.ax.set_ylim([max(self.offsets_ppm), min(self.offsets_ppm)])

        self.UpdateCESTFrame()


class DiffusionFit(wx.Frame):
    def __init__(self, title, parent=None):
        self.main_frame = parent
        # Get the monitor size and set the window size to 85% of the monitor size
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        self.display_index = wx.Display.GetFromWindow(parent)
        self.display_index_current = self.display_index
        self.width = int(1.0 * sizes[self.display_index][0])
        self.height = int(0.875 * sizes[self.display_index][1])
        wx.Frame.__init__(
            self, parent=parent, title=title, size=(self.width, self.height)
        )
        self.panel_diffusion = wx.Panel(self, -1)
        self.main_diffusion_sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_diffusion_sizer)

        self.fig_diffusion = Figure()
        self.fig_diffusion.tight_layout()
        self.canvas_diffusion = FigCanvas(self, -1, self.fig_diffusion)
        self.main_diffusion_sizer.Add(self.canvas_diffusion, 10, flag=wx.GROW)
        self.toolbar_diffusion = NavigationToolbar(self.canvas_diffusion)
        self.main_diffusion_sizer.Add(self.toolbar_diffusion, 0, wx.EXPAND)

        self.sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.titlecolor = "black"

        self.initial_values()
        self.make_diffusion_sizer()
        self.plot_diffusion_data()
        self.Show()
        self.Centre()
        # Bind method to check/resize the window when the frame is moved
        self.Bind(wx.EVT_MOVE, self.OnMoveFrame)

        # Bind method to resize the window when the frame is resized
        self.Bind(wx.EVT_SIZE, self.OnSizeFrame)

    def OnMoveFrame(self, event):
        # Get the new default display if the frame is moved
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        display_index = wx.Display.GetFromWindow(self)
        if display_index != self.display_index_current:
            self.display_index_current = display_index
            self.width = int(1.0 * sizes[display_index][0])
            self.height = int(0.875 * sizes[display_index][1])
            self.SetSize((self.width, self.height))
            self.canvas_diffusion.SetSize(
                (
                    self.width * 0.0104,
                    (self.height - self.diffusion_sizer_total.GetMinSize()[1] - 100)
                    * 0.0104,
                )
            )
            self.fig_diffusion.set_size_inches(
                self.width * 0.0104,
                (self.height - self.diffusion_sizer_total.GetMinSize()[1] - 100)
                * 0.0104,
            )
            self.UpdateDiffusionFrame()
        event.Skip()

    def OnSizeFrame(self, event):
        # Get the new frame size
        self.width, self.height = self.GetSize()
        self.SetSize((self.width, self.height))
        self.canvas_diffusion.SetSize(
            (
                self.width * 0.0104,
                (self.height - self.diffusion_sizer_total.GetMinSize()[1] - 100)
                * 0.0104,
            )
        )
        self.fig_diffusion.set_size_inches(
            self.width * 0.0104,
            (self.height - self.diffusion_sizer_total.GetMinSize()[1] - 100) * 0.0104,
        )
        self.UpdateDiffusionFrame()
        event.Skip()

    def UpdateDiffusionFrame(self):
        self.canvas_diffusion.draw()
        self.canvas_diffusion.Refresh()
        self.canvas_diffusion.Update()
        self.panel_diffusion.Refresh()
        self.panel_diffusion.Update()

    # The place where initial global variables are defined
    def initial_values(self):
        self.spectrometer = "Bruker"
        self.nucleus_type = "1H"
        self.bipolar_gradients = False
        self.little_delta = 1000
        self.big_delta = 0.1
        self.integral_factor = 1.0
        self.max_gradient = 53.0
        self.DAC_conversion = 0.002
        self.gamma_dictionary = {}
        self.gamma_dictionary["1H"] = 2.67522e4  # rad s-1 G-1
        self.gamma_dictionary["13C"] = 0.672828e4  # rad s-1 G-1
        self.gamma_dictionary["15N"] = -0.27116e4  # rad s-1 G-1
        self.gamma_dictionary["19F"] = 2.51815e4  # rad s-1 G-1
        self.gamma = self.gamma_dictionary["1H"]
        self.whole_plot = False  # Default to having only the diffusion data in a single plot with no diffusion coefficient subplots
        self.monoexponential_fit = False

        # Input ppms and y data from the main frame
        self.x_data = self.main_frame.new_x_ppms
        self.y_data = self.main_frame.nmrdata.data.T

        # Initially have noise region selection set to false
        self.noise_region_selection = False

        # Create an array to store the min and max ppm values for the selected regions
        self.selected_regions_of_interest = []

        self.AddROI = False

        self.ROI_color = []  # Empty array to store the colors of the ROIs
        self.deleted_ROI_number = (
            0  # Parameter to store the number of ROI's which have been deleted
        )

        self.deleted_slices = (
            []
        )  # Array to hold the indexes of the slices which have been deleted

    def make_diffusion_sizer(self):
        self.diffusion_sizer = wx.BoxSizer(wx.HORIZONTAL)
        # Check box for Varian/Bruker data (default to Bruker)
        self.diffusion_data_type = wx.BoxSizer(wx.HORIZONTAL)
        self.diffusion_data_type_label = wx.StaticBox(self, -1, "Spectrometer Type")
        self.diffusion_data_type_sizer = wx.StaticBoxSizer(
            self.diffusion_data_type_label, wx.HORIZONTAL
        )
        self.diffusion_data_type_sizer.AddSpacer(5)

        # Make a radio button
        self.diffusion_data_type_radio = wx.RadioBox(
            self, -1, choices=["Bruker", "Varian"], style=wx.RA_HORIZONTAL
        )
        self.diffusion_data_type_radio.Bind(wx.EVT_RADIOBOX, self.OnDiffusionDataType)
        if self.spectrometer == "Bruker":
            self.diffusion_data_type_radio.SetSelection(0)
        else:
            self.diffusion_data_type_radio.SetSelection(1)

        self.diffusion_data_type_sizer.Add(self.diffusion_data_type_radio)
        self.diffusion_data_type_sizer.AddSpacer(5)

        # Box to put in all the delay parameters
        self.experimental_parameters_label = wx.StaticBox(self, -1, "Delay Parameters")
        self.experimental_parameters_sizer = wx.StaticBoxSizer(
            self.experimental_parameters_label, wx.HORIZONTAL
        )
        self.experimental_parameters_sizer.AddSpacer(5)

        # Then can have a find gradient percentages button (will use difframp file for Bruker and procpar for Varian)
        # Can have a TextCtrl for the little delta, big delta with a button which can search for these values in the acqus/procpar file and fill them in automatically
        # If that fails can pop up with a window where the gradient percentages used can be entered manually
        # Can then have a TextCtrl for the little delta, big delta with a button which can search for these values in the acqus/procpar file and fill them in automatically
        # If that fails the user can enter them manually
        # Checkbox for whether bipolar gradients were used in the experiment or not (default to no)
        self.bipolar_gradients_checkbox = wx.CheckBox(self, label="Bipolar Gradients")
        self.bipolar_gradients_checkbox.SetValue(False)
        self.bipolar_gradients_checkbox.Bind(wx.EVT_CHECKBOX, self.OnBipolarGradients)
        self.experimental_parameters_sizer.Add(self.bipolar_gradients_checkbox)
        self.little_delta_label = wx.StaticText(self, -1, "δ (μs):")
        self.experimental_parameters_sizer.AddSpacer(5)
        self.experimental_parameters_sizer.Add(self.little_delta_label)
        self.experimental_parameters_sizer.AddSpacer(5)
        self.little_delta_box = wx.TextCtrl(
            self, -1, str(self.little_delta), size=(50, -1)
        )
        self.experimental_parameters_sizer.Add(self.little_delta_box)
        self.big_delta_label = wx.StaticText(self, -1, "Δ (s):")
        self.experimental_parameters_sizer.AddSpacer(5)
        self.experimental_parameters_sizer.Add(self.big_delta_label)
        self.experimental_parameters_sizer.AddSpacer(5)
        self.big_delta_box = wx.TextCtrl(self, -1, str(self.big_delta), size=(50, -1))
        self.experimental_parameters_sizer.Add(self.big_delta_box)
        self.find_parameters_button = wx.Button(self, -1, "Find Parameters")
        self.find_parameters_button.Bind(wx.EVT_BUTTON, self.find_parameters)
        self.experimental_parameters_sizer.AddSpacer(5)
        self.experimental_parameters_sizer.Add(self.find_parameters_button)
        self.experimental_parameters_sizer.AddSpacer(5)

        # Box to put in all the gradient parameters
        self.gradient_parameters_label = wx.StaticBox(self, -1, "Gradient Parameters")
        self.gradient_parameters_sizer = wx.StaticBoxSizer(
            self.gradient_parameters_label, wx.HORIZONTAL
        )
        self.gradient_parameters_sizer.AddSpacer(5)

        if self.spectrometer == "Bruker":
            # Have a box to put in the gradient integral factor (default to 1)
            self.integral_factor_label = wx.StaticText(
                self, -1, "Gradient Integral Factor:"
            )
            self.gradient_parameters_sizer.Add(self.integral_factor_label)
            self.gradient_parameters_sizer.AddSpacer(5)
            self.integral_factor_box = wx.TextCtrl(
                self, -1, str(self.integral_factor), size=(30, -1)
            )
            self.gradient_parameters_sizer.Add(self.integral_factor_box)
            self.gradient_parameters_sizer.AddSpacer(5)

            # Have a box where the user can insert the max spectrometer gradient (default to 53G/cm for Bruker)
            self.max_gradient_label = wx.StaticText(self, -1, "Max Gradient (G/cm):")
            self.gradient_parameters_sizer.Add(self.max_gradient_label)
            self.gradient_parameters_sizer.AddSpacer(5)
            self.max_gradient_box = wx.TextCtrl(
                self, -1, str(self.max_gradient), size=(30, -1)
            )
            self.gradient_parameters_sizer.Add(self.max_gradient_box)
            self.gradient_parameters_sizer.AddSpacer(5)

            self.find_gradient_percentages_button = wx.Button(
                self, -1, "Find Gradient Percentages"
            )
            self.find_gradient_percentages_button.Bind(
                wx.EVT_BUTTON, self.find_gradient_percentages
            )
            self.gradient_parameters_sizer.Add(self.find_gradient_percentages_button)
            self.gradient_parameters_sizer.AddSpacer(5)

        else:
            self.max_gradient = 60.0
            # Have a box to put in the gradient integral factor (default to 1)
            self.integral_factor_label = wx.StaticText(
                self, -1, "Gradient Integral Factor:"
            )
            self.gradient_parameters_sizer.Add(self.integral_factor_label)
            self.gradient_parameters_sizer.AddSpacer(5)
            self.integral_factor_box = wx.TextCtrl(
                self, -1, str(self.integral_factor), size=(30, -1)
            )
            self.gradient_parameters_sizer.Add(self.integral_factor_box)
            self.gradient_parameters_sizer.AddSpacer(5)

            # Have a box where the user can insert the max spectrometer gradient (default to 53G/cm for Bruker)
            self.max_gradient_label = wx.StaticText(self, -1, "Max Gradient (G/cm):")
            self.gradient_parameters_sizer.Add(self.max_gradient_label)
            self.gradient_parameters_sizer.AddSpacer(5)
            self.max_gradient_box = wx.TextCtrl(
                self, -1, str(self.max_gradient), size=(30, -1)
            )
            self.gradient_parameters_sizer.Add(self.max_gradient_box)
            self.gradient_parameters_sizer.AddSpacer(5)

            # Have a box for DAC-G/cm conversion (default = 0.002)
            self.dac_conversion_label = wx.StaticText(
                self, -1, "DAC to G/cm Conversion:"
            )
            self.gradient_parameters_sizer.Add(self.dac_conversion_label)
            self.gradient_parameters_sizer.AddSpacer(5)
            self.dac_conversion_box = wx.TextCtrl(
                self, -1, str(self.DAC_conversion), size=(50, -1)
            )
            self.gradient_parameters_sizer.Add(self.dac_conversion_box)
            self.gradient_parameters_sizer.AddSpacer(5)

            self.find_gradient_percentages_button = wx.Button(
                self, -1, "Find Gradients"
            )
            self.find_gradient_percentages_button.Bind(
                wx.EVT_BUTTON, self.find_gradient_percentages
            )
            self.gradient_parameters_sizer.Add(self.find_gradient_percentages_button)
            self.gradient_parameters_sizer.AddSpacer(5)

        # Create a button to open a textbox window where a user can input the gradient values manually
        self.input_gradients_text = wx.Button(self, -1, "Input Manually")
        self.input_gradients_text.Bind(wx.EVT_BUTTON, self.input_gradients_text_button)
        self.gradient_parameters_sizer.Add(self.input_gradients_text)
        self.gradient_parameters_sizer.AddSpacer(5)

        # Add all sizers to first row of buttons
        self.diffusion_sizer.AddSpacer(5)
        self.diffusion_sizer.Add(self.diffusion_data_type_sizer)
        self.diffusion_sizer.AddSpacer(5)
        self.diffusion_sizer.Add(self.experimental_parameters_sizer)
        self.diffusion_sizer.AddSpacer(5)
        self.diffusion_sizer.Add(self.gradient_parameters_sizer)
        self.diffusion_sizer.AddSpacer(5)

        self.diffusion_sizer_total = wx.BoxSizer(wx.VERTICAL)
        self.diffusion_sizer_total.Add(self.diffusion_sizer)

        # Can have a drop down menu for the nuclei used in the experiment (default to 1H, other options include 19F, 13C, 15N)
        self.nucleus_label = wx.StaticBox(self, -1, "Nucleus:")
        self.nucleus_sizer = wx.StaticBoxSizer(self.nucleus_label, wx.HORIZONTAL)
        self.nucleus_sizer.AddSpacer(5)
        self.nucleus_choices = ["1H", "19F", "13C", "15N"]
        self.nucleus_dropdown = wx.Choice(self, -1, choices=self.nucleus_choices)
        self.nucleus_dropdown.SetSelection(0)
        self.nucleus_sizer.Add(self.nucleus_dropdown)
        self.nucleus_sizer.AddSpacer(5)
        self.nucleus_dropdown.Bind(wx.EVT_CHOICE, self.OnNucleusChoice)

        # Then have button which will allow a user to drag over a section where they wish to estimate the noise level
        # This can then be plotted as a shaded region on the plot
        self.noise_label = wx.StaticBox(self, -1, "Noise Region")
        self.noise_sizer = wx.StaticBoxSizer(self.noise_label, wx.HORIZONTAL)
        self.select_noise_button = wx.Button(self, -1, "Select Noise Region")
        self.select_noise_button.Bind(wx.EVT_BUTTON, self.OnSelectNoise)
        self.noise_sizer.AddSpacer(5)
        self.noise_sizer.Add(self.select_noise_button)
        self.noise_sizer.AddSpacer(5)

        # Then have a TextCtrl for the minimum SNR for the diffusion coefficient to be estimated (default to 10)
        self.noise_factor = 10
        self.noise_factor_label = wx.StaticText(self, -1, "Minimum SNR:")
        self.noise_factor_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.noise_sizer.AddSpacer(10)
        self.noise_sizer.Add(self.noise_factor_label)
        self.noise_sizer.AddSpacer(5)
        self.noise_factor_box = wx.TextCtrl(
            self, -1, str(self.noise_factor), size=(50, -1)
        )
        self.noise_sizer.Add(self.noise_factor_box)
        self.noise_sizer.AddSpacer(5)

        # Need to have a fitting sizer which will contain all the fitting buttons
        self.fitting_label = wx.StaticBox(self, -1, "Fitting")
        self.fitting_sizer = wx.StaticBoxSizer(self.fitting_label, wx.HORIZONTAL)
        self.fitting_sizer.AddSpacer(5)

        # Can then have a button which will fit the Stejskal Tanner equation at all ppms across the whole spectrum that are higher than the noise level
        self.whole_spectrum_fitting_button = wx.Button(self, -1, "Fit Whole Spectrum")
        self.whole_spectrum_fitting_button.Bind(
            wx.EVT_BUTTON, self.OnWholeSpectrumFitting
        )
        self.fitting_sizer.Add(self.whole_spectrum_fitting_button)
        self.fitting_sizer.AddSpacer(5)
        # This can then be plotted
        # In addition, for each ppm can get a plot of I/I0 for all points which is also plotted next to this

        # Then can have a button saying select region of interest. The diffusion coefficient in this region can be estimated along with an error (from the standard deviation of the points)
        # Can plot this distribution of diffusion coefficients and it should resemble a Gaussian distribution

        # Have a button to add a new region of interest
        self.add_region_button = wx.Button(self, -1, "Add ROI")
        self.add_region_button.Bind(wx.EVT_BUTTON, self.OnAddROI)
        self.fitting_sizer.Add(self.add_region_button)
        self.fitting_sizer.AddSpacer(5)

        # Have a button to delete a region of interest
        self.delete_region_button = wx.Button(self, -1, "Delete ROI")
        self.delete_region_button.Bind(wx.EVT_BUTTON, self.OnDeleteROI)
        self.fitting_sizer.Add(self.delete_region_button)
        self.fitting_sizer.AddSpacer(5)

        # Can then have a button which will fit the Stejskal Tanner equation to the mean values of the points above the noise in the region of interest
        self.region_fitting_button = wx.Button(self, -1, "Fit")
        self.region_fitting_button.Bind(wx.EVT_BUTTON, self.OnRegionFitting)
        self.fitting_sizer.Add(self.region_fitting_button)
        self.fitting_sizer.AddSpacer(5)

        # Have a button which will perform a biexponential fit on the data in the region of interest
        self.biexponential_fitting_button = wx.Button(self, -1, "Biexponential Fit")
        self.biexponential_fitting_button.Bind(
            wx.EVT_BUTTON, self.OnBiexponentialFitting
        )
        self.fitting_sizer.Add(self.biexponential_fitting_button)
        self.fitting_sizer.AddSpacer(5)

        # Have a box containing other functions such as a button to delete a slice from the plot and repeat the fitting
        self.other_functions_label = wx.StaticBox(self, -1, "Other Functions")
        self.other_functions_sizer = wx.StaticBoxSizer(
            self.other_functions_label, wx.HORIZONTAL
        )
        self.other_functions_sizer.AddSpacer(5)
        self.delete_slice_button = wx.Button(self, -1, "Delete Slice")
        self.delete_slice_button.Bind(wx.EVT_BUTTON, self.OnDeleteSlice)
        self.other_functions_sizer.Add(self.delete_slice_button)
        self.other_functions_sizer.AddSpacer(5)

        # Creating a sizer for changing the y axis limits in the spectrum
        self.intensity_label = wx.StaticBox(self, -1, "Y Axis Zoom (%):")
        self.intensity_sizer = wx.StaticBoxSizer(self.intensity_label, wx.VERTICAL)
        width = 100
        self.intensity_slider = FloatSlider(
            self, id=-1, value=0, minval=-1, maxval=10, res=0.01, size=(width, height)
        )
        self.intensity_slider.Bind(wx.EVT_SLIDER, self.OnIntensityScrollDiffusion)
        self.intensity_sizer.AddSpacer(5)
        self.intensity_sizer.Add(self.intensity_slider)

        self.diffusion_fitting_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.diffusion_fitting_sizer.AddSpacer(5)
        self.diffusion_fitting_sizer.Add(self.nucleus_sizer)
        self.diffusion_fitting_sizer.AddSpacer(5)
        self.diffusion_fitting_sizer.Add(self.noise_sizer)
        self.diffusion_fitting_sizer.AddSpacer(5)
        self.diffusion_fitting_sizer.Add(self.fitting_sizer)
        self.diffusion_fitting_sizer.AddSpacer(5)
        self.diffusion_fitting_sizer.Add(self.other_functions_sizer)
        self.diffusion_fitting_sizer.AddSpacer(5)
        self.diffusion_fitting_sizer.Add(self.intensity_sizer)

        self.diffusion_sizer_total.AddSpacer(5)
        self.diffusion_sizer_total.Add(self.diffusion_fitting_sizer)
        self.diffusion_sizer_total.AddSpacer(5)

        self.main_diffusion_sizer.Add(self.diffusion_sizer_total)

    def plot_diffusion_data(self):

        self.ax_diffusion = self.fig_diffusion.add_subplot(111)
        count = 1
        self.slice_plots = []
        for i, data in enumerate(self.y_data):
            (line,) = self.ax_diffusion.plot(
                self.x_data, data, linewidth=0.5, label=str(count)
            )
            self.slice_plots.append(line)
            count += 1
        self.ax_diffusion.set_xlim([self.x_data[0], self.x_data[-1]])

        self.ax_diffusion.set_xlabel(self.main_frame.nmrdata.axislabels[1])
        legend = self.ax_diffusion.legend(title="Slice Number")
        legend.get_title().set_color(self.titlecolor)
        self.ax_diffusion.set_ylabel("Intensity")

        self.noise_region = self.ax_diffusion.axvspan(
            min(self.x_data), min(self.x_data), alpha=0.2, color="gray"
        )

    def OnDiffusionDataType(self, event):
        # Find out whether Bruker or Varian data is being used
        if self.diffusion_data_type_radio.GetSelection() == 0:
            self.spectrometer = "Bruker"
        else:
            self.spectrometer = "Varian"
        self.diffusion_sizer.Clear(True)
        self.diffusion_fitting_sizer.Clear(True)
        self.make_diffusion_sizer()
        self.Refresh()
        self.Layout()

    def OnIntensityScrollDiffusion(self, event):
        # Function to change the y axis limits
        intensity_percent = 10 ** float(self.intensity_slider.GetValue())

        self.ax_diffusion.set_ylim(
            -(np.max(self.y_data) / 8) / (intensity_percent / 100),
            np.max(self.y_data) / (intensity_percent / 100),
        )
        self.UpdateDiffusionFrame()

    def OnBipolarGradients(self, event):
        if self.bipolar_gradients_checkbox.GetValue() == True:
            self.bipolar_gradients = True
        else:
            self.bipolar_gradients = False

    def find_parameters(self, event):
        if self.spectrometer == "Bruker":
            # Search through acqus file to get the little delta (p30) and big delta (d20) values used
            try:
                file = open("acqus", "r")
            except:
                # Give an error message saying unable to find acqus file
                msg = wx.MessageDialog(
                    self,
                    "Unable to find acqus file. Please input delays manually",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                return

            try:
                delays_total = []
                durations_total = []
                add_delays = False
                add_durations = False
                for line in file:
                    if add_delays == True:
                        if "##" in line:
                            add_delays = False
                            continue
                        else:
                            delays = line.split("\n")[0].split()
                            for delay in delays:
                                delays_total.append(float(delay))
                    if "##$D=" in line:
                        add_delays = True

                    if add_durations == True:
                        if "##" in line:
                            add_durations = False
                            continue
                        else:
                            durations = line.split("\n")[0].split()
                            for duration in durations:
                                durations_total.append(float(duration))
                    if "##$P=" in line:
                        add_durations = True
                self.big_delta = delays_total[20]
                if self.bipolar_gradients == True:
                    self.small_delta = durations_total[30] * 2
                else:
                    self.small_delta = durations_total[30]
            except:
                # Give an error message saying unable to find delays in the acqus file (./acqus)
                msg = wx.MessageDialog(
                    self,
                    "Unable to find delays in the acqus file (./acqus). Please input delays manually",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                return

            file.close()

        else:
            # Search through Varian procpar file to find out the little delta and big delta values used
            try:
                self.dic, self.data = ng.varian.read("./")
            except:
                # Give an error message saying unable to find procpar file
                msg = wx.MessageDialog(
                    self,
                    "Unable to find procpar file. Please input delays manually",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                return

            try:
                # Find big delta and small delta
                self.big_delta = float(self.dic["procpar"]["BigT"]["values"][0])
                self.small_delta = float(self.dic["procpar"]["gt1"]["values"][0]) * 1e6
                if self.bipolar_gradients == True:
                    self.small_delta = self.small_delta * 2
            except:
                # Give an error message saying unable to find delays in the procpar file (./procpar)
                msg = wx.MessageDialog(
                    self,
                    "Unable to find delays in the procpar file (./procpar). Please input delays manually",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                return

        # Set the little delta and big delta values in the GUI to the found values
        self.little_delta_box.SetValue(str(self.small_delta))
        self.big_delta_box.SetValue(str(self.big_delta))

    def find_gradient_percentages(self, event):
        if self.spectrometer == "Bruker":
            self.max_gradient = float(self.max_gradient_box.GetValue())
            # Search through the difframp file to get the gradient percentages used
            try:
                with open("./lists/gp/Difframp", "r") as file:
                    self.gradients_percent = []
                    skip_line = True
                    for line in file:
                        if skip_line == True:
                            if "##XYDATA= (X++(Y..Y))" not in line:
                                pass
                            else:
                                skip_line = False
                        else:
                            if "##END=" not in line:
                                self.gradients_percent.append(
                                    float(line.split()[0]) * 100
                                )

                self.gradients = (
                    np.array(self.gradients_percent) / 100
                ) * self.max_gradient

                if len(self.y_data) != len(self.gradients):
                    for i, deleted_slice in enumerate(self.deleted_slices):
                        self.gradients = np.delete(self.gradients, deleted_slice, 0)
                        self.gradients_percent = np.delete(
                            self.gradients_percent, deleted_slice, 0
                        )

                if len(self.y_data) != len(self.gradients):
                    # Give an error message saying unable to find gradient percentages in the difframp file (./lists/gp/Difframp)
                    msg = wx.MessageDialog(
                        self,
                        "Number of gradients in difframp is not equal to the data size",
                        "Error",
                        wx.OK | wx.ICON_ERROR,
                    )
                    msg.ShowModal()
                    msg.Destroy()
                    # Bring up a window where the user can enter the gradient percentages manually (TextCtrl for min/max gradient percentages and then a radiobox for linear, squared, exponential distribution)
                    # Can then press okay and will produce gradient percentages and gradients manually
                    self.gradients_percent = []
                    self.gradients = []
                    self.gradients_input_manual = DiffusionGradientManualInput(
                        title="Manual Gradient Input",
                        parent=self,
                        spectrometer=self.spectrometer,
                    )

                else:
                    # Give a pop out window showing the gradient percentages and values used
                    gradient_percent_string = "Gradient Percentages (%): "
                    gradient_string = "Gradients (G/cm): "
                    for i, gradient_percent in enumerate(self.gradients_percent):
                        gradient_percent_string = (
                            gradient_percent_string
                            + "{:.2f}, ".format(gradient_percent)
                        )
                        gradient_string = gradient_string + "{:.2f}, ".format(
                            self.gradients[i]
                        )

                    gradient_percent_string = gradient_percent_string[:-2]
                    gradient_string = gradient_string[:-2]

                    msg = wx.MessageDialog(
                        self,
                        gradient_percent_string + "\n" + gradient_string,
                        "Gradient Percentages and Values",
                        wx.OK | wx.ICON_INFORMATION,
                    )
                    msg.ShowModal()
                    msg.Destroy()

            except:
                try:
                    # Try to read in gradients from (./Difframp) file, older versions of topspin save the file here
                    with open("./Difframp", "r") as file:
                        self.gradients_percent = []
                        skip_line = True
                        for line in file:
                            if skip_line == True:
                                if "##XYDATA= (X++(Y..Y))" not in line:
                                    pass
                                else:
                                    skip_line = False
                            else:
                                if "##END=" not in line:
                                    self.gradients_percent.append(
                                        float(line.split()[0]) * 100
                                    )

                    self.gradients = (
                        np.array(self.gradients_percent) / 100
                    ) * self.max_gradient

                    if len(self.y_data) != len(self.gradients):
                        for i, deleted_slice in enumerate(self.deleted_slices):
                            self.gradients = np.delete(self.gradients, deleted_slice, 0)
                            self.gradients_percent = np.delete(
                                self.gradients_percent, deleted_slice, 0
                            )

                    if len(self.y_data) != len(self.gradients):
                        # Give an error message saying unable to find gradient percentages in the difframp file (./lists/gp/Difframp)
                        msg = wx.MessageDialog(
                            self,
                            "Number of gradients in difframp is not equal to the data size",
                            "Error",
                            wx.OK | wx.ICON_ERROR,
                        )
                        msg.ShowModal()
                        msg.Destroy()
                        # Bring up a window where the user can enter the gradient percentages manually (TextCtrl for min/max gradient percentages and then a radiobox for linear, squared, exponential distribution)
                        # Can then press okay and will produce gradient percentages and gradients manually
                        self.gradients_percent = []
                        self.gradients = []
                        self.gradients_input_manual = DiffusionGradientManualInput(
                            title="Manual Gradient Input",
                            parent=self,
                            spectrometer=self.spectrometer,
                        )

                    else:
                        # Give a pop out window showing the gradient percentages and values used
                        gradient_percent_string = "Gradient Percentages (%): "
                        gradient_string = "Gradients (G/cm): "
                        for i, gradient_percent in enumerate(self.gradients_percent):
                            gradient_percent_string = (
                                gradient_percent_string
                                + "{:.2f}, ".format(gradient_percent)
                            )
                            gradient_string = gradient_string + "{:.2f}, ".format(
                                self.gradients[i]
                            )

                        gradient_percent_string = gradient_percent_string[:-2]
                        gradient_string = gradient_string[:-2]

                        msg = wx.MessageDialog(
                            self,
                            gradient_percent_string + "\n" + gradient_string,
                            "Gradient Percentages and Values",
                            wx.OK | wx.ICON_INFORMATION,
                        )
                        msg.ShowModal()
                        msg.Destroy()

                except:

                    # Give an error message saying unable to find gradient percentages in the difframp file (./lists/gp/Difframp)
                    msg = wx.MessageDialog(
                        self,
                        "Unable to find gradient percentages in the difframp file (./lists/gp/Difframp or ./Difframp). Please input gradients manually",
                        "Error",
                        wx.OK | wx.ICON_ERROR,
                    )
                    msg.ShowModal()
                    msg.Destroy()
                    # Bring up a window where the user can enter the gradient percentages manually (TextCtrl for min/max gradient percentages and then a radiobox for linear, squared, exponential distribution)
                    # Can then press okay and will produce gradient percentages and gradients manually
                    self.gradients_percent = []
                    self.gradients = []
                    self.gradients_input_manual = DiffusionGradientManualInput(
                        title="Manual Gradient Input",
                        parent=self,
                        spectrometer=self.spectrometer,
                    )

        else:
            self.DAC_conversion = float(self.dac_conversion_box.GetValue())
            self.max_gradient = float(self.max_gradient_box.GetValue())
            self.gradient_list = []

            # Search through the procpar file to get the gradient percentages used
            try:
                self.dic, self.data = ng.varian.read("./")

                # Get the gradient strength parameters
                gradient_name = self.dic["procpar"]["array"]["values"][0]

                # separate the data for each gradient strength
                for i, gradient in enumerate(
                    self.dic["procpar"][gradient_name]["values"]
                ):
                    self.gradient_list.append(float(gradient))
                self.gradients = np.array(self.gradient_list) * self.DAC_conversion
                self.gradients_percent = (
                    np.array(self.gradient_list)
                    * self.DAC_conversion
                    / self.max_gradient
                    * 100
                )

                # Give a pop out window showing the gradient percentages and values used
                gradient_percent_string = "Gradient Percentages (%): "
                gradient_string = "Gradients (G/cm): "
                for i, gradient_percent in enumerate(self.gradients_percent):
                    gradient_percent_string = (
                        gradient_percent_string + "{:.2f}, ".format(gradient_percent)
                    )
                    gradient_string = gradient_string + "{:.2f}, ".format(
                        self.gradients[i]
                    )

                gradient_percent_string = gradient_percent_string[:-2]
                gradient_string = gradient_string[:-2]

                msg = wx.MessageDialog(
                    self,
                    gradient_percent_string + "\n" + gradient_string,
                    "Gradient Percentages and Values",
                    wx.OK | wx.ICON_INFORMATION,
                )
                msg.ShowModal()
                msg.Destroy()

                if len(self.y_data) != len(self.gradients):
                    for i, deleted_slice in enumerate(self.deleted_slices):
                        self.gradients = np.delete(self.gradients, deleted_slice, 0)
                        self.gradients_percent = np.delete(
                            self.gradients_percent, deleted_slice, 0
                        )

            except:
                # Give an error message saying unable to find procpar file
                msg = wx.MessageDialog(
                    self,
                    "Unable to find procpar file. Please input gradients manually",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                # Bring up a window where the user can enter the gradient percentages manually (TextCtrl for min/max gradient percentages and then a radiobox for linear, squared, exponential distribution)
                # Can then press okay and will produce gradient percentages and gradients manually
                self.gradients_percent = []
                self.gradients = []
                self.gradients_input_manual = DiffusionGradientManualInput(
                    title="Manual Gradient Input",
                    parent=self,
                    spectrometer=self.spectrometer,
                )

    def OnNucleusChoice(self, event):
        self.nucleus_type = self.nucleus_choices[self.nucleus_dropdown.GetSelection()]
        self.gamma = self.gamma_dictionary[self.nucleus_type]

    def OnSelectNoise(self, event):
        # self.noise_region.set_xy([[min(self.x_data),0],[min(self.x_data),1],[min(self.x_data),1],[min(self.x_data),0]])
        # if(self.whole_plot==True):
        #     self.noise_region_2.set_xy([[min(self.x_data),0],[min(self.x_data),1],[min(self.x_data),1],[min(self.x_data),0]])
        #     self.noise_region_3.set_xy([[min(self.x_data),0],[min(self.x_data),1],[min(self.x_data),1],[min(self.x_data),0]])

        self.UpdateDiffusionFrame()

        self.press = False
        self.move = False
        if self.whole_plot == True:
            self.canvas_diffusion = FigCanvas(self, -1, self.fig_diffusion)
        self.noise_select_press = self.canvas_diffusion.mpl_connect(
            "button_press_event", self.OnPress
        )
        self.noise_select_release = self.canvas_diffusion.mpl_connect(
            "button_release_event", self.OnReleaseNoise
        )
        self.noise_select_move = self.canvas_diffusion.mpl_connect(
            "motion_notify_event", self.OnMove
        )

    def OnPress(self, event):
        if self.whole_plot == False:
            if event.inaxes == self.ax_diffusion:
                self.press = True
                self.x0 = event.xdata
        else:
            if (
                event.inaxes == self.ax_diffusion
                or event.inaxes == self.ax_diffusion_whole_fit
                or event.inaxes == self.ax_diffusion_I0_whole_fit
                or event.inaxes == self.diffusion_coefficient_plot
            ):
                self.press = True
                self.x0 = event.xdata

    def OnMove(self, event):
        if self.whole_plot == False:
            if event.inaxes == self.ax_diffusion:
                self.move_noise(event)

        else:
            if (
                event.inaxes == self.ax_diffusion
                or event.inaxes == self.ax_diffusion_whole_fit
                or event.inaxes == self.ax_diffusion_I0_whole_fit
                or event.inaxes == self.diffusion_coefficient_plot
            ):
                self.move_noise(event)

    def move_noise(self, event):
        if self.press:
            self.move = True
            self.x1 = event.xdata
            if self.x1 > self.x0:
                xmax = self.x1
                xmin = self.x0
            else:
                xmax = self.x0
                xmin = self.x1
            self.noise_region.set_x(xmin)
            self.noise_region.set_width(xmax - xmin)
            # self.noise_region.set(xy=[[xmin,0],[xmin,1],[xmax,1],[xmax,0]])  # no longer works in recent matplotlib versions
            if self.whole_plot == True:
                # self.noise_region_2.set(xy=[[xmin,0],[xmin,1],[xmax,1],[xmax,0]])   # no longer works in recent matplotlib versions
                # self.noise_region_3.set(xy=[[xmin,0],[xmin,1],[xmax,1],[xmax,0]])   # no longer works in recent matplotlib versions

                self.noise_region_2.set_x(xmin)
                self.noise_region_2.set_width(xmax - xmin)

                self.noise_region_3.set_x(xmin)
                self.noise_region_3.set_width(xmax - xmin)

            self.UpdateDiffusionFrame()

    def release_noise(self, event):
        if self.press:
            self.x2 = event.xdata
            if self.x2 > self.x0:
                xmax = self.x2
                xmin = self.x0
            else:
                xmax = self.x0
                xmin = self.x2
            self.noise_x_initial = xmin
            self.noise_x_final = xmax

            self.noise_region.set_x(xmin)
            self.noise_region.set_width(xmax - xmin)
            # self.noise_region.set(xy=[[xmin,0],[xmin,1],[xmax,1],[xmax,0]])  # no longer works in recent matplotlib versions
            if self.whole_plot == True:
                # self.noise_region_2.set(xy=[[xmin,0],[xmin,1],[xmax,1],[xmax,0]])  # no longer works in recent matplotlib versions
                # self.noise_region_3.set(xy=[[xmin,0],[xmin,1],[xmax,1],[xmax,0]])  # no longer works in recent matplotlib versions

                self.noise_region_2.set_x(xmin)
                self.noise_region_2.set_width(xmax - xmin)

                self.noise_region_3.set_x(xmin)
                self.noise_region_3.set_width(xmax - xmin)

            self.UpdateDiffusionFrame()
        self.press = False
        self.move = False
        self.canvas_diffusion.mpl_disconnect(self.noise_select_press)
        self.canvas_diffusion.mpl_disconnect(self.noise_select_move)
        self.canvas_diffusion.mpl_disconnect(self.noise_select_release)

        # Turn on noise region selection flag
        self.noise_region_selection = True

    def OnReleaseNoise(self, event):
        if self.whole_plot == False:
            if event.inaxes == self.ax_diffusion:
                self.release_noise(event)

        else:
            if (
                event.inaxes == self.ax_diffusion
                or event.inaxes == self.ax_diffusion_whole_fit
                or event.inaxes == self.ax_diffusion_I0_whole_fit
                or event.inaxes == self.diffusion_coefficient_plot
            ):
                self.release_noise(event)

    def OnWholeSpectrumFitting(self, event):
        # Initially check that the noise region has been selected
        try:
            self.noise_x_initial
            self.noise_x_final
        except:
            # Give an error message saying noise region has not been selected
            msg = wx.MessageDialog(
                self,
                "Please select a noise region before fitting",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return

        # Check that the minimum SNR has been entered and is a value greater than 0
        try:
            self.noise_factor = float(self.noise_factor_box.GetValue())
        except:
            # Give an error message saying noise factor has not been entered
            msg = wx.MessageDialog(
                self,
                "Please enter a minimum SNR before fitting",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return

        if self.noise_factor <= 0:
            # Give an error message saying noise factor has not been entered correctly
            msg = wx.MessageDialog(
                self,
                "Please enter a minimum SNR greater than 0 before fitting",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return

        # Check that the gradients have been found
        try:
            self.gradients
        except:
            # Give an error message saying gradients have not been found
            msg = wx.MessageDialog(
                self,
                "Please input the gradients before fitting",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return

        # Find out the standard deviation of the noise region
        noise_region = self.y_data[
            :,
            np.where(
                (self.x_data >= self.noise_x_initial)
                & (self.x_data <= self.noise_x_final)
            )[0],
        ]
        noise_region_std = np.std(noise_region)

        # Find out the ppms which have intensity of the minimum intensity slice greater than the noise level in all slices
        self.ppms_above_noise_slices = []
        for i, data in enumerate(self.y_data):
            self.ppms_above_noise_slice = []
            for j, intensity in enumerate(data):
                if intensity > self.noise_factor * noise_region_std:
                    self.ppms_above_noise_slice.append(self.x_data[j])

            self.ppms_above_noise_slices.append(self.ppms_above_noise_slice)

        # Find out the ppms which are in all slices and remove ones which are not
        self.ppms_above_noise = self.ppms_above_noise_slices[0]
        for i, ppm in enumerate(self.ppms_above_noise):
            for j, ppms_slice in enumerate(self.ppms_above_noise_slices):
                if ppm not in ppms_slice:
                    del self.ppms_above_noise[i]
                    break

        # Get the indices of all the ppms which are above the noise level
        self.ppms_above_noise_indices = []
        for i, ppm in enumerate(self.ppms_above_noise):
            self.ppms_above_noise_indices.append(np.where(self.x_data == ppm)[0][0])

        # Remove all the y data points which are below the noise level
        self.SelectDataAboveThreshold()

        # Separate data into point by point
        self.SeparateDataIntoPointByPoint()

        # For all the ppms that have intensity above the noise threshold, fit the Stejskal Tanner equation to the data
        self.fitted_I0_global = []
        self.fitted_D_global = []

        self.little_delta = float(self.little_delta_box.GetValue())
        self.big_delta = float(self.big_delta_box.GetValue())

        for i, ppm in enumerate(self.ppms_above_noise):
            self.y_vals = np.real(self.y_data_point_by_point[i])

            # Start at a few different initial diffusion coefficients so that don't get stuck in local minima
            fits = []
            chi_squareds = []
            for j, D_initial in enumerate(10 ** np.linspace(-5, -10, 10)):
                fit = self.leastsq_global([np.max(self.y_vals), D_initial])
                fits.append(fit)
                chi_squareds.append(np.sum(self.chi_global(fit) ** 2))

            fit = fits[np.argmin(chi_squareds)]
            self.fitted_I0_global.append(fit[0])
            self.fitted_D_global.append(fit[1])

        self.PlotWholeSpectrumFitting()

    def PlotWholeSpectrumFitting(self):
        self.fig_diffusion.clear()
        self.fig_diffusion.tight_layout()

        gs = gridspec.GridSpec(2, 2)

        self.ax_diffusion = self.fig_diffusion.add_subplot(gs[0, :])
        self.ax_diffusion_whole_fit = self.fig_diffusion.add_subplot(
            gs[1, 0], sharex=self.ax_diffusion, sharey=self.ax_diffusion
        )
        self.ax_diffusion_I0_whole_fit = self.fig_diffusion.add_subplot(
            gs[1, 1], sharex=self.ax_diffusion
        )

        count = 1
        self.slice_plots = []
        for i, data in enumerate(self.y_data):
            (line,) = self.ax_diffusion.plot(
                self.x_data, data, linewidth=0.5, label=str(count)
            )
            self.slice_plots.append(line)
            count += 1
        self.ax_diffusion.set_xlim([self.x_data[0], self.x_data[-1]])

        self.ax_diffusion.set_xlabel(self.main_frame.nmrdata.axislabels[1])
        legend = self.ax_diffusion.legend(title="Slice Number")
        legend.get_title().set_color(self.titlecolor)
        self.ax_diffusion.set_ylabel("Intensity")
        self.ax_diffusion.set_title("Diffusion Data", color=self.titlecolor)

        self.noise_region = self.ax_diffusion.axvspan(
            self.noise_x_initial, self.noise_x_final, alpha=0.2, color="gray"
        )

        # Plot the fitted diffusion coefficients and use a twiny to also plot the initial slice of the spectrum
        self.ax_diffusion_whole_fit.plot(self.x_data, self.y_data[0])
        self.ax_diffusion_whole_fit.set_xlabel(self.main_frame.nmrdata.axislabels[1])
        self.ax_diffusion_whole_fit.set_yticks([])
        self.diffusion_coefficient_plot = self.ax_diffusion_whole_fit.twinx()
        self.diffusion_coefficient_plot.set_ylabel(r"Diffusion Coefficient (cm$^2$/s)")
        self.diffusion_coefficient_plot.set_xlabel(
            self.main_frame.nmrdata.axislabels[1]
        )
        self.diffusion_coefficient_plot.set_xlim([self.x_data[0], self.x_data[-1]])
        self.diffusion_coefficient_plot.scatter(
            self.ppms_above_noise, self.fitted_D_global, color="tab:red", s=0.5
        )
        self.diffusion_coefficient_plot.yaxis.tick_left()
        self.diffusion_coefficient_plot.yaxis.set_label_position("left")
        self.diffusion_coefficient_plot.set_title(
            "Diffusion Coefficient vs PPM", color=self.titlecolor
        )
        self.diffusion_coefficient_plot.set_yscale("log")
        self.noise_region_2 = self.ax_diffusion_whole_fit.axvspan(
            self.noise_x_initial, self.noise_x_final, alpha=0.2, color="gray"
        )

        # Plot I/I0 for every chosen ppm across the spectrum for all slices
        for i, selected_y_data in enumerate(self.y_data_above_noise):
            self.ax_diffusion_I0_whole_fit.scatter(
                self.ppms_above_noise,
                np.array(selected_y_data) / self.fitted_I0_global,
                s=0.5,
            )

        self.ax_diffusion_I0_whole_fit.set_xlabel(self.main_frame.nmrdata.axislabels[1])
        self.ax_diffusion_I0_whole_fit.set_ylabel(r"I/I$_0$")
        self.ax_diffusion_I0_whole_fit.set_xlim([self.x_data[0], self.x_data[-1]])
        self.ax_diffusion_I0_whole_fit.set_title(
            r"I/I$_0$ vs PPM", color=self.titlecolor
        )
        self.noise_region_3 = self.ax_diffusion_I0_whole_fit.axvspan(
            self.noise_x_initial, self.noise_x_final, alpha=0.2, color="gray"
        )

        if self.whole_plot == True:
            if len(self.ROI_regions) == 0:
                self.ROI_regions = []
                self.ROI_regions_2 = []
                self.ROI_regions_3 = []
            else:
                ROI_regions = []
                ROI_regions_2 = []
                ROI_regions_3 = []
                for i, ROI_region in enumerate(self.ROI_regions):
                    bottom_left = ROI_region.get_xy()
                    width = ROI_region.get_width()
                    ROI_regions.append(
                        self.ax_diffusion.axvspan(
                            bottom_left[0],
                            bottom_left[0] + width,
                            alpha=0.2,
                            color=self.ROI_color[i],
                        )
                    )
                    ROI_regions_2.append(
                        self.ax_diffusion_whole_fit.axvspan(
                            bottom_left[0],
                            bottom_left[0] + width,
                            alpha=0.2,
                            color=self.ROI_color[i],
                        )
                    )
                    ROI_regions_3.append(
                        self.ax_diffusion_I0_whole_fit.axvspan(
                            bottom_left[0],
                            bottom_left[0] + width,
                            alpha=0.2,
                            color=self.ROI_color[i],
                        )
                    )
                self.ROI_regions = ROI_regions
                self.ROI_regions_2 = ROI_regions_2
                self.ROI_regions_3 = ROI_regions_3

        else:
            self.ROI_regions = []
            self.ROI_regions_2 = []
            self.ROI_regions_3 = []

        # Turn on whole plot mode once whole plot mode has been completed
        self.whole_plot = True

        self.UpdateDiffusionFrame()

    def SelectDataAboveThreshold(self):
        # Remove all the y data points which are below the noise level
        self.y_data_above_noise = []
        for i, data in enumerate(self.y_data):
            self.y_data_above_noise_slice = []
            for index in self.ppms_above_noise_indices:
                self.y_data_above_noise_slice.append(data[index])

            self.y_data_above_noise.append(self.y_data_above_noise_slice)

    def SeparateDataIntoPointByPoint(self):
        # Separate the data into arrays of intensities for each ppm that has intensity above noise threshold in all slices
        self.y_data_point_by_point = []
        for i, ppm in enumerate(self.ppms_above_noise):
            y_data = []
            for j, data in enumerate(self.y_data_above_noise):
                y_data.append(data[i])

            y_data = np.array(y_data)
            self.y_data_point_by_point.append(y_data)

        self.y_data_point_by_point = np.array(self.y_data_point_by_point)

    def StejsktalTanner(self, p0):
        I0, D = p0
        return I0 * np.exp(
            -((self.gamma**2) * (self.gradients**2) * (self.little_delta * 1e-6) ** 2)
            * (self.big_delta - (self.little_delta * 1e-6) / 3)
            * D
        )

    def chi_global(self, p0):
        return self.y_vals - self.StejsktalTanner(p0)

    def leastsq_global(self, p0):
        fit = leastsq(self.chi_global, p0)
        return fit[0]

    def OnAddROI(self, event):
        # Check that the full spectrum has been fitted first
        if self.whole_plot != True:
            # Give an error message saying full spectrum has not been fitted
            msg = wx.MessageDialog(
                self,
                "Please fit the whole spectrum before selecting a region of interest",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return

        # Check if have pressed AddROI but have not selected a region of interest then delete the previous value
        if self.AddROI == True:
            if self.ROI_regions[-1].get_xy()[0][0] == self.x_data[0]:
                del self.ROI_regions[-1]
                del self.ROI_regions_2[-1]
                del self.ROI_regions_3[-1]

        self.AddROI == True

        # Add new region plots with the default values (min ppm values)
        self.ROI_color.append(
            self.main_frame.colours[
                len(self.selected_regions_of_interest) + self.deleted_ROI_number
            ]
        )
        self.ROI_regions.append(
            self.ax_diffusion.axvspan(
                self.x_data[0], self.x_data[0], alpha=0.2, color=self.ROI_color[-1]
            )
        )
        self.ROI_regions_2.append(
            self.ax_diffusion_whole_fit.axvspan(
                self.x_data[0], self.x_data[0], alpha=0.2, color=self.ROI_color[-1]
            )
        )
        self.ROI_regions_3.append(
            self.ax_diffusion_I0_whole_fit.axvspan(
                self.x_data[0], self.x_data[0], alpha=0.2, color=self.ROI_color[-1]
            )
        )

        self.UpdateDiffusionFrame()

        self.canvas_diffusion.mpl_disconnect(self.noise_select_press)
        self.canvas_diffusion.mpl_disconnect(self.noise_select_move)
        self.canvas_diffusion.mpl_disconnect(self.noise_select_release)

        self.press = False
        self.move = False

        self.select_ROI_press = self.canvas_diffusion.mpl_connect(
            "button_press_event", self.OnPressROI
        )
        self.select_ROI_release = self.canvas_diffusion.mpl_connect(
            "button_release_event", self.OnReleaseROI
        )
        self.select_ROI_move = self.canvas_diffusion.mpl_connect(
            "motion_notify_event", self.OnMoveROI
        )

    def OnDeleteROI(self, event):
        # Check that a region of interest has been added first
        if len(self.selected_regions_of_interest) == 0:
            # Give an error message saying no regions of interest have been added
            msg = wx.MessageDialog(
                self,
                "No regions of interest have been added",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return

        # When mouse is over a region of interest highlight that region (make alpha=0.75)
        # When mouse is moved away from region of interest make alpha=0.2 again
        # When mouse is clicked delete that region of interest

        self.canvas_diffusion.mpl_disconnect(self.noise_select_press)
        self.canvas_diffusion.mpl_disconnect(self.noise_select_move)
        self.canvas_diffusion.mpl_disconnect(self.noise_select_release)

        self.canvas_diffusion.mpl_disconnect(self.select_ROI_press)
        self.canvas_diffusion.mpl_disconnect(self.select_ROI_move)
        self.canvas_diffusion.mpl_disconnect(self.select_ROI_release)

        self.delete_ROI_press = self.canvas_diffusion.mpl_connect(
            "button_press_event", self.OnPressDeleteROI
        )
        self.delete_ROI_highlight = self.canvas_diffusion.mpl_connect(
            "motion_notify_event", self.OnHighlightROI
        )

    def OnPressDeleteROI(self, event):
        if (
            event.inaxes == self.ax_diffusion
            or event.inaxes == self.ax_diffusion_whole_fit
            or event.inaxes == self.ax_diffusion_I0_whole_fit
            or event.inaxes == self.diffusion_coefficient_plot
        ):
            # Find out the highlighted slices
            if len(self.selected_regions_of_interest) == 1:
                for i in self.highlighted_regions:
                    self.selected_regions_of_interest = []
                    self.ROI_regions[i].set_alpha(0.2)
                    self.ROI_regions_2[i].set_alpha(0.2)
                    self.ROI_regions_3[i].set_alpha(0.2)

                    self.ROI_regions[i].set_x(self.x_data[0])
                    self.ROI_regions[i].set_width(0)
                    self.ROI_regions_2[i].set_x(self.x_data[0])
                    self.ROI_regions_2[i].set_width(0)
                    self.ROI_regions_3[i].set_x(self.x_data[0])
                    self.ROI_regions_3[i].set_width(0)
                    # self.noise_region.set_width(xmax-xmin)
                    # self.ROI_regions[i].set_xy([[self.x_data[0],0],[self.x_data[0],1],[self.x_data[0],1],[self.x_data[0],0]])
                    # self.ROI_regions_2[i].set_xy([[self.x_data[0],0],[self.x_data[0],1],[self.x_data[0],1],[self.x_data[0],0]])
                    # self.ROI_regions_3[i].set_xy([[self.x_data[0],0],[self.x_data[0],1],[self.x_data[0],1],[self.x_data[0],0]])
                    del self.ROI_regions[i]
                    del self.ROI_regions_2[i]
                    del self.ROI_regions_3[i]
                    del self.ROI_color[i]
                    self.deleted_ROI_number += 1
            else:
                for i in self.highlighted_regions:
                    del self.selected_regions_of_interest[i]
                    self.ROI_regions[i].set_alpha(0.2)
                    self.ROI_regions_2[i].set_alpha(0.2)
                    self.ROI_regions_3[i].set_alpha(0.2)
                    self.ROI_regions[i].set_x(self.x_data[0])
                    self.ROI_regions[i].set_width(0)
                    self.ROI_regions_2[i].set_x(self.x_data[0])
                    self.ROI_regions_2[i].set_width(0)
                    self.ROI_regions_3[i].set_x(self.x_data[0])
                    self.ROI_regions_3[i].set_width(0)
                    # self.ROI_regions[i].set_xy([[self.x_data[0],0],[self.x_data[0],1],[self.x_data[0],1],[self.x_data[0],0]])
                    # self.ROI_regions_2[i].set_xy([[self.x_data[0],0],[self.x_data[0],1],[self.x_data[0],1],[self.x_data[0],0]])
                    # self.ROI_regions_3[i].set_xy([[self.x_data[0],0],[self.x_data[0],1],[self.x_data[0],1],[self.x_data[0],0]])
                    del self.ROI_regions[i]
                    del self.ROI_regions_2[i]
                    del self.ROI_regions_3[i]
                    del self.ROI_color[i]
                    self.deleted_ROI_number += 1

            self.UpdateDiffusionFrame()

            # Disconnect highlight and press events
            self.canvas_diffusion.mpl_disconnect(self.delete_ROI_press)
            self.canvas_diffusion.mpl_disconnect(self.delete_ROI_highlight)

            if self.monoexponential_fit == True:
                self.OnRegionFitting(event)

    def OnHighlightROI(self, event):
        if (
            event.inaxes == self.ax_diffusion
            or event.inaxes == self.ax_diffusion_whole_fit
            or event.inaxes == self.ax_diffusion_I0_whole_fit
            or event.inaxes == self.diffusion_coefficient_plot
        ):
            x0 = event.xdata
            self.highlight_ROI(x0)

    def highlight_ROI(self, x0):
        self.highlighted_regions = []
        # Check if x0 is within any of the regions of interest
        if len(self.selected_regions_of_interest) == 1:
            region = self.selected_regions_of_interest[0]
            if x0 >= region[0] and x0 <= region[1]:
                self.ROI_regions[0].set_alpha(0.75)
                self.ROI_regions_2[0].set_alpha(0.75)
                self.ROI_regions_3[0].set_alpha(0.75)
                self.highlighted_regions.append(0)

        else:
            for i, region in enumerate(self.selected_regions_of_interest):
                if x0 >= region[0] and x0 <= region[1]:
                    self.ROI_regions[i].set_alpha(0.75)
                    self.ROI_regions_2[i].set_alpha(0.75)
                    self.ROI_regions_3[i].set_alpha(0.75)
                    self.highlighted_regions.append(i)
            for i, region in enumerate(self.selected_regions_of_interest):
                if i not in self.highlighted_regions:
                    self.ROI_regions[i].set_alpha(0.2)
                    self.ROI_regions_2[i].set_alpha(0.2)
                    self.ROI_regions_3[i].set_alpha(0.2)

        self.UpdateDiffusionFrame()

    def OnPressROI(self, event):
        if (
            event.inaxes == self.ax_diffusion
            or event.inaxes == self.ax_diffusion_whole_fit
            or event.inaxes == self.ax_diffusion_I0_whole_fit
            or event.inaxes == self.diffusion_coefficient_plot
        ):
            self.press = True
            self.x0 = event.xdata

    def OnMoveROI(self, event):

        if (
            event.inaxes == self.ax_diffusion
            or event.inaxes == self.ax_diffusion_whole_fit
            or event.inaxes == self.ax_diffusion_I0_whole_fit
            or event.inaxes == self.diffusion_coefficient_plot
        ):
            self.move_ROI(event)

    def move_ROI(self, event):
        if self.press:
            self.move = True
            self.x1 = event.xdata
            if self.x1 > self.x0:
                xmax = self.x1
                xmin = self.x0
            else:
                xmax = self.x0
                xmin = self.x1

            self.ROI_regions[-1].set_x(xmin)
            self.ROI_regions[-1].set_width(xmax - xmin)
            self.ROI_regions_2[-1].set_x(xmin)
            self.ROI_regions_2[-1].set_width(xmax - xmin)
            self.ROI_regions_3[-1].set_x(xmin)
            self.ROI_regions_3[-1].set_width(xmax - xmin)
            # self.ROI_regions[-1].set_xy([[xmin,0],[xmin,1],[xmax,1],[xmax,0]])
            # self.ROI_regions_2[-1].set_xy([[xmin,0],[xmin,1],[xmax,1],[xmax,0]])
            # self.ROI_regions_3[-1].set_xy([[xmin,0],[xmin,1],[xmax,1],[xmax,0]])
            self.UpdateDiffusionFrame()

    def release_ROI(self, event):
        if self.press:
            self.x2 = event.xdata
            if self.x2 > self.x0:
                xmax = self.x2
                xmin = self.x0
            else:
                xmax = self.x0
                xmin = self.x2
            self.ROI_x_initial = xmin
            self.ROI_x_final = xmax

            self.ROI_regions[-1].set_x(xmin)
            self.ROI_regions[-1].set_width(xmax - xmin)
            self.ROI_regions_2[-1].set_x(xmin)
            self.ROI_regions_2[-1].set_width(xmax - xmin)
            self.ROI_regions_3[-1].set_x(xmin)
            self.ROI_regions_3[-1].set_width(xmax - xmin)
            # self.ROI_regions[-1].set_xy([[xmin,0],[xmin,1],[xmax,1],[xmax,0]])
            # self.ROI_regions_2[-1].set_xy([[xmin,0],[xmin,1],[xmax,1],[xmax,0]])
            # self.ROI_regions_3[-1].set_xy([[xmin,0],[xmin,1],[xmax,1],[xmax,0]])

            # Add the min and max ppm values to the array of selected regions of interest
            self.selected_regions_of_interest.append([xmin, xmax])

            self.UpdateDiffusionFrame()

            self.press = False
            self.move = False
            self.canvas_diffusion.mpl_disconnect(self.select_ROI_press)
            self.canvas_diffusion.mpl_disconnect(self.select_ROI_move)
            self.canvas_diffusion.mpl_disconnect(self.select_ROI_release)

    def OnReleaseROI(self, event):
        if (
            event.inaxes == self.ax_diffusion
            or event.inaxes == self.ax_diffusion_whole_fit
            or event.inaxes == self.ax_diffusion_I0_whole_fit
            or event.inaxes == self.diffusion_coefficient_plot
        ):
            self.release_ROI(event)

    def OnRegionFitting(self, event):
        # Remove the deleted slices from the gradient and gradient percentages lists

        self.ppms_in_ROI_total = []
        self.ppms_in_ROI_indices_total = []
        self.average_y_data_in_ROI_above_noise_total = []
        self.error_y_data_in_ROI_above_noise_total = []
        self.error_I_I0_in_ROI_total = []
        self.error_log_I_I0_in_ROI_total = []
        self.error_I_I0_in_ROI_total = []
        self.I0_average_in_ROI_total = []
        self.fitted_D_ROI_total = []
        self.fitted_I0_ROI_total = []
        self.mean_fitted_D_ROI_total = []
        self.mean_fitted_I0_ROI_total = []
        self.fitted_I0_total = []
        self.fitted_D_total = []

        for i, region in enumerate(self.selected_regions_of_interest):
            self.ROI_x_initial = region[0]
            self.ROI_x_final = region[1]
            # Get the indices of the ppms which are in the ROI and have intensity above the noise in all slices
            self.ppms_in_ROI = []
            self.ppms_in_ROI_indices = []
            for i, ppm in enumerate(self.ppms_above_noise):
                if ppm >= self.ROI_x_initial and ppm <= self.ROI_x_final:
                    self.ppms_in_ROI.append(ppm)
                    self.ppms_in_ROI_indices.append(i)

            self.average_y_data_in_ROI_above_noise = []
            self.error_y_data_in_ROI_above_noise = []
            self.error_I_I0_in_ROI = []

            for i, data in enumerate(self.y_data_above_noise):
                self.y_data_in_ROI_above_noise_slice = []
                self.I_I0_in_ROI_slice = []
                for index in self.ppms_in_ROI_indices:
                    self.y_data_in_ROI_above_noise_slice.append(np.real(data[index]))
                    self.I_I0_in_ROI_slice.append(
                        np.real(data[index] / self.fitted_I0_global[index])
                    )

                self.average_y_data_in_ROI_above_noise.append(
                    np.mean(np.array(self.y_data_in_ROI_above_noise_slice))
                )
                self.error_y_data_in_ROI_above_noise.append(
                    np.std(np.array(self.y_data_in_ROI_above_noise_slice))
                )
                self.error_I_I0_in_ROI.append(np.std(np.array(self.I_I0_in_ROI_slice)))

            self.average_y_data_in_ROI_above_noise = np.array(
                self.average_y_data_in_ROI_above_noise
            )
            self.error_y_data_in_ROI_above_noise = np.array(
                self.error_y_data_in_ROI_above_noise
            )
            self.error_I_I0_in_ROI = np.array(self.error_I_I0_in_ROI)

            # Also need the error in I/I0 for each slice in the ROI
            self.error_log_I_I0_in_ROI = []
            self.error_I_I0_in_ROI = []
            self.I0_average_in_ROI = []
            for i, slice_I0 in enumerate(self.y_data_above_noise):
                self.I_I0_in_ROI_slice = []
                self.I0_average_in_ROI_slice = []
                for index in self.ppms_in_ROI_indices:
                    self.I_I0_in_ROI_slice.append(
                        np.real(slice_I0[index] / self.fitted_I0_global[index])
                    )
                    self.I0_average_in_ROI_slice.append(
                        np.real(self.fitted_I0_global[index])
                    )

                self.error_log_I_I0_in_ROI.append(
                    np.std(np.log(np.array(self.I_I0_in_ROI_slice)))
                )
                self.error_I_I0_in_ROI.append(np.std(np.array(self.I_I0_in_ROI_slice)))
                self.I0_average_in_ROI.append(
                    np.mean(np.array(self.I0_average_in_ROI_slice))
                )

            self.error_log_I_I0_in_ROI = np.array(self.error_log_I_I0_in_ROI)
            self.I0_average_in_ROI = np.array(self.I0_average_in_ROI)
            self.error_I_I0_in_ROI = np.array(self.error_I_I0_in_ROI)

            self.fitted_D_ROI = []
            self.fitted_I0_ROI = []

            for index in self.ppms_in_ROI_indices:
                self.fitted_D_ROI.append(self.fitted_D_global[index])
                self.fitted_I0_ROI.append(self.fitted_I0_global[index])

            self.mean_fitted_D_ROI = np.mean(np.array(self.fitted_D_ROI))
            self.mean_fitted_I0_ROI = np.mean(np.array(self.fitted_I0_ROI))

            self.ppms_in_ROI_total.append(self.ppms_in_ROI)
            self.ppms_in_ROI_indices_total.append(self.ppms_in_ROI_indices)
            self.average_y_data_in_ROI_above_noise_total.append(
                self.average_y_data_in_ROI_above_noise
            )
            self.error_y_data_in_ROI_above_noise_total.append(
                self.error_y_data_in_ROI_above_noise
            )
            self.error_I_I0_in_ROI_total.append(self.error_I_I0_in_ROI)
            self.error_log_I_I0_in_ROI_total.append(self.error_log_I_I0_in_ROI)
            self.error_I_I0_in_ROI_total.append(self.error_I_I0_in_ROI)
            self.I0_average_in_ROI_total.append(self.I0_average_in_ROI)
            self.fitted_D_ROI_total.append(self.fitted_D_ROI)
            self.fitted_I0_ROI_total.append(self.fitted_I0_ROI)
            self.mean_fitted_D_ROI_total.append(self.mean_fitted_D_ROI)
            self.mean_fitted_I0_ROI_total.append(self.mean_fitted_I0_ROI)

            # Fit the Stejskal Tanner equation to the data for all points in the ROI, use the standard deviation of all I/I0 values as the error
            self.fitted_I0, self.fitted_D = self.leastsq_ROI(
                [np.max(self.average_y_data_in_ROI_above_noise), 1e-9]
            )

            self.fitted_I0_total.append(self.fitted_I0)
            self.fitted_D_total.append(self.fitted_D)

        self.monoexponential_fit = True

        self.PlotRegionFitting()

    def OnBiexponentialFitting(self, event):
        if self.monoexponential_fit != True:
            # Give an error message to say please perform monoexponential fitting first
            msg = wx.MessageDialog(
                self,
                "Please perform monoexponential fitting first",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        elif len(self.selected_regions_of_interest) > 1:
            # Give an error message saying that biexponential fitting is only supported while one region of interest is present
            msg = wx.MessageDialog(
                self,
                "Biexponential fitting is only supported while one region of interest is present",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        else:
            # Perform a biexponential fit
            # Loop over various initial guesses for the biexponential fit
            chi_squared = 100000
            D1_array = 10 ** np.linspace(-5, -10, 5)
            D2_array = 10 ** np.linspace(-5, -10, 5)
            f1_array = np.linspace(0.1, 0.9, 5)
            for i, D1 in enumerate(D1_array):
                for j, D2 in enumerate(D2_array):
                    for k, f1 in enumerate(f1_array):
                        p0 = [
                            np.max(self.average_y_data_in_ROI_above_noise),
                            D1,
                            D2,
                            f1,
                        ]
                        fit = leastsq(self.chi_biexponential_ROI, p0)
                        if (
                            np.sum(self.chi_biexponential_ROI(fit[0]) ** 2)
                            < chi_squared
                        ):
                            chi_squared = np.sum(
                                self.chi_biexponential_ROI(fit[0]) ** 2
                            )
                            best_fit = fit
            fit = best_fit
            I0_ROI = np.abs(fit[0][0])
            D_ROI = np.abs(fit[0][1])
            D2_ROI = np.abs(fit[0][2])
            f1_ROI = np.abs(fit[0][3])
            xvals = np.linspace(0, 1, 100)
            gradient_vals = self.gradients
            self.gradients = np.sqrt(xvals) * self.max_gradient
            # Put the diffusion coefficient into the correct format in the legend
            self.dval1, self.dpow1 = "{:.3e}".format(D_ROI).split("e-")
            if self.dpow1[0] == "0":
                self.dpow1 = self.dpow1[1:]
            self.dval2, self.dpow2 = "{:.3e}".format(D2_ROI).split("e-")
            if self.dpow2[0] == "0":
                self.dpow2 = self.dpow2[1:]

            # Plot the biexponential fit
            self.ax_diffusion_fit.plot(
                xvals,
                self.StejsktalTannerBiexponential([I0_ROI, D_ROI, D2_ROI, f1_ROI])
                / I0_ROI,
                color="tab:red",
                linestyle="--",
                label=r"D$_1$ = "
                + self.dval1
                + r"$\times$10$^{-"
                + r"{}".format(self.dpow1)
                + r"}$ cm$^2$/s, "
                + r"D$_2$ = "
                + self.dval2
                + r"$\times$10$^{-"
                + r"{}".format(self.dpow2)
                + r"}$ cm$^2$/s",
            )
            legend = self.ax_diffusion_fit.legend(fontsize=8)
            legend.get_title().set_color(self.titlecolor)
            self.gradients = gradient_vals
            self.UpdateDiffusionFrame()

    def PlotRegionFitting(self):
        # Generate 3 extra plots for the region fitting (I/I0 vs gradient^2 with fitted curve, log(I/I0) vs gradient^2 with fitted curve, histogram of diffusion coefficients within ROI)
        self.fig_diffusion.clear()
        self.fig_diffusion.tight_layout()

        gs = gridspec.GridSpec(2, 3)

        self.ax_diffusion = self.fig_diffusion.add_subplot(gs[0, 0:2])
        self.ax_diffusion_whole_fit = self.fig_diffusion.add_subplot(
            gs[1, 0], sharex=self.ax_diffusion
        )
        self.ax_diffusion_I0_whole_fit = self.fig_diffusion.add_subplot(
            gs[1, 1], sharex=self.ax_diffusion
        )
        self.ax_diffusion_fit = self.fig_diffusion.add_subplot(gs[0, 2])
        self.ax_diffusion_histogram = self.fig_diffusion.add_subplot(gs[1, 2])

        matplotlib.rcParams.update({"font.size": 8})

        count = 1
        self.slice_plots = []
        for i, data in enumerate(self.y_data):
            (line,) = self.ax_diffusion.plot(
                self.x_data, data, linewidth=0.5, label=str(count)
            )
            self.slice_plots.append(line)
            count += 1
        self.ax_diffusion.set_xlim([self.x_data[0], self.x_data[-1]])

        self.ax_diffusion.set_xlabel(self.main_frame.nmrdata.axislabels[1])
        legend = self.ax_diffusion.legend(title="Slice Number", fontsize=8)
        legend.get_title().set_color(self.titlecolor)
        self.ax_diffusion.set_ylabel("Intensity", fontsize=8)
        self.ax_diffusion.set_title(
            "Diffusion Data", color=self.titlecolor, fontsize=10
        )
        self.ax_diffusion.tick_params(axis="both", which="major", labelsize=8)

        self.noise_region = self.ax_diffusion.axvspan(
            self.noise_x_initial, self.noise_x_final, alpha=0.2, color="gray"
        )

        # Plot the fitted diffusion coefficients and use a twiny to also plot the initial slice of the spectrum
        self.ax_diffusion_whole_fit.plot(self.x_data, self.y_data[0])
        self.ax_diffusion_whole_fit.set_xlabel(
            self.main_frame.nmrdata.axislabels[1], fontsize=8
        )
        self.ax_diffusion_whole_fit.tick_params(axis="both", which="major", labelsize=8)
        self.ax_diffusion_whole_fit.set_yticks([])
        self.diffusion_coefficient_plot = self.ax_diffusion_whole_fit.twinx()
        self.diffusion_coefficient_plot.set_ylabel(r"Diffusion Coefficient (cm$^2$/s)")
        self.diffusion_coefficient_plot.set_xlabel(
            self.main_frame.nmrdata.axislabels[1], fontsize=8
        )
        self.diffusion_coefficient_plot.set_xlim([self.x_data[0], self.x_data[-1]])
        self.diffusion_coefficient_plot.scatter(
            self.ppms_above_noise, self.fitted_D_global, color="tab:red", s=0.5
        )
        self.diffusion_coefficient_plot.set_yscale("log")
        self.diffusion_coefficient_plot.yaxis.tick_left()
        self.diffusion_coefficient_plot.yaxis.set_label_position("left")

        self.diffusion_coefficient_plot.set_title(
            "Diffusion Coefficient vs PPM", color=self.titlecolor, fontsize=10
        )
        self.diffusion_coefficient_plot.tick_params(
            axis="both", which="major", labelsize=8
        )
        self.noise_region_2 = self.ax_diffusion_whole_fit.axvspan(
            self.noise_x_initial, self.noise_x_final, alpha=0.2, color="gray"
        )

        # Plot I/I0 for every chosen ppm across the spectrum for all slices
        for i, selected_y_data in enumerate(self.y_data_above_noise):
            self.ax_diffusion_I0_whole_fit.scatter(
                self.ppms_above_noise,
                np.array(selected_y_data) / self.fitted_I0_global,
                s=0.5,
            )

        self.ax_diffusion_I0_whole_fit.set_xlabel(
            self.main_frame.nmrdata.axislabels[1], fontsize=8
        )
        self.ax_diffusion_I0_whole_fit.set_ylabel(r"I/I$_0$", fontsize=8)
        self.ax_diffusion_I0_whole_fit.set_xlim([self.x_data[0], self.x_data[-1]])
        self.ax_diffusion_I0_whole_fit.set_title(
            r"I/I$_0$ vs PPM", color=self.titlecolor, fontsize=10
        )
        self.ax_diffusion_I0_whole_fit.tick_params(
            axis="both", which="major", labelsize=8
        )
        self.noise_region_3 = self.ax_diffusion_I0_whole_fit.axvspan(
            self.noise_x_initial, self.noise_x_final, alpha=0.2, color="gray"
        )

        # Plot the ROI regions on all plots
        self.ROI_regions = []
        self.ROI_regions_2 = []
        self.ROI_regions_3 = []
        for i, region in enumerate(self.selected_regions_of_interest):
            self.ROI_x_initial = region[0]
            self.ROI_x_final = region[1]
            color = self.ROI_color[i]
            self.ROI_regions.append(
                self.ax_diffusion.axvspan(
                    self.ROI_x_initial, self.ROI_x_final, alpha=0.2, color=color
                )
            )
            self.ROI_regions_2.append(
                self.ax_diffusion_whole_fit.axvspan(
                    self.ROI_x_initial, self.ROI_x_final, alpha=0.2, color=color
                )
            )
            self.ROI_regions_3.append(
                self.ax_diffusion_I0_whole_fit.axvspan(
                    self.ROI_x_initial, self.ROI_x_final, alpha=0.2, color=color
                )
            )

        self.fitted_gaussian_parameters = []
        # Plot a histogram of the diffusion coefficients in the ROI
        # try:
        for i, fitted_D_ROI in enumerate(self.fitted_D_ROI_total):
            if int(len(fitted_D_ROI)) > 0:
                self.ax_diffusion_histogram.hist(
                    fitted_D_ROI,
                    bins=int(len(fitted_D_ROI)),
                    color=self.main_frame.colours[i],
                    edgecolor=self.ROI_color[i],
                    alpha=0.25,
                )
                self.ax_diffusion_histogram.set_xlabel(
                    r"Diffusion Coefficient (cm$^2$/s)", fontsize=8
                )
                self.ax_diffusion_histogram.set_ylabel("Frequency Density", fontsize=8)
                self.ax_diffusion_histogram.set_title(
                    "Histogram of Diffusion Coefficients",
                    color=self.titlecolor,
                    fontsize=10,
                )
                # Get the bin size of the histogram
                self.bin_size = self.ax_diffusion_histogram.patches[0].get_width()
                self.bin_centers = np.arange(
                    min(fitted_D_ROI) + self.bin_size / 2,
                    max(fitted_D_ROI),
                    self.bin_size,
                )
                self.bin_centers = np.array(self.bin_centers)
                self.bin_centers = self.bin_centers[
                    np.where(self.bin_centers <= max(fitted_D_ROI))
                ]
                self.bin_centers = self.bin_centers[
                    np.where(self.bin_centers >= min(fitted_D_ROI))
                ]
                self.bin_centers = np.array(self.bin_centers)
                # Get the frequency densities of the histogram in each bin
                self.frequency_density = []
                for j, bin_center in enumerate(self.bin_centers):
                    self.frequency_density.append(
                        len(
                            np.where(
                                (fitted_D_ROI >= bin_center - self.bin_size / 2)
                                & (fitted_D_ROI < bin_center + self.bin_size / 2)
                            )[0]
                        )
                    )
                # Fit a gaussian to the histogram of diffusion coefficients, this will be the error in the diffusion coefficient
                self.fitted_D_ROI = np.array(fitted_D_ROI)
                result = self.leastsq_gaussian_ROI(
                    [1, np.mean(fitted_D_ROI), np.std(fitted_D_ROI)]
                )

                if result[0] != "Failed":
                    A, mu, sigma = result
                    self.fitted_gaussian_parameters.append([A, mu, sigma])
                    if np.abs(sigma) > 1.25 * np.std(fitted_D_ROI):
                        self.ax_diffusion_histogram.plot(
                            self.bin_centers,
                            self.gaussian_ROI(self.bin_centers, A, mu, sigma),
                            label=r"$\sigma_{gauss}$ = "
                            + "{:.3e}".format(np.abs(sigma))
                            + r", $\sigma_{stdev}$ = "
                            + "{:.3e}".format(np.std(fitted_D_ROI)),
                            color=self.ROI_color[i],
                        )
                    else:
                        self.ax_diffusion_histogram.plot(
                            self.bin_centers,
                            self.gaussian_ROI(self.bin_centers, A, mu, sigma),
                            label=r"$\sigma_{gauss}$ = "
                            + "{:.3e}".format(np.abs(sigma)),
                            color=self.ROI_color[i],
                        )
                    legend = self.ax_diffusion_histogram.legend(fontsize=8)
                    legend.get_title().set_color(self.titlecolor)
                else:
                    # Gaussian fit failed - setting sigma to the standard deviation of the diffusion coefficients and A to max frequency density
                    sigma = np.std(fitted_D_ROI)
                    A = max(self.frequency_density)
                    self.fitted_gaussian_parameters.append(
                        [A, np.mean(fitted_D_ROI), sigma]
                    )

                    # # Give an error message saying that one of the ROI windows is too small. Please increase the size of the ROI window and try again
                    msg = wx.MessageDialog(
                        self,
                        "Gaussian error fit did not work, error set to standard deviation. If desired, increase the size of the ROI window and try again",
                        "Error",
                        wx.OK | wx.ICON_ERROR,
                    )
                    msg.ShowModal()
                    msg.Destroy()
            else:
                # Give an error message saying that one of the ROI windows is too small. Please delete the ROI and try again
                msg = wx.MessageDialog(
                    self,
                    "One of the ROI windows is too small. Please delete the ROI and try again",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                return

        # excpt:
        #     pass
        # # Give an error message saying that one of the ROI windows is too small. Please increase the size of the ROI window and try again
        # msg = wx.MessageDialog(self, 'One of the ROI windows is too small. Please increase the size of the ROI window and try again', 'Error', wx.OK | wx.ICON_ERROR)
        # msg.ShowModal()
        # msg.Destroy()
        # return

        for i, region in enumerate(self.selected_regions_of_interest):
            self.ROI_x_initial = region[0]
            self.ROI_x_final = region[1]
            self.mean_fitted_D_ROI = self.mean_fitted_D_ROI_total[i]
            self.average_y_data_in_ROI_above_noise = (
                self.average_y_data_in_ROI_above_noise_total[i]
            )
            self.error_y_data_in_ROI_above_noise = (
                self.error_y_data_in_ROI_above_noise_total[i]
            )
            self.error_I_I0_in_ROI = self.error_I_I0_in_ROI_total[i]
            self.I0_average_in_ROI = self.I0_average_in_ROI_total[i]
            self.fitted_D_ROI = self.fitted_D_ROI_total[i]
            self.fitted_I0_ROI = self.fitted_I0_ROI_total[i]
            self.mean_fitted_I0_ROI = self.mean_fitted_I0_ROI_total[i]
            self.fitted_I0 = self.fitted_I0_total[i]
            self.fitted_D = self.fitted_D_total[i]
            self.error_log_I_I0_in_ROI = self.error_log_I_I0_in_ROI_total[i]
            self.error_I_I0_in_ROI = self.error_I_I0_in_ROI_total[i]

            # Put the diffusion coefficient into the correct format in the legend
            self.dval, self.dpow = "{:.3e}".format(self.mean_fitted_D_ROI).split("e-")
            if self.dpow[0] == "0":
                self.dpow = self.dpow[1:]

            if (self.fitted_gaussian_parameters[i][2]) > 1.25 * np.std(
                self.fitted_D_ROI
            ):
                self.errorval, self.errorpow = "{:.3e}".format(
                    np.abs(np.std(self.fitted_D_ROI))
                ).split("e-")
            else:
                self.errorval, self.errorpow = "{:.3e}".format(
                    np.abs(self.fitted_gaussian_parameters[i][2])
                ).split("e-")
            if self.errorpow[0] == "0":
                self.errorpow = self.errorpow[1:]
            difference = int(self.errorpow) - int(self.dpow)
            if difference > 0:
                # add zeros to the front of self.errorval
                self.errorval = (
                    "0."
                    + "0" * (difference - 1)
                    + self.errorval.split(".")[0]
                    + self.errorval.split(".")[1]
                )

            # Plot the fitted curve for the ROI data
            self.ax_diffusion_fit.errorbar(
                (np.array(self.gradients_percent) / 100) ** 2,
                self.average_y_data_in_ROI_above_noise / self.I0_average_in_ROI,
                yerr=self.error_I_I0_in_ROI,
                fmt="o",
                markersize=1,
                capsize=2,
                color=self.ROI_color[i],
            )

            xvals = np.linspace(0, 1, 100)
            gradient_vals = self.gradients
            self.gradients = np.sqrt(xvals) * self.max_gradient
            self.ax_diffusion_fit.plot(
                xvals,
                self.StejsktalTanner([self.mean_fitted_I0_ROI, self.mean_fitted_D_ROI])
                / self.mean_fitted_I0_ROI,
                label=r"D = ({}".format(self.dval)
                + r"$\pm$"
                + r"{})".format(self.errorval)
                + r"$\times$10$^{-"
                + r"{}".format(self.dpow)
                + r"}$ cm$^2$/s",
                color=self.ROI_color[i],
            )

            self.gradients = gradient_vals
        self.ax_diffusion_fit.set_xlabel(r"(G/G$_{max}$)$^2$", fontsize=8)
        self.ax_diffusion_fit.set_ylabel(r"I/I$_0$", fontsize=8)
        self.ax_diffusion_fit.set_title(
            "Fitted Stejskal Tanner", color=self.titlecolor, fontsize=10
        )
        legend = self.ax_diffusion_fit.legend(fontsize=8)
        legend.get_title().set_color(self.titlecolor)

        self.fig_diffusion.tight_layout()

        self.UpdateDiffusionFrame()

    def StejsktalTannerBiexponential(self, p0):
        I0, D1, D2, f1 = p0
        # Ensure all values are positive
        D1 = np.abs(D1)
        D2 = np.abs(D2)
        f1 = np.abs(f1)
        I0 = np.abs(I0)
        return I0 * (
            f1
            * np.exp(
                -(
                    (self.gamma**2)
                    * (self.gradients**2)
                    * (self.little_delta * 1e-6) ** 2
                )
                * (self.big_delta - (self.little_delta * 1e-6) / 3)
                * D1
            )
            + (1 - f1)
            * np.exp(
                -(
                    (self.gamma**2)
                    * (self.gradients**2)
                    * (self.little_delta * 1e-6) ** 2
                )
                * (self.big_delta - (self.little_delta * 1e-6) / 3)
                * D2
            )
        )

    def chi_biexponential_ROI(self, p0):
        return (
            self.average_y_data_in_ROI_above_noise
            - self.StejsktalTannerBiexponential(p0)
        ) / self.error_y_data_in_ROI_above_noise

    def leastsq_biexponential(self, p0):
        fit = leastsq(self.chi_biexponential_ROI, p0)
        return fit[0]

    def leastsq_ROI(self, p0):
        fit = leastsq(self.chi_ROI, p0)
        return fit[0]

    def chi_ROI(self, p0):
        return (
            self.average_y_data_in_ROI_above_noise - self.StejsktalTanner(p0)
        ) / self.error_y_data_in_ROI_above_noise

    def gaussian_ROI(self, x, A, mu, sigma):
        return A * np.exp(-((x - mu) ** 2) / (2 * sigma**2))

    def chi_gaussian_ROI(self, p0):
        return self.frequency_density - self.gaussian_ROI(
            self.bin_centers, p0[0], p0[1], p0[2]
        )

    def leastsq_gaussian_ROI(self, p0):
        try:
            fit = leastsq(self.chi_gaussian_ROI, p0)
            return fit[0]
        except:
            return ["Failed"]

    def OnDeleteSlice(self, event):

        # Check to see if the gradients have already been inputted
        try:
            self.gradients
        except:
            # Give an error message saying that the gradients must be inputted first
            msg = wx.MessageDialog(
                self,
                "Please input the gradients first before deleting a slice",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        # Check to see that the number of slices is four or more
        if len(self.y_data) < 4:
            # Give an error message saying that there must be at least four slices
            msg = wx.MessageDialog(
                self,
                "There must be at least four slices in the data to perform diffusion data fitting",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        self.delete_slice = True
        # Bring up a dialog box to ask which slice to delete
        self.delete_slice_index = 0
        self.delete_slice_dialog = DeleteSliceDialog("Delete Slice", self)

    # Function to continue the deletion of a slice after the dialog box has been closed and completed by the user
    def continue_deletion(self):
        # Check to see if the full spectrum fitting has been performed
        if self.whole_plot != True:
            # Delete the correct slice in the y data
            self.y_data = np.delete(self.y_data, self.delete_slice_index, axis=0)

            # Delete the correct value in gradients
            self.gradients = np.delete(self.gradients, self.delete_slice_index)
            self.gradients_percent = np.delete(
                self.gradients_percent, self.delete_slice_index
            )

            # Redo the plotting
            self.fig_diffusion.clear()
            self.fig_diffusion.tight_layout()
            self.plot_diffusion_data()

            # If the noise region has already been selected, then redo the plotting of this
            if self.noise_region_selection == True:
                self.noise_region = self.ax_diffusion.axvspan(
                    self.noise_x_initial, self.noise_x_final, alpha=0.2, color="gray"
                )
            self.UpdateDiffusionFrame()
        elif self.whole_plot == True and self.monoexponential_fit != True:
            # Delete the correct slice in the y data
            self.y_data = np.delete(self.y_data, self.delete_slice_index, axis=0)

            # Delete the correct value in gradients
            self.gradients = np.delete(self.gradients, self.delete_slice_index)
            self.gradients_percent = np.delete(
                self.gradients_percent, self.delete_slice_index
            )

            self.OnWholeSpectrumFitting(event=None)

        elif self.whole_plot == True and self.monoexponential_fit == True:
            # Delete the correct slice in the y data
            self.y_data = np.delete(self.y_data, self.delete_slice_index, axis=0)

            # Delete the correct value in gradients
            self.gradients = np.delete(self.gradients, self.delete_slice_index)
            self.gradients_percent = np.delete(
                self.gradients_percent, self.delete_slice_index
            )

            self.OnWholeSpectrumFitting(event=None)
            self.OnRegionFitting(event=None)

    def input_gradients_text_button(self, event):
        self.input_gradients_text_dialog = DiffusionGradientManualInput(
            "Input Gradients", self, self.spectrometer
        )


# Need this box to delete any slice from the data, as well as remove it from the list of gradients
class DeleteSliceDialog(wx.Frame):
    def __init__(self, title, parent):
        self.main_frame = parent
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        width = int(0.2 * self.monitorWidth)
        height = int(0.1 * self.monitorHeight)
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.panel_delete_slice = wx.Panel(self, -1)
        self.main_delete_slice = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_delete_slice)

        self.make_delete_slice_sizer()
        self.Show()

    def make_delete_slice_sizer(self):
        # Make a sizer to hold the text box and button
        self.delete_slice_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.delete_slice_sizer.AddSpacer(5)
        # ComboBox for slice number
        self.slice_number_label = wx.StaticText(self, -1, "Slice Number:")
        self.delete_slice_sizer.Add(self.slice_number_label)
        self.delete_slice_sizer.AddSpacer(5)
        self.slice_number_choices = []
        for i in range(len(self.main_frame.y_data)):
            self.slice_number_choices.append(str(i + 1))
        self.slice_number_combobox = wx.ComboBox(
            self, -1, choices=self.slice_number_choices, style=wx.CB_READONLY
        )
        self.delete_slice_sizer.Add(self.slice_number_combobox)
        self.delete_slice_sizer.AddSpacer(5)
        # Have a button to confirm the deletion
        self.confirm_button = wx.Button(self, -1, "Delete")
        self.confirm_button.Bind(wx.EVT_BUTTON, self.OnConfirmDelete)
        self.delete_slice_sizer.Add(self.confirm_button)
        self.delete_slice_sizer.AddSpacer(5)

        self.main_delete_slice.AddSpacer(5)
        self.main_delete_slice.Add(self.delete_slice_sizer)
        self.main_delete_slice.AddSpacer(5)

    def OnConfirmDelete(self, event):
        self.main_frame.delete_slice_index = (
            int(self.slice_number_combobox.GetValue()) - 1
        )
        self.main_frame.deleted_slices.append(self.main_frame.delete_slice_index)
        self.main_frame.continue_deletion()
        self.Destroy()


class GradientsManualTextInput(wx.Frame):
    def __init__(self, title, parent):
        self.main_frame = parent
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        width = int(0.3 * self.monitorWidth)
        height = int(0.5 * self.monitorHeight)
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.panel_gradients_text_input = wx.Panel(self, -1)
        self.main_gradients_text_input = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_gradients_text_input)

        self.make_manual_gradients_text_input_sizer()
        self.Show()

    def make_manual_gradients_text_input_sizer(self):
        if self.main_frame.spectrometer == "Bruker":
            self.input_gradient_text_label = wx.StaticBox(
                self, -1, "Input gradient percentages (one per line)"
            )
            self.input_gradient_text_sizer = wx.StaticBoxSizer(
                self.input_gradient_text_label, wx.VERTICAL
            )
        else:
            self.input_gradient_text_label = wx.StaticBox(
                self, -1, "Input gradient DAC values (one per line)"
            )
            self.input_gradient_text_sizer = wx.StaticBoxSizer(
                self.input_gradient_text_label, wx.VERTICAL
            )

        try:
            file = open("gradients.txt")
            gradients_percent = []
            for line in file.readlines():
                line = line.split("\n")[0]
                gradients_percent.append(line)

            label = ""
            for gradient in gradients_percent:
                label = label + gradient + "\n"
            file.close()
        except:
            label = ""

        self.gradient_box = wx.TextCtrl(
            self, -1, value=label, size=(250, 400), style=wx.TE_MULTILINE
        )
        self.input_gradient_text_sizer.AddSpacer(3)
        self.input_gradient_text_sizer.Add(self.gradient_box)
        self.input_gradient_text_sizer.AddSpacer(10)

        self.save_gradients_button = wx.Button(self, -1, "Save gradients")

        self.save_gradients_button.Bind(wx.EVT_BUTTON, self.OnSaveGradients)

        self.input_gradient_text_sizer.Add(self.save_gradients_button)

        self.main_gradients_text_input.Add(self.input_gradient_text_sizer)

    def OnSaveGradients(self, event):
        # Remove all extra empty lines at the end of the text box
        if self.gradient_box.GetValue().split("\n")[-1] == "":
            old_value = self.gradient_box.GetValue()
            new_value = ""
            for i, line in enumerate(old_value.split("\n")):
                if line != "":
                    if i == len(old_value.split("\n")) - 1:
                        new_value = new_value + line.rstrip()
                    else:
                        new_value = new_value + line.rstrip() + "\n"
            # If the last element in the string contains a newline character, remove it
            if new_value[-1] == "\n":
                new_value = new_value[:-1]
            self.gradient_box.SetValue(new_value)
            self.panel_gradients_text_input.Layout()

        if self.gradient_box.GetValue() == "":
            # Give an error message saying that the gradients must be inputted first
            msg = wx.MessageDialog(
                self,
                "Please input the gradients first before saving",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        # Ensure all text are floats/integers
        for line in self.gradient_box.GetValue().split("\n"):
            if line != "":
                try:
                    float(line)
                except:
                    # Give an error message saying that all lines must be floats
                    error = (
                        "All gradient values must be numbers: "
                        + line
                        + " is not a number"
                    )
                    msg = wx.MessageDialog(self, error, "Error", wx.OK | wx.ICON_ERROR)
                    msg.ShowModal()
                    msg.Destroy()
                    return

        # Ensure there are no negative values
        for line in self.gradient_box.GetValue().split("\n"):
            if line != "":
                if float(line) < 0:
                    # Give an error message saying that there must be no negative values
                    msg = wx.MessageDialog(
                        self,
                        "There must be no negative values",
                        "Error",
                        wx.OK | wx.ICON_ERROR,
                    )
                    msg.ShowModal()
                    msg.Destroy()
                    return
        # If Bruker, ensure that all gradient percentages are between 0 and 100
        if self.main_frame.spectrometer == "Bruker":
            for line in self.gradient_box.GetValue().split("\n"):
                if line != "":
                    if float(line) > 100:
                        # Give an error message saying that there must be no values greater than 100
                        error = (
                            "There must be no gradient percentages greater than 100: "
                            + line
                            + " is greater than 100"
                        )
                        msg = wx.MessageDialog(
                            self, error, "Error", wx.OK | wx.ICON_ERROR
                        )
                        msg.ShowModal()
                        msg.Destroy()
                        return
        # If Varian, ensure that all gradient DAC values are between 0 and 30000
        else:
            for line in self.gradient_box.GetValue().split("\n"):
                if line != "":
                    if float(line) > 30000:
                        # Give an error message saying that there must be no values greater than 30000
                        error = (
                            "There must be no DAC gradient values greater than 30000: "
                            + line
                            + " is greater than 30000"
                        )
                        msg = wx.MessageDialog(
                            self, error, "Error", wx.OK | wx.ICON_ERROR
                        )
                        msg.ShowModal()
                        msg.Destroy()
                        return
        # Ensure that the number of gradients entered is the same as the number of slices
        if len(self.gradient_box.GetValue().split("\n")) != len(self.main_frame.y_data):
            # Give an error message saying that the number of gradients must be the same as the number of slices
            error = (
                "The number of gradients must be the same as the number of slices in the data: "
                + str(len(self.gradient_box.GetValue().split("\n")))
                + " gradients entered, "
                + str(len(self.main_frame.y_data))
                + " slices in the data"
            )
            msg = wx.MessageDialog(self, error, "Error", wx.OK | wx.ICON_ERROR)
            msg.ShowModal()
            msg.Destroy()
            return
        file = open("gradients.txt", "w")
        file.write(self.gradient_box.GetValue())
        file.close()
        if self.main_frame.spectrometer == "Bruker":
            self.main_frame.gradients_percent = []
            self.main_frame.gradients = []
            for line in self.gradient_box.GetValue().split("\n"):
                if line != "":
                    self.main_frame.gradients_percent.append(float(line))

            self.main_frame.gradients_percent = np.array(
                self.main_frame.gradients_percent
            )
            self.main_frame.gradients = (
                (self.main_frame.gradients_percent / 100)
                * self.main_frame.max_gradient
                * self.main_frame.integral_factor
            )

        else:
            self.gradients_DAC = []
            for line in self.gradient_box.GetValue().split("\n"):
                if line != "":
                    self.gradients_DAC.append(float(line))
            self.gradients_DAC = np.array(self.gradients_DAC)

            self.gradients = (
                np.array(self.gradients_DAC) * self.main_frame.DAC_conversion
            )
            self.gradients_percent = (
                np.array(self.gradients_DAC)
                * self.main_frame.DAC_conversion
                / self.main_frame.max_gradient
                * 100
            )
            self.main_frame.gradients_percent = self.gradients_percent
            self.main_frame.gradients = self.gradients


class DiffusionGradientManualInput(wx.Frame):
    def __init__(self, title, parent, spectrometer="Bruker"):
        self.main_frame = parent
        self.spectrometer = spectrometer
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        width = int(0.3 * self.monitorWidth)
        height = int(0.35 * self.monitorHeight)
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.panel_gradient_input = wx.Panel(self, -1)

        # Define initial default values
        self.number_of_gradients = 5
        self.min_gradient_percent = 10.0
        self.max_gradient_percent = 90.0
        self.min_gradient_DAC = 10000
        self.max_gradient_DAC = 30000
        self.gradient_distribution = "squared"

        self.main_sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.make_manual_gradient_input_sizer()
        self.make_manual_gradients_text_input_sizer()

        self.main_sizer.AddSpacer(5)
        self.main_sizer.Add(self.gradient_input_sizer)
        self.main_sizer.AddSpacer(20)
        self.main_sizer.Add(self.input_gradient_text_sizer)
        self.main_sizer.AddSpacer(5)

        self.SetSizer(self.main_sizer)
        self.Show()

    def make_manual_gradient_input_sizer(self):
        self.gradient_input_sizer_label = wx.StaticBox(self, -1, "Calculate Gradients")
        self.gradient_input_sizer = wx.StaticBoxSizer(
            self.gradient_input_sizer_label, wx.VERTICAL
        )
        self.gradient_input_sizer.AddSpacer(5)
        # TextCtrl for number of gradient values used
        self.number_of_gradients_label = wx.StaticText(self, -1, "Number of Gradients:")
        self.gradient_input_sizer.Add(self.number_of_gradients_label)
        self.gradient_input_sizer.AddSpacer(5)
        self.number_of_gradients = int(len(self.main_frame.y_data))
        self.number_of_gradients_box = wx.TextCtrl(
            self, -1, str(self.number_of_gradients), size=(50, -1)
        )
        self.gradient_input_sizer.Add(self.number_of_gradients_box)
        self.gradient_input_sizer.AddSpacer(5)

        if self.spectrometer == "Bruker":
            # TextCtrl for minimum gradient percentage
            self.min_gradient_label = wx.StaticText(self, -1, "Min Gradient (%):")
            self.gradient_input_sizer.Add(self.min_gradient_label)
            self.gradient_input_sizer.AddSpacer(5)
            self.min_gradient_box = wx.TextCtrl(
                self, -1, str(self.min_gradient_percent), size=(50, -1)
            )
            self.gradient_input_sizer.Add(self.min_gradient_box)
            self.gradient_input_sizer.AddSpacer(5)

            # TextCtrl for maximum gradient percentage
            self.max_gradient_label = wx.StaticText(self, -1, "Max Gradient (%):")
            self.gradient_input_sizer.Add(self.max_gradient_label)
            self.gradient_input_sizer.AddSpacer(5)
            self.max_gradient_box = wx.TextCtrl(
                self, -1, str(self.max_gradient_percent), size=(50, -1)
            )
            self.gradient_input_sizer.Add(self.max_gradient_box)
            self.gradient_input_sizer.AddSpacer(5)
        else:
            # TextCtrl for minimum gradient DAC value
            self.min_gradient_label = wx.StaticText(self, -1, "Min Gradient (DAC):")
            self.gradient_input_sizer.Add(self.min_gradient_label)
            self.gradient_input_sizer.AddSpacer(5)
            self.min_gradient_box = wx.TextCtrl(
                self, -1, str(self.min_gradient_DAC), size=(50, -1)
            )
            self.gradient_input_sizer.Add(self.min_gradient_box)
            self.gradient_input_sizer.AddSpacer(5)

            # TextCtrl for maximum gradient DAC value
            self.max_gradient_label = wx.StaticText(self, -1, "Max Gradient (DAC):")
            self.gradient_input_sizer.Add(self.max_gradient_label)
            self.gradient_input_sizer.AddSpacer(5)
            self.max_gradient_box = wx.TextCtrl(
                self, -1, str(self.max_gradient_DAC), size=(50, -1)
            )
            self.gradient_input_sizer.Add(self.max_gradient_box)
            self.gradient_input_sizer.AddSpacer(5)

        # RadioBox for gradient distribution
        self.gradient_distribution_label = wx.StaticText(self, -1, "Gradient Spacing:")
        self.gradient_input_sizer.Add(self.gradient_distribution_label)
        self.gradient_input_sizer.AddSpacer(5)
        self.gradient_distribution_choices = ["linear", "squared", "exponential"]
        self.gradient_distribution_radiobox = wx.RadioBox(
            self, -1, choices=self.gradient_distribution_choices, style=wx.RA_VERTICAL
        )
        self.gradient_distribution_radiobox.SetSelection(1)
        self.gradient_input_sizer.Add(self.gradient_distribution_radiobox)
        self.gradient_input_sizer.AddSpacer(5)

        # Button to confirm input
        self.confirm_button = wx.Button(self, -1, "Calculate")
        if self.spectrometer == "Bruker":
            self.confirm_button.Bind(wx.EVT_BUTTON, self.OnCalculateGradientsBruker)
        else:
            self.confirm_button.Bind(wx.EVT_BUTTON, self.OnCalculateGradientsVarian)
        self.gradient_input_sizer.Add(self.confirm_button)
        self.gradient_input_sizer.AddSpacer(5)

    def OnCalculateGradientsBruker(self, event):
        # Check number of gradients is an integer
        if self.number_of_gradients_box.GetValue().isdigit() == False:
            msg = wx.MessageDialog(
                self,
                "Number of gradients must be an integer",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return

        if (
            int(self.number_of_gradients_box.GetValue())
            - len(self.main_frame.deleted_slices)
        ) != len(self.main_frame.y_data):
            error = (
                "Number of gradients must be equal to the number of slices: "
                + str(len(self.main_frame.y_data))
                + " slices in the data, "
                + str(int(self.number_of_gradients_box.GetValue()))
                + " gradients entered"
            )
            msg = wx.MessageDialog(self, error, "Error", wx.OK | wx.ICON_ERROR)
            msg.ShowModal()
            msg.Destroy()
            return

        # Check min and max gradient percentages are between 0 and 100
        if (
            float(self.min_gradient_box.GetValue()) < 0
            or float(self.min_gradient_box.GetValue()) > 100
        ):
            msg = wx.MessageDialog(
                self,
                "Min gradient percentage must be between 0 and 100",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        if (
            float(self.max_gradient_box.GetValue()) < 0
            or float(self.max_gradient_box.GetValue()) > 100
        ):
            msg = wx.MessageDialog(
                self,
                "Max gradient percentage must be between 0 and 100",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        # Check min gradient percentage is less than max gradient percentage
        if float(self.min_gradient_box.GetValue()) > float(
            self.max_gradient_box.GetValue()
        ):
            msg = wx.MessageDialog(
                self,
                "Min gradient percentage must be less than max gradient percentage",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        self.number_of_gradients = int(self.number_of_gradients_box.GetValue())
        self.min_gradient_percent = float(self.min_gradient_box.GetValue())
        self.max_gradient_percent = float(self.max_gradient_box.GetValue())
        self.gradient_distribution = self.gradient_distribution_choices[
            self.gradient_distribution_radiobox.GetSelection()
        ]
        self.CalculateGradientsBruker()
        self.gradients = (
            self.gradients_percent
            * self.main_frame.max_gradient
            * self.main_frame.integral_factor
            / 100
        )
        self.main_frame.gradients_percent = self.gradients_percent
        self.main_frame.gradients = self.gradients
        string_of_gradients = ""
        for i, gradient in enumerate(self.gradients_percent):
            if i == len(self.gradients_percent) - 1:
                string_of_gradients = string_of_gradients + str(gradient)
            else:
                string_of_gradients = string_of_gradients + str(gradient) + "\n"

        # Save gradients to text file
        file = open("gradients.txt", "w")
        file.write(string_of_gradients)
        file.close()

        # Write the gradients to the text box
        self.gradient_box.SetValue(string_of_gradients)

        if len(self.main_frame.y_data) != len(self.main_frame.gradients):
            for i, deleted_slice in enumerate(self.main_frame.deleted_slices):
                self.main_frame.gradients = np.delete(
                    self.main_frame.gradients, deleted_slice, axis=0
                )
                self.main_frame.gradients_percent = np.delete(
                    self.main_frame.gradients_percent, deleted_slice, axis=0
                )

    def CalculateGradientsBruker(self):
        if self.gradient_distribution == "linear":
            self.gradients_percent = np.linspace(
                self.min_gradient_percent,
                self.max_gradient_percent,
                self.number_of_gradients,
            )
        elif self.gradient_distribution == "squared":
            self.gradients_percent = np.sqrt(
                np.linspace(
                    self.min_gradient_percent**2,
                    self.max_gradient_percent**2,
                    self.number_of_gradients,
                )
            )
        else:
            self.gradients_percent = np.log(
                np.linspace(
                    np.exp(self.min_gradient_percent),
                    np.exp(self.max_gradient_percent),
                    self.number_of_gradients,
                )
            )

    def OnCalculateGradientsVarian(self, event):
        # Check number of gradients is an integer
        if self.number_of_gradients_box.GetValue().isdigit() == False:
            msg = wx.MessageDialog(
                self,
                "Number of gradients must be an integer",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return

        if (
            int(self.number_of_gradients_box.GetValue())
            - len(self.main_frame.deleted_slices)
        ) != len(self.main_frame.y_data):
            error = (
                "Number of gradients must be equal to the number of slices: "
                + str(len(self.main_frame.y_data))
                + " slices in the data, "
                + str(int(self.number_of_gradients_box.GetValue()))
                + " gradients entered"
            )
            msg = wx.MessageDialog(self, error, "Error", wx.OK | wx.ICON_ERROR)
            msg.ShowModal()
            msg.Destroy()
            return

        # Check min gradient DAC is a number
        if self.min_gradient_box.GetValue().isdigit() == False:
            msg = wx.MessageDialog(
                self,
                "Min gradient DAC value must be a number",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        # Check min gradient DAC is a number above 0
        if float(self.min_gradient_box.GetValue()) <= 0:
            msg = wx.MessageDialog(
                self,
                "Min gradient DAC value must be above 0",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        # Check max gradient DAC is a number
        if self.max_gradient_box.GetValue().isdigit() == False:
            msg = wx.MessageDialog(
                self,
                "Max gradient DAC value must be a number",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        # Check max gradient DAC is a number above 0
        if float(self.max_gradient_box.GetValue()) <= 0:
            msg = wx.MessageDialog(
                self,
                "Max gradient DAC value must be above 0",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return

        # Check min gradient percentage is less than max gradient percentage
        if float(self.min_gradient_box.GetValue()) > float(
            self.max_gradient_box.GetValue()
        ):
            msg = wx.MessageDialog(
                self,
                "Min gradient DAC value must be less than max gradient DAC value",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        self.number_of_gradients = int(self.number_of_gradients_box.GetValue())
        self.min_gradient_DAC = float(self.min_gradient_box.GetValue())
        self.max_gradient_DAC = float(self.max_gradient_box.GetValue())
        self.gradient_distribution = self.gradient_distribution_choices[
            self.gradient_distribution_radiobox.GetSelection()
        ]
        self.CalculateGradientsVarian()
        self.gradients = np.array(self.gradients_DAC) * self.main_frame.DAC_conversion
        self.gradients_percent = (
            np.array(self.gradients_DAC)
            * self.main_frame.DAC_conversion
            / self.main_frame.max_gradient
            * 100
        )
        self.main_frame.gradients_percent = self.gradients_percent
        self.main_frame.gradients = self.gradients
        if len(self.main_frame.y_data) != len(self.main_frame.gradients):
            for i, deleted_slice in enumerate(self.main_frame.deleted_slices):
                self.main_frame.gradients = np.delete(
                    self.main_frame.gradients, deleted_slice, axis=0
                )
                self.main_frame.gradients_percent = np.delete(
                    self.main_frame.gradients_percent, deleted_slice, axis=0
                )

        string_of_gradients = ""
        for i, gradient in enumerate(self.gradients_DAC):
            if i == len(self.gradients_DAC) - 1:
                string_of_gradients = string_of_gradients + str(gradient)
            else:
                string_of_gradients = string_of_gradients + str(gradient) + "\n"

        # Save gradients to text file
        file = open("gradients.txt", "w")
        file.write(string_of_gradients)
        file.close()

        # Write the gradients to the text box
        self.gradient_box.SetValue(string_of_gradients)

    def CalculateGradientsVarian(self):
        if self.gradient_distribution == "linear":
            self.gradients_DAC = np.linspace(
                self.min_gradient_DAC, self.max_gradient_DAC, self.number_of_gradients
            )
        elif self.gradient_distribution == "squared":
            self.gradients_DAC = np.sqrt(
                np.linspace(
                    self.min_gradient_DAC**2,
                    self.max_gradient_DAC**2,
                    self.number_of_gradients,
                )
            )
        else:
            self.gradients_DAC = np.log(
                np.linspace(
                    np.exp(self.min_gradient_DAC),
                    np.exp(self.max_gradient_DAC),
                    self.number_of_gradients,
                )
            )

    def make_manual_gradients_text_input_sizer(self):
        if self.main_frame.spectrometer == "Bruker":
            self.input_gradient_text_label = wx.StaticBox(
                self, -1, "Input gradient percentages (one per line)"
            )
            self.input_gradient_text_sizer = wx.StaticBoxSizer(
                self.input_gradient_text_label, wx.VERTICAL
            )
        else:
            self.input_gradient_text_label = wx.StaticBox(
                self, -1, "Input gradient DAC values (one per line)"
            )
            self.input_gradient_text_sizer = wx.StaticBoxSizer(
                self.input_gradient_text_label, wx.VERTICAL
            )

        try:
            file = open("gradients.txt")
            gradients_percent = []
            for line in file.readlines():
                line = line.split("\n")[0]
                gradients_percent.append(line)

            label = ""
            for gradient in gradients_percent:
                label = label + gradient + "\n"
            file.close()
        except:
            label = ""

        self.gradient_box = wx.TextCtrl(
            self, -1, value=label, size=(250, 235), style=wx.TE_MULTILINE
        )
        self.input_gradient_text_sizer.AddSpacer(3)
        self.input_gradient_text_sizer.Add(self.gradient_box)
        self.input_gradient_text_sizer.AddSpacer(10)

        self.save_gradients_button = wx.Button(self, -1, "Save gradients")

        self.save_gradients_button.Bind(wx.EVT_BUTTON, self.OnSaveGradients)

        self.input_gradient_text_sizer.Add(self.save_gradients_button)

    def OnSaveGradients(self, event):
        # Remove all extra empty lines at the end of the text box
        if self.gradient_box.GetValue().split("\n")[-1] == "":
            old_value = self.gradient_box.GetValue()
            new_value = ""
            for i, line in enumerate(old_value.split("\n")):
                if line != "":
                    if i == len(old_value.split("\n")) - 1:
                        new_value = new_value + line.rstrip()
                    else:
                        new_value = new_value + line.rstrip() + "\n"
            # If the last element in the string contains a newline character, remove it
            if new_value[-1] == "\n":
                new_value = new_value[:-1]
            self.gradient_box.SetValue(new_value)
            self.panel_gradient_input.Layout()

        if self.gradient_box.GetValue() == "":
            # Give an error message saying that the gradients must be inputted first
            msg = wx.MessageDialog(
                self,
                "Please input the gradients first before saving",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        # Ensure all text are floats/integers
        for line in self.gradient_box.GetValue().split("\n"):
            if line != "":
                try:
                    float(line)
                except:
                    # Give an error message saying that all lines must be floats
                    error = (
                        "All gradient values must be numbers: "
                        + line
                        + " is not a number"
                    )
                    msg = wx.MessageDialog(self, error, "Error", wx.OK | wx.ICON_ERROR)
                    msg.ShowModal()
                    msg.Destroy()
                    return

        # Ensure there are no negative values
        for line in self.gradient_box.GetValue().split("\n"):
            if line != "":
                if float(line) < 0:
                    # Give an error message saying that there must be no negative values
                    msg = wx.MessageDialog(
                        self,
                        "There must be no negative values",
                        "Error",
                        wx.OK | wx.ICON_ERROR,
                    )
                    msg.ShowModal()
                    msg.Destroy()
                    return
        # If Bruker, ensure that all gradient percentages are between 0 and 100
        if self.main_frame.spectrometer == "Bruker":
            for line in self.gradient_box.GetValue().split("\n"):
                if line != "":
                    if float(line) > 100:
                        # Give an error message saying that there must be no values greater than 100
                        error = (
                            "There must be no gradient percentages greater than 100: "
                            + line
                            + " is greater than 100"
                        )
                        msg = wx.MessageDialog(
                            self, error, "Error", wx.OK | wx.ICON_ERROR
                        )
                        msg.ShowModal()
                        msg.Destroy()
                        return
        # If Varian, ensure that all gradient DAC values are between 0 and 30000
        else:
            for line in self.gradient_box.GetValue().split("\n"):
                if line != "":
                    if float(line) > 30000:
                        # Give an error message saying that there must be no values greater than 30000
                        error = (
                            "There must be no DAC gradient values greater than 30000: "
                            + line
                            + " is greater than 30000"
                        )
                        msg = wx.MessageDialog(
                            self, error, "Error", wx.OK | wx.ICON_ERROR
                        )
                        msg.ShowModal()
                        msg.Destroy()
                        return
        # Ensure that the number of gradients entered is the same as the number of slices
        if len(self.gradient_box.GetValue().split("\n")) != len(self.main_frame.y_data):
            # Give an error message saying that the number of gradients must be the same as the number of slices
            error = (
                "The number of gradients must be the same as the number of slices in the data: "
                + str(len(self.gradient_box.GetValue().split("\n")))
                + " gradients entered, "
                + str(len(self.main_frame.y_data))
                + " slices in the data"
            )
            msg = wx.MessageDialog(self, error, "Error", wx.OK | wx.ICON_ERROR)
            msg.ShowModal()
            msg.Destroy()
            return
        file = open("gradients.txt", "w")
        file.write(self.gradient_box.GetValue())
        file.close()
        if self.main_frame.spectrometer == "Bruker":
            self.main_frame.gradients_percent = []
            self.main_frame.gradients = []
            for line in self.gradient_box.GetValue().split("\n"):
                if line != "":
                    self.main_frame.gradients_percent.append(float(line))

            self.main_frame.gradients_percent = np.array(
                self.main_frame.gradients_percent
            )
            self.main_frame.gradients = (
                (self.main_frame.gradients_percent / 100)
                * self.main_frame.max_gradient
                * self.main_frame.integral_factor
            )

        else:
            self.gradients_DAC = []
            for line in self.gradient_box.GetValue().split("\n"):
                if line != "":
                    self.gradients_DAC.append(float(line))
            self.gradients_DAC = np.array(self.gradients_DAC)

            self.gradients = (
                np.array(self.gradients_DAC) * self.main_frame.DAC_conversion
            )
            self.gradients_percent = (
                np.array(self.gradients_DAC)
                * self.main_frame.DAC_conversion
                / self.main_frame.max_gradient
                * 100
            )
            self.main_frame.gradients_percent = self.gradients_percent
            self.main_frame.gradients = self.gradients


class RelaxFit(wx.Frame):
    def __init__(self, title, parent=None):
        self.main_frame = parent
        # Get the monitor size and set the window size to 85% of the monitor size
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        self.display_index = wx.Display.GetFromWindow(parent)
        self.display_index_current = self.display_index
        self.width = int(1.0 * sizes[self.display_index][0])
        self.height = int(0.875 * sizes[self.display_index][1])
        wx.Frame.__init__(
            self, parent=parent, title=title, size=(self.width, self.height)
        )
        self.panel_relax = wx.Panel(self, -1)
        self.main_relax_sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_relax_sizer)

        self.fig_relax = Figure()
        self.fig_relax.tight_layout()
        self.canvas_relax = FigCanvas(self, -1, self.fig_relax)
        self.main_relax_sizer.Add(self.canvas_relax, 10, flag=wx.GROW)
        self.toolbar_relax = NavigationToolbar(self.canvas_relax)
        self.main_relax_sizer.Add(self.toolbar_relax, 0, wx.EXPAND)

        self.sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.titlecolor = "black"

        self.initial_values()
        self.make_relax_sizer()
        self.plot_relax_data()
        self.Show()

        # Bind method to check/resize the window when the frame is moved
        self.Bind(wx.EVT_MOVE, self.OnMoveFrame)

        # Bind method to resize the window when the frame is resized
        self.Bind(wx.EVT_SIZE, self.OnSizeFrame)

    def OnMoveFrame(self, event):
        # Get the new default display if the frame is moved
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        display_index = wx.Display.GetFromWindow(self)
        if display_index != self.display_index_current:
            self.display_index_current = display_index
            self.width = int(1.0 * sizes[display_index][0])
            self.height = int(0.875 * sizes[display_index][1])
            self.SetSize((self.width, self.height))
            self.canvas_relax.SetSize(
                (
                    self.width * 0.0104,
                    (self.height - self.relax_sizer.GetMinSize()[1] - 100) * 0.0104,
                )
            )
            self.fig_relax.set_size_inches(
                self.width * 0.0104,
                (self.height - self.relax_sizer.GetMinSize()[1] - 100) * 0.0104,
            )
            self.UpdateRelaxFrame()
        event.Skip()

    def OnSizeFrame(self, event):
        # Get the new frame size
        self.width, self.height = self.GetSize()
        self.SetSize((self.width, self.height))
        self.canvas_relax.SetSize(
            (
                self.width * 0.0104,
                (self.height - self.relax_sizer.GetMinSize()[1] - 100) * 0.0104,
            )
        )
        self.fig_relax.set_size_inches(
            self.width * 0.0104,
            (self.height - self.relax_sizer.GetMinSize()[1] - 100) * 0.0104,
        )
        self.UpdateRelaxFrame()
        event.Skip()

    def UpdateRelaxFrame(self):
        self.canvas_relax.draw()
        self.canvas_relax.Refresh()
        self.canvas_relax.Update()
        self.panel_relax.Refresh()
        self.panel_relax.Update()

    # The place where initial global variables are defined
    def initial_values(self):

        self.whole_plot = False  # Default to having only the diffusion data in a single plot with no diffusion coefficient subplots
        self.monoexponential_fit = False

        # Input ppms and y data from the main frame
        self.x_data = self.main_frame.new_x_ppms
        self.y_data = self.main_frame.nmrdata.data.T

        # Initially have noise region selection set to false
        self.noise_region_selection = False

        # Create an array to store the min and max ppm values for the selected regions
        self.selected_regions_of_interest = []

        self.AddROI = False

        self.ROI_color = []  # Empty array to store the colors of the ROIs
        self.deleted_ROI_number = (
            0  # Parameter to store the number of ROI's which have been deleted
        )

        self.deleted_slices = (
            []
        )  # Array to hold the indexes of the slices which have been deleted

    def make_relax_sizer(self):
        self.relax_sizer = wx.BoxSizer(wx.HORIZONTAL)
        # Check box for Varian/Bruker data (default to Bruker)

        # Create a button that opens a file for a user to input the delay times
        self.delays_label = wx.StaticBox(self, -1, "Delay Times")
        self.delays_sizer_total = wx.StaticBoxSizer(self.delays_label, wx.VERTICAL)
        self.delays_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.delay_times_button = wx.Button(self, -1, "Input Delay Times")
        self.delay_times_button.Bind(wx.EVT_BUTTON, self.OnInputDelayTimes)
        self.delays_sizer.AddSpacer(5)
        self.delays_sizer.Add(self.delay_times_button)
        self.delays_sizer.AddSpacer(5)
        self.relax_sizer.AddSpacer(5)
        self.delays_sizer_total.AddSpacer(4)
        self.delays_sizer_total.Add(self.delays_sizer)
        self.delays_sizer_total.AddSpacer(3)
        self.relax_sizer.Add(self.delays_sizer_total)

        # Then have button which will allow a user to drag over a section where they wish to estimate the noise level
        # This can then be plotted as a shaded region on the plot
        self.noise_label = wx.StaticBox(self, -1, "Noise Region")
        self.noise_sizer_total = wx.StaticBoxSizer(self.noise_label, wx.VERTICAL)
        self.noise_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.select_noise_button = wx.Button(self, -1, "Select Noise Region")
        self.select_noise_button.Bind(wx.EVT_BUTTON, self.OnSelectNoise)
        self.noise_sizer.AddSpacer(5)
        self.noise_sizer.Add(self.select_noise_button)
        self.noise_sizer.AddSpacer(5)

        # Then have a TextCtrl for the minimum SNR for the relaxation coefficient to be estimated (default to 10)
        self.noise_factor = 10
        self.noise_factor_label = wx.StaticText(self, -1, "Minimum SNR:")
        self.noise_factor_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.noise_sizer.AddSpacer(10)
        self.noise_sizer.Add(self.noise_factor_label)
        self.noise_sizer.AddSpacer(5)
        self.noise_factor_box = wx.TextCtrl(
            self, -1, str(self.noise_factor), size=(50, -1)
        )
        self.noise_sizer.Add(self.noise_factor_box)
        self.noise_sizer.AddSpacer(5)
        self.noise_sizer_total.AddSpacer(3)
        self.noise_sizer_total.Add(self.noise_sizer)
        self.noise_sizer_total.AddSpacer(2)

        # Have radio box where a user can choose to fit the data to obtain R1 or R2 relaxation rates
        self.fitting_type_label = wx.StaticBox(self, -1, "Relaxation Type:")
        self.fitting_type_sizer = wx.StaticBoxSizer(
            self.fitting_type_label, wx.HORIZONTAL
        )
        self.R1_fit = False
        self.R2_fit = True
        self.choices = ["R\u2081", "R\u2082"]
        self.R1R2_radiobox = wx.RadioBox(
            self, -1, choices=self.choices, style=wx.RA_HORIZONTAL
        )
        self.R1R2_radiobox.SetSelection(1)
        self.R1R2_radiobox.Bind(wx.EVT_RADIOBOX, self.OnFitSelection)

        self.fitting_type_sizer.AddSpacer(5)
        self.fitting_type_sizer.Add(self.R1R2_radiobox)
        self.fitting_type_sizer.AddSpacer(5)

        # Need to have a fitting sizer which will contain all the fitting buttons
        self.fitting_label = wx.StaticBox(self, -1, "Fitting")
        self.fitting_sizer_total = wx.StaticBoxSizer(self.fitting_label, wx.VERTICAL)
        self.fitting_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.fitting_sizer.AddSpacer(5)

        # Can then have a button which will fit the relaxation equation at all ppms across the whole spectrum that are higher than the noise level
        self.whole_spectrum_fitting_button = wx.Button(self, -1, "Fit Whole Spectrum")
        self.whole_spectrum_fitting_button.Bind(
            wx.EVT_BUTTON, self.OnWholeSpectrumFitting
        )
        self.fitting_sizer.Add(self.whole_spectrum_fitting_button)
        self.fitting_sizer.AddSpacer(5)
        # This can then be plotted
        # In addition, for each ppm can get a plot of I/I0 for all points which is also plotted next to this

        # Then can have a button saying select region of interest. The relaxation coefficient in this region can be estimated along with an error (from the standard deviation of the points)
        # Can plot this distribution of relaxation coefficients and it should resemble a Gaussian distribution

        # Have a button to add a new region of interest
        self.add_region_button = wx.Button(self, -1, "Add ROI")
        self.add_region_button.Bind(wx.EVT_BUTTON, self.OnAddROI)
        self.fitting_sizer.Add(self.add_region_button)
        self.fitting_sizer.AddSpacer(5)

        # Have a button to delete a region of interest
        self.delete_region_button = wx.Button(self, -1, "Delete ROI")
        self.delete_region_button.Bind(wx.EVT_BUTTON, self.OnDeleteROI)
        self.fitting_sizer.Add(self.delete_region_button)
        self.fitting_sizer.AddSpacer(5)

        # Can then have a button which will fit the Relaxation equation to the mean values of the points above the noise in the region of interest
        self.region_fitting_button = wx.Button(self, -1, "Fit")
        self.region_fitting_button.Bind(wx.EVT_BUTTON, self.OnRegionFitting)
        self.fitting_sizer.Add(self.region_fitting_button)
        self.fitting_sizer.AddSpacer(5)

        # Have a button which will perform a biexponential fit on the data in the region of interest
        self.biexponential_fitting_button = wx.Button(self, -1, "Biexponential Fit")
        self.biexponential_fitting_button.Bind(
            wx.EVT_BUTTON, self.OnBiexponentialFitting
        )
        self.fitting_sizer.Add(self.biexponential_fitting_button)
        self.fitting_sizer.AddSpacer(5)

        self.fitting_sizer_total.AddSpacer(4)
        self.fitting_sizer_total.Add(self.fitting_sizer)
        self.fitting_sizer_total.AddSpacer(4)

        # Have a button for printing fitted values
        self.print_fitted_values_button = wx.Button(self, -1, "Print Fit")
        self.print_fitted_values_button.Bind(wx.EVT_BUTTON, self.OnPrintFittedValues)
        self.fitting_sizer.Add(self.print_fitted_values_button)
        self.fitting_sizer.AddSpacer(5)

        # Have a box containing other functions such as a button to delete a slice from the plot and repeat the fitting
        self.other_functions_label = wx.StaticBox(self, -1, "Other Functions")
        self.other_functions_sizer = wx.StaticBoxSizer(
            self.other_functions_label, wx.VERTICAL
        )
        self.other_functions_sizer.AddSpacer(4)
        self.delete_slice_button = wx.Button(self, -1, "Delete Slice")
        self.delete_slice_button.Bind(wx.EVT_BUTTON, self.OnDeleteSlice)
        self.other_functions_sizer.Add(self.delete_slice_button)
        self.other_functions_sizer.AddSpacer(4)

        # Creating a sizer for changing the y axis limits in the spectrum
        self.intensity_label = wx.StaticBox(self, -1, "Y Axis Zoom (%):")
        self.intensity_sizer = wx.StaticBoxSizer(self.intensity_label, wx.VERTICAL)
        width = 100
        self.intensity_slider = FloatSlider(
            self, id=-1, value=0, minval=-1, maxval=10, res=0.01, size=(width, height)
        )
        self.intensity_slider.Bind(wx.EVT_SLIDER, self.OnIntensityScrollRelax)
        self.intensity_sizer.AddSpacer(5)
        self.intensity_sizer.Add(self.intensity_slider)

        self.relax_sizer.AddSpacer(5)

        self.relax_sizer.Add(self.noise_sizer_total)
        self.relax_sizer.AddSpacer(5)
        self.relax_sizer.Add(self.fitting_type_sizer)
        self.relax_sizer.AddSpacer(5)
        self.relax_sizer.Add(self.fitting_sizer_total)
        self.relax_sizer.AddSpacer(5)
        self.relax_sizer.Add(self.other_functions_sizer)
        self.relax_sizer.AddSpacer(5)
        self.relax_sizer.Add(self.intensity_sizer)

        self.main_relax_sizer.AddSpacer(5)
        self.main_relax_sizer.Add(self.relax_sizer)
        self.main_relax_sizer.AddSpacer(5)

        # See if delays.txt file containing delays exists
        try:
            file = open("delays.txt", "r")
            self.delays = np.loadtxt("delays.txt")
            file.close()
        except:
            pass

    def OnFitSelection(self, event):
        if self.R1R2_radiobox.GetSelection() == 0:
            self.R1_fit = True
            self.R2_fit = False
        else:
            self.R1_fit = False
            self.R2_fit = True

    def OnPrintFittedValues(self, event):
        if self.whole_plot == False:
            # Give a message saying please perform whole spectrum and ROI fitting first
            msg = wx.MessageDialog(
                self,
                "Please perform whole spectrum and ROI fitting first",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        if self.whole_plot == True:
            if self.monoexponential_fit == False:
                # Give a message saying please perform ROI fitting first
                msg = wx.MessageDialog(
                    self,
                    "Please perform ROI fitting first",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                return
            else:
                # For each exponential decay, print the raw data and errors
                for i, ROI in enumerate(self.selected_regions_of_interest):
                    self.average_y_data_in_ROI_above_noise = (
                        self.average_y_data_in_ROI_above_noise_total[i]
                    )
                    self.error_y_data_in_ROI_above_noise = (
                        self.error_y_data_in_ROI_above_noise_total[i]
                    )
                    self.error_I_I0_in_ROI = self.error_I_I0_in_ROI_total[i]
                    self.I0_average_in_ROI = self.I0_average_in_ROI_total[i]

                    print(
                        self.average_y_data_in_ROI_above_noise / self.I0_average_in_ROI
                    )
                    print(self.error_I_I0_in_ROI)

    def plot_relax_data(self):

        self.ax_relax = self.fig_relax.add_subplot(111)
        count = 1
        self.slice_plots = []
        for i, data in enumerate(self.y_data):
            (line,) = self.ax_relax.plot(
                self.x_data, data, linewidth=0.5, label=str(count)
            )
            self.slice_plots.append(line)
            count += 1
        self.ax_relax.set_xlim([self.x_data[0], self.x_data[-1]])

        self.ax_relax.set_xlabel(self.main_frame.nmrdata.axislabels[1])
        legend = self.ax_relax.legend(title="Slice Number")
        legend.get_title().set_color(self.titlecolor)
        self.ax_relax.set_ylabel("Intensity")

        self.noise_region = self.ax_relax.axvspan(
            min(self.x_data), min(self.x_data), alpha=0.2, color="gray"
        )

    def OnInputDelayTimes(self, event):
        delays = DelaysManualInput(title="Delays Manual Input", parent=self)

    def OnIntensityScrollRelax(self, event):
        # Function to change the y axis limits
        intensity_percent = 10 ** float(self.intensity_slider.GetValue())

        self.ax_relax.set_ylim(
            -(np.max(self.y_data) / 8) / (intensity_percent / 100),
            np.max(self.y_data) / (intensity_percent / 100),
        )
        self.UpdateRelaxFrame()

    def OnSelectNoise(self, event):

        # self.noise_region.set_xy([[min(self.x_data),0],[min(self.x_data),1],[min(self.x_data),1],[min(self.x_data),0]])
        # self.noise_region.set_xy([min(self.x_data),0])
        # if(self.whole_plot==True):
        #     self.noise_region_2.set_xy([[min(self.x_data),0],[min(self.x_data),1],[min(self.x_data),1],[min(self.x_data),0]])
        #     self.noise_region_3.set_xy([[min(self.x_data),0],[min(self.x_data),1],[min(self.x_data),1],[min(self.x_data),0]])

        self.UpdateRelaxFrame()

        self.press = False
        self.move = False
        if self.whole_plot == True:
            self.canvas_relax = FigCanvas(self, -1, self.fig_relax)
        self.noise_select_press = self.canvas_relax.mpl_connect(
            "button_press_event", self.OnPress
        )
        self.noise_select_release = self.canvas_relax.mpl_connect(
            "button_release_event", self.OnReleaseNoise
        )
        self.noise_select_move = self.canvas_relax.mpl_connect(
            "motion_notify_event", self.OnMove
        )

    def OnPress(self, event):
        if self.whole_plot == False:
            if event.inaxes == self.ax_relax:
                self.press = True
                self.x0 = event.xdata
        else:
            if (
                event.inaxes == self.ax_relax
                or event.inaxes == self.ax_relax_whole_fit
                or event.inaxes == self.ax_relax_I0_whole_fit
                or event.inaxes == self.relax_coefficient_plot
            ):
                self.press = True
                self.x0 = event.xdata

    def OnMove(self, event):
        if self.whole_plot == False:
            if event.inaxes == self.ax_relax:
                self.move_noise(event)

        else:
            if (
                event.inaxes == self.ax_relax
                or event.inaxes == self.ax_relax_whole_fit
                or event.inaxes == self.ax_relax_I0_whole_fit
                or event.inaxes == self.relax_coefficient_plot
            ):
                self.move_noise(event)

    def move_noise(self, event):
        if self.press:
            self.move = True
            self.x1 = event.xdata
            if self.x1 > self.x0:
                xmax = self.x1
                xmin = self.x0
            else:
                xmax = self.x0
                xmin = self.x1
            self.noise_region.set_x(xmin)
            self.noise_region.set_width(xmax - xmin)
            # self.noise_region.set(xy=[[xmin,0],[xmin,1],[xmax,1],[xmax,0]])
            if self.whole_plot == True:
                # self.noise_region_2.set(xy=[[xmin,0],[xmin,1],[xmax,1],[xmax,0]])
                # self.noise_region_3.set(xy=[[xmin,0],[xmin,1],[xmax,1],[xmax,0]])
                self.noise_region_2.set_x(xmin)
                self.noise_region_2.set_width(xmax - xmin)
                self.noise_region_3.set_x(xmin)
                self.noise_region_3.set_width(xmax - xmin)

            self.UpdateRelaxFrame()

    def release_noise(self, event):
        if self.press:
            self.x2 = event.xdata
            if self.x2 > self.x0:
                xmax = self.x2
                xmin = self.x0
            else:
                xmax = self.x0
                xmin = self.x2
            self.noise_x_initial = xmin
            self.noise_x_final = xmax

            self.noise_region.set_x(xmin)
            self.noise_region.set_width(xmax - xmin)
            # self.noise_region.set(xy=[[xmin,0],[xmin,1],[xmax,1],[xmax,0]])
            if self.whole_plot == True:
                self.noise_region_2.set_x(xmin)
                self.noise_region_2.set_width(xmax - xmin)
                self.noise_region_3.set_x(xmin)
                self.noise_region_3.set_width(xmax - xmin)
                # self.noise_region_2.set(xy=[[xmin,0],[xmin,1],[xmax,1],[xmax,0]])
                # self.noise_region_3.set(xy=[[xmin,0],[xmin,1],[xmax,1],[xmax,0]])

            self.UpdateRelaxFrame()
        self.press = False
        self.move = False
        self.canvas_relax.mpl_disconnect(self.noise_select_press)
        self.canvas_relax.mpl_disconnect(self.noise_select_move)
        self.canvas_relax.mpl_disconnect(self.noise_select_release)

        # Turn on noise region selection flag
        self.noise_region_selection = True

    def OnReleaseNoise(self, event):
        if self.whole_plot == False:
            if event.inaxes == self.ax_relax:
                self.release_noise(event)

        else:
            if (
                event.inaxes == self.ax_relax
                or event.inaxes == self.ax_relax_whole_fit
                or event.inaxes == self.ax_relax_I0_whole_fit
                or event.inaxes == self.relax_coefficient_plot
            ):
                self.release_noise(event)

    def OnWholeSpectrumFitting(self, event):
        # Initially check that the noise region has been selected
        try:
            self.noise_x_initial
            self.noise_x_final
        except:
            # Give an error message saying noise region has not been selected
            msg = wx.MessageDialog(
                self,
                "Please select a noise region before fitting",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return

        # Check that the minimum SNR has been entered and is a value greater than 0
        try:
            self.noise_factor = float(self.noise_factor_box.GetValue())
        except:
            # Give an error message saying noise factor has not been entered
            msg = wx.MessageDialog(
                self,
                "Please enter a minimum SNR before fitting",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return

        if self.noise_factor <= 0:
            # Give an error message saying noise factor has not been entered correctly
            msg = wx.MessageDialog(
                self,
                "Please enter a minimum SNR greater than 0 before fitting",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return

        # Check that the delays have been found
        try:
            self.delays
        except:
            # Give an error message saying gradients have not been found
            msg = wx.MessageDialog(
                self,
                "Please input the delays before fitting",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return

        # Find out the standard deviation of the noise region
        noise_region = self.y_data[
            :,
            np.where(
                (self.x_data >= self.noise_x_initial)
                & (self.x_data <= self.noise_x_final)
            )[0],
        ]
        noise_region_std = np.std(noise_region)

        # Find out the ppms which have intensity of the minimum intensity slice greater than the noise level in all slices
        self.ppms_above_noise_slices = []
        for i, data in enumerate(self.y_data):
            self.ppms_above_noise_slice = []
            for j, intensity in enumerate(data):
                if intensity > self.noise_factor * noise_region_std:
                    self.ppms_above_noise_slice.append(self.x_data[j])

            self.ppms_above_noise_slices.append(self.ppms_above_noise_slice)

        # Find out the ppms which are in all slices and remove ones which are not
        self.ppms_above_noise = self.ppms_above_noise_slices[0]
        for i, ppm in enumerate(self.ppms_above_noise):
            for j, ppms_slice in enumerate(self.ppms_above_noise_slices):
                if ppm not in ppms_slice:
                    del self.ppms_above_noise[i]
                    break

        # Get the indices of all the ppms which are above the noise level
        self.ppms_above_noise_indices = []
        for i, ppm in enumerate(self.ppms_above_noise):
            self.ppms_above_noise_indices.append(np.where(self.x_data == ppm)[0][0])

        # Remove all the y data points which are below the noise level
        self.SelectDataAboveThreshold()

        # Separate data into point by point
        self.SeparateDataIntoPointByPoint()

        # For all the ppms that have intensity above the noise threshold, fit the Relaxation equation to the data
        self.fitted_I0_global = []
        self.fitted_relax_global = []

        for i, ppm in enumerate(self.ppms_above_noise):
            self.y_vals = np.real(self.y_data_point_by_point[i])

            # Start at a few different initial relaxation values coefficients so that don't get stuck in local minima
            fits = []
            chi_squareds = []
            for j, relaxation_initial in enumerate(np.linspace(1, 100, 10)):
                fit = self.leastsq_global([np.max(self.y_vals), relaxation_initial])
                fits.append(fit)
                chi_squareds.append(np.sum(self.chi_global(fit) ** 2))

            fit = fits[np.argmin(chi_squareds)]
            self.fitted_I0_global.append(fit[0])
            self.fitted_relax_global.append(fit[1])

        self.PlotWholeSpectrumFitting()

    def PlotWholeSpectrumFitting(self):
        self.fig_relax.clear()
        self.fig_relax.tight_layout()

        gs = gridspec.GridSpec(2, 2)

        self.ax_relax = self.fig_relax.add_subplot(gs[0, :])
        self.ax_relax_whole_fit = self.fig_relax.add_subplot(
            gs[1, 0], sharex=self.ax_relax, sharey=self.ax_relax
        )
        self.ax_relax_I0_whole_fit = self.fig_relax.add_subplot(
            gs[1, 1], sharex=self.ax_relax
        )

        count = 1
        self.slice_plots = []
        for i, data in enumerate(self.y_data):
            (line,) = self.ax_relax.plot(
                self.x_data, data, linewidth=0.5, label=str(count)
            )
            self.slice_plots.append(line)
            count += 1
        self.ax_relax.set_xlim([self.x_data[0], self.x_data[-1]])

        self.ax_relax.set_xlabel(self.main_frame.nmrdata.axislabels[1])
        legend = self.ax_relax.legend(title="Slice Number")
        legend.get_title().set_color(self.titlecolor)
        self.ax_relax.set_ylabel("Intensity")
        if self.R1_fit == True:
            self.ax_relax.set_title("R1 Data", color=self.titlecolor)
        else:
            self.ax_relax.set_title("R2 Data", color=self.titlecolor)

        self.noise_region = self.ax_relax.axvspan(
            self.noise_x_initial, self.noise_x_final, alpha=0.2, color="gray"
        )

        # Plot the fitted relaxation coefficients and use a twiny to also plot the initial slice of the spectrum
        self.ax_relax_whole_fit.plot(self.x_data, self.y_data[0])
        self.ax_relax_whole_fit.set_xlabel(self.main_frame.nmrdata.axislabels[1])
        self.ax_relax_whole_fit.set_yticks([])
        self.relax_coefficient_plot = self.ax_relax_whole_fit.twinx()
        if self.R1_fit == True:
            self.relax_coefficient_plot.set_ylabel(r"R$_1$ Coefficient (s$^{-1}$)")
        else:
            self.relax_coefficient_plot.set_ylabel(r"R$_2$ Coefficient (s$^{-1}$)")
        self.relax_coefficient_plot.set_xlabel(self.main_frame.nmrdata.axislabels[1])
        self.relax_coefficient_plot.set_xlim([self.x_data[0], self.x_data[-1]])
        self.relax_coefficient_plot.scatter(
            self.ppms_above_noise, self.fitted_relax_global, color="tab:red", s=0.5
        )
        self.relax_coefficient_plot.yaxis.tick_left()
        self.relax_coefficient_plot.yaxis.set_label_position("left")
        if self.R1_fit == True:
            self.relax_coefficient_plot.set_title("R1 vs PPM", color=self.titlecolor)
        else:
            self.relax_coefficient_plot.set_title("R2 vs PPM", color=self.titlecolor)
        self.noise_region_2 = self.ax_relax_whole_fit.axvspan(
            self.noise_x_initial, self.noise_x_final, alpha=0.2, color="gray"
        )

        # Plot I/I0 for every chosen ppm across the spectrum for all slices
        for i, selected_y_data in enumerate(self.y_data_above_noise):
            self.ax_relax_I0_whole_fit.scatter(
                self.ppms_above_noise,
                np.array(selected_y_data) / self.fitted_I0_global,
                s=0.5,
            )

        self.ax_relax_I0_whole_fit.set_xlabel(self.main_frame.nmrdata.axislabels[1])
        self.ax_relax_I0_whole_fit.set_ylabel(r"I/I$_0$")
        self.ax_relax_I0_whole_fit.set_xlim([self.x_data[0], self.x_data[-1]])
        self.ax_relax_I0_whole_fit.set_title(r"I/I$_0$ vs PPM", color=self.titlecolor)
        self.noise_region_3 = self.ax_relax_I0_whole_fit.axvspan(
            self.noise_x_initial, self.noise_x_final, alpha=0.2, color="gray"
        )

        if self.whole_plot == True:
            if len(self.ROI_regions) == 0:
                self.ROI_regions = []
                self.ROI_regions_2 = []
                self.ROI_regions_3 = []
            else:
                ROI_regions = []
                ROI_regions_2 = []
                ROI_regions_3 = []
                for i, ROI_region in enumerate(self.ROI_regions):
                    bottom_left = ROI_region.get_xy()
                    width = ROI_region.get_width()
                    ROI_regions.append(
                        self.ax_relax.axvspan(
                            bottom_left[0],
                            bottom_left[0] + width,
                            alpha=0.2,
                            color=self.ROI_color[i],
                        )
                    )
                    ROI_regions_2.append(
                        self.ax_relax_whole_fit.axvspan(
                            bottom_left[0],
                            bottom_left[0] + width,
                            alpha=0.2,
                            color=self.ROI_color[i],
                        )
                    )
                    ROI_regions_3.append(
                        self.ax_relax_I0_whole_fit.axvspan(
                            bottom_left[0],
                            bottom_left[0] + width,
                            alpha=0.2,
                            color=self.ROI_color[i],
                        )
                    )
                self.ROI_regions = ROI_regions
                self.ROI_regions_2 = ROI_regions_2
                self.ROI_regions_3 = ROI_regions_3

        else:
            self.ROI_regions = []
            self.ROI_regions_2 = []
            self.ROI_regions_3 = []

        # Turn on whole plot mode once whole plot mode has been completed
        self.whole_plot = True

        self.UpdateRelaxFrame()

    def SelectDataAboveThreshold(self):
        # Remove all the y data points which are below the noise level
        self.y_data_above_noise = []
        for i, data in enumerate(self.y_data):
            self.y_data_above_noise_slice = []
            for index in self.ppms_above_noise_indices:
                self.y_data_above_noise_slice.append(data[index])

            self.y_data_above_noise.append(self.y_data_above_noise_slice)

    def SeparateDataIntoPointByPoint(self):
        # Separate the data into arrays of intensities for each ppm that has intensity above noise threshold in all slices
        self.y_data_point_by_point = []
        for i, ppm in enumerate(self.ppms_above_noise):
            y_data = []
            for j, data in enumerate(self.y_data_above_noise):
                y_data.append(data[i])

            y_data = np.array(y_data)
            self.y_data_point_by_point.append(y_data)

        self.y_data_point_by_point = np.array(self.y_data_point_by_point)

    def T2_RelaxationEquation(self, p0):
        I0, R = p0
        return I0 * np.exp(-self.delays * R)

    def T1_RelaxationEquation(self, p0):
        I0, R = p0
        # Check to see if the first slice is positively or negatively phased
        if np.max(self.x_data[0]) != np.abs(self.x_data[0]):
            return I0 * (1 - 2 * np.exp(-self.delays * R))
        else:
            return I0 * (2 * np.exp(-self.delays * R) - 1)

    def chi_global(self, p0):
        if self.R1_fit == True:
            return self.y_vals - self.T1_RelaxationEquation(p0)
        else:
            return self.y_vals - self.T2_RelaxationEquation(p0)

    def leastsq_global(self, p0):
        fit = leastsq(self.chi_global, p0)
        return fit[0]

    def OnAddROI(self, event):
        # Check that the full spectrum has been fitted first
        if self.whole_plot != True:
            # Give an error message saying full spectrum has not been fitted
            msg = wx.MessageDialog(
                self,
                "Please fit the whole spectrum before selecting a region of interest",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return

        # Check if have pressed AddROI but have not selected a region of interest then delete the previous value
        if self.AddROI == True:
            if self.ROI_regions[-1].get_xy()[0][0] == self.x_data[0]:
                del self.ROI_regions[-1]
                del self.ROI_regions_2[-1]
                del self.ROI_regions_3[-1]

        self.AddROI == True

        # Add new region plots with the default values (min ppm values)
        self.ROI_color.append(
            self.main_frame.colours[
                len(self.selected_regions_of_interest) + self.deleted_ROI_number
            ]
        )
        self.ROI_regions.append(
            self.ax_relax.axvspan(
                self.x_data[0], self.x_data[0], alpha=0.2, color=self.ROI_color[-1]
            )
        )
        self.ROI_regions_2.append(
            self.ax_relax_whole_fit.axvspan(
                self.x_data[0], self.x_data[0], alpha=0.2, color=self.ROI_color[-1]
            )
        )
        self.ROI_regions_3.append(
            self.ax_relax_I0_whole_fit.axvspan(
                self.x_data[0], self.x_data[0], alpha=0.2, color=self.ROI_color[-1]
            )
        )

        self.UpdateRelaxFrame()

        self.canvas_relax.mpl_disconnect(self.noise_select_press)
        self.canvas_relax.mpl_disconnect(self.noise_select_move)
        self.canvas_relax.mpl_disconnect(self.noise_select_release)

        self.press = False
        self.move = False

        self.select_ROI_press = self.canvas_relax.mpl_connect(
            "button_press_event", self.OnPressROI
        )
        self.select_ROI_release = self.canvas_relax.mpl_connect(
            "button_release_event", self.OnReleaseROI
        )
        self.select_ROI_move = self.canvas_relax.mpl_connect(
            "motion_notify_event", self.OnMoveROI
        )

    def OnDeleteROI(self, event):
        # Check that a region of interest has been added first
        if len(self.selected_regions_of_interest) == 0:
            # Give an error message saying no regions of interest have been added
            msg = wx.MessageDialog(
                self,
                "No regions of interest have been added",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return

        # When mouse is over a region of interest highlight that region (make alpha=0.75)
        # When mouse is moved away from region of interest make alpha=0.2 again
        # When mouse is clicked delete that region of interest

        self.canvas_relax.mpl_disconnect(self.noise_select_press)
        self.canvas_relax.mpl_disconnect(self.noise_select_move)
        self.canvas_relax.mpl_disconnect(self.noise_select_release)

        self.canvas_relax.mpl_disconnect(self.select_ROI_press)
        self.canvas_relax.mpl_disconnect(self.select_ROI_move)
        self.canvas_relax.mpl_disconnect(self.select_ROI_release)

        self.delete_ROI_press = self.canvas_relax.mpl_connect(
            "button_press_event", self.OnPressDeleteROI
        )
        self.delete_ROI_highlight = self.canvas_relax.mpl_connect(
            "motion_notify_event", self.OnHighlightROI
        )

    def OnPressDeleteROI(self, event):
        if (
            event.inaxes == self.ax_relax
            or event.inaxes == self.ax_relax_whole_fit
            or event.inaxes == self.ax_relax_I0_whole_fit
            or event.inaxes == self.relax_coefficient_plot
        ):
            # Find out the highlighted slices
            if len(self.selected_regions_of_interest) == 1:
                for i in self.highlighted_regions:
                    self.selected_regions_of_interest = []
                    self.ROI_regions[i].set_alpha(0.2)
                    self.ROI_regions_2[i].set_alpha(0.2)
                    self.ROI_regions_3[i].set_alpha(0.2)
                    self.ROI_regions[i].set_x(self.x_data[0])
                    self.ROI_regions[i].set_width(0)
                    self.ROI_regions_2[i].set_x(self.x_data[0])
                    self.ROI_regions_2[i].set_width(0)
                    self.ROI_regions_3[i].set_x(self.x_data[0])
                    self.ROI_regions_3[i].set_width(0)
                    # self.ROI_regions[i].set_xy([[self.x_data[0],0],[self.x_data[0],1],[self.x_data[0],1],[self.x_data[0],0]])
                    # self.ROI_regions_2[i].set_xy([[self.x_data[0],0],[self.x_data[0],1],[self.x_data[0],1],[self.x_data[0],0]])
                    # self.ROI_regions_3[i].set_xy([[self.x_data[0],0],[self.x_data[0],1],[self.x_data[0],1],[self.x_data[0],0]])
                    del self.ROI_regions[i]
                    del self.ROI_regions_2[i]
                    del self.ROI_regions_3[i]
                    del self.ROI_color[i]
                    self.deleted_ROI_number += 1
            else:
                for i in self.highlighted_regions:
                    del self.selected_regions_of_interest[i]
                    self.ROI_regions[i].set_alpha(0.2)
                    self.ROI_regions_2[i].set_alpha(0.2)
                    self.ROI_regions_3[i].set_alpha(0.2)
                    self.ROI_regions[i].set_x(self.x_data[0])
                    self.ROI_regions[i].set_width(0)
                    self.ROI_regions_2[i].set_x(self.x_data[0])
                    self.ROI_regions_2[i].set_width(0)
                    self.ROI_regions_3[i].set_x(self.x_data[0])
                    self.ROI_regions_3[i].set_width(0)
                    # self.ROI_regions[i].set_xy([[self.x_data[0],0],[self.x_data[0],1],[self.x_data[0],1],[self.x_data[0],0]])
                    # self.ROI_regions_2[i].set_xy([[self.x_data[0],0],[self.x_data[0],1],[self.x_data[0],1],[self.x_data[0],0]])
                    # self.ROI_regions_3[i].set_xy([[self.x_data[0],0],[self.x_data[0],1],[self.x_data[0],1],[self.x_data[0],0]])
                    del self.ROI_regions[i]
                    del self.ROI_regions_2[i]
                    del self.ROI_regions_3[i]
                    del self.ROI_color[i]
                    self.deleted_ROI_number += 1

            self.UpdateRelaxFrame()

            # Disconnect highlight and press events
            self.canvas_relax.mpl_disconnect(self.delete_ROI_press)
            self.canvas_relax.mpl_disconnect(self.delete_ROI_highlight)

            if self.monoexponential_fit == True:
                self.OnRegionFitting(event)

    def OnHighlightROI(self, event):
        if (
            event.inaxes == self.ax_relax
            or event.inaxes == self.ax_relax_whole_fit
            or event.inaxes == self.ax_relax_I0_whole_fit
            or event.inaxes == self.relax_coefficient_plot
        ):
            x0 = event.xdata
            self.highlight_ROI(x0)

    def highlight_ROI(self, x0):
        self.highlighted_regions = []
        # Check if x0 is within any of the regions of interest
        if len(self.selected_regions_of_interest) == 1:
            region = self.selected_regions_of_interest[0]
            if x0 >= region[0] and x0 <= region[1]:
                self.ROI_regions[0].set_alpha(0.75)
                self.ROI_regions_2[0].set_alpha(0.75)
                self.ROI_regions_3[0].set_alpha(0.75)
                self.highlighted_regions.append(0)

        else:
            for i, region in enumerate(self.selected_regions_of_interest):
                if x0 >= region[0] and x0 <= region[1]:
                    self.ROI_regions[i].set_alpha(0.75)
                    self.ROI_regions_2[i].set_alpha(0.75)
                    self.ROI_regions_3[i].set_alpha(0.75)
                    self.highlighted_regions.append(i)
            for i, region in enumerate(self.selected_regions_of_interest):
                if i not in self.highlighted_regions:
                    self.ROI_regions[i].set_alpha(0.2)
                    self.ROI_regions_2[i].set_alpha(0.2)
                    self.ROI_regions_3[i].set_alpha(0.2)

        self.UpdateRelaxFrame()

    def OnPressROI(self, event):
        if (
            event.inaxes == self.ax_relax
            or event.inaxes == self.ax_relax_whole_fit
            or event.inaxes == self.ax_relax_I0_whole_fit
            or event.inaxes == self.relax_coefficient_plot
        ):
            self.press = True
            self.x0 = event.xdata

    def OnMoveROI(self, event):

        if (
            event.inaxes == self.ax_relax
            or event.inaxes == self.ax_relax_whole_fit
            or event.inaxes == self.ax_relax_I0_whole_fit
            or event.inaxes == self.relax_coefficient_plot
        ):
            self.move_ROI(event)

    def move_ROI(self, event):
        if self.press:
            self.move = True
            self.x1 = event.xdata
            if self.x1 > self.x0:
                xmax = self.x1
                xmin = self.x0
            else:
                xmax = self.x0
                xmin = self.x1

            self.ROI_regions[-1].set_x(xmin)
            self.ROI_regions[-1].set_width(xmax - xmin)
            self.ROI_regions_2[-1].set_x(xmin)
            self.ROI_regions_2[-1].set_width(xmax - xmin)
            self.ROI_regions_3[-1].set_x(xmin)
            self.ROI_regions_3[-1].set_width(xmax - xmin)
            # self.ROI_regions[-1].set_xy([[xmin,0],[xmin,1],[xmax,1],[xmax,0]])
            # self.ROI_regions_2[-1].set_xy([[xmin,0],[xmin,1],[xmax,1],[xmax,0]])
            # self.ROI_regions_3[-1].set_xy([[xmin,0],[xmin,1],[xmax,1],[xmax,0]])
            self.UpdateRelaxFrame()

    def release_ROI(self, event):
        if self.press:
            self.x2 = event.xdata
            if self.x2 > self.x0:
                xmax = self.x2
                xmin = self.x0
            else:
                xmax = self.x0
                xmin = self.x2
            self.ROI_x_initial = xmin
            self.ROI_x_final = xmax
            self.ROI_regions[-1].set_x(xmin)
            self.ROI_regions[-1].set_width(xmax - xmin)
            self.ROI_regions_2[-1].set_x(xmin)
            self.ROI_regions_2[-1].set_width(xmax - xmin)
            self.ROI_regions_3[-1].set_x(xmin)
            self.ROI_regions_3[-1].set_width(xmax - xmin)
            # self.ROI_regions[-1].set_xy([[xmin,0],[xmin,1],[xmax,1],[xmax,0]])
            # self.ROI_regions_2[-1].set_xy([[xmin,0],[xmin,1],[xmax,1],[xmax,0]])
            # self.ROI_regions_3[-1].set_xy([[xmin,0],[xmin,1],[xmax,1],[xmax,0]])

            # Add the min and max ppm values to the array of selected regions of interest
            self.selected_regions_of_interest.append([xmin, xmax])

            self.UpdateRelaxFrame()

            self.press = False
            self.move = False
            self.canvas_relax.mpl_disconnect(self.select_ROI_press)
            self.canvas_relax.mpl_disconnect(self.select_ROI_move)
            self.canvas_relax.mpl_disconnect(self.select_ROI_release)

    def OnReleaseROI(self, event):
        if (
            event.inaxes == self.ax_relax
            or event.inaxes == self.ax_relax_whole_fit
            or event.inaxes == self.ax_relax_I0_whole_fit
            or event.inaxes == self.relax_coefficient_plot
        ):
            self.release_ROI(event)

    def OnRegionFitting(self, event):
        # Remove the deleted slices from the gradient and gradient percentages lists

        self.ppms_in_ROI_total = []
        self.ppms_in_ROI_indices_total = []
        self.average_y_data_in_ROI_above_noise_total = []
        self.error_y_data_in_ROI_above_noise_total = []
        self.error_I_I0_in_ROI_total = []
        self.error_log_I_I0_in_ROI_total = []
        self.error_I_I0_in_ROI_total = []
        self.I0_average_in_ROI_total = []
        self.fitted_relax_ROI_total = []
        self.fitted_I0_ROI_total = []
        self.mean_fitted_relax_ROI_total = []
        self.mean_fitted_I0_ROI_total = []
        self.fitted_I0_total = []
        self.fitted_relax_total = []

        for i, region in enumerate(self.selected_regions_of_interest):
            self.ROI_x_initial = region[0]
            self.ROI_x_final = region[1]
            # Get the indices of the ppms which are in the ROI and have intensity above the noise in all slices
            self.ppms_in_ROI = []
            self.ppms_in_ROI_indices = []
            for i, ppm in enumerate(self.ppms_above_noise):
                if ppm >= self.ROI_x_initial and ppm <= self.ROI_x_final:
                    self.ppms_in_ROI.append(ppm)
                    self.ppms_in_ROI_indices.append(i)

            self.average_y_data_in_ROI_above_noise = []
            self.error_y_data_in_ROI_above_noise = []
            self.error_I_I0_in_ROI = []

            for i, data in enumerate(self.y_data_above_noise):
                self.y_data_in_ROI_above_noise_slice = []
                self.I_I0_in_ROI_slice = []
                for index in self.ppms_in_ROI_indices:
                    self.y_data_in_ROI_above_noise_slice.append(np.real(data[index]))
                    self.I_I0_in_ROI_slice.append(
                        np.real(data[index] / self.fitted_I0_global[index])
                    )

                self.average_y_data_in_ROI_above_noise.append(
                    np.mean(np.array(self.y_data_in_ROI_above_noise_slice))
                )
                self.error_y_data_in_ROI_above_noise.append(
                    np.std(np.array(self.y_data_in_ROI_above_noise_slice))
                )
                self.error_I_I0_in_ROI.append(np.std(np.array(self.I_I0_in_ROI_slice)))

            self.average_y_data_in_ROI_above_noise = np.array(
                self.average_y_data_in_ROI_above_noise
            )
            self.error_y_data_in_ROI_above_noise = np.array(
                self.error_y_data_in_ROI_above_noise
            )
            self.error_I_I0_in_ROI = np.array(self.error_I_I0_in_ROI)

            # Also need the error in I/I0 for each slice in the ROI
            self.error_log_I_I0_in_ROI = []
            self.error_I_I0_in_ROI = []
            self.I0_average_in_ROI = []
            for i, slice_I0 in enumerate(self.y_data_above_noise):
                self.I_I0_in_ROI_slice = []
                self.I0_average_in_ROI_slice = []
                for index in self.ppms_in_ROI_indices:
                    self.I_I0_in_ROI_slice.append(
                        np.real(slice_I0[index] / self.fitted_I0_global[index])
                    )
                    self.I0_average_in_ROI_slice.append(
                        np.real(self.fitted_I0_global[index])
                    )

                self.error_log_I_I0_in_ROI.append(
                    np.std(np.log(np.array(self.I_I0_in_ROI_slice)))
                )
                self.error_I_I0_in_ROI.append(np.std(np.array(self.I_I0_in_ROI_slice)))
                self.I0_average_in_ROI.append(
                    np.mean(np.array(self.I0_average_in_ROI_slice))
                )

            self.error_log_I_I0_in_ROI = np.array(self.error_log_I_I0_in_ROI)
            self.I0_average_in_ROI = np.array(self.I0_average_in_ROI)
            self.error_I_I0_in_ROI = np.array(self.error_I_I0_in_ROI)

            self.fitted_relax_ROI = []
            self.fitted_I0_ROI = []

            for index in self.ppms_in_ROI_indices:
                self.fitted_relax_ROI.append(np.real(self.fitted_relax_global[index]))
                self.fitted_I0_ROI.append(np.real(self.fitted_I0_global[index]))

            self.mean_fitted_relax_ROI = np.mean(np.array(self.fitted_relax_ROI))
            self.mean_fitted_I0_ROI = np.mean(np.array(self.fitted_I0_ROI))

            self.ppms_in_ROI_total.append(self.ppms_in_ROI)
            self.ppms_in_ROI_indices_total.append(self.ppms_in_ROI_indices)
            self.average_y_data_in_ROI_above_noise_total.append(
                self.average_y_data_in_ROI_above_noise
            )
            self.error_y_data_in_ROI_above_noise_total.append(
                self.error_y_data_in_ROI_above_noise
            )
            self.error_I_I0_in_ROI_total.append(self.error_I_I0_in_ROI)
            self.error_log_I_I0_in_ROI_total.append(self.error_log_I_I0_in_ROI)
            self.error_I_I0_in_ROI_total.append(self.error_I_I0_in_ROI)
            self.I0_average_in_ROI_total.append(self.I0_average_in_ROI)
            self.fitted_relax_ROI_total.append(self.fitted_relax_ROI)
            self.fitted_I0_ROI_total.append(self.fitted_I0_ROI)
            self.mean_fitted_relax_ROI_total.append(self.mean_fitted_relax_ROI)
            self.mean_fitted_I0_ROI_total.append(self.mean_fitted_I0_ROI)

            # Fit the relaxation equation to the data for all points in the ROI, use the standard deviation of all I/I0 values as the error
            self.fitted_I0, self.fitted_D = self.leastsq_ROI(
                [np.max(self.average_y_data_in_ROI_above_noise), 1e-9]
            )

            self.fitted_I0_total.append(self.fitted_I0)
            self.fitted_relax_total.append(self.fitted_D)

        self.monoexponential_fit = True

        self.PlotRegionFitting()

    def OnBiexponentialFitting(self, event):
        if self.monoexponential_fit != True:
            # Give an error message to say please perform monoexponential fitting first
            msg = wx.MessageDialog(
                self,
                "Please perform monoexponential fitting first",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        elif len(self.selected_regions_of_interest) > 1:
            # Give an error message saying that biexponential fitting is only supported while one region of interest is present
            msg = wx.MessageDialog(
                self,
                "Biexponential fitting is only supported while one region of interest is present",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        else:
            # Perform a biexponential fit
            # Loop over various initial guesses for the biexponential fit
            chi_squared = 100000
            relax_1_array = np.linspace(1, 100, 5)
            relax_2_array = np.linspace(1, 100, 5)
            f1_array = np.linspace(0.1, 0.9, 5)
            for i, relax_1 in enumerate(relax_1_array):
                for j, relax_2 in enumerate(relax_2_array):
                    for k, f1 in enumerate(f1_array):
                        p0 = [
                            np.max(self.average_y_data_in_ROI_above_noise),
                            relax_1,
                            relax_2,
                            f1,
                        ]
                        fit = leastsq(self.chi_biexponential_ROI, p0)
                        if (
                            np.sum(self.chi_biexponential_ROI(fit[0]) ** 2)
                            < chi_squared
                        ):
                            chi_squared = np.sum(
                                self.chi_biexponential_ROI(fit[0]) ** 2
                            )
                            best_fit = fit
            fit = best_fit
            I0_ROI = np.abs(fit[0][0])
            relax_1_ROI = np.abs(fit[0][1])
            relax_2_ROI = np.abs(fit[0][2])
            f1_ROI = np.abs(fit[0][3])
            xvals = np.linspace(min(self.delays), max(self.delays), 100)

            # Plot the biexponential fit
            delays = self.delays
            self.delays = xvals
            if self.R1_fit == True:
                self.ax_relax_fit.plot(
                    xvals * 1000,
                    self.T1_Biexponential([I0_ROI, relax_1_ROI, relax_2_ROI, f1_ROI])
                    / I0_ROI,
                    color="tab:red",
                    linestyle="--",
                    label=r"R$_{1}$ = "
                    + str(np.round(1 / relax_1_ROI, 2))
                    + r" s$^{-1}$, R$_{2}$ = "
                    + str(np.round(1 / relax_2_ROI, 2))
                    + r" s$^{-1}$, f = "
                    + str(np.round(f1_ROI, 2)),
                )
            else:
                self.ax_relax_fit.plot(
                    xvals * 1000,
                    self.T2_Biexponential([I0_ROI, relax_1_ROI, relax_2_ROI, f1_ROI])
                    / I0_ROI,
                    color="tab:red",
                    linestyle="--",
                    label=r"R$_{2}$ = "
                    + str(np.round(relax_1_ROI, 2))
                    + r" s$^{-1}$, R$_{2}$ = "
                    + str(np.round(relax_2_ROI, 2))
                    + r" s$^{-1}$, f = "
                    + str(np.round(f1_ROI, 2)),
                )
            legend = self.ax_relax_fit.legend(fontsize=8)
            legend.get_title().set_color(self.titlecolor)
            self.delays = delays
            self.UpdateRelaxFrame()

    def PlotRegionFitting(self):
        # Generate 3 extra plots for the region fitting (I/I0 vs gradient^2 with fitted curve, log(I/I0) vs gradient^2 with fitted curve, histogram of T2 coefficients within ROI)
        self.fig_relax.clear()
        self.fig_relax.tight_layout()

        gs = gridspec.GridSpec(2, 3)

        self.ax_relax = self.fig_relax.add_subplot(gs[0, 0:2])
        self.ax_relax_whole_fit = self.fig_relax.add_subplot(
            gs[1, 0], sharex=self.ax_relax
        )
        self.ax_relax_I0_whole_fit = self.fig_relax.add_subplot(
            gs[1, 1], sharex=self.ax_relax
        )
        self.ax_relax_fit = self.fig_relax.add_subplot(gs[0, 2])
        self.ax_relax_histogram = self.fig_relax.add_subplot(gs[1, 2])

        matplotlib.rcParams.update({"font.size": 8})

        count = 1
        self.slice_plots = []
        for i, data in enumerate(self.y_data):
            (line,) = self.ax_relax.plot(
                self.x_data, data, linewidth=0.5, label=str(count)
            )
            self.slice_plots.append(line)
            count += 1
        self.ax_relax.set_xlim([self.x_data[0], self.x_data[-1]])

        self.ax_relax.set_xlabel(self.main_frame.nmrdata.axislabels[1])
        legend = self.ax_relax.legend(title="Slice Number", fontsize=8)
        legend.get_title().set_color(self.titlecolor)
        self.ax_relax.set_ylabel("Intensity", fontsize=8)
        if self.R1_fit == True:
            self.ax_relax.set_title("R1 Data", color=self.titlecolor, fontsize=10)
        else:
            self.ax_relax.set_title("T2 Data", color=self.titlecolor, fontsize=10)
        self.ax_relax.tick_params(axis="both", which="major", labelsize=8)

        self.noise_region = self.ax_relax.axvspan(
            self.noise_x_initial, self.noise_x_final, alpha=0.2, color="gray"
        )

        # Plot the fitted relaxation coefficients and use a twiny to also plot the initial slice of the spectrum
        self.ax_relax_whole_fit.plot(self.x_data, self.y_data[0])
        self.ax_relax_whole_fit.set_xlabel(
            self.main_frame.nmrdata.axislabels[1], fontsize=8
        )
        self.ax_relax_whole_fit.tick_params(axis="both", which="major", labelsize=8)
        self.ax_relax_whole_fit.set_yticks([])
        self.relax_coefficient_plot = self.ax_relax_whole_fit.twinx()
        if self.R1_fit == True:
            self.relax_coefficient_plot.set_ylabel(r"R1 Coefficient ($s^{-1}$)")
        else:
            self.relax_coefficient_plot.set_ylabel(r"R2 Coefficient ($s^{-1}$)")
        self.relax_coefficient_plot.set_xlabel(
            self.main_frame.nmrdata.axislabels[1], fontsize=8
        )
        self.relax_coefficient_plot.set_xlim([self.x_data[0], self.x_data[-1]])
        self.relax_coefficient_plot.scatter(
            self.ppms_above_noise, self.fitted_relax_global, color="tab:red", s=0.5
        )
        self.relax_coefficient_plot.yaxis.tick_left()
        self.relax_coefficient_plot.yaxis.set_label_position("left")

        if self.R1_fit == True:
            self.relax_coefficient_plot.set_title(
                "R1 vs PPM", color=self.titlecolor, fontsize=10
            )
        else:
            self.relax_coefficient_plot.set_title(
                "R2 vs PPM", color=self.titlecolor, fontsize=10
            )
        self.relax_coefficient_plot.tick_params(axis="both", which="major", labelsize=8)
        self.noise_region_2 = self.ax_relax_whole_fit.axvspan(
            self.noise_x_initial, self.noise_x_final, alpha=0.2, color="gray"
        )

        # Plot I/I0 for every chosen ppm across the spectrum for all slices
        for i, selected_y_data in enumerate(self.y_data_above_noise):
            self.ax_relax_I0_whole_fit.scatter(
                self.ppms_above_noise,
                np.array(selected_y_data) / self.fitted_I0_global,
                s=0.5,
            )

        self.ax_relax_I0_whole_fit.set_xlabel(
            self.main_frame.nmrdata.axislabels[1], fontsize=8
        )
        self.ax_relax_I0_whole_fit.set_ylabel(r"I/I$_0$", fontsize=8)
        self.ax_relax_I0_whole_fit.set_xlim([self.x_data[0], self.x_data[-1]])
        self.ax_relax_I0_whole_fit.set_title(
            r"I/I$_0$ vs PPM", color=self.titlecolor, fontsize=10
        )
        self.ax_relax_I0_whole_fit.tick_params(axis="both", which="major", labelsize=8)
        self.noise_region_3 = self.ax_relax_I0_whole_fit.axvspan(
            self.noise_x_initial, self.noise_x_final, alpha=0.2, color="gray"
        )

        # Plot the ROI regions on all plots
        self.ROI_regions = []
        self.ROI_regions_2 = []
        self.ROI_regions_3 = []
        for i, region in enumerate(self.selected_regions_of_interest):
            self.ROI_x_initial = region[0]
            self.ROI_x_final = region[1]
            color = self.ROI_color[i]
            self.ROI_regions.append(
                self.ax_relax.axvspan(
                    self.ROI_x_initial, self.ROI_x_final, alpha=0.2, color=color
                )
            )
            self.ROI_regions_2.append(
                self.ax_relax_whole_fit.axvspan(
                    self.ROI_x_initial, self.ROI_x_final, alpha=0.2, color=color
                )
            )
            self.ROI_regions_3.append(
                self.ax_relax_I0_whole_fit.axvspan(
                    self.ROI_x_initial, self.ROI_x_final, alpha=0.2, color=color
                )
            )

        self.fitted_gaussian_parameters = []
        sigmas = []
        # Plot a histogram of the relaxation coefficients in the ROI
        for i, fitted_relax_ROI in enumerate(self.fitted_relax_ROI_total):
            if len(fitted_relax_ROI) > 0:
                self.ax_relax_histogram.hist(
                    fitted_relax_ROI,
                    bins=int(len(fitted_relax_ROI)),
                    color=self.main_frame.colours[i],
                    edgecolor=self.ROI_color[i],
                    alpha=0.25,
                )
                if self.R1_fit == True:
                    self.ax_relax_histogram.set_xlabel(
                        r"R1 Value (s$^{-1}$)", fontsize=8
                    )
                    self.ax_relax_histogram.set_title(
                        "Histogram of R1 Values", color=self.titlecolor, fontsize=10
                    )
                else:
                    self.ax_relax_histogram.set_xlabel(
                        r"R2 Value (s$^{-1}$)", fontsize=8
                    )
                    self.ax_relax_histogram.set_title(
                        "Histogram of R2 Values", color=self.titlecolor, fontsize=10
                    )
                self.ax_relax_histogram.set_ylabel("Frequency Density", fontsize=8)

                # Get the bin size of the histogram
                self.bin_size = self.ax_relax_histogram.patches[0].get_width()
                self.bin_centers = np.arange(
                    min(fitted_relax_ROI) + self.bin_size / 2,
                    max(fitted_relax_ROI),
                    self.bin_size,
                )
                self.bin_centers = np.array(self.bin_centers)
                self.bin_centers = self.bin_centers[
                    np.where(self.bin_centers <= max(fitted_relax_ROI))
                ]
                self.bin_centers = self.bin_centers[
                    np.where(self.bin_centers >= min(fitted_relax_ROI))
                ]
                self.bin_centers = np.array(self.bin_centers)

                # Get the frequency densities of the histogram in each bin
                self.frequency_density = []
                for j, bin_center in enumerate(self.bin_centers):
                    self.frequency_density.append(
                        len(
                            np.where(
                                (fitted_relax_ROI >= bin_center - self.bin_size / 2)
                                & (fitted_relax_ROI < bin_center + self.bin_size / 2)
                            )[0]
                        )
                    )

                # Fit a gaussian to the histogram of relaxation coefficients, this will be the error in the relaxation coefficient
                self.fitted_relax_ROI = np.array(fitted_relax_ROI)
                result = self.leastsq_gaussian_ROI(
                    [1, np.mean(fitted_relax_ROI), np.std(fitted_relax_ROI)]
                )
                if result[0] != "Failed":
                    A, mu, sigma = self.leastsq_gaussian_ROI(
                        [1, np.mean(fitted_relax_ROI), np.std(fitted_relax_ROI)]
                    )
                    if sigma > max(fitted_relax_ROI) - min(fitted_relax_ROI):
                        sigma = np.std(fitted_relax_ROI)
                        self.fitted_gaussian_parameters.append([A, mu, sigma])
                        self.ax_relax_histogram.plot(
                            self.bin_centers,
                            self.gaussian_ROI(self.bin_centers, A, mu, sigma),
                            label=r"$\sigma$ = "
                            + "{:.3e}".format(sigma)
                            + r"s$^{-1}$ (std)",
                            color=self.ROI_color[i],
                        )
                    else:
                        self.fitted_gaussian_parameters.append([A, mu, sigma])
                        self.ax_relax_histogram.plot(
                            self.bin_centers,
                            self.gaussian_ROI(self.bin_centers, A, mu, sigma),
                            label=r"$\sigma$ = "
                            + "{:.3e}".format(sigma)
                            + r"s$^{-1}$ (gauss fit)",
                            color=self.ROI_color[i],
                        )
                    sigmas.append(np.abs(sigma))
                    legend = self.ax_relax_histogram.legend(fontsize=8)
                    legend.get_title().set_color(self.titlecolor)
                else:
                    # Gaussian fit failed - setting sigma to the standard deviation of the diffusion coefficients and A to max frequency density
                    sigma = np.std(fitted_relax_ROI)
                    A = max(self.frequency_density)
                    self.fitted_gaussian_parameters.append(
                        [A, np.mean(fitted_relax_ROI), sigma]
                    )
                    # # Give an error message saying that one of the ROI windows is too small. Please increase the size of the ROI window and try again
                    msg = wx.MessageDialog(
                        self,
                        "Gaussian error fit did not work, error set to standard deviation. If desired, increase the size of the ROI window and try again",
                        "Error",
                        wx.OK | wx.ICON_ERROR,
                    )
                    msg.ShowModal()
                    msg.Destroy()
                    sigmas.append(np.abs(sigma))
                    legend = self.ax_relax_histogram.legend(fontsize=8)
                    legend.get_title().set_color(self.titlecolor)

            else:
                # Give an error message saying that one of the ROI windows is too small. Please delete the ROI and try again
                msg = wx.MessageDialog(
                    self,
                    "One of the ROI windows is too small. Please delete the ROI and try again",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                return

        for i, region in enumerate(self.selected_regions_of_interest):
            self.ROI_x_initial = region[0]
            self.ROI_x_final = region[1]
            self.mean_fitted_relax_ROI = self.mean_fitted_relax_ROI_total[i]
            self.average_y_data_in_ROI_above_noise = (
                self.average_y_data_in_ROI_above_noise_total[i]
            )
            self.error_y_data_in_ROI_above_noise = (
                self.error_y_data_in_ROI_above_noise_total[i]
            )
            self.error_I_I0_in_ROI = self.error_I_I0_in_ROI_total[i]
            self.I0_average_in_ROI = self.I0_average_in_ROI_total[i]
            self.fitted_relax_ROI = self.fitted_relax_ROI_total[i]
            self.fitted_I0_ROI = self.fitted_I0_ROI_total[i]
            self.mean_fitted_I0_ROI = self.mean_fitted_I0_ROI_total[i]
            self.fitted_I0 = self.fitted_I0_total[i]
            self.fitted_D = self.fitted_relax_total[i]
            self.error_log_I_I0_in_ROI = self.error_log_I_I0_in_ROI_total[i]
            self.error_I_I0_in_ROI = self.error_I_I0_in_ROI_total[i]

            # Plot the fitted curve for the ROI data
            self.ax_relax_fit.errorbar(
                np.array(self.delays) * 1000,
                self.average_y_data_in_ROI_above_noise / self.I0_average_in_ROI,
                yerr=self.error_I_I0_in_ROI,
                fmt="o",
                markersize=1,
                capsize=2,
                color=self.ROI_color[i],
            )

            delays = self.delays
            xvals = np.linspace(np.min(self.delays), np.max(self.delays), 100)
            self.delays = xvals
            if self.R1_fit == True:
                self.ax_relax_fit.plot(
                    xvals * 1000,
                    self.T1_RelaxationEquation(
                        [self.mean_fitted_I0_ROI, self.mean_fitted_relax_ROI]
                    )
                    / self.mean_fitted_I0_ROI,
                    label=r"R$_{1}$ = "
                    + str(round(self.mean_fitted_relax_ROI, 2))
                    + "+/-"
                    + str(round(sigmas[i], 2)),
                    color=self.ROI_color[i],
                )
            else:
                self.ax_relax_fit.plot(
                    xvals * 1000,
                    self.T2_RelaxationEquation(
                        [self.mean_fitted_I0_ROI, self.mean_fitted_relax_ROI]
                    )
                    / self.mean_fitted_I0_ROI,
                    label=r"R$_{2}$ = "
                    + str(round(self.mean_fitted_relax_ROI, 2))
                    + "+/-"
                    + str(round(sigmas[i], 2)),
                    color=self.ROI_color[i],
                )
            self.delays = delays

        self.ax_relax_fit.set_xlabel(r"Delays (ms)", fontsize=8)
        self.ax_relax_fit.set_ylabel(r"I/I$_0$", fontsize=8)
        if self.R1_fit == True:
            self.ax_relax_fit.set_title(
                "Fitted R1 Relaxation", color=self.titlecolor, fontsize=10
            )
        else:
            self.ax_relax_fit.set_title(
                "Fitted R2 Relaxation", color=self.titlecolor, fontsize=10
            )
        legend = self.ax_relax_fit.legend(fontsize=8)
        legend.get_title().set_color(self.titlecolor)

        self.fig_relax.tight_layout()

        self.UpdateRelaxFrame()

    def T2_Biexponential(self, p0):
        I0, R2_1, R2_2, f1 = p0
        # Ensure all values are positive
        R2_1 = np.abs(R2_1)
        R2_2 = np.abs(R2_2)
        f1 = np.abs(f1)
        I0 = np.abs(I0)
        return I0 * (
            f1 * np.exp(-self.delays * R2_1) + (1 - f1) * np.exp(-self.delays * R2_2)
        )

    def T1_Biexponential(self, p0):
        I0, R1_1, R1_2, f1 = p0
        # Ensure all values are positive
        R1_1 = np.abs(R1_1)
        R1_2 = np.abs(R1_2)
        f1 = np.abs(f1)
        I0 = np.abs(I0)
        if np.max(self.x_data[0]) != np.max(np.abs(self.x_data[0])):
            return I0 * (
                f1 * (1 - 2 * np.exp(-self.delays * R1_1))
                + (1 - f1) * (1 - 2 * np.exp(-self.delays * R1_2))
            )
        else:
            return I0 * (
                f1 * (2 * np.exp(-self.delays * R1_1) - 1)
                + (1 - f1) * (2 * np.exp(-self.delays * R1_2) - 1)
            )

    def chi_biexponential_ROI(self, p0):
        if self.R1_fit == True:
            return (
                self.average_y_data_in_ROI_above_noise - self.T1_Biexponential(p0)
            ) / self.error_y_data_in_ROI_above_noise
        else:
            return (
                self.average_y_data_in_ROI_above_noise - self.T2_Biexponential(p0)
            ) / self.error_y_data_in_ROI_above_noise

    def leastsq_biexponential(self, p0):
        fit = leastsq(self.chi_biexponential_ROI, p0)
        return fit[0]

    def leastsq_ROI(self, p0):
        fit = leastsq(self.chi_ROI, p0)
        return fit[0]

    def chi_ROI(self, p0):
        return (
            self.average_y_data_in_ROI_above_noise - self.T2_RelaxationEquation(p0)
        ) / self.error_y_data_in_ROI_above_noise

    def gaussian_ROI(self, x, A, mu, sigma):
        return A * np.exp(-((x - mu) ** 2) / (2 * sigma**2))

    def chi_gaussian_ROI(self, p0):
        return self.frequency_density - self.gaussian_ROI(
            self.bin_centers, p0[0], p0[1], p0[2]
        )

    def leastsq_gaussian_ROI(self, p0):
        try:
            fit = leastsq(self.chi_gaussian_ROI, p0)
            return fit[0]
        except:
            return ["Failed"]

    def OnDeleteSlice(self, event):

        # Check to see if the gradients have already been inputted
        try:
            self.delays
        except:
            # Give an error message saying that the gradients must be inputted first
            msg = wx.MessageDialog(
                self,
                "Please input the delays first before deleting a slice",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        # Check to see that the number of slices is four or more
        if len(self.y_data) < 4:
            # Give an error message saying that there must be at least four slices
            msg = wx.MessageDialog(
                self,
                "There must be at least three slices in the data to perform relaxation data fitting so cannot delete a slice",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        self.delete_slice = True
        # Bring up a dialog box to ask which slice to delete
        self.delete_slice_index = 0
        self.delete_slice_dialog = DeleteSliceDialog("Delete Slice", self)

    # Function to continue the deletion of a slice after the dialog box has been closed and completed by the user
    def continue_deletion(self):
        # Check to see if the full spectrum fitting has been performed
        if self.whole_plot != True:
            # Delete the correct slice in the y data
            self.y_data = np.delete(self.y_data, self.delete_slice_index, axis=0)

            # Delete the correct value in gradients
            self.delays = np.delete(self.delays, self.delete_slice_index)

            # Redo the plotting
            self.fig_relax.clear()
            self.fig_relax.tight_layout()
            self.plot_relax_data()

            # If the noise region has already been selected, then redo the plotting of this
            if self.noise_region_selection == True:
                self.noise_region = self.ax_relax.axvspan(
                    self.noise_x_initial, self.noise_x_final, alpha=0.2, color="gray"
                )
            self.UpdateRelaxFrame()
        elif self.whole_plot == True and self.monoexponential_fit != True:
            # Delete the correct slice in the y data
            self.y_data = np.delete(self.y_data, self.delete_slice_index, axis=0)

            # Delete the correct value in gradients
            self.delays = np.delete(self.delays, self.delete_slice_index)

            self.OnWholeSpectrumFitting(event=None)

        elif self.whole_plot == True and self.monoexponential_fit == True:
            # Delete the correct slice in the y data
            self.y_data = np.delete(self.y_data, self.delete_slice_index, axis=0)

            # Delete the correct value in gradients
            self.delays = np.delete(self.delays, self.delete_slice_index)
            self.OnWholeSpectrumFitting(event=None)
            self.OnRegionFitting(event=None)


# Need this box to delete any slice from the data, as well as remove it from the list of gradients
class DeleteSliceDialog(wx.Frame):
    def __init__(self, title, parent):
        self.main_frame = parent
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        width = int(0.2 * self.monitorWidth)
        height = int(0.1 * self.monitorHeight)
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.panel_delete_slice = wx.Panel(self, -1)
        self.main_delete_slice = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_delete_slice)

        self.make_delete_slice_sizer()
        self.Show()

    def make_delete_slice_sizer(self):
        # Make a sizer to hold the text box and button
        self.delete_slice_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.delete_slice_sizer.AddSpacer(5)
        # ComboBox for slice number
        self.slice_number_label = wx.StaticText(self, -1, "Slice Number:")
        self.delete_slice_sizer.Add(self.slice_number_label)
        self.delete_slice_sizer.AddSpacer(5)
        self.slice_number_choices = []
        for i in range(len(self.main_frame.y_data)):
            self.slice_number_choices.append(str(i + 1))
        self.slice_number_combobox = wx.ComboBox(
            self, -1, choices=self.slice_number_choices, style=wx.CB_READONLY
        )
        self.delete_slice_sizer.Add(self.slice_number_combobox)
        self.delete_slice_sizer.AddSpacer(5)
        # Have a button to confirm the deletion
        self.confirm_button = wx.Button(self, -1, "Delete")
        self.confirm_button.Bind(wx.EVT_BUTTON, self.OnConfirmDelete)
        self.delete_slice_sizer.Add(self.confirm_button)
        self.delete_slice_sizer.AddSpacer(5)

        self.main_delete_slice.AddSpacer(5)
        self.main_delete_slice.Add(self.delete_slice_sizer)
        self.main_delete_slice.AddSpacer(5)

    def OnConfirmDelete(self, event):
        self.main_frame.delete_slice_index = (
            int(self.slice_number_combobox.GetValue()) - 1
        )
        self.main_frame.deleted_slices.append(self.main_frame.delete_slice_index)
        self.main_frame.continue_deletion()
        self.Destroy()


class DelaysManualInput(wx.Frame):
    def __init__(self, title, parent):
        self.main_frame = parent
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        width = int(0.3 * self.monitorWidth)
        height = int(0.5 * self.monitorHeight)
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.panel_delays_input = wx.Panel(self, -1)
        self.main_delays_input = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_delays_input)

        self.make_manual_delays_input_sizer()
        self.Show()

    def make_manual_delays_input_sizer(self):
        self.input_delays_label = wx.StaticBox(
            self, -1, "Input delays in seconds (one delay per line)"
        )
        self.input_delays_sizer = wx.StaticBoxSizer(
            self.input_delays_label, wx.VERTICAL
        )

        try:
            file = open("delays.txt")
            delays = []
            for line in file.readlines():
                line = line.split("\n")[0]
                delays.append(line)

            label = ""
            for delay in delays:
                label = label + delay + "\n"
            file.close()
        except:
            label = ""

        self.delay_box = wx.TextCtrl(
            self, -1, value=label, size=(250, 400), style=wx.TE_MULTILINE
        )
        self.input_delays_sizer.AddSpacer(3)
        self.input_delays_sizer.Add(self.delay_box)
        self.input_delays_sizer.AddSpacer(10)

        self.save_delays_button = wx.Button(self, -1, "Save delays")

        self.save_delays_button.Bind(wx.EVT_BUTTON, self.OnSaveDelays)

        self.input_delays_sizer.Add(self.save_delays_button)

        self.main_delays_input.Add(self.input_delays_sizer)

    def OnSaveDelays(self, event):
        # Ensure that there are not any blank lines in the delays
        new_lines = ""
        for i, line in enumerate(self.delay_box.GetValue().split("\n")):
            if line != "":
                new_lines = new_lines + line.rstrip() + "\n"

        if new_lines[-1] == "\n":
            new_lines = new_lines[:-1]
        self.delay_box.SetValue(new_lines)

        # Check that all the delays are numbers

        for delay in self.delay_box.GetValue().split("\n"):
            try:
                float(delay)
            except:
                error = (
                    "Please ensure that all delays entrered are numbers: "
                    + delay
                    + " is not a number"
                )
                msg = wx.MessageDialog(self, error, "Error", wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                return

        # Ensure that all delays are positive
        for delay in self.delay_box.GetValue().split("\n"):
            if float(delay) < 0:
                error = (
                    "Please ensure that all delays entrered are positive: "
                    + delay
                    + " is negative"
                )
                msg = wx.MessageDialog(self, error, "Error", wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                return

        # Ensure that the number of delays is the same as the number of slices
        if len(self.delay_box.GetValue().split("\n")) != len(self.main_frame.y_data):
            error = (
                "Please ensure that the number of delays is the same as the number of slices. There are "
                + str(len(self.main_frame.y_data))
                + " slices, but "
                + str(len(self.delay_box.GetValue().split("\n")))
                + " delays were given."
            )
            msg = wx.MessageDialog(self, error, "Error", wx.OK | wx.ICON_ERROR)
            msg.ShowModal()
            msg.Destroy()
            return

        file = open("delays.txt", "w")
        file.write(self.delay_box.GetValue())
        file.close()
        self.main_frame.delays = []
        for line in self.delay_box.GetValue().split("\n"):
            if line != "":
                self.main_frame.delays.append(float(line))

        self.main_frame.delays = np.array(self.main_frame.delays)


class ReadSession:
    def __init__(self, parent, session_file):
        self.main_frame = parent
        self.session_file = session_file
        # try:
        self.read_session_file()
        # except:
        #     # Give an error message saying that the session file could not be read (exiting)
        #     msg = wx.MessageDialog(self.main_frame, 'The session file could not be read. Please check the file paths in the .session file to ensure they are correct.', 'Error', wx.OK | wx.ICON_ERROR)
        #     msg.ShowModal()
        #     msg.Destroy()
        #     exit()

    def read_session_file(self):
        file = open(self.session_file, "r")
        lines = file.readlines()
        file.close()
        # Check to see what the first line is
        if len(lines[0].split("\n")[0].split()) == 1:
            # This is a session containing just one window
            # Check to see if the window is 1D, 2D or 3D
            plot_type = lines[0].split("\n")[0].split()[0]
            if plot_type == "1D" or plot_type == "1D stack":
                # This is a 1D window
                if plot_type == "1D":
                    stack = False
                else:
                    stack = True
                # Check to see if multiplot mode is on
                if lines[1].split("\n")[0].split(":")[1].split()[0] == "True":
                    multiplot = True
                    # Get the file path of the original data
                    file_path_original = lines[2].split("\n")[0].split("file_path:")[1]
                    self.main_frame.nmrdata = GetData(file_path_original)
                    self.main_frame.viewer = OneDViewer(
                        parent=self.main_frame, nmrdata=self.main_frame.nmrdata
                    )
                    self.main_frame.viewer.stack = stack
                    self.main_frame.main_sizer.Add(self.main_frame.viewer, 1, wx.EXPAND)
                    title = lines[3].split("\n")[0].split(":")[1]
                    p0_coarse = float(lines[4].split("\n")[0].split(":")[1])
                    p0_fine = float(lines[5].split("\n")[0].split(":")[1])
                    p1_coarse = float(lines[6].split("\n")[0].split(":")[1])
                    p1_fine = float(lines[7].split("\n")[0].split(":")[1])
                    colour = int(lines[8].split("\n")[0].split(":")[1])
                    linewidth = float(lines[9].split("\n")[0].split(":")[1])
                    reference_range = int(lines[10].split("\n")[0].split(":")[1])
                    reference_value = float(lines[11].split("\n")[0].split(":")[1])
                    vertical_range = int(lines[12].split("\n")[0].split(":")[1])
                    vertical_value = float(lines[13].split("\n")[0].split(":")[1])
                    multiply_range = int(lines[14].split("\n")[0].split(":")[1])
                    multiply_value = float(lines[15].split("\n")[0].split(":")[1])
                    pivot_point = float(lines[16].split("\n")[0].split(":")[1])
                    pivot_x = float(lines[17].split("\n")[0].split(":")[1])
                    pivot_visible = lines[18].split("\n")[0].split(":")[1]
                    self.choices = [title]
                    self.main_frame.viewer.line1.set_label(title)
                    self.main_frame.viewer.P0_slider.SetValue(p0_coarse)
                    self.main_frame.viewer.P1_slider.SetValue(p1_coarse)
                    self.main_frame.viewer.P0_slider_fine.SetValue(p0_fine)
                    self.main_frame.viewer.P1_slider_fine.SetValue(p1_fine)
                    self.main_frame.viewer.index = colour
                    self.main_frame.viewer.colour_chooser.SetSelection(colour)
                    self.main_frame.viewer.OnColourChoice1D(event=None)
                    self.main_frame.viewer.linewidth_slider.SetValue(linewidth)
                    self.main_frame.viewer.OnLinewidthScroll1D(event=None)
                    self.main_frame.viewer.reference_range_chooser.SetSelection(
                        reference_range
                    )
                    self.main_frame.viewer.OnReferenceCombo(event=None)
                    self.main_frame.viewer.reference_slider.SetValue(reference_value)
                    self.main_frame.viewer.OnReferenceScroll1D(event=None)
                    self.main_frame.viewer.vertical_range_chooser.SetSelection(
                        vertical_range
                    )
                    self.main_frame.viewer.OnVerticalCombo(event=None)
                    self.main_frame.viewer.vertical_slider.SetValue(vertical_value)
                    self.main_frame.viewer.OnVerticalScroll1D(event=None)
                    self.main_frame.viewer.multiply_range_chooser.SetSelection(
                        multiply_range
                    )
                    self.main_frame.viewer.OnMultiplyCombo(event=None)
                    self.main_frame.viewer.multiply_slider.SetValue(multiply_value)
                    self.main_frame.viewer.OnMultiplyScroll1D(event=None)
                    self.main_frame.viewer.pivot_line.set_xdata([pivot_point])
                    self.main_frame.viewer.pivot_x = pivot_x
                    if pivot_visible == "True":
                        self.main_frame.viewer.pivot_line.set_visible(True)
                    else:
                        self.main_frame.viewer.pivot_line.set_visible(False)
                    self.main_frame.viewer.OnSliderScroll1D(event=None)
                    # Read all the values into the values dictionary for the first plot
                    self.main_frame.viewer.values_dictionary[0] = {}
                    self.main_frame.viewer.values_dictionary[0][
                        "path"
                    ] = file_path_original
                    self.main_frame.viewer.values_dictionary[0]["title"] = title
                    self.main_frame.viewer.values_dictionary[0]["linewidth"] = linewidth
                    self.main_frame.viewer.values_dictionary[0]["color index"] = colour
                    self.main_frame.viewer.values_dictionary[0][
                        "original_ppms"
                    ] = self.main_frame.viewer.ppm_original
                    self.main_frame.viewer.values_dictionary[0][
                        "original_data"
                    ] = self.main_frame.viewer.nmrdata.data
                    self.main_frame.viewer.values_dictionary[0][
                        "dictionary"
                    ] = self.main_frame.viewer.nmrdata.dic
                    self.main_frame.viewer.values_dictionary[0][
                        "move left/right"
                    ] = reference_value
                    self.main_frame.viewer.values_dictionary[0][
                        "move left/right range index"
                    ] = reference_range
                    self.main_frame.viewer.values_dictionary[0][
                        "move up/down"
                    ] = vertical_value
                    self.main_frame.viewer.values_dictionary[0][
                        "move up/down range index"
                    ] = vertical_range
                    self.main_frame.viewer.values_dictionary[0][
                        "multiply value"
                    ] = multiply_value
                    self.main_frame.viewer.values_dictionary[0][
                        "multiply range index"
                    ] = multiply_range
                    self.main_frame.viewer.values_dictionary[0]["p0 Coarse"] = p0_coarse
                    self.main_frame.viewer.values_dictionary[0]["p1 Coarse"] = p1_coarse
                    self.main_frame.viewer.values_dictionary[0]["p0 Fine"] = p0_fine
                    self.main_frame.viewer.values_dictionary[0]["p1 Fine"] = p1_fine
                    self.main_frame.viewer.multiplot_mode = True
                    self.main_frame.viewer.files.first_drop = False
                    count = 1
                    # Loop over the rest of the lines to get the file paths of the other data
                    for i, line in enumerate(lines):
                        if i < 19:
                            continue
                        if line.split("\n")[0].split(":")[0] == "file_path":
                            file_path = line.split("\n")[0].split("file_path:")[1]
                            self.main_frame.viewer.values_dictionary[count] = {}
                            self.main_frame.viewer.values_dictionary[count][
                                "path"
                            ] = file_path
                        if line.split("\n")[0].split(":")[0] == "title":
                            title = line.split("\n")[0].split(":")[1]
                            self.main_frame.viewer.values_dictionary[count][
                                "title"
                            ] = title
                        if line.split("\n")[0].split(":")[0] == "linewidth":
                            linewidth = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "linewidth"
                            ] = linewidth
                        if line.split("\n")[0].split(":")[0] == "p0_coarse":
                            p0_coarse = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "p0 Coarse"
                            ] = p0_coarse
                        if line.split("\n")[0].split(":")[0] == "p1_coarse":
                            p1_coarse = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "p1 Coarse"
                            ] = p1_coarse
                        if line.split("\n")[0].split(":")[0] == "p0_fine":
                            p0_fine = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "p0 Fine"
                            ] = p0_fine
                        if line.split("\n")[0].split(":")[0] == "p1_fine":
                            p1_fine = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "p1 Fine"
                            ] = p1_fine
                        if line.split("\n")[0].split(":")[0] == "colour":
                            colour = int(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "color index"
                            ] = colour
                        if line.split("\n")[0].split(":")[0] == "reference_range":
                            reference_range = int(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "move left/right range index"
                            ] = reference_range
                        if line.split("\n")[0].split(":")[0] == "reference_value":
                            reference_value = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "move left/right"
                            ] = reference_value
                        if line.split("\n")[0].split(":")[0] == "vertical_range":
                            vertical_range = int(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "move up/down range index"
                            ] = vertical_range
                        if line.split("\n")[0].split(":")[0] == "vertical_value":
                            vertical_value = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "move up/down"
                            ] = vertical_value
                        if line.split("\n")[0].split(":")[0] == "multiply_range":
                            multiply_range = int(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "multiply range index"
                            ] = multiply_range
                        if line.split("\n")[0].split(":")[0] == "multiply_value":
                            multiply_value = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "multiply value"
                            ] = multiply_value
                        if line.split("\n")[0].split(":")[0] == "pivot_visible":
                            # Have put all the saved parameters into the dictionary
                            # Now need to add this plot to the canvas
                            self.add_saved_plot1D(count)
                            count += 1
                    self.main_frame.viewer.files.choices = self.choices

                    # Loop through all the plots and perform OnSelectPlot function
                    for i in range(len(self.main_frame.viewer.files.choices)):
                        self.main_frame.viewer.plot_combobox.SetSelection(i)
                        self.main_frame.viewer.OnSelectPlot(event=None)

                    self.main_frame.viewer.plot_combobox.SetSelection(0)
                    self.main_frame.viewer.OnSelectPlot(event=None)

                else:
                    multiplot = False
                    # Try to read in the data
                    for i, line in enumerate(lines):
                        if i < 2:
                            continue
                        else:
                            if line.split("\n")[0].split(":")[0] == "file_path":
                                file_path = line.split("\n")[0].split("file_path:")[0]
                            elif line.split("\n")[0].split(":")[0] == "p0_coarse":
                                p0_coarse = float(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "p0_fine":
                                p0_fine = float(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "p1_coarse":
                                p1_coarse = float(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "p1_fine":
                                p1_fine = float(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "colour":
                                colour = int(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "linewidth":
                                linewidth = float(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "reference_range":
                                reference_range = int(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "reference_value":
                                reference_value = float(
                                    line.split("\n")[0].split(":")[1]
                                )
                            elif line.split("\n")[0].split(":")[0] == "vertical_range":
                                vertical_range = int(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "vertical_value":
                                vertical_value = float(
                                    line.split("\n")[0].split(":")[1]
                                )
                            elif line.split("\n")[0].split(":")[0] == "multiply_range":
                                multiply_range = int(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "multiply_value":
                                multiply_value = float(
                                    line.split("\n")[0].split(":")[1]
                                )
                            elif line.split("\n")[0].split(":")[0] == "pivot_point":
                                pivot_point = float(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "pivot_x":
                                pivot_x = float(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "pivot_visible":
                                pivot_visible = line.split("\n")[0].split(":")[1]
                    self.main_frame.nmrdata = GetData(file_path)
                    self.main_frame.viewer = OneDViewer(
                        parent=self.main_frame, nmrdata=self.main_frame.nmrdata
                    )
                    self.main_frame.main_sizer.Add(self.main_frame.viewer, 1, wx.EXPAND)
                    # Update the plot with all the saved parameters
                    self.main_frame.viewer.P0_slider.SetValue(p0_coarse)
                    self.main_frame.viewer.P1_slider.SetValue(p1_coarse)
                    self.main_frame.viewer.P0_slider_fine.SetValue(p0_fine)
                    self.main_frame.viewer.P1_slider_fine.SetValue(p1_fine)
                    self.main_frame.viewer.OnSliderScroll1D(event=None)
                    self.main_frame.viewer.colour_chooser.SetSelection(colour)
                    self.main_frame.viewer.OnColourChoice1D(event=None)
                    self.main_frame.viewer.linewidth_slider.SetValue(linewidth)
                    self.main_frame.viewer.OnLinewidthScroll1D(event=None)
                    self.main_frame.viewer.reference_range_chooser.SetSelection(
                        reference_range
                    )
                    self.main_frame.viewer.OnReferenceCombo(event=None)
                    self.main_frame.viewer.reference_slider.SetValue(reference_value)
                    self.main_frame.viewer.OnReferenceScroll1D(event=None)
                    self.main_frame.viewer.vertical_range_chooser.SetSelection(
                        vertical_range
                    )
                    self.main_frame.viewer.OnVerticalCombo(event=None)
                    self.main_frame.viewer.vertical_slider.SetValue(vertical_value)
                    self.main_frame.viewer.OnVerticalScroll1D(event=None)
                    self.main_frame.viewer.multiply_range_chooser.SetSelection(
                        multiply_range
                    )
                    self.main_frame.viewer.OnMultiplyCombo(event=None)
                    self.main_frame.viewer.multiply_slider.SetValue(multiply_value)
                    self.main_frame.viewer.OnMultiplyScroll1D(event=None)
                    self.main_frame.viewer.pivot_line.set_xdata([pivot_point])
                    self.main_frame.viewer.pivot_x = pivot_x
                    if pivot_visible == "True":
                        self.main_frame.viewer.pivot_line.set_visible(True)
                    else:
                        self.main_frame.viewer.pivot_line.set_visible(False)
                    self.main_frame.viewer.OnSliderScroll1D(event=None)
            elif lines[0].split("\n")[0].split()[0] == "2D":
                # This is a 2D window
                if lines[1].split("\n")[0].split(":")[1] == "True":
                    self.multiplot_mode = True
                    transposed2D = lines[2].split("\n")[0].split(":")[1]

                    # Get the file path of the original data
                    file_path_original = lines[3].split("\n")[0].split("file_path:")[1]
                    self.main_frame.nmrdata = GetData(file_path_original)
                    self.main_frame.viewer = TwoDViewer(
                        parent=self.main_frame, nmrdata=self.main_frame.nmrdata
                    )
                    self.main_frame.main_sizer.Add(self.main_frame.viewer, 1, wx.EXPAND)
                    title = lines[4].split("\n")[0].split(":")[1]
                    p0_coarse = float(lines[5].split("\n")[0].split(":")[1])
                    p0_fine = float(lines[6].split("\n")[0].split(":")[1])
                    p1_coarse = float(lines[7].split("\n")[0].split(":")[1])
                    p1_fine = float(lines[8].split("\n")[0].split(":")[1])
                    move_x = float(lines[9].split("\n")[0].split(":")[1])
                    move_y = float(lines[10].split("\n")[0].split(":")[1])
                    move_x_index = int(lines[11].split("\n")[0].split(":")[1])
                    move_y_index = int(lines[12].split("\n")[0].split(":")[1])
                    contour_linewidth = float(lines[13].split("\n")[0].split(":")[1])
                    multiply_factor = float(lines[14].split("\n")[0].split(":")[1])
                    contour_levels = int(lines[15].split("\n")[0].split(":")[1])
                    transposed = lines[16].split("\n")[0].split(":")[1]
                    self.choices = [title]
                    self.main_frame.viewer.contour_width_slider.SetValue(
                        contour_linewidth
                    )
                    self.main_frame.viewer.contour_levels_slider.SetValue(
                        contour_levels
                    )
                    self.main_frame.viewer.multiply_slider.SetValue(multiply_factor)
                    self.main_frame.viewer.reference_range_chooserX.SetSelection(
                        move_x_index
                    )
                    self.main_frame.viewer.reference_range_chooserY.SetSelection(
                        move_y_index
                    )
                    self.main_frame.viewer.move_x_slider.SetValue(move_x)
                    self.main_frame.viewer.move_y_slider.SetValue(move_y)
                    self.main_frame.viewer.P0_slider.SetValue(p0_coarse)
                    self.main_frame.viewer.P1_slider.SetValue(p1_coarse)
                    self.main_frame.viewer.P0_slider_fine.SetValue(p0_fine)
                    self.main_frame.viewer.P1_slider_fine.SetValue(p1_fine)

                    self.main_frame.viewer.files.first_drop = False
                    self.main_frame.viewer.OnSliderScroll2D(event=None)
                    # Read all the values into the values dictionary for the first plot
                    self.main_frame.viewer.values_dictionary[0] = {}
                    self.main_frame.viewer.values_dictionary[0][
                        "path"
                    ] = file_path_original
                    self.main_frame.viewer.values_dictionary[0]["title"] = title
                    self.main_frame.viewer.values_dictionary[0][
                        "contour linewidth"
                    ] = contour_linewidth
                    self.main_frame.viewer.values_dictionary[0][
                        "contour levels"
                    ] = contour_levels
                    self.main_frame.viewer.values_dictionary[0]["move x"] = move_x
                    self.main_frame.viewer.values_dictionary[0][
                        "move x range index"
                    ] = move_x_index
                    self.main_frame.viewer.values_dictionary[0]["move y"] = move_y
                    self.main_frame.viewer.values_dictionary[0][
                        "move y range index"
                    ] = move_y_index
                    self.main_frame.viewer.values_dictionary[0][
                        "multiply factor"
                    ] = multiply_factor
                    self.main_frame.viewer.values_dictionary[0]["p0 Coarse"] = p0_coarse
                    self.main_frame.viewer.values_dictionary[0]["p1 Coarse"] = p1_coarse
                    self.main_frame.viewer.values_dictionary[0]["p0 Fine"] = p0_fine
                    self.main_frame.viewer.values_dictionary[0]["p1 Fine"] = p1_fine
                    self.main_frame.viewer.values_dictionary[0][
                        "transposed"
                    ] = transposed
                    self.main_frame.viewer.values_dictionary[0][
                        "original_x_ppms"
                    ] = self.main_frame.viewer.ppms_0
                    self.main_frame.viewer.values_dictionary[0][
                        "original_y_ppms"
                    ] = self.main_frame.viewer.ppms_1
                    self.main_frame.viewer.values_dictionary[0]["new_x_ppms"] = (
                        self.main_frame.viewer.ppms_0
                        + np.ones(len(self.main_frame.viewer.ppms_0)) * move_x
                    )
                    self.main_frame.viewer.values_dictionary[0]["new_y_ppms"] = (
                        self.main_frame.viewer.ppms_1
                        + np.ones(len(self.main_frame.viewer.ppms_1)) * move_y
                    )
                    self.main_frame.viewer.values_dictionary[0][
                        "z_data"
                    ] = self.main_frame.viewer.nmrdata.data
                    self.main_frame.viewer.values_dictionary[0][
                        "uc0"
                    ] = self.main_frame.viewer.uc0
                    self.main_frame.viewer.values_dictionary[0][
                        "uc1"
                    ] = self.main_frame.viewer.uc1
                    self.main_frame.viewer.values_dictionary[0]["linewidth 1D"] = 1.0
                    self.main_frame.viewer.multiplot_mode = True
                    count = 1
                    self.main_frame.viewer.values_dictionary[count] = {}
                    # Loop over the rest of the lines to get the file paths of the other data
                    for i, line in enumerate(lines):
                        if i < 17:
                            continue
                        if "file_path:" in line.split("\n")[0]:
                            file_path = line.split("\n")[0].split("file_path:")[1]
                            self.main_frame.viewer.values_dictionary[count][
                                "path"
                            ] = file_path
                        elif line.split("\n")[0].split(":")[0] == "title":
                            title = line.split("\n")[0].split(":")[1]
                            self.main_frame.viewer.values_dictionary[count][
                                "title"
                            ] = title
                        elif line.split("\n")[0].split(":")[0] == "contour linewidth":
                            contour_linewidth = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "contour linewidth"
                            ] = contour_linewidth
                        elif line.split("\n")[0].split(":")[0] == "move x":
                            move_x = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "move x"
                            ] = move_x
                        elif line.split("\n")[0].split(":")[0] == "move x range index":
                            move_x_index = int(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "move x range index"
                            ] = move_x_index
                        elif line.split("\n")[0].split(":")[0] == "move y":
                            move_y = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "move y"
                            ] = move_y
                        elif line.split("\n")[0].split(":")[0] == "move y range index":
                            move_y_index = int(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "move y range index"
                            ] = move_y_index
                        elif line.split("\n")[0].split(":")[0] == "multiply factor":
                            multiply_factor = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "multiply factor"
                            ] = multiply_factor
                        elif line.split("\n")[0].split(":")[0] == "p0 Coarse":
                            p0_coarse = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "p0 Coarse"
                            ] = p0_coarse
                        elif line.split("\n")[0].split(":")[0] == "p1 Coarse":
                            p1_coarse = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "p1 Coarse"
                            ] = p1_coarse
                        elif line.split("\n")[0].split(":")[0] == "p0 Fine":
                            p0_fine = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "p0 Fine"
                            ] = p0_fine
                        elif line.split("\n")[0].split(":")[0] == "p1 Fine":
                            p1_fine = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "p1 Fine"
                            ] = p1_fine
                        elif line.split("\n")[0].split(":")[0] == "contour levels":
                            contour_levels = int(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "contour levels"
                            ] = contour_levels
                        elif line.split("\n")[0].split(":")[0] == "transposed":
                            self.main_frame.viewer.values_dictionary[count][
                                "transposed"
                            ] = line.split("\n")[0].split(":")[1]
                            self.main_frame.viewer.values_dictionary[count][
                                "linewidth 1D"
                            ] = 1.0
                            self.add_saved_plot2D(count)
                            count += 1
                            self.main_frame.viewer.values_dictionary[count] = {}

                    # Search thrugh values dictionary and remove empty entries
                    keys = list(self.main_frame.viewer.values_dictionary.keys())
                    for key in keys:
                        if self.main_frame.viewer.values_dictionary[key] == {}:
                            del self.main_frame.viewer.values_dictionary[key]

                    self.plot_overlaid_2D()

            elif lines[0].split("\n")[0].split()[0] == "3D":
                # This is a 3D window
                pass
            else:
                # Give a popout saying that the session file is not formatted correctly
                msg = wx.MessageDialog(
                    self,
                    "The session file is not formatted correctly",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                return

    def add_saved_plot1D(self, count):
        # Add the saved plot to the canvas
        # Read in the data
        dic, data = ng.pipe.read(
            self.main_frame.viewer.values_dictionary[count]["path"]
        )
        self.main_frame.viewer.values_dictionary[count]["original_data"] = data
        self.main_frame.viewer.values_dictionary[count]["dictionary"] = dic
        # Make the uc object
        uc0 = ng.pipe.make_uc(dic, data)

        # Get the ppm scale
        ppm_scale = uc0.ppm_scale()
        self.main_frame.viewer.values_dictionary[count]["original_ppms"] = ppm_scale

        data = data * self.main_frame.viewer.values_dictionary[count][
            "multiply value"
        ] + self.main_frame.viewer.values_dictionary[count]["move up/down"] * np.ones(
            len(data)
        )

        self.choices.append(self.main_frame.viewer.values_dictionary[count]["title"])
        self.main_frame.viewer.plot_combobox.Clear()
        self.main_frame.viewer.plot_combobox.AppendItems(self.choices)
        self.main_frame.viewer.plot_combobox.SetSelection(0)
        xlim, ylim = (
            self.main_frame.viewer.ax.get_xlim(),
            self.main_frame.viewer.ax.get_ylim(),
        )
        self.main_frame.viewer.extra_plots.append(
            self.main_frame.viewer.ax.plot(
                uc0.ppm_scale(),
                data,
                color=self.main_frame.viewer.colours[
                    self.main_frame.viewer.values_dictionary[count]["color index"]
                ],
                label=self.choices[-1],
                linewidth=self.main_frame.viewer.values_dictionary[count]["linewidth"],
            )
        )

        self.main_frame.viewer.ax.legend()
        self.main_frame.viewer.ax.set_xlim(xlim)
        self.main_frame.viewer.ax.set_ylim(ylim)
        self.main_frame.viewer.OnSelectPlot(wx.EVT_COMBOBOX)

    def add_saved_plot2D(self, count):
        # Add the saved plot to the canvas
        # Read in the data
        dic, data = ng.pipe.read(
            self.main_frame.viewer.values_dictionary[count]["path"]
        )
        self.main_frame.viewer.values_dictionary[count]["original_data"] = data
        self.main_frame.viewer.values_dictionary[count]["dictionary"] = dic
        # Make the uc object
        uc0 = ng.pipe.make_uc(dic, data, dim=0)
        uc1 = ng.pipe.make_uc(dic, data, dim=1)
        ppm0 = uc0.ppm_scale()
        ppm1 = uc1.ppm_scale()
        x, y = np.meshgrid(ppm1, ppm0)

        self.main_frame.viewer.values_dictionary[count]["original_x_ppms"] = ppm0
        self.main_frame.viewer.values_dictionary[count]["original_y_ppms"] = ppm1
        self.main_frame.viewer.values_dictionary[count]["new_x_ppms"] = ppm0
        self.main_frame.viewer.values_dictionary[count]["new_y_ppms"] = ppm1
        self.main_frame.viewer.values_dictionary[count]["z_data"] = data
        self.main_frame.viewer.values_dictionary[count]["contour linewidth"] = 1.0
        self.main_frame.viewer.values_dictionary[count]["linewidth 1D"] = 1.0
        self.main_frame.viewer.values_dictionary[count]["uc0"] = uc0
        self.main_frame.viewer.values_dictionary[count]["uc1"] = uc1

        length = len(self.main_frame.viewer.values_dictionary.keys())

        # If transpose is false, then the x-axis is the first axis and the y-axis is the second axis
        self.main_frame.viewer.values_dictionary[count]["new_x_ppms_old"] = (
            self.main_frame.viewer.values_dictionary[count]["new_x_ppms"]
        )
        self.main_frame.viewer.values_dictionary[count]["new_y_ppms_old"] = (
            self.main_frame.viewer.values_dictionary[count]["new_y_ppms"]
        )
        transposed = self.main_frame.viewer.values_dictionary[count]["transposed"]

        if transposed == "True":
            self.main_frame.viewer.values_dictionary[count]["new_x_ppms_old"] = (
                self.main_frame.viewer.values_dictionary[count]["new_x_ppms"]
            )
            self.main_frame.viewer.values_dictionary[count]["new_y_ppms_old"] = (
                self.main_frame.viewer.values_dictionary[count]["new_y_ppms"]
            )
            self.main_frame.viewer.values_dictionary[count]["new_x_ppms"] = (
                self.main_frame.viewer.values_dictionary[count]["new_y_ppms_old"]
            )
            self.main_frame.viewer.values_dictionary[count]["new_y_ppms"] = (
                self.main_frame.viewer.values_dictionary[count]["new_x_ppms_old"]
            )
            self.main_frame.viewer.values_dictionary[count]["original_x_ppms"] = (
                self.main_frame.viewer.values_dictionary[count]["new_x_ppms"]
            )
            self.main_frame.viewer.values_dictionary[count]["original_y_ppms"] = (
                self.main_frame.viewer.values_dictionary[count]["original_y_ppms"]
            )
            uc0 = self.main_frame.viewer.values_dictionary[count]["uc1"]
            uc1 = self.main_frame.viewer.values_dictionary[count]["uc0"]
            self.main_frame.viewer.values_dictionary[count]["uc0"] = uc0
            self.main_frame.viewer.values_dictionary[count]["uc1"] = uc1
            self.main_frame.viewer.values_dictionary[count]["z_data"] = (
                self.main_frame.viewer.values_dictionary[count]["z_data"].T
            )

    def plot_overlaid_2D(self):
        xlim, ylim = (
            self.main_frame.viewer.ax.get_xlim(),
            self.main_frame.viewer.ax.get_ylim(),
        )
        xlabel = self.main_frame.viewer.ax.get_xlabel()
        ylabel = self.main_frame.viewer.ax.get_ylabel()
        self.main_frame.viewer.ax.clear()
        self.main_frame.viewer.axes1D.clear()
        self.main_frame.viewer.axes1D_2.clear()
        self.main_frame.viewer.axes1D.set_yticks([])
        self.main_frame.viewer.axes1D_2.set_xticks([])
        self.main_frame.viewer.twoD_spectra = []
        self.main_frame.viewer.twoD_slices_horizontal = []
        self.main_frame.viewer.twoD_slices_vertical = []

        for i in range(len(self.main_frame.viewer.values_dictionary)):
            multiply_factor = self.main_frame.viewer.values_dictionary[i][
                "multiply factor"
            ]

            x, y = np.meshgrid(
                self.main_frame.viewer.values_dictionary[i]["new_y_ppms"],
                self.main_frame.viewer.values_dictionary[i]["new_x_ppms"],
            )
            self.main_frame.viewer.twoD_spectra.append(
                self.main_frame.viewer.ax.contour(
                    y,
                    x,
                    self.main_frame.viewer.values_dictionary[i]["z_data"]
                    * multiply_factor,
                    colors=self.main_frame.viewer.twoD_colours[i],
                    levels=self.main_frame.viewer.cl,
                    linewidths=self.main_frame.viewer.values_dictionary[i][
                        "contour linewidth"
                    ],
                )
            )

            if self.main_frame.viewer.transposed2D == False:
                self.main_frame.viewer.twoD_slices_horizontal.append(
                    self.main_frame.viewer.axes1D.plot(
                        self.main_frame.viewer.values_dictionary[i]["new_x_ppms"],
                        self.main_frame.viewer.values_dictionary[i]["z_data"][:, 1]
                        * multiply_factor,
                        color=self.main_frame.viewer.twoD_label_colours[i],
                        linewidth=self.main_frame.viewer.values_dictionary[i][
                            "linewidth 1D"
                        ],
                    )
                )
                self.main_frame.viewer.twoD_slices_vertical.append(
                    self.main_frame.viewer.axes1D_2.plot(
                        self.main_frame.viewer.values_dictionary[i]["new_y_ppms"],
                        self.main_frame.viewer.values_dictionary[i]["z_data"][1, :]
                        * multiply_factor,
                        color=self.main_frame.viewer.twoD_label_colours[i],
                        linewidth=self.main_frame.viewer.values_dictionary[i][
                            "linewidth 1D"
                        ],
                    )
                )

            else:
                if i == 0:
                    self.main_frame.viewer.twoD_slices_horizontal.append(
                        self.main_frame.viewer.axes1D.plot(
                            self.main_frame.viewer.values_dictionary[i]["new_x_ppms"],
                            self.main_frame.viewer.values_dictionary[i]["z_data"].T[
                                1, :
                            ]
                            * multiply_factor,
                            color=self.main_frame.viewer.twoD_label_colours[i],
                            linewidth=self.main_frame.viewer.values_dictionary[i][
                                "linewidth 1D"
                            ],
                        )
                    )
                    self.main_frame.viewer.twoD_slices_vertical.append(
                        self.main_frame.viewer.axes1D_2.plot(
                            self.main_frame.viewer.values_dictionary[i]["new_y_ppms"],
                            self.main_frame.viewer.values_dictionary[i]["z_data"].T[
                                :, 1
                            ]
                            * multiply_factor,
                            color=self.main_frame.viewer.twoD_label_colours[i],
                            linewidth=self.main_frame.viewer.values_dictionary[i][
                                "linewidth 1D"
                            ],
                        )
                    )
                else:
                    self.main_frame.viewer.twoD_slices_horizontal.append(
                        self.main_frame.viewer.axes1D.plot(
                            self.main_frame.viewer.values_dictionary[i]["new_x_ppms"],
                            self.main_frame.viewer.values_dictionary[i]["z_data"].T[
                                :, 1
                            ]
                            * multiply_factor,
                            color=self.main_frame.viewer.twoD_label_colours[i],
                            linewidth=self.main_frame.viewer.values_dictionary[i][
                                "linewidth 1D"
                            ],
                        )
                    )
                    self.main_frame.viewer.twoD_slices_vertical.append(
                        self.main_frame.viewer.axes1D_2.plot(
                            self.main_frame.viewer.values_dictionary[i]["new_y_ppms"],
                            self.main_frame.viewer.values_dictionary[i]["z_data"].T[
                                1, :
                            ]
                            * multiply_factor,
                            color=self.main_frame.viewer.twoD_label_colours[i],
                            linewidth=self.main_frame.viewer.values_dictionary[i][
                                "linewidth 1D"
                            ],
                        )
                    )
        self.main_frame.viewer.line_h = self.main_frame.viewer.ax.axhline(
            y=self.main_frame.viewer.values_dictionary[i]["new_x_ppms"][1],
            color="black",
            lw=1.5,
        )
        self.main_frame.viewer.line_v = self.main_frame.viewer.ax.axvline(
            x=self.main_frame.viewer.values_dictionary[i]["new_y_ppms"][1],
            color="black",
            lw=1.5,
        )
        self.main_frame.viewer.line_h.set_visible(False)
        self.main_frame.viewer.line_v.set_visible(False)
        self.custom_labels = []

        for i in range(len(self.main_frame.viewer.values_dictionary)):
            self.custom_labels.append(
                self.main_frame.viewer.values_dictionary[i]["title"]
            )
            self.main_frame.viewer.files.custom_lines.append(
                Line2D(
                    [0], [0], color=self.main_frame.viewer.twoD_label_colours[i], lw=1.5
                )
            )

        # Set all vertical and horizontal slices to invisible initia
        for i in range(len(self.main_frame.viewer.twoD_slices_horizontal)):
            self.main_frame.viewer.twoD_slices_horizontal[i][0].set_visible(False)
            self.main_frame.viewer.twoD_slices_vertical[i][0].set_visible(False)

        self.main_frame.viewer.custom_labels = self.custom_labels
        self.main_frame.viewer.ax.legend(
            self.main_frame.viewer.files.custom_lines, self.custom_labels
        )
        self.main_frame.viewer.ax.set_xlim(xlim)
        self.main_frame.viewer.ax.set_ylim(ylim)
        self.main_frame.viewer.ax.set_xlabel(xlabel)
        self.main_frame.viewer.ax.set_ylabel(ylabel)

        # Add labels of the extra plots to the select plot box
        self.main_frame.viewer.plot_combobox.Clear()
        self.main_frame.viewer.plot_combobox.AppendItems(self.custom_labels)
        self.main_frame.viewer.plot_combobox.SetSelection(0)

        self.main_frame.viewer.UpdateFrame()


def main():
    app = wx.App()
    frame = MyApp()
    app.MainLoop()


if __name__ == "__main__":
    main()
