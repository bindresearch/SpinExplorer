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


print("-------------------------------------------------------------")
print("                         SpinProcess                         ")
print("-------------------------------------------------------------")
print("                (version 1.2) 20th June 2025                 ")
print(" (c) 2025 James Eaton, Andrew Baldwin (University of Oxford) ")
print("                        MIT License                          ")
print("-------------------------------------------------------------")
print("                     Processing NMR Data                     ")
print("-------------------------------------------------------------")
print(" Documentation at:")
print(" https://github.com/james-eaton-1/SpinExplorer")
print("-------------------------------------------------------------")
print("")


import sys

import wx
import wx.lib.agw.hyperlink as hl

# Import relevant modules
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("WXAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import (
    NavigationToolbar2WxAgg as NavigationToolbar,
)
import nmrglue as ng
import subprocess
import os

matplotlib.rcParams["font.sans-serif"] = "Arial"
matplotlib.rcParams["font.family"] = "sans-serif"

# Suppress complex warning from numpy
import warnings

# warnings.simplefilter("ignore", np.ComplexWarning)  # For old numpy versions
warnings.simplefilter("ignore", np.exceptions.ComplexWarning)  # For new numpy versions


# Importing SpinProcess modules
from SpinExplorer.SpinProcess.ReadingData.read_fid import ReadFID
from SpinExplorer.SpinProcess.FormattingGUI.notebook import NotebookProcess


# James Eaton, 10/06/2025, University of Oxford
# This program is designed to allow the user to process NMR FID data that has been converted to nmrPipe format.


# Read the FID data from the nmrPipe file


class SpinProcess(wx.Frame):
    def __init__(
        self, original_frame=None, file_parser=False, path="", cwd="", reprocess=False
    ):
        # Get the monitor size and set the window size to 85% of the monitor size
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 0.7 * self.monitorWidth
        self.height = 0.75 * self.monitorHeight

        # Read the NMR data in the current directory
        self.nmr_data = ReadFID(self)

        # Initially set the reprocessing flag to False
        self.reprocess = reprocess
        self.original_frame = original_frame
        self.file_parser = file_parser
        self.path = path
        self.cwd = cwd

        # Create the main window
        self.main_window = wx.Frame.__init__(
            self, None, title="SpinProcess", size=(self.width, self.height)
        )

        self.notebook = NotebookProcess(self, self.nmr_data)
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.main_sizer.AddSpacer(10)
        self.main_sizer.Add(self.notebook, 1, wx.EXPAND)
        self.notebook.create_buttons(parent=self)

        self.SetSizerAndFit(self.main_sizer)
        # self.SetWindowStyle(wx.STAY_ON_TOP)
        self.Show()
        self.Centre()

        self.Bind(wx.EVT_CLOSE, self.OnClose)

    def OnClose(self, event):
        self.Destroy()
        sys.exit()

    def change_frame_size(self, width, height):
        self.SetSize(width, height)

        # Centre the window on the screen
        self.Centre()


class OneDFrame(wx.Panel):

    def __init__(self, parent):
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 0.7 * self.monitorWidth
        self.height = 0.75 * self.monitorHeight
        self.parent = parent
        wx.Panel.__init__(self, parent, id=wx.ID_ANY, size=(self.width, self.height))
        # Create panel for processing dimension 1 of the data
        self.nmr_data = parent.nmr_data
        self.set_variables()
        self.create_canvas()
        self.create_menu_bar()

    def set_variables(self):

        # See if NMR processing file (nmrproc.com) can be found, if it can try to load the variables from it
        if os.path.exists("processing_parameters.txt"):
            found_nmrproc_com = True
        else:
            found_nmrproc_com = False

        self.set_initial_solvent_suppression_variables()
        self.set_initial_linear_prediction_variables()
        self.set_initial_apodization_variables()
        self.set_initial_zero_filling_variables()
        self.set_initial_fourier_transform_variables()
        self.set_initial_phasing_variables()
        self.set_initial_extraction_variables()
        self.set_initial_baseline_correction_variables()

        self.parent.load_variables = False
        if found_nmrproc_com == False:
            pass
        else:
            # Ask the user if they want to load the variables from the nmrproc.com file
            dlg = wx.MessageDialog(
                self,
                "A file containing NMR processing parameters has been found (processing_parameters.txt). Do you want to load the variables from it?",
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

    def create_canvas(self):

        pass

    def create_menu_bar(self):
        # Create the main sizer
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_sizer)

        self.sizer_1 = wx.BoxSizer(wx.VERTICAL)
        self.sizer_1.AddSpacer(10)

        # Create all the sizers
        self.create_solvent_suppression_sizer(parent=self)
        self.create_linear_prediction_sizer(parent=self)
        self.create_apodization_sizer(parent=self)
        self.create_zero_filling_sizer(parent=self)
        self.create_fourier_transform_sizer(parent=self)
        self.create_phase_correction_sizer(parent=self)
        self.create_extraction_sizer(parent=self)
        self.create_baseline_correction_sizer(parent=self)

        self.main_sizer.Add(self.sizer_1, 0, wx.EXPAND)

        self.SetSizerAndFit(self.main_sizer)
        self.Layout()

        # Get the size of the main sizer and set the window size to 1.05 times the size of the main sizer
        self.width, self.height = self.main_sizer.GetSize()
        self.parent.parent.change_frame_size(
            int(self.width * 1.05), int(self.height * 1.25)
        )


class TwoDFrame(wx.Panel):
    def __init__(self, parent, oneDFrame):
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 0.7 * self.monitorWidth
        self.height = 0.75 * self.monitorHeight
        self.parent = parent
        wx.Panel.__init__(self, parent, id=wx.ID_ANY, size=(self.width, self.height))

        self.oneDFrame = oneDFrame

        # Create panel for processing dimension 1 of the data
        self.nmr_data = parent.nmr_data
        self.set_variables_dim2()
        self.create_canvas_dim2()
        self.create_menu_bar_dim2()

    def set_variables_dim2(self):
        self.set_initial_linear_prediction_variables_dim2()
        self.set_initial_apodization_variables_dim2()
        self.set_initial_zero_filling_variables_dim2()
        self.set_initial_fourier_transform_variables_dim2()
        self.set_initial_phasing_variables_dim2()
        self.set_initial_extraction_variables_dim2()
        self.set_initial_baseline_correction_variables_dim2()

        if self.parent.load_variables == True:
            try:
                self.load_variables_from_nmrproc_com_2D()
            except:
                pass

    def load_variables_from_nmrproc_com_2D(self):
        # Open processing_parameters.txt file and load the variables from it
        file = open("processing_parameters.txt", "r")
        lines = file.readlines()
        file.close()

        include_line = False
        for line in lines:
            if "Dimension 2" in line:
                include_line = True
                continue
            if include_line == False:
                continue
            if include_line == True and "Dimension 3" in line:
                include_line = False
                break
            if include_line == True:
                line = line.split("\n")[0]
                if line.split(":")[0] == "Linear Prediction":
                    self.linear_prediction_radio_box_dim2_selection = int(
                        line.split(": ")[1]
                    )
                if line.split(":")[0] == "Linear Prediction Options Selection":
                    self.linear_prediction_dim2_options_selection = int(
                        line.split(": ")[1]
                    )
                if line.split(":")[0] == "Linear Prediction Coefficients Selection":
                    self.linear_prediction_dim2_coefficients_selection = int(
                        line.split(": ")[1]
                    )
                if line.split(":")[0] == "NUS file":
                    self.nuslist_name_dim2 = line.split(": ")[1]
                if line.split(":")[0] == "NUS CPU":
                    self.smile_nus_cpu_textcontrol_dim2 = int(line.split(": ")[1])
                if line.split(":")[0] == "NUS Iterations":
                    self.smile_nus_iterations_textcontrol_dim2 = int(
                        line.split(": ")[1]
                    )
                if line.split(":")[0] == "Apodization":
                    if "True" in line:
                        self.apodization_dim2_checkbox_value = True
                    else:
                        self.apodization_dim2_checkbox_value = False
                if line.split(":")[0] == "Apodization Combobox Selection":
                    self.apodization_dim2_combobox_selection = int(line.split(": ")[1])
                    self.apodization_dim2_combobox_selection_old = int(
                        line.split(": ")[1]
                    )
                if line.split(":")[0] == "Exponential Line Broadening":
                    self.exponential_line_broadening_dim2 = float(line.split(": ")[1])
                if line.split(":")[0] == "Apodization First Point Scaling":
                    self.apodization_first_point_scaling_dim2 = float(
                        line.split(": ")[1]
                    )
                if line.split(":")[0] == "G1":
                    self.g1_dim2 = float(line.split(": ")[1])
                if line.split(":")[0] == "G2":
                    self.g2_dim2 = float(line.split(": ")[1])
                if line.split(":")[0] == "G3":
                    self.g3_dim2 = float(line.split(": ")[1])
                if line.split(":")[0] == "Offset":
                    self.offset_dim2 = float(line.split(": ")[1])
                if line.split(":")[0] == "End":
                    self.end_dim2 = float(line.split(": ")[1])
                if line.split(":")[0] == "Power":
                    self.power_dim2 = float(line.split(": ")[1])
                if line.split(":")[0] == "A":
                    self.a_dim2 = float(line.split(": ")[1])
                if line.split(":")[0] == "B":
                    self.b_dim2 = float(line.split(": ")[1])
                if line.split(":")[0] == "T1":
                    self.t1_dim2 = float(line.split(": ")[1])
                if line.split(":")[0] == "T2":
                    self.t2_dim2 = float(line.split(": ")[1])
                if line.split(":")[0] == "Loc":
                    self.loc_dim2 = float(line.split(": ")[1])

                if line.split(":")[0] == "Zero Filling":
                    if "True" in line:
                        self.zero_filling_checkbox_dim2_value = True
                    else:
                        self.zero_filling_checkbox_dim2_value = False
                if line.split(":")[0] == "Zero Filling Combobox Selection":
                    self.zero_filling_dim2_combobox_selection = int(line.split(": ")[1])
                if line.split(":")[0] == "Zero Filling Value Doubling Times":
                    self.zero_filling_dim2_value_doubling_times = int(
                        line.split(": ")[1]
                    )
                if line.split(":")[0] == "Zero Filling Value Zeros to Add":
                    self.zero_filling_dim2_value_zeros_to_add = int(line.split(": ")[1])
                if line.split(":")[0] == "Zero Filling Value Final Data Size":
                    self.zero_filling_dim2_value_final_data_size = int(
                        line.split(": ")[1]
                    )
                if line.split(":")[0] == "Zero Filling Round Checkbox":
                    if "True" in line:
                        self.zero_filling_round_checkbox_dim2_value = True
                    else:
                        self.zero_filling_round_checkbox_dim2_value = False
                if line.split(":")[0] == "Fourier Transform":
                    if "True" in line:
                        self.fourier_transform_checkbox_dim2_value = True
                    else:
                        self.fourier_transform_checkbox_dim2_value = False
                if line.split(":")[0] == "Fourier Transform Method Selection":
                    self.ft_method_selection_dim2 = int(line.split(": ")[1])
                if line.split(":")[0] == "Phase Correction":
                    if "True" in line:
                        self.phase_correction_checkbox_dim2_value = True
                    else:
                        self.phase_correction_checkbox_dim2_value = False
                if line.split(":")[0] == "Phase Correction P0":
                    self.p0_total_dim2 = float(line.split(": ")[1])
                if line.split(":")[0] == "Phase Correction P1":
                    self.p1_total_dim2 = float(line.split(": ")[1])
                if line.split(":")[0] == "F1180":
                    if "True" in line:
                        self.f1180_dim2 = True
                    else:
                        self.f1180_dim2 = False
                if line.split(":")[0] == "Extraction":
                    if "True" in line:
                        self.extraction_checkbox_dim2_value = True
                    else:
                        self.extraction_checkbox_dim2_value = False
                if line.split(":")[0] == "Extraction PPM Start":
                    self.extraction_ppm_start_dim2 = float(line.split(": ")[1])
                if line.split(":")[0] == "Extraction PPM End":
                    self.extraction_ppm_end_dim2 = float(line.split(": ")[1])
                if line.split(":")[0] == "Baseline Correction":
                    if "True" in line:
                        self.baseline_correction_checkbox_dim2_value = True
                    else:
                        self.baseline_correction_checkbox_dim2_value = False
                if line.split(":")[0] == "Baseline Correction Radio Box Selection":
                    self.baseline_correction_radio_box_selection_dim2 = int(
                        line.split(": ")[1]
                    )
                if line.split(":")[0] == "Baseline Correction Nodes":
                    self.baseline_correction_nodes_dim2 = int(line.split(": ")[1])
                if line.split(":")[0] == "Baseline Correction Node List":
                    self.baseline_correction_node_list_dim2 = line.split(": ")[1]
                if line.split(":")[0] == "Baseline Correction Polynomial Order":
                    self.baseline_correction_polynomial_order_dim2 = int(
                        line.split(": ")[1]
                    )

    def set_initial_linear_prediction_variables_dim2(self):
        self.linear_prediction_radio_box_dim2_selection = 1
        self.linear_prediction_dim2_checkbox_value = False
        self.linear_prediction_dim2_options_selection = 0
        self.linear_prediction_dim2_coefficients_selection = 0
        self.linear_prediction_selection = 0

        # Check to see if the nuslist file exists in the current directory using os.path.isfile('nuslist')
        if os.path.isfile("nuslist"):
            self.nuslist_name_dim2 = "nuslist"
        else:
            self.nuslist_name_dim2 = ""

        self.number_of_nus_CPU_dim2 = 1
        self.nus_iterations_dim2 = 50
        self.smile_data_extension_number_dim2 = (
            0  # int(self.nmr_data.number_of_points[1]*1.5)
        )

    def set_initial_apodization_variables_dim2(self):
        self.apodization_dim2_checkbox_value = True
        self.apodization_dim2_combobox_selection = 1
        self.apodization_dim2_combobox_selection_old = 1

        # Initial values for exponential apodization
        self.exponential_line_broadening_dim2 = 0.5
        self.apodization_first_point_scaling_dim2 = 0.5

        # Initial values for Lorentz to Gauss apodization
        self.g1_dim2 = 0.33
        self.g2_dim2 = 1
        self.g3_dim2 = 0.0

        # Initial values for Sinebell apodization
        self.offset_dim2 = 0.5
        self.end_dim2 = 0.98
        self.power_dim2 = 1.0

        # Initial values for Gauss Broadening apodization
        self.a_dim2 = 1.0
        self.b_dim2 = 1.0

        # Initial values for Trapezoid apodization
        self.t1_dim2 = int((self.nmr_data.number_of_points[1] / 2) / 4)
        self.t2_dim2 = int((self.nmr_data.number_of_points[1] / 2) / 4)

        # Initial values for Triangle apodization
        self.loc_dim2 = 0.5

    def set_initial_zero_filling_variables_dim2(self):
        self.zero_filling_dim2_checkbox_value = True
        self.zero_filling_dim2_combobox_selection = 0
        self.zero_filling_dim2_combobox_selection_old = 0
        self.zero_filling_dim2_value_doubling_times = 1
        self.zero_filling_dim2_value_zeros_to_add = 0
        self.zero_filling_dim2_value_final_data_size = 0
        self.zero_filling_dim2_round_checkbox_value = True

    def set_initial_fourier_transform_variables_dim2(self):
        self.fourier_transform_dim2_checkbox_value = True
        self.ft_method_selection_dim2 = (
            0  # Initially use the 'auto' method of FT as default
        )

    def set_initial_phasing_variables_dim2(self):
        self.phasing_dim2_checkbox_value = True
        self.p0_total_dim2 = 0.0
        self.p1_total_dim2 = 0.0
        self.p0_total_dim2_old = 0.0
        self.p1_total_dim2_old = 0.0
        self.phasing_from_smile = False
        self.f1180_dim2 = False

    def set_initial_extraction_variables_dim2(self):
        self.extraction_checkbox_value_dim2 = False
        self.extraction_start_dim2 = "0.0"
        self.extraction_end_dim2 = "0.0"

    def set_initial_baseline_correction_variables_dim2(self):
        self.baseline_correction_checkbox_value_dim2 = False
        self.baseline_correction_radio_box_selection_dim2 = 0
        self.node_list_dim2 = "0,5,95,100"
        self.node_width_dim2 = "2"
        self.polynomial_order_dim2 = "4"

    def create_canvas_dim2(self):

        pass

    def create_menu_bar_dim2(self):
        # Create the main sizer
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_sizer)

        # Create a sizer for the processing options for the first dimension
        self.sizer_2 = wx.BoxSizer(wx.VERTICAL)
        self.sizer_2.AddSpacer(10)

        # Create all the sizers (allow a checkbox at the top for SMILE NUS reconstruction which will change the possible options)
        # For NUS reconstruction using SMILE need to have exact phasing paramaters (first process without NUS and then calculate phase in indirect dimension, then process again using SMILE containing exact phasing parameters)
        self.create_linear_prediction_sizer_dim2(parent=self)
        self.create_apodization_sizer_dim2(parent=self)
        self.create_zero_filling_sizer_dim2(parent=self)
        self.create_fourier_transform_sizer_dim2(parent=self)
        self.create_phase_correction_sizer_dim2(parent=self)
        self.create_extraction_sizer_dim2(parent=self)
        self.create_baseline_correction_sizer_dim2(parent=self)

        self.main_sizer.Add(self.sizer_2, 0, wx.EXPAND)

        self.SetSizerAndFit(self.main_sizer)
        self.Layout()

        # Get the size of the main sizer and set the window size to 1.05 times the size of the main sizer
        self.width, self.height = self.main_sizer.GetSize()
        self.parent.parent.change_frame_size(
            int(self.width * 1.05), int(self.height * 1.25)
        )

    def create_apodization_sizer_dim2(self, parent):

        # Create a box for apodization options
        self.apodization_box_dim2 = wx.StaticBox(parent, -1, "Apodization")
        self.apodization_sizer_dim2 = wx.StaticBoxSizer(
            self.apodization_box_dim2, wx.HORIZONTAL
        )
        self.apodization_checkbox_dim2 = wx.CheckBox(parent, -1, "Apply apodization")
        self.apodization_checkbox_dim2.Bind(
            wx.EVT_CHECKBOX, self.on_apodization_checkbox_dim2
        )
        self.apodization_checkbox_dim2.SetValue(self.apodization_dim2_checkbox_value)
        self.apodization_sizer_dim2.Add(
            self.apodization_checkbox_dim2, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.apodization_sizer_dim2.AddSpacer(10)
        # Have a combobox for apodization options
        self.apodization_options_dim2 = [
            "None",
            "Exponential",
            "Lorentz to Gauss",
            "Sinebell",
            "Gauss Broadening",
            "Trapazoid",
            "Triangle",
        ]
        self.apodization_combobox_dim2 = wx.ComboBox(
            parent, -1, choices=self.apodization_options_dim2, style=wx.CB_READONLY
        )
        self.apodization_combobox_dim2.SetSelection(
            self.apodization_dim2_combobox_selection
        )
        self.apodization_combobox_dim2.Bind(
            wx.EVT_COMBOBOX, self.on_apodization_combobox_dim2
        )
        self.apodization_sizer_dim2.Add(
            self.apodization_combobox_dim2, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.apodization_sizer_dim2.AddSpacer(10)
        if self.apodization_dim2_combobox_selection == 1:
            # Have a textcontrol for the line broadening
            self.apodization_line_broadening_label = wx.StaticText(
                parent, -1, "Line Broadening (Hz):"
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_line_broadening_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_line_broadening_textcontrol_dim2 = wx.TextCtrl(
                parent, -1, str(self.exponential_line_broadening_dim2), size=(30, 20)
            )
            self.apodization_line_broadening_textcontrol_dim2.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_line_broadening_textcontrol_dim2,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(
                parent, -1, "First Point Scaling:"
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_first_point_textcontrol_dim2 = wx.TextCtrl(
                parent,
                -1,
                str(self.apodization_first_point_scaling_dim2),
                size=(30, 20),
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_first_point_textcontrol_dim2,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            self.apodization_sizer_dim2.AddSpacer(10)
        elif self.apodization_dim2_combobox_selection == 2:
            # Have a textcontrol for the g1 value
            self.apodization_g1_label = wx.StaticText(
                parent, -1, "Inverse Lorentzian (Hz):"
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_g1_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_g1_textcontrol_dim2 = wx.TextCtrl(
                parent, -1, str(self.g1_dim2), size=(40, 20)
            )
            self.apodization_g1_textcontrol_dim2.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_g1_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the g2 value
            self.apodization_g2_label = wx.StaticText(
                parent, -1, "Gaussian Broadening (Hz):"
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_g2_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_g2_textcontrol_dim2 = wx.TextCtrl(
                parent, -1, str(self.g2_dim2), size=(40, 20)
            )
            self.apodization_g2_textcontrol_dim2.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_g2_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the g3 value
            self.apodization_g3_label = wx.StaticText(parent, -1, "Gaussian Shift:")
            self.apodization_sizer_dim2.Add(
                self.apodization_g3_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_g3_textcontrol_dim2 = wx.TextCtrl(
                parent, -1, str(self.g3_dim2), size=(40, 20)
            )
            self.apodization_g3_textcontrol_dim2.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_g3_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(
                parent, -1, "First Point Scaling:"
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_first_point_textcontrol_dim2 = wx.TextCtrl(
                self, -1, str(self.apodization_first_point_scaling_dim2), size=(30, 20)
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_first_point_textcontrol_dim2,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            self.apodization_sizer_dim2.AddSpacer(10)
        elif self.apodization_dim2_combobox_selection == 3:
            # Have a textcontrol for the offset value
            self.apodization_offset_label = wx.StaticText(
                parent, -1, "Offset (\u03c0):"
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_offset_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_offset_textcontrol_dim2 = wx.TextCtrl(
                parent, -1, str(self.offset_dim2), size=(40, 20)
            )
            self.apodization_offset_textcontrol_dim2.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_offset_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the end value
            self.apodization_end_label = wx.StaticText(parent, -1, "End (\u03c0):")
            self.apodization_sizer_dim2.Add(
                self.apodization_end_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_end_textcontrol_dim2 = wx.TextCtrl(
                parent, -1, str(self.end_dim2), size=(40, 20)
            )
            self.apodization_end_textcontrol_dim2.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_end_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the power value
            self.apodization_power_label = wx.StaticText(parent, -1, "Power:")
            self.apodization_sizer_dim2.Add(
                self.apodization_power_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_power_textcontrol_dim2 = wx.TextCtrl(
                parent, -1, str(self.power_dim2), size=(30, 20)
            )
            self.apodization_power_textcontrol_dim2.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_power_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(
                parent, -1, "First Point Scaling:"
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_first_point_textcontrol_dim2 = wx.TextCtrl(
                parent,
                -1,
                str(self.apodization_first_point_scaling_dim2),
                size=(30, 20),
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_first_point_textcontrol_dim2,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            self.apodization_sizer_dim2.AddSpacer(10)
        elif self.apodization_dim2_combobox_selection == 4:
            # Have a textcontrol for the a value
            self.apodization_a_label = wx.StaticText(
                parent, -1, "Line Broadening (Hz):"
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_a_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_a_textcontrol_dim2 = wx.TextCtrl(
                parent, -1, str(self.a_dim2), size=(40, 20)
            )
            self.apodization_a_textcontrol_dim2.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_a_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the b value
            self.apodization_b_label = wx.StaticText(
                parent, -1, "Gaussian Broadening (Hz):"
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_b_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_b_textcontrol_dim2 = wx.TextCtrl(
                parent, -1, str(self.b_dim2), size=(40, 20)
            )
            self.apodization_b_textcontrol_dim2.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_b_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(
                parent, -1, "First Point Scaling:"
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_first_point_textcontrol_dim2 = wx.TextCtrl(
                self, -1, str(self.apodization_first_point_scaling_dim2), size=(30, 20)
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_first_point_textcontrol_dim2,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            self.apodization_sizer_dim2.AddSpacer(10)
        elif self.apodization_dim2_combobox_selection == 5:
            # Have a textcontrol for the t1 value
            self.apodization_t1_label = wx.StaticText(parent, -1, "Ramp up points:")
            self.apodization_sizer_dim2.Add(
                self.apodization_t1_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_t1_textcontrol_dim2 = wx.TextCtrl(
                parent, -1, str(self.t1_dim2), size=(50, 20)
            )
            self.apodization_t1_textcontrol_dim2.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_t1_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the t2 value
            self.apodization_t2_label = wx.StaticText(parent, -1, "Ramp down points:")
            self.apodization_sizer_dim2.Add(
                self.apodization_t2_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_t2_textcontrol_dim2 = wx.TextCtrl(
                parent, -1, str(self.t2_dim2), size=(50, 20)
            )
            self.apodization_t2_textcontrol_dim2.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_t2_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(
                parent, -1, "First Point Scaling:"
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_first_point_textcontrol_dim2 = wx.TextCtrl(
                parent,
                -1,
                str(self.apodization_first_point_scaling_dim2),
                size=(30, 20),
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_first_point_textcontrol_dim2,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            self.apodization_sizer_dim2.AddSpacer(10)
        elif self.apodization_dim2_combobox_selection == 6:
            # Have a textcontrol for the loc value
            self.apodization_loc_label = wx.StaticText(
                parent, -1, "Location of maximum:"
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_loc_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_loc_textcontrol_dim2 = wx.TextCtrl(
                parent, -1, str(self.loc_dim2), size=(40, 20)
            )
            self.apodization_loc_textcontrol_dim2.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_loc_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(
                parent, -1, "First Point Scaling:"
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_first_point_textcontrol_dim2 = wx.TextCtrl(
                parent,
                -1,
                str(self.apodization_first_point_scaling_dim2),
                size=(30, 20),
            )
            self.apodization_sizer_dim2.Add(
                self.apodization_first_point_textcontrol_dim2,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            self.apodization_sizer_dim2.AddSpacer(10)

        # Have a button for information on currently selected apodization containing unicode i in a circle
        self.apodization_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))
        self.apodization_info.Bind(wx.EVT_BUTTON, self.oneDFrame.on_apodization_info)
        self.apodization_sizer_dim2.Add(
            self.apodization_info, 0, wx.ALIGN_CENTER_VERTICAL
        )

        # Have a mini plots of the apodization function along with the FID first slice
        self.plot_window_function_dim2()

        self.sizer_2.Add(self.apodization_sizer_dim2)
        self.sizer_2.AddSpacer(10)

    def on_apodization_checkbox_dim2(self, event):
        # Get the selection from the checkbox
        self.apodization_dim2_checkbox_selection = (
            self.apodization_checkbox_dim2.GetValue()
        )

    def on_apodization_textcontrol_dim2(self, event):
        # If the user presses enter, update the plot
        keycode = event.GetKeyCode()
        if keycode == wx.WXK_RETURN:
            self.update_window_function_plot_dim2()
        event.Skip()

    def plot_window_function_dim2(self):
        self.apodization_plot_sizer_dim2 = wx.BoxSizer(wx.VERTICAL)
        self.apodization_plot_dim2 = Figure(figsize=(1, 0.5), facecolor="#e6e6e7")
        self.apodization_plot_ax_dim2 = self.apodization_plot_dim2.add_subplot(111)
        # self.apodization_plot_ax.set_axis_off()

        self.apodization_plot_ax_dim2.set_xticks([])
        self.apodization_plot_ax_dim2.set_yticks([])

        # If the apodization function is None, make remove the axes of the plot
        if self.apodization_dim2_combobox_selection == 0:
            self.apodization_plot_ax_dim2.spines["top"].set_visible(False)
            self.apodization_plot_ax_dim2.spines["right"].set_visible(False)
            self.apodization_plot_ax_dim2.spines["bottom"].set_visible(False)
            self.apodization_plot_ax_dim2.spines["left"].set_visible(False)

        # If have a pseudo axis and the pseudo axis is the 2nd dimension
        if self.nmr_data.pseudo_axis == True and self.nmr_data.index == 1:
            x = np.linspace(
                0,
                (self.nmr_data.number_of_points[2] / 2)
                / self.nmr_data.spectral_width[2],
                int(self.nmr_data.number_of_points[2] / 2),
            )
            if self.apodization_dim2_combobox_selection == 1:
                # Exponential window function
                (self.line1,) = self.apodization_plot_ax_dim2.plot(
                    x,
                    np.exp(-(np.pi * x * self.exponential_line_broadening_dim2)),
                    color="#1f77b4",
                )
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(
                    0,
                    (self.nmr_data.number_of_points[2] / 2)
                    / self.nmr_data.spectral_width[2],
                )

            elif self.apodization_dim2_combobox_selection == 2:
                # Lorentz to Gauss window function
                e = (
                    np.pi
                    * (self.nmr_data.number_of_points[2] / 2)
                    / self.nmr_data.spectral_width[2]
                    * self.g1_dim2
                )
                g = (
                    0.6
                    * np.pi
                    * self.g2_dim2
                    * (
                        self.g3_dim2
                        * (
                            (self.nmr_data.number_of_points[2] / 2)
                            / self.nmr_data.spectral_width[2]
                            - 1
                        )
                        - x
                    )
                )
                func = np.exp(e - g * g)
                (self.line1,) = self.apodization_plot_ax_dim2.plot(
                    x, func, color="#1f77b4"
                )
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(
                    0,
                    (self.nmr_data.number_of_points[2] / 2)
                    / self.nmr_data.spectral_width[2],
                )
            elif self.apodization_dim2_combobox_selection == 3:
                # Sinebell window function
                func = (
                    np.sin(
                        (
                            np.pi * self.offset_dim2
                            + np.pi * (self.end_dim2 - self.offset_dim2) * x
                        )
                        / (
                            (
                                (
                                    (self.nmr_data.number_of_points[2] / 2)
                                    / self.nmr_data.spectral_width[2]
                                )
                            )
                        )
                    )
                    ** self.power_dim2
                )
                (self.line1,) = self.apodization_plot_ax_dim2.plot(
                    x, func, color="#1f77b4"
                )
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(
                    0,
                    (self.nmr_data.number_of_points[2] / 2)
                    / self.nmr_data.spectral_width[2],
                )
            elif self.apodization_dim2_combobox_selection == 4:
                # Gauss broadening window function
                func = np.exp(-self.a_dim2 * (x**2) - self.b_dim2 * x)
                (self.line1,) = self.apodization_plot_ax_dim2.plot(
                    x, func, color="#1f77b4"
                )
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(
                    0,
                    self.nmr_data.number_of_points[2] / self.nmr_data.spectral_width[2],
                )
            elif self.apodization_dim2_combobox_selection == 5:
                # Trapazoid window function
                func = np.concatenate(
                    (
                        np.linspace(0, 1, int(self.t1_dim2)),
                        np.ones(
                            int(self.nmr_data.number_of_points[2] / 2)
                            - int(self.t1_dim2)
                            - int(self.t1_dim2)
                        ),
                        np.linspace(1, 0, int(self.t2_dim2)),
                    )
                )
                (self.line1,) = self.apodization_plot_ax_dim2.plot(
                    x, func, color="#1f77b4"
                )
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(
                    0,
                    self.nmr_data.number_of_points[2] / self.nmr_data.spectral_width[2],
                )
            elif self.apodization_dim2_combobox_selection == 6:
                # Triangle window function
                func = np.concatenate(
                    (
                        np.linspace(
                            0,
                            1,
                            int(
                                self.loc_dim2 * (self.nmr_data.number_of_points[2] / 2)
                            ),
                        ),
                        np.linspace(
                            1,
                            0,
                            int(
                                (1 - self.loc_dim2)
                                * (self.nmr_data.number_of_points[2] / 2)
                            ),
                        ),
                    )
                )
                (self.line1,) = self.apodization_plot_ax_dim2.plot(
                    x, func, color="#1f77b4"
                )

                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(
                    0,
                    (self.nmr_data.number_of_points[2] / 2)
                    / self.nmr_data.spectral_width[2],
                )

            self.apodization_plot_ax_dim2.set_xlim(
                0, self.nmr_data.number_of_points[2] / self.nmr_data.spectral_width[2]
            )

        else:
            x = np.linspace(
                0,
                (self.nmr_data.number_of_points[1] / 2)
                / self.nmr_data.spectral_width[1],
                int(self.nmr_data.number_of_points[1] / 2),
            )
            if self.apodization_dim2_combobox_selection == 1:
                # Exponential window function
                (self.line1,) = self.apodization_plot_ax_dim2.plot(
                    x,
                    np.exp(-(np.pi * x * self.exponential_line_broadening_dim2)),
                    color="#1f77b4",
                )
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(
                    0,
                    (self.nmr_data.number_of_points[1] / 2)
                    / self.nmr_data.spectral_width[1],
                )

            elif self.apodization_dim2_combobox_selection == 2:
                # Lorentz to Gauss window function
                e = (
                    np.pi
                    * (self.nmr_data.number_of_points[1] / 2)
                    / self.nmr_data.spectral_width[1]
                    * self.g1_dim2
                )
                g = (
                    0.6
                    * np.pi
                    * self.g2_dim2
                    * (
                        self.g3_dim2
                        * (
                            (self.nmr_data.number_of_points[1] / 2)
                            / self.nmr_data.spectral_width[1]
                            - 1
                        )
                        - x
                    )
                )
                func = np.exp(e - g * g)
                (self.line1,) = self.apodization_plot_ax_dim2.plot(
                    x, func, color="#1f77b4"
                )
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(
                    0,
                    (self.nmr_data.number_of_points[1] / 2)
                    / self.nmr_data.spectral_width[1],
                )
            elif self.apodization_dim2_combobox_selection == 3:
                # Sinebell window function
                func = (
                    np.sin(
                        (
                            np.pi * self.offset_dim2
                            + np.pi * (self.end_dim2 - self.offset_dim2) * x
                        )
                        / (
                            (
                                (
                                    (self.nmr_data.number_of_points[1] / 2)
                                    / self.nmr_data.spectral_width[1]
                                )
                            )
                        )
                    )
                    ** self.power_dim2
                )
                (self.line1,) = self.apodization_plot_ax_dim2.plot(
                    x, func, color="#1f77b4"
                )
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(
                    0,
                    (self.nmr_data.number_of_points[1] / 2)
                    / self.nmr_data.spectral_width[1],
                )
            elif self.apodization_dim2_combobox_selection == 4:
                # Gauss broadening window function
                func = np.exp(-self.a_dim2 * (x**2) - self.b_dim2 * x)
                (self.line1,) = self.apodization_plot_ax_dim2.plot(
                    x, func, color="#1f77b4"
                )
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(
                    0,
                    self.nmr_data.number_of_points[1] / self.nmr_data.spectral_width[1],
                )
            elif self.apodization_dim2_combobox_selection == 5:
                # Trapazoid window function
                func = np.concatenate(
                    (
                        np.linspace(0, 1, int(self.t1_dim2)),
                        np.ones(
                            int(self.nmr_data.number_of_points[1] / 2)
                            - int(self.t1_dim2)
                            - int(self.t1_dim2)
                        ),
                        np.linspace(1, 0, int(self.t2_dim2)),
                    )
                )
                (self.line1,) = self.apodization_plot_ax_dim2.plot(
                    x, func, color="#1f77b4"
                )
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(
                    0,
                    self.nmr_data.number_of_points[1] / self.nmr_data.spectral_width[1],
                )
            elif self.apodization_dim2_combobox_selection == 6:
                # Triangle window function
                func = np.concatenate(
                    (
                        np.linspace(
                            0,
                            1,
                            int(
                                self.loc_dim2 * (self.nmr_data.number_of_points[1] / 2)
                            ),
                        ),
                        np.linspace(
                            1,
                            0,
                            int(
                                (1 - self.loc_dim2)
                                * (self.nmr_data.number_of_points[1] / 2)
                            ),
                        ),
                    )
                )
                (self.line1,) = self.apodization_plot_ax_dim2.plot(
                    x, func, color="#1f77b4"
                )

                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(
                    0,
                    (self.nmr_data.number_of_points[1] / 2)
                    / self.nmr_data.spectral_width[1],
                )

            self.apodization_plot_ax_dim2.set_xlim(
                0, self.nmr_data.number_of_points[1] / self.nmr_data.spectral_width[1]
            )

        self.apodization_plot_canvas = FigCanvas(self, -1, self.apodization_plot_dim2)
        self.apodization_plot_sizer_dim2.Add(self.apodization_plot_canvas, 0, wx.EXPAND)

        self.apodization_sizer_dim2.Add(
            self.apodization_plot_sizer_dim2, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.apodization_sizer_dim2.AddSpacer(10)

    def update_window_function_plot_dim2(self):
        if self.nmr_data.pseudo_axis == True and self.nmr_data.index == 1:
            x = np.linspace(
                0,
                (self.nmr_data.number_of_points[2] / 2)
                / self.nmr_data.spectral_width[2],
                int(self.nmr_data.number_of_points[2] / 2),
            )
        else:
            x = np.linspace(
                0,
                (self.nmr_data.number_of_points[1] / 2)
                / self.nmr_data.spectral_width[1],
                int(self.nmr_data.number_of_points[1] / 2),
            )
        try:
            c = float(self.apodization_first_point_textcontrol_dim2.GetValue())
            self.apodization_first_point_scaling_dim2 = c
        except:
            # Give a popout window saying that the values are not valid
            msg = wx.MessageDialog(
                self,
                "The value entered for apodization first point scaling is not valid (use 0.5 or 1.0)",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            self.apodization_first_point_textcontrol_dim2.SetValue(
                str(self.apodization_first_point_scaling_dim2)
            )
            return
        if c != 0.5 and c != 1.0:
            msg = wx.MessageDialog(
                self,
                "The value entered for apodization first point scaling is not valid (use 0.5 or 1.0)",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            self.apodization_first_point_textcontrol_dim2.SetValue(
                str(self.apodization_first_point_scaling_dim2)
            )
            return
        self.apodization_first_point_scaling = c
        if self.apodization_dim2_combobox_selection == 1:
            try:
                em = float(self.apodization_line_broadening_textcontrol_dim2.GetValue())
            except:
                # Give a popout window saying that the values are not valid
                msg = wx.MessageDialog(
                    self,
                    "The values entered are not valid",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_line_broadening_textcontrol_dim2.SetValue(
                    str(self.exponential_line_broadening_dim2)
                )
                return
            self.exponential_line_broadening_dim2 = em

            self.line1.set_ydata(
                np.exp(-(np.pi * x * self.exponential_line_broadening_dim2))
            )
        elif self.apodization_dim2_combobox_selection == 2:
            try:
                g1 = float(self.apodization_g1_textcontrol_dim2.GetValue())
                g2 = float(self.apodization_g2_textcontrol_dim2.GetValue())
                g3 = float(self.apodization_g3_textcontrol_dim2.GetValue())
            except:
                # Give a popout window saying that the values are not valid
                msg = wx.MessageDialog(
                    self,
                    "The values entered are not valid",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_g1_textcontrol_dim2.SetValue(str(self.g1_dim2))
                self.apodization_g2_textcontrol_dim2.SetValue(str(self.g2_dim2))
                self.apodization_g3_textcontrol_dim2.SetValue(str(self.g3_dim2))
                return
            # Check to see if g3 is between 0 and 1
            if g3 < 0 or g3 > 1:
                msg = wx.MessageDialog(
                    self,
                    "Gaussian shift must be between 0 and 1",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_g3_textcontrol_dim2.SetValue(str(self.g3_dim2))
                return
            self.g1_dim2 = g1
            self.g2_dim2 = g2
            self.g3_dim2 = g3
            if self.nmr_data.pseudo_axis == True and self.nmr_data.index == 1:
                e = (
                    np.pi
                    * (self.nmr_data.number_of_points[2] / 2)
                    / self.nmr_data.spectral_width[2]
                    * self.g1_dim2
                )
                g = (
                    0.6
                    * np.pi
                    * self.g2_dim2
                    * (
                        self.g3_dim2
                        * (
                            (self.nmr_data.number_of_points[2] / 2)
                            / self.nmr_data.spectral_width[2]
                            - 1
                        )
                        - x
                    )
                )
                func = np.exp(e - g * g)
                self.line1.set_ydata(func)

                self.apodization_plot_ax_dim2.set_xlim(
                    0,
                    (self.nmr_data.number_of_points[1] / 2)
                    / self.nmr_data.spectral_width[1],
                )

            else:
                e = (
                    np.pi
                    * (self.nmr_data.number_of_points[1] / 2)
                    / self.nmr_data.spectral_width[1]
                    * self.g1_dim2
                )
                g = (
                    0.6
                    * np.pi
                    * self.g2_dim2
                    * (
                        self.g3_dim2
                        * (
                            (self.nmr_data.number_of_points[1] / 2)
                            / self.nmr_data.spectral_width[1]
                            - 1
                        )
                        - x
                    )
                )
                func = np.exp(e - g * g)
                self.line1.set_ydata(func)

                self.apodization_plot_ax_dim2.set_xlim(
                    0,
                    (self.nmr_data.number_of_points[1] / 2)
                    / self.nmr_data.spectral_width[1],
                )

        elif self.apodization_dim2_combobox_selection == 3:
            try:
                offset = float(self.apodization_offset_textcontrol_dim2.GetValue())
                end = float(self.apodization_end_textcontrol_dim2.GetValue())
                power = float(self.apodization_power_textcontrol_dim2.GetValue())
                power = int(power)
            except:
                msg = wx.MessageDialog(
                    self,
                    "The values entered are not valid",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_offset_textcontrol_dim2.SetValue(str(self.offset_dim2))
                self.apodization_end_textcontrol_dim2.SetValue(str(self.end_dim2))
                self.apodization_power_textcontrol_dim2.SetValue(str(self.power_dim2))
                return
            # Check that offset and end are between 0 and 1
            if offset < 0 or offset > 1:
                msg = wx.MessageDialog(
                    self,
                    "Offset values must be between 0 and 1",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_offset_textcontrol_dim2.SetValue(str(self.offset_dim2))
                return
            if end < 0 or end > 1:
                msg = wx.MessageDialog(
                    self,
                    "End values must be between 0 and 1",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_end_textcontrol_dim2.SetValue(str(self.end_dim2))
                return
            # Check that power is greater than 0
            if power < 0:
                msg = wx.MessageDialog(
                    self, "Power must be greater than 0", "Error", wx.OK | wx.ICON_ERROR
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_power_textcontrol_dim2.SetValue(str(self.power_dim2))
                return
            self.offset_dim2 = offset
            self.end_dim2 = end
            self.power_dim2 = power
            if self.nmr_data.pseudo_axis == True and self.nmr_data.index == 1:
                func = (
                    np.sin(
                        (
                            np.pi * self.offset_dim2
                            + np.pi * (self.end_dim2 - self.offset_dim2) * x
                        )
                        / (
                            (
                                (
                                    (self.nmr_data.number_of_points[2] / 2)
                                    / self.nmr_data.spectral_width[2]
                                )
                            )
                        )
                    )
                    ** self.power_dim2
                )
            else:
                func = (
                    np.sin(
                        (
                            np.pi * self.offset_dim2
                            + np.pi * (self.end_dim2 - self.offset_dim2) * x
                        )
                        / (
                            (
                                (
                                    (self.nmr_data.number_of_points[1] / 2)
                                    / self.nmr_data.spectral_width[1]
                                )
                            )
                        )
                    )
                    ** self.power_dim2
                )
            self.line1.set_ydata(func)
        elif self.apodization_dim2_combobox_selection == 4:
            try:
                a = float(self.apodization_a_textcontrol_dim2.GetValue())
                b = float(self.apodization_b_textcontrol_dim2.GetValue())
            except:
                msg = wx.MessageDialog(
                    self,
                    "The values entered are not valid",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_a_textcontrol_dim2.SetValue(str(self.a_dim2))
                self.apodization_b_textcontrol_dim2.SetValue(str(self.b_dim2))
                return
            self.a_dim2 = a
            self.b_dim2 = b
            func = np.exp(-self.a_dim2 * (x**2) - self.b_dim2 * x)
            self.line1.set_ydata(func)
        elif self.apodization_dim2_combobox_selection == 5:
            try:
                t1 = float(self.apodization_t1_textcontrol_dim2.GetValue())
                t2 = float(self.apodization_t2_textcontrol_dim2.GetValue())
            except:
                msg = wx.MessageDialog(
                    self,
                    "The values entered are not valid",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_t1_textcontrol_dim2.SetValue(str(self.t1_dim2))
                self.apodization_t2_textcontrol_dim2.SetValue(str(self.t2_dim2))
                return
            # Ensure that t1 and t2 are greater than 0
            if t1 < 0 or t2 < 0:
                msg = wx.MessageDialog(
                    self,
                    "Ramp up and ramp down points must be greater than 0",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_t1_textcontrol_dim2.SetValue(str(self.t1_dim2))
                self.apodization_t2_textcontrol_dim2.SetValue(str(self.t2_dim2))
                return
            # Ensure that t1 + t2 is less than the number of points
            if self.nmr_data.pseudo_axis == True and self.nmr_data.index == 1:
                if t1 + t2 > self.nmr_data.number_of_points[2]:
                    message = (
                        "Ramp up and ramp down points must be less than the number of points ("
                        + str(self.nmr_data.number_of_points[2])
                        + ")"
                    )
                    msg = wx.MessageDialog(
                        self, message, "Error", wx.OK | wx.ICON_ERROR
                    )
                    msg.ShowModal()
                    msg.Destroy()
                    self.apodization_t1_textcontrol_dim2.SetValue(str(self.t1_dim2))
                    self.apodization_t2_textcontrol_dim2.SetValue(str(self.t2_dim2))
                    return

            else:
                if t1 + t2 > self.nmr_data.number_of_points[1]:
                    message = (
                        "Ramp up and ramp down points must be less than the number of points ("
                        + str(self.nmr_data.number_of_points[1])
                        + ")"
                    )
                    msg = wx.MessageDialog(
                        self, message, "Error", wx.OK | wx.ICON_ERROR
                    )
                    msg.ShowModal()
                    msg.Destroy()
                    self.apodization_t1_textcontrol_dim2.SetValue(str(self.t1_dim2))
                    self.apodization_t2_textcontrol_dim2.SetValue(str(self.t2_dim2))
                    return
            self.t1_dim2 = t1
            self.t2_dim2 = t2
            if self.nmr_data.pseudo_axis == True and self.nmr_data.index == 1:
                func = np.concatenate(
                    (
                        np.linspace(0, 1, int(self.t1_dim2)),
                        np.ones(
                            int(self.nmr_data.number_of_points[2] / 2)
                            - int(self.t1_dim2)
                            - int(self.t2_dim2)
                        ),
                        np.linspace(1, 0, int(self.t2_dim2)),
                    )
                )
            else:
                func = np.concatenate(
                    (
                        np.linspace(0, 1, int(self.t1_dim2)),
                        np.ones(
                            int(self.nmr_data.number_of_points[1] / 2)
                            - int(self.t1_dim2)
                            - int(self.t2_dim2)
                        ),
                        np.linspace(1, 0, int(self.t2_dim2)),
                    )
                )
            self.line1.set_ydata(func)
        elif self.apodization_dim2_combobox_selection == 6:
            try:
                loc = float(self.apodization_loc_textcontrol_dim2.GetValue())
            except:
                msg = wx.MessageDialog(
                    self,
                    "The values entered are not valid",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_loc_textcontrol_dim2.SetValue(str(self.loc_dim2))
                return
            # Ensure that loc is between 0 and 1
            if loc < 0 or loc > 1:
                msg = wx.MessageDialog(
                    self,
                    "Location of maximum must be between 0 and 1",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_loc_textcontrol_dim2.SetValue(str(self.loc_dim2))
                return
            self.loc_dim2 = loc
            if self.nmr_data.pseudo_axis == True and self.nmr_data.index == 1:
                func = np.concatenate(
                    (
                        np.linspace(
                            0,
                            1,
                            int(
                                self.loc_dim2 * (self.nmr_data.number_of_points[2] / 2)
                            ),
                        ),
                        np.linspace(
                            1,
                            0,
                            int(self.nmr_data.number_of_points[2] / 2)
                            - int(self.loc_dim2 * self.nmr_data.number_of_points[2]),
                        ),
                    )
                )
            else:
                func = np.concatenate(
                    (
                        np.linspace(
                            0,
                            1,
                            int(
                                self.loc_dim2 * (self.nmr_data.number_of_points[1] / 2)
                            ),
                        ),
                        np.linspace(
                            1,
                            0,
                            int(self.nmr_data.number_of_points[1] / 2)
                            - int(self.loc_dim2 * self.nmr_data.number_of_points[1]),
                        ),
                    )
                )
            self.line1.set_ydata(func)

        self.apodization_plot_canvas.draw()

    def on_apodization_combobox_dim2(self, event):
        self.apodization_dim2_combobox_selection = (
            self.apodization_combobox_dim2.GetSelection()
        )

        # Destroy the combobox and textcontrols for the previous apodization function
        # self.apodization_sizer.Detach(self.apodization_combobox)
        # self.apodization_combobox.Destroy()

        # # Remove the zf sizer
        self.zero_filling_sizer_dim2.Clear(delete_windows=True)

        if self.apodization_dim2_combobox_selection_old == 1:
            # Remove the previous textcontrols

            self.apodization_sizer_dim2.Detach(self.apodization_line_broadening_label)
            self.apodization_line_broadening_label.Destroy()
            self.apodization_sizer_dim2.Detach(
                self.apodization_line_broadening_textcontrol_dim2
            )
            self.apodization_line_broadening_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(
                self.apodization_first_point_textcontrol_dim2
            )
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)

        elif self.apodization_dim2_combobox_selection_old == 2:
            self.apodization_sizer_dim2.Detach(self.apodization_g1_label)
            self.apodization_g1_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g1_textcontrol_dim2)
            self.apodization_g1_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g2_label)
            self.apodization_g2_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g2_textcontrol_dim2)
            self.apodization_g2_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g3_label)
            self.apodization_g3_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g3_textcontrol_dim2)
            self.apodization_g3_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(
                self.apodization_first_point_textcontrol_dim2
            )
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)

        elif self.apodization_dim2_combobox_selection_old == 3:
            self.apodization_sizer_dim2.Detach(self.apodization_offset_label)
            self.apodization_offset_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_offset_textcontrol_dim2)
            self.apodization_offset_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_end_label)
            self.apodization_end_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_end_textcontrol_dim2)
            self.apodization_end_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_power_label)
            self.apodization_power_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_power_textcontrol_dim2)
            self.apodization_power_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(
                self.apodization_first_point_textcontrol_dim2
            )
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        elif self.apodization_dim2_combobox_selection_old == 4:
            self.apodization_sizer_dim2.Detach(self.apodization_a_label)
            self.apodization_a_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_a_textcontrol_dim2)
            self.apodization_a_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_b_label)
            self.apodization_b_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_b_textcontrol_dim2)
            self.apodization_b_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(
                self.apodization_first_point_textcontrol_dim2
            )
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        elif self.apodization_dim2_combobox_selection_old == 5:
            self.apodization_sizer_dim2.Detach(self.apodization_t1_label)
            self.apodization_t1_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_t1_textcontrol_dim2)
            self.apodization_t1_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_t2_label)
            self.apodization_t2_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_t2_textcontrol_dim2)
            self.apodization_t2_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(
                self.apodization_first_point_textcontrol_dim2
            )
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        elif self.apodization_dim2_combobox_selection_old == 6:
            self.apodization_sizer_dim2.Detach(self.apodization_loc_label)
            self.apodization_loc_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_loc_textcontrol_dim2)
            self.apodization_loc_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(
                self.apodization_first_point_textcontrol_dim2
            )
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        elif self.apodization_dim2_combobox_selection_old == 0:
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)

        self.apodization_sizer_dim2.Detach(self.apodization_checkbox_dim2)
        self.apodization_checkbox_dim2.Destroy()
        self.apodization_sizer_dim2.Detach(self.apodization_combobox_dim2)
        self.apodization_combobox_dim2.Hide()

        # Delete the current apodization sizer and then create a new one

        self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
        self.apodization_plot_ax_dim2.clear()
        self.apodization_plot_ax_dim2.clear()
        self.apodization_plot_sizer_dim2.Clear(True)

        self.sizer_2.Remove(self.apodization_sizer_dim2)
        # self.apodization_sizer.Clear(delete_windows=True)

        # Remove the linear prediction sizers
        self.linear_prediction_sizer_dim2.Clear(delete_windows=True)
        # self.sizer_1.Remove(self.linear_prediction_sizer)

        # self.sizer_1.Remove(self.solvent_suppression_sizer)

        self.sizer_2.Clear(delete_windows=True)

        self.create_menu_bar_dim2()
        self.Refresh()
        self.Update()
        self.Layout()

        self.apodization_dim2_combobox_selection_old = (
            self.apodization_dim2_combobox_selection
        )

    def create_zero_filling_sizer_dim2(self, parent):
        # Create a box for zero filling options
        self.zero_filling_box_dim2 = wx.StaticBox(parent, -1, "Zero Filling")
        self.zero_filling_sizer_dim2 = wx.StaticBoxSizer(
            self.zero_filling_box_dim2, wx.HORIZONTAL
        )
        self.zero_filling_checkbox_dim2 = wx.CheckBox(parent, -1, "Apply zero filling")
        self.zero_filling_checkbox_dim2.SetValue(self.zero_filling_dim2_checkbox_value)
        self.zero_filling_checkbox_dim2.Bind(
            wx.EVT_CHECKBOX, self.on_zero_filling_checkbox_dim2
        )
        self.zero_filling_sizer_dim2.Add(
            self.zero_filling_checkbox_dim2, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.zero_filling_sizer_dim2.AddSpacer(10)
        # Have a combobox for zero filling options
        self.zf_options_label = wx.StaticText(parent, -1, "Options:")
        self.zero_filling_sizer_dim2.Add(
            self.zf_options_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.zero_filling_sizer_dim2.AddSpacer(5)
        self.zero_filling_options_dim2 = [
            "Doubling spectrum size",
            "Adding Zeros",
            "Final data size",
        ]
        self.zero_filling_combobox_dim2 = wx.ComboBox(
            parent, -1, choices=self.zero_filling_options_dim2, style=wx.CB_READONLY
        )
        self.zero_filling_combobox_dim2.Bind(
            wx.EVT_COMBOBOX, self.on_zero_filling_combobox_dim2
        )
        self.zero_filling_combobox_dim2.SetSelection(
            self.zero_filling_dim2_combobox_selection
        )
        self.zero_filling_sizer_dim2.Add(
            self.zero_filling_combobox_dim2, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.zero_filling_sizer_dim2.AddSpacer(10)
        if self.zero_filling_dim2_combobox_selection == 0:
            # Have a textcontrol for the doubling number/number of zeros/final data size
            self.zf_value_label = wx.StaticText(parent, -1, "Doubling number:")
            self.zero_filling_sizer_dim2.Add(
                self.zf_value_label, 0, wx.ALIGN_CENTER_VERTICAL
            )

            self.zero_filling_textcontrol_dim2 = wx.TextCtrl(
                parent,
                -1,
                str(self.zero_filling_dim2_value_doubling_times),
                size=(40, 20),
            )
            self.zero_filling_sizer_dim2.AddSpacer(5)
            self.zero_filling_sizer_dim2.Add(
                self.zero_filling_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.zero_filling_textcontrol_dim2.Bind(
                wx.EVT_TEXT, self.on_zero_filling_doubling_number_dim2
            )
            self.zero_filling_sizer_dim2.AddSpacer(20)
        elif self.zero_filling_dim2_combobox_selection == 1:
            # Have a textcontrol for the doubling number/number of zeros/final data size
            self.zf_value_label = wx.StaticText(parent, -1, "Number of zeros to add:")
            self.zero_filling_sizer_dim2.Add(
                self.zf_value_label, 0, wx.ALIGN_CENTER_VERTICAL
            )

            self.zero_filling_textcontrol_dim2 = wx.TextCtrl(
                parent,
                -1,
                str(self.zero_filling_dim2_value_zeros_to_add),
                size=(40, 20),
            )
            self.zero_filling_sizer_dim2.AddSpacer(5)
            self.zero_filling_sizer_dim2.Add(
                self.zero_filling_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.zero_filling_textcontrol_dim2.Bind(
                wx.EVT_TEXT, self.on_zero_filling_zeros_to_add_dim2
            )
            self.zero_filling_sizer_dim2.AddSpacer(20)
        elif self.zero_filling_dim2_combobox_selection == 2:
            # Have a textcontrol for the doubling number/number of zeros/final data size
            self.zf_value_label = wx.StaticText(parent, -1, "Final data size:")
            self.zero_filling_sizer_dim2.Add(
                self.zf_value_label, 0, wx.ALIGN_CENTER_VERTICAL
            )

            self.zero_filling_textcontrol_dim2 = wx.TextCtrl(
                parent,
                -1,
                str(self.zero_filling_dim2_value_final_data_size),
                size=(40, 20),
            )
            self.zero_filling_textcontrol_dim2.Bind(
                wx.EVT_TEXT, self.on_zero_filling_final_size_dim2
            )
            self.zero_filling_sizer_dim2.AddSpacer(5)
            self.zero_filling_sizer_dim2.Add(
                self.zero_filling_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.zero_filling_sizer_dim2.AddSpacer(20)

        # Have a checkbox for rounding to the nearest power of 2
        self.zero_filling_round_checkbox_dim2 = wx.CheckBox(
            parent, -1, "Round to nearest power of 2"
        )
        self.zero_filling_round_checkbox_dim2.SetValue(
            self.zero_filling_dim2_round_checkbox_value
        )
        self.zero_filling_round_checkbox_dim2.Bind(
            wx.EVT_CHECKBOX, self.on_zero_filling_round_checkbox_dim2
        )
        self.zero_filling_sizer_dim2.Add(
            self.zero_filling_round_checkbox_dim2, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.zero_filling_sizer_dim2.AddSpacer(10)

        # Have a button showing information on zero filling
        self.zero_filling_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))
        self.zero_filling_info.Bind(wx.EVT_BUTTON, self.oneDFrame.on_zero_fill_info)
        self.zero_filling_sizer_dim2.Add(
            self.zero_filling_info, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.zero_filling_sizer_dim2.AddSpacer(10)

        self.sizer_2.Add(self.zero_filling_sizer_dim2)
        self.sizer_2.AddSpacer(10)

    def on_zero_filling_checkbox_dim2(self, event):
        self.zero_filling_dim2_checkbox_value = (
            self.zero_filling_checkbox_dim2.GetValue()
        )

    def on_zero_filling_round_checkbox_dim2(self, event):
        self.zero_filling_dim2_round_checkbox_value = (
            self.zero_filling_round_checkbox_dim2.GetValue()
        )

    def on_zero_filling_final_size_dim2(self, event):
        self.zero_filling_dim2_value_final_data_size = (
            self.zero_filling_textcontrol_dim2.GetValue()
        )

    def on_zero_filling_zeros_to_add_dim2(self, event):
        self.zero_filling_dim2_value_zeros_to_add = (
            self.zero_filling_textcontrol_dim2.GetValue()
        )

    def on_zero_filling_doubling_number_dim2(self, event):
        self.zero_filling_dim2_value_doubling_times = (
            self.zero_filling_textcontrol_dim2.GetValue()
        )

    def on_zero_filling_combobox_dim2(self, event):
        self.zero_filling_dim2_combobox_selection = (
            self.zero_filling_combobox_dim2.GetSelection()
        )
        # # # Remove the zf sizer
        self.zero_filling_sizer_dim2.Clear()
        self.zero_filling_sizer_dim2.Detach(self.zero_filling_checkbox_dim2)
        self.zero_filling_checkbox_dim2.Destroy()
        self.zero_filling_sizer_dim2.Detach(self.zf_options_label)
        self.zf_options_label.Destroy()
        self.zero_filling_sizer_dim2.Detach(self.zero_filling_info)
        self.zero_filling_info.Destroy()
        self.zero_filling_sizer_dim2.Detach(self.zf_value_label)
        self.zf_value_label.Destroy()
        self.zero_filling_sizer_dim2.Detach(self.zero_filling_round_checkbox_dim2)
        self.zero_filling_round_checkbox_dim2.Destroy()
        self.zero_filling_sizer_dim2.Detach(self.zero_filling_textcontrol_dim2)
        self.zero_filling_textcontrol_dim2.Destroy()

        self.zero_filling_sizer_dim2.Detach(self.zero_filling_combobox_dim2)
        self.zero_filling_combobox_dim2.Hide()

        if self.apodization_dim2_combobox_selection_old == 1:
            # Remove the previous textcontrols

            self.apodization_sizer_dim2.Detach(self.apodization_line_broadening_label)
            self.apodization_line_broadening_label.Destroy()
            self.apodization_sizer_dim2.Detach(
                self.apodization_line_broadening_textcontrol_dim2
            )
            self.apodization_line_broadening_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(
                self.apodization_first_point_textcontrol_dim2
            )
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)

        elif self.apodization_dim2_combobox_selection_old == 2:
            self.apodization_sizer_dim2.Detach(self.apodization_g1_label)
            self.apodization_g1_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g1_textcontrol_dim2)
            self.apodization_g1_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g2_label)
            self.apodization_g2_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g2_textcontrol_dim2)
            self.apodization_g2_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g3_label)
            self.apodization_g3_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g3_textcontrol_dim2)
            self.apodization_g3_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(
                self.apodization_first_point_textcontrol_dim2
            )
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)

        elif self.apodization_dim2_combobox_selection_old == 3:
            self.apodization_sizer_dim2.Detach(self.apodization_offset_label)
            self.apodization_offset_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_offset_textcontrol_dim2)
            self.apodization_offset_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_end_label)
            self.apodization_end_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_end_textcontrol_dim2)
            self.apodization_end_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_power_label)
            self.apodization_power_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_power_textcontrol_dim2)
            self.apodization_power_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(
                self.apodization_first_point_textcontrol_dim2
            )
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        elif self.apodization_dim2_combobox_selection_old == 4:
            self.apodization_sizer_dim2.Detach(self.apodization_a_label)
            self.apodization_a_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_a_textcontrol_dim2)
            self.apodization_a_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_b_label)
            self.apodization_b_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_b_textcontrol_dim2)
            self.apodization_b_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(
                self.apodization_first_point_textcontrol_dim2
            )
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        elif self.apodization_dim2_combobox_selection_old == 5:
            self.apodization_sizer_dim2.Detach(self.apodization_t1_label)
            self.apodization_t1_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_t1_textcontrol_dim2)
            self.apodization_t1_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_t2_label)
            self.apodization_t2_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_t2_textcontrol_dim2)
            self.apodization_t2_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(
                self.apodization_first_point_textcontrol_dim2
            )
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        elif self.apodization_dim2_combobox_selection_old == 6:
            self.apodization_sizer_dim2.Detach(self.apodization_loc_label)
            self.apodization_loc_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_loc_textcontrol_dim2)
            self.apodization_loc_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(
                self.apodization_first_point_textcontrol_dim2
            )
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        elif self.apodization_dim2_combobox_selection_old == 0:
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)

        self.apodization_sizer_dim2.Detach(self.apodization_checkbox_dim2)
        self.apodization_checkbox_dim2.Destroy()
        self.apodization_sizer_dim2.Detach(self.apodization_combobox_dim2)
        self.apodization_combobox_dim2.Hide()

        # Delete the current apodization sizer and then create a new one

        self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
        self.apodization_plot_ax_dim2.clear()
        self.apodization_plot_ax_dim2.clear()
        self.apodization_plot_sizer_dim2.Clear(True)

        self.sizer_2.Remove(self.apodization_sizer_dim2)
        # self.apodization_sizer.Clear(delete_windows=True)

        # Remove the linear prediction sizers
        self.linear_prediction_sizer_dim2.Clear(delete_windows=True)
        # self.sizer_1.Remove(self.linear_prediction_sizer)

        # self.sizer_1.Remove(self.solvent_suppression_sizer)

        self.sizer_2.Clear(delete_windows=True)

        self.create_menu_bar_dim2()
        self.Refresh()
        self.Update()
        self.Layout()

    def create_fourier_transform_sizer_dim2(self, parent):
        # Create a box for fourier transform options
        self.fourier_transform_box = wx.StaticBox(parent, -1, "Fourier Transform")
        self.fourier_transform_sizer_dim2 = wx.StaticBoxSizer(
            self.fourier_transform_box, wx.HORIZONTAL
        )
        self.fourier_transform_checkbox_dim2 = wx.CheckBox(
            parent, -1, "Apply fourier transform"
        )
        self.fourier_transform_checkbox_dim2.Bind(
            wx.EVT_CHECKBOX, self.on_fourier_transform_checkbox_dim2
        )
        self.fourier_transform_checkbox_dim2.SetValue(
            self.fourier_transform_dim2_checkbox_value
        )
        self.fourier_transform_sizer_dim2.Add(
            self.fourier_transform_checkbox_dim2, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.fourier_transform_sizer_dim2.AddSpacer(10)
        # Have a button for advanced options for fourier transform
        self.fourier_transform_advanced_options_dim2 = wx.Button(
            parent, -1, "Advanced Options"
        )
        self.fourier_transform_advanced_options_dim2.Bind(
            wx.EVT_BUTTON, self.on_fourier_transform_advanced_options_dim2
        )
        self.fourier_transform_sizer_dim2.Add(
            self.fourier_transform_advanced_options_dim2, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.fourier_transform_sizer_dim2.AddSpacer(10)

        # Have a button showing information on fourier transform
        self.fourier_transform_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))
        self.fourier_transform_info.Bind(
            wx.EVT_BUTTON, self.oneDFrame.on_fourier_transform_info
        )
        self.fourier_transform_sizer_dim2.Add(
            self.fourier_transform_info, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.sizer_2.Add(self.fourier_transform_sizer_dim2)
        self.sizer_2.AddSpacer(10)

    def on_fourier_transform_checkbox_dim2(self, event):
        self.fourier_transform_dim2_checkbox_value = (
            self.fourier_transform_checkbox_dim2.GetValue()
        )

    def on_fourier_transform_advanced_options_dim2(self, event):
        # Create a frame with a set of advanced options for the fourier transform implementation
        self.fourier_transform_advanced_options_window_dim2 = wx.Frame(
            self,
            -1,
            "Fourier Transform Advanced Options (Dimension 2)",
            size=(700, 300),
        )

        self.fourier_transform_advanced_options_window_sizer_dim2 = wx.BoxSizer(
            wx.VERTICAL
        )
        self.fourier_transform_advanced_options_window_dim2.SetSizer(
            self.fourier_transform_advanced_options_window_sizer_dim2
        )

        # Create a sizer for the fourier transform advanced options
        self.ft_label = wx.StaticBox(
            self.fourier_transform_advanced_options_window_dim2,
            -1,
            "Fourier Transform Method:",
        )
        self.fourier_transform_advanced_options_sizer_dim2 = wx.StaticBoxSizer(
            self.ft_label, wx.VERTICAL
        )

        # Have a radiobox for auto, real, inverse, sign alternation
        self.fourier_transform_advanced_options_sizer_dim2.AddSpacer(10)
        self.fourier_transform_auto_real_inverse_sign_alternation_radio_box_dim2 = (
            wx.RadioBox(
                self.fourier_transform_advanced_options_window_dim2,
                -1,
                choices=["Auto", "Real", "Inverse", "Sign Alternation", "Negative"],
                style=wx.RA_SPECIFY_COLS,
            )
        )
        self.fourier_transform_auto_real_inverse_sign_alternation_radio_box_dim2.SetSelection(
            self.ft_method_selection_dim2
        )
        self.fourier_transform_advanced_options_sizer_dim2.Add(
            self.fourier_transform_auto_real_inverse_sign_alternation_radio_box_dim2,
            0,
            wx.ALIGN_CENTER_HORIZONTAL,
        )
        self.fourier_transform_advanced_options_sizer_dim2.AddSpacer(10)

        self.ft_method_text = "Auto: The auto method will automatically select the best method for the fourier transform of the FID. \n\nReal: The Fourier Transform will be applied to the real part of the FID only. \n\nInverse: The inverse Fourier Transform will be applied to the FID. \n\nSign Alternation: The sign alternation method will be applied to the FID. \n\n"

        self.ft_method_info = wx.StaticText(
            self.fourier_transform_advanced_options_window_dim2, -1, self.ft_method_text
        )
        self.fourier_transform_advanced_options_sizer_dim2.Add(
            self.ft_method_info, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.fourier_transform_advanced_options_sizer_dim2.AddSpacer(10)

        # Have a save and close button
        self.fourier_transform_advanced_options_save_button_dim2 = wx.Button(
            self.fourier_transform_advanced_options_window_dim2, -1, "Save and Close"
        )
        self.fourier_transform_advanced_options_save_button_dim2.Bind(
            wx.EVT_BUTTON, self.on_fourier_transform_advanced_options_save_dim2
        )
        self.fourier_transform_advanced_options_sizer_dim2.Add(
            self.fourier_transform_advanced_options_save_button_dim2,
            0,
            wx.ALIGN_CENTER_HORIZONTAL,
        )

        self.fourier_transform_advanced_options_window_sizer_dim2.Add(
            self.fourier_transform_advanced_options_sizer_dim2,
            0,
            wx.ALIGN_CENTER_HORIZONTAL,
        )

        self.fourier_transform_advanced_options_window_dim2.Show()

    def on_fourier_transform_advanced_options_save_dim2(self, event):
        # Save the current selection and close the window
        self.ft_method_selection_dim2 = (
            self.fourier_transform_auto_real_inverse_sign_alternation_radio_box_dim2.GetSelection()
        )
        self.fourier_transform_advanced_options_window_dim2.Close()

    def create_extraction_sizer_dim2(self, parent):
        # A box for extraction of data between two ppm values
        self.extraction_box_dim2 = wx.StaticBox(parent, -1, "Extraction")
        self.extraction_sizer_dim2 = wx.StaticBoxSizer(
            self.extraction_box_dim2, wx.HORIZONTAL
        )
        self.extraction_checkbox_dim2 = wx.CheckBox(
            parent, -1, "Include data extraction"
        )
        self.extraction_checkbox_dim2.Bind(
            wx.EVT_CHECKBOX, self.on_extraction_checkbox_dim2
        )
        self.extraction_checkbox_dim2.SetValue(self.extraction_checkbox_value_dim2)
        self.extraction_sizer_dim2.Add(
            self.extraction_checkbox_dim2, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.extraction_sizer_dim2.AddSpacer(10)
        # Have a textcontrol for the ppm start value
        self.extraction_ppm_start_label = wx.StaticText(
            parent, -1, "Start chemical shift (ppm):"
        )
        self.extraction_sizer_dim2.Add(
            self.extraction_ppm_start_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.extraction_ppm_start_textcontrol_dim2 = wx.TextCtrl(
            parent, -1, str(self.extraction_start_dim2), size=(40, 20)
        )
        self.extraction_ppm_start_textcontrol_dim2.Bind(
            wx.EVT_TEXT, self.on_extraction_dim2
        )
        self.extraction_sizer_dim2.Add(
            self.extraction_ppm_start_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.extraction_sizer_dim2.AddSpacer(10)
        # Have a textcontrol for the ppm end value
        self.extraction_ppm_end_label = wx.StaticText(
            parent, -1, "End chemical shift (ppm):"
        )
        self.extraction_sizer_dim2.Add(
            self.extraction_ppm_end_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.extraction_ppm_end_textcontrol_dim2 = wx.TextCtrl(
            parent, -1, str(self.extraction_end_dim2), size=(40, 20)
        )
        self.extraction_ppm_end_textcontrol_dim2.Bind(
            wx.EVT_TEXT, self.on_extraction_dim2
        )
        self.extraction_sizer_dim2.Add(
            self.extraction_ppm_end_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.extraction_sizer_dim2.AddSpacer(10)
        # Have a button showing information on extraction
        self.extraction_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))
        self.extraction_info.Bind(wx.EVT_BUTTON, self.oneDFrame.on_extraction_info)
        self.extraction_sizer_dim2.Add(
            self.extraction_info, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.sizer_2.Add(self.extraction_sizer_dim2)
        self.sizer_2.AddSpacer(10)

    def on_extraction_checkbox_dim2(self, event):
        self.extraction_checkbox_value_dim2 = self.extraction_checkbox_dim2.GetValue()

    def on_extraction_dim2(self, event):
        self.extraction_start_dim2 = (
            self.extraction_ppm_start_textcontrol_dim2.GetValue()
        )
        self.extraction_end_dim2 = self.extraction_ppm_end_textcontrol_dim2.GetValue()

    def create_baseline_correction_sizer_dim2(self, parent):
        # Create a box for baseline correction options (linear/polynomial)
        self.baseline_correction_box_dim2 = wx.StaticBox(
            parent, -1, "Baseline Correction"
        )
        self.baseline_correction_sizer_dim2 = wx.StaticBoxSizer(
            self.baseline_correction_box_dim2, wx.HORIZONTAL
        )
        self.baseline_correction_checkbox_dim2 = wx.CheckBox(
            parent, -1, "Apply baseline correction"
        )
        self.baseline_correction_checkbox_dim2.Bind(
            wx.EVT_CHECKBOX, self.on_baseline_correction_checkbox_dim2
        )
        self.baseline_correction_checkbox_dim2.SetValue(
            self.baseline_correction_checkbox_value_dim2
        )
        self.baseline_correction_sizer_dim2.Add(
            self.baseline_correction_checkbox_dim2, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.baseline_correction_sizer_dim2.AddSpacer(10)
        # Have a radio box for linear or polynomial baseline correction
        self.baseline_correction_radio_box_dim2 = wx.RadioBox(
            parent, -1, "Baseline Correction Method", choices=["Linear", "Polynomial"]
        )
        # Bind the radio box to a function that will update the baseline correction options
        self.baseline_correction_radio_box_dim2.Bind(
            wx.EVT_RADIOBOX, self.on_baseline_correction_radio_box_dim2
        )
        self.baseline_correction_radio_box_dim2.SetSelection(
            self.baseline_correction_radio_box_selection_dim2
        )
        self.baseline_correction_sizer_dim2.Add(
            self.baseline_correction_radio_box_dim2, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.baseline_correction_sizer_dim2.AddSpacer(10)

        # If linear baseline correction is selected, have a textcontrol for the node values to use
        self.baseline_correction_nodes_label = wx.StaticText(
            parent, -1, "Node width (pts):"
        )
        self.baseline_correction_sizer_dim2.Add(
            self.baseline_correction_nodes_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.baseline_correction_nodes_textcontrol_dim2 = wx.TextCtrl(
            parent, -1, self.node_width_dim2, size=(30, 20)
        )
        self.baseline_correction_nodes_textcontrol_dim2.Bind(
            wx.EVT_TEXT, self.on_baseline_correction_textcontrol_dim2
        )
        self.baseline_correction_sizer_dim2.Add(
            self.baseline_correction_nodes_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.baseline_correction_sizer_dim2.AddSpacer(10)
        # Have a textcontrol for the node list (percentages)
        self.baseline_correction_node_list_label = wx.StaticText(
            parent, -1, "Node list (%):"
        )
        self.baseline_correction_sizer_dim2.Add(
            self.baseline_correction_node_list_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.baseline_correction_node_list_textcontrol_dim2 = wx.TextCtrl(
            parent, -1, self.node_list_dim2, size=(100, 20)
        )
        self.baseline_correction_node_list_textcontrol_dim2.Bind(
            wx.EVT_TEXT, self.on_baseline_correction_textcontrol_dim2
        )
        self.baseline_correction_sizer_dim2.Add(
            self.baseline_correction_node_list_textcontrol_dim2,
            0,
            wx.ALIGN_CENTER_VERTICAL,
        )
        self.baseline_correction_sizer_dim2.AddSpacer(10)
        # If polynomial baseline correction is selected, have a textcontrol for the polynomial order

        self.baseline_correction_polynomial_order_label = wx.StaticText(
            parent, -1, "Polynomial order:"
        )
        self.baseline_correction_sizer_dim2.Add(
            self.baseline_correction_polynomial_order_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.baseline_correction_polynomial_order_textcontrol_dim2 = wx.TextCtrl(
            parent, -1, self.polynomial_order_dim2, size=(30, 20)
        )
        self.baseline_correction_polynomial_order_textcontrol_dim2.Bind(
            wx.EVT_TEXT, self.on_baseline_correction_textcontrol_dim2
        )
        self.baseline_correction_sizer_dim2.Add(
            self.baseline_correction_polynomial_order_textcontrol_dim2,
            0,
            wx.ALIGN_CENTER_VERTICAL,
        )
        self.baseline_correction_sizer_dim2.AddSpacer(10)

        if self.baseline_correction_radio_box_selection_dim2 == 0:
            self.baseline_correction_polynomial_order_label.Hide()
            self.baseline_correction_polynomial_order_textcontrol_dim2.Hide()

        # Have a button showing information on baseline correction
        self.baseline_correction_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))

        self.baseline_correction_info.Bind(
            wx.EVT_BUTTON, self.oneDFrame.on_baseline_correction_info
        )
        self.baseline_correction_sizer_dim2.Add(
            self.baseline_correction_info, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.sizer_2.Add(self.baseline_correction_sizer_dim2)
        self.sizer_2.AddSpacer(10)

    def on_baseline_correction_checkbox_dim2(self, event):
        self.baseline_correction_checkbox_value_dim2 = (
            self.baseline_correction_checkbox_dim2.GetValue()
        )

    def on_baseline_correction_radio_box_dim2(self, event):
        # If the user selects linear or polynomial baseline correction, update the options
        self.baseline_correction_radio_box_selection_dim2 = (
            self.baseline_correction_radio_box_dim2.GetSelection()
        )

        if self.baseline_correction_radio_box_selection_dim2 == 0:
            # Remove the polynomial order textcontrol
            self.baseline_correction_sizer_dim2.Hide(
                self.baseline_correction_polynomial_order_label
            )
            self.baseline_correction_sizer_dim2.Hide(
                self.baseline_correction_polynomial_order_textcontrol_dim2
            )
            self.baseline_correction_sizer_dim2.Layout()
        elif self.baseline_correction_radio_box_selection_dim2 == 1:
            # Add the polynomial order textcontrol
            self.baseline_correction_sizer_dim2.Show(
                self.baseline_correction_polynomial_order_label
            )
            self.baseline_correction_sizer_dim2.Show(
                self.baseline_correction_polynomial_order_textcontrol_dim2
            )
            self.baseline_correction_sizer_dim2.Layout()

    def on_baseline_correction_textcontrol_dim2(self, event):
        # If the node width or node list textcontrols are changed, update the node width and node list
        self.node_width_dim2 = (
            self.baseline_correction_nodes_textcontrol_dim2.GetValue()
        )
        self.node_list_dim2 = (
            self.baseline_correction_node_list_textcontrol_dim2.GetValue()
        )
        self.polynomial_order_dim2 = (
            self.baseline_correction_polynomial_order_textcontrol_dim2.GetValue()
        )


class ThreeDFrame(wx.Panel):
    def __init__(self, parent, oneDFrame):
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 0.7 * self.monitorWidth
        self.height = 0.75 * self.monitorHeight
        self.parent = parent
        wx.Panel.__init__(self, parent, id=wx.ID_ANY, size=(self.width, self.height))

        self.oneDFrame = oneDFrame
        # Create panel for processing dimension 1 of the data
        self.nmr_data = parent.nmr_data
        self.set_variables_dim3()
        self.create_canvas_dim3()
        self.create_menu_bar_dim3()

    def set_variables_dim3(self):

        self.set_initial_linear_prediction_variables_dim3()
        self.set_initial_apodization_variables_dim3()
        self.set_initial_zero_filling_variables_dim3()
        self.set_initial_fourier_transform_variables_dim3()
        self.set_initial_phasing_variables_dim3()
        self.set_initial_extraction_variables_dim3()
        self.set_initial_baseline_correction_variables_dim3()

        if self.parent.load_variables == True:
            try:
                self.load_variables_from_nmrproc_com_3D()
            except:
                pass

    def load_variables_from_nmrproc_com_3D(self):
        # Open processing_parameters.txt file and load the variables from it
        file = open("processing_parameters.txt", "r")
        lines = file.readlines()
        file.close()

        include_line = False
        for line in lines:
            if "Dimension 3" in line:
                include_line = True
                continue
            if include_line == False:
                continue
            if include_line == True and "Dimension 4" in line:
                include_line = False
                break
            if include_line == True:
                line = line.split("\n")[0]
                if line.split(":")[0] == "Linear Prediction":
                    self.linear_prediction_radio_box_dim3_selection = int(
                        line.split(": ")[1]
                    )
                if line.split(":")[0] == "Linear Prediction Options Selection":
                    self.linear_prediction_dim3_options_selection = int(
                        line.split(": ")[1]
                    )
                if line.split(":")[0] == "Linear Prediction Coefficients Selection":
                    self.linear_prediction_dim3_coefficients_selection = int(
                        line.split(": ")[1]
                    )
                if line.split(":")[0] == "NUS file":
                    self.nuslist_name_dim3 = line.split(": ")[1]
                if line.split(":")[0] == "NUS CPU":
                    self.smile_nus_cpu_textcontrol_dim3 = int(line.split(": ")[1])
                if line.split(":")[0] == "NUS Iterations":
                    self.smile_nus_iterations_textcontrol_dim3 = int(
                        line.split(": ")[1]
                    )
                if line.split(":")[0] == "Apodization":
                    if "True" in line:
                        self.apodization_dim3_checkbox_value = True
                    else:
                        self.apodization_dim3_checkbox_value = False
                if line.split(":")[0] == "Apodization Combobox Selection":
                    self.apodization_dim3_combobox_selection = int(line.split(": ")[1])
                    self.apodization_dim3_combobox_selection_old = int(
                        line.split(": ")[1]
                    )
                if line.split(":")[0] == "Exponential Line Broadening":
                    self.exponential_line_broadening_dim3 = float(line.split(": ")[1])
                if line.split(":")[0] == "Apodization First Point Scaling":
                    self.apodization_first_point_scaling_dim3 = float(
                        line.split(": ")[1]
                    )
                if line.split(":")[0] == "G1":
                    self.g1_dim3 = float(line.split(": ")[1])
                if line.split(":")[0] == "G2":
                    self.g2_dim3 = float(line.split(": ")[1])
                if line.split(":")[0] == "G3":
                    self.g3_dim3 = float(line.split(": ")[1])
                if line.split(":")[0] == "Offset":
                    self.offset_dim3 = float(line.split(": ")[1])
                if line.split(":")[0] == "End":
                    self.end_dim3 = float(line.split(": ")[1])
                if line.split(":")[0] == "Power":
                    self.power_dim3 = float(line.split(": ")[1])
                if line.split(":")[0] == "A":
                    self.a_dim3 = float(line.split(": ")[1])
                if line.split(":")[0] == "B":
                    self.b_dim3 = float(line.split(": ")[1])
                if line.split(":")[0] == "T1":
                    self.t1_dim3 = float(line.split(": ")[1])
                if line.split(":")[0] == "T2":
                    self.t2_dim3 = float(line.split(": ")[1])
                if line.split(":")[0] == "Loc":
                    self.loc_dim3 = float(line.split(": ")[1])

                if line.split(":")[0] == "Zero Filling":
                    if "True" in line:
                        self.zero_filling_checkbox_dim3_value = True
                    else:
                        self.zero_filling_checkbox_dim3_value = False
                if line.split(":")[0] == "Zero Filling Combobox Selection":
                    self.zero_filling_dim3_combobox_selection = int(line.split(": ")[1])
                if line.split(":")[0] == "Zero Filling Value Doubling Times":
                    self.zero_filling_dim3_value_doubling_times = int(
                        line.split(": ")[1]
                    )
                if line.split(":")[0] == "Zero Filling Value Zeros to Add":
                    self.zero_filling_dim3_value_zeros_to_add = int(line.split(": ")[1])
                if line.split(":")[0] == "Zero Filling Value Final Data Size":
                    self.zero_filling_dim3_value_final_data_size = int(
                        line.split(": ")[1]
                    )
                if line.split(":")[0] == "Zero Filling Round Checkbox":
                    if "True" in line:
                        self.zero_filling_round_checkbox_dim3_value = True
                    else:
                        self.zero_filling_round_checkbox_dim3_value = False
                if line.split(":")[0] == "Fourier Transform":
                    if "True" in line:
                        self.fourier_transform_checkbox_dim3_value = True
                    else:
                        self.fourier_transform_checkbox_dim3_value = False
                if line.split(":")[0] == "Fourier Transform Method Selection":
                    self.ft_method_selection_dim3 = int(line.split(": ")[1])
                if line.split(":")[0] == "Phase Correction":
                    if "True" in line:
                        self.phase_correction_checkbox_dim3_value = True
                    else:
                        self.phase_correction_checkbox_dim3_value = False
                if line.split(":")[0] == "Phase Correction P0":
                    self.p0_total_dim3 = float(line.split(": ")[1])
                if line.split(":")[0] == "Phase Correction P1":
                    self.p1_total_dim3 = float(line.split(": ")[1])
                if line.split(":")[0] == "F1180":
                    if "True" in line:
                        self.f1180_dim3 = True
                    else:
                        self.f1180_dim3 = False
                if line.split(":")[0] == "Extraction":
                    if "True" in line:
                        self.extraction_checkbox_dim3_value = True
                    else:
                        self.extraction_checkbox_dim3_value = False
                if line.split(":")[0] == "Extraction PPM Start":
                    self.extraction_ppm_start_dim3 = float(line.split(": ")[1])
                if line.split(":")[0] == "Extraction PPM End":
                    self.extraction_ppm_end_dim3 = float(line.split(": ")[1])
                if line.split(":")[0] == "Baseline Correction":
                    if "True" in line:
                        self.baseline_correction_checkbox_dim3_value = True
                    else:
                        self.baseline_correction_checkbox_dim3_value = False
                if line.split(":")[0] == "Baseline Correction Radio Box Selection":
                    self.baseline_correction_radio_box_selection_dim3 = int(
                        line.split(": ")[1]
                    )
                if line.split(":")[0] == "Baseline Correction Nodes":
                    self.baseline_correction_nodes_dim3 = int(line.split(": ")[1])
                if line.split(":")[0] == "Baseline Correction Node List":
                    self.baseline_correction_node_list_dim3 = line.split(": ")[1]
                if line.split(":")[0] == "Baseline Correction Polynomial Order":
                    self.baseline_correction_polynomial_order_dim3 = int(
                        line.split(": ")[1]
                    )

    def set_initial_linear_prediction_variables_dim3(self):
        self.linear_prediction_radio_box_dim3_selection = 1
        self.linear_prediction_dim3_checkbox_value = False
        self.linear_prediction_dim3_options_selection = 0
        self.linear_prediction_dim3_coefficients_selection = 0
        self.linear_prediction_selection = 0

        # Check to see if the nuslist file exists in the current directory using os.path.isfile('nuslist')
        if os.path.isfile("nuslist"):
            self.nuslist_name_dim3 = "nuslist"
        else:
            self.nuslist_name_dim3 = ""

        self.number_of_nus_CPU_dim3 = 1
        self.nus_iterations_dim3 = 50
        self.nus_data_extension_dim3 = 0  # int(self.nmr_data.number_of_points[2]*1.5)

    def set_initial_apodization_variables_dim3(self):
        self.apodization_dim3_checkbox_value = True
        self.apodization_dim3_combobox_selection = 1
        self.apodization_dim3_combobox_selection_old = 1

        # Initial values for exponential apodization
        self.exponential_line_broadening_dim3 = 0.5
        self.apodization_first_point_scaling_dim3 = 0.5

        # Initial values for Lorentz to Gauss apodization
        self.g1_dim3 = 0.33
        self.g2_dim3 = 1
        self.g3_dim3 = 0.0

        # Initial values for Sinebell apodization
        self.offset_dim3 = 0.5
        self.end_dim3 = 0.98
        self.power_dim3 = 1.0

        # Initial values for Gauss Broadening apodization
        self.a_dim3 = 1.0
        self.b_dim3 = 1.0

        # Initial values for Trapezoid apodization
        self.t1_dim3 = int((self.nmr_data.number_of_points[2] / 2) / 4)
        self.t2_dim3 = int((self.nmr_data.number_of_points[2] / 2) / 4)

        # Initial values for Triangle apodization
        self.loc_dim3 = 0.5

    def set_initial_zero_filling_variables_dim3(self):
        self.zero_filling_dim3_checkbox_value = True
        self.zero_filling_dim3_combobox_selection = 0
        self.zero_filling_dim3_combobox_selection_old = 0
        self.zero_filling_dim3_value_doubling_times = 1
        self.zero_filling_dim3_value_zeros_to_add = 0
        self.zero_filling_dim3_value_final_data_size = 0
        self.zero_filling_dim3_round_checkbox_value = True

    def set_initial_fourier_transform_variables_dim3(self):
        self.fourier_transform_dim3_checkbox_value = True
        self.ft_method_selection_dim3 = (
            0  # Initially use the 'auto' method of FT as default
        )

    def set_initial_phasing_variables_dim3(self):
        self.phasing_dim3_checkbox_value = True
        self.p0_total_dim3 = 0.0
        self.p1_total_dim3 = 0.0
        self.p0_total_dim3_old = 0.0
        self.p1_total_dim3_old = 0.0
        self.phasing_from_smile = False
        self.f1180_dim3 = False

    def set_initial_extraction_variables_dim3(self):
        self.extraction_checkbox_value_dim3 = False
        self.extraction_start_dim3 = "0.0"
        self.extraction_end_dim3 = "0.0"

    def set_initial_baseline_correction_variables_dim3(self):
        self.baseline_correction_checkbox_value_dim3 = False
        self.baseline_correction_radio_box_selection_dim3 = 0
        self.node_list_dim3 = "0,5,95,100"
        self.node_width_dim3 = "2"
        self.polynomial_order_dim3 = "4"

    def create_canvas_dim3(self):

        pass

    def create_menu_bar_dim3(self):
        # Create the main sizer
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_sizer)

        # Create a sizer for the processing options for the first dimension
        self.sizer_2 = wx.BoxSizer(wx.VERTICAL)
        self.sizer_2.AddSpacer(10)

        # Create all the sizers (allow a checkbox at the top for SMILE NUS reconstruction which will change the possible options)
        # For NUS reconstruction using SMILE need to have exact phasing paramaters (first process without NUS and then calculate phase in indirect dimension, then process again using SMILE containing exact phasing parameters)
        self.create_linear_prediction_sizer_dim3(parent=self)
        self.create_apodization_sizer_dim3(parent=self)
        self.create_zero_filling_sizer_dim3(parent=self)
        self.create_fourier_transform_sizer_dim3(parent=self)
        self.create_phase_correction_sizer_dim3(parent=self)
        self.create_extraction_sizer_dim3(parent=self)
        self.create_baseline_correction_sizer_dim3(parent=self)

        self.main_sizer.Add(self.sizer_2, 0, wx.EXPAND)

        self.SetSizerAndFit(self.main_sizer)
        self.Layout()

        # Get the size of the main sizer and set the window size to 1.05 times the size of the main sizer
        self.width, self.height = self.main_sizer.GetSize()
        self.parent.parent.change_frame_size(
            int(self.width * 1.05), int(self.height * 1.25)
        )

    def create_linear_prediction_sizer_dim3(self, parent):
        # Create a sizer for the linear prediction options
        self.linear_prediction_sizer_dim3_label = wx.StaticBox(
            self, -1, "Linear Prediction/SMILE NUS Reconstruction"
        )
        self.linear_prediction_sizer_dim3 = wx.StaticBoxSizer(
            self.linear_prediction_sizer_dim3_label, wx.HORIZONTAL
        )
        self.linear_prediction_sizer_dim3.AddSpacer(10)

        # Have a radiobox for None, Linear Prediction and SMILE NUS Reconstruction
        self.linear_prediction_radio_box_dim3 = wx.RadioBox(
            parent,
            -1,
            "",
            choices=["None", "Linear Prediction", "SMILE NUS Reconstruction"],
            style=wx.RA_SPECIFY_ROWS,
        )
        self.linear_prediction_radio_box_dim3.Bind(
            wx.EVT_RADIOBOX, self.on_linear_prediction_radio_box_dim3
        )
        self.linear_prediction_radio_box_dim3.SetSelection(
            self.linear_prediction_radio_box_dim3_selection
        )

        self.linear_prediction_sizer_dim3.Add(
            self.linear_prediction_radio_box_dim3, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.linear_prediction_sizer_dim3.AddSpacer(10)

        if self.linear_prediction_radio_box_dim3.GetSelection() == 1:
            # Have a combobox for linear prediction options
            self.linear_prediction_options_text = wx.StaticText(
                parent, -1, "Add Predicted Points:"
            )
            self.linear_prediction_sizer_dim3.Add(
                self.linear_prediction_options_text, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_dim3.AddSpacer(5)
            self.linear_prediction_options = ["After FID", "Before FID"]
            self.linear_prediction_combobox_dim3 = wx.ComboBox(
                parent, -1, choices=self.linear_prediction_options, style=wx.CB_READONLY
            )
            self.linear_prediction_combobox_dim3.SetSelection(
                self.linear_prediction_dim3_options_selection
            )
            self.linear_prediction_combobox_dim3.Bind(
                wx.EVT_COMBOBOX, self.on_linear_prediction_combobox_dim3
            )
            self.linear_prediction_sizer_dim3.Add(
                self.linear_prediction_combobox_dim3, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_dim3.AddSpacer(10)
            # Have a combobox of predicted coefficient options
            self.linear_prediction_coefficients_text = wx.StaticText(
                parent, -1, "Predicted Coefficients:"
            )
            self.linear_prediction_sizer_dim3.Add(
                self.linear_prediction_coefficients_text, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_dim3.AddSpacer(5)
            self.linear_prediction_coefficients_options = [
                "Forward",
                "Backward",
                "Both",
            ]
            self.linear_prediction_coefficients_combobox_dim3 = wx.ComboBox(
                parent,
                -1,
                choices=self.linear_prediction_coefficients_options,
                style=wx.CB_READONLY,
            )
            self.linear_prediction_coefficients_combobox_dim3.SetSelection(
                self.linear_prediction_dim3_coefficients_selection
            )
            self.linear_prediction_coefficients_combobox_dim3.Bind(
                wx.EVT_COMBOBOX, self.on_linear_prediction_combobox_coefficients_dim3
            )
            self.linear_prediction_sizer_dim3.Add(
                self.linear_prediction_coefficients_combobox_dim3,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            self.linear_prediction_sizer_dim3.AddSpacer(10)
        elif self.linear_prediction_radio_box_dim3.GetSelection() == 2:
            # Have a set of options for SMILE NUS processing

            # NUS file
            self.smile_nus_file_text = wx.StaticText(parent, -1, "NUS File:")
            self.linear_prediction_sizer_dim3.Add(
                self.smile_nus_file_text, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_dim3.AddSpacer(5)

            self.smile_nus_file_textcontrol_dim3 = wx.TextCtrl(
                parent, -1, self.nuslist_name_dim3, size=(100, 20)
            )
            self.smile_nus_file_textcontrol_dim3.Bind(
                wx.EVT_TEXT, self.on_smile_nus_file_textcontrol_dim3
            )

            self.linear_prediction_sizer_dim3.Add(
                self.smile_nus_file_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_dim3.AddSpacer(10)

            # # Zero order phase correction
            # self.smile_nus_p0_text = wx.StaticText(parent, -1, 'Zero Order Phase Correction (p0):')
            # self.linear_prediction_sizer_dim3.Add(self.smile_nus_p0_text, 0, wx.ALIGN_CENTER_VERTICAL)
            # self.linear_prediction_sizer_dim3.AddSpacer(5)
            # self.smile_nus_p0_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.p0_total_dim3), size=(50, 20))
            # self.smile_nus_p0_textcontrol_dim3.Bind(wx.EVT_TEXT, self.on_smile_nus_p0_textcontrol_dim3)
            # self.linear_prediction_sizer_dim3.Add(self.smile_nus_p0_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            # self.linear_prediction_sizer_dim3.AddSpacer(10)

            # # First order phase correction
            # self.smile_nus_p1_text = wx.StaticText(parent, -1, 'First Order Phase Correction (p1):')
            # self.linear_prediction_sizer_dim3.Add(self.smile_nus_p1_text, 0, wx.ALIGN_CENTER_VERTICAL)
            # self.linear_prediction_sizer_dim3.AddSpacer(5)
            # self.smile_nus_p1_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.p1_total_dim3), size=(50, 20))
            # self.smile_nus_p1_textcontrol_dim3.Bind(wx.EVT_TEXT, self.on_smile_nus_p1_textcontrol_dim3)
            # self.linear_prediction_sizer_dim3.Add(self.smile_nus_p1_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            # self.linear_prediction_sizer_dim3.AddSpacer(10)

            # Have a data extension textcontrol
            self.smile_nus_data_extension_text = wx.StaticText(
                parent, -1, "Data Extension:"
            )
            self.linear_prediction_sizer_dim3.Add(
                self.smile_nus_data_extension_text, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_dim3.AddSpacer(5)
            self.smile_nus_data_extension_textcontrol_dim3 = wx.TextCtrl(
                parent, -1, str(self.nus_data_extension_dim3), size=(50, 20)
            )
            self.smile_nus_data_extension_textcontrol_dim3.Bind(
                wx.EVT_TEXT, self.on_smile_nus_data_extension_textcontrol_dim3
            )
            self.linear_prediction_sizer_dim3.Add(
                self.smile_nus_data_extension_textcontrol_dim3,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            self.linear_prediction_sizer_dim3.AddSpacer(10)

            # Number of CPU's
            self.smile_nus_cpu_text = wx.StaticText(parent, -1, "Number of CPU's:")
            self.linear_prediction_sizer_dim3.Add(
                self.smile_nus_cpu_text, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_dim3.AddSpacer(5)
            self.smile_nus_cpu_textcontrol_dim3 = wx.TextCtrl(
                parent, -1, str(self.number_of_nus_CPU_dim3), size=(30, 20)
            )
            self.smile_nus_cpu_textcontrol_dim3.Bind(
                wx.EVT_TEXT, self.on_smile_nus_cpu_textcontrol_dim3
            )
            self.linear_prediction_sizer_dim3.Add(
                self.smile_nus_cpu_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_dim3.AddSpacer(10)

            # Number of iterations
            self.smile_nus_iterations_text = wx.StaticText(
                parent, -1, "Number of Iterations:"
            )
            self.linear_prediction_sizer_dim3.Add(
                self.smile_nus_iterations_text, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.linear_prediction_sizer_dim3.AddSpacer(5)
            self.smile_nus_iterations_textcontrol_dim3 = wx.TextCtrl(
                parent, -1, str(self.nus_iterations_dim3), size=(30, 20)
            )
            self.smile_nus_iterations_textcontrol_dim3.Bind(
                wx.EVT_TEXT, self.on_smile_nus_iterations_textcontrol_dim3
            )
            self.linear_prediction_sizer_dim3.Add(
                self.smile_nus_iterations_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
            )

        # Have a button showing information on linear prediction
        self.linear_prediction_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))
        self.linear_prediction_info.Bind(
            wx.EVT_BUTTON, self.on_linear_prediction_info_dim3
        )
        self.linear_prediction_sizer_dim3.AddSpacer(10)
        self.linear_prediction_sizer_dim3.Add(
            self.linear_prediction_info, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.sizer_2.Add(self.linear_prediction_sizer_dim3)
        self.sizer_2.AddSpacer(10)

    def on_linear_prediction_combobox_dim3(self, event):
        # Get the selection from the combobox and update the linear prediction options
        self.linear_prediction_dim3_options_selection = (
            self.linear_prediction_combobox_dim3.GetSelection()
        )

    def on_linear_prediction_combobox_coefficients_dim3(self, event):
        # Get the selection from the combobox and update the linear prediction options
        self.linear_prediction_dim3_coefficients_selection = (
            self.linear_prediction_coefficients_combobox_dim3.GetSelection()
        )

    def on_smile_nus_file_textcontrol_dim3(self, event):
        # Get the value from the textcontrol
        self.nuslist_name_dim3 = self.smile_nus_file_textcontrol_dim3.GetValue()

    def on_smile_nus_p0_textcontrol_dim3(self, event):
        self.p0_total_dim3 = self.smile_nus_p0_textcontrol_dim3.GetValue()
        self.phasing_from_smile = True
        # Update the phasing values in the phasing section too
        self.phase_correction_p0_textcontrol_dim3.SetValue(str(self.p0_total_dim3))
        self.phasing_from_smile = False

    def on_smile_nus_p1_textcontrol_dim3(self, event):
        self.p1_total_dim3 = self.smile_nus_p1_textcontrol_dim3.GetValue()
        self.phasing_from_smile = True
        # Update the phasing values in the phasing section too
        self.phase_correction_p1_textcontrol_dim3.SetValue(str(self.p1_total_dim3))
        self.phasing_from_smile = False

    def on_smile_nus_data_extension_textcontrol_dim3(self, event):
        if self.smile_nus_data_extension_textcontrol_dim3.GetValue() != "":
            try:
                self.nus_data_extension_dim3 = int(
                    self.smile_nus_data_extension_textcontrol_dim3.GetValue()
                )
            except:
                self.nus_data_extension_dim3 = int(
                    self.nmr_data.number_of_points[2] * 1.5
                )
                # Give an error message
                error_message = wx.MessageDialog(
                    self,
                    "Data extension value needs to be an integer. Resetting to original value.",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                error_message.ShowModal()
                error_message.Destroy()
        else:
            self.nus_data_extension_dim3 = ""

    def on_smile_nus_cpu_textcontrol_dim3(self, event):
        if self.smile_nus_cpu_textcontrol_dim3.GetValue() != "":
            try:
                self.number_of_nus_CPU_dim3 = int(
                    self.smile_nus_cpu_textcontrol_dim3.GetValue()
                )
            except:
                self.number_of_nus_CPU_dim3 = 1
                # Give an error message
                error_message = wx.MessageDialog(
                    self,
                    "Number of CPU's needs to be an integer. Resetting to original value.",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                error_message.ShowModal()
                error_message.Destroy()
        else:
            self.number_of_nus_CPU_dim3 = ""

    def on_smile_nus_iterations_textcontrol_dim3(self, event):
        if self.smile_nus_iterations_textcontrol_dim3.GetValue() != "":
            try:
                self.nus_iterations_dim3 = int(
                    self.smile_nus_iterations_textcontrol_dim3.GetValue()
                )
            except:
                self.nus_iterations_dim3 = 800
                # Give an error message
                error_message = wx.MessageDialog(
                    self,
                    "Number of iterations needs to be an integer. Resetting to original value.",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                error_message.ShowModal()
                error_message.Destroy()
        else:
            self.nus_iterations_dim3 = ""

    def on_linear_prediction_info_dim3(self, event):
        # Create a popout window with information about linear prediction

        # Create a new frame
        self.linear_prediction_info_frame = wx.Frame(
            self, -1, "Linear Prediction / SMILE Information", size=(500, 500)
        )
        colour = "BLUE"

        # Create a sizer to hold the box
        self.linear_prediction_info_sizer_window = wx.BoxSizer(wx.VERTICAL)
        self.linear_prediction_info_sizer_window.AddSpacer(10)

        # Create a sizer to hold the text
        self.linear_prediction_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.linear_prediction_info_sizer.AddSpacer(10)

        # Create a text box with the information
        # Linear prediction information
        linear_prediction_information = "Linear prediction is a method used to increase the resolution of NMR spectra. It is used to predict the points of truncated FIDs (especially in indirect dimensions) and increase signal resolution.\n\n The linear prediction coefficients can be predicted using the forward FID data, backward data or an average of both directions. Then these can be used to add predicted points either before or after the current FID.\n\n Note that advanced options such as  -pred (number of predicted points) and -ord (number of predicted coefficients) can be implemented by manually added them to the nmrproc.com file."

        self.linear_prediction_info_text = wx.StaticText(
            self.linear_prediction_info_frame,
            -1,
            linear_prediction_information,
            size=(450, 200),
            style=wx.ALIGN_CENTER_HORIZONTAL,
        )

        # Add the text to the sizer
        self.linear_prediction_info_sizer.Add(
            self.linear_prediction_info_text, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.linear_prediction_info_sizer.AddSpacer(10)

        # Add a url to the nmrPipe help page
        url = "http://www.nmrscience.com/ref/nmrpipe/lp.html"
        self.linear_prediction_info_url = hl.HyperLinkCtrl(
            self.linear_prediction_info_frame,
            -1,
            "NMRPipe Help Page for Linear Prediction",
            URL=url,
        )
        self.linear_prediction_info_url.SetColours(colour, colour, colour)
        self.linear_prediction_info_url.SetUnderlines(False, False, False)
        self.linear_prediction_info_url.UpdateLink()

        # Add url to the sizer
        self.linear_prediction_info_sizer.Add(
            self.linear_prediction_info_url, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.linear_prediction_info_sizer.AddSpacer(10)

        # Add the sizer to the window sizer
        self.linear_prediction_info_sizer_window.Add(
            self.linear_prediction_info_sizer, 0, wx.ALIGN_CENTER
        )
        self.linear_prediction_info_sizer_window.AddSpacer(10)

        # Have text to explain SMILE NUS reconstruction
        smile_nus_text = "SMILE NUS reconstruction is a method used to reconstruct non-uniformly sampled data. The NUS file is a list of points that have been sampled in the FID.\nThe number of CPU's (default=1) is the number of cores that will be used to perform the reconstruction and the number of iterations can be changed to improve the accuracy (default=800).\n Furthermore, in order for accurate SMILE reconstruction, the correct zero (p0) and first (p1) order phase correction values need to be inputted."

        self.smile_nus_text = wx.StaticText(
            self.linear_prediction_info_frame,
            -1,
            smile_nus_text,
            size=(450, 200),
            style=wx.ALIGN_CENTER_HORIZONTAL,
        )
        self.linear_prediction_info_sizer_window.Add(
            self.smile_nus_text, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.linear_prediction_info_sizer_window.AddSpacer(10)

        # Add the window sizer to the frame
        self.linear_prediction_info_frame.SetSizer(
            self.linear_prediction_info_sizer_window
        )

        # Show the frame
        self.linear_prediction_info_frame.Show()

    def on_linear_prediction_radio_box_dim3(self, event):
        # Get the selection from the radio box and update the linear prediction options
        self.linear_prediction_radio_box_dim3_selection = (
            self.linear_prediction_radio_box_dim3.GetSelection()
        )

        # Remove all the old sizers and replot

        if self.apodization_dim3_combobox_selection_old == 1:
            # Remove the previous textcontrols

            self.apodization_sizer_dim3.Detach(self.apodization_line_broadening_label)
            self.apodization_line_broadening_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_line_broadening_textcontrol_dim3
            )
            self.apodization_line_broadening_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_first_point_textcontrol_dim3
            )
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)

        elif self.apodization_dim3_combobox_selection_old == 2:
            self.apodization_sizer_dim3.Detach(self.apodization_g1_label)
            self.apodization_g1_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g1_textcontrol_dim3)
            self.apodization_g1_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g2_label)
            self.apodization_g2_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g2_textcontrol_dim3)
            self.apodization_g2_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g3_label)
            self.apodization_g3_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g3_textcontrol_dim3)
            self.apodization_g3_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_first_point_textcontrol_dim3
            )
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)

        elif self.apodization_dim3_combobox_selection_old == 3:
            self.apodization_sizer_dim3.Detach(self.apodization_offset_label)
            self.apodization_offset_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_offset_textcontrol)
            self.apodization_offset_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_end_label)
            self.apodization_end_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_end_textcontrol_dim3)
            self.apodization_end_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_power_label)
            self.apodization_power_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_power_textcontrol_dim3)
            self.apodization_power_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_first_point_textcontrol_dim3
            )
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif self.apodization_dim3_combobox_selection_old == 4:
            self.apodization_sizer_dim3.Detach(self.apodization_a_label)
            self.apodization_a_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_a_textcontrol_dim3)
            self.apodization_a_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_b_label)
            self.apodization_b_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_b_textcontrol_dim3)
            self.apodization_b_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_first_point_textcontrol_dim3
            )
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif self.apodization_dim3_combobox_selection_old == 5:
            self.apodization_sizer_dim3.Detach(self.apodization_t1_label)
            self.apodization_t1_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_t1_textcontrol_dim3)
            self.apodization_t1_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_t2_label)
            self.apodization_t2_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_t2_textcontrol_dim3)
            self.apodization_t2_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_first_point_textcontrol_dim3
            )
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif self.apodization_dim3_combobox_selection_old == 6:
            self.apodization_sizer_dim3.Detach(self.apodization_loc_label)
            self.apodization_loc_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_loc_textcontrol_dim3)
            self.apodization_loc_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_first_point_textcontrol_dim3
            )
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif self.apodization_dim3_combobox_selection_old == 0:
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)

        self.apodization_sizer_dim3.Detach(self.apodization_checkbox_dim3)
        self.apodization_checkbox_dim3.Destroy()
        self.apodization_sizer_dim3.Detach(self.apodization_combobox_dim3)
        self.apodization_combobox_dim3.Hide()

        # Delete the current apodization sizer and then create a new one

        self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
        self.apodization_plot_ax_dim3.clear()
        self.apodization_plot_ax_dim3.clear()
        self.apodization_plot_sizer_dim3.Clear(True)

        self.sizer_2.Remove(self.apodization_sizer_dim3)
        # self.apodization_sizer.Clear(delete_windows=True)

        # Remove the linear prediction sizers
        self.linear_prediction_sizer_dim3.Clear(delete_windows=True)
        # self.sizer_1.Remove(self.linear_prediction_sizer)

        # self.sizer_1.Remove(self.solvent_suppression_sizer)

        self.sizer_2.Clear(delete_windows=True)

        self.create_menu_bar_dim3()
        self.Refresh()
        self.Update()
        self.Layout()

    def create_apodization_sizer_dim3(self, parent):

        # Create a box for apodization options
        self.apodization_box_dim3 = wx.StaticBox(parent, -1, "Apodization")
        self.apodization_sizer_dim3 = wx.StaticBoxSizer(
            self.apodization_box_dim3, wx.HORIZONTAL
        )
        self.apodization_checkbox_dim3 = wx.CheckBox(parent, -1, "Apply apodization")
        self.apodization_checkbox_dim3.Bind(
            wx.EVT_CHECKBOX, self.on_apodization_checkbox_dim3
        )
        self.apodization_checkbox_dim3.SetValue(self.apodization_dim3_checkbox_value)
        self.apodization_sizer_dim3.Add(
            self.apodization_checkbox_dim3, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.apodization_sizer_dim3.AddSpacer(10)
        # Have a combobox for apodization options
        self.apodization_options_dim3 = [
            "None",
            "Exponential",
            "Lorentz to Gauss",
            "Sinebell",
            "Gauss Broadening",
            "Trapazoid",
            "Triangle",
        ]
        self.apodization_combobox_dim3 = wx.ComboBox(
            parent, -1, choices=self.apodization_options_dim3, style=wx.CB_READONLY
        )
        self.apodization_combobox_dim3.SetSelection(
            self.apodization_dim3_combobox_selection
        )
        self.apodization_combobox_dim3.Bind(
            wx.EVT_COMBOBOX, self.on_apodization_combobox_dim3
        )
        self.apodization_sizer_dim3.Add(
            self.apodization_combobox_dim3, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.apodization_sizer_dim3.AddSpacer(10)
        if self.apodization_dim3_combobox_selection == 1:
            # Have a textcontrol for the line broadening
            self.apodization_line_broadening_label = wx.StaticText(
                parent, -1, "Line Broadening (Hz):"
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_line_broadening_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_line_broadening_textcontrol_dim3 = wx.TextCtrl(
                parent, -1, str(self.exponential_line_broadening_dim3), size=(30, 20)
            )
            self.apodization_line_broadening_textcontrol_dim3.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_line_broadening_textcontrol_dim3,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(
                parent, -1, "First Point Scaling:"
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_first_point_textcontrol_dim3 = wx.TextCtrl(
                parent,
                -1,
                str(self.apodization_first_point_scaling_dim3),
                size=(30, 20),
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_first_point_textcontrol_dim3,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            self.apodization_sizer_dim3.AddSpacer(10)
        elif self.apodization_dim3_combobox_selection == 2:
            # Have a textcontrol for the g1 value
            self.apodization_g1_label = wx.StaticText(
                parent, -1, "Inverse Lorentzian (Hz):"
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_g1_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_g1_textcontrol_dim3 = wx.TextCtrl(
                parent, -1, str(self.g1_dim3), size=(40, 20)
            )
            self.apodization_g1_textcontrol_dim3.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_g1_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the g2 value
            self.apodization_g2_label = wx.StaticText(
                parent, -1, "Gaussian Broadening (Hz):"
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_g2_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_g2_textcontrol_dim3 = wx.TextCtrl(
                parent, -1, str(self.g2_dim3), size=(40, 20)
            )
            self.apodization_g2_textcontrol_dim3.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_g2_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the g3 value
            self.apodization_g3_label = wx.StaticText(parent, -1, "Gaussian Shift:")
            self.apodization_sizer_dim3.Add(
                self.apodization_g3_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_g3_textcontrol_dim3 = wx.TextCtrl(
                parent, -1, str(self.g3_dim3), size=(40, 20)
            )
            self.apodization_g3_textcontrol_dim3.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_g3_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(
                parent, -1, "First Point Scaling:"
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_first_point_textcontrol_dim3 = wx.TextCtrl(
                self, -1, str(self.apodization_first_point_scaling_dim3), size=(30, 20)
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_first_point_textcontrol_dim3,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            self.apodization_sizer_dim3.AddSpacer(10)
        elif self.apodization_dim3_combobox_selection == 3:
            # Have a textcontrol for the offset value
            self.apodization_offset_label = wx.StaticText(
                parent, -1, "Offset (\u03c0):"
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_offset_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_offset_textcontrol_dim3 = wx.TextCtrl(
                parent, -1, str(self.offset_dim3), size=(40, 20)
            )
            self.apodization_offset_textcontrol_dim3.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_offset_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the end value
            self.apodization_end_label = wx.StaticText(parent, -1, "End (\u03c0):")
            self.apodization_sizer_dim3.Add(
                self.apodization_end_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_end_textcontrol_dim3 = wx.TextCtrl(
                parent, -1, str(self.end_dim3), size=(40, 20)
            )
            self.apodization_end_textcontrol_dim3.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_end_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the power value
            self.apodization_power_label = wx.StaticText(parent, -1, "Power:")
            self.apodization_sizer_dim3.Add(
                self.apodization_power_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_power_textcontrol_dim3 = wx.TextCtrl(
                parent, -1, str(self.power_dim3), size=(30, 20)
            )
            self.apodization_power_textcontrol_dim3.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_power_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(
                parent, -1, "First Point Scaling:"
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_first_point_textcontrol_dim3 = wx.TextCtrl(
                parent,
                -1,
                str(self.apodization_first_point_scaling_dim3),
                size=(30, 20),
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_first_point_textcontrol_dim3,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            self.apodization_sizer_dim3.AddSpacer(10)
        elif self.apodization_dim3_combobox_selection == 4:
            # Have a textcontrol for the a value
            self.apodization_a_label = wx.StaticText(
                parent, -1, "Line Broadening (Hz):"
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_a_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_a_textcontrol_dim3 = wx.TextCtrl(
                parent, -1, str(self.a_dim3), size=(40, 20)
            )
            self.apodization_a_textcontrol_dim3.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_a_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the b value
            self.apodization_b_label = wx.StaticText(
                parent, -1, "Gaussian Broadening (Hz):"
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_b_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_b_textcontrol_dim3 = wx.TextCtrl(
                parent, -1, str(self.b_dim3), size=(40, 20)
            )
            self.apodization_b_textcontrol_dim3.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_b_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(
                parent, -1, "First Point Scaling:"
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_first_point_textcontrol_dim3 = wx.TextCtrl(
                self, -1, str(self.apodization_first_point_scaling_dim3), size=(30, 20)
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_first_point_textcontrol_dim3,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            self.apodization_sizer_dim3.AddSpacer(10)
        elif self.apodization_dim3_combobox_selection == 5:
            # Have a textcontrol for the t1 value
            self.apodization_t1_label = wx.StaticText(parent, -1, "Ramp up points:")
            self.apodization_sizer_dim3.Add(
                self.apodization_t1_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_t1_textcontrol_dim3 = wx.TextCtrl(
                parent, -1, str(self.t1_dim3), size=(50, 20)
            )
            self.apodization_t1_textcontrol_dim3.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_t1_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the t2 value
            self.apodization_t2_label = wx.StaticText(parent, -1, "Ramp down points:")
            self.apodization_sizer_dim3.Add(
                self.apodization_t2_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_t2_textcontrol_dim3 = wx.TextCtrl(
                parent, -1, str(self.t2_dim3), size=(50, 20)
            )
            self.apodization_t2_textcontrol_dim3.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_t2_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(
                parent, -1, "First Point Scaling:"
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_first_point_textcontrol_dim3 = wx.TextCtrl(
                parent,
                -1,
                str(self.apodization_first_point_scaling_dim3),
                size=(30, 20),
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_first_point_textcontrol_dim3,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            self.apodization_sizer_dim3.AddSpacer(10)
        elif self.apodization_dim3_combobox_selection == 6:
            # Have a textcontrol for the loc value
            self.apodization_loc_label = wx.StaticText(
                parent, -1, "Location of maximum:"
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_loc_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_loc_textcontrol_dim3 = wx.TextCtrl(
                parent, -1, str(self.loc_dim3), size=(40, 20)
            )
            self.apodization_loc_textcontrol_dim3.Bind(
                wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_loc_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(
                parent, -1, "First Point Scaling:"
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.apodization_first_point_textcontrol_dim3 = wx.TextCtrl(
                parent,
                -1,
                str(self.apodization_first_point_scaling_dim3),
                size=(30, 20),
            )
            self.apodization_sizer_dim3.Add(
                self.apodization_first_point_textcontrol_dim3,
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            self.apodization_sizer_dim3.AddSpacer(10)

        # Have a button for information on currently selected apodization containing unicode i in a circle
        self.apodization_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))
        self.apodization_info.Bind(wx.EVT_BUTTON, self.oneDFrame.on_apodization_info)
        self.apodization_sizer_dim3.Add(
            self.apodization_info, 0, wx.ALIGN_CENTER_VERTICAL
        )

        # Have a mini plots of the apodization function along with the FID first slice
        self.plot_window_function_dim3()

        self.sizer_2.Add(self.apodization_sizer_dim3)
        self.sizer_2.AddSpacer(10)

    def on_apodization_checkbox_dim3(self, event):
        # Get the selection from the checkbox
        self.apodization_dim3_checkbox_selection = (
            self.apodization_checkbox_dim3.GetValue()
        )

    def on_apodization_textcontrol_dim3(self, event):
        # If the user presses enter, update the plot
        keycode = event.GetKeyCode()
        if keycode == wx.WXK_RETURN:
            self.update_window_function_plot_dim3()
        event.Skip()

    def plot_window_function_dim3(self):
        self.apodization_plot_sizer_dim3 = wx.BoxSizer(wx.VERTICAL)
        self.apodization_plot_dim3 = Figure(figsize=(1, 0.5), facecolor="#e6e6e7")
        self.apodization_plot_ax_dim3 = self.apodization_plot_dim3.add_subplot(111)
        # self.apodization_plot_ax.set_axis_off()

        self.apodization_plot_ax_dim3.set_xticks([])
        self.apodization_plot_ax_dim3.set_yticks([])

        # If the apodization function is None, make remove the axes of the plot
        if self.apodization_dim3_combobox_selection == 0:
            self.apodization_plot_ax_dim3.spines["top"].set_visible(False)
            self.apodization_plot_ax_dim3.spines["right"].set_visible(False)
            self.apodization_plot_ax_dim3.spines["bottom"].set_visible(False)
            self.apodization_plot_ax_dim3.spines["left"].set_visible(False)

        x = np.linspace(
            0,
            (self.nmr_data.number_of_points[2] / 2) / self.nmr_data.spectral_width[2],
            int(self.nmr_data.number_of_points[2] / 2),
        )
        if self.apodization_dim3_combobox_selection == 1:
            # Exponential window function
            (self.line1,) = self.apodization_plot_ax_dim3.plot(
                x,
                np.exp(-(np.pi * x * self.exponential_line_broadening_dim3)),
                color="#1f77b4",
            )
            self.apodization_plot_ax_dim3.set_ylim(0, 1.5)
            self.apodization_plot_ax_dim3.set_xlim(
                0,
                (self.nmr_data.number_of_points[2] / 2)
                / self.nmr_data.spectral_width[2],
            )

        elif self.apodization_dim3_combobox_selection == 2:
            # Lorentz to Gauss window function
            e = (
                np.pi
                * (self.nmr_data.number_of_points[2] / 2)
                / self.nmr_data.spectral_width[2]
                * self.g1_dim3
            )
            g = (
                0.6
                * np.pi
                * self.g2_dim3
                * (
                    self.g3_dim3
                    * (
                        (self.nmr_data.number_of_points[2] / 2)
                        / self.nmr_data.spectral_width[2]
                        - 1
                    )
                    - x
                )
            )
            func = np.exp(e - g * g)
            (self.line1,) = self.apodization_plot_ax_dim3.plot(x, func, color="#1f77b4")
            self.apodization_plot_ax_dim3.set_ylim(0, 1.5)
            self.apodization_plot_ax_dim3.set_xlim(
                0,
                (self.nmr_data.number_of_points[2] / 2)
                / self.nmr_data.spectral_width[2],
            )
        elif self.apodization_dim3_combobox_selection == 3:
            # Sinebell window function
            func = (
                np.sin(
                    (
                        np.pi * self.offset_dim3
                        + np.pi * (self.end_dim3 - self.offset_dim3) * x
                    )
                    / (
                        (
                            (
                                (self.nmr_data.number_of_points[2] / 2)
                                / self.nmr_data.spectral_width[2]
                            )
                        )
                    )
                )
                ** self.power_dim3
            )
            (self.line1,) = self.apodization_plot_ax_dim3.plot(x, func, color="#1f77b4")
            self.apodization_plot_ax_dim3.set_ylim(0, 1.5)
            self.apodization_plot_ax_dim3.set_xlim(
                0,
                (self.nmr_data.number_of_points[2] / 2)
                / self.nmr_data.spectral_width[2],
            )
        elif self.apodization_dim3_combobox_selection == 4:
            # Gauss broadening window function
            func = np.exp(-self.a_dim3 * (x**2) - self.b_dim3 * x)
            (self.line1,) = self.apodization_plot_ax_dim3.plot(x, func, color="#1f77b4")
            self.apodization_plot_ax_dim3.set_ylim(0, 1.5)
            self.apodization_plot_ax_dim3.set_xlim(
                0,
                (self.nmr_data.number_of_points[2] / 2)
                / self.nmr_data.spectral_width[2],
            )
        elif self.apodization_dim3_combobox_selection == 5:
            # Trapazoid window function
            func = np.concatenate(
                (
                    np.linspace(0, 1, int(self.t1_dim3)),
                    np.ones(
                        int(self.nmr_data.number_of_points[2] / 2)
                        - int(self.t1_dim3)
                        - int(self.t1_dim3)
                    ),
                    np.linspace(1, 0, int(self.t2_dim3)),
                )
            )
            (self.line1,) = self.apodization_plot_ax_dim3.plot(x, func, color="#1f77b4")
            self.apodization_plot_ax_dim3.set_ylim(0, 1.5)
            self.apodization_plot_ax_dim3.set_xlim(
                0,
                (self.nmr_data.number_of_points[2] / 2)
                / self.nmr_data.spectral_width[2],
            )
        elif self.apodization_dim3_combobox_selection == 6:
            # Triangle window function
            func = np.concatenate(
                (
                    np.linspace(
                        0,
                        1,
                        int(self.loc_dim3 * (self.nmr_data.number_of_points[2] / 2)),
                    ),
                    np.linspace(
                        1,
                        0,
                        int(
                            (1 - self.loc_dim3)
                            * (self.nmr_data.number_of_points[2] / 2)
                        ),
                    ),
                )
            )
            (self.line1,) = self.apodization_plot_ax_dim3.plot(x, func, color="#1f77b4")

            self.apodization_plot_ax_dim3.set_ylim(0, 1.5)
            self.apodization_plot_ax_dim3.set_xlim(
                0,
                (self.nmr_data.number_of_points[2] / 2)
                / self.nmr_data.spectral_width[2],
            )

        self.apodization_plot_ax_dim3.set_xlim(
            0, (self.nmr_data.number_of_points[2] / 2) / self.nmr_data.spectral_width[2]
        )

        self.apodization_plot_canvas = FigCanvas(self, -1, self.apodization_plot_dim3)
        self.apodization_plot_sizer_dim3.Add(self.apodization_plot_canvas, 0, wx.EXPAND)

        self.apodization_sizer_dim3.Add(
            self.apodization_plot_sizer_dim3, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.apodization_sizer_dim3.AddSpacer(10)

    def update_window_function_plot_dim3(self):
        x = np.linspace(
            0,
            (self.nmr_data.number_of_points[2] / 2) / self.nmr_data.spectral_width[1],
            int(self.nmr_data.number_of_points[2] / 2),
        )
        try:
            c = float(self.apodization_first_point_textcontrol_dim3.GetValue())
            self.apodization_first_point_scaling_dim3 = c
        except:
            # Give a popout window saying that the values are not valid
            msg = wx.MessageDialog(
                self,
                "The value entered for apodization first point scaling is not valid (use 0.5 or 1.0)",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            self.apodization_first_point_textcontrol_dim3.SetValue(
                str(self.apodization_first_point_scaling_dim3)
            )
            return
        if c != 0.5 and c != 1.0:
            msg = wx.MessageDialog(
                self,
                "The value entered for apodization first point scaling is not valid (use 0.5 or 1.0)",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            self.apodization_first_point_textcontrol_dim3.SetValue(
                str(self.apodization_first_point_scaling_dim3)
            )
            return
        self.apodization_first_point_scaling = c
        if self.apodization_dim3_combobox_selection == 1:
            try:
                em = float(self.apodization_line_broadening_textcontrol_dim3.GetValue())
            except:
                # Give a popout window saying that the values are not valid
                msg = wx.MessageDialog(
                    self,
                    "The values entered are not valid",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_line_broadening_textcontrol_dim3.SetValue(
                    str(self.exponential_line_broadening_dim3)
                )
                return
            self.exponential_line_broadening_dim3 = em

            self.line1.set_ydata(
                np.exp(-(np.pi * x * self.exponential_line_broadening_dim3))
            )
        elif self.apodization_dim3_combobox_selection == 2:
            try:
                g1 = float(self.apodization_g1_textcontrol_dim3.GetValue())
                g2 = float(self.apodization_g2_textcontrol_dim3.GetValue())
                g3 = float(self.apodization_g3_textcontrol_dim3.GetValue())
            except:
                # Give a popout window saying that the values are not valid
                msg = wx.MessageDialog(
                    self,
                    "The values entered are not valid",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_g1_textcontrol_dim3.SetValue(str(self.g1_dim3))
                self.apodization_g2_textcontrol_dim3.SetValue(str(self.g2_dim3))
                self.apodization_g3_textcontrol_dim3.SetValue(str(self.g3_dim3))
                return
            # Check to see if g3 is between 0 and 1
            if g3 < 0 or g3 > 1:
                msg = wx.MessageDialog(
                    self,
                    "Gaussian shift must be between 0 and 1",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_g3_textcontrol_dim3.SetValue(str(self.g3_dim3))
                return
            self.g1_dim3 = g1
            self.g2_dim3 = g2
            self.g3_dim3 = g3
            e = (
                np.pi
                * (self.nmr_data.number_of_points[2] / 2)
                / self.nmr_data.spectral_width[2]
                * self.g1_dim3
            )
            g = (
                0.6
                * np.pi
                * self.g2_dim3
                * (
                    self.g3_dim3
                    * (
                        (self.nmr_data.number_of_points[2] / 2)
                        / self.nmr_data.spectral_width[2]
                        - 1
                    )
                    - x
                )
            )
            func = np.exp(e - g * g)
            self.line1.set_ydata(func)

            self.apodization_plot_ax_dim3.set_xlim(
                0,
                (self.nmr_data.number_of_points[2] / 2)
                / self.nmr_data.spectral_width[2],
            )

        elif self.apodization_dim3_combobox_selection == 3:
            try:
                offset = float(self.apodization_offset_textcontrol_dim3.GetValue())
                end = float(self.apodization_end_textcontrol_dim3.GetValue())
                power = float(self.apodization_power_textcontrol_dim3.GetValue())
                power = int(power)
            except:
                msg = wx.MessageDialog(
                    self,
                    "The values entered are not valid",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_offset_textcontrol_dim3.SetValue(str(self.offset_dim3))
                self.apodization_end_textcontrol_dim3.SetValue(str(self.end_dim3))
                self.apodization_power_textcontrol_dim3.SetValue(str(self.power_dim3))
                return
            # Check that offset and end are between 0 and 1
            if offset < 0 or offset > 1:
                msg = wx.MessageDialog(
                    self,
                    "Offset values must be between 0 and 1",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_offset_textcontrol_dim3.SetValue(str(self.offset_dim3))
                return
            if end < 0 or end > 1:
                msg = wx.MessageDialog(
                    self,
                    "End values must be between 0 and 1",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_end_textcontrol_dim3.SetValue(str(self.end_dim3))
                return
            # Check that power is greater than 0
            if power < 0:
                msg = wx.MessageDialog(
                    self, "Power must be greater than 0", "Error", wx.OK | wx.ICON_ERROR
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_power_textcontrol_dim3.SetValue(str(self.power_dim3))
                return
            self.offset_dim3 = offset
            self.end_dim3 = end
            self.power_dim3 = power
            func = (
                np.sin(
                    (
                        np.pi * self.offset_dim3
                        + np.pi * (self.end_dim3 - self.offset_dim3) * x
                    )
                    / (
                        (
                            (
                                (self.nmr_data.number_of_points[2] / 2)
                                / self.nmr_data.spectral_width[2]
                            )
                        )
                    )
                )
                ** self.power_dim3
            )
            self.line1.set_ydata(func)
        elif self.apodization_dim3_combobox_selection == 4:
            try:
                a = float(self.apodization_a_textcontrol_dim3.GetValue())
                b = float(self.apodization_b_textcontrol_dim3.GetValue())
            except:
                msg = wx.MessageDialog(
                    self,
                    "The values entered are not valid",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_a_textcontrol_dim3.SetValue(str(self.a_dim3))
                self.apodization_b_textcontrol_dim3.SetValue(str(self.b_dim3))
                return
            self.a_dim3 = a
            self.b_dim3 = b
            func = np.exp(-self.a_dim3 * (x**2) - self.b_dim3 * x)
            self.line1.set_ydata(func)
        elif self.apodization_dim3_combobox_selection == 5:
            try:
                t1 = float(self.apodization_t1_textcontrol_dim3.GetValue())
                t2 = float(self.apodization_t2_textcontrol_dim3.GetValue())
            except:
                msg = wx.MessageDialog(
                    self,
                    "The values entered are not valid",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_t1_textcontrol_dim3.SetValue(str(self.t1_dim3))
                self.apodization_t2_textcontrol_dim3.SetValue(str(self.t2_dim3))
                return
            # Ensure that t1 and t2 are greater than 0
            if t1 < 0 or t2 < 0:
                msg = wx.MessageDialog(
                    self,
                    "Ramp up and ramp down points must be greater than 0",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_t1_textcontrol_dim3.SetValue(str(self.t1_dim3))
                self.apodization_t2_textcontrol_dim3.SetValue(str(self.t2_dim3))
                return
            # Ensure that t1 + t2 is less than the number of points
            if t1 + t2 > self.nmr_data.number_of_points[2]:
                message = (
                    "Ramp up and ramp down points must be less than the number of points ("
                    + str(self.nmr_data.number_of_points[2])
                    + ")"
                )
                msg = wx.MessageDialog(self, message, "Error", wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_t1_textcontrol_dim3.SetValue(str(self.t1_dim3))
                self.apodization_t2_textcontrol_dim3.SetValue(str(self.t2_dim3))
                return
            self.t1_dim3 = t1
            self.t2_dim3 = t2
            func = np.concatenate(
                (
                    np.linspace(0, 1, int(self.t1_dim3)),
                    np.ones(
                        int(self.nmr_data.number_of_points[2] / 2)
                        - int(self.t1_dim3)
                        - int(self.t2_dim3)
                    ),
                    np.linspace(1, 0, int(self.t2_dim3)),
                )
            )
            self.line1.set_ydata(func)
        elif self.apodization_dim3_combobox_selection == 6:
            try:
                loc = float(self.apodization_loc_textcontrol_dim3.GetValue())
            except:
                msg = wx.MessageDialog(
                    self,
                    "The values entered are not valid",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_loc_textcontrol_dim3.SetValue(str(self.loc_dim3))
                return
            # Ensure that loc is between 0 and 1
            if loc < 0 or loc > 1:
                msg = wx.MessageDialog(
                    self,
                    "Location of maximum must be between 0 and 1",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                self.apodization_loc_textcontrol_dim3.SetValue(str(self.loc_dim3))
                return
            self.loc_dim3 = loc
            func = np.concatenate(
                (
                    np.linspace(
                        0,
                        1,
                        int(self.loc_dim3 * (self.nmr_data.number_of_points[2] / 2)),
                    ),
                    np.linspace(
                        1,
                        0,
                        int(self.nmr_data.number_of_points[2] / 2)
                        - int(self.loc_dim3 * (self.nmr_data.number_of_points[2]) / 2),
                    ),
                )
            )
            self.line1.set_ydata(func)

        self.apodization_plot_canvas.draw()

    def on_apodization_combobox_dim3(self, event):
        self.apodization_dim3_combobox_selection = (
            self.apodization_combobox_dim3.GetSelection()
        )

        # Destroy the combobox and textcontrols for the previous apodization function
        # self.apodization_sizer.Detach(self.apodization_combobox)
        # self.apodization_combobox.Destroy()

        # # Remove the zf sizer
        self.zero_filling_sizer_dim3.Clear(delete_windows=True)

        if self.apodization_dim3_combobox_selection_old == 1:
            # Remove the previous textcontrols

            self.apodization_sizer_dim3.Detach(self.apodization_line_broadening_label)
            self.apodization_line_broadening_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_line_broadening_textcontrol_dim3
            )
            self.apodization_line_broadening_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_first_point_textcontrol_dim3
            )
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)

        elif self.apodization_dim3_combobox_selection_old == 2:
            self.apodization_sizer_dim3.Detach(self.apodization_g1_label)
            self.apodization_g1_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g1_textcontrol_dim3)
            self.apodization_g1_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g2_label)
            self.apodization_g2_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g2_textcontrol_dim3)
            self.apodization_g2_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g3_label)
            self.apodization_g3_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g3_textcontrol_dim3)
            self.apodization_g3_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_first_point_textcontrol_dim3
            )
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)

        elif self.apodization_dim3_combobox_selection_old == 3:
            self.apodization_sizer_dim3.Detach(self.apodization_offset_label)
            self.apodization_offset_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_offset_textcontrol_dim3)
            self.apodization_offset_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_end_label)
            self.apodization_end_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_end_textcontrol_dim3)
            self.apodization_end_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_power_label)
            self.apodization_power_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_power_textcontrol_dim3)
            self.apodization_power_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_first_point_textcontrol_dim3
            )
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif self.apodization_dim3_combobox_selection_old == 4:
            self.apodization_sizer_dim3.Detach(self.apodization_a_label)
            self.apodization_a_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_a_textcontrol_dim3)
            self.apodization_a_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_b_label)
            self.apodization_b_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_b_textcontrol_dim3)
            self.apodization_b_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_first_point_textcontrol_dim3
            )
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif self.apodization_dim3_combobox_selection_old == 5:
            self.apodization_sizer_dim3.Detach(self.apodization_t1_label)
            self.apodization_t1_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_t1_textcontrol_dim3)
            self.apodization_t1_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_t2_label)
            self.apodization_t2_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_t2_textcontrol_dim3)
            self.apodization_t2_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_first_point_textcontrol_dim3
            )
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif self.apodization_dim3_combobox_selection_old == 6:
            self.apodization_sizer_dim3.Detach(self.apodization_loc_label)
            self.apodization_loc_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_loc_textcontrol_dim3)
            self.apodization_loc_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_first_point_textcontrol_dim3
            )
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif self.apodization_dim3_combobox_selection_old == 0:
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)

        self.apodization_sizer_dim3.Detach(self.apodization_checkbox_dim3)
        self.apodization_checkbox_dim3.Destroy()
        self.apodization_sizer_dim3.Detach(self.apodization_combobox_dim3)
        self.apodization_combobox_dim3.Hide()

        # Delete the current apodization sizer and then create a new one

        self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
        self.apodization_plot_ax_dim3.clear()
        self.apodization_plot_ax_dim3.clear()
        self.apodization_plot_sizer_dim3.Clear(True)

        self.sizer_2.Remove(self.apodization_sizer_dim3)
        # self.apodization_sizer.Clear(delete_windows=True)

        # Remove the linear prediction sizers
        self.linear_prediction_sizer_dim3.Clear(delete_windows=True)
        # self.sizer_1.Remove(self.linear_prediction_sizer)

        # self.sizer_1.Remove(self.solvent_suppression_sizer)

        self.sizer_2.Clear(delete_windows=True)

        self.create_menu_bar_dim3()
        self.Refresh()
        self.Update()
        self.Layout()

        self.apodization_dim3_combobox_selection_old = (
            self.apodization_dim3_combobox_selection
        )

    def create_zero_filling_sizer_dim3(self, parent):
        # Create a box for zero filling options
        self.zero_filling_box_dim3 = wx.StaticBox(parent, -1, "Zero Filling")
        self.zero_filling_sizer_dim3 = wx.StaticBoxSizer(
            self.zero_filling_box_dim3, wx.HORIZONTAL
        )
        self.zero_filling_checkbox_dim3 = wx.CheckBox(parent, -1, "Apply zero filling")
        self.zero_filling_checkbox_dim3.SetValue(self.zero_filling_dim3_checkbox_value)
        self.zero_filling_checkbox_dim3.Bind(
            wx.EVT_CHECKBOX, self.on_zero_filling_checkbox_dim3
        )
        self.zero_filling_sizer_dim3.Add(
            self.zero_filling_checkbox_dim3, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.zero_filling_sizer_dim3.AddSpacer(10)
        # Have a combobox for zero filling options
        self.zf_options_label = wx.StaticText(parent, -1, "Options:")
        self.zero_filling_sizer_dim3.Add(
            self.zf_options_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.zero_filling_sizer_dim3.AddSpacer(5)
        self.zero_filling_options_dim3 = [
            "Doubling spectrum size",
            "Adding Zeros",
            "Final data size",
        ]
        self.zero_filling_combobox_dim3 = wx.ComboBox(
            parent, -1, choices=self.zero_filling_options_dim3, style=wx.CB_READONLY
        )
        self.zero_filling_combobox_dim3.Bind(
            wx.EVT_COMBOBOX, self.on_zero_filling_combobox_dim3
        )
        self.zero_filling_combobox_dim3.SetSelection(
            self.zero_filling_dim3_combobox_selection
        )
        self.zero_filling_sizer_dim3.Add(
            self.zero_filling_combobox_dim3, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.zero_filling_sizer_dim3.AddSpacer(10)
        if self.zero_filling_dim3_combobox_selection == 0:
            # Have a textcontrol for the doubling number/number of zeros/final data size
            self.zf_value_label = wx.StaticText(parent, -1, "Doubling number:")
            self.zero_filling_sizer_dim3.Add(
                self.zf_value_label, 0, wx.ALIGN_CENTER_VERTICAL
            )

            self.zero_filling_textcontrol_dim3 = wx.TextCtrl(
                parent,
                -1,
                str(self.zero_filling_dim3_value_doubling_times),
                size=(40, 20),
            )
            self.zero_filling_sizer_dim3.AddSpacer(5)
            self.zero_filling_sizer_dim3.Add(
                self.zero_filling_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.zero_filling_textcontrol_dim3.Bind(
                wx.EVT_TEXT, self.on_zero_filling_doubling_number_dim3
            )
            self.zero_filling_sizer_dim3.AddSpacer(20)
        elif self.zero_filling_dim3_combobox_selection == 1:
            # Have a textcontrol for the doubling number/number of zeros/final data size
            self.zf_value_label = wx.StaticText(parent, -1, "Number of zeros to add:")
            self.zero_filling_sizer_dim3.Add(
                self.zf_value_label, 0, wx.ALIGN_CENTER_VERTICAL
            )

            self.zero_filling_textcontrol_dim3 = wx.TextCtrl(
                parent,
                -1,
                str(self.zero_filling_dim3_value_zeros_to_add),
                size=(40, 20),
            )
            self.zero_filling_sizer_dim3.AddSpacer(5)
            self.zero_filling_sizer_dim3.Add(
                self.zero_filling_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.zero_filling_textcontrol_dim3.Bind(
                wx.EVT_TEXT, self.on_zero_filling_zeros_to_add_dim3
            )
            self.zero_filling_sizer_dim3.AddSpacer(20)
        elif self.zero_filling_dim3_combobox_selection == 2:
            # Have a textcontrol for the doubling number/number of zeros/final data size
            self.zf_value_label = wx.StaticText(parent, -1, "Final data size:")
            self.zero_filling_sizer_dim3.Add(
                self.zf_value_label, 0, wx.ALIGN_CENTER_VERTICAL
            )

            self.zero_filling_textcontrol_dim3 = wx.TextCtrl(
                parent,
                -1,
                str(self.zero_filling_dim3_value_final_data_size),
                size=(40, 20),
            )
            self.zero_filling_textcontrol_dim3.Bind(
                wx.EVT_TEXT, self.on_zero_filling_final_size_dim3
            )
            self.zero_filling_sizer_dim3.AddSpacer(5)
            self.zero_filling_sizer_dim3.Add(
                self.zero_filling_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
            )
            self.zero_filling_sizer_dim3.AddSpacer(20)

        # Have a checkbox for rounding to the nearest power of 2
        self.zero_filling_round_checkbox_dim3 = wx.CheckBox(
            parent, -1, "Round to nearest power of 2"
        )
        self.zero_filling_round_checkbox_dim3.SetValue(
            self.zero_filling_dim3_round_checkbox_value
        )
        self.zero_filling_round_checkbox_dim3.Bind(
            wx.EVT_CHECKBOX, self.on_zero_filling_round_checkbox_dim3
        )
        self.zero_filling_sizer_dim3.Add(
            self.zero_filling_round_checkbox_dim3, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.zero_filling_sizer_dim3.AddSpacer(10)

        # Have a button showing information on zero filling
        self.zero_filling_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))
        self.zero_filling_info.Bind(wx.EVT_BUTTON, self.oneDFrame.on_zero_fill_info)
        self.zero_filling_sizer_dim3.Add(
            self.zero_filling_info, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.zero_filling_sizer_dim3.AddSpacer(10)

        self.sizer_2.Add(self.zero_filling_sizer_dim3)
        self.sizer_2.AddSpacer(10)

    def on_zero_filling_checkbox_dim3(self, event):
        self.zero_filling_dim3_checkbox_value = (
            self.zero_filling_checkbox_dim3.GetValue()
        )

    def on_zero_filling_round_checkbox_dim3(self, event):
        self.zero_filling_dim3_round_checkbox_value = (
            self.zero_filling_round_checkbox_dim3.GetValue()
        )

    def on_zero_filling_final_size_dim3(self, event):
        self.zero_filling_dim3_value_final_data_size = (
            self.zero_filling_textcontrol_dim3.GetValue()
        )

    def on_zero_filling_zeros_to_add_dim3(self, event):
        self.zero_filling_dim3_value_zeros_to_add = (
            self.zero_filling_textcontrol_dim3.GetValue()
        )

    def on_zero_filling_doubling_number_dim3(self, event):
        self.zero_filling_dim3_value_doubling_times = (
            self.zero_filling_textcontrol_dim3.GetValue()
        )

    def on_zero_filling_combobox_dim3(self, event):
        self.zero_filling_dim3_combobox_selection = (
            self.zero_filling_combobox_dim3.GetSelection()
        )
        # # # Remove the zf sizer
        self.zero_filling_sizer_dim3.Clear()
        self.zero_filling_sizer_dim3.Detach(self.zero_filling_checkbox_dim3)
        self.zero_filling_checkbox_dim3.Destroy()
        self.zero_filling_sizer_dim3.Detach(self.zf_options_label)
        self.zf_options_label.Destroy()
        self.zero_filling_sizer_dim3.Detach(self.zero_filling_info)
        self.zero_filling_info.Destroy()
        self.zero_filling_sizer_dim3.Detach(self.zf_value_label)
        self.zf_value_label.Destroy()
        self.zero_filling_sizer_dim3.Detach(self.zero_filling_round_checkbox_dim3)
        self.zero_filling_round_checkbox_dim3.Destroy()
        self.zero_filling_sizer_dim3.Detach(self.zero_filling_textcontrol_dim3)
        self.zero_filling_textcontrol_dim3.Destroy()

        self.zero_filling_sizer_dim3.Detach(self.zero_filling_combobox_dim3)
        self.zero_filling_combobox_dim3.Hide()

        if self.apodization_dim3_combobox_selection_old == 1:
            # Remove the previous textcontrols

            self.apodization_sizer_dim3.Detach(self.apodization_line_broadening_label)
            self.apodization_line_broadening_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_line_broadening_textcontrol_dim3
            )
            self.apodization_line_broadening_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_first_point_textcontrol_dim3
            )
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)

        elif self.apodization_dim3_combobox_selection_old == 2:
            self.apodization_sizer_dim3.Detach(self.apodization_g1_label)
            self.apodization_g1_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g1_textcontrol_dim3)
            self.apodization_g1_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g2_label)
            self.apodization_g2_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g2_textcontrol_dim3)
            self.apodization_g2_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g3_label)
            self.apodization_g3_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g3_textcontrol_dim3)
            self.apodization_g3_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_first_point_textcontrol_dim3
            )
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)

        elif self.apodization_dim3_combobox_selection_old == 3:
            self.apodization_sizer_dim3.Detach(self.apodization_offset_label)
            self.apodization_offset_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_offset_textcontrol_dim3)
            self.apodization_offset_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_end_label)
            self.apodization_end_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_end_textcontrol_dim3)
            self.apodization_end_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_power_label)
            self.apodization_power_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_power_textcontrol_dim3)
            self.apodization_power_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_first_point_textcontrol_dim3
            )
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif self.apodization_dim3_combobox_selection_old == 4:
            self.apodization_sizer_dim3.Detach(self.apodization_a_label)
            self.apodization_a_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_a_textcontrol_dim3)
            self.apodization_a_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_b_label)
            self.apodization_b_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_b_textcontrol_dim3)
            self.apodization_b_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_first_point_textcontrol_dim3
            )
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif self.apodization_dim3_combobox_selection_old == 5:
            self.apodization_sizer_dim3.Detach(self.apodization_t1_label)
            self.apodization_t1_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_t1_textcontrol_dim3)
            self.apodization_t1_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_t2_label)
            self.apodization_t2_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_t2_textcontrol_dim3)
            self.apodization_t2_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_first_point_textcontrol_dim3
            )
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif self.apodization_dim3_combobox_selection_old == 6:
            self.apodization_sizer_dim3.Detach(self.apodization_loc_label)
            self.apodization_loc_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_loc_textcontrol_dim3)
            self.apodization_loc_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(
                self.apodization_first_point_textcontrol_dim3
            )
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif self.apodization_dim3_combobox_selection_old == 0:
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)

        self.apodization_sizer_dim3.Detach(self.apodization_checkbox_dim3)
        self.apodization_checkbox_dim3.Destroy()
        self.apodization_sizer_dim3.Detach(self.apodization_combobox_dim3)
        self.apodization_combobox_dim3.Hide()

        # Delete the current apodization sizer and then create a new one

        self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
        self.apodization_plot_ax_dim3.clear()
        self.apodization_plot_ax_dim3.clear()
        self.apodization_plot_sizer_dim3.Clear(True)

        self.sizer_2.Remove(self.apodization_sizer_dim3)
        # self.apodization_sizer.Clear(delete_windows=True)

        # Remove the linear prediction sizers
        self.linear_prediction_sizer_dim3.Clear(delete_windows=True)
        # self.sizer_1.Remove(self.linear_prediction_sizer)

        # self.sizer_1.Remove(self.solvent_suppression_sizer)

        self.sizer_2.Clear(delete_windows=True)

        self.create_menu_bar_dim3()
        self.Refresh()
        self.Update()
        self.Layout()

    def create_fourier_transform_sizer_dim3(self, parent):
        # Create a box for fourier transform options
        self.fourier_transform_box = wx.StaticBox(parent, -1, "Fourier Transform")
        self.fourier_transform_sizer_dim3 = wx.StaticBoxSizer(
            self.fourier_transform_box, wx.HORIZONTAL
        )
        self.fourier_transform_checkbox_dim3 = wx.CheckBox(
            parent, -1, "Apply fourier transform"
        )
        self.fourier_transform_checkbox_dim3.Bind(
            wx.EVT_CHECKBOX, self.on_fourier_transform_checkbox_dim3
        )
        self.fourier_transform_checkbox_dim3.SetValue(
            self.fourier_transform_dim3_checkbox_value
        )
        self.fourier_transform_sizer_dim3.Add(
            self.fourier_transform_checkbox_dim3, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.fourier_transform_sizer_dim3.AddSpacer(10)
        # Have a button for advanced options for fourier transform
        self.fourier_transform_advanced_options_dim3 = wx.Button(
            parent, -1, "Advanced Options"
        )
        self.fourier_transform_advanced_options_dim3.Bind(
            wx.EVT_BUTTON, self.on_fourier_transform_advanced_options_dim3
        )
        self.fourier_transform_sizer_dim3.Add(
            self.fourier_transform_advanced_options_dim3, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.fourier_transform_sizer_dim3.AddSpacer(10)

        # Have a button showing information on fourier transform
        self.fourier_transform_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))
        self.fourier_transform_info.Bind(
            wx.EVT_BUTTON, self.oneDFrame.on_fourier_transform_info
        )
        self.fourier_transform_sizer_dim3.Add(
            self.fourier_transform_info, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.sizer_2.Add(self.fourier_transform_sizer_dim3)
        self.sizer_2.AddSpacer(10)

    def on_fourier_transform_checkbox_dim3(self, event):
        self.fourier_transform_dim3_checkbox_value = (
            self.fourier_transform_checkbox_dim3.GetValue()
        )

    def on_fourier_transform_advanced_options_dim3(self, event):
        # Create a frame with a set of advanced options for the fourier transform implementation
        self.fourier_transform_advanced_options_window_dim3 = wx.Frame(
            self,
            -1,
            "Fourier Transform Advanced Options (Dimension 2)",
            size=(700, 300),
        )

        self.fourier_transform_advanced_options_window_sizer_dim3 = wx.BoxSizer(
            wx.VERTICAL
        )
        self.fourier_transform_advanced_options_window_dim3.SetSizer(
            self.fourier_transform_advanced_options_window_sizer_dim3
        )

        # Create a sizer for the fourier transform advanced options
        self.ft_label = wx.StaticBox(
            self.fourier_transform_advanced_options_window_dim3,
            -1,
            "Fourier Transform Method:",
        )
        self.fourier_transform_advanced_options_sizer_dim3 = wx.StaticBoxSizer(
            self.ft_label, wx.VERTICAL
        )

        # Have a radiobox for auto, real, inverse, sign alternation
        self.fourier_transform_advanced_options_sizer_dim3.AddSpacer(10)
        self.fourier_transform_auto_real_inverse_sign_alternation_radio_box_dim3 = (
            wx.RadioBox(
                self.fourier_transform_advanced_options_window_dim3,
                -1,
                choices=["Auto", "Real", "Inverse", "Sign Alternation", "Negative"],
                style=wx.RA_SPECIFY_COLS,
            )
        )
        self.fourier_transform_auto_real_inverse_sign_alternation_radio_box_dim3.SetSelection(
            self.ft_method_selection_dim3
        )
        self.fourier_transform_advanced_options_sizer_dim3.Add(
            self.fourier_transform_auto_real_inverse_sign_alternation_radio_box_dim3,
            0,
            wx.ALIGN_CENTER_HORIZONTAL,
        )
        self.fourier_transform_advanced_options_sizer_dim3.AddSpacer(10)

        self.ft_method_text = "Auto: The auto method will automatically select the best method for the fourier transform of the FID. \n\nReal: The Fourier Transform will be applied to the real part of the FID only. \n\nInverse: The inverse Fourier Transform will be applied to the FID. \n\nSign Alternation: The sign alternation method will be applied to the FID. \n\n"

        self.ft_method_info = wx.StaticText(
            self.fourier_transform_advanced_options_window_dim3, -1, self.ft_method_text
        )
        self.fourier_transform_advanced_options_sizer_dim3.Add(
            self.ft_method_info, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.fourier_transform_advanced_options_sizer_dim3.AddSpacer(10)

        # Have a save and close button
        self.fourier_transform_advanced_options_save_button_dim3 = wx.Button(
            self.fourier_transform_advanced_options_window_dim3, -1, "Save and Close"
        )
        self.fourier_transform_advanced_options_save_button_dim3.Bind(
            wx.EVT_BUTTON, self.on_fourier_transform_advanced_options_save_dim3
        )
        self.fourier_transform_advanced_options_sizer_dim3.Add(
            self.fourier_transform_advanced_options_save_button_dim3,
            0,
            wx.ALIGN_CENTER_HORIZONTAL,
        )

        self.fourier_transform_advanced_options_window_sizer_dim3.Add(
            self.fourier_transform_advanced_options_sizer_dim3,
            0,
            wx.ALIGN_CENTER_HORIZONTAL,
        )

        self.fourier_transform_advanced_options_window_dim3.Show()

    def on_fourier_transform_advanced_options_save_dim3(self, event):
        # Save the current selection and close the window
        self.ft_method_selection_dim3 = (
            self.fourier_transform_auto_real_inverse_sign_alternation_radio_box_dim3.GetSelection()
        )
        self.fourier_transform_advanced_options_window_dim3.Close()

    def create_phase_correction_sizer_dim3(self, parent):
        # Create a box for phase correction options
        self.phase_correction_box_dim3 = wx.StaticBox(parent, -1, "Phase Correction")
        self.phase_correction_sizer_dim3 = wx.StaticBoxSizer(
            self.phase_correction_box_dim3, wx.HORIZONTAL
        )
        self.phase_correction_checkbox_dim3 = wx.CheckBox(
            parent, -1, "Apply phase correction"
        )
        self.phase_correction_checkbox_dim3.Bind(
            wx.EVT_CHECKBOX, self.on_phase_correction_checkbox_dim3
        )
        self.phase_correction_checkbox_dim3.SetValue(self.phasing_dim3_checkbox_value)
        self.phase_correction_sizer_dim3.Add(
            self.phase_correction_checkbox_dim3, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.phase_correction_sizer_dim3.AddSpacer(10)
        # Have a textcontrol for p0 and p1 values
        self.phase_correction_p0_label = wx.StaticText(
            parent, -1, "Zero order correction (p0):"
        )
        self.phase_correction_sizer_dim3.Add(
            self.phase_correction_p0_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.phase_correction_p0_textcontrol_dim3 = wx.TextCtrl(
            parent, -1, str(self.p0_total_dim3), size=(50, 20)
        )
        self.phase_correction_p0_textcontrol_dim3.Bind(
            wx.EVT_TEXT, self.on_phase_correction_p0_dim3
        )
        self.phase_correction_sizer_dim3.Add(
            self.phase_correction_p0_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.phase_correction_sizer_dim3.AddSpacer(10)
        self.phase_correction_p1_label = wx.StaticText(
            parent, -1, "First order correction (p1):"
        )
        self.phase_correction_sizer_dim3.Add(
            self.phase_correction_p1_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.phase_correction_p1_textcontrol_dim3 = wx.TextCtrl(
            parent, -1, str(self.p1_total_dim3), size=(50, 20)
        )
        self.phase_correction_p1_textcontrol_dim3.Bind(
            wx.EVT_TEXT, self.on_phase_correction_p1_dim3
        )
        self.phase_correction_sizer_dim3.Add(
            self.phase_correction_p1_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.phase_correction_sizer_dim3.AddSpacer(10)

        # Have a checkbox for f1180
        self.phase_correction_f1180_button_dim3 = wx.CheckBox(parent, -1, "F1180")
        self.phase_correction_f1180_button_dim3.Bind(
            wx.EVT_CHECKBOX, self.on_phase_correction_f1180_dim3
        )
        self.phase_correction_f1180_button_dim3.SetValue(self.f1180_dim3)
        self.phase_correction_sizer_dim3.Add(
            self.phase_correction_f1180_button_dim3, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.phase_correction_sizer_dim3.AddSpacer(10)

        # Have a button showing information on phase correction
        self.phase_correction_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))
        self.phase_correction_info.Bind(
            wx.EVT_BUTTON, self.on_phase_correction_info_dim3
        )
        self.phase_correction_sizer_dim3.Add(
            self.phase_correction_info, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.sizer_2.Add(self.phase_correction_sizer_dim3)
        self.sizer_2.AddSpacer(10)

    def on_phase_correction_checkbox_dim3(self, event):
        self.phasing_dim3_checkbox_value = (
            self.phase_correction_checkbox_dim3.GetValue()
        )

    def on_phase_correction_p0_dim3(self, event):
        self.p0_total_dim3 = self.phase_correction_p0_textcontrol_dim3.GetValue()

    def on_phase_correction_p1_dim3(self, event):
        self.p1_total_dim3 = self.phase_correction_p1_textcontrol_dim3.GetValue()
        try:
            if np.abs(float(self.p1_total_dim3)) > 45:
                self.apodization_first_point_scaling_dim3 = 1.0
                self.apodization_first_point_textcontrol_dim3.SetValue(
                    str(self.apodization_first_point_scaling_dim3)
                )
            else:
                self.apodization_first_point_scaling_dim3 = 0.0
                self.apodization_first_point_textcontrol_dim3.SetValue(
                    str(self.apodization_first_point_scaling_dim3)
                )
        except:
            pass

    def on_phase_correction_f1180_dim3(self, event):
        if self.phase_correction_f1180_button_dim3.GetValue() == True:
            self.c_old = self.apodization_first_point_scaling_dim3
            self.p0_total_dim3_old = self.p0_total_dim3
            self.p1_total_dim3_old = self.p1_total_dim3
            # Apply -90 p0 and 180 p1 to the phase correction textcontrols
            self.p0_total_dim3 = -90.0
            self.p1_total_dim3 = 180.0
            self.phase_correction_p0_textcontrol_dim3.SetValue(str(self.p0_total_dim3))
            self.phase_correction_p1_textcontrol_dim3.SetValue(str(self.p1_total_dim3))
            # Disable the phase correction textcontrols
            self.phase_correction_p0_textcontrol_dim3.Disable()
            self.phase_correction_p1_textcontrol_dim3.Disable()

            self.apodization_first_point_scaling_dim3 = 1.0
            self.apodization_first_point_textcontrol_dim3.SetValue(
                str(self.apodization_first_point_scaling_dim3)
            )
        else:
            self.p0_total_dim3 = self.p0_total_dim3_old
            self.p1_total_dim3 = self.p1_total_dim3_old
            self.phase_correction_p0_textcontrol_dim3.SetValue(str(self.p0_total_dim3))
            self.phase_correction_p1_textcontrol_dim3.SetValue(str(self.p1_total_dim3))
            self.phase_correction_p0_textcontrol_dim3.Enable()
            self.phase_correction_p1_textcontrol_dim3.Enable()
            self.apodization_first_point_scaling_dim3 = self.c_old
            self.apodization_first_point_textcontrol_dim3.SetValue(
                str(self.apodization_first_point_scaling_dim3)
            )

    def on_phase_correction_info_dim3(self, event):
        phase_correction_text = "Phase correction is a method to correct for phase errors in the FID. Zero order phase correction (p0) is used to correct a phase offset that is applied equally across the spectrum. However, a first order phase correction (p1) is used to correct the phasing in a spectrum where peaks in different locations of the spectrum require a different phasing value. For the indirect dimension, it is often the case that the acquisition is delayed by an exact time so that the resulting spectrum can be phased using the phase values of: p0=-90, p1=180. This is often termed F1180. \n Further information can be found using the link below."

        # Create a popup window with the information
        self.phase_correction_info_window = wx.Frame(
            self, -1, "Phase Correction Information", size=(450, 300)
        )
        colour = "BLUE"

        self.phase_correction_info_window_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.phase_correction_info_window.SetSizer(
            self.phase_correction_info_window_sizer
        )

        self.phase_correction_info_window_sizer.AddSpacer(10)

        # Create a sizer for the phase correction information
        self.phase_correction_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.phase_correction_info_sizer.AddSpacer(10)
        self.phase_correction_info_sizer.Add(
            wx.StaticText(
                self.phase_correction_info_window,
                -1,
                phase_correction_text,
                size=(400, 200),
            ),
            0,
            wx.ALIGN_CENTER,
        )
        self.phase_correction_info_sizer.AddSpacer(10)

        # Have a hyperlink to the phase correction information
        self.phase_correction_info_hyperlink = hl.HyperLinkCtrl(
            self.phase_correction_info_window,
            -1,
            "NMRPipe Help Page for Phase Correction",
            URL="http://www.nmrscience.com/ref/nmrpipe/ps.html",
        )
        self.phase_correction_info_hyperlink.SetColours(colour, colour, colour)
        self.phase_correction_info_hyperlink.SetUnderlines(False, False, False)
        self.phase_correction_info_hyperlink.SetBold(False)
        self.phase_correction_info_hyperlink.UpdateLink()
        self.phase_correction_info_sizer.Add(
            self.phase_correction_info_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.phase_correction_info_sizer.AddSpacer(10)

        self.phase_correction_info_window_sizer.Add(
            self.phase_correction_info_sizer, 0, wx.ALIGN_CENTER
        )

        self.phase_correction_info_window.Show()

    def create_extraction_sizer_dim3(self, parent):
        # A box for extraction of data between two ppm values
        self.extraction_box_dim3 = wx.StaticBox(parent, -1, "Extraction")
        self.extraction_sizer_dim3 = wx.StaticBoxSizer(
            self.extraction_box_dim3, wx.HORIZONTAL
        )
        self.extraction_checkbox_dim3 = wx.CheckBox(
            parent, -1, "Include data extraction"
        )
        self.extraction_checkbox_dim3.Bind(
            wx.EVT_CHECKBOX, self.on_extraction_checkbox_dim3
        )
        self.extraction_checkbox_dim3.SetValue(self.extraction_checkbox_value_dim3)
        self.extraction_sizer_dim3.Add(
            self.extraction_checkbox_dim3, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.extraction_sizer_dim3.AddSpacer(10)
        # Have a textcontrol for the ppm start value
        self.extraction_ppm_start_label = wx.StaticText(
            parent, -1, "Start chemical shift (ppm):"
        )
        self.extraction_sizer_dim3.Add(
            self.extraction_ppm_start_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.extraction_ppm_start_textcontrol_dim3 = wx.TextCtrl(
            parent, -1, str(self.extraction_start_dim3), size=(40, 20)
        )
        self.extraction_ppm_start_textcontrol_dim3.Bind(
            wx.EVT_TEXT, self.on_extraction_dim3
        )
        self.extraction_sizer_dim3.Add(
            self.extraction_ppm_start_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.extraction_sizer_dim3.AddSpacer(10)
        # Have a textcontrol for the ppm end value
        self.extraction_ppm_end_label = wx.StaticText(
            parent, -1, "End chemical shift (ppm):"
        )
        self.extraction_sizer_dim3.Add(
            self.extraction_ppm_end_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.extraction_ppm_end_textcontrol_dim3 = wx.TextCtrl(
            parent, -1, str(self.extraction_end_dim3), size=(40, 20)
        )
        self.extraction_ppm_end_textcontrol_dim3.Bind(
            wx.EVT_TEXT, self.on_extraction_dim3
        )
        self.extraction_sizer_dim3.Add(
            self.extraction_ppm_end_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.extraction_sizer_dim3.AddSpacer(10)
        # Have a button showing information on extraction
        self.extraction_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))
        self.extraction_info.Bind(wx.EVT_BUTTON, self.oneDFrame.on_extraction_info)
        self.extraction_sizer_dim3.Add(
            self.extraction_info, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.sizer_2.Add(self.extraction_sizer_dim3)
        self.sizer_2.AddSpacer(10)

    def on_extraction_checkbox_dim3(self, event):
        self.extraction_checkbox_value_dim3 = self.extraction_checkbox_dim3.GetValue()

    def on_extraction_dim3(self, event):
        self.extraction_start_dim3 = (
            self.extraction_ppm_start_textcontrol_dim3.GetValue()
        )
        self.extraction_end_dim3 = self.extraction_ppm_end_textcontrol_dim3.GetValue()

    def create_baseline_correction_sizer_dim3(self, parent):
        # Create a box for baseline correction options (linear/polynomial)
        self.baseline_correction_box_dim3 = wx.StaticBox(
            parent, -1, "Baseline Correction"
        )
        self.baseline_correction_sizer_dim3 = wx.StaticBoxSizer(
            self.baseline_correction_box_dim3, wx.HORIZONTAL
        )
        self.baseline_correction_checkbox_dim3 = wx.CheckBox(
            parent, -1, "Apply baseline correction"
        )
        self.baseline_correction_checkbox_dim3.Bind(
            wx.EVT_CHECKBOX, self.on_baseline_correction_checkbox_dim3
        )
        self.baseline_correction_checkbox_dim3.SetValue(
            self.baseline_correction_checkbox_value_dim3
        )
        self.baseline_correction_sizer_dim3.Add(
            self.baseline_correction_checkbox_dim3, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.baseline_correction_sizer_dim3.AddSpacer(10)
        # Have a radio box for linear or polynomial baseline correction
        self.baseline_correction_radio_box_dim3 = wx.RadioBox(
            parent, -1, "Baseline Correction Method", choices=["Linear", "Polynomial"]
        )
        # Bind the radio box to a function that will update the baseline correction options
        self.baseline_correction_radio_box_dim3.Bind(
            wx.EVT_RADIOBOX, self.on_baseline_correction_radio_box_dim3
        )
        self.baseline_correction_radio_box_dim3.SetSelection(
            self.baseline_correction_radio_box_selection_dim3
        )
        self.baseline_correction_sizer_dim3.Add(
            self.baseline_correction_radio_box_dim3, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.baseline_correction_sizer_dim3.AddSpacer(10)

        # If linear baseline correction is selected, have a textcontrol for the node values to use
        self.baseline_correction_nodes_label = wx.StaticText(
            parent, -1, "Node width (pts):"
        )
        self.baseline_correction_sizer_dim3.Add(
            self.baseline_correction_nodes_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.baseline_correction_nodes_textcontrol_dim3 = wx.TextCtrl(
            parent, -1, self.node_width_dim3, size=(30, 20)
        )
        self.baseline_correction_nodes_textcontrol_dim3.Bind(
            wx.EVT_TEXT, self.on_baseline_correction_textcontrol_dim3
        )
        self.baseline_correction_sizer_dim3.Add(
            self.baseline_correction_nodes_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.baseline_correction_sizer_dim3.AddSpacer(10)
        # Have a textcontrol for the node list (percentages)
        self.baseline_correction_node_list_label = wx.StaticText(
            parent, -1, "Node list (%):"
        )
        self.baseline_correction_sizer_dim3.Add(
            self.baseline_correction_node_list_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.baseline_correction_node_list_textcontrol_dim3 = wx.TextCtrl(
            parent, -1, self.node_list_dim3, size=(100, 20)
        )
        self.baseline_correction_node_list_textcontrol_dim3.Bind(
            wx.EVT_TEXT, self.on_baseline_correction_textcontrol_dim3
        )
        self.baseline_correction_sizer_dim3.Add(
            self.baseline_correction_node_list_textcontrol_dim3,
            0,
            wx.ALIGN_CENTER_VERTICAL,
        )
        self.baseline_correction_sizer_dim3.AddSpacer(10)
        # If polynomial baseline correction is selected, have a textcontrol for the polynomial order

        self.baseline_correction_polynomial_order_label = wx.StaticText(
            parent, -1, "Polynomial order:"
        )
        self.baseline_correction_sizer_dim3.Add(
            self.baseline_correction_polynomial_order_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.baseline_correction_polynomial_order_textcontrol_dim3 = wx.TextCtrl(
            parent, -1, self.polynomial_order_dim3, size=(30, 20)
        )
        self.baseline_correction_polynomial_order_textcontrol_dim3.Bind(
            wx.EVT_TEXT, self.on_baseline_correction_textcontrol_dim3
        )
        self.baseline_correction_sizer_dim3.Add(
            self.baseline_correction_polynomial_order_textcontrol_dim3,
            0,
            wx.ALIGN_CENTER_VERTICAL,
        )
        self.baseline_correction_sizer_dim3.AddSpacer(10)

        if self.baseline_correction_radio_box_selection_dim3 == 0:
            self.baseline_correction_polynomial_order_label.Hide()
            self.baseline_correction_polynomial_order_textcontrol_dim3.Hide()

        # Have a button showing information on baseline correction
        self.baseline_correction_info = wx.Button(parent, -1, "\u24d8", size=(25, 32))

        self.baseline_correction_info.Bind(
            wx.EVT_BUTTON, self.oneDFrame.on_baseline_correction_info
        )
        self.baseline_correction_sizer_dim3.Add(
            self.baseline_correction_info, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.sizer_2.Add(self.baseline_correction_sizer_dim3)
        self.sizer_2.AddSpacer(10)

    def on_baseline_correction_checkbox_dim3(self, event):
        self.baseline_correction_checkbox_value_dim3 = (
            self.baseline_correction_checkbox_dim3.GetValue()
        )

    def on_baseline_correction_radio_box_dim3(self, event):
        # If the user selects linear or polynomial baseline correction, update the options
        self.baseline_correction_radio_box_selection_dim3 = (
            self.baseline_correction_radio_box_dim3.GetSelection()
        )

        if self.baseline_correction_radio_box_selection_dim3 == 0:
            # Remove the polynomial order textcontrol
            self.baseline_correction_sizer_dim3.Hide(
                self.baseline_correction_polynomial_order_label
            )
            self.baseline_correction_sizer_dim3.Hide(
                self.baseline_correction_polynomial_order_textcontrol_dim3
            )
            self.baseline_correction_sizer_dim3.Layout()
        elif self.baseline_correction_radio_box_selection_dim3 == 1:
            # Add the polynomial order textcontrol
            self.baseline_correction_sizer_dim3.Show(
                self.baseline_correction_polynomial_order_label
            )
            self.baseline_correction_sizer_dim3.Show(
                self.baseline_correction_polynomial_order_textcontrol_dim3
            )
            self.baseline_correction_sizer_dim3.Layout()

    def on_baseline_correction_textcontrol_dim3(self, event):
        # If the node width or node list textcontrols are changed, update the node width and node list
        self.node_width_dim3 = (
            self.baseline_correction_nodes_textcontrol_dim3.GetValue()
        )
        self.node_list_dim3 = (
            self.baseline_correction_node_list_textcontrol_dim3.GetValue()
        )
        self.polynomial_order_dim3 = (
            self.baseline_correction_polynomial_order_textcontrol_dim3.GetValue()
        )
