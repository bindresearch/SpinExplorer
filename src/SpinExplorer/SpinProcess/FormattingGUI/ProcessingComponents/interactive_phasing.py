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

import sys
import wx
import numpy as np
import nmrglue as ng
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import (
    NavigationToolbar2WxAgg as NavigationToolbar,
)

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

# Interactive phasing window


class InteractivePhasingFrame(wx.Frame):
    def __init__(self, main_frame, nmr_spectrum, ppms, nmr_d):
        # Get the monitor size and set the window size to 85% of the monitor size
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 1.0 * self.monitorWidth
        self.height = 0.85 * self.monitorHeight
        self.phasing_frame = wx.Frame.__init__(
            self,
            None,
            wx.ID_ANY,
            "Interactive Phasing",
            wx.DefaultPosition,
            size=(int(self.width), int(self.height)),
        )

        self.main_frame = main_frame

        self.nmr_spectrum = nmr_spectrum
        self.ppms = ppms

        self.total_P0 = 0.0
        self.total_P1 = 0.0

        try:
            if len(self.nmr_spectrum[0]) > 1:
                try:
                    len(self.nmr_spectrum[0][0])
                    self.nmr_spectrum = self.nmr_spectrum[0][0]
                except:
                    max_value = 0
                    max_index = 0
                    for index, slice in enumerate(self.nmr_spectrum):
                        if np.max(np.abs(slice)) > max_value:
                            max_value = np.max(np.abs(slice))
                            max_index = index
                    self.nmr_spectrum = self.nmr_spectrum[max_index]

        except:
            pass

        self.nmr_d = nmr_d

        self.create_canvas()

    def create_canvas(self):

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.sizer)
        self.fig = Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.toolbar = NavigationToolbar(self.canvas)

        self.sizer.Add(self.toolbar, 0, wx.EXPAND)
        self.sizer.AddSpacer(10)

        # Suppress complex warning from numpy
        import warnings

        # warnings.simplefilter("ignore", np.ComplexWarning)  # For old numpy versions
        warnings.simplefilter(
            "ignore", np.exceptions.ComplexWarning
        )  # For new numpy versions

        self.create_sliders()
        self.draw_figure_1D_phasing()
        self.Layout()
        self.Show()

    def create_sliders(self):
        from SpinExplorer.SpinView.SpinView import FloatSlider

        # Create the phasing 1D sizer
        self.phasing_label = wx.StaticBox(self, -1, "Phasing:")
        self.phasing_sizer = wx.StaticBoxSizer(self.phasing_label, wx.VERTICAL)
        self.P0_label = wx.StaticText(self, label="P0 (Coarse):")
        self.P1_label = wx.StaticText(self, label="P1 (Coarse):")
        self.P0_slider = FloatSlider(
            self, id=-1, value=0.0, minval=-180, maxval=180, res=0.1, size=(300, height)
        )
        self.P1_slider = FloatSlider(
            self, id=-1, value=0.0, minval=-180, maxval=180, res=0.1, size=(300, height)
        )
        self.P0_slider.Bind(wx.EVT_SLIDER, self.OnSliderScroll1D)
        self.P1_slider.Bind(wx.EVT_SLIDER, self.OnSliderScroll1D)
        self.P0_label_fine = wx.StaticText(self, label="P0 (Fine):     ")
        self.P1_label_fine = wx.StaticText(self, label="P1 (Fine):     ")
        self.P0_slider_fine = FloatSlider(
            self, id=-1, value=0.0, minval=-10, maxval=10, res=0.01, size=(300, height)
        )
        self.P1_slider_fine = FloatSlider(
            self, id=-1, value=0.0, minval=-10, maxval=10, res=0.01, size=(300, height)
        )
        self.P0_slider_fine.Bind(wx.EVT_SLIDER, self.OnSliderScroll1D)
        self.P1_slider_fine.Bind(wx.EVT_SLIDER, self.OnSliderScroll1D)
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
        self.phasing_combined.AddSpacer(160)
        self.phasing_combined.Add(self.P0_total_value)
        self.phasing_combined.AddSpacer(170)
        self.phasing_combined.Add(self.P1_total)
        self.phasing_combined.AddSpacer(160)
        self.phasing_combined.Add(self.P1_total_value)
        self.phasing_sizer.AddSpacer(5)
        self.phasing_sizer.Add(self.sizer_coarse)
        self.phasing_sizer.AddSpacer(10)
        self.phasing_sizer.Add(self.sizer_fine)
        self.phasing_sizer.AddSpacer(10)
        self.phasing_sizer.Add(self.phasing_combined)

        # Add a button to set the pivot point for phasing
        self.pivot_button = wx.Button(self, label="Set Pivot Point")
        self.pivot_button.Bind(wx.EVT_BUTTON, self.OnPivotButton)
        self.pivot_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.pivot_sizer.AddSpacer(500)
        self.pivot_sizer.Add(self.pivot_button)

        self.pivot_x_default = 0

        # Add a button to remove the pivot point
        self.remove_pivot_button = wx.Button(self, label="Remove Pivot Point")
        self.remove_pivot_button.Bind(wx.EVT_BUTTON, self.OnRemovePivotButton)
        self.pivot_sizer.AddSpacer(20)
        self.pivot_sizer.Add(self.remove_pivot_button)

        # Add the pivot point buttons to the phasing sizer
        self.phasing_sizer.AddSpacer(10)
        self.phasing_sizer.Add(self.pivot_sizer)

        # Create a sizer for changing the y axis limits in the spectrum
        self.zoom_label = wx.StaticBox(self, -1, "Y Axis Zoom (%):")
        self.zoom_sizer = wx.StaticBoxSizer(self.zoom_label, wx.VERTICAL)
        self.intensity_slider = FloatSlider(
            self, id=-1, value=0, minval=-1, maxval=10, res=0.01, size=(300, height)
        )
        self.intensity_slider.Bind(wx.EVT_SLIDER, self.OnIntensityScroll1D)
        self.zoom_sizer.AddSpacer(5)
        self.zoom_sizer.Add(self.intensity_slider)

        # Have a save and close button
        self.save_button = wx.Button(self, label="Save and Close")
        self.save_button.Bind(wx.EVT_BUTTON, self.OnSavePhasing)

        # Add all the sizers to the main sizer
        self.sizer1 = wx.BoxSizer(wx.VERTICAL)
        self.sizer2 = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer2.Add(self.phasing_sizer, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer2.AddSpacer(20)
        self.sizer2.Add(self.zoom_sizer, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer1.Add(self.sizer2, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.sizer1.AddSpacer(20)
        self.sizer1.Add(self.save_button, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.sizer1.AddSpacer(20)
        self.sizer.Add(self.sizer1, 0, wx.ALIGN_CENTER_HORIZONTAL)

    def UpdateFrame(self):
        self.canvas.draw()
        self.canvas.Refresh()
        self.canvas.Update()

    def OnSavePhasing(self, event):
        # Function to save the phasing values and close the window
        self.main_frame.phase_correction_p0_textcontrol.SetValue(str(self.total_P0))
        self.main_frame.phase_correction_p1_textcontrol.SetValue(str(self.total_P1))
        self.Close()

    def OnIntensityScroll1D(self, event):
        # Function to change the y axis limits
        intensity_percent = 10 ** float(self.intensity_slider.GetValue())
        self.ax.set_ylim(
            -(np.max(self.data) / 8) / (intensity_percent / 100),
            np.max(self.data) / (intensity_percent / 100),
        )
        self.UpdateFrame()

    def draw_figure_1D_phasing(self):
        # Function to plot the 1D spectrum
        self.ax = self.fig.add_subplot(111)

        self.data = self.nmr_spectrum

        (self.line1,) = self.ax.plot(self.ppms, self.data, linewidth=0.5)
        self.pivot_x = self.pivot_x_default
        self.ax.set_xlabel("Chemical shift (ppm)")
        self.ax.set_ylabel("Intensity")
        self.ax.set_xlim(max(self.ppms), min(self.ppms))
        self.line1.set_color("tab:blue")
        self.pivot_line = self.ax.axvline(
            self.pivot_x_default, color="black", linestyle="--"
        )
        self.pivot_line.set_visible(False)
        self.UpdateFrame()

    def OnPivotButton(self, event):
        # Get the user to select a pivot point for phasing by clicking on the spectrum
        # Give a message box to tell the user to click on the spectrum where they want the pivot point to be
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

        # Find the index of the point closest to the pivot point
        self.pivot_index = np.abs(self.ppms - self.pivot_x).argmin()
        self.pivot_x = self.pivot_index
        self.canvas.mpl_disconnect(self.pivot_press)
        self.pivot_line.set_visible(True)
        self.OnSliderScroll1D(wx.EVT_SCROLL)

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

    def OnSliderScroll1D(self, event):
        # Get all the slider values for P0 and P1 (coarse and fine), put the combined coarse and fine values on the screen
        self.total_P0 = self.P0_slider.GetValue() + self.P0_slider_fine.GetValue()
        self.total_P1 = self.P1_slider.GetValue() + self.P1_slider_fine.GetValue()
        self.P0_total_value.SetLabel("{:.2f}".format(self.total_P0))
        self.P1_total_value.SetLabel("{:.2f}".format(self.total_P1))
        self.phase1D()

    def phase1D(self):
        # Function to phase the data using the combined course/fine phasing values and plot
        imaginary_data = ng.process.proc_base.ht(
            self.nmr_spectrum, self.nmr_spectrum.shape[0]
        )
        self.data = imaginary_data * np.exp(
            1j
            * (
                self.total_P0 * np.pi / 180
                + self.total_P1
                * (np.pi / 180)
                * (
                    np.arange(-self.pivot_x, -self.pivot_x + self.nmr_spectrum.shape[0])
                    / self.nmr_spectrum.shape[0]
                )
            )
        ) + np.ones(len(self.nmr_spectrum))
        self.line1.set_ydata(self.data + np.ones(len(self.data)))
        self.UpdateFrame()
