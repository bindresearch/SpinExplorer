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


import sys
import wx
import numpy as np
import nmrglue as ng
from matplotlib.figure import Figure
from SpinExplorer.SpinView.config import reference_range_values, multiply_range_values, vertical_range_values
from SpinExplorer.SpinView.UI_objects.UI_tools import FloatSlider, PhasingSliderRange
from SpinExplorer.SpinView.Viewers.loading import GetData
from SpinExplorer.SpinView.Viewers.oned_view import OneDViewer
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


# Referencing window
class ReferenceSpectrum(wx.Frame):
    def __init__(self, main_frame, file):
        # Get the monitor size and set the window size to 85% of the monitor size
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 1.0 * self.monitorWidth
        self.height = 0.85 * self.monitorHeight
        self.display_index = wx.Display.GetFromWindow(main_frame)
        self.display_index_current = self.display_index
        self.reference_frame = wx.Frame.__init__(
            self,
            None,
            wx.ID_ANY,
            "Referencing",
            wx.DefaultPosition,
            size=(int(self.width), int(self.height)),
        )

        """
        Allow a user to select a file containing a reference 1D dataset
        containing DSS. This assumes that the reference is the same
        sample at the same temperature recorded with preferably the
        same o1 value.
        """

        self.main_frame = main_frame

        try:
            # Read the data
            self.nmrdata = GetData(self, file=file)
            if self.nmrdata.dim == 1:
                self.viewer = OneDViewer(parent=self, nmrdata=self.nmrdata)
            else:
                # The data is not a 1D reference spectrum
                dlg = wx.MessageDialog(
                        None,
                        "The selected spectrum was not a 1D spectrum. Please select a 1D reference spectrum and try again.",
                        "Error",
                        wx.OK | wx.ICON_INFORMATION,
                )
                self.Raise()
                self.SetFocus()
                result = dlg.ShowModal()
                self.Destroy()
                return

        except:
            # The data is not a 1D reference spectrum
            dlg = wx.MessageDialog(
                    None,
                    "The selected file could not be read correctly. Please select a 1D reference spectrum in nmrPipe format (.ft) and try again.",
                    "Error",
                    wx.OK | wx.ICON_INFORMATION,
            )
            self.Raise()
            self.SetFocus()
            result = dlg.ShowModal()
            self.Destroy()
            return


        dlg = wx.MessageDialog(
                None,
                "The carrier frequency (e.g. O1) of the reference spectrum and the loaded spectrum must be the same for the referencing to work correctly. Would you like to continue?",
                "Referencing",
                wx.YES_NO | wx.ICON_INFORMATION,
        )
        self.Raise()
        self.SetFocus()
        result = dlg.ShowModal()
        if(result == wx.ID_NO):
            self.Destroy()
            return

        self.set_initial_variables()
        self.create_canvas()

        # Resize to ensure that the canvas gets the correct DPI of the current display
        w, h = self.GetSize()
        self.SetSize(w + 1, h)
        self.SetSize(w, h)
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
            self.canvas.SetSize(
                (
                    self.width * 0.0104,
                    (self.height - self.sizer.GetMinSize()[1] - 100)
                    * 0.0104,
                )
            )
            self.fig.set_size_inches(
                self.width * 0.0104,
                (self.height - self.sizer.GetMinSize()[1] - 100)
                * 0.0104,
            )
        # Resize to ensure that the canvas gets the correct DPI of the current display
        w, h = self.GetSize()
        self.SetSize(w + 1, h)
        self.SetSize(w, h)
        self.UpdateFrame()
        event.Skip()

    def OnSizeFrame(self, event):
        # Get the new frame size
        self.width, self.height = self.GetSize()
        self.SetSize((self.width, self.height))
        self.canvas.SetSize(
            (
                self.width * 0.0104,
                (self.height - self.sizer.GetMinSize()[1] - 100)
                * 0.0104,
            )
        )
        self.fig.set_size_inches(
            self.width * 0.0104,
            (self.height - self.sizer.GetMinSize()[1] - 100) * 0.0104,
        )
        self.UpdateFrame()
        event.Skip()


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

        self.create_buttons()
        self.draw_figure_1D_referencing()
        self.Layout()
        self.Show()


    def set_initial_variables(self):
        self.reference_select = False
        self.reference_value = 0.0
        self.reference_selected_value = 0.0

    def create_buttons(self):

        # Create the phasing 1D sizer
        self.referencing_label = wx.StaticBox(self, -1, "Referencing options:")
        self.referencing_sizer = wx.StaticBoxSizer(self.referencing_label, wx.HORIZONTAL)

        # Have a toggle button allowing a user to add and update the location of a reference peak
        self.reference_toggle_button = wx.ToggleButton(self.referencing_label, label="Select reference point")
        self.reference_toggle_button.Bind(wx.EVT_TOGGLEBUTTON, self.OnUpdateReference)

        # Have a textcontrol where users can input the desired chemical shift of the reference point
        self.reference_text = wx.StaticText(self.referencing_label, label="Desired chemical shift of reference point (ppm):")
        self.reference_input = wx.TextCtrl(self.referencing_label, -1, value=str(self.reference_value), size=(100,20))

        # Have a button to apply the updated referencing
        self.apply_referencing_button = wx.Button(self.referencing_label, label='Apply referencing')
        self.apply_referencing_button.Bind(wx.EVT_BUTTON, self.OnApplyReferencing)


        self.referencing_sizer.AddSpacer(10)
        self.referencing_sizer.Add(self.reference_toggle_button)
        self.referencing_sizer.AddSpacer(10)
        self.referencing_sizer.Add(self.reference_text)
        self.referencing_sizer.AddSpacer(5)
        self.referencing_sizer.Add(self.reference_input)
        self.referencing_sizer.AddSpacer(10)
        self.referencing_sizer.Add(self.apply_referencing_button)
        self.referencing_sizer.AddSpacer(10)


        # Create a sizer for changing the y axis limits in the spectrum
        self.zoom_label = wx.StaticBox(self, -1, "Y Axis Zoom (%):")
        self.zoom_sizer = wx.StaticBoxSizer(self.zoom_label, wx.VERTICAL)
        self.intensity_slider = FloatSlider(
            self.zoom_label, id=-1, value=0, minval=-1, maxval=10, res=0.01, size=(300, height)
        )
        self.intensity_slider.Bind(wx.EVT_SLIDER, self.OnIntensityScroll1D)
        self.zoom_sizer.AddSpacer(5)
        self.zoom_sizer.Add(self.intensity_slider)

        # Add all the sizers to the main sizer
        self.sizer1 = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer1.Add(self.referencing_sizer, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer1.AddSpacer(20)
        self.sizer1.Add(self.zoom_sizer, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer.Add(self.sizer1, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.sizer.AddSpacer(10)

    def UpdateFrame(self):
        self.canvas.draw()
        self.canvas.Refresh()
        self.canvas.Update()


    def OnUpdateReference(self, event):
        if(self.reference_select==False):
            self.reference_connect = self.fig.canvas.mpl_connect(
                    "button_press_event", self.on_click_reference
                )
            self.reference_toggle_button.SetValue(True)
            self.reference_select=True
            self.reference_line.set_visible(True)
            self.UpdateFrame()
        else:
            self.reference_connect = self.fig.canvas.mpl_connect(
                    "button_press_event", self.on_click_reference
                )
            self.fig.canvas.mpl_disconnect(self.reference_connect)
            self.reference_toggle_button.SetValue(False)
            self.reference_select=False
            self.reference_line.set_visible(False)
            self.UpdateFrame()
        
    def on_click_reference(self, event):
        """
        Edit the value of the self.reference_position to
        the selected position and 
        """

        x, y = self.ax.transData.inverted().transform((event.x, event.y))
        if(x==None):
            return
        self.reference_selected_value = float(x)
        self.reference_line.set_xdata([self.reference_selected_value])
        self.UpdateFrame()



    def OnApplyReferencing(self, event):
        """
        Find the central chemical shift of the data after setting the chemical shift
        of the selected reference point to the desired reference value
        """
        # Previous carrier
        previous_carrier = float(self.dic['FDF2CAR'])

        reference_correction = self.reference_selected_value - self.reference_value

        new_carrier = previous_carrier - reference_correction

        self.main_frame.shared_format.on_save_reference(new_carrier)

        self.Destroy()



    def OnIntensityScroll1D(self, event):
        # Function to change the y axis limits
        intensity_percent = 10 ** float(self.intensity_slider.GetValue())
        self.ax.set_ylim(
            -(np.max(self.data) / 8) / (intensity_percent / 100),
            np.max(self.data) / (intensity_percent / 100),
        )
        self.UpdateFrame()

    def on_mouse_wheel(self, event):
        
        toolbar = self.fig.canvas.toolbar
        if toolbar:
            toolbar.push_current() # logs position in toolbar so commands back, forward, home work
        mx, my = event.GetPosition()

        scale = self.fig.canvas.GetDPIScaleFactor()
        mx *= scale
        my *= scale

        h = self.fig.canvas.GetSize().height * scale
        my = h - my

        zoom = 1.1 if event.GetWheelRotation() < 0 else 1/1.1

        renderer = self.fig.canvas.get_renderer()

        for ax in self.fig.axes:
            bbox = ax.get_window_extent(renderer=renderer)

            if not bbox.contains(mx, my):
                continue

            inv = ax.transData.inverted()
            x, y = inv.transform((mx, my))

            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

            ax.set_xlim([x + (v - x) * zoom for v in xlim])
            ax.set_ylim([y + (v - y) * zoom for v in ylim])

        self.fig.canvas.draw_idle()

    def on_key_1d(self, event):
        # navigator options
        if event.key == "z":
            self.toolbar.zoom()
        if event.key == "p":
            self.toolbar.pan()
        if event.key == "q":
            self.toolbar.home()
        if event.key == "b":
            self.toolbar.back()
        if event.key == "f":
            self.toolbar.forward()



    def draw_figure_1D_referencing(self):
        # Function to plot the 1D spectrum
        self.ax = self.fig.add_subplot(111)

        self.key_press_connect = self.fig.canvas.mpl_connect(
            "key_press_event", self.on_key_1d
        )
        self.mouse_wheel_connect = self.fig.canvas.Bind(wx.EVT_MOUSEWHEEL, self.on_mouse_wheel)

        self.data = self.nmrdata.data
        self.dic = self.nmrdata.dic
        self.ppms = ng.pipe.make_uc(self.dic, self.data, dim=-1).ppm_scale()

        (self.line1,) = self.ax.plot(self.ppms, self.data, linewidth=0.5)
        self.ax.set_xlabel("Chemical shift (ppm)")
        self.ax.set_ylabel("Intensity")
        self.ax.set_xlim(max(self.ppms), min(self.ppms))
        self.line1.set_color("tab:blue")
        self.reference_line = self.ax.axvline(
            self.reference_selected_value, color="black", linestyle="--"
        )
        self.reference_line.set_visible(False)
        self.UpdateFrame()

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