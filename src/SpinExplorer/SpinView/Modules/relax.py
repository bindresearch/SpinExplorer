import wx # type: ignore
import numpy as np
import sys
import os
import json
import math
import matplotlib 
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_wxagg import (
    NavigationToolbar2WxAgg as NavigationToolbar,
)
from scipy.optimize import leastsq # type: ignore

from SpinExplorer.SpinView.UI_objects.UI_tools import FloatSlider
from SpinExplorer.SpinView.Viewers.overlays import DeleteSliceDialog
from SpinExplorer.SpinView.Viewers.module_utils import DelaysManualInput, InputROI

if sys.platform == "linux":
    platform = "linux"
    height = 30
elif sys.platform == "darwin":
    platform = "mac"
    height = 16
else:
    platform = "windows"
    height = 30


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
        self.title = title
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
        self.relax_sizer = wx.BoxSizer(wx.VERTICAL)
        
        self.row1 = wx.BoxSizer(wx.HORIZONTAL)
        self.row2 = wx.BoxSizer(wx.HORIZONTAL)

        # Create a button that opens a file for a user to input the delay times
        self.delays_label = wx.StaticBox(self, -1, "Delay Times")
        self.delays_sizer_total = wx.StaticBoxSizer(self.delays_label, wx.VERTICAL)
        self.delays_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.delay_times_button = wx.Button(self.delays_label, -1, "Input Delay Times")
        self.delay_times_button.Bind(wx.EVT_BUTTON, self.OnInputDelayTimes)
        self.delays_sizer.AddSpacer(5)
        self.delays_sizer.Add(self.delay_times_button)
        self.delays_sizer.AddSpacer(5)
        
        self.delays_sizer_total.AddSpacer(4)
        self.delays_sizer_total.Add(self.delays_sizer)
        self.delays_sizer_total.AddSpacer(3)


        # Then have button which will allow a user to drag over a section where they wish to estimate the noise level
        # This can then be plotted as a shaded region on the plot
        self.noise_label = wx.StaticBox(self, -1, "Noise Region")
        self.noise_sizer_total = wx.StaticBoxSizer(self.noise_label, wx.VERTICAL)
        self.noise_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.select_noise_button = wx.Button(self.noise_label, -1, "Select Noise Region")
        self.select_noise_button.Bind(wx.EVT_BUTTON, self.OnSelectNoise)
        self.noise_sizer.AddSpacer(5)
        self.noise_sizer.Add(self.select_noise_button)
        self.noise_sizer.AddSpacer(5)

        # Then have a TextCtrl for the minimum SNR for the relaxation coefficient to be estimated (default to 10)
        self.noise_factor = 10
        self.noise_factor_label = wx.StaticText(self.noise_label, -1, "Minimum SNR:")
        self.noise_factor_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.noise_sizer.AddSpacer(10)
        self.noise_sizer.Add(self.noise_factor_label)
        self.noise_sizer.AddSpacer(5)
        self.noise_factor_box = wx.TextCtrl(
            self.noise_label, -1, str(self.noise_factor), size=(50, -1)
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
            self.fitting_type_label, -1, choices=self.choices, style=wx.RA_HORIZONTAL
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
        self.whole_spectrum_fitting_button = wx.Button(self.fitting_label, -1, "Fit Whole Spectrum")
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
        self.add_region_button = wx.Button(self.fitting_label, -1, "Add ROI")
        self.add_region_button.Bind(wx.EVT_BUTTON, self.OnAddROI)
        self.fitting_sizer.Add(self.add_region_button)
        self.fitting_sizer.AddSpacer(5)

        self.input_region_button = wx.Button(self.fitting_label, -1, "Input ROI")
        self.input_region_button.Bind(wx.EVT_BUTTON, self.OnInputROI)
        self.fitting_sizer.Add(self.input_region_button)
        self.fitting_sizer.AddSpacer(5)

        # Have a button to delete a region of interest
        self.delete_region_button = wx.Button(self.fitting_label, -1, "Delete ROI")
        self.delete_region_button.Bind(wx.EVT_BUTTON, self.OnDeleteROI)
        self.fitting_sizer.Add(self.delete_region_button)
        self.fitting_sizer.AddSpacer(5)

        # Can then have a button which will fit the Relaxation equation to the mean values of the points above the noise in the region of interest
        self.region_fitting_button = wx.Button(self.fitting_label, -1, "Fit")
        self.region_fitting_button.Bind(wx.EVT_BUTTON, self.OnRegionFitting)
        self.fitting_sizer.Add(self.region_fitting_button)
        self.fitting_sizer.AddSpacer(5)

        # Have a button which will perform a biexponential fit on the data in the region of interest
        self.biexponential_fitting_button = wx.Button(self.fitting_label, -1, "Biexponential Fit")
        self.biexponential_fitting_button.Bind(
            wx.EVT_BUTTON, self.OnBiexponentialFitting
        )
        self.fitting_sizer.Add(self.biexponential_fitting_button)
        self.fitting_sizer.AddSpacer(5)

        self.save_fitting_button = wx.Button(self.fitting_label, -1, "Save Fit")
        self.save_fitting_button.Bind(wx.EVT_BUTTON, self.OnSaveFitting)
        self.fitting_sizer.Add(self.save_fitting_button)
        self.fitting_sizer.AddSpacer(5)

        self.fitting_sizer_total.AddSpacer(4)
        self.fitting_sizer_total.Add(self.fitting_sizer)
        self.fitting_sizer_total.AddSpacer(4)

        # # Have a button for printing fitted values
        # self.print_fitted_values_button = wx.Button(self, -1, "Print Fit")
        # self.print_fitted_values_button.Bind(wx.EVT_BUTTON, self.OnPrintFittedValues)
        # self.fitting_sizer.Add(self.print_fitted_values_button)
        # self.fitting_sizer.AddSpacer(5)

        # Have a box containing other functions such as a button to delete a slice from the plot and repeat the fitting
        self.other_functions_label = wx.StaticBox(self, -1, "Other Functions")
        self.other_functions_sizer = wx.StaticBoxSizer(
            self.other_functions_label, wx.VERTICAL
        )
        self.other_functions_sizer.AddSpacer(4)
        self.delete_slice_button = wx.Button(self.other_functions_label, -1, "Delete Slice")
        self.delete_slice_button.Bind(wx.EVT_BUTTON, self.OnDeleteSlice)
        self.other_functions_sizer.Add(self.delete_slice_button)
        self.other_functions_sizer.AddSpacer(4)

        # Creating a sizer for changing the y axis limits in the spectrum
        self.intensity_label = wx.StaticBox(self, -1, "Y Axis Zoom (%):")
        self.intensity_sizer = wx.StaticBoxSizer(self.intensity_label, wx.VERTICAL)
        width = 100
        self.intensity_slider = FloatSlider(
            self.intensity_label, id=-1, value=0, minval=-1, maxval=10, res=0.01, size=(width, height)
        )
        self.intensity_slider.Bind(wx.EVT_SLIDER, self.OnIntensityScrollRelax)
        self.intensity_sizer.AddSpacer(5)
        self.intensity_sizer.Add(self.intensity_slider)

        
        self.row1.AddSpacer(5)
        self.row1.Add(self.delays_sizer_total)
        self.row1.AddSpacer(5)
        self.row1.Add(self.noise_sizer_total)
        self.row1.AddSpacer(5)
        self.row1.Add(self.fitting_type_sizer)
        self.row1.AddSpacer(5)
        self.row1.Add(self.other_functions_sizer)
        self.row1.AddSpacer(5)
        self.row1.Add(self.intensity_sizer)

        self.row2.Add(self.fitting_sizer_total)
        
        self.relax_sizer.Add(self.row1, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.relax_sizer.AddSpacer(5)
        self.relax_sizer.Add(self.row2, 0, wx.ALIGN_CENTER_HORIZONTAL)


        self.main_relax_sizer.AddSpacer(5)
        self.main_relax_sizer.Add(self.relax_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
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
        """
        Replaced by save fit button
        """
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
        legend = self.ax_relax.legend(title="Slice Number", ncol=math.ceil(len(self.slice_plots)/8))
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


        self.press = False
        self.move = False

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

        self.check_noise = self.check_noise_points_relax()

        # If self.check_noise = False, the noise region is too small and covers less than 2 points
        if(self.check_noise == False):
            # Give an error message saying noise region has not been selected
            msg = wx.MessageDialog(
                self,
                "The selected noise region does not cover enough points. Please re-select the noise region covering more points and try again.",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return

        # Turn on noise region selection flag
        self.noise_region_selection = True

    def check_noise_points_relax(self):
        """
        Checking that the noise regions selected contains data points
        """
        noise_region = self.y_data[
            :,
            np.where(
                (self.x_data >= self.noise_x_initial)
                & (self.x_data <= self.noise_x_final)
            )[0],
        ][0]

        if(len(noise_region)<2):
            # There are less than 2 points in the current selected noise region
            return False

        return True


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
        legend = self.ax_relax.legend(title="Slice Number", ncol=math.ceil(len(self.slice_plots)/8))
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
    
    def OnInputROI(self, event):
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
        


        # Getting the user to input the new ROI values
        self.user_input_region()


    def user_input_region(self):
        """
        Have a poput mini window where a user can insert their desired chemical shift range
        for the region of interest
        """

        input_roi = InputROI(title='Input ROI', parent=self)


    def add_user_input_region(self, xmin, xmax):
        # Adding an ROI region

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
        

        self.ROI_regions[-1].set_x(xmin)
        self.ROI_regions[-1].set_width(xmax - xmin)
        self.ROI_regions_2[-1].set_x(xmin)
        self.ROI_regions_2[-1].set_width(xmax - xmin)
        self.ROI_regions_3[-1].set_x(xmin)
        self.ROI_regions_3[-1].set_width(xmax - xmin)
        # self.ROI_regions[-1].set_xy([[xmin,0],[xmin,1],[xmax,1],[xmax,0]])
        # self.ROI_regions_2[-1].set_xy([[xmin,0],[xmin,1],[xmax,1],[xmax,0]])
        # self.ROI_regions_3[-1].set_xy([[xmin,0],[xmin,1],[xmax,1],[xmax,0]])
        # Check that the number of points in the selected region of interest is greater than 2
        check_ROI = self.check_ROI(xmin, xmax)
        if(check_ROI == True):
            # Add the min and max ppm values to the array of selected regions of interest
            self.selected_regions_of_interest.append([xmin, xmax])
            self.UpdateRelaxFrame()
            return True
        else:
            # Deleting this selected ROI because it contained less than 2 points
            self.DeleteSmallROI()
            self.UpdateRelaxFrame()
            return False

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

    def DeleteSmallROI(self):
        """
        Deleting a ROI which did not contain enough points when selected by a user
        """
        del self.ROI_regions[-1]
        del self.ROI_regions_2[-1]
        del self.ROI_regions_3[-1]
        del self.ROI_color[-1]
        self.deleted_ROI_number += 1
        self.PlotWholeSpectrumFitting()

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
                    self.PlotWholeSpectrumFitting()
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
                    self.PlotWholeSpectrumFitting()

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

            # Check that the number of points in the selected region of interest is greater than 2
            check_ROI = self.check_ROI(xmin, xmax)

            if(check_ROI == True):
                # Add the min and max ppm values to the array of selected regions of interest
                self.selected_regions_of_interest.append([xmin, xmax])
            else:
                # Deleting this selected ROI because it contained less than 2 points
                self.DeleteSmallROI()


            self.UpdateRelaxFrame()

            self.press = False
            self.move = False
            self.canvas_relax.mpl_disconnect(self.select_ROI_press)
            self.canvas_relax.mpl_disconnect(self.select_ROI_move)
            self.canvas_relax.mpl_disconnect(self.select_ROI_release)

    def check_ROI(self, xmin, xmax):
        """
        Checking that the ROI selected contains data points
        """

        ppms_in_ROI = []
        for i, ppm in enumerate(self.ppms_above_noise):
            if ppm >= xmin and ppm <= xmax:
                ppms_in_ROI.append(ppm)

        if(len(ppms_in_ROI)<2):
            # There are less than 2 points in the current selected noise region

            # Give an error message saying noise region has not been selected
            msg = wx.MessageDialog(
                self,
                "The selected region of interest (ROI) does not cover enough points. Please re-select the ROI covering more points and try again.",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()

            return False

        return True

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
        legend = self.ax_relax.legend(title="Slice Number", fontsize=8, ncol=math.ceil(len(self.slice_plots)/8))
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
                    color=self.ROI_color[i],
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
                    + str(round(sigmas[i], 2)) + r' $s^{-1}$',
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
                    + str(round(sigmas[i], 2))  + r' $s^{-1}$',
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
        


    def OnSaveFitting(self, event):
        """
        Save each of the current fits along with the data to a csv file
        (include ppm ranges used and the data title and path etc)

        Metadata JSON 
        - Data title
        - Global fit noise region (ppm)
        - Global fit noise value
        - Global fit minimum S/N
        
        Region of interest csv files 
        - Delays, Intensity, Intensity error, relaxation rates

        """

        # Make sure that there are regions of interest that have fits performed on them
        check = self.check_regions()

        if(check==True):
            # Ask the user to provide a path of a directory to save the fit data and metadata to
            self.ask_user_path()
        else:
            # Outputting error message saving that no regions of interest could be found
            dlg = wx.MessageDialog(
            self,
            "No regions of interest (ROI) could be found. Please add a region of interest, perform a fit, and try again.",
            "Error",
            wx.OK,
            )
            self.Raise()
            self.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
        
        return

    def check_regions(self):
        """
        Check that regions of interest have been selected and that they have fits associated
        with them

        Returns
        -------
        True - check passed
        False - check failed
        """
        if(len(self.selected_regions_of_interest)>0):
            return True
        else:
            return False
        
    def ask_user_path(self):
        """
        Ask the user where they would like to save the folder which will contain the outputs
        """

        dlg = wx.FileDialog (self, "Input directory name to create and save output files to",
            style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
        )
        dlg.SetDirectory(os.getcwd())
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            os.makedirs(path, exist_ok=True)
            self.save_metadata(path)
            self.save_data(path)
            dlg.Destroy()
        else:
            return

    def save_metadata(self, path):
        """
        Saving the fit metadata to the path\metadata.json
        """

        noise_region_text = str(self.noise_x_initial) + '-' + str(self.noise_x_final)

        metadata = {"Data title": self.title, 'Global fitting noise region (ppm)': noise_region_text, 'Global fit minimum S/N': self.noise_factor}
        
        # Save the metadata
        with open(os.path.join(path, 'metadata.json'), 'w') as f:
            json.dump(metadata, f)


    def save_data(self, path):
        """
        For each region of interest, saving the fit as a csv
        """
        
        for i, region in enumerate(self.selected_regions_of_interest):
            self.region_min = region[0]
            self.region_max = region[1]
            self.mean_fitted_relax_ROI = self.mean_fitted_relax_ROI_total[i]
            self.mean_fitted_I0_ROI = self.mean_fitted_I0_ROI_total[i]


            self.mean_fitted_I0_ROI, self.mean_fitted_relax_ROI


            if (self.fitted_gaussian_parameters[i][2]) > 1.25 * np.std(
                self.fitted_relax_ROI_total[i]
            ):
                r_error = np.abs(np.std(self.fitted_relax_ROI_total[i]))

            else:
                r_error = np.abs(self.fitted_gaussian_parameters[i][2])

            I0_error = np.abs(np.std(self.fitted_I0_ROI))

            self.average_y_data_in_ROI_above_noise = self.average_y_data_in_ROI_above_noise_total[i]
            
            self.error_y_data_in_ROI_above_noise = self.error_y_data_in_ROI_above_noise_total[i]

            self.error_I_I0_in_ROI = self.error_I_I0_in_ROI_total[i]
            self.I0_average_in_ROI = self.I0_average_in_ROI_total[i]



            # The data will be saved in the following manner

            # Fit results
            # region_min =
            # region_max = 
            # R2/R1 (s^-1) =
            # I0 =  

            # Delay (ms), I, I error, I/I0, I/I0 error


            with open(os.path.join(path, 'Fit_ROI_'+str(i+1)+'.csv'), 'w') as file:
                file.write('Region min (ppm),'+ str(self.region_min)+'\n')
                file.write('Region max (ppm),'+ str(self.region_max)+'\n')
                if(self.R1_fit==True):
                    relax_name = 'R1'
                else:
                    relax_name = 'R2'
                file.write(relax_name+' (s^-1),{:.3e}\n'.format(self.mean_fitted_relax_ROI))
                file.write(relax_name+' error (s^-1),{:.3e}\n'.format(r_error))
                file.write('I0,{:.3e}\n'.format(self.mean_fitted_I0_ROI))
                file.write('I0 error,{:.3e}\n'.format(I0_error))

                file.write('\n\n')


                file.write('Delay (s), I, I error, I/I0, I/I0 error\n')
                for j in range(len(self.delays)):
                    file.write('{:.3f},{:.3f},{:.3e},{:.3e},{:.3f}\n'.format(self.delays[j], self.average_y_data_in_ROI_above_noise[j], self.error_y_data_in_ROI_above_noise[j], self.average_y_data_in_ROI_above_noise[j] / self.I0_average_in_ROI[j], self.error_I_I0_in_ROI[j]))



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
