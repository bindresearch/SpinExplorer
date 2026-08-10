import wx # type: ignore
import numpy as np 
import nmrglue as ng # type: ignore
import sys
import os
import math
import json
import matplotlib
matplotlib.use("wxAgg")
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import (
    NavigationToolbar2WxAgg as NavigationToolbar,
)
import matplotlib.gridspec as gridspec
from scipy.optimize import leastsq # type: ignore
from SpinExplorer.SpinView.UI_objects.UI_tools import FloatSlider
from SpinExplorer.SpinView.Viewers.overlays import DeleteSliceDialog
from SpinExplorer.SpinView.Viewers.module_utils import DiffusionGradientManualInput, InputROI

if sys.platform == "linux":
    platform = "linux"
    height = 30
elif sys.platform == "darwin":
    platform = "mac"
    height = 16
else:
    platform = "windows"
    height = 30

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
        self.title=title
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
        self.gradient_integral_factor = 1.0
        self.max_gradient = 53.0
        self.DAC_conversion = 0.002
        self.gamma_dictionary = {}
        self.gamma_dictionary["1H"] = 2.67522e4  # rad s-1 G-1
        self.gamma_dictionary["2H"] = 0.41065e4  # rad s-1 G-1
        self.gamma_dictionary["13C"] = 0.672828e4  # rad s-1 G-1
        self.gamma_dictionary["15N"] = -0.27116e4  # rad s-1 G-1
        self.gamma_dictionary["19F"] = 2.51815e4  # rad s-1 G-1
        self.gamma_dictionary["31P"]= 1.08291e4 # rad s-1 G-1
        self.gamma_dictionary["23Na"] = 0.70761e4 # rad s-1 G-1
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

        self.read_parameters()


    def read_parameters(self):
        """
        Checking to see if pulseprogram is in the current directory. If it is,
        check to see if the gradients are bipolar gradients, what the gradient shape
        is (for the gradient integral factor) and what the little_delta and big_delta
        values are. Then reading difflist and difframp and use the gradient integral
        factor to work out Gmax in G/cm
        """
        try:
            from pathlib import Path
            self.folder = Path.cwd()
    
            pulseprogram_path = self.folder / 'pulseprogram'
            
            if not pulseprogram_path.exists():
                pulseprogram_path = self.folder / 'pulseprogram.precomp'

            with open(pulseprogram_path, 'r') as file:
                sequence_line = file.readlines()[0]
                if('stebpesgp1s' in sequence_line or 'stebpgp1s' in sequence_line):
                    self.bipolar_gradients = True

            # Update little and big delta values based on if bipolar gradients or not
            self.find_parameters(wx.EVT_BUTTON, initial_find=True)

            # Find out GP6 gradient shape
            with open(self.folder / 'acqus', 'r') as acqus_file:
                lines = acqus_file.readlines()
                found_GPNAM = False
                for line in lines:
                    if(found_GPNAM==True):
                        gradient = line.split()[5]
                        if(gradient == '<SMSQ10.100>'):
                            self.gradient_integral_factor = 0.9
                        break

                    if('GPNAM' in line):
                        found_GPNAM = True



            self.find_max_gradient_estimate()


        except:
            pass

            



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
            self.diffusion_data_type_label, -1, choices=["Bruker", "Varian"], style=wx.RA_HORIZONTAL
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
        self.bipolar_gradients_checkbox = wx.CheckBox(self.experimental_parameters_label, label="Bipolar Gradients")
        self.bipolar_gradients_checkbox.SetValue(self.bipolar_gradients)
        self.bipolar_gradients_checkbox.Bind(wx.EVT_CHECKBOX, self.OnBipolarGradients)
        self.experimental_parameters_sizer.Add(self.bipolar_gradients_checkbox)
        self.little_delta_label = wx.StaticText(self.experimental_parameters_label, -1, "δ (μs):")
        self.experimental_parameters_sizer.AddSpacer(5)
        self.experimental_parameters_sizer.Add(self.little_delta_label)
        self.experimental_parameters_sizer.AddSpacer(5)
        self.little_delta_box = wx.TextCtrl(
            self.experimental_parameters_label, -1, str(self.little_delta), size=(50, -1)
        )
        self.experimental_parameters_sizer.Add(self.little_delta_box)
        self.big_delta_label = wx.StaticText(self.experimental_parameters_label, -1, "Δ (s):")
        self.experimental_parameters_sizer.AddSpacer(5)
        self.experimental_parameters_sizer.Add(self.big_delta_label)
        self.experimental_parameters_sizer.AddSpacer(5)
        self.big_delta_box = wx.TextCtrl(self.experimental_parameters_label, -1, str(self.big_delta), size=(50, -1))
        self.experimental_parameters_sizer.Add(self.big_delta_box)
        self.find_parameters_button = wx.Button(self.experimental_parameters_label, -1, "Find Parameters")
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
                self.gradient_parameters_label, -1, "Gradient Integral Factor:"
            )
            self.gradient_parameters_sizer.Add(self.integral_factor_label)
            self.gradient_parameters_sizer.AddSpacer(5)
            self.integral_factor_box = wx.TextCtrl(
                self.gradient_parameters_label, -1, str(self.gradient_integral_factor), size=(30, -1)
            )
            self.gradient_parameters_sizer.Add(self.integral_factor_box)
            self.gradient_parameters_sizer.AddSpacer(5)

            # Have a box where the user can insert the max spectrometer gradient (default to 53G/cm for Bruker)
            self.max_gradient_label = wx.StaticText(self.gradient_parameters_label, -1, "Max Gradient (G/cm):")
            self.gradient_parameters_sizer.Add(self.max_gradient_label)
            self.gradient_parameters_sizer.AddSpacer(5)
            self.max_gradient_box = wx.TextCtrl(
                self.gradient_parameters_label, -1, str(self.max_gradient), size=(30, -1)
            )
            self.gradient_parameters_sizer.Add(self.max_gradient_box)
            self.gradient_parameters_sizer.AddSpacer(5)

            self.find_gradient_percentages_button = wx.Button(
                self.gradient_parameters_label, -1, "Find Gradient Percentages"
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
                self.gradient_parameters_label, -1, "Gradient Integral Factor:"
            )
            self.gradient_parameters_sizer.Add(self.integral_factor_label)
            self.gradient_parameters_sizer.AddSpacer(5)
            self.integral_factor_box = wx.TextCtrl(
                self.gradient_parameters_label, -1, str(self.gradient_integral_factor), size=(30, -1)
            )
            self.gradient_parameters_sizer.Add(self.integral_factor_box)
            self.gradient_parameters_sizer.AddSpacer(5)

            # Have a box where the user can insert the max spectrometer gradient (default to 53G/cm for Bruker)
            self.max_gradient_label = wx.StaticText(self.gradient_parameters_label, -1, "Max Gradient (G/cm):")
            self.gradient_parameters_sizer.Add(self.max_gradient_label)
            self.gradient_parameters_sizer.AddSpacer(5)
            self.max_gradient_box = wx.TextCtrl(
                self.gradient_parameters_label, -1, str(self.max_gradient), size=(30, -1)
            )
            self.gradient_parameters_sizer.Add(self.max_gradient_box)
            self.gradient_parameters_sizer.AddSpacer(5)

            # Have a box for DAC-G/cm conversion (default = 0.002)
            self.dac_conversion_label = wx.StaticText(
                self.gradient_parameters_label, -1, "DAC to G/cm Conversion:"
            )
            self.gradient_parameters_sizer.Add(self.dac_conversion_label)
            self.gradient_parameters_sizer.AddSpacer(5)
            self.dac_conversion_box = wx.TextCtrl(
                self.gradient_parameters_label, -1, str(self.DAC_conversion), size=(50, -1)
            )
            self.gradient_parameters_sizer.Add(self.dac_conversion_box)
            self.gradient_parameters_sizer.AddSpacer(5)

            self.find_gradient_percentages_button = wx.Button(
                self.gradient_parameters_label, -1, "Find Gradients"
            )
            self.find_gradient_percentages_button.Bind(
                wx.EVT_BUTTON, self.find_gradient_percentages
            )
            self.gradient_parameters_sizer.Add(self.find_gradient_percentages_button)
            self.gradient_parameters_sizer.AddSpacer(5)

        # Create a button to open a textbox window where a user can input the gradient values manually
        self.input_gradients_text = wx.Button(self.gradient_parameters_label, -1, "Input Manually")
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
        self.nucleus_choices = ["1H", "2H", "13C", "15N","19F", "23Na","31P"]
        self.nucleus_dropdown = wx.Choice(self.nucleus_label, -1, choices=self.nucleus_choices)
        self.nucleus_dropdown.SetSelection(0)
        self.nucleus_sizer.Add(self.nucleus_dropdown)
        self.nucleus_sizer.AddSpacer(5)
        self.nucleus_dropdown.Bind(wx.EVT_CHOICE, self.OnNucleusChoice)

        # Then have button which will allow a user to drag over a section where they wish to estimate the noise level
        # This can then be plotted as a shaded region on the plot
        self.noise_label = wx.StaticBox(self, -1, "Noise Region")
        self.noise_sizer = wx.StaticBoxSizer(self.noise_label, wx.HORIZONTAL)
        self.select_noise_button = wx.Button(self.noise_label, -1, "Select Noise Region")
        self.select_noise_button.Bind(wx.EVT_BUTTON, self.OnSelectNoise)
        self.noise_sizer.AddSpacer(5)
        self.noise_sizer.Add(self.select_noise_button)
        self.noise_sizer.AddSpacer(5)

        # Then have a TextCtrl for the minimum SNR for the diffusion coefficient to be estimated (default to 10)
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

        # Need to have a fitting sizer which will contain all the fitting buttons
        self.fitting_label = wx.StaticBox(self, -1, "Fitting")
        self.fitting_sizer = wx.StaticBoxSizer(self.fitting_label, wx.HORIZONTAL)
        self.fitting_sizer.AddSpacer(5)

        # Can then have a button which will fit the Stejskal Tanner equation at all ppms across the whole spectrum that are higher than the noise level
        self.whole_spectrum_fitting_button = wx.Button(self.fitting_label, -1, "Fit Whole Spectrum")
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

        # Can then have a button which will fit the Stejskal Tanner equation to the mean values of the points above the noise in the region of interest
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

        # Have a box containing other functions such as a button to delete a slice from the plot and repeat the fitting
        self.other_functions_label = wx.StaticBox(self, -1, "Other Functions")
        self.other_functions_sizer = wx.StaticBoxSizer(
            self.other_functions_label, wx.HORIZONTAL
        )
        self.other_functions_sizer.AddSpacer(5)
        self.delete_slice_button = wx.Button(self.other_functions_label, -1, "Delete Slice")
        self.delete_slice_button.Bind(wx.EVT_BUTTON, self.OnDeleteSlice)
        self.other_functions_sizer.Add(self.delete_slice_button)
        self.other_functions_sizer.AddSpacer(5)

        # Creating a sizer for changing the y axis limits in the spectrum
        self.intensity_label = wx.StaticBox(self, -1, "Y Axis Zoom (%):")
        self.intensity_sizer = wx.StaticBoxSizer(self.intensity_label, wx.VERTICAL)
        width = 100
        self.intensity_slider = FloatSlider(
            self.intensity_label, id=-1, value=0, minval=-1, maxval=10, res=0.01, size=(width, height)
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
        legend = self.ax_diffusion.legend(title="Slice Number", ncol=math.ceil(len(self.slice_plots)/8))
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

    def find_parameters(self, event, initial_find=False):
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
                    self.little_delta = durations_total[30] * 2
                else:
                    self.small_delta = durations_total[30]
                    self.little_delta = durations_total[30]
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
        if(initial_find==False):
            self.little_delta_box.SetValue(str(self.small_delta))
            self.big_delta_box.SetValue(str(self.big_delta))


    def error_message(self, parameter, value):
        msg = wx.MessageDialog(
            self,
            "The following parameter ({}) in the main diffusion window with value {} cannot be converted to a float. Please input a float point number for this value and try again.".format(parameter, value),
            "Error",
            wx.OK | wx.ICON_ERROR,
        )
        msg.ShowModal()
        msg.Destroy()
        return
    

    def find_max_gradient_estimate(self):
        """
        Reading difflist to get gradient values and difframp to get gradient percentages
        and working backwards to get the max spectrometer gradient strength (using the
        gradient integral factor)
        """

        try:
            with open("./lists/gp/Difframp", "r") as file:
                gradients_percent = []
                skip_line = True
                for line in file:
                    if skip_line == True:
                        if "##XYDATA= (X++(Y..Y))" not in line:
                            pass
                        else:
                            skip_line = False
                    else:
                        if "##END=" not in line:
                            gradients_percent.append(
                                float(line.split()[0]) * 100
                            )
        except:
            try:
                with open("./Difframp", "r") as file:
                    gradients_percent = []
                    skip_line = True
                    for line in file:
                        if skip_line == True:
                            if "##XYDATA= (X++(Y..Y))" not in line:
                                pass
                            else:
                                skip_line = False
                        else:
                            if "##END=" not in line:
                                gradients_percent.append(
                                    float(line.split()[0]) * 100
                                )
            except:
                return False
            
        
        try:
            with open("./lists/gp/difflist", "r") as file:
                gradients = []
                for line in file:
                    line=line.split('\n')[0]
                    gradients.append(float(line))
        except:
            try:
                with open("./difflist", "r") as file:
                    gradients = []
                    for line in file:
                        line=line.split('\n')[0]
                        gradients.append(float(line))
            except:
                return False
            

        self.max_gradient = ((gradients[-1]/gradients_percent[-1])/self.gradient_integral_factor)*100



    def find_gradient_percentages(self, event):

        if self.spectrometer == "Bruker":
            # Check the max gradient and integral factors are valid numbers
            try:
                self.max_gradient = float(self.max_gradient_box.GetValue())
            except:
                self.error_message('max gradient', self.max_gradient_box.GetValue())
                return
            try:
                self.gradient_integral_factor = float(self.integral_factor_box.GetValue())
            except:
                self.error_message('integral factor', self.integral_factor_box.GetValue())
                return
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
                ) * self.max_gradient * self.gradient_integral_factor

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
                    ) * self.max_gradient * self.gradient_integral_factor

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
            # Check the max gradient and integral factors are valid numbers
            try:
                self.DAC_conversion = float(self.dac_conversion_box.GetValue())
            except:
                self.error_message('DAC conversion', self.dac_conversion_box.GetValue())
                return
            try:
                self.gradient_integral_factor = float(self.integral_factor_box.GetValue())
            except:
                self.error_message('integral factor', self.integral_factor_box.GetValue())
                return
            try:
                self.max_gradient = float(self.max_gradient_box.GetValue())
            except:
                self.error_message('max gradient', self.max_gradient_box.GetValue())
                return
            
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

                # Multiplying gradients by the gradient integral factor
                self.gradients = self.gradients * self.gradient_integral_factor

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

        self.press = False
        self.move = False

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

        # Check to see that there are points within the noise region
        self.check_noise = self.check_noise_points()


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

    def check_noise_points(self):
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
        legend = self.ax_diffusion.legend(title="Slice Number", ncol=math.ceil(len(self.slice_plots)/8))
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
            self.UpdateDiffusionFrame()
            return True
        else:
            # Deleting this selected ROI because it contained less than 2 points
            self.DeleteSmallROI()
            self.UpdateDiffusionFrame()
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

            # Check that the number of points in the selected region of interest is greater than 2
            check_ROI = self.check_ROI(xmin, xmax)

            if(check_ROI == True):
                # Add the min and max ppm values to the array of selected regions of interest
                self.selected_regions_of_interest.append([xmin, xmax])
            else:
                # Deleting this selected ROI because it contained less than 2 points
                self.DeleteSmallROI()


            self.UpdateDiffusionFrame()

            self.press = False
            self.move = False
            self.canvas_diffusion.mpl_disconnect(self.select_ROI_press)
            self.canvas_diffusion.mpl_disconnect(self.select_ROI_move)
            self.canvas_diffusion.mpl_disconnect(self.select_ROI_release)



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
        # Generate 2 extra plots for the region fitting (I/I0 vs gradient^2 with fitted curve, histogram of diffusion coefficients within ROI)
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
        legend = self.ax_diffusion.legend(title="Slice Number", fontsize=8, ncol=math.ceil(len(self.slice_plots)/8))
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
                    color=self.ROI_color[i],
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
                (np.array(self.gradients_percent) / 100 * self.gradient_integral_factor) ** 2,
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
        


    def OnSaveFitting(self, event):
        """
        Save each of the current fits along with the data to a csv file
        (include ppm ranges used and the data title and path etc)

        Metadata JSON 
        - Data title
        - Gmax (G/cm)
        - Nucleus
        - Gyromagnetic ratio
        - Big Delta (s)
        - Little Delta (ms)
        - Global fit noise region (ppm)
        - Global fit noise value
        - Global fit minimum S/N
        
        Region of interest csv files 
        - (ppm range), G/Gmax^2, Intensity, Intensity error, diffusion coefficient, diffusion coefficient error

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
        dlg = wx.DirDialog (self, "Choose a location to save the output files to:",
            defaultPath=os.getcwd()
        )
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            # Choose a name for the output directory
            dlg2 = wx.TextEntryDialog(
                self,
                "Enter an output directory name:",
                "Output directory",
                ""
            )
            if dlg2.ShowModal() == wx.ID_OK:
                name = dlg2.GetValue()
                final_path = os.path.join(path, name)
                os.makedirs(final_path, exist_ok=True)
                self.save_metadata(final_path)
                self.save_data(final_path)
            dlg.Destroy()
            dlg2.Destroy()
        else:
            dlg.Destroy()
            return

    def save_metadata(self, path):
        """
        Saving the fit metadata to the path\metadata.json
        """

        noise_region_text = str(self.noise_x_initial) + '-' + str(self.noise_x_final)

        metadata = {"Data title": self.title, "G_max (G/cm)": self.max_gradient, "Gradient integral factor":self.gradient_integral_factor,'Nucleus':self.nucleus_choices[self.nucleus_dropdown.GetSelection()], 'Gyromagnetic ratio (rad s^-1 G^-1)':self.gamma, 'Big Delta Δ (s)': self.big_delta, 'Little delta δ (ms)':self.little_delta, 'Global fitting noise region (ppm)': noise_region_text, 'Global fit minimum S/N': self.noise_factor}
        
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
            self.mean_fitted_D_ROI = self.mean_fitted_D_ROI_total[i]
            self.mean_fitted_I0_ROI = self.mean_fitted_I0_ROI_total[i]


            if (self.fitted_gaussian_parameters[i][2]) > 1.25 * np.std(
                self.fitted_D_ROI_total[i]
            ):
                d_error = np.abs(np.std(self.fitted_D_ROI_total[i]))

            else:
                d_error = np.abs(self.fitted_gaussian_parameters[i][2])

            I0_error = np.abs(np.std(self.fitted_I0_ROI_total[i]))

            self.average_y_data_in_ROI_above_noise = self.average_y_data_in_ROI_above_noise_total[i]
            
            self.error_y_data_in_ROI_above_noise = self.error_y_data_in_ROI_above_noise_total[i]

            self.error_I_I0_in_ROI = self.error_I_I0_in_ROI_total[i]
            self.I0_average_in_ROI = self.I0_average_in_ROI_total[i]

            grad_percent = np.array(self.gradients_percent) / 100 * self.gradient_integral_factor
            G_Gmax2 = (np.array(self.gradients_percent) / 100 * self.gradient_integral_factor)**2
            grad_vals = grad_percent*self.max_gradient


            # The data will be saved in the following manner

            # Fit results
            # region_min =
            # region_max = 
            # diffusion_coefficient (cm^2/s) =
            # I0 =  

            # G, (G/G_max)^2, I, I error, I/I0, I/I0 error


            with open(os.path.join(path, 'Fit_ROI_'+str(i+1)+'.csv'), 'w') as file:
                file.write('Region min (ppm),'+ str(self.region_min)+'\n')
                file.write('Region max (ppm),'+ str(self.region_max)+'\n')
                file.write('Diffusion coefficient (cm^2/s),{:.3e}\n'.format(self.mean_fitted_D_ROI))
                file.write('Diffusion coefficient error (cm^2/s),{:.3e}\n'.format(d_error))
                file.write('I0,{:.3e}\n'.format(self.mean_fitted_I0_ROI))
                file.write('I0 error,{:.3e}\n'.format(I0_error))

                file.write('\n\n')


                file.write('G (G/cm), (G/G_max)^2, I, I error, I/I0, I/I0 error\n')
                for j in range(len(grad_vals)):
                    file.write('{:.3f},{:.3f},{:.3e},{:.3e},{:.3f},{:.3f}\n'.format(grad_vals[j], G_Gmax2[j], self.average_y_data_in_ROI_above_noise[j], self.error_y_data_in_ROI_above_noise[j], self.average_y_data_in_ROI_above_noise[j] / self.I0_average_in_ROI[j], self.error_I_I0_in_ROI[j]))




            



        


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
        self.input_gradients_text_dialog = DiffusionGradientManualInput("Input Gradients", self, self.spectrometer)