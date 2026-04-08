import wx # type: ignore
import numpy as np 
import nmrglue as ng # type: ignore
import sys
import os
import matplotlib
matplotlib.use("wxAgg")
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import (
    NavigationToolbar2WxAgg as NavigationToolbar,
)
import matplotlib.patches as patches
from scipy.interpolate import make_interp_spline # type: ignore

from SpinExplorer.SpinView.UI_objects.UI_tools import FloatSlider, PhasingSliderRange
from SpinExplorer.SpinView.Viewers.oned_view import OneDViewer

from SpinExplorer.SpinView.Modules.usta import uSTA_Dialog
from SpinExplorer.SpinView.Modules.diffusion import DiffusionFit
from SpinExplorer.SpinView.Modules.cest import CESTFrame, CESTOrder_Dialog
from SpinExplorer.SpinView.Modules.relax import RelaxFit

from SpinExplorer.SpinView.Peaks.peaks import PeakListWindow2D
from SpinExplorer.SpinView.Viewers.overlays import FileDrop, ReadProjection
from SpinExplorer.SpinView.IO import GetData
from SpinExplorer.SpinView.config import height, platform, colours, twoD_colours
from SpinExplorer.SpinView.config import reference_range_values, multiply_range_values


# A class to create a panel for viewing 2D NMR spectra
class TwoDViewer(wx.Panel):
    def __init__(self, parent, nmrdata, threeDprojection=False, fid_viewer=False):
        # Get the monitor size and set the window size to 85% of the monitor size
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        self.display_index = wx.Display.GetFromWindow(parent)
        self.width = int(1.0 * sizes[self.display_index][0])
        self.height = int(0.875 * sizes[self.display_index][1])
        self.parent = parent
        self.threeDprojection = threeDprojection
        self.fid_viewer=fid_viewer
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
        self.colours = colours
        self.colour_value = self.colours[0]

        # Initial 1D slice colour for 2D/3D spectra is set to navy
        self.colour_slice = "navy"

        # List of cmap colours for when overlaying multiple spectra
        self.cmap = "#e41a1c"
        self.cmap_neg = "#377eb8"
        self.twoD_colours = twoD_colours
        self.twoD_label_colours = twoD_colours

        self.twoD_slices_horizontal = []
        self.twoD_slices_vertical = []

        # Range of the sliders to for moving spectra left/right/up/down
        self.reference_range_values = reference_range_values
        self.reference_range = float(self.reference_range_values[0])
        self.reference_rangeX = float(self.reference_range_values[0])
        self.reference_rangeY = float(self.reference_range_values[0])

        # Range of the sliders to for moving spectra up/down in 1D spectra
        self.vertical_range_values = reference_range_values

        # Range of the sliders to for multiplying 1D spectra
        self.multiply_range_values = multiply_range_values

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

        self.start_point = None
        self.rect = None

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
        height1 = 25
        if self.threeDprojection == False:
            self.label_button = wx.Button(
                self, label="Change Labels", size=(width, height1)
            )
            self.label_button.Bind(wx.EVT_BUTTON, self.OnLabelButton)

            self.save_session_button = wx.Button(
                self, label="Save Session", size=(width, height1)
            )
            self.save_session_button.Bind(wx.EVT_BUTTON, self.OnSaveSessionButton2D)

        self.reset_button = wx.Button(self, label="Reset", size=(width, height1))
        self.reset_button.Bind(wx.EVT_BUTTON, self.OnResetButton2D)

        # Create a button to transpose the given NMR spectrum
        self.transpose_button = wx.Button(self, label="Transpose", size=(width, height1))
        self.transpose_button.Bind(wx.EVT_BUTTON, self.OnTransposeButton)

        # Create a button to stack the slices of the given NMR spectrum
        self.stack_button = wx.Button(self, label="Stack Slices", size=(width, height1))
        self.stack_button.Bind(wx.EVT_BUTTON, self.OnStackButton)

        # Create a button for Re-Processing
        self.reprocess_button = wx.Button(
            self, label="Re-Process", size=(width, height1)
        )
        self.reprocess_button.Bind(wx.EVT_BUTTON, self.OnReprocessButton)

        # Create a button to fit the diffusion data of the given NMR spectrum
        self.fit_diffusion_button = wx.Button(
            self, label="Fit Diffusion", size=(width, height1)
        )
        self.fit_diffusion_button.Bind(wx.EVT_BUTTON, self.OnFitDiffusionButton)

        # Create a button to fit the relaxation data of the given NMR spectrum
        self.fit_relax_button = wx.Button(
            self, label="Fit Relaxation", size=(width, height1)
        )
        self.fit_relax_button.Bind(wx.EVT_BUTTON, self.OnFitRelaxButton)

        # Create a button which will open a CESTView panel to analyse pseudo2D CEST data
        self.CEST_button = wx.Button(self, label="CEST Analysis", size=(width, height1))
        self.CEST_button.Bind(wx.EVT_BUTTON, self.OnCESTButton)

        # Create a button which will make the correct files in order to perform uSTA analysis
        self.uSTA_button = wx.Button(self, label="uSTA", size=(width, height1))
        self.uSTA_button.Bind(wx.EVT_BUTTON, self.OnuSTAButton)

        # Create a button to toggle the main sizer between shown and hidden
        self.toggle_button = wx.Button(self, label="Hide Options", size=(width, height1))
        self.toggle_button.Bind(wx.EVT_BUTTON, self.OnHideButton)

        self.peaklist_button = wx.Button(self, label="Read Peaks", size=(width, height1))
        self.peaklist_button.Bind(wx.EVT_BUTTON, self.OnReadPeaks)

        self.calc_intensity_button = wx.Button(
            self, label="Find Intensity", size=(width, height1)
        )
        self.calc_intensity_button.Bind(wx.EVT_BUTTON, self.OnCalculateIntensity2D)

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
            self.hide_sizer.AddSpacer(5)
            self.hide_sizer.Add(self.peaklist_button)
            self.hide_sizer.AddSpacer(5)
            self.hide_sizer.Add(self.calc_intensity_button)

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
            self.hide_sizer.AddSpacer(5)
            self.hide_sizer.Add(self.peaklist_button)
            self.hide_sizer.AddSpacer(5)
            self.hide_sizer.Add(self.calc_intensity_button)

        # Create a sizer to phase the data
        self.phasing_label = wx.StaticBox(self, -1, "Phasing:")
        self.phasing_sizer = wx.StaticBoxSizer(self.phasing_label, wx.VERTICAL)
        self.phasing_sizer1 = wx.BoxSizer(wx.HORIZONTAL)
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
        self.P0_total_value = wx.TextCtrl(self, value = "0", 
                                    size = (50,height), style = wx.TE_PROCESS_ENTER)
        self.P0_total_value.Bind(wx.EVT_TEXT_ENTER, self.P0_text_change)

        self.P1_total_value = wx.TextCtrl(self, value = "0", 
                                    size = (50,height), style = wx.TE_PROCESS_ENTER)
        self.P1_total_value.Bind(wx.EVT_TEXT_ENTER, self.P1_text_change)

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

        # Adding a button to change the range of the coarse and fine sliders (default to +/-180 and +/-10 degrees)
        self.update_phasing_range = wx.Button(self, label="Change slider range")
        self.update_phasing_range.Bind(wx.EVT_BUTTON, self.OnSliderRange2D)

        # Add a button to set the pivot point for phasing
        self.pivot_button = wx.Button(self, label="Set Pivot Point")
        self.pivot_button.Bind(wx.EVT_BUTTON, self.OnPivotButton2D)
        self.pivot_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.pivot_sizer.Add(self.update_phasing_range)
        self.pivot_sizer.AddSpacer(20)
        self.pivot_sizer.Add(self.pivot_button)

        # Add a button to remove the pivot point
        self.remove_pivot_button = wx.Button(self, label="Remove Pivot Point")
        self.remove_pivot_button.Bind(wx.EVT_BUTTON, self.OnRemovePivotButton2D)
        self.pivot_sizer.AddSpacer(20)
        self.pivot_sizer.Add(self.remove_pivot_button)


        self.phasing_sizer1.Add(self.P0_label_sizer, wx.ALIGN_TOP)
        self.phasing_sizer1.AddSpacer(10)
        self.phasing_sizer1.Add(self.P0_slider_sizer, wx.ALIGN_TOP)
        self.phasing_sizer1.AddSpacer(50)
        self.phasing_sizer1.Add(self.P1_label_sizer, wx.ALIGN_TOP)
        self.phasing_sizer1.AddSpacer(10)
        self.phasing_sizer1.Add(self.P1_slider_sizer, wx.ALIGN_TOP)

        self.phasing_sizer.Add(self.phasing_sizer1)
        self.phasing_sizer.AddSpacer(5)
        self.phasing_sizer.Add(self.pivot_sizer, wx.ALIGN_CENTER_HORIZONTAL, 1)


        # Create a sizer for changing the contour levels of the spectrum
        self.contour_label = wx.StaticBox(self, -1, "Contour Start = max(data)/x")
        self.contour_sizer = wx.StaticBoxSizer(self.contour_label, wx.VERTICAL)
        self.csizer = wx.BoxSizer(wx.HORIZONTAL)
        self.x_val = 10.00
        self.contour2_label = wx.StaticText(self, label="x:")
        self.contour_slider = FloatSlider(
            self, id=-1, value=1, minval=0, maxval=3, res=0.001, size=(200, height)
        )
        self.contour_slider.Bind(wx.EVT_SLIDER, self.OnMinContour2D)
        self.csizer.Add(self.contour2_label)
        self.csizer.AddSpacer(5)
        self.csizer.Add(self.contour_slider)
        self.contour_sizer.AddSpacer(5)
        self.contour_sizer.Add(self.csizer)
        self.contour_sizer.AddSpacer(5)
        self.contour_value_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.contour_value_sizer.AddSpacer(75)
        self.contour_value_label = wx.TextCtrl(
            self, value="10", size=(50, 20), style=wx.TE_PROCESS_ENTER
        )
        self.contour_value_label.Bind(wx.EVT_TEXT_ENTER, self.OnTextContour2D)
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
        self.multiply_label = wx.StaticBox(self, -1, "Multiply 2D Data by " + "n" + ":")
        self.multiply_sizer = wx.StaticBoxSizer(self.multiply_label, wx.VERTICAL)
        self.multiply_inner_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.multiply_ranges = multiply_range_values
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
        self.multiply_value_label = wx.TextCtrl(self, value = "1.0", 
                                    size = (50,height), style = wx.TE_PROCESS_ENTER)
        self.multiply_value_label.Bind(wx.EVT_TEXT_ENTER, self.multiply_text_change)
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

        self.reference_range_chooserX.SetSelection(0)
        self.reference_range_chooserY.SetSelection(0)

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
        self.move_x_value_label = wx.TextCtrl(self, value = "0.0", 
                                    size = (50,height), style = wx.TE_PROCESS_ENTER)
        self.move_x_value_label.Bind(wx.EVT_TEXT_ENTER, self.move_xtext_change)
        self.move_values.Add(self.move_x_value_label)
        self.move_values.AddSpacer(50)
        self.move_values.Add(wx.StaticText(self, label="Y Movement (ppm):"))
        self.move_values.AddSpacer(10)
        self.move_y_value_label = wx.TextCtrl(self, value = "0.0", 
                                    size = (50,height), style = wx.TE_PROCESS_ENTER)
        self.move_y_value_label.Bind(wx.EVT_TEXT_ENTER, self.move_ytext_change)
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
        self.buttons_sizer.AddSpacer(5)
        if self.threeDprojection == False:
            self.buttons_sizer.Add(self.fit_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.buttons_sizer.AddSpacer(5)
        self.buttons_sizer.Add(self.hide_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.bottom_right_sizer.Add(self.buttons_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)

        self.bottom_right_sizer.AddSpacer(5)

        self.bottom_sizer.Add(self.bottom_right_sizer)

        self.slice_mode = None

    def OnSliderRange2D(self, event):
        """
        Creating a popout where a user can update the slider range
        """
        self.slider_range_window = PhasingSliderRange("Phasing slider ranges",self)


    def OnTextContour2D(self, event):
        """
        First check that the input is valid.
        Then Update the Slider value
        """
        try:
            self.x_val = float(self.contour_value_label.GetValue())
            self.OnMinContour2D(event, textcontrol=True)
        except:
            self.contour_value_label.SetValue(self.x_val)

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
        from SpinExplorer.SpinProcess.SpinProcess import SpinProcess

        reprocessing_frame = SpinProcess(self, reprocess=True)
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

    def OnCalculateIntensity2D(self, event):
        """
        If the find peaks window exits. Check that none of the toggled buttons
        are selected. If they are turn them off.

        Then ask the user to drag over a region to find the max intensity
        """

        for window in wx.GetTopLevelWindows():
            if isinstance(window, wx.Frame) and window.GetTitle() == "Peak Lists":
                self.peaklist_frame.turn_off_togglebuttons()

        dlg = wx.MessageDialog(
            None,
            "Drag over a selected region of the spectrum. The max intensity, mean intensity, integral, and standard deviation of the selected region are then outputted.",
            "Find Intensity",
            wx.OK,
        )
        dlg.ShowModal()
        dlg.Destroy()

        # If drag, finds new peaks
        self.calculate_press = self.fig.canvas.mpl_connect(
            "button_press_event", self.on_press_calculate
        )
        self.calculate_release = self.fig.canvas.mpl_connect(
            "button_release_event", self.on_release_calculate
        )
        self.calculate_motion = self.fig.canvas.mpl_connect(
            "motion_notify_event", self.on_motion_calculate
        )

    def on_press_calculate(self, event):
        """
        This is activated when the mouse is clicked when calculate intensity
        is clicked
        """
        x, y = self.ax.transData.inverted().transform((event.x, event.y))

        if x != None and y != None:
            self.start_point = (x, y)

            # Create the rectangle
            self.rect = patches.Rectangle(
                self.start_point, 0, 0, linewidth=1, edgecolor="red", facecolor="none"
            )
            self.ax.add_patch(self.rect)
            self.fig.canvas.draw()
            self.UpdateFrame()

    def on_motion_calculate(self, event):
        """
        This is activated when the mouse is moved when calculate intensity
        was clicked
        """
        if not self.start_point:
            return

        x, y = self.ax.transData.inverted().transform((event.x, event.y))

        if x != None and y != None:

            # Update rectangle size
            x0, y0 = self.start_point
            x1, y1 = x, y
            width = x1 - x0
            height = y1 - y0

            self.rect.set_width(width)
            self.rect.set_height(height)
            self.rect.set_xy((x0, y0))
            self.canvas.draw_idle()
            self.UpdateFrame()

    def on_release_calculate(self, event):
        """
        This is activated when the mouse is released when calculate
        intensity is clicked
        """
        if not self.start_point:
            return
        x, y = self.ax.transData.inverted().transform((event.x, event.y))

        if x != None and y != None:
            x0, y0 = self.start_point
            x1, y1 = x, y
            xmin, xmax = sorted([x0, x1])
            ymin, ymax = sorted([y0, y1])

        self.fig.canvas.mpl_disconnect(self.calculate_motion)
        self.fig.canvas.mpl_disconnect(self.calculate_press)
        self.fig.canvas.mpl_disconnect(self.calculate_release)

        self.intensity2Dpopout([xmin, xmax], [ymin, ymax])

        # Cleanup
        self.start_point = None
        self.rect.set_visible(False)
        self.rect = None
        self.canvas.draw()
        self.UpdateFrame()

    def intensity2Dpopout(self, xlim, ylim):
        """
        xlim and ylim are the limits of the rectangle dragged by the user.
        This function finds the desired intensity outputs to show the user.
        """
        if self.multiplot_mode == False:
            data = self.nmrdata.data * self.multiply_factor
            xppms = self.new_x_ppms
            yppms = self.new_y_ppms
        else:
            data = self.values_dictionary[self.active_plot_index]["z_data"] * self.values_dictionary[self.active_plot_index]['multiply factor']
            xppms = self.values_dictionary[self.active_plot_index]["new_x_ppms"]
            yppms = self.values_dictionary[self.active_plot_index]["new_y_ppms"]


        if (
            xlim[0] < np.min(xppms)
            or xlim[1] > np.max(xppms)
            or ylim[0] < np.min(yppms)
            or ylim[1] > np.max(yppms)
        ):
            # Output a return to say that the selected rectangle goes outside the chemical shift range of the plot
            dlg = wx.MessageDialog(
                None,
                "The selected rectangle goes outside the chemical shift range of the plot. Please try again.",
                "Find Intensity",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return

        # Find the nearest index of the limits to those in xppms and yppms

        self.x_index_initial = np.abs(xppms - xlim[0]).argmin()
        self.x_index_final = np.abs(xppms - xlim[1]).argmin()
        self.y_index_initial = np.abs(yppms - ylim[0]).argmin()
        self.y_index_final = np.abs(yppms - ylim[1]).argmin()

        data_selected = data[
            self.x_index_final : self.x_index_initial,
            self.y_index_final : self.y_index_initial,
        ]

        max_value = np.max(data_selected)
        mean_value = np.mean(data_selected)
        integral = np.sum(data_selected)
        stdev = np.std(data_selected)

        wx.MessageBox(
            "Maximum Intensity:\n{:E}\nMean Intensity:{:E}\nIntegral:\n{:E}\nStandard deviation:\n{:E}".format(
                max_value,
                mean_value,
                integral,
                stdev,
            ),
            "Find Intensity",
            wx.OK | wx.ICON_INFORMATION,
        )

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
                self.contour_value_label.SetValue("10")
                self.contour_levels_slider.SetValue(20)
                self.move_x_slider.SetValue(0)
                self.move_y_slider.SetValue(0)
                self.move_x_value_label.SetLabel("0.00")
                self.move_y_value_label.SetLabel("0.00")
                self.multiply_slider.SetValue(1.0)
                self.multiply_value_label.SetLabel("1.0")
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
                    self.contour_value_label.SetValue("10")
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
                    self.contour_value_label.SetValue("10")
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
        self.dlg = wx.Dialog(None, title="Change Labels")
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
        self.dlg.Show()

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


        self.UpdateFrame()
        self.dlg.Destroy()

    def OnuSTAButton(self, event):
        # Producing an output to the user telling them this will produce a .data file required for uSTA analysis. Full uSTA implementation is in development.
        dlg = wx.MessageDialog(
            self,
            "This will produce a .data file required for downstream saturation transfer analysis using uSTA from the Baldwin lab (https://usta.chem.ox.ac.uk). Do you want to continue?",
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

    def OnReadPeaks(self, event):
        """
        Have a popout window where a user can choose a peaklist (in Sparky/tabular format)
        This can also suggest if it is advisable to transpose the peaklist to match
        the data.

        Supported file formats are:
        (test).ft2.list files (Sparky/tabular format)
        CCPN peak table output 
        """

        for window in wx.GetTopLevelWindows():
            if isinstance(window, wx.Frame) and window.GetTitle() == "Peak Lists":
                # The window already exists (return)
                return

        self.peaklist_frame = PeakListWindow2D(title="Peak Lists", parent=self)

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
        self.mouse_wheel_connect = self.fig.canvas.Bind(wx.EVT_MOUSEWHEEL, self.on_mouse_wheel)

        contour_start = np.max(self.nmrdata.data) / 10  # contour level start value
        self.contour_num = 20  # number of contour levels
        self.contour_factor = 1.20  # scaling factor between contour levels
        # calculate contour levels
        self.cl = contour_start * self.contour_factor ** np.arange(self.contour_num)
        self.cl_neg = -contour_start * self.contour_factor ** np.flip(
            np.arange(self.contour_num)
        )

        if(self.fid_viewer==True):
            self.nmrdata.dic, self.nmrdata.data = ng.pipe_proc.tp(self.nmrdata.dic,self.nmrdata.data, auto=True)
            self.nmrdata.dic, self.nmrdata.data = ng.pipe_proc.di(self.nmrdata.dic,self.nmrdata.data)
            self.nmrdata.dic, self.nmrdata.data = ng.pipe_proc.tp(self.nmrdata.dic,self.nmrdata.data, auto=True)

        # Get ppm values for x and y axis
        if self.nmrdata.file != ".":
            self.uc0 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=0)
            self.uc1 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data, dim=1)
        else:
            udic = ng.bruker.guess_udic(self.nmrdata.dic, self.nmrdata.data)
            self.uc0 = ng.fileiobase.uc_from_udic(udic, dim=0)
            self.uc1 = ng.fileiobase.uc_from_udic(udic, dim=1)
        if(self.fid_viewer==False):
            self.ppms_0 = self.uc0.ppm_scale()
            self.ppms_1 = self.uc1.ppm_scale()
        else:
            self.ppms_0 = np.arange(0, len(self.uc0.ppm_scale()),1)
            self.ppms_1 = np.arange(0, len(self.uc1.ppm_scale()),1)
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
        if(self.fid_viewer==False):
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
        self.OnMinContour2D(wx.EVT_SCROLL, textcontrol=True)

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
            if(self.fid_viewer==False):
                self.ax.set_xlim([max(self.new_x_ppms), min(self.new_x_ppms)])
                self.ax.set_ylim([max(self.new_y_ppms), min(self.new_y_ppms)])
            else:
                self.ax.set_xlim([min(self.new_x_ppms), max(self.new_x_ppms)])
                self.ax.set_ylim([min(self.new_y_ppms), max(self.new_y_ppms)])
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

            for window in wx.GetTopLevelWindows():
                if isinstance(window, wx.Frame) and window.GetTitle() == "Peak Lists":
                    # Need to swap the order of the peaklist dimensions
                    try:
                        for (
                            peaklist_name,
                            dictionary,
                        ) in self.peaklist_frame.peak_list_dictionary.items():
                            shift1 = dictionary["shift1"]
                            shift2 = dictionary["shift2"]
                            self.peaklist_frame.peak_list_dictionary[peaklist_name][
                                "shift1"
                            ] = shift2
                            self.peaklist_frame.peak_list_dictionary[peaklist_name][
                                "shift2"
                            ] = shift1
                            axis_names_old = self.peaklist_frame.names[peaklist_name]
                            self.peaklist_frame.names[peaklist_name] = [axis_names_old[1], axis_names_old[0]]

                    except:
                        pass
                else:
                    pass

            self.OnMinContour2D(wx.EVT_SCROLL, textcontrol=True)
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

            if(self.fid_viewer==False):
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

            self.OnMinContour2D(wx.EVT_SCROLL, textcontrol=True)
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
                self.stacks = Stack2D(
                    title="Stacked Slices - " + self.parent.title, parent=self
                )
            else:
                self.stacks = Stack2D(
                    title="Stacked Slices - " + self.parent.title, parent=self
                )
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
            self.diffusion = DiffusionFit(
                title="Diffusion Fit - " + self.parent.title, parent=self
            )
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
            # See if the number of slices is even, if it is not even then the data is likely not CEST data
            if(len(self.nmrdata.data.T)%2 != 0):
                self.continue_window = wx.MessageDialog(
                        self,
                        "There is not an even number of slices in the indirect dimension. The data should be aranged as an array of on-resonance/off-resonance for each saturation frequency.",
                        "Warning",
                        wx.OK | wx.ICON_WARNING,
                    )
                self.continue_window.ShowModal()
                self.continue_window.Destroy()
                return


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
            title="CEST - " + self.parent.title,
            parent=self,
            CESTArrayOrder=self.CESTArrayOrder,
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
            self.RelaxFit = RelaxFit(
                title="Relaxation Fit - " + self.parent.title, parent=self
            )
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

    def OnMinContour2D(self, event, textcontrol=False):
        # Function to update the contour levels when the user changes the number of contour levels
        if textcontrol == False:
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
                zorder=1,
            )
            self.contour1_neg = self.ax.contour(
                self.Y,
                self.X,
                self.nmrdata.data * self.multiply_factor,
                self.cl_neg,
                colors=self.cmap_neg,
                linewidths=self.linewidth,
                zorder=1,
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
                    zorder=1,
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

        if textcontrol == False:
            self.contour_value_label.SetValue(
                "{:.2f}".format(10 ** float(self.contour_slider.GetValue()))
            )

        for window in wx.GetTopLevelWindows():
            self.peaklist_colours = [
                "black",
                "gray",
                "saddlebrown",
                "purple",
                "purple",
                "blue",
                "red",
                "orange",
            ]
            if isinstance(window, wx.Frame) and window.GetTitle() == "Peak Lists":
                # Plot Peaklists
                count = 0
                self.points = []
                self.annotations = []
                for (
                    peaklist_name,
                    dictionary,
                ) in self.peaklist_frame.peak_list_dictionary.items():

                    if (
                        self.peaklist_frame.select_peak_button.GetValue() == True
                        or self.peaklist_frame.select_peaks_button.GetValue() == True
                    ):
                        if peaklist_name == self.peaklist_frame.selected_peaklist:
                            if "N/A" in self.peaklist_frame.selected_peak_indexes:
                                cs = self.peaklist_colours[count]
                            else:
                                cs = []
                                for i, peak in enumerate(dictionary["peak_name"]):
                                    if i in self.peaklist_frame.selected_peak_indexes:
                                        cs.append("darkviolet")
                                    else:
                                        cs.append(self.peaklist_colours[count])
                        else:
                            cs = self.peaklist_colours[count]
                    else:
                        cs = self.peaklist_colours[count]

                    shift1 = dictionary["shift1"]
                    shift2 = dictionary["shift2"]
                    self.points.append(
                        self.ax.scatter(
                            shift1,
                            shift2,
                            s=5,
                            marker="o",
                            c=cs,
                            picker=5,
                            zorder=2,
                        )
                    )
                    count += 1

                    # Annotation for hover
                    self.annotations.append(
                        self.ax.annotate(
                            "",
                            xy=(0, 0),
                            xytext=(15, 15),
                            textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"),
                        )
                    )
                    self.annotations[-1].set_visible(False)
                    # adjust_text(self.annotations, ax=self.ax)

                    # Connect event
                    self.hover_connect = self.canvas.mpl_connect(
                        "motion_notify_event", self.on_hover
                    )

            else:
                pass

        self.OnSliderScroll2D(wx.EVT_SCROLL)
        self.OnIntensityScroll2D(wx.EVT_SCROLL)

        self.UpdateFrame()

    def on_hover(self, event):

        if event.inaxes != None:
            # Calculate distance from mouse to each point

            for i, points in enumerate(self.points):
                cont, ind = points.contains(event)

                if cont:
                    # Show annotation
                    index = ind["ind"][0]  # first index found
                    peaklist_name = self.peaklist_frame.peak_list_choices[i]
                    dictionary = self.peaklist_frame.peak_list_dictionary[peaklist_name]

                    peakname = (
                        dictionary["peak_name"][index] + " (" + peaklist_name + ")"
                    )
                    x = dictionary["shift1"][index]
                    y = dictionary["shift2"][index]
                    self.annotations[i].xy = (x, y)

                    text = peakname
                    self.annotations[i].set_text(text)
                    self.annotations[i].set_color(self.peaklist_colours[i])
                    self.annotations[i].set_position((36, i * 36))
                    self.annotations[i].set_visible(True)
                    self.canvas.draw_idle()
                else:
                    if self.annotations[i].get_visible():
                        self.annotations[i].set_visible(False)

            # adjust_text(self.annotations, ax=self.ax, time_lim=5)
            self.canvas.draw_idle()

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
        self.OnMinContour2D(event, textcontrol=True)

    def multiply_text_change(self, event):
        self.multiply_factor = float(self.multiply_value_label.GetValue())
        self.multiply_slider.SetValue(self.multiply_factor)
        self.ApplyMultiplication(event)

    def OnMultiplyScroll2D(self, event):
        self.multiply_factor = float(self.multiply_slider.GetValue())
        self.multiply_value_label.SetLabel(
            "{:.2f}".format(float(self.multiply_slider.GetValue()))
        )
        self.ApplyMultiplication(event)

    def ApplyMultiplication(self, event):
        if self.multiplot_mode == True:
            if self.select_all_checkbox.IsChecked() == False:
                self.values_dictionary[self.active_plot_index][
                    "multiply factor"
                ] = self.multiply_factor
            else:
                for i in range(len(self.twoD_slices_horizontal)):
                    self.values_dictionary[i]["multiply factor"] = self.multiply_factor
        self.OnMinContour2D(event, textcontrol=True)

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
        self.OnMinContour2D(event, textcontrol=True)

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

        self.OnMinContour2D(event, textcontrol=True)

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

    def move_xtext_change(self, event):
        self.x_movement = float(self.move_x_value_label.GetValue())
        self.move_x_slider.SetValue(self.x_movement)
        self.MoveX()

    def OnMoveX(self, event):
        # update x-axis
        self.x_movement = float(self.move_x_slider.GetValue())
        self.move_x_value_label.SetLabel("{:.4f}".format(self.x_movement))
        self.MoveX()
    
    def MoveX(self):
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
            self.OnMinContour2D(wx.EVT_SCROLL, textcontrol=True)
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
            self.OnMinContour2D(wx.EVT_SCROLL, textcontrol=True)
            self.UpdateFrame()


    def move_ytext_change(self, event):
        self.y_movement = float(self.move_y_value_label.GetValue())
        self.move_y_slider.SetValue(self.y_movement)
        self.MoveY()

    def OnMoveY(self, event):
        # update y-axis
        self.y_movement = float(self.move_y_slider.GetValue())
        self.move_y_value_label.SetLabel("{:.4f}".format(self.y_movement))
        self.MoveY()

    def MoveY(self):
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
            self.OnMinContour2D(wx.EVT_SCROLL, textcontrol=True)
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
            self.OnMinContour2D(wx.EVT_SCROLL, textcontrol=True)
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


    def on_key_2d(self, event):
        
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
                        if(self.fid_viewer==False):
                            data = self.nmrdata.data[
                                    :, self.uc1(str(self.new_y_ppms[1]) + "ppm")
                                ]
                        else:
                            data = self.nmrdata.data[
                                    :, int(self.new_y_ppms[1])
                                ]
                        (self.line1,) = self.axes1D.plot(
                            self.new_x_ppms,
                            data
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
                        if(self.fid_viewer==False):
                            data = self.nmrdata.data[
                                self.uc0(str(self.new_x_ppms[1]) + "ppm"), :
                            ]
                        else:
                            data = self.nmrdata.data[
                                int(self.new_x_ppms[1]), :
                            ]
                        (self.line3,) = self.axes1D_2.plot(
                            data
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
                    if(self.fid_viewer==False):
                        data = self.nmrdata.data[
                            :, self.uc1(str(self.y1 - self.y_movement) + "ppm")
                        ]
                    else:
                        data = self.nmrdata.data[
                            :, int(self.y1 - self.y_movement)
                        ]
                    self.line1.set_ydata(data*self.multiply_factor)
                    self.line2.set_ydata([self.y1])
                    self.line1.set_xdata(self.ppms_0 + self.x_movement)
                    self.OnSliderScroll2D(wx.EVT_SCROLL)
                    self.UpdateFrame()
                if self.line3.get_visible() == True:
                    if(self.fid_viewer == False):
                        data = self.nmrdata.data[
                            self.uc0(str(self.x1 - self.x_movement) + "ppm"), :
                        ]
                    else:
                        data = self.nmrdata.data[
                            int(self.x1 - self.x_movement), :
                        ]
                    self.line3.set_xdata(
                        data
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
    
    def P0_text_change(self, event):
        self.total_P0 = float(self.P0_total_value.GetValue())
        self.P0_slider.SetValue(self.total_P0)
        self.P0_slider_fine.SetValue(0.0)
        self.total_P1 = self.P1_slider.GetValue() + self.P1_slider_fine.GetValue()
        self.phase2D()

    def P1_text_change(self, event):
        self.total_P1 = float(self.P1_total_value.GetValue())
        self.P1_slider.SetValue(self.total_P1)
        self.P1_slider_fine.SetValue(0.0)
        self.total_P0 = self.P0_slider.GetValue() + self.P0_slider_fine.GetValue()
        self.phase2D()

    def phase2D(self):
        # Phase the 2D data with the combined coarse/fine phasing values and plot the result
        if self.multiplot_mode == False:
            try:
                if self.line1.get_visible() == True:
                    if(self.fid_viewer==False):
                        data = (
                            self.nmrdata.data[
                                :, self.uc1(str(self.y1 - self.y_movement) + "ppm")
                            ]
                            * self.multiply_factor
                        )
                    else:
                        data = (
                            self.nmrdata.data[
                                :, int(self.y1 - self.y_movement)
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
                    if(self.fid_viewer==False):
                        data = (
                            self.nmrdata.data[
                                self.uc0(str(self.x1 - self.x_movement) + "ppm"), :
                            ]
                            * self.multiply_factor
                        )
                    else:
                        data = (
                            self.nmrdata.data[
                                int(self.x1 - self.x_movement), :
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
            nmr_data_0 = GetData(self, file=self.main_frame.parent.nmrdata.file)
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
