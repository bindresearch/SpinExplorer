import wx # type: ignore
import numpy as np
import nmrglue as ng  # type: ignore
import matplotlib
import sys
import os
matplotlib.use("wxAgg")
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import (
    NavigationToolbar2WxAgg as NavigationToolbar,
)
from scipy.interpolate import make_interp_spline # type: ignore

from SpinExplorer.SpinView.UI_objects.UI_tools import FloatSlider, PhasingSliderRange
from SpinExplorer.SpinView.Viewers.overlays import FileDrop
from SpinExplorer.SpinView.config import *

from abc import ABC, abstractmethod

class spectrum_window(wx.Frame, ABC):
    
    def __init__(self, parent, nmrdata, uc0 = None)->None:
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        self.display_index = wx.Display.GetFromWindow(parent)
        self.width = int(1.0*sizes[self.display_index][0])
        self.height = int(0.875*sizes[self.display_index][1])
        self.parent = parent
        self.uc0_initial = uc0
        self.stack = False
        self.uc0 = uc0
        super().__init__(parent, id=wx.ID_ANY, size = (self.width,self.height))
        self.nmrdata = nmrdata
    
    def make_static_vbox(self, label, *items, spacing=15):
        """
        Creates a StaticBoxSizer with optional leading/trailing spacer.
        Items can be widgets or (widget, flag, proportion) tuples.
        """
        static_box = wx.StaticBox(self, -1, label)
        sizer = wx.StaticBoxSizer(static_box, wx.VERTICAL)
        sizer.AddSpacer(spacing)
        for i, item in enumerate(items):
            if isinstance(item, tuple):
                widget, flag, proportion = item
                sizer.Add(widget, proportion, flag)
            else:
                sizer.Add(item)
            if i < len(items) - 1:
                sizer.AddSpacer(spacing)
        sizer.AddSpacer(spacing)
        return sizer
    
    def make_static_hbox(self, label, *items, spacing=5, outer=0):
        static_box = wx.StaticBox(self, -1, label)
        sizer = wx.StaticBoxSizer(static_box, wx.HORIZONTAL)
        if outer:
            sizer.AddSpacer(outer)
        for i, item in enumerate(items):
            if isinstance(item, tuple):
                widget, flag, proportion = item
                sizer.Add(widget, proportion, flag)
            else:
                sizer.Add(item)
            if i < len(items) - 1:
                sizer.AddSpacer(spacing)
        if outer:
            sizer.AddSpacer(outer)
        return sizer

    def make_vbox(self, *items, spacing=10):
        """
        Items can be widgets or (widget, flag, proportion) tuples.
        e.g. make_vbox(label, (slider, wx.ALIGN_CENTER_HORIZONTAL, 0))
        """
        sizer = wx.BoxSizer(wx.VERTICAL)
        for i, item in enumerate(items):
            if isinstance(item, tuple):
                widget, flag, proportion = item
                sizer.Add(widget, proportion, flag)
            else:
                sizer.Add(item)
            if i < len(items) - 1:
                sizer.AddSpacer(spacing)
        return sizer

    def make_hbox(self, *items, spacing=10, outer=5):
        """
        Items can be widgets or (widget, flag, proportion) tuples.
        e.g. make_hbox((sizer, wx.ALIGN_TOP, 0), ...)
        """
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.AddSpacer(outer)
        for i, item in enumerate(items):
            if isinstance(item, tuple):
                widget, flag, proportion = item
                sizer.Add(widget, proportion, flag)
            else:
                sizer.Add(item, 0, wx.ALIGN_TOP)
            if i < len(items) - 1:
                sizer.AddSpacer(spacing)
        sizer.AddSpacer(outer)
        return sizer
    
    def create_canvas(self)->None:
        self.panel = wx.Panel(self)
        self.fig = Figure()
        self.canvas = FigCanvas(self,-1,self.fig)
        self.toolbar = NavigationToolbar(self.canvas)
    
    
    def set_initial_variables(self):
        # Colours for 1D lines
        self.colours = colours
        self.colour_value = self.colours[0]

        # Initial 1D slice colour for 2D/3D spectra is set to navy
        self.colour_slice = "navy"


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


        self.linewidth = 1.0
        self.linewidth1D = 1.5

        self.x_difference = 0
        self.y_difference = 0


        # Default options for pivot point for P1 phasing
        self.pivot_x_default = 0
        self.pivot_x = self.pivot_x_default

        self.pivot_y_default = 0
        self.pivot_y = self.pivot_y_default


        # Suppress complex warning from numpy
        import warnings

        # warnings.simplefilter("ignore", np.ComplexWarning)  # For old numpy versions
        warnings.simplefilter(
            "ignore", np.exceptions.ComplexWarning
        )  # For new numpy versions
    
    
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


        self.create_phasing_panel()
        self.create_yaxis_limits_panel()
        self.move_spec_horizontal_panel()
        self.move_spec_vertical_panel()
        self.choose_colour_panel()
        self.change_linewidth_panel()
        self.multiply_spec_panel()

        
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

        # self.load_session_button =  wx.Button(self, label="Load Session", size=(130, 30))
        # self.load_session_button.Bind(wx.EVT_BUTTON, self.OnLoadSession)

        # Button stack
        self.button_sizers = self.make_vbox(
            self.max_button,
            self.baseline,
            self.reset_button,
            self.subtract_button,
            self.reprocess_button,
            # self.load_session_button,
            self.save_button,
            self.save_session_button,
            self.hide_button,
            spacing=5,
        )

        # Intensity/reference/vertical stack
        platform_spacer = 15 if platform == "linux" else 10
        self.intensity_reference_sizer = self.make_vbox(
            self.contour_sizer,
            self.reference_total,
            self.vertical_sizer,
            spacing=platform_spacer,
        )

        # Top-left: plot selection and display options
        spacer1 = int(self.parent.width / 100)
        self.top_left_sizer = self.make_hbox(
            self.select_plot_sizer,
            self.colour_sizer,
            self.linewidth_sizer,
            self.multiply_total,
            spacing=spacer1,
            outer=0,
        )

        # Left side: top options + phasing
        self.left_sizer = self.make_vbox(
            self.top_left_sizer,
            self.phasing_sizer,
            spacing=5,
        )

        # Bottom-right: intensity/reference/vertical
        self.bottom_right_sizer = self.make_hbox(
            self.intensity_reference_sizer,
            spacing=5,
            outer=0,
        )

        # Full bottom: left + right + buttons
        self.bottom_sizer = self.make_hbox(
            self.left_sizer,
            self.bottom_right_sizer,
            self.button_sizers,
            spacing=10,
            outer=5,
        )


    def create_phasing_panel(self):
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
        self.P0_slider.Bind(wx.EVT_SLIDER, self.OnSliderScroll)
        self.P1_slider.Bind(wx.EVT_SLIDER, self.OnSliderScroll)
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
        self.P0_slider_fine.Bind(wx.EVT_SLIDER, self.OnSliderScroll)
        self.P1_slider_fine.Bind(wx.EVT_SLIDER, self.OnSliderScroll)
        self.P0_total = wx.StaticText(self, label="P0 (Total):", size=(70, height))
        self.P1_total = wx.StaticText(self, label="P1 (Total):", size=(70, height))
        self.P0_total_value = wx.TextCtrl(self, value = "0", 
                                    size = (70, height), style = wx.TE_PROCESS_ENTER)
        self.P0_total_value.Bind(wx.EVT_TEXT_ENTER, self.P0_text_change)

        self.P1_total_value = wx.TextCtrl(self, value = "0", 
                                    size = (70,height), style = wx.TE_PROCESS_ENTER)
        self.P1_total_value.Bind(wx.EVT_TEXT_ENTER, self.P1_text_change)

        # Adding a button to change the range of the coarse and fine sliders (default to +/-180 and +/-10 degrees)
        self.update_phasing_range = wx.Button(self, label="Change slider range")
        self.update_phasing_range.Bind(wx.EVT_BUTTON, self.OnSliderRange)

        # Adding a button to set the pivot point for phasing
        self.pivot_button = wx.Button(self, label="Set Pivot Point")
        self.pivot_button.Bind(wx.EVT_BUTTON, self.OnPivotButton)
        self.pivot_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.pivot_sizer.Add(self.update_phasing_range)
        self.pivot_sizer.AddSpacer(10)
        self.pivot_sizer.Add(self.pivot_button)

        # Adding a button to remove the pivot point
        self.remove_pivot_button = wx.Button(self, label="Remove Pivot Point")
        self.remove_pivot_button.Bind(wx.EVT_BUTTON, self.OnRemovePivotButton)
        self.pivot_sizer.AddSpacer(10)
        self.pivot_sizer.Add(self.remove_pivot_button)

        # P0 and P1 label/slider sizers
        for prefix, idx in (('P0','0'), ('P1','1')):
            labels_sizer = self.make_vbox(
                getattr(self, f'{prefix}_label'),
                getattr(self, f'{prefix}_label_fine'),
                getattr(self, f'{prefix}_total'),
            )
            sliders_sizer = self.make_vbox(
                (getattr(self, f'{prefix}_slider'),       wx.ALIGN_CENTER_HORIZONTAL, 0),
                (getattr(self, f'{prefix}_slider_fine'),  wx.ALIGN_CENTER_HORIZONTAL, 0),
                (getattr(self, f'{prefix}_total_value'),  wx.ALIGN_CENTER_HORIZONTAL, 5),
            )
            setattr(self, f'p{idx}_sizer_labels',  labels_sizer)
            setattr(self, f'p{idx}_sizer_sliders', sliders_sizer)

        # Combine into horizontal phasing sizer
        self.phasing_sizer1 = self.make_hbox(
            self.p0_sizer_labels,
            self.p0_sizer_sliders,
            self.p1_sizer_labels,
            self.p1_sizer_sliders,
        )


        self.phasing_sizer.Add(self.phasing_sizer1)
        self.phasing_sizer.AddSpacer(10)
        self.phasing_sizer.Add(self.pivot_sizer, wx.ALIGN_CENTER, 1)


    def create_yaxis_limits_panel(self):
        width = int(self.parent.width / 4.5)

        self.intensity_slider = FloatSlider(
            self, id=-1, value=0, minval=-1, maxval=10, res=0.01, size=(width, height)
        )
        self.intensity_slider.Bind(wx.EVT_SLIDER, self.OnIntensityScroll1D)

        # Sizer assembly
        self.contour_sizer = self.make_static_vbox(
            "Y Axis Zoom (%):",
            self.intensity_slider,
            spacing=5,)

    def move_spec_horizontal_panel(self):
        # Layout dimensions
        total_zoom_width = self.contour_sizer.GetMinSize()[0] - 15
        if total_zoom_width < 150:
            width = int(total_zoom_width * 0.4)
            slider_width = int(total_zoom_width * 0.6)
        else:
            width = 55
            slider_width = total_zoom_width - 70  # 15 padding + 55 width

        # Widgets
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

        self.reference_value_label = wx.TextCtrl(
            self, value="0.0", size=(70, height), style=wx.TE_PROCESS_ENTER
        )
        self.reference_value_label.Bind(wx.EVT_TEXT_ENTER, self.OnReferenceText)

        self.reference_range_chooser = wx.ComboBox(
            self,
            value=self.reference_range_values[0],
            choices=self.reference_range_values,
            size=(width, height),
        )
        self.reference_range_chooser.Bind(wx.EVT_COMBOBOX, self.OnReferenceCombo)
        self.reference_range_chooser.SetSelection(0)

        self.reference_range_text = wx.StaticText(self, label="Range")

        # Sizer assembly
        self.reference_sizer = self.make_vbox(
            (self.reference_slider,      wx.ALIGN_CENTER_HORIZONTAL, 5),
            (self.reference_value_label, wx.ALIGN_CENTER_HORIZONTAL, 5),
            spacing=5,
        )

        self.reference_sizer2 = self.make_vbox(
            (self.reference_range_chooser, wx.ALIGN_CENTER_HORIZONTAL, 5),
            (self.reference_range_text,    wx.ALIGN_CENTER_HORIZONTAL, 5),
            spacing=5,
        )

        self.reference_sizer_full = self.make_hbox(
            self.reference_sizer,
            self.reference_sizer2,
            spacing=5,
            outer=0,
        )

        self.reference_total = self.make_static_vbox(
            "Move \u2190/\u2192 (ppm):",
            self.reference_sizer_full,
            spacing=0,
        )

        # Pad to match contour sizer width if needed
        padding = self.contour_sizer.GetMinSize()[0] - self.reference_total.GetMinSize()[0]
        if padding > 0:
            self.reference_sizer_full.AddSpacer(padding)


    def move_spec_vertical_panel(self):
        # Data range
        self.vertical_range = int(max(self.nmrdata.data))
        vertical_scale = self.vertical_range * float(self.vertical_range_values[0]) / 100

        # Widgets
        self.vertical_slider = FloatSlider(
            self,
            id=-1,
            value=0,
            minval=-vertical_scale,
            maxval=vertical_scale,
            res=vertical_scale / 100,
            size=(slider_width, height),
        )
        self.vertical_slider.Bind(wx.EVT_SLIDER, self.OnVerticalScroll1D)

        self.vertical_range_chooser = wx.ComboBox(
            self,
            value=self.vertical_range_values[0],
            choices=self.vertical_range_values,
            size=(width, height),
        )
        self.vertical_range_chooser.Bind(wx.EVT_COMBOBOX, self.OnVerticalCombo)
        self.vertical_range_chooser.SetSelection(0)

        self.vertical_range_label = wx.StaticText(self, label="Range")

        # Sizer assembly
        self.vertical_sizer2 = self.make_vbox(
            (self.vertical_range_chooser, wx.ALIGN_CENTER_HORIZONTAL, 5),
            (self.vertical_range_label,   wx.ALIGN_CENTER_HORIZONTAL, 5),
            spacing=5,
        )

        self.vertical_sizer = self.make_static_hbox(
            "Move \u2191/\u2193 (%):",
            self.vertical_slider,
            self.vertical_sizer2,
            spacing=5,
        )

    def choose_colour_panel(self):
        # Creating a combobox to change the colour of the 1D spectrum
        self.colour_label = wx.StaticBox(self, -1, "1D Line Colour")
        self.colour_sizer = wx.StaticBoxSizer(self.colour_label, wx.VERTICAL)
        self.options = colour_options
        self.colour_chooser = wx.ComboBox(
            self, value=colour_options[0], choices=colour_options, size=(100, height))
        self.colour_chooser.Bind(wx.EVT_COMBOBOX, self.OnColourChoice1D)
        self.colour_chooser.SetSelection(0)
        self.colour_sizer = self.make_static_vbox("1D Line Colour", self.colour_chooser)
    
    def change_linewidth_panel(self):

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

        
    def multiply_spec_panel(self):
        # Layout dimensions
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
            range_width  = int(leftover_width * 0.4)
            slider_width = int(leftover_width * 0.6)
        else:
            slider_width = leftover_width - 100
            range_width  = 100

        multiply_max = float(self.multiply_range_values[0])

        # Widgets
        self.multiply_value = 1

        self.multiply_slider = FloatSlider(
            self,
            id=-1,
            value=1,
            minval=0.1,
            maxval=multiply_max,
            res=multiply_max / 1000,
            size=(slider_width, height),
        )
        self.multiply_slider.Bind(wx.EVT_SLIDER, self.OnMultiplyScroll1D)

        self.multiply_range_chooser = wx.ComboBox(
            self,
            value=self.multiply_range_values[0],
            choices=self.multiply_range_values,
            size=(range_width, height),
        )
        self.multiply_range_chooser.Bind(wx.EVT_COMBOBOX, self.OnMultiplyCombo)
        self.multiply_range_chooser.SetSelection(0)

        self.multiply_label_value = wx.TextCtrl(
            self, value="1.000", size=(70, height), style=wx.TE_PROCESS_ENTER
        )
        self.multiply_label_value.Bind(wx.EVT_TEXT_ENTER, self.OnMultiplyText)

        self.multiply_combobox_label = wx.StaticText(self, label="Range")

        # Sizer assembly
        self.multiply_sizer_column1 = self.make_vbox(
            (self.multiply_slider,      wx.ALIGN_CENTER_HORIZONTAL, 5),
            (self.multiply_label_value, wx.ALIGN_CENTER_HORIZONTAL, 5),
            spacing=5,
        )

        self.multiply_sizer_column2 = self.make_vbox(
            (self.multiply_range_chooser,  wx.ALIGN_CENTER_HORIZONTAL, 5),
            (self.multiply_combobox_label, wx.ALIGN_CENTER_HORIZONTAL, 5),
            spacing=5,
        )

        self.multiply_sizer = self.make_hbox(
            self.multiply_sizer_column1,
            self.multiply_sizer_column2,
            spacing=5,
            outer=0,
        )

        self.multiply_total = self.make_static_vbox(
            "Multiplication Factor:",
            self.multiply_sizer,
            spacing=5,
        )

        # Pad to match linewidth sizer height if needed
        padding = self.linewidth_sizer.GetMinSize()[1] - self.multiply_total.GetMinSize()[1]
        if padding > 0:
            self.multiply_total.AddSpacer(padding)