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

# Frame for One-Dimensional NMR Spectra
class OneDViewer(wx.Panel):
    def __init__(self, parent, nmrdata, uc0=None, fid_viewer=False, title=''):
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
        self.fid_viewer=fid_viewer
        self.title=title
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
        self.main_sizer.Add(self.canvas, 1, wx.EXPAND)
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
        self.colours = colours
        self.colour_value = self.colours[0]

        # Range of the sliders to for moving spectra left/right/up/down
        self.reference_range_values = reference_range_values
        self.reference_range = float(self.reference_range_values[0])
        self.reference_rangeX = float(self.reference_range_values[0])
        self.reference_rangeY = float(self.reference_range_values[0])

        # Range of the sliders to for moving spectra up/down in 1D spectra
        self.vertical_range_values = vertical_range_values

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

        # A list to host spectra which are currently hidden
        self.hidden_list = []

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
            self.select_plot_label, choices=["Main Plot"], style=wx.CB_READONLY
        )
        self.plot_combobox.Bind(wx.EVT_COMBOBOX, self.OnSelectPlot)
        self.select_plot_sizer.Add(self.plot_combobox, 0, wx.ALL, 5)
        # Checkbox where can select all plots to be edited at the same time
        self.select_all_checkbox = wx.CheckBox(self.select_plot_label, label="Select All")
        self.hide_select_row = wx.BoxSizer(wx.HORIZONTAL)

        self.hide_checkbox = wx.CheckBox(self.select_plot_label, label="Hide")
        self.hide_checkbox.Bind(wx.EVT_CHECKBOX, self.OnHideSpectrum)

        self.hide_select_row.Add(self.hide_checkbox)
        self.hide_select_row.AddSpacer(5)
        self.hide_select_row.Add(self.select_all_checkbox)

        self.select_plot_sizer.Add(
            self.hide_select_row, 0, wx.ALIGN_CENTER_HORIZONTAL, 5
        )

        # Creating the phasing 1D sizer
        self.phasing_label = wx.StaticBox(self, -1, "Phasing:")
        self.phasing_sizer = wx.StaticBoxSizer(self.phasing_label, wx.VERTICAL)
        self.P0_label = wx.StaticText(self.phasing_label, label="P0 (Coarse):", size=(70, height))
        self.P1_label = wx.StaticText(self.phasing_label, label="P1 (Coarse):", size=(70, height))
        self.P0_slider = FloatSlider(
            self.phasing_label,
            id=-1,
            value=0,
            minval=-180,
            maxval=180,
            res=0.1,
            size=(int(self.parent.width / 5), height),
        )
        self.P1_slider = FloatSlider(
            self.phasing_label,
            id=-1,
            value=0,
            minval=-180,
            maxval=180,
            res=0.1,
            size=(int(self.parent.width / 5), height),
        )
        self.P0_slider.Bind(wx.EVT_SLIDER, self.OnSliderScroll1D)
        self.P1_slider.Bind(wx.EVT_SLIDER, self.OnSliderScroll1D)
        self.P0_label_fine = wx.StaticText(self.phasing_label, label="P0 (Fine):", size=(70, height))
        self.P1_label_fine = wx.StaticText(self.phasing_label, label="P1 (Fine):", size=(70, height))
        self.P0_slider_fine = FloatSlider(
            self.phasing_label,
            id=-1,
            value=0,
            minval=-10,
            maxval=10,
            res=0.01,
            size=(int(self.parent.width / 5), height),
        )
        self.P1_slider_fine = FloatSlider(
            self.phasing_label,
            id=-1,
            value=0,
            minval=-10,
            maxval=10,
            res=0.01,
            size=(int(self.parent.width / 5), height),
        )
        self.P0_slider_fine.Bind(wx.EVT_SLIDER, self.OnSliderScroll1D)
        self.P1_slider_fine.Bind(wx.EVT_SLIDER, self.OnSliderScroll1D)
        self.P0_total = wx.StaticText(self.phasing_label, label="P0 (Total):", size=(70, height))
        self.P1_total = wx.StaticText(self.phasing_label, label="P1 (Total):", size=(70, height))
        self.P0_total_value = wx.TextCtrl(self.phasing_label, value = "0", 
                                    size = (70, height), style = wx.TE_PROCESS_ENTER)
        self.P0_total_value.Bind(wx.EVT_TEXT_ENTER, self.P0_text_change)

        self.P1_total_value = wx.TextCtrl(self.phasing_label, value = "0", 
                                    size = (70,height), style = wx.TE_PROCESS_ENTER)
        self.P1_total_value.Bind(wx.EVT_TEXT_ENTER, self.P1_text_change)

        # Adding a button to change the range of the coarse and fine sliders (default to +/-180 and +/-10 degrees)
        self.update_phasing_range = wx.Button(self.phasing_label, label="Change slider range")
        self.update_phasing_range.Bind(wx.EVT_BUTTON, self.OnSliderRange1D)

        # Adding a button to set the pivot point for phasing
        self.pivot_button = wx.Button(self.phasing_label, label="Set Pivot Point")
        self.pivot_button.Bind(wx.EVT_BUTTON, self.OnPivotButton)
        self.pivot_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.pivot_sizer.Add(self.update_phasing_range)
        self.pivot_sizer.AddSpacer(10)
        self.pivot_sizer.Add(self.pivot_button)

        # Adding a button to remove the pivot point
        self.remove_pivot_button = wx.Button(self.phasing_label, label="Remove Pivot Point")
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
            self.contour_label, id=-1, value=0, minval=-1, maxval=10, res=0.01, size=(width, height)
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
            self.reference_label,
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
            self.reference_label,
            value=self.reference_range_values[0],
            choices=self.reference_range_values,
            size=(width, height),
        )
        self.reference_range_chooser.Bind(wx.EVT_COMBOBOX, self.OnReferenceCombo)
        self.reference_range_chooser.SetSelection(0)
        self.reference_sizer2.Add(
            self.reference_range_chooser, wx.ALIGN_CENTER_HORIZONTAL, 5
        )
        self.reference_value_label = wx.TextCtrl(self.reference_label, value = "0.0", 
                                    size = (70, height), style = wx.TE_PROCESS_ENTER)
        self.reference_value_label.Bind(wx.EVT_TEXT_ENTER, self.OnReferenceText)
        self.reference_sizer.AddSpacer(5)
        self.reference_sizer.Add(
            self.reference_value_label, wx.ALIGN_CENTER_HORIZONTAL, 5
        )
        self.reference_range_text = wx.StaticText(self.reference_label, label="Range")
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
            self.vertical_label,
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
            self.vertical_label,
            value=self.vertical_range_values[0],
            choices=self.vertical_range_values,
            size=(width, height),
        )
        self.vertical_range_chooser.Bind(wx.EVT_COMBOBOX, self.OnVerticalCombo)
        self.vertical_range_chooser.SetSelection(0)
        self.vertical_sizer2.Add(
            self.vertical_range_chooser, wx.ALIGN_CENTER_HORIZONTAL, 5
        )
        self.vertical_range_label = wx.StaticText(self.vertical_label, label="Range")
        self.vertical_sizer2.AddSpacer(5)
        self.vertical_sizer2.Add(
            self.vertical_range_label, wx.ALIGN_CENTER_HORIZONTAL, 5
        )

        self.vertical_sizer.Add(self.vertical_sizer2)

        # Creating a combobox to change the colour of the 1D spectrum
        self.colour_label = wx.StaticBox(self, -1, "Colour")
        self.colour_sizer = wx.StaticBoxSizer(self.colour_label, wx.VERTICAL)
        self.options = colour_options
        self.colour_chooser = wx.ComboBox(
            self.colour_label, value=self.options[0], choices=self.options, size=(75, height)
        )
        self.colour_chooser.Bind(wx.EVT_COMBOBOX, self.OnColourChoice1D)
        self.colour_chooser.SetSelection(0)
        spacer = 15
        self.colour_sizer.AddSpacer(spacer)
        self.colour_sizer.Add(self.colour_chooser)
        self.colour_sizer.AddSpacer(spacer)

        # Creating a slider to change the linewidth of the 1D spectrum
        self.linewidth_label = wx.StaticBox(self, -1, "Linewidth")
        self.linewidth_sizer = wx.StaticBoxSizer(self.linewidth_label, wx.VERTICAL)
        self.linewidth_slider = FloatSlider(
            self.linewidth_label, id=-1, value=0.5, minval=0.1, maxval=2, res=0.1, size=(50, height)
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
            self.multiply_label,
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
            self.multiply_label,
            value=self.multiply_range_values[0],
            choices=self.multiply_range_values,
            size=(range_width, height),
        )
        self.multiply_range_chooser.Bind(wx.EVT_COMBOBOX, self.OnMultiplyCombo)
        self.multiply_range_chooser.SetSelection(0)
        self.multiply_sizer_column2.Add(self.multiply_range_chooser, border=5)
        self.multiply_label_value = wx.TextCtrl(self.multiply_label, value = "1.000", 
                                    size = (70, height), style = wx.TE_PROCESS_ENTER)
        self.multiply_label_value.Bind(wx.EVT_TEXT_ENTER, self.OnMultiplyText)

        self.multiply_combobox_label = wx.StaticText(self.multiply_label, label="Range")
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

        # self.load_session_button =  wx.Button(self, label="Load Session", size=(130, 30))
        # self.load_session_button.Bind(wx.EVT_BUTTON, self.OnLoadSession)

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
        # self.button_sizers.Add(self.load_session_button)
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

    def OnHideSpectrum(self, event):
        if(self.multiplot_mode==False):
            if(self.hide_checkbox.GetValue()==True):
                self.hidden_list.append('Main Plot')
            else:
                self.hidden_list = []
        else:
            if(self.hide_checkbox.GetValue()==True):
                if(self.active_plot_index==0):
                    self.hidden_list.append(self.line1.get_label())
                else:
                    self.hidden_list.append(self.extra_plots[self.active_plot_index-1][0].get_label())
            else:
                if(self.active_plot_index==0):
                    self.hidden_list.remove(self.line1.get_label())
                else:
                    self.hidden_list.remove(self.extra_plots[self.active_plot_index-1][0].get_label())
                if(self.hidden_list==None):
                    self.hidden_list=[]
        
        self.hide_spectra()


    def hide_spectra(self):
        if(self.multiplot_mode==False):
            if(self.hide_checkbox.GetValue()==True):
                self.line1.set_visible(False)
            else:
                self.line1.set_visible(True)
        else:
            visible_lines = []
            if(self.line1.get_label() in self.hidden_list):
                self.line1.set_visible(False)
            else:
                self.line1.set_visible(True)
                visible_lines.append(self.line1)
            
            for i in range(len(self.extra_plots)):
                if(self.extra_plots[i][0].get_label() in self.hidden_list):
                    self.extra_plots[i][0].set_visible(False)
                else:
                    self.extra_plots[i][0].set_visible(True)
                    visible_lines.append(self.extra_plots[i][0])
            
            self.ax.legend(visible_lines, [l.get_label() for l in visible_lines])
        self.UpdateFrame()
    
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

    
    def OnLoadSession(self, event):
        self.parent.find_sessions(ask_user=False)
        session=self.parent.session_file
        self.parent.Destroy()
        from SpinExplorer.SpinView.SpinView import SpinView
        frame = SpinView(session_file=session)
        
       

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

        hidden=False
        if(self.active_plot_index==0):
            if(self.line1.get_label() in self.hidden_list):
                hidden=True
        elif(self.extra_plots[self.active_plot_index-1][0].get_label() in self.hidden_list):
            hidden=True
        else:
            pass
        self.hide_checkbox.SetValue(hidden)


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

        self.hide_spectra()

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

                reprocessing_frame = SpinProcess(self, path=path, cwd=self.parent.cwd, reprocess=True)
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

    def OnMultiplyText(self, event):
        # Function to multiply the 1D spectrum by a constant value
        self.multiply_value = float(self.multiply_label_value.GetValue())
        self.multiply_slider.SetValue(self.multiply_value)
        self.ApplyMultiplication()

    def OnMultiplyScroll1D(self, event):
        # Function to multiply the 1D spectrum by a constant value
        self.multiply_value = float(self.multiply_slider.GetValue())
        self.multiply_label_value.SetValue("{:.3f}".format(self.multiply_value))
        self.ApplyMultiplication()

    def ApplyMultiplication(self):
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


    def OnReferenceText(self, event):
        reference_value = float(self.reference_value_label.GetValue())
        self.reference_slider.SetValue(reference_value)
        self.ApplyReference(reference_value)


    def OnReferenceScroll1D(self, event):
        # Function to move the spectrum left/right in the ppm scale when the slider position is changed
        reference_value = float(self.reference_slider.GetValue())
        self.reference_value_label.SetValue("{:.4f}".format(reference_value))
        self.ApplyReference(reference_value)

    def ApplyReference(self, reference_value):
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
        


    def draw_figure_1D(self):
        # Function to plot the 1D spectrum
        self.ax = self.fig.add_subplot(111)
        self.key_press_connect = self.fig.canvas.mpl_connect(
            "key_press_event", self.on_key_1d
        )
        self.mouse_wheel_connect = self.fig.canvas.Bind(wx.EVT_MOUSEWHEEL, self.on_mouse_wheel)
        # Get ppm values for x axis
        if self.uc0 == None:
            if self.nmrdata.file != ".":
                self.uc0 = ng.pipe.make_uc(self.nmrdata.dic, self.nmrdata.data)
            else:
                udic = ng.bruker.guess_udic(self.nmrdata.dic, self.nmrdata.data)
                self.uc0 = ng.fileiobase.uc_from_udic(udic)


        if(self.fid_viewer==False):
            if(self.nmrdata.dic['FDF2FTFLAG']==1):
                # Fourier transformed axis
                self.ppm_original = self.uc0.ppm_scale()
                if('(ppm)' not in self.nmrdata.axislabels[0]):
                    self.nmrdata.axislabels[0] += ' (ppm)'
            else:
                # Data has not been Fourier transformed
                self.ppm_original = np.arange(0, len(self.nmrdata.data),1)
                if('(points)' not in self.nmrdata.axislabels[0]):
                    self.nmrdata.axislabels[0] += ' (points)'
        else:
            # showing the fid axis as points
            self.ppm_original = np.arange(0, len(self.nmrdata.data),1)
            if('(points)' not in self.nmrdata.axislabels[0]):
                self.nmrdata.axislabels[0] += ' (points)'
        self.ppms = self.ppm_original
        self.data = self.nmrdata.data
        (self.line1,) = self.ax.plot(self.ppms, self.data, linewidth=0.5)
        self.ax.set_xlabel(self.nmrdata.axislabels[0])
        self.ax.set_ylabel("Intensity")
        if('(points)' not in self.nmrdata.axislabels[0]):
            self.ax.set_xlim(max(self.ppms), min(self.ppms))
        else:
            self.ax.set_xlim(min(self.ppms), max(self.ppms))
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


    def OnSliderRange1D(self, event):
        """
        Creating a popout where a user can update the slider range
        """
        self.slider_range_window = PhasingSliderRange("Phasing slider ranges",self)


    def OnSliderScroll1D(self, event):
        # Get all the slider values for P0 and P1 (coarse and fine), put the combined coarse and fine values on the screen
        self.total_P0 = self.P0_slider.GetValue() + self.P0_slider_fine.GetValue()
        self.total_P1 = self.P1_slider.GetValue() + self.P1_slider_fine.GetValue()
        self.P0_total_value.SetValue("{:.2f}".format(self.total_P0))
        self.P1_total_value.SetValue("{:.2f}".format(self.total_P1))
        self.phase1D()

    def P0_text_change(self, event):
        self.total_P0 = float(self.P0_total_value.GetValue())
        self.P0_slider.SetValue(self.total_P0)
        self.P0_slider_fine.SetValue(0.0)
        self.total_P1 = self.P1_slider.GetValue() + self.P1_slider_fine.GetValue()
        self.phase1D()

    def P1_text_change(self, event):
        self.total_P1 = float(self.P1_total_value.GetValue())
        self.P1_slider.SetValue(self.total_P1)
        self.P1_slider_fine.SetValue(0.0)
        self.total_P0 = self.P0_slider.GetValue() + self.P0_slider_fine.GetValue()
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
