import wx # type: ignore
import numpy as np 
import nmrglue as ng # type: ignore
import sys 
import os
import copy
import matplotlib
matplotlib.use("wxAgg")
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import (
    NavigationToolbar2WxAgg as NavigationToolbar,
)
from SpinExplorer.SpinView.UI_objects.UI_tools import FloatSlider
from SpinExplorer.SpinView.Viewers.overlays import ReadProjection
from SpinExplorer.SpinView.Viewers.twod_view import TwoDViewer

from SpinExplorer.SpinView.Peaks.peaks import PeakListWindow3D
from SpinExplorer.SpinView.config import *


class ThreeDViewer(wx.Panel):
    def __init__(self, parent, nmrdata, fid_viewer=False):
        # Get the monitor size and set the window size to 85% of the monitor size
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        self.display_index = wx.Display.GetFromWindow(parent)
        self.width = int(1.0 * sizes[self.display_index][0])
        self.height = int(0.875 * sizes[self.display_index][1])
        self.parent = parent
        wx.Panel.__init__(self, parent, id=wx.ID_ANY, size=(self.width, self.height))
        self.nmrdata = nmrdata
        self.mouse_wheel_mode = ScrollMode.ZOOM

        self.fid_viewer = fid_viewer
        self.set_initial_variables_3D()
        self.create_button_panel_3D()
        self.create_hidden_button_panel_3D()
        self.create_canvas_3D()
        self.add_to_main_sizer_3D()


        if(self.fid_viewer==True):
            self.nmrdata.dic, self.nmrdata.data = ng.pipe_proc.tp(self.nmrdata.dic, self.nmrdata.data, auto=True)
            self.nmrdata.dic, self.nmrdata.data = ng.pipe_proc.di(self.nmrdata.dic,self.nmrdata.data)
            self.nmrdata.dic, self.nmrdata.data = self.transpose_3d(self.nmrdata.dic, self.nmrdata.data)
            self.nmrdata.dic, self.nmrdata.data = self.zero_transpose_3d(self.nmrdata.dic, self.nmrdata.data, unpack_complex=True)
            self.nmrdata.dic, self.nmrdata.data = ng.pipe_proc.di(self.nmrdata.dic,self.nmrdata.data)
            self.nmrdata.dic, self.nmrdata.data = self.zero_transpose_3d(self.nmrdata.dic, self.nmrdata.data)
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
        self.colours = colours
        self.colour_value = self.colours[0]

        # Initial 1D slice colour for 2D/3D spectra is set to navy
        self.colour_slice = "navy"

        # List of cmap colours for when overlaying multiple spectra
        self.cmap = "#e41a1c"
        self.cmap_neg = "#377eb8"
        self.twoD_colours = twoD_colours
        self.twoD_label_colours = self.twoD_colours

        self.twoD_slices_horizontal = []
        self.twoD_slices_vertical = []

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
        self.x_val = 10.00
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
        self.contour_val_box.AddSpacer(75)
        self.contour_val = wx.TextCtrl(
            self, value="10", size=(50, 20), style=wx.TE_PROCESS_ENTER
        )
        self.contour_val.Bind(wx.EVT_TEXT_ENTER, self.OnTextContour3D)
        self.contour_val_box.Add(self.contour_val)
        self.contour_sizer.Add(self.contour_val_box)

        # self.contour_sizer.AddSpacer(5)
        # self.contour_val_box = wx.BoxSizer(wx.HORIZONTAL)
        # self.contour_val_box.AddSpacer(125)
        # self.contour_val = wx.StaticText(self, label="10")
        # self.contour_val_box.Add(self.contour_val)
        # self.contour_sizer.Add(self.contour_val_box)

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
        from SpinExplorer.SpinProcess.SpinProcess import SpinProcess

        reprocessing_frame = SpinProcess(self, reprocess=True)
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

        frame = SpinBore(
            title="SpinBore - " + self.parent.title, projection=projection, parent=self
        )
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
            if(self.fid_viewer==False):
                self.ppms_0 = self.uc0.ppm_scale()
                self.ppms_1 = self.uc1.ppm_scale()
                self.ppms_2 = self.uc2.ppm_scale()
            else:
                self.ppms_0 = np.arange(0, len(self.uc0.ppm_scale()),1)
                self.ppms_1 = np.arange(0, len(self.uc1.ppm_scale()),1)
                self.ppms_2 = np.arange(0, len(self.uc2.ppm_scale()),1)


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
            if(self.fid_viewer==False):
                self.ppms_0 = self.uc0.ppm_scale()
                self.ppms_1 = self.uc1.ppm_scale()
                self.ppms_2 = self.uc2.ppm_scale()
            else:
                self.ppms_0 = np.arange(0, len(self.uc0.ppm_scale()),1)
                self.ppms_1 = np.arange(0, len(self.uc1.ppm_scale()),1)
                self.ppms_2 = np.arange(0, len(self.uc2.ppm_scale()),1)

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
            if(self.fid_viewer==False):
                self.ppms_0 = self.uc0.ppm_scale()
                self.ppms_1 = self.uc1.ppm_scale()
                self.ppms_2 = self.uc2.ppm_scale()
            else:
                self.ppms_0 = np.arange(0, len(self.uc0.ppm_scale()),1)
                self.ppms_1 = np.arange(0, len(self.uc1.ppm_scale()),1)
                self.ppms_2 = np.arange(0, len(self.uc2.ppm_scale()),1)

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
            if(self.fid_viewer==False):
                self.ppms_0 = self.uc0.ppm_scale()
                self.ppms_1 = self.uc1.ppm_scale()
                self.ppms_2 = self.uc2.ppm_scale()
            else:
                self.ppms_0 = np.arange(0, len(self.uc0.ppm_scale()),1)
                self.ppms_1 = np.arange(0, len(self.uc1.ppm_scale()),1)
                self.ppms_2 = np.arange(0, len(self.uc2.ppm_scale()),1)

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
            if(self.fid_viewer==False):
                self.ppms_0 = self.uc0.ppm_scale()
                self.ppms_1 = self.uc1.ppm_scale()
                self.ppms_2 = self.uc2.ppm_scale()
            else:
                self.ppms_0 = np.arange(0, len(self.uc0.ppm_scale()),1)
                self.ppms_1 = np.arange(0, len(self.uc1.ppm_scale()),1)
                self.ppms_2 = np.arange(0, len(self.uc2.ppm_scale()),1)

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

        if(self.fid_viewer==False):
            self.ax.set_xlim(max(self.ppms_0), min(self.ppms_0))
            self.ax.set_ylim(max(self.ppms_1), min(self.ppms_1))
        else:
            self.ax.set_xlim(min(self.ppms_0), max(self.ppms_0))
            self.ax.set_ylim(min(self.ppms_1), max(self.ppms_1))
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


    """
    Obtained from nmrglue followed by customisation
    
    Copyright Notice and Statement for the nmrglue Project
    Copyright (c) 2010-2015 Jonathan J. Helmus
    All rights reserved.
    """

    def transpose_3d(
        self, dic, data, hyper=False, nohyper=False, auto=False, nohdr=False
    ):
        """
        Transpose data (2D).

        Parameters
        ----------
        dic : dict
            Dictionary of NMRPipe parameters.
        data : ndarray
            Array of NMR data.
        hyper : bool
            True to perform hypercomplex transpose.
        nohyper : bool
            True to suppress hypercomplex transpose.
        auto : bool
            True to choose transpose mode automatically.
        nohdr : bool
            True to not update the transpose parameters in ndic.

        Returns
        -------
        ndic : dict
            Dictionary of updated NMRPipe parameters.
        ndata : ndarray
            Array of NMR data which has been transposed.

        """
        # XXX test if works with TPPI
        if nohyper:
            hyper = False

        fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
        fn2 = "FDF" + str(int(dic["FDDIMORDER"][1]))  # F1, F2, etc

        if auto:
            if (dic[fn + "QUADFLAG"] != 1) and (dic[fn2 + "QUADFLAG"] != 1):
                hyper = True
            else:
                hyper = False

        if hyper:  # Hypercomplex transpose need type recast
            data = np.array(ng.proc_base.tp_hyper(data), dtype="complex64")
        else:
            data = np.transpose(data, axes=(0, 2, 1))
            if dic[fn2 + "QUADFLAG"] != 1 and nohyper is False:
                # unpack complex as needed
                data = np.array(ng.proc_base.c2ri(data), dtype="complex64")

        # update the dimensionality and order
        dic["FDSLICECOUNT"] = data.shape[-2]
        if (data.dtype == "float32") and (nohyper is True):
            # when nohyper is True and the new last dimension was complex
            # prior to transposing then FDSIZE is set as if the dimension was
            # converted to complex data, that is half the actual size.
            dic["FDSIZE"] = data.shape[-1] / 2
        else:
            dic["FDSIZE"] = data.shape[-1]

        dic["FDSPECNUM"] = dic["FDSLICECOUNT"]
        dic["FDDIMORDER1"], dic["FDDIMORDER2"] = (
            dic["FDDIMORDER2"],
            dic["FDDIMORDER1"],
        )
        dic["FDDIMORDER"] = [
            dic["FDDIMORDER1"],
            dic["FDDIMORDER2"],
            dic["FDDIMORDER3"],
            dic["FDDIMORDER4"],
        ]

        if dic["FD2DPHASE"] == 0:
            dic["FDF1QUADFLAG"], dic["FDF2QUADFLAG"] = (
                dic["FDF2QUADFLAG"],
                dic["FDF1QUADFLAG"],
            )

        if nohdr is not True:
            dic["FDTRANSPOSED"] = (dic["FDTRANSPOSED"] + 1) % 2

        dic = ng.pipe_proc.clean_minmax(dic)
        return dic, data
    

    """
    Obtained from nmrglue followed by customisation
    
    Copyright Notice and Statement for the nmrglue Project
    Copyright (c) 2010-2015 Jonathan J. Helmus
    All rights reserved.
    """


    def zero_transpose_3d(self, dic, data, unpack_complex=False):
        """
        Transpose NMRPipe-style data from (X, Y, Z) to (Z, Y, X),
        including correct updates to the NMRPipe dictionary.

        Parameters:
            dic (dict): NMRPipe dictionary
            data (ndarray): NMRPipe data, assumed shape (X, Y, Z)

        Returns:
            new_dic (dict): Transposed dictionary
            new_data (ndarray): Transposed data, shape (Z, Y, X)
        """
        # Transpose data from (X, Y, Z) to (Z, Y, X)
        new_data = np.transpose(data, axes=(2, 1, 0))

        fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
        fn3 = "FDF" + str(int(dic["FDDIMORDER"][2]))  # F1, F2, etc

        # Create new dictionary
        new_dic = dic.copy()

        # for i, new_i in enumerate(
        #     new_axis_order
        # ):  # i = new dim index, new_i = old dim index
        #     for key in dic:
        #         if key.startswith(axis_keys[new_i]):
        #             # e.g., FDF1SW -> FDF1SW, becomes FDF3SW when i == 0 (Z)
        #             suffix = key[len(axis_keys[new_i]) :]  # e.g. SW, ORIG
        #             new_key = axis_keys[i] + suffix  # FDF1SW, FDF2SW, etc.
        #             new_dic[new_key] = dic[key]

        # swapping the FDDIMORDER1 and FDDIMORDER3 values
        order1 = dic["FDDIMORDER1"]
        order3 = dic["FDDIMORDER3"]
        new_dic["FDDIMORDER1"] = order3
        new_dic["FDDIMORDER3"] = order1
        new_dic["FDDIMORDER"][0] = order3
        new_dic["FDDIMORDER"][2] = order1

        new_dic["FDDIMORDER"] = [
            new_dic["FDDIMORDER1"],
            new_dic["FDDIMORDER2"],
            new_dic["FDDIMORDER3"],
            new_dic["FDDIMORDER4"],
        ]

        new_dic["FDSLICECOUNT"] = new_data.shape[-2]
        new_dic["FDSPECNUM"] = new_dic["FDSLICECOUNT"]
        new_dic["FDSIZE"] = new_data.shape[-1]

        if(unpack_complex==True):
            if dic[fn3 + "QUADFLAG"] != 1:
                # unpack complex as needed
                new_data = np.array(ng.proc_base.c2ri(new_data), dtype="complex64")
                new_dic[fn3 + "SIZE"] = int(new_dic[fn3 + "SIZE"] / 2)

        return new_dic, new_data

    def draw_figure_3D(self):
        self.ax = self.fig.add_subplot(111)
        self.axes1D = self.ax.twinx()
        self.axes1D_2 = self.ax.twiny()

        self.fig.canvas.mpl_connect("key_press_event", self.on_key_3d)
        self.fig.canvas.mpl_connect("button_press_event", self.on_click_3d)
        self.mouse_wheel_connect = self.fig.canvas.Bind(wx.EVT_MOUSEWHEEL, self.on_mouse_wheel)

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

        if(self.fid_viewer==False):
            self.ppms_0 = self.uc0.ppm_scale()
            self.ppms_1 = self.uc1.ppm_scale()
            self.ppms_2 = self.uc2.ppm_scale()
        else:
            self.ppms_0 = np.arange(0, len(self.uc0.ppm_scale()),1)
            self.ppms_1 = np.arange(0, len(self.uc1.ppm_scale()),1)
            self.ppms_2 = np.arange(0, len(self.uc2.ppm_scale()),1)

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

        if(self.fid_viewer==False):
            self.ax.set_xlim(max(self.ppms_0), min(self.ppms_0))
            self.ax.set_ylim(max(self.ppms_1), min(self.ppms_1))
        else:
            self.ax.set_xlim(min(self.ppms_0), max(self.ppms_0))
            self.ax.set_ylim(min(self.ppms_1), max(self.ppms_1))

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
        self.projection_window = ProjectionFrame(
            parent=self, title="Projections - " + self.parent.title
        )

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
            self.plot3D_window = Plot3DFrame(
                parent=self, title="3D Plot - " + self.parent.title
            )
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
                parent=self,
                title="Waterfall Plot - " + self.parent.title,
                visible=visible,
            )
            self.waterfall_warning.Destroy()
        else:
            self.waterfall_warning.Destroy()
            return

    def OnTextContour3D(self,event):

        try:
            self.contour_slider.SetValue(np.log10(float(self.contour_val.GetValue())))
            self.x_val = self.contour_val.GetValue()
            self.OnMinContour3D(event)
        except:
            self.OnMinContour3D(event)

    def OnLinewidthScroll3D(self, event):
        self.contour_linewidth = float(self.linewidth_slider.GetValue())
        self.OnMinContour3D(event)

    def mouse_wheel_zoom(self, event):
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

    def on_mouse_wheel(self, event):

        toolbar = self.fig.canvas.toolbar
        if toolbar:
            toolbar.push_current() # logs position in toolbar so commands back, forward, home work

        if self.mouse_wheel_mode == ScrollMode.ZOOM:
            self.mouse_wheel_zoom(event)
        
        if self.mouse_wheel_mode == ScrollMode.CONTOUR:
            delta = 0.1 if event.GetWheelRotation() > 0 else -0.1
            current = float(self.contour_slider.GetValue())
            self.contour_slider.SetValue(current + delta)
            self.x_val = 10 ** float(self.contour_slider.GetValue())
            self.contour_val.SetValue(
                "{:.2f}".format(10 ** float(self.contour_slider.GetValue()))
            )
            self.OnMinContour3D(event)
        
        if self.mouse_wheel_mode == ScrollMode.PLANE:
            delta = 0.1 if event.GetWheelRotation() > 0 else -0.1 
            current = float(self.z_slider.GetValue())
            self.z_slider.SetValue(current+delta)
            self.OnZScroll3D(event)

    def on_key_3d(self, event):
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
        
        if event.key == "c":
            self.mouse_wheel_mode = cycle_scroll_mode(self.mouse_wheel_mode, THREED_SCROLL_MODES)

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
                    if(self.fid_viewer==False):
                        data = self.nmrdata.data[z_index][
                            :, self.uc1(str(self.ppms_1[1]) + "ppm")
                        ]
                    else:
                        data = self.nmrdata.data[z_index][
                            :, int(self.ppms_1[1])
                        ]
                    (self.line1,) = self.axes1D.plot(
                        self.ppms_0,
                        data,
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
                    if(self.fid_viewer==False):
                        data = self.nmrdata.data[z_index][
                            self.uc0(str(self.ppms_0[1]) + "ppm"), :
                        ]
                    else:
                        data = self.nmrdata.data[z_index][
                            int(self.ppms_0[1]), :
                        ]
                    (self.line3,) = self.axes1D_2.plot(
                        data,
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
                if(self.fid_viewer==False):
                    data = self.nmrdata.data[z_index][
                        :, self.uc1(str(self.y1 - self.y_movement) + "ppm")
                    ]
                else:
                    data = self.nmrdata.data[z_index][
                        :, int(self.y1 - self.y_movement)
                    ]
                self.line1.set_ydata(
                    data
                )
                self.line2.set_ydata([self.y1])
                self.line1.set_xdata(self.new_x_ppms)
                self.OnSliderScroll3D(None)
            if self.line3.get_visible() == True:
                if(self.fid_viewer==False):
                    data = self.nmrdata.data[z_index][
                        self.uc0(str(self.x1 - self.x_movement) + "ppm"), :
                    ]
                else:
                    data = self.nmrdata.data[z_index][
                        int(self.x1 - self.x_movement), :
                    ]
                self.line3.set_xdata(
                    data
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
            if(self.fid_viewer==False):
                data = self.nmrdata.data[z_index][:, self.uc1(str(self.y1) + "ppm")]
            else:
                data = self.nmrdata.data[z_index][:, int(self.y1)]
            complex_data = ng.process.proc_base.ht(data, self.nmrdata.data.shape[1])
            self.phased_data = ng.process.proc_base.ps(
                complex_data, p0=self.total_P0, p1=self.total_P1
            )
            self.line1.set_ydata(self.phased_data)
        if self.line3.get_visible() == True:
            if(self.fid_viewer==False):
                data = self.nmrdata.data[z_index][self.uc0(str(self.x1) + "ppm"), :]
            else:
                data = self.nmrdata.data[z_index][int(self.x1), :]
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

    # def OnContourSlider(self, event):
    #     contour_val = 10 ** float(self.contour_slider.GetValue())
    #     contour_start = (
    #         np.max(self.main_frame.nmrdata.data) / contour_val
    #     )  # contour level start value
    #     self.contour_num = 20  # number of contour levels
    #     self.contour_factor = 1.2  # scaling factor between contour levels
    #     # calculate contour levels
    #     self.cl = contour_start * self.contour_factor ** np.arange(self.contour_num)
    #     self.cl_neg = -contour_start * self.contour_factor ** np.flip(
    #         np.arange(self.contour_num)
    #     )

    #     xlim = self.ax.get_xlim3d()
    #     ylim = self.ax.get_ylim3d()
    #     zlim = self.ax.get_zlim3d()
    #     xlabel = self.ax.get_xlabel()
    #     ylabel = self.ax.get_ylabel()
    #     zlabel = self.ax.get_zlabel()
    #     view = self.ax.azim

    #     self.ax.clear()

    #     for i in range(self.main_frame.nmrdata.data.shape[0]):
    #         x = self.main_frame.ppms_1
    #         y = self.main_frame.ppms_0
    #         z = self.main_frame.nmrdata.data[i]
    #         x, y = np.meshgrid(x, y)
    #         self.ax.contour(
    #             x,
    #             y,
    #             z,
    #             zdir="z",
    #             offset=self.main_frame.ppms_2[i],
    #             levels=self.cl,
    #             colors="red",
    #             linewidths=0.5,
    #         )

    #     self.ax.set_zlim3d(zlim[0], zlim[1])
    #     self.ax.set_xlim3d(xlim[0], xlim[1])
    #     self.ax.set_ylim3d(ylim[0], ylim[1])
    #     self.ax.set_xlabel(xlabel)
    #     self.ax.set_ylabel(ylabel)
    #     self.ax.set_zlabel(zlabel)
    #     self.ax.view_init(azim=view)

    #     self.Update3DFrame()

    # def DrawContours3D(self):
    #     """Core drawing function — call this from any trigger."""
    #     contour_start = np.max(self.main_frame.nmrdata.data) / self.x_val
    #     self.contour_num = 20
    #     self.contour_factor = 1.2
    #     self.cl = contour_start * self.contour_factor ** np.arange(self.contour_num)
    #     self.cl_neg = -contour_start * self.contour_factor ** np.flip(np.arange(self.contour_num))

    #     xlim, ylim, zlim = self.ax.get_xlim3d(), self.ax.get_ylim3d(), self.ax.get_zlim3d()
    #     xlabel, ylabel, zlabel = self.ax.get_xlabel(), self.ax.get_ylabel(), self.ax.get_zlabel()
    #     view = self.ax.azim

    #     self.ax.clear()

    #     for i in range(self.main_frame.nmrdata.data.shape[0]):
    #         x, y = np.meshgrid(self.main_frame.ppms_1, self.main_frame.ppms_0)
    #         self.ax.contour(
    #             x, y,
    #             self.main_frame.nmrdata.data[i],
    #             zdir="z",
    #             offset=self.main_frame.ppms_2[i],
    #             levels=self.cl,
    #             colors="red",
    #             linewidths=0.5,
    #         )

    #     self.ax.set_xlim3d(xlim[0], xlim[1])
    #     self.ax.set_ylim3d(ylim[0], ylim[1])
    #     self.ax.set_zlim3d(zlim[0], zlim[1])
    #     self.ax.set_xlabel(xlabel)
    #     self.ax.set_ylabel(ylabel)
    #     self.ax.set_zlabel(zlabel)
    #     self.ax.view_init(azim=view)

    #     self.Update3DFrame()


    # def OnContourSlider(self, event):
    #     """Triggered by slider."""
    #     self.x_val = 10 ** float(self.contour_slider.GetValue())
    #     self.contour_value_label.SetValue("{:.2f}".format(self.x_val))
    #     self.DrawContours3D()


    # def OnContourText(self, event):
    #     """Triggered by text control."""
    #     self.x_val = float(self.contour_value_label.GetValue())
    #     self.DrawContours3D()


    # def OnMouseWheel(self, event):
    #     match self.mouse_wheel_mode:
    #         case ScrollMode.ZOOM:
    #             self.OnZoom(event)
    #         case ScrollMode.CONTOUR:
    #             delta = 0.1 if event.GetWheelRotation() > 0 else -0.1
    #             current = float(self.contour_slider.GetValue())
    #             self.contour_slider.SetValue(current + delta)
    #             self.x_val = 10 ** float(self.contour_slider.GetValue())
    #             self.contour_value_label.SetValue("{:.2f}".format(self.x_val))
    #             self.DrawContours3D()
    #         case ScrollMode.PLANE:
    #             self.OnChangePlane(event)


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
        # Checking if the projection data needs transposing to match the main frame
        self.check_for_transpose()

        self.sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.titlecolor = "black"

        self.selected_bore_peaks = []
        self.free_protein = 'Free'
        self.nucleus = 'H'

        # Plot a cross for where the user has clicked on the 2D projection
        self.plot_cross = True

        # A flag to describe whether the strip plot orientation is reversed. This is necessary when overlaying peaks on top of the data
        self.alternative_orientation = False

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
            self, -1, choices=["1H", "13C", "15N"], style=wx.CB_READONLY
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


        self.read_peaks_button = wx.Button(self, -1, "Read 3D Peaks")
        self.read_peaks_button.Bind(wx.EVT_BUTTON, self.OnReadPeaksBore)

        self.bore_sizer.Add(self.bore_sizer_2D)
        self.bore_sizer.AddSpacer(10)
        self.bore_sizer.Add(self.bore_sizer_1D)
        self.bore_sizer_row2 = wx.BoxSizer(wx.HORIZONTAL)
        self.bore_sizer_row2.AddSpacer(10)
        self.bore_sizer_row2.Add(self.bore_sizer_strip)
        self.bore_sizer_row2.AddSpacer(10)
        self.bore_sizer_row2.Add(self.read_peaks_button, 0, wx.ALIGN_CENTER_VERTICAL)

        self.main_bore_sizer.Add(self.bore_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.main_bore_sizer.Add(self.bore_sizer_row2, 0, wx.ALIGN_CENTER_HORIZONTAL)


    def OnReadPeaksBore(self, event):
        """
        Open up the Peak Lists 3D frame
        """

        self.peak_lists3D = PeakListWindow3D(
            title="3D Peak List - " + self.main_frame.parent.title, parent=self
        )
        self.peak_lists3D.Show()


    def check_for_transpose(self):
        """
        Check to see if the projection needs transposing in order to
        match that in the main spectrum x/y axes etc.

        For example: if the current 3D is showing (15N, 15N_1), 1H and the
        15N_1.15N.dat projection file is being read then the projection data
        needs to be transposed after reading (without setting transposed equal
        to true)
        """

        main_frame_axes = self.main_frame.orientation_chooser.GetValue()[1:].split(')')[0].split(',')
        projection_axes = self.nmrdata.axislabels

        self.swap_labels = False

        if(main_frame_axes!=projection_axes):
            # transpose the data
            self.nmrdata.dic, self.nmrdata.data = ng.pipe_proc.tp(copy.deepcopy(self.nmrdata.dic), copy.deepcopy(self.nmrdata.data))
            self.swap_labels = True


    def plot_bore_data(self):
        # Make a figure containing 2 plots, one large 2D contour plot and a vertical smaller plot showing the bore down a selected 2D coordinate
        self.ax_bore, self.ax_bore_2, self.ax_bore_3 = self.fig_bore.subplots(
            1, 3, gridspec_kw={"width_ratios": [3, 1, 1]}
        )

        self.ax_bore_2.sharey(self.ax_bore_3)

        self.click_press_connect = self.fig_bore.canvas.mpl_connect(
            "button_press_event", self.on_click_bore
        )

        self.cmap = "#e41a1c"
        self.cmap_neg = "#377eb8"
        self.transposed2D = False

        contour_start = np.max(self.nmrdata.data) / 10  # contour level start value
        self.contour_num = 50  # number of contour levels
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
        if(self.swap_labels==False):
            self.ax_bore.set_xlabel(self.nmrdata.axislabels[1])
            self.ax_bore.set_ylabel(self.nmrdata.axislabels[0])
        else:
            self.ax_bore.set_xlabel(self.nmrdata.axislabels[0])
            self.ax_bore.set_ylabel(self.nmrdata.axislabels[1])
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
        self.ax_bore_3.set_xlim(max(self.ppms_0), min(self.ppms_0))
        self.ax_bore_3.set_ylim(max(self.ppms_2), min(self.ppms_2))
        if(self.swap_labels==False):
            self.ax_bore_3.set_xlabel(self.nmrdata.axislabels[1])
        else:
            self.ax_bore_3.set_xlabel(self.nmrdata.axislabels[0])
        if self.Xstrip.shape != self.bore_data_strip1.shape:
            self.alternative_orientation = True
            self.Xstrip, self.Ystrip = np.meshgrid(self.ppms_1, self.ppms_2)
            self.ax_bore_3.set_xlim(max(self.ppms_1), min(self.ppms_1))
            self.ax_bore_3.set_ylim(max(self.ppms_2), min(self.ppms_2))
            if(self.swap_labels==False):
                self.ax_bore_3.set_xlabel(self.nmrdata.axislabels[0])
            else:
                self.ax_bore_3.set_xlabel(self.nmrdata.axislabels[1])
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

        self.ax_bore_3.set_title("Strip Plot")

        self.original_limits = [
            (ax.get_xlim(), ax.get_ylim()) for ax in [self.ax_bore_2, self.ax_bore_3]
        ]

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

                title = self.ax_bore_3.get_title()
                xlim3, ylim3 = self.ax_bore_3.get_xlim(), self.ax_bore_3.get_ylim()
                self.ax_bore_3.clear()

                self.Xstrip, self.Ystrip = np.meshgrid(self.ppms_0, self.ppms_2)
                if(self.swap_labels==False):
                    self.ax_bore_3.set_xlabel(self.nmrdata.axislabels[1])
                else:
                    self.ax_bore_3.set_xlabel(self.nmrdata.axislabels[0])
                if self.Xstrip.shape != self.bore_data_strip1.shape:
                    self.Xstrip, self.Ystrip = np.meshgrid(self.ppms_1, self.ppms_2)
                    self.ax_bore_3.set_xlim(max(self.ppms_1), min(self.ppms_1))
                    self.ax_bore_3.set_ylim(max(self.ppms_2), min(self.ppms_2))
                    if(self.swap_labels==False):
                        self.ax_bore_3.set_xlabel(self.nmrdata.axislabels[0])
                    else:
                        self.ax_bore_3.set_xlabel(self.nmrdata.axislabels[1])
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
                self.ax_bore_3.set_xlim(xlim3)
                self.ax_bore_3.set_title(title)
                self.ax_bore_3.set_ylim(ylim3)

                # for window in wx.GetTopLevelWindows():
                #     if (
                #         isinstance(window, wx.Frame)
                #         and window.GetTitle().split()[0] == "3D"
                #         and window.GetTitle().split()[1] == "Peak"
                #     ):
                #         if(self.selected_bore_peaks!=[]):
                #             # Plot these bore peaks
                #             xvals = []
                #             yvals = []
                #             names = []
                #             for index in self.selected_bore_peaks:
                #                 names.append(self.peak_lists3D.peak_list_dictionary[self.peak_lists3D.peak_list_choices[0]]['peak_names'][index])
                #                 xvals.append(self.peak_lists3D.peak_list_dictionary[self.peak_lists3D.peak_list_choices[0]]['shift1'][index])
                #                 yvals.append(self.peak_lists3D.peak_list_dictionary[self.peak_lists3D.peak_list_choices[0]]['shift3'][index])

                #             self.scatter_strip = self.ax_bore_3.scatter(xvals, yvals, s=5,
                #             marker="o",
                #             picker=5,
                #             zorder=2)

                #     # Annotation for hover
                #             self.annotations_strip = self.ax_bore_3.annotate(
                #                 "",
                #                 xy=(0, 0),
                #                 xytext=(15, 15),
                #                 textcoords="offset points",
                #                 bbox=dict(boxstyle="round", fc="w"),
                #                 arrowprops=dict(arrowstyle="->"))
                #             self.annotations[-1].set_visible(False)

                #             # Connect event
                #             self.hover_connect_strip = self.canvas_bore.mpl_connect(
                #                 "motion_notify_event", self.on_hover_strip
                #             )

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
                if(self.swap_labels==False):
                    self.ax_bore_3.set_xlabel(self.nmrdata.axislabels[1])
                else:
                    self.ax_bore_3.set_xlabel(self.nmrdata.axislabels[0])
                if self.Xstrip.shape != self.bore_data_strip1.shape:
                    self.Xstrip, self.Ystrip = np.meshgrid(self.ppms_1, self.ppms_2)
                    self.ax_bore_3.set_xlim(max(self.ppms_1), min(self.ppms_1))
                    self.ax_bore_3.set_ylim(max(self.ppms_2), min(self.ppms_2))
                    if(self.swap_labels==False):
                        self.ax_bore_3.set_xlabel(self.nmrdata.axislabels[0])
                    else:
                        self.ax_bore_3.set_xlabel(self.nmrdata.axislabels[1])

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
                self.ax_bore_3.set_xlim(xlim3)
                self.ax_bore_3.set_title(title)
                self.ax_bore_3.set_ylim(ylim3)

            self.OverlayBore()
            self.ax_bore_2.set_ylim(self.original_limits[0][1])
            self.ax_bore_3.set_xlim(self.original_limits[1][0])
            self.ax_bore_3.set_ylim(self.original_limits[1][1])
        # try:
        #     self.overlay_peaklist()
        # except:
        #     pass
        # self.canvas_bore.draw_idle()
        self.UpdateBoreFrame()

    def overlay_peaklist(self):
        """
        This will show all peaks down the bore that have been selected in
        the peaklist
        """

        for window in wx.GetTopLevelWindows():
            if (
                isinstance(window, wx.Frame)
                and window.GetTitle().split()[0] == "3D"
                and window.GetTitle().split()[1] == "Peak"
            ):
                if len(self.selected_bore_peaks) > 0:
                    # Plot these bore peaks
                    xvals = []
                    yvals = []
                    names = []
                    colors = []
                    for index in self.selected_bore_peaks:
                        names.append(
                            self.peak_lists3D.peak_list_dictionary[
                                self.peak_lists3D.peak_list_choices[0]
                            ]["peak_name"][index]
                        )
                        if('N/A' not in self.peak_lists3D.selected_peak_indexes):
                            if(index in self.peak_lists3D.selected_peak_indexes):
                                colors.append('darkviolet')
                            else:
                                colors.append('k')
                        else:
                            colors.append('k')
                

                        s = self.peak_lists3D.bore_xdim 
                        if(self.alternative_orientation == True):
                            if(s == 'shift1'):
                                s = 'shift2'
                            else:
                                s = 'shift1'
                        
                        xvals.append(
                            self.peak_lists3D.peak_list_dictionary[
                                self.peak_lists3D.peak_list_choices[0]
                            ][s][index]
                        )
                        yvals.append(
                            self.peak_lists3D.peak_list_dictionary[
                                self.peak_lists3D.peak_list_choices[0]
                            ]["shift3"][index]
                        )



                    self.scatter_strip = self.ax_bore_3.scatter(
                        xvals, yvals, s=5, marker="o", picker=5, zorder=2, color=colors
                    )

                    # Annotation for hover
                    self.annotations_strip = self.ax_bore_3.annotate(
                        "",
                        xy=(0, 0),
                        xytext=(15, 15),
                        textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"),
                    )
                    self.annotations[-1].set_visible(False)

                    # Connect event
                    self.hover_connect_strip = self.canvas_bore.mpl_connect(
                        "motion_notify_event", self.on_hover_strip
                    )

    def OverlayBore(self):
        if self.include_overlay == True:
            intensity_percent = 10 ** float(self.bore_intensity_slider.GetValue())
            # Plot current amino acid selection
            if self.bore_free_protein_combo_box.GetValue() == "Free":
                data = self.bmrb_free[self.amino_acid][self.nucleus][0]
                labels = self.bmrb_free[self.amino_acid][self.nucleus][1]
            else:
                data = self.bmrb_protein[self.amino_acid][self.nucleus][0]
                labels = self.bmrb_protein[self.amino_acid][self.nucleus][1]

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

        try:
            xvalue = self.line3.get_ydata()
        except:
            xvalue = "1"

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
        self.ax_bore_3.set_xlabel(xlabel)
        self.ax_bore_3.set_title(title)

        if xvalue != "1":
            if self.transposed2D == False:
                xvalue = self.bore_initial[0]
            else:
                xvalue = self.bore_initial[1]

            self.line3 = self.ax_bore_3.axvline(x=xvalue, color="black", linewidth=0.5)

        self.overlay_peaklist()

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
        if(self.plot_cross == True):
            (self.cross,) = self.ax_bore.plot(
                cross_data[0], cross_data[1], marker="X", color="k"
            )

        self.ax_bore.set_xlim(xlim)
        self.ax_bore.set_ylim(ylim)
        if(self.swap_labels==False):
            self.ax_bore.set_xlabel(self.nmrdata.axislabels[1])
            self.ax_bore.set_ylabel(self.nmrdata.axislabels[0])
        else:
            self.ax_bore.set_xlabel(self.nmrdata.axislabels[0])
            self.ax_bore.set_ylabel(self.nmrdata.axislabels[1])

        self.add_peaklist()

        self.OnBoreSliderStripPlot(wx.EVT_BUTTON)

        self.UpdateBoreFrame()

    def transpose_peaklist(self):
        """
        Transpose the projection peaklist when the 2D plot transpose button
        is pressed
        """

        for window in wx.GetTopLevelWindows():
            if (
                isinstance(window, wx.Frame)
                and window.GetTitle().split()[0] == "3D"
                and window.GetTitle().split()[1] == "Peak"
            ):
                for (
                    peaklist_name,
                    dictionary,
                ) in self.peak_lists3D.peak_list_dictionary.items():
                    shift1 = dictionary["shift1"]
                    shift2 = dictionary["shift2"]
                    self.peak_lists3D.peak_list_dictionary[peaklist_name][
                        "shift1"
                    ] = shift2
                    self.peak_lists3D.peak_list_dictionary[peaklist_name][
                        "shift2"
                    ] = shift1

                if(self.peak_lists3D.bore_xdim == 'shift1'):
                    self.peak_lists3D.bore_xdim = 'shift2'
                else:
                    self.peak_lists3D.bore_xdim = 'shift1'
                

    def add_peaklist(self):
        """
        Adding a peaklist to the projection plot
        """

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
            if (
                isinstance(window, wx.Frame)
                and window.GetTitle().split()[0] == "3D"
                and window.GetTitle().split()[1] == "Peak"
            ):
                # Plot Peaklists
                count = 0
                self.points = []
                self.annotations = []
                for (
                    peaklist_name,
                    dictionary,
                ) in self.peak_lists3D.peak_list_dictionary.items():
                    
                    if (
                        self.peak_lists3D.select_peak_button.GetValue() == True
                    ):
                                cs = []
                                for i, peak in enumerate(dictionary["peak_name"]):
                                    if i in self.peak_lists3D.selected_peak_indexes:
                                        cs.append("darkviolet")
                                    else:
                                        cs.append(self.peaklist_colours[count])
                    else:
                        cs = self.peaklist_colours[count]

                    # cs = self.peaklist_colours[count]

                    shift1 = dictionary["shift1"]
                    shift2 = dictionary["shift2"]
                    self.points.append(
                        self.ax_bore.scatter(
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
                        self.ax_bore.annotate(
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
                    self.hover_connect = self.canvas_bore.mpl_connect(
                        "motion_notify_event", self.on_hover
                    )
                    self.click_connect = self.canvas_bore.mpl_connect(
                        "button_press_event", self.on_pick
                    )

            else:
                pass

            self.UpdateBoreFrame()

    def on_hover(self, event):

        if event.inaxes is self.ax_bore:
            # Calculate distance from mouse to each point

            cont, ind = self.points[0].contains(event)
            if cont:
                # Show annotation
                index = ind["ind"][0]  # first index found
                peaklist_name = self.peak_lists3D.peak_list_choices[0]
                dictionary = self.peak_lists3D.peak_list_dictionary[peaklist_name]
                peakname = dictionary["peak_name"][index] + " (" + peaklist_name + ")"
                x = dictionary["shift1"][index]
                y = dictionary["shift2"][index]
                self.annotations[0].xy = (x, y)
                text = peakname
                self.annotations[0].set_text(text)
                self.annotations[0].set_color(self.peaklist_colours[0])
                self.annotations[0].set_position((36, 0))
                self.annotations[0].set_visible(True)
                self.canvas_bore.draw_idle()
            else:
                if self.annotations[0].get_visible():
                    self.annotations[0].set_visible(False)

            # adjust_text(self.annotations, ax=self.ax, time_lim=5)
            self.canvas_bore.draw_idle()

    def on_pick(self, event):
        if event.inaxes is self.ax_bore:
            contains, details = self.points[0].contains(event)
            try:
                self.scatter_strip.clear()
                self.annotations_strip.clear()
            except:
                pass
            if contains:
                indices = details["ind"]
                self.selected_bore_peaks = indices
                self.overlay_peaklist()
                self.UpdateBoreFrame()
            else:
                self.selected_bore_peaks = []

    def on_hover_strip(self, event):

        if event.inaxes is self.ax_bore_3:
            # Calculate distance from mouse to each point

            cont, ind = self.scatter_strip.contains(event)
            if cont:
                # Show annotation
                index = ind["ind"][0]  # first index found
                peaklist_name = self.peak_lists3D.peak_list_choices[0]
                dictionary = self.peak_lists3D.peak_list_dictionary[peaklist_name]
                peakname = (
                    dictionary["peak_name"][self.selected_bore_peaks[index]]
                    + " ("
                    + peaklist_name
                    + ")"
                )

                s = self.peak_lists3D.bore_xdim 
                if(self.alternative_orientation == True):
                    if(s == 'shift1'):
                        s = 'shift2'
                    else:
                        s = 'shift1'
                
                x = dictionary[s][self.selected_bore_peaks[index]]
                y = dictionary["shift3"][self.selected_bore_peaks[index]]
                self.annotations_strip.xy = (x, y)
                text = peakname
                self.annotations_strip.set_text(text)
                self.annotations_strip.set_color(self.peaklist_colours[0])
                self.annotations_strip.set_position((36, 0))
                self.annotations_strip.set_visible(True)
                self.canvas_bore.draw_idle()
            else:
                if self.annotations_strip.get_visible():
                    self.annotations_strip.set_visible(False)

            # adjust_text(self.annotations, ax=self.ax, time_lim=5)
            self.canvas_bore.draw_idle()

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
        if(self.plot_cross == True):
            (self.cross,) = self.ax_bore.plot(
                cross_data[1], cross_data[0], marker="X", color="k"
            )
        self.transpose_peaklist()
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
        self.OnFreeProteinSelection(event, silent=True)
        if self.bore_combo_box.GetSelection() == 0:
            self.nucleus = "H"
        elif self.bore_combo_box.GetSelection() == 1:
            self.nucleus = "C"
        else:
            if(self.free_protein == 'Protein'):
                self.nucleus = "N"
            else:
                # Nitrogen shifts are not available for the free amino acids
                message = 'Nitrogen shifts are not available for the free amino acids'
                dlg = wx.MessageBox(message, "Nucleus selection", wx.OK)
                
                self.bore_combo_box.SetSelection(1)
                self.nucleus = "C"
        
        self.OnFreeProteinSelection(event, silent=False)
        

    def OnFreeProteinSelection(self, event, silent=False):
        if self.bore_free_protein_combo_box.GetSelection() == 0:
            self.free_protein = "Free"
        else:
            self.free_protein = "Protein"

        if(silent==False):
            self.OverlayBore()

    def read_bmrb(self):

        self.read_free()
        self.read_protein()
        # else:
        #     # Give a warning saying that the 'bmrb_free.txt' file is missing
        #     message = "The file 'bmrb_free.txt' is missing from the SpinView directory. Unable to show free amino acid overlays."
        #     dlg = wx.MessageDialog(None, message, 'File Missing', wx.OK | wx.ICON_WARNING)
        #     dlg.ShowModal()
        #     dlg.Destroy()
        #     self.include_overlay = False
        #     return

    def read_protein(self):
        self.bmrb_protein = {}
        data = """
ALA,H,8.194
ALA,HA,4.237
ALA,HB,1.357
ALA,CO,177.815
ALA,CA,53.138
ALA,CB,18.956
ALA,N,123.380
ARG,H,8.235
ARG,HA,4.285
ARG,HB2,1.795
ARG,HB3,1.765
ARG,HD2,3.118
ARG,HD3,3.103
ARG,HE,7.355
ARG,HG2,1.567
ARG,HG3,1.547
ARG,HH11,6.895
ARG,HH12,6.845
ARG,HH21,6.807
ARG,HH22,6.796
ARG,CO,176.491
ARG,CA,56.765
ARG,CB,30.623
ARG,CD,43.163
ARG,CG,27.216
ARG,CZ,159.864
ARG,N,120.904
ARG,NE,84.593
ARG,NH1,74.089
ARG,NH2,72.837
ASN,H,8.322
ASN,HA,4.658
ASN,HB2,2.803
ASN,HB3,2.749
ASN,HD21,7.328
ASN,HD22,7.148
ASN,CO,175.290
ASN,CA,53.517
ASN,CB,38.681
ASN,CG,176.714
ASN,N,118.974
ASN,ND2,112.749
ASP,H,8.295
ASP,HA,4.581
ASP,HB2,2.711
ASP,HB3,2.661
ASP,HD2,6.567
ASP,CO,176.427
ASP,CA,54.664
ASP,CB,40.863
ASP,CG,179.232
ASP,N,120.744
CYS,H,8.380
CYS,HA,4.655
CYS,HB2,2.956
CYS,HB3,2.893
CYS,HG,2.086
CYS,CO,174.820
CYS,CA,58.035
CYS,CB,33.441
CYS,N,120.092
GLN,H,8.219
GLN,HA,4.256
GLN,HB2,2.045
GLN,HB3,2.016
GLN,HE21,7.229
GLN,HE22,7.036
GLN,HG2,2.316
GLN,HG3,2.296
GLN,CO,176.345
GLN,CA,56.532
GLN,CB,29.144
GLN,CD,179.705
GLN,CG,33.785
GLN,N,120.091
GLN,NE2,111.870
GLU,H,8.330
GLU,HA,4.238
GLU,HB2,2.021
GLU,HB3,1.998
GLU,HE2,4.132
GLU,HG2,2.268
GLU,HG3,2.250
GLU,CO,176.930
GLU,CA,57.299
GLU,CB,29.946
GLU,CD,182.205
GLU,CG,36.110
GLU,N,120.802
GLY,H,8.328
GLY,H1,8.525
GLY,HA2,3.958
GLY,HA3,3.894
GLY,CO,173.899
GLY,CA,45.347
GLY,N,109.547
HIS,H,8.244
HIS,HA,4.599
HIS,HB2,3.104
HIS,HB3,3.048
HIS,HD1,8.907
HIS,HD2,7.001
HIS,HE1,7.951
HIS,HE2,9.628
HIS,CO,175.249
HIS,CA,56.469
HIS,CB,30.251
HIS,CD2,120.290
HIS,CE1,137.590
HIS,CG,132.265
HIS,N,119.731
HIS,ND1,193.190
HIS,NE2,183.264
ILE,H,8.256
ILE,HA,4.152
ILE,HB,1.788
ILE,HG12,1.275
ILE,HG13,1.204
ILE,HD,0.683
ILE,HG,0.778
ILE,CO,175.971
ILE,CA,61.676
ILE,CB,38.527
ILE,CD1,13.390
ILE,CG1,27.757
ILE,CG2,17.519
ILE,N,121.416
LEU,H,8.214
LEU,H1,7.340
LEU,HA,4.293
LEU,HB2,1.611
LEU,HB3,1.529
LEU,HG,1.512
LEU,HD1,0.755
LEU,HD2,0.736
LEU,CO,177.106
LEU,CA,55.674
LEU,CB,42.205
LEU,CD1,24.621
LEU,CD2,24.098
LEU,CG,26.782
LEU,N,121.845
LYS,H,8.177
LYS,HA,4.253
LYS,HB2,1.779
LYS,HB3,1.753
LYS,HD2,1.607
LYS,HD3,1.600
LYS,HE2,2.914
LYS,HE3,2.908
LYS,HG2,1.369
LYS,HG3,1.354
LYS,CO,176.725
LYS,CA,56.949
LYS,CB,32.737
LYS,CD,28.967
LYS,CE,41.896
LYS,CG,24.896
LYS,N,121.123
LYS,NZ,33.117
LYS,QZ,7.413
MET,H,8.251
MET,HA,4.388
MET,HB2,2.030
MET,HB3,1.995
MET,HG2,2.421
MET,HG3,2.396
MET,HE,1.891
MET,CO,176.256
MET,CA,56.129
MET,CB,32.904
MET,CE,17.121
MET,CG,32.029
MET,N,120.146
PHE,H,8.328
PHE,HA,4.604
PHE,HB2,3.002
PHE,HB3,2.942
PHE,HD1,7.058
PHE,HD2,7.062
PHE,HE1,7.088
PHE,HE2,7.085
PHE,HZ,6.998
PHE,CO,175.486
PHE,CA,58.118
PHE,CB,39.870
PHE,CD1,131.583
PHE,CD2,131.574
PHE,CE1,130.726
PHE,CE2,130.760
PHE,CG,138.248
PHE,CZ,129.252
PHE,N,120.363
PRO,H,8.519
PRO,HA,4.386
PRO,HB2,2.078
PRO,HB3,2.004
PRO,HD2,3.651
PRO,HD3,3.620
PRO,HG2,1.928
PRO,HG3,1.906
PRO,CO,176.770
PRO,CA,63.334
PRO,CB,31.834
PRO,CD,50.342
PRO,CG,27.200
PRO,N,135.632
SER,H,8.275
SER,HA,4.466
SER,HB2,3.870
SER,HB3,3.847
SER,HG,5.336
SER,CO,174.631
SER,CA,58.679
SER,CB,63.791
SER,N,116.317
THR,H,8.224
THR,HA,4.446
THR,HB,4.166
THR,HG1,5.061
THR,HG,1.138
THR,CO,174.554
THR,CA,62.203
THR,CB,69.694
THR,CG2,21.550
THR,N,115.329
TRP,H,8.260
TRP,HA,4.651
TRP,HB2,3.187
TRP,HB3,3.124
TRP,HD1,7.138
TRP,HE1,10.088
TRP,HE3,7.321
TRP,HH2,6.985
TRP,HZ2,7.289
TRP,HZ3,6.881
TRP,CO,176.244
TRP,CA,57.748
TRP,CB,29.885
TRP,CD1,126.575
TRP,CD2,127.266
TRP,CE2,137.632
TRP,CE3,120.500
TRP,CG,110.962
TRP,CH2,123.805
TRP,CZ2,114.261
TRP,CZ3,121.366
TRP,N,121.554
TRP,NE1,129.262
TYR,H,8.279
TYR,HA,4.599
TYR,HB2,2.909
TYR,HB3,2.846
TYR,HD1,6.941
TYR,HD2,6.937
TYR,HE1,6.702
TYR,HE2,6.704
TYR,HH,9.079
TYR,CO,175.541
TYR,CA,58.161
TYR,CB,39.227
TYR,CD1,132.730
TYR,CD2,132.720
TYR,CE1,117.947
TYR,CE2,117.922
TYR,CG,129.665
TYR,CZ,156.898
TYR,N,120.446
VAL,H,8.265
VAL,HA,4.157
VAL,HB,1.986
VAL,HG1,0.823
VAL,HG2,0.806
VAL,CO,175.717
VAL,CA,62.529
VAL,CB,32.668
VAL,CG1,21.490
VAL,CG2,21.303
VAL,N,121.079
"""
        for line in data.splitlines():
            line = line.split("\n")[0].split(",")
            if len(line) != 3:
                continue
            if line[0] not in self.bmrb_protein.keys():
                self.bmrb_protein[line[0]] = {}
            if "H" not in self.bmrb_protein[line[0]].keys():
                self.bmrb_protein[line[0]]["H"] = []
            if "C" not in self.bmrb_protein[line[0]].keys():
                self.bmrb_protein[line[0]]["C"] = []
            if "N" not in self.bmrb_protein[line[0]].keys():
                self.bmrb_protein[line[0]]["N"] = []
            if line[1][0] == "H":
                self.bmrb_protein[line[0]]["H"].append([line[1], float(line[2])])
            elif line[1][0] == "C":
                self.bmrb_protein[line[0]]["C"].append([line[1], float(line[2])])
            elif line[1][0] == "N":
                self.bmrb_protein[line[0]]["N"].append([line[1], float(line[2])])
        for key in self.bmrb_protein.keys():
            H_values = []
            H_labels = []
            N_values = []
            N_labels = []
            C_values = []
            C_labels = []
            for i in range(len(self.bmrb_protein[key]["H"])):
                H_values.append(self.bmrb_protein[key]["H"][i][1])
                H_labels.append(self.bmrb_protein[key]["H"][i][0])
            for i in range(len(self.bmrb_protein[key]["N"])):
                N_values.append(self.bmrb_protein[key]["N"][i][1])
                N_labels.append(self.bmrb_protein[key]["N"][i][0])
            for i in range(len(self.bmrb_protein[key]["C"])):
                C_values.append(self.bmrb_protein[key]["C"][i][1])
                C_labels.append(self.bmrb_protein[key]["C"][i][0])
            self.bmrb_protein[key]["H"] = [H_values, H_labels]
            self.bmrb_protein[key]["C"] = [C_values, C_labels]
            self.bmrb_protein[key]["N"] = [N_values, N_labels]

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