import wx # type: ignore
import os 
import nmrglue as ng # type: ignore
import numpy as np
import sys
from matplotlib.lines import Line2D
from SpinExplorer.SpinView.config import height, platform, colours, twoD_colours


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
        self.color_list = colours
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
        self.color_list = twoD_colours
        self.custom_lines = []
        self.custom_labels = []
        self.parent.active_plot_index = 0

    def OnDropFiles(self, x, y, filenames):

        if len(filenames)==1 and ".session" in filenames[0]:
            # Loading a new session
            # pos = self.GetPosition()
            # size = self.GetSize()
            self.parent.parent.Destroy()
            from SpinExplorer.SpinView.SpinView import SpinView
            frame = SpinView(session_file=filenames[0])
            
            
        else:
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
                            self.color_list = colours
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

                            self.parent.OnMinContour2D(wx.EVT_BUTTON,textcontrol=True)

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
