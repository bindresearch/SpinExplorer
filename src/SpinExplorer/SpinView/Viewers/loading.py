import wx # type: ignore
import numpy as np
import nmrglue as ng # type: ignore
import os

from matplotlib.lines import Line2D

from SpinExplorer.SpinView.Viewers.oned_view import OneDViewer
from SpinExplorer.SpinView.Viewers.twod_view import TwoDViewer

from SpinExplorer.SpinView.IO import ChooseFile, GetData

class ReadSession:
    def __init__(self, parent, session_file):
        self.main_frame = parent
        self.session_file = session_file
        # try:
        self.read_session_file()
        # except:
        #     # Give an error message saying that the session file could not be read (exiting)
        #     msg = wx.MessageDialog(self.main_frame, 'The session file could not be read. Please check the file paths in the .session file to ensure they are correct.', 'Error', wx.OK | wx.ICON_ERROR)
        #     msg.ShowModal()
        #     msg.Destroy()
        #     exit()

    def read_session_file(self):
        file = open(self.session_file, "r")
        lines = file.readlines()
        file.close()
        # Check to see what the first line is
        if len(lines[0].split("\n")[0].split()) == 1:
            # This is a session containing just one window
            # Check to see if the window is 1D, 2D or 3D
            plot_type = lines[0].split("\n")[0].split()[0]
            if plot_type == "1D" or plot_type == "1D stack":
                # This is a 1D window
                if plot_type == "1D":
                    stack = False
                else:
                    stack = True
                # Check to see if multiplot mode is on
                if lines[1].split("\n")[0].split(":")[1].split()[0] == "True":
                    multiplot = True
                    # Get the file path of the original data
                    file_path_original = lines[2].split("\n")[0].split("file_path:")[1]
                    self.main_frame.nmrdata = GetData(self, file_path_original)
                    self.main_frame.viewer = OneDViewer(
                        parent=self.main_frame, nmrdata=self.main_frame.nmrdata
                    )
                    self.main_frame.viewer.stack = stack
                    self.main_frame.main_sizer.Add(self.main_frame.viewer, 1, wx.EXPAND)
                    title = lines[3].split("\n")[0].split(":")[1]
                    p0_coarse = float(lines[4].split("\n")[0].split(":")[1])
                    p0_fine = float(lines[5].split("\n")[0].split(":")[1])
                    p1_coarse = float(lines[6].split("\n")[0].split(":")[1])
                    p1_fine = float(lines[7].split("\n")[0].split(":")[1])
                    colour = int(lines[8].split("\n")[0].split(":")[1])
                    linewidth = float(lines[9].split("\n")[0].split(":")[1])
                    reference_range = int(lines[10].split("\n")[0].split(":")[1])
                    reference_value = float(lines[11].split("\n")[0].split(":")[1])
                    vertical_range = int(lines[12].split("\n")[0].split(":")[1])
                    vertical_value = float(lines[13].split("\n")[0].split(":")[1])
                    multiply_range = int(lines[14].split("\n")[0].split(":")[1])
                    multiply_value = float(lines[15].split("\n")[0].split(":")[1])
                    pivot_point = float(lines[16].split("\n")[0].split(":")[1])
                    pivot_x = float(lines[17].split("\n")[0].split(":")[1])
                    pivot_visible = lines[18].split("\n")[0].split(":")[1]
                    self.choices = [title]
                    self.main_frame.viewer.line1.set_label(title)
                    self.main_frame.viewer.P0_slider.SetValue(p0_coarse)
                    self.main_frame.viewer.P1_slider.SetValue(p1_coarse)
                    self.main_frame.viewer.P0_slider_fine.SetValue(p0_fine)
                    self.main_frame.viewer.P1_slider_fine.SetValue(p1_fine)
                    self.main_frame.viewer.index = colour
                    self.main_frame.viewer.colour_chooser.SetSelection(colour)
                    self.main_frame.viewer.OnColourChoice1D(event=None)
                    self.main_frame.viewer.linewidth_slider.SetValue(linewidth)
                    self.main_frame.viewer.OnLinewidthScroll1D(event=None)
                    self.main_frame.viewer.reference_range_chooser.SetSelection(
                        reference_range
                    )
                    self.main_frame.viewer.OnReferenceCombo(event=None)
                    self.main_frame.viewer.reference_slider.SetValue(reference_value)
                    self.main_frame.viewer.OnReferenceScroll1D(event=None)
                    self.main_frame.viewer.vertical_range_chooser.SetSelection(
                        vertical_range
                    )
                    self.main_frame.viewer.OnVerticalCombo(event=None)
                    self.main_frame.viewer.vertical_slider.SetValue(vertical_value)
                    self.main_frame.viewer.OnVerticalScroll1D(event=None)
                    self.main_frame.viewer.multiply_range_chooser.SetSelection(
                        multiply_range
                    )
                    self.main_frame.viewer.OnMultiplyCombo(event=None)
                    self.main_frame.viewer.multiply_slider.SetValue(multiply_value)
                    self.main_frame.viewer.OnMultiplyScroll1D(event=None)
                    self.main_frame.viewer.pivot_line.set_xdata([pivot_point])
                    self.main_frame.viewer.pivot_x = pivot_x
                    if pivot_visible == "True":
                        self.main_frame.viewer.pivot_line.set_visible(True)
                    else:
                        self.main_frame.viewer.pivot_line.set_visible(False)
                    self.main_frame.viewer.OnSliderScroll1D(event=None)
                    # Read all the values into the values dictionary for the first plot
                    self.main_frame.viewer.values_dictionary[0] = {}
                    self.main_frame.viewer.values_dictionary[0][
                        "path"
                    ] = file_path_original
                    self.main_frame.viewer.values_dictionary[0]["title"] = title
                    self.main_frame.viewer.values_dictionary[0]["linewidth"] = linewidth
                    self.main_frame.viewer.values_dictionary[0]["color index"] = colour
                    self.main_frame.viewer.values_dictionary[0][
                        "original_ppms"
                    ] = self.main_frame.viewer.ppm_original
                    self.main_frame.viewer.values_dictionary[0][
                        "original_data"
                    ] = self.main_frame.viewer.nmrdata.data
                    self.main_frame.viewer.values_dictionary[0][
                        "dictionary"
                    ] = self.main_frame.viewer.nmrdata.dic
                    self.main_frame.viewer.values_dictionary[0][
                        "move left/right"
                    ] = reference_value
                    self.main_frame.viewer.values_dictionary[0][
                        "move left/right range index"
                    ] = reference_range
                    self.main_frame.viewer.values_dictionary[0][
                        "move up/down"
                    ] = vertical_value
                    self.main_frame.viewer.values_dictionary[0][
                        "move up/down range index"
                    ] = vertical_range
                    self.main_frame.viewer.values_dictionary[0][
                        "multiply value"
                    ] = multiply_value
                    self.main_frame.viewer.values_dictionary[0][
                        "multiply range index"
                    ] = multiply_range
                    self.main_frame.viewer.values_dictionary[0]["p0 Coarse"] = p0_coarse
                    self.main_frame.viewer.values_dictionary[0]["p1 Coarse"] = p1_coarse
                    self.main_frame.viewer.values_dictionary[0]["p0 Fine"] = p0_fine
                    self.main_frame.viewer.values_dictionary[0]["p1 Fine"] = p1_fine
                    self.main_frame.viewer.multiplot_mode = True
                    self.main_frame.viewer.files.first_drop = False
                    count = 1
                    # Loop over the rest of the lines to get the file paths of the other data

                    for i, line in enumerate(lines):
                        if i < 19:
                            continue
                        if line.split("\n")[0].split(":")[0] == "file_path":
                            file_path = line.split("\n")[0].split("file_path:")[1]
                            self.main_frame.viewer.values_dictionary[count] = {}
                            self.main_frame.viewer.values_dictionary[count][
                                "path"
                            ] = file_path
                        if line.split("\n")[0].split(":")[0] == "title":
                            title = line.split("\n")[0].split(":")[1]
                            self.main_frame.viewer.values_dictionary[count][
                                "title"
                            ] = title
                        if line.split("\n")[0].split(":")[0] == "linewidth":
                            linewidth = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "linewidth"
                            ] = linewidth
                        if line.split("\n")[0].split(":")[0] == "p0_coarse":
                            p0_coarse = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "p0 Coarse"
                            ] = p0_coarse
                        if line.split("\n")[0].split(":")[0] == "p1_coarse":
                            p1_coarse = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "p1 Coarse"
                            ] = p1_coarse
                        if line.split("\n")[0].split(":")[0] == "p0_fine":
                            p0_fine = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "p0 Fine"
                            ] = p0_fine
                        if line.split("\n")[0].split(":")[0] == "p1_fine":
                            p1_fine = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "p1 Fine"
                            ] = p1_fine
                        if line.split("\n")[0].split(":")[0] == "colour":
                            colour = int(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "color index"
                            ] = colour
                        if line.split("\n")[0].split(":")[0] == "reference_range":
                            reference_range = int(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "move left/right range index"
                            ] = reference_range
                        if line.split("\n")[0].split(":")[0] == "reference_value":
                            reference_value = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "move left/right"
                            ] = reference_value
                        if line.split("\n")[0].split(":")[0] == "vertical_range":
                            vertical_range = int(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "move up/down range index"
                            ] = vertical_range
                        if line.split("\n")[0].split(":")[0] == "vertical_value":
                            vertical_value = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "move up/down"
                            ] = vertical_value
                        if line.split("\n")[0].split(":")[0] == "multiply_range":
                            multiply_range = int(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "multiply range index"
                            ] = multiply_range
                        if line.split("\n")[0].split(":")[0] == "multiply_value":
                            multiply_value = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "multiply value"
                            ] = multiply_value
                        if line.split("\n")[0].split(":")[0] == "pivot_visible":
                            # Have put all the saved parameters into the dictionary
                            # Now need to add this plot to the canvas
                            self.add_saved_plot1D(count)
                            count += 1
                    self.main_frame.viewer.files.choices = self.choices      
                    
                    # Loop through all the plots and perform OnSelectPlot function
                    for i in range(len(self.main_frame.viewer.files.choices)):
                        self.main_frame.viewer.plot_combobox.SetSelection(i)
                        self.main_frame.viewer.OnSelectPlot(event=None)

                    self.main_frame.viewer.plot_combobox.SetSelection(0)
                    self.main_frame.viewer.OnSelectPlot(event=None)

          


                else:
                    multiplot = False
                    # Try to read in the data
                    for i, line in enumerate(lines):
                        if i < 2:
                            continue
                        else:
                            if line.split("\n")[0].split(":")[0] == "file_path":
                                file_path = line.split("\n")[0].split("file_path:")[0]
                            elif line.split("\n")[0].split(":")[0] == "p0_coarse":
                                p0_coarse = float(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "p0_fine":
                                p0_fine = float(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "p1_coarse":
                                p1_coarse = float(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "p1_fine":
                                p1_fine = float(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "colour":
                                colour = int(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "linewidth":
                                linewidth = float(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "reference_range":
                                reference_range = int(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "reference_value":
                                reference_value = float(
                                    line.split("\n")[0].split(":")[1]
                                )
                            elif line.split("\n")[0].split(":")[0] == "vertical_range":
                                vertical_range = int(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "vertical_value":
                                vertical_value = float(
                                    line.split("\n")[0].split(":")[1]
                                )
                            elif line.split("\n")[0].split(":")[0] == "multiply_range":
                                multiply_range = int(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "multiply_value":
                                multiply_value = float(
                                    line.split("\n")[0].split(":")[1]
                                )
                            elif line.split("\n")[0].split(":")[0] == "pivot_point":
                                pivot_point = float(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "pivot_x":
                                pivot_x = float(line.split("\n")[0].split(":")[1])
                            elif line.split("\n")[0].split(":")[0] == "pivot_visible":
                                pivot_visible = line.split("\n")[0].split(":")[1]
                    self.main_frame.nmrdata = GetData(self, file_path)
                    self.main_frame.viewer = OneDViewer(
                        parent=self.main_frame, nmrdata=self.main_frame.nmrdata
                    )
                    self.main_frame.main_sizer.Add(self.main_frame.viewer, 1, wx.EXPAND)
                    # Update the plot with all the saved parameters
                    self.main_frame.viewer.P0_slider.SetValue(p0_coarse)
                    self.main_frame.viewer.P1_slider.SetValue(p1_coarse)
                    self.main_frame.viewer.P0_slider_fine.SetValue(p0_fine)
                    self.main_frame.viewer.P1_slider_fine.SetValue(p1_fine)
                    self.main_frame.viewer.OnSliderScroll1D(event=None)
                    self.main_frame.viewer.colour_chooser.SetSelection(colour)
                    self.main_frame.viewer.OnColourChoice1D(event=None)
                    self.main_frame.viewer.linewidth_slider.SetValue(linewidth)
                    self.main_frame.viewer.OnLinewidthScroll1D(event=None)
                    self.main_frame.viewer.reference_range_chooser.SetSelection(
                        reference_range
                    )
                    self.main_frame.viewer.OnReferenceCombo(event=None)
                    self.main_frame.viewer.reference_slider.SetValue(reference_value)
                    self.main_frame.viewer.OnReferenceScroll1D(event=None)
                    self.main_frame.viewer.vertical_range_chooser.SetSelection(
                        vertical_range
                    )
                    self.main_frame.viewer.OnVerticalCombo(event=None)
                    self.main_frame.viewer.vertical_slider.SetValue(vertical_value)
                    self.main_frame.viewer.OnVerticalScroll1D(event=None)
                    self.main_frame.viewer.multiply_range_chooser.SetSelection(
                        multiply_range
                    )
                    self.main_frame.viewer.OnMultiplyCombo(event=None)
                    self.main_frame.viewer.multiply_slider.SetValue(multiply_value)
                    self.main_frame.viewer.OnMultiplyScroll1D(event=None)
                    self.main_frame.viewer.pivot_line.set_xdata([pivot_point])
                    self.main_frame.viewer.pivot_x = pivot_x
                    if pivot_visible == "True":
                        self.main_frame.viewer.pivot_line.set_visible(True)
                    else:
                        self.main_frame.viewer.pivot_line.set_visible(False)
                    self.main_frame.viewer.OnSliderScroll1D(event=None)
            elif lines[0].split("\n")[0].split()[0] == "2D":
                # This is a 2D window
                if lines[1].split("\n")[0].split(":")[1] == "True":
                    self.multiplot_mode = True
                    transposed2D = lines[2].split("\n")[0].split(":")[1]


                    # Get the file path of the original data
                    file_path_original = lines[3].split("\n")[0].split("file_path:")[1]
                    self.main_frame.nmrdata = GetData(self, file_path_original)
                    self.main_frame.viewer = TwoDViewer(
                        parent=self.main_frame, nmrdata=self.main_frame.nmrdata
                    )

                    self.main_frame.main_sizer.Add(self.main_frame.viewer, 1, wx.EXPAND)
                    title = lines[4].split("\n")[0].split(":")[1]
                    p0_coarse = float(lines[5].split("\n")[0].split(":")[1])
                    p0_fine = float(lines[6].split("\n")[0].split(":")[1])
                    p1_coarse = float(lines[7].split("\n")[0].split(":")[1])
                    p1_fine = float(lines[8].split("\n")[0].split(":")[1])
                    if(transposed2D=="False"):
                        move_x = float(lines[9].split("\n")[0].split(":")[1])
                        move_y = float(lines[10].split("\n")[0].split(":")[1])
                        move_x_index = int(lines[11].split("\n")[0].split(":")[1])
                        move_y_index = int(lines[12].split("\n")[0].split(":")[1])
                    else:
                        move_y = float(lines[9].split("\n")[0].split(":")[1])
                        move_x = float(lines[10].split("\n")[0].split(":")[1])
                        move_y_index = int(lines[11].split("\n")[0].split(":")[1])
                        move_x_index = int(lines[12].split("\n")[0].split(":")[1])

                    contour_linewidth = float(lines[13].split("\n")[0].split(":")[1])
                    multiply_factor = float(lines[14].split("\n")[0].split(":")[1])
                    contour_levels = int(lines[15].split("\n")[0].split(":")[1])
                    transposed = lines[16].split("\n")[0].split(":")[1]
                    self.choices = [title]
                    self.main_frame.viewer.contour_width_slider.SetValue(
                        contour_linewidth
                    )
                    self.main_frame.viewer.contour_levels_slider.SetValue(
                        contour_levels
                    )
                    self.main_frame.viewer.multiply_slider.SetValue(multiply_factor)
                    self.main_frame.viewer.reference_range_chooserX.SetSelection(
                        move_x_index
                    )
                    self.main_frame.viewer.reference_range_chooserY.SetSelection(
                        move_y_index
                    )
                    self.main_frame.viewer.move_x_slider.SetValue(move_x)
                    self.main_frame.viewer.move_y_slider.SetValue(move_y)
                    self.main_frame.viewer.P0_slider.SetValue(p0_coarse)
                    self.main_frame.viewer.P1_slider.SetValue(p1_coarse)
                    self.main_frame.viewer.P0_slider_fine.SetValue(p0_fine)
                    self.main_frame.viewer.P1_slider_fine.SetValue(p1_fine)

                    self.main_frame.viewer.files.first_drop = False
                    self.main_frame.viewer.OnSliderScroll2D(event=None)
                    # Read all the values into the values dictionary for the first plot
                    self.main_frame.viewer.values_dictionary[0] = {}
                    self.main_frame.viewer.values_dictionary[0][
                        "path"
                    ] = file_path_original
                    self.main_frame.viewer.values_dictionary[0]["title"] = title
                    self.main_frame.viewer.values_dictionary[0][
                        "contour linewidth"
                    ] = contour_linewidth
                    self.main_frame.viewer.values_dictionary[0][
                        "contour levels"
                    ] = contour_levels
                    self.main_frame.viewer.values_dictionary[0]["move x"] = move_x
                    self.main_frame.viewer.values_dictionary[0][
                        "move x range index"
                    ] = move_x_index
                    self.main_frame.viewer.values_dictionary[0]["move y"] = move_y
                    self.main_frame.viewer.values_dictionary[0][
                        "move y range index"
                    ] = move_y_index
                    self.main_frame.viewer.values_dictionary[0][
                        "multiply factor"
                    ] = multiply_factor
                    self.main_frame.viewer.values_dictionary[0]["p0 Coarse"] = p0_coarse
                    self.main_frame.viewer.values_dictionary[0]["p1 Coarse"] = p1_coarse
                    self.main_frame.viewer.values_dictionary[0]["p0 Fine"] = p0_fine
                    self.main_frame.viewer.values_dictionary[0]["p1 Fine"] = p1_fine
                    self.main_frame.viewer.values_dictionary[0][
                        "transposed"
                    ] = transposed
                    self.main_frame.viewer.values_dictionary[0][
                        "original_x_ppms"
                    ] = self.main_frame.viewer.ppms_0
                    self.main_frame.viewer.values_dictionary[0][
                        "original_y_ppms"
                    ] = self.main_frame.viewer.ppms_1
                    self.main_frame.viewer.values_dictionary[0]["new_x_ppms"] = (
                        self.main_frame.viewer.ppms_0
                        + np.ones(len(self.main_frame.viewer.ppms_0)) * move_x
                    )
                    self.main_frame.viewer.values_dictionary[0]["new_y_ppms"] = (
                        self.main_frame.viewer.ppms_1
                        + np.ones(len(self.main_frame.viewer.ppms_1)) * move_y
                    )
                    self.main_frame.viewer.values_dictionary[0][
                        "z_data"
                    ] = self.main_frame.viewer.nmrdata.data
                    self.main_frame.viewer.values_dictionary[0][
                        "uc0"
                    ] = self.main_frame.viewer.uc0
                    self.main_frame.viewer.values_dictionary[0][
                        "uc1"
                    ] = self.main_frame.viewer.uc1
                    self.main_frame.viewer.values_dictionary[0]["linewidth 1D"] = 1.0
                    self.main_frame.viewer.multiplot_mode = True
                    count = 1
                    self.main_frame.viewer.values_dictionary[count] = {}
                    # Loop over the rest of the lines to get the file paths of the other data
                    peaklists = False
                    peaklist_line = 0
                    for i, line in enumerate(lines):
                        if i < 17:
                            continue
                        elif('Peaklist path' in line):
                            peaklists = True
                            peaklist_line = i
                            break
                        else:
                            pass
                        if "file_path:" in line.split("\n")[0]:
                            file_path = line.split("\n")[0].split("file_path:")[1]
                            self.main_frame.viewer.values_dictionary[count][
                                "path"
                            ] = file_path
                        elif line.split("\n")[0].split(":")[0] == "title":
                            title = line.split("\n")[0].split(":")[1]
                            self.main_frame.viewer.values_dictionary[count][
                                "title"
                            ] = title
                        elif line.split("\n")[0].split(":")[0] == "contour linewidth":
                            contour_linewidth = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "contour linewidth"
                            ] = contour_linewidth
                        elif line.split("\n")[0].split(":")[0] == "move x":
                            if(transposed2D=='False'):
                                move_x = float(line.split("\n")[0].split(":")[1])
                                self.main_frame.viewer.values_dictionary[count][
                                    "move x"
                                ] = move_x
                            else:
                                move_y = float(line.split("\n")[0].split(":")[1])
                                self.main_frame.viewer.values_dictionary[count][
                                    "move y"
                                ] = move_y
                        elif line.split("\n")[0].split(":")[0] == "move x range index":
                            if(transposed2D=='False'):
                                move_x_index = int(line.split("\n")[0].split(":")[1])
                                self.main_frame.viewer.values_dictionary[count][
                                    "move x range index"
                                ] = move_x_index
                            else:
                                move_y_index = int(line.split("\n")[0].split(":")[1])
                                self.main_frame.viewer.values_dictionary[count][
                                    "move y range index"
                                ] = move_y_index
                        elif line.split("\n")[0].split(":")[0] == "move y":
                            if(transposed2D=='False'):
                                move_y = float(line.split("\n")[0].split(":")[1])
                                self.main_frame.viewer.values_dictionary[count][
                                    "move y"
                                ] = move_y
                            else:
                                move_x = float(line.split("\n")[0].split(":")[1])
                                self.main_frame.viewer.values_dictionary[count][
                                    "move x"
                                ] = move_x
                        elif line.split("\n")[0].split(":")[0] == "move y range index":
                            if(transposed2D=='False'):
                                move_y_index = int(line.split("\n")[0].split(":")[1])
                                self.main_frame.viewer.values_dictionary[count][
                                    "move y range index"
                                ] = move_y_index
                            else:
                                move_x_index = int(line.split("\n")[0].split(":")[1])
                                self.main_frame.viewer.values_dictionary[count][
                                    "move x range index"
                                ] = move_x_index
                        elif line.split("\n")[0].split(":")[0] == "multiply factor":
                            multiply_factor = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "multiply factor"
                            ] = multiply_factor
                        elif line.split("\n")[0].split(":")[0] == "p0 Coarse":
                            p0_coarse = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "p0 Coarse"
                            ] = p0_coarse
                        elif line.split("\n")[0].split(":")[0] == "p1 Coarse":
                            p1_coarse = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "p1 Coarse"
                            ] = p1_coarse
                        elif line.split("\n")[0].split(":")[0] == "p0 Fine":
                            p0_fine = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "p0 Fine"
                            ] = p0_fine
                        elif line.split("\n")[0].split(":")[0] == "p1 Fine":
                            p1_fine = float(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "p1 Fine"
                            ] = p1_fine
                        elif line.split("\n")[0].split(":")[0] == "contour levels":
                            contour_levels = int(line.split("\n")[0].split(":")[1])
                            self.main_frame.viewer.values_dictionary[count][
                                "contour levels"
                            ] = contour_levels
                        elif line.split("\n")[0].split(":")[0] == "transposed":
                            self.main_frame.viewer.values_dictionary[count][
                                "transposed"
                            ] = line.split("\n")[0].split(":")[1]
                            self.main_frame.viewer.values_dictionary[count][
                                "linewidth 1D"
                            ] = 1.0
                            self.add_saved_plot2D(count)
                            count += 1
                            self.main_frame.viewer.values_dictionary[count] = {}


                    # Search thrugh values dictionary and remove empty entries
                    keys = list(self.main_frame.viewer.values_dictionary.keys())
                    for key in keys:
                        if self.main_frame.viewer.values_dictionary[key] == {}:
                            del self.main_frame.viewer.values_dictionary[key]

                    self.plot_overlaid_2D()
                    if(peaklists==True):
                        self.main_frame.viewer.OnReadPeaks(wx.EVT_BUTTON)
                    for l in lines[peaklist_line+1:]:
                        peaklist_file=l.split('\n')[0]
                        self.main_frame.viewer.peaklist_frame.OnAddPeakList(wx.EVT_BUTTON, file=peaklist_file)
                    
                else:
                    self.multiplot_mode = False
                    transposed2D = lines[2].split("\n")[0].split(":")[1]
                    # Get the file path of the original data
                    file_path_original = lines[3].split("\n")[0].split("file_path:")[1]
                    self.main_frame.nmrdata = GetData(self, file_path_original)
                    self.main_frame.viewer = TwoDViewer(
                        parent=self.main_frame, nmrdata=self.main_frame.nmrdata
                    )
                    self.main_frame.main_sizer.Add(self.main_frame.viewer, 1, wx.EXPAND)
                    p0_coarse = float(lines[4].split("\n")[0].split(":")[1])
                    p0_fine = float(lines[5].split("\n")[0].split(":")[1])
                    p1_coarse = float(lines[6].split("\n")[0].split(":")[1])
                    p1_fine = float(lines[7].split("\n")[0].split(":")[1])
                    if(transposed2D==False):
                        move_x = float(lines[8].split("\n")[0].split(":")[1])
                        move_y = float(lines[9].split("\n")[0].split(":")[1])
                        move_x_index = int(lines[10].split("\n")[0].split(":")[1])
                        move_y_index = int(lines[11].split("\n")[0].split(":")[1])
                    else:
                        move_y = float(lines[8].split("\n")[0].split(":")[1])
                        move_x = float(lines[9].split("\n")[0].split(":")[1])
                        move_y_index = int(lines[10].split("\n")[0].split(":")[1])
                        move_x_index = int(lines[11].split("\n")[0].split(":")[1])
                    contour_linewidth = float(lines[12].split("\n")[0].split(":")[1])
                    multiply_factor = float(lines[13].split("\n")[0].split(":")[1])
                    contour_levels = int(lines[14].split("\n")[0].split(":")[1])
                    transposed = lines[15].split("\n")[0].split(":")[1]
                    self.main_frame.viewer.contour_width_slider.SetValue(
                        contour_linewidth
                    )
                    self.main_frame.viewer.contour_levels_slider.SetValue(
                        contour_levels
                    )
                    self.main_frame.viewer.multiply_slider.SetValue(multiply_factor)
                    self.main_frame.viewer.reference_range_chooserX.SetSelection(
                        move_x_index
                    )
                    self.main_frame.viewer.reference_range_chooserY.SetSelection(
                        move_y_index
                    )
                    self.main_frame.viewer.move_x_slider.SetValue(move_x)
                    self.main_frame.viewer.move_y_slider.SetValue(move_y)
                    self.main_frame.viewer.P0_slider.SetValue(p0_coarse)
                    self.main_frame.viewer.P1_slider.SetValue(p1_coarse)
                    self.main_frame.viewer.P0_slider_fine.SetValue(p0_fine)
                    self.main_frame.viewer.P1_slider_fine.SetValue(p1_fine)
                    self.main_frame.viewer.OnSliderScroll2D(event=None)
                    
                    peaklists = False
                    peaklist_line = 0
                    # Loop over the rest of the lines to get the file paths of the other data
                    for i, line in enumerate(lines):
                        if('Peaklist paths' in line):
                            peaklists = True
                            peaklist_line = i
                    if(peaklists==True):
                        self.main_frame.viewer.OnReadPeaks(wx.EVT_BUTTON)
                    for l in lines[peaklist_line+1:]:
                        peaklist_file=l.split('\n')[0]
                        self.main_frame.viewer.peaklist_frame.OnAddPeakList(wx.EVT_BUTTON, file=peaklist_file)


            elif lines[0].split("\n")[0].split()[0] == "3D":
                # This is a 3D window
                pass
            else:
                # Give a popout saying that the session file is not formatted correctly
                msg = wx.MessageDialog(
                    self,
                    "The session file is not formatted correctly",
                    "Error",
                    wx.OK | wx.ICON_ERROR,
                )
                msg.ShowModal()
                msg.Destroy()
                return

    def add_saved_plot1D(self, count):
        # Add the saved plot to the canvas
        # Read in the data
        dic, data = ng.pipe.read(
            self.main_frame.viewer.values_dictionary[count]["path"]
        )
        self.main_frame.viewer.values_dictionary[count]["original_data"] = data
        self.main_frame.viewer.values_dictionary[count]["dictionary"] = dic
        # Make the uc object
        uc0 = ng.pipe.make_uc(dic, data)

        # Get the ppm scale
        ppm_scale = uc0.ppm_scale()
        self.main_frame.viewer.values_dictionary[count]["original_ppms"] = ppm_scale

        data = data * self.main_frame.viewer.values_dictionary[count][
            "multiply value"
        ] + self.main_frame.viewer.values_dictionary[count]["move up/down"] * np.ones(
            len(data)
        )

        self.choices.append(self.main_frame.viewer.values_dictionary[count]["title"])
        self.main_frame.viewer.plot_combobox.Clear()
        self.main_frame.viewer.plot_combobox.AppendItems(self.choices)
        self.main_frame.viewer.plot_combobox.SetSelection(0)
        xlim, ylim = (
            self.main_frame.viewer.ax.get_xlim(),
            self.main_frame.viewer.ax.get_ylim(),
        )
        self.main_frame.viewer.extra_plots.append(
            self.main_frame.viewer.ax.plot(
                uc0.ppm_scale(),
                data,
                color=self.main_frame.viewer.colours[
                    self.main_frame.viewer.values_dictionary[count]["color index"]
                ],
                label=self.choices[-1],
                linewidth=self.main_frame.viewer.values_dictionary[count]["linewidth"],
            )
        )

        self.main_frame.viewer.ax.legend()
        self.main_frame.viewer.ax.set_xlim(xlim)
        self.main_frame.viewer.ax.set_ylim(ylim)
        self.main_frame.viewer.OnSelectPlot(wx.EVT_COMBOBOX)

    def add_saved_plot2D(self, count):
        # Add the saved plot to the canvas
        # Read in the data
        dic, data = ng.pipe.read(
            self.main_frame.viewer.values_dictionary[count]["path"]
        )
        self.main_frame.viewer.values_dictionary[count]["original_data"] = data
        self.main_frame.viewer.values_dictionary[count]["dictionary"] = dic
        # Make the uc object
        uc0 = ng.pipe.make_uc(dic, data, dim=0)
        uc1 = ng.pipe.make_uc(dic, data, dim=1)
        ppm0 = uc0.ppm_scale()
        ppm1 = uc1.ppm_scale()
        x, y = np.meshgrid(ppm1, ppm0)

        self.main_frame.viewer.values_dictionary[count]["original_x_ppms"] = ppm0
        self.main_frame.viewer.values_dictionary[count]["original_y_ppms"] = ppm1
        self.main_frame.viewer.values_dictionary[count]["new_x_ppms"] = ppm0
        self.main_frame.viewer.values_dictionary[count]["new_y_ppms"] = ppm1
        self.main_frame.viewer.values_dictionary[count]["z_data"] = data
        self.main_frame.viewer.values_dictionary[count]["contour linewidth"] = 1.0
        self.main_frame.viewer.values_dictionary[count]["linewidth 1D"] = 1.0
        self.main_frame.viewer.values_dictionary[count]["uc0"] = uc0
        self.main_frame.viewer.values_dictionary[count]["uc1"] = uc1

        length = len(self.main_frame.viewer.values_dictionary.keys())

        # If transpose is false, then the x-axis is the first axis and the y-axis is the second axis
        self.main_frame.viewer.values_dictionary[count]["new_x_ppms_old"] = (
            self.main_frame.viewer.values_dictionary[count]["new_x_ppms"]
        )
        self.main_frame.viewer.values_dictionary[count]["new_y_ppms_old"] = (
            self.main_frame.viewer.values_dictionary[count]["new_y_ppms"]
        )
        transposed = self.main_frame.viewer.values_dictionary[count]["transposed"]

        if transposed == "True":
            self.main_frame.viewer.values_dictionary[count]["new_x_ppms_old"] = (
                self.main_frame.viewer.values_dictionary[count]["new_x_ppms"]
            )
            self.main_frame.viewer.values_dictionary[count]["new_y_ppms_old"] = (
                self.main_frame.viewer.values_dictionary[count]["new_y_ppms"]
            )
            self.main_frame.viewer.values_dictionary[count]["new_x_ppms"] = (
                self.main_frame.viewer.values_dictionary[count]["new_y_ppms_old"]
            )
            self.main_frame.viewer.values_dictionary[count]["new_y_ppms"] = (
                self.main_frame.viewer.values_dictionary[count]["new_x_ppms_old"]
            )
            self.main_frame.viewer.values_dictionary[count]["original_x_ppms"] = (
                self.main_frame.viewer.values_dictionary[count]["new_x_ppms"]
            )
            self.main_frame.viewer.values_dictionary[count]["original_y_ppms"] = (
                self.main_frame.viewer.values_dictionary[count]["original_y_ppms"]
            )
            uc0 = self.main_frame.viewer.values_dictionary[count]["uc1"]
            uc1 = self.main_frame.viewer.values_dictionary[count]["uc0"]
            self.main_frame.viewer.values_dictionary[count]["uc0"] = uc0
            self.main_frame.viewer.values_dictionary[count]["uc1"] = uc1
            self.main_frame.viewer.values_dictionary[count]["z_data"] = (
                self.main_frame.viewer.values_dictionary[count]["z_data"].T
            )

    def plot_overlaid_2D(self):
        xlim, ylim = (
            self.main_frame.viewer.ax.get_xlim(),
            self.main_frame.viewer.ax.get_ylim(),
        )
        xlabel = self.main_frame.viewer.ax.get_xlabel()
        ylabel = self.main_frame.viewer.ax.get_ylabel()
        self.main_frame.viewer.ax.clear()
        self.main_frame.viewer.axes1D.clear()
        self.main_frame.viewer.axes1D_2.clear()
        self.main_frame.viewer.axes1D.set_yticks([])
        self.main_frame.viewer.axes1D_2.set_xticks([])
        self.main_frame.viewer.twoD_spectra = []
        self.main_frame.viewer.twoD_slices_horizontal = []
        self.main_frame.viewer.twoD_slices_vertical = []

        for i in range(len(self.main_frame.viewer.values_dictionary)):
            multiply_factor = self.main_frame.viewer.values_dictionary[i][
                "multiply factor"
            ]

            x, y = np.meshgrid(
                self.main_frame.viewer.values_dictionary[i]["new_y_ppms"],
                self.main_frame.viewer.values_dictionary[i]["new_x_ppms"],
            )
            self.main_frame.viewer.twoD_spectra.append(
                self.main_frame.viewer.ax.contour(
                    y,
                    x,
                    self.main_frame.viewer.values_dictionary[i]["z_data"]
                    * multiply_factor,
                    colors=self.main_frame.viewer.twoD_colours[i],
                    levels=self.main_frame.viewer.cl,
                    linewidths=self.main_frame.viewer.values_dictionary[i][
                        "contour linewidth"
                    ],
                )
            )

            if self.main_frame.viewer.transposed2D == False:
                self.main_frame.viewer.twoD_slices_horizontal.append(
                    self.main_frame.viewer.axes1D.plot(
                        self.main_frame.viewer.values_dictionary[i]["new_x_ppms"],
                        self.main_frame.viewer.values_dictionary[i]["z_data"][:, 1]
                        * multiply_factor,
                        color=self.main_frame.viewer.twoD_label_colours[i],
                        linewidth=self.main_frame.viewer.values_dictionary[i][
                            "linewidth 1D"
                        ],
                    )
                )
                self.main_frame.viewer.twoD_slices_vertical.append(
                    self.main_frame.viewer.axes1D_2.plot(
                        self.main_frame.viewer.values_dictionary[i]["new_y_ppms"],
                        self.main_frame.viewer.values_dictionary[i]["z_data"][1, :]
                        * multiply_factor,
                        color=self.main_frame.viewer.twoD_label_colours[i],
                        linewidth=self.main_frame.viewer.values_dictionary[i][
                            "linewidth 1D"
                        ],
                    )
                )

            else:
                if i == 0:
                    self.main_frame.viewer.twoD_slices_horizontal.append(
                        self.main_frame.viewer.axes1D.plot(
                            self.main_frame.viewer.values_dictionary[i]["new_x_ppms"],
                            self.main_frame.viewer.values_dictionary[i]["z_data"].T[
                                1, :
                            ]
                            * multiply_factor,
                            color=self.main_frame.viewer.twoD_label_colours[i],
                            linewidth=self.main_frame.viewer.values_dictionary[i][
                                "linewidth 1D"
                            ],
                        )
                    )
                    self.main_frame.viewer.twoD_slices_vertical.append(
                        self.main_frame.viewer.axes1D_2.plot(
                            self.main_frame.viewer.values_dictionary[i]["new_y_ppms"],
                            self.main_frame.viewer.values_dictionary[i]["z_data"].T[
                                :, 1
                            ]
                            * multiply_factor,
                            color=self.main_frame.viewer.twoD_label_colours[i],
                            linewidth=self.main_frame.viewer.values_dictionary[i][
                                "linewidth 1D"
                            ],
                        )
                    )
                else:
                    self.main_frame.viewer.twoD_slices_horizontal.append(
                        self.main_frame.viewer.axes1D.plot(
                            self.main_frame.viewer.values_dictionary[i]["new_x_ppms"],
                            self.main_frame.viewer.values_dictionary[i]["z_data"].T[
                                :, 1
                            ]
                            * multiply_factor,
                            color=self.main_frame.viewer.twoD_label_colours[i],
                            linewidth=self.main_frame.viewer.values_dictionary[i][
                                "linewidth 1D"
                            ],
                        )
                    )
                    self.main_frame.viewer.twoD_slices_vertical.append(
                        self.main_frame.viewer.axes1D_2.plot(
                            self.main_frame.viewer.values_dictionary[i]["new_y_ppms"],
                            self.main_frame.viewer.values_dictionary[i]["z_data"].T[
                                1, :
                            ]
                            * multiply_factor,
                            color=self.main_frame.viewer.twoD_label_colours[i],
                            linewidth=self.main_frame.viewer.values_dictionary[i][
                                "linewidth 1D"
                            ],
                        )
                    )
        self.main_frame.viewer.line_h = self.main_frame.viewer.ax.axhline(
            y=self.main_frame.viewer.values_dictionary[i]["new_x_ppms"][1],
            color="black",
            lw=1.5,
        )
        self.main_frame.viewer.line_v = self.main_frame.viewer.ax.axvline(
            x=self.main_frame.viewer.values_dictionary[i]["new_y_ppms"][1],
            color="black",
            lw=1.5,
        )
        self.main_frame.viewer.line_h.set_visible(False)
        self.main_frame.viewer.line_v.set_visible(False)
        self.custom_labels = []

        for i in range(len(self.main_frame.viewer.values_dictionary)):
            self.custom_labels.append(
                self.main_frame.viewer.values_dictionary[i]["title"]
            )
            self.main_frame.viewer.files.custom_lines.append(
                Line2D(
                    [0], [0], color=self.main_frame.viewer.twoD_label_colours[i], lw=1.5
                )
            )

        # Set all vertical and horizontal slices to invisible initia
        for i in range(len(self.main_frame.viewer.twoD_slices_horizontal)):
            self.main_frame.viewer.twoD_slices_horizontal[i][0].set_visible(False)
            self.main_frame.viewer.twoD_slices_vertical[i][0].set_visible(False)

        self.main_frame.viewer.ax.set_xlim(xlim)
        self.main_frame.viewer.ax.set_ylim(ylim)
        self.main_frame.viewer.ax.set_xlabel(xlabel)
        self.main_frame.viewer.ax.set_ylabel(ylabel)

        # Add labels of the extra plots to the select plot box
        self.main_frame.viewer.plot_combobox.Clear()
        self.main_frame.viewer.plot_combobox.AppendItems(self.custom_labels)
        self.main_frame.viewer.plot_combobox.SetSelection(self.main_frame.viewer.plot_combobox.GetCount()-1)
        self.main_frame.viewer.files.custom_labels = self.custom_labels
        self.main_frame.viewer.ax.legend(
            self.main_frame.viewer.files.custom_lines, self.custom_labels
        )

        self.main_frame.viewer.OnSelectPlot2D(wx.EVT_COMBOBOX)
        self.main_frame.viewer.UpdateFrame()