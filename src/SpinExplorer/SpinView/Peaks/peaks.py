import wx 
import os 
import pathlib
import numpy as np
import pandas as pd
import nmrglue as ng
import copy 
import re
import wx.grid as gridlib 
import matplotlib.patches as patches
from matplotlib.backend_bases import MouseEvent as MPLMouseEvent
import matplotlib
matplotlib.use("wxAgg")
from SpinExplorer.SpinView.Peaks import bore_candidates
from SpinExplorer.SpinView.Peaks.fit_peaks import fit_peaks
from SpinExplorer.SpinView.Peaks.fit_peaks import fit_peaks_2D_window
from SpinExplorer.SpinView.Peaks.analysis import analysis_frame

# Shown in front of the name of the peak mode which is selected
PEAK_MODE_MARK = "\u25cf "

# The linewidths of a peak, found by fitting, in each dimension. They are held
# both in Hz and in ppm
PEAK_LINEWIDTH_KEYS = [
    "linewidth1_hz",
    "linewidth1_ppm",
    "linewidth2_hz",
    "linewidth2_ppm",
]

# The two ways a reference plane can be used when picking peaks in a 3D
TRACE_METHOD = "Pick down the bore of each peak"
TEMPLATE_METHOD = "Pick in 3D and match"


# How many peaks each experiment gives down the bore dimension of a 3D, used to
# fill in the expected number when an experiment is chosen. Adding an experiment
# here is all that is needed to offer it
BORE_EXPERIMENTS = [
    ("None", None),
    ("HNCO", 1),
    ("HNCACO", 2),
    ("HNCANH", 2),
    ("HNCOCA", 1),
    ("HNCOCANH", 1),
]


# Everything which is held for each peak of a 2D peaklist, so that the lists
# stay the same length as peaks are added and removed
PEAK_ENTRY_KEYS = [
    "peak_name",
    "shift1",
    "shift2",
    "intensity",
] + PEAK_LINEWIDTH_KEYS


def find_axis_column_label(name, number: int) -> str:
    """
    How a chemical shift column is named in the peaklist table and in the header
    of a saved peaklist. The name of the axis is used when it is known, so that
    the columns say which dimension they hold.
    """
    name = str(name).strip()

    # The axis labels of a spectrum carry their units, which are added here
    for units in [" (ppm)", " (points)", "(ppm)", "(points)"]:
        if name.endswith(units) == True:
            name = name[: -len(units)].strip()

    if name == "" or name.lower() == "none":
        return "Shift {} (ppm)".format(number)

    return "{} (ppm)".format(name)


def find_axis_name(label) -> str:
    """
    The name of the axis a column label was made from, which is what a peaklist
    holds for each of its dimensions.
    """
    label = str(label).strip()

    for units in [" (ppm)", " (points)", "(ppm)", "(points)"]:
        if label.endswith(units) == True:
            label = label[: -len(units)].strip()

    if label.startswith("Shift") == True:
        return ""

    return label


def split_peaklist_header(line) -> list:
    """
    The column names of a peaklist header. The columns are separated by tabs, so
    that a column name made of more than one word stays in one piece.
    """
    if "\t" in line:
        return [name.strip() for name in line.split("\t")]

    return line.split()


def peaklist_file_is_empty(peaklist_file) -> bool:
    """
    Whether a peaklist file holds nothing at all, as a file which has only just
    been created does. Such a file is a peaklist which peaks have not been added
    to yet rather than a file which could not be read.
    """
    try:
        with open(peaklist_file) as file:
            for line in file:
                if line.strip() != "":
                    return False
    except OSError:
        return False

    return True


class PeakModeButtons:
    """
    Shared behaviour for the buttons which put the spectrum into a peak
    picking mode (adding, selecting and moving peaks).

    Only one of these modes can be used at a time, and none of them can be
    used while the pan or zoom tools are on, as those take the mouse clicks
    for themselves. Selecting a mode therefore turns the navigation tools
    off, and selecting pan or zoom turns the peak modes off.
    """

    # The buttons which turn on a peak mode, with the value which says
    # whether that mode is currently active. A window only has some of these.
    peak_mode_buttons = [
        ("add_peaks_button", "active_add", "Adding peaks"),
        ("select_peak_button", "active_select_peak", "Selecting a peak"),
        ("select_peaks_button", "active_select_peaks", "Selecting a peak group"),
        ("move_peaks_button", "active_move", "Moving peaks"),
        ("move_peaks_bore_button", "active_movez", "Moving peaks in z"),
    ]

    def find_peak_mode_buttons(self):
        """
        The peak mode buttons which this window has, as (button, value name,
        description).
        """
        found = []
        for button_name, flag_name, description in self.peak_mode_buttons:
            button = getattr(self, button_name, None)
            if button == None or hasattr(self, flag_name) == False:
                continue
            found.append((button, flag_name, description))

        return found

    def find_active_peak_mode(self):
        """
        The description of the peak mode which is currently on, or None when
        the spectrum is not in a peak mode.
        """
        for button, flag_name, description in self.find_peak_mode_buttons():
            if getattr(self, flag_name) == True:
                return description

        return None

    def find_peak_mode_label(self, button):
        """
        The name of a mode button without its mark. The button is also made
        wide enough to hold its name with the mark in front of it, so that
        marking it does not push the name onto a second line.
        """
        if button in self.peak_mode_labels:
            return self.peak_mode_labels[button]

        label = button.GetLabel()
        self.peak_mode_labels[button] = label

        button.SetLabel(PEAK_MODE_MARK + label)
        button.SetMinSize(button.GetBestSize())
        button.SetLabel(label)

        return label

    def update_mode_buttons(self):
        """
        Show which peak mode is selected. The button of the active mode is
        pressed in and marked, so that it is obvious what clicking on the
        spectrum will do.
        """
        changed = False
        for button, flag_name, description in self.find_peak_mode_buttons():
            active = getattr(self, flag_name) == True

            try:
                label = self.find_peak_mode_label(button)
                button.SetValue(active)
                button.SetLabel(PEAK_MODE_MARK + label if active else label)
                changed = True
            except RuntimeError:
                # The button has been destroyed
                continue

        if changed == True:
            try:
                self.Layout()
            except RuntimeError:
                pass

    def find_dataset_choices(self):
        """
        The spectra which a peaklist can belong to, named as they are in the
        spectrum window.
        """
        try:
            choices = list(self.main_frame.plot_combobox.GetItems())
        except (RuntimeError, AttributeError):
            choices = []

        if len(choices) == 0:
            choices = ["Main Plot"]

        return choices

    def find_peaklist_dataset(self, peaklist=None) -> int:
        """
        The spectrum which a peaklist belongs to, as an index into the list of
        spectra. A peaklist belongs to the spectrum which was selected when it
        was loaded until it is given a different one.
        """
        if peaklist == None:
            peaklist = self.current_peaklist_box.GetValue()

        index = self.peaklist_datasets.get(peaklist)

        if index == None or index >= len(self.find_dataset_choices()):
            return int(getattr(self.main_frame, "active_plot_index", 0))

        return int(index)

    def set_peaklist_dataset(self, peaklist, index):
        """
        Record which spectrum a peaklist belongs to.
        """
        if peaklist == "" or peaklist == None:
            return

        self.peaklist_datasets[peaklist] = int(index)

    def update_dataset_box(self):
        """
        Show which spectrum the selected peaklist belongs to.
        """
        dataset_box = getattr(self, "dataset_box", None)
        if dataset_box == None:
            return

        try:
            choices = self.find_dataset_choices()
            if list(dataset_box.GetItems()) != choices:
                dataset_box.SetItems(choices)

            peaklist = self.current_peaklist_box.GetValue()
            if peaklist == "":
                dataset_box.SetSelection(wx.NOT_FOUND)
                dataset_box.Enable(False)
                return

            dataset_box.Enable(True)
            dataset_box.SetSelection(self.find_peaklist_dataset(peaklist))
        except (RuntimeError, AttributeError):
            pass

    def update_peaklist_intensities(self, peaklist=None) -> bool:
        """
        Read the intensity of every peak of a peaklist from the spectrum the
        peaklist belongs to, returning whether the intensities were updated.
        This is needed when a peaklist is given a different spectrum, as the
        intensities it holds were read from the spectrum it belonged to before.
        """
        if peaklist == None:
            peaklist = self.current_peaklist_box.GetValue()

        if peaklist not in self.peak_list_dictionary:
            return False

        # Only the 2D window reads intensities from the spectrum
        find_new_intensity = getattr(self, "find_new_intensity", None)
        if find_new_intensity == None:
            return False

        dictionary = self.peak_list_dictionary[peaklist]
        intensities = []

        try:
            for i, peak_name in enumerate(dictionary["peak_name"]):
                intensities.append(
                    find_new_intensity(
                        dictionary["shift1"][i], dictionary["shift2"][i], peaklist
                    )
                )
        except (KeyError, IndexError, AttributeError, TypeError):
            # The spectrum cannot be read, the intensities are left as they were
            return False

        dictionary["intensity"] = intensities

        return True

    def find_remembered_peaklists(self) -> list:
        """
        The peaklists which were loaded the last time this peaks window was
        open, as [file, the spectrum it belongs to]. They are held by the
        spectrum window, which stays open while the peaks window is closed and
        opened again.
        """
        remembered = getattr(self.main_frame, "remembered_peaklists", None)
        if remembered == None:
            remembered = []
            self.main_frame.remembered_peaklists = remembered

        return remembered

    def find_remembered_path(self, peaklist_file) -> str:
        """
        A peaklist file as it is remembered, which is its whole path.
        """
        try:
            return str(pathlib.Path(peaklist_file).absolute())
        except TypeError:
            return str(peaklist_file)

    def remember_peaklist(self, peaklist_file, dataset=None):
        """
        Remember a peaklist so that it is loaded again when this window is
        opened next time, along with the spectrum it belongs to.
        """
        path = self.find_remembered_path(peaklist_file)
        remembered = self.find_remembered_peaklists()

        for held in remembered:
            if held[0] == path:
                held[1] = dataset
                return

        remembered.append([path, dataset])

    def forget_peaklist(self, peaklist_file):
        """
        Stop remembering a peaklist, so that it does not come back when this
        window is opened next time.
        """
        path = self.find_remembered_path(peaklist_file)
        remembered = self.find_remembered_peaklists()

        for held in list(remembered):
            if held[0] == path or str(held[0]).endswith(str(peaklist_file)) == True:
                remembered.remove(held)

    def load_remembered_peaklists(self):
        """
        Load the peaklists which were shown the last time this window was open,
        so that closing and opening it does not lose them.
        """
        for path, dataset in list(self.find_remembered_peaklists()):
            if os.path.exists(path) == False:
                # The file has gone, so it is not remembered any more
                self.forget_peaklist(path)
                continue

            try:
                self.AddPeaklist(path)
                if dataset != None:
                    self.set_peaklist_dataset(
                        self.current_peaklist_box.GetValue(), dataset
                    )
            except Exception:
                # A peaklist which cannot be read is not worth holding on to
                self.forget_peaklist(path)

        update_dataset_box = getattr(self, "update_dataset_box", None)
        if update_dataset_box != None:
            update_dataset_box()

    def zoom_to_region(self, axes, xlimits, ylimits):
        """
        Zoom a plot onto a region, keeping each axis running in the direction it
        already runs in. Chemical shift axes are drawn with the largest shift
        first, whereas an axis showing points runs the other way, so the
        direction is taken from the plot rather than assumed.
        """
        xmin, xmax = sorted([float(limit) for limit in xlimits])
        ymin, ymax = sorted([float(limit) for limit in ylimits])

        if axes.get_xlim()[0] > axes.get_xlim()[1]:
            axes.set_xlim([xmax, xmin])
        else:
            axes.set_xlim([xmin, xmax])

        if axes.get_ylim()[0] > axes.get_ylim()[1]:
            axes.set_ylim([ymax, ymin])
        else:
            axes.set_ylim([ymin, ymax])

    def OnRefreshIntensities(self, event):
        """
        Read the intensity of every peak of the current peaklist from the
        spectrum, leaving the positions of the peaks alone. This is used when a
        peaklist picked on one spectrum is loaded onto another one.
        """
        peaklist = self.find_current_peaklist()
        if peaklist == None:
            self.no_peaklist_message("Refreshing intensities")
            return

        # Keep the peaklist as it was, so that the intensities can be undone
        previous = getattr(self, "previous_peaklists", None)
        if previous != None:
            if len(previous) > 10:
                previous.pop(0)
            previous.append(copy.deepcopy(self.peak_list_dictionary))

        number_of_peaks = len(self.peak_list_dictionary[peaklist]["peak_name"])

        if self.update_peaklist_intensities(peaklist) == False:
            dlg = wx.MessageDialog(
                self,
                "The intensities could not be read from the spectrum, so the "
                "peaklist has been left as it was.",
                "Refreshing intensities",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return

        self.AddToTable()
        self.redraw_spectrum()

        dlg = wx.MessageDialog(
            self,
            "The intensities of the {} peaks of {} have been read from the "
            "spectrum. The positions of the peaks have not been changed.".format(
                number_of_peaks, peaklist
            ),
            "Refreshing intensities",
            wx.OK,
        )
        dlg.ShowModal()
        dlg.Destroy()

    def redraw_spectrum(self):
        """
        Draw the spectrum and its peaks again, whichever window this is.
        """
        for name, arguments in [
            ("OnMinContour2D", {"textcontrol": True}),
            ("OnBoreSlider", {}),
        ]:
            redraw = getattr(self.main_frame, name, None)
            if redraw != None:
                try:
                    redraw(wx.EVT_BUTTON, **arguments)
                except (RuntimeError, AttributeError, TypeError):
                    pass
                return

    def OnDatasetSelection(self, event):
        """
        The user has chosen which spectrum the selected peaklist belongs to.
        The intensities of its peaks are read again from the new spectrum.
        """
        peaklist = self.current_peaklist_box.GetValue()
        selection = self.dataset_box.GetSelection()

        if selection < 0 or selection == self.find_peaklist_dataset(peaklist):
            # The peaklist already belongs to this spectrum
            return

        self.set_peaklist_dataset(peaklist, selection)

        # The spectrum a peaklist belongs to is remembered along with it, so
        # that it is put back on the same spectrum next time
        find_peaklist_file = getattr(self, "find_peaklist_file", None)
        if find_peaklist_file != None:
            peaklist_file = find_peaklist_file()
            if peaklist_file != "":
                self.remember_peaklist(peaklist_file, selection)

        if self.update_peaklist_intensities(peaklist) == True:
            self.AddToTable()

    def find_current_peaklist(self):
        """
        The peaklist which is currently selected, or None if the selection
        does not name one of the peaklists which are loaded.
        """
        try:
            peaklist = self.current_peaklist_box.GetValue()
        except (RuntimeError, AttributeError):
            return None

        if peaklist in self.peak_list_dictionary:
            return peaklist

        return None

    def no_peaklist_message(self, title="Peak lists"):
        """
        Tell the user that the selected peaklist cannot be used.
        """
        dlg = wx.MessageDialog(
            None,
            "No peaklist is selected. Please choose one of the loaded peaklists and try again.",
            title,
            wx.OK | wx.ICON_WARNING,
        )
        dlg.ShowModal()
        dlg.Destroy()

    def find_navigation_toolbar(self):
        """
        The pan/zoom toolbar of the spectrum this peak window belongs to. The
        bore window names its toolbar differently from the 2D window.
        """
        main_frame = getattr(self, "main_frame", None)

        for name in ["toolbar", "toolbar_bore"]:
            toolbar = getattr(main_frame, name, None)
            if toolbar != None:
                return toolbar

        return None

    def find_navigation_mode(self) -> str:
        """
        Whether the spectrum is currently in pan or zoom mode.
        """
        toolbar = self.find_navigation_toolbar()
        if toolbar == None:
            return ""

        return str(getattr(toolbar, "mode", "")).lower()

    def turn_off_navigation(self):
        """
        Turn off the pan and zoom tools so that clicking on the spectrum
        picks peaks rather than moving the spectrum.
        """
        toolbar = self.find_navigation_toolbar()
        if toolbar == None:
            return

        mode = self.find_navigation_mode()

        # The tools are turned off by selecting them again, which is noticed
        # by the watcher below, so it is told to ignore this one
        self.changing_peak_mode = True
        try:
            if "pan" in mode:
                toolbar.pan()
            elif "zoom" in mode:
                toolbar.zoom()
        except Exception:
            pass
        finally:
            self.changing_peak_mode = False

    def on_navigation_selected(self):
        """
        Pan or zoom has been selected, so the peak modes are turned off.
        """
        if getattr(self, "changing_peak_mode", False) == True:
            return

        if self.find_navigation_mode() in ["", "none"]:
            # The tool has been turned off rather than selected
            return

        if self.find_active_peak_mode() == None:
            return

        self.turn_off_togglebuttons()
        self.update_mode_buttons()

    def watch_navigation_toolbar(self):
        """
        Follow the pan and zoom tools so that selecting one of them turns off
        any peak mode. The tools are followed by wrapping them, which covers
        both their buttons and their keyboard shortcuts.
        """
        toolbar = self.find_navigation_toolbar()
        if toolbar == None:
            return

        # The window the tools report to, so that a peak window which has been
        # closed and opened again is the one which is told
        toolbar.peak_modes_window = self

        if getattr(toolbar, "peak_modes_watched", False) == True:
            return

        def watched(tool):
            def use_tool(*args, **kwargs):
                result = tool(*args, **kwargs)
                window = getattr(toolbar, "peak_modes_window", None)
                if window != None:
                    try:
                        window.on_navigation_selected()
                    except RuntimeError:
                        # The peak window has been closed
                        toolbar.peak_modes_window = None
                return result

            return use_tool

        toolbar.pan = watched(toolbar.pan)
        toolbar.zoom = watched(toolbar.zoom)

        # The buttons of the toolbar are connected to the functions above when
        # the toolbar is made, so they are connected again here to make them
        # use the followed versions rather than the originals
        for name, tool in [("Pan", "pan"), ("Zoom", "zoom")]:
            tool_id = getattr(toolbar, "wx_ids", {}).get(name)
            if tool_id == None:
                continue
            toolbar.Bind(wx.EVT_TOOL, getattr(toolbar, tool), id=tool_id)

        toolbar.peak_modes_watched = True

    def follow_peak_modes(self):
        """
        Make the peak mode buttons turn the pan and zoom tools off, and keep
        the buttons showing which mode is selected. This is done by wrapping
        the functions the buttons and their shortcuts use, so that it applies
        however the mode was changed.
        """
        self.peak_mode_labels = {}
        self.changing_peak_mode = False

        handlers = [
            "OnAddPeaks",
            "OnSelectPeak",
            "OnSelectPeaks",
            "OnMovePeaks",
            "OnMovePeak",
            "OnMovePeakz",
            "OnRemovePeaks",
            "OnFindPeaks",
        ]

        for name in handlers:
            handler = getattr(self, name, None)
            if handler == None:
                continue
            setattr(self, name, self.with_peak_mode(handler))

    def with_peak_mode(self, handler):
        """
        Turn the navigation tools off whenever a peak mode is turned on, and
        show which mode is selected once the mode has changed.
        """

        def use_handler(event, *args, **kwargs):
            result = handler(event, *args, **kwargs)

            if self.find_active_peak_mode() != None:
                self.turn_off_navigation()

            self.update_mode_buttons()

            return result

        return use_handler


class PeakListWindow2D(PeakModeButtons, wx.Frame):
    def __init__(self, title, parent):
        """
        This class contains all the information relating to loading in
        peaklists. For now the peak lists loaded will not be saved in a
        session but hopefully that can be added in the future
        """
        self.main_frame = parent
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        width = 900
        height = 600
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.panel_peaklist = wx.Panel(self, -1)
        self.main_peaklist_sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_peaklist_sizer)


        self.set_initial_values()
        self.follow_peak_modes()
        self.make_peaklist_window()
        self.watch_navigation_toolbar()
        self.update_mode_buttons()
        self.update_dataset_box()

        # The spectrum draws the peaks of this window, which is not the one it
        # holds until this window has been made
        self.main_frame.peaklist_frame = self

        # The peaklists which were shown last time are loaded again
        self.load_remembered_peaklists()

        self.Show()

        self.Bind(wx.EVT_CLOSE, self.OnClose)

    def set_initial_values(self):
        """
        Setting initial values such as the peak list colour choices
        """
        self.peak_list_choices = [""]
        self.initial_peak_list_colours = ["black", "gray", "saddlebrown"]
        self.selected_colour = "darkviolet"
        self.peak_list_dictionary = {}
        self.selected_peakname = ""
        self.selected_peaklist = ""
        self.selected_peak_indexes = ""

        # Set when a peak has been picked out to be removed
        self.remove_peak = False

        # Flags showing whether a given button is active or not
        self.active_add = False
        self.active_select_peak = False
        self.active_select_peaks = False
        self.active_remove = False
        self.active_move = False
        self.active_find = False

        # The spectrum each peaklist belongs to, by peaklist name
        self.peaklist_datasets = {}

        self.rect = None
        self.start_point = None
        self.start_point_move = None

        self.old_key = None
        self.old_num = None

        # list to store the state of a peaklist at given time points
        self.previous_peaklists = []

        # list to store peaklists which are hidden
        self.hidden_peaklists = []

        # Initially, the selected area is set to an empty list which will get populated
        self.selected_area = []

        self.names = {}

        # A list to hold the paths of loaded peaklists
        self.peaklist_paths = []

    def make_peaklist_window(self):
        """
        This window will have the following:
        - a button to add peaklists
        - a selection of buttons associated with picking peaks using nmrglue
        - buttons to toggle add peak(s), select peak, select region, remove peak(s), move peak(s), find peak
        """

        self.row1_label = wx.StaticBox(self, -1, "Loading Peaklists:")
        self.row1 = wx.StaticBoxSizer(self.row1_label, wx.HORIZONTAL)

        self.add_peaklist_button = wx.Button(self.row1_label, label="Add peaklist")
        self.add_peaklist_button.Bind(wx.EVT_BUTTON, self.OnAddPeakList)

        self.peaklist_selection_text = wx.StaticText(self.row1_label, -1, "Selected Peaklist:")

        self.current_peaklist_box = wx.ComboBox(
            self.row1_label, choices=self.peak_list_choices, size=(250, 20),
            style=wx.CB_READONLY
        )
        
        self.current_peaklist_box.Bind(wx.EVT_COMBOBOX, self.OnPeakListSelection)

        self.remove_peaklist_button = wx.Button(self.row1_label, label="Remove Peaklist")
        self.remove_peaklist_button.Bind(wx.EVT_BUTTON, self.OnRemovePeakList)

        self.dataset_text = wx.StaticText(self.row1_label, -1, "Spectrum:")
        self.dataset_box = wx.ComboBox(
            self.row1_label, choices=self.find_dataset_choices(), size=(200, 20),
            style=wx.CB_READONLY
        )
        self.dataset_box.Bind(wx.EVT_COMBOBOX, self.OnDatasetSelection)



        self.row2_label = wx.StaticBox(
            self, -1, "Manipulate Peaklists: (shorcuts for Mac - Command+key in brackets)"
        )
        self.row2 = wx.StaticBoxSizer(self.row2_label, wx.VERTICAL)

        self.other_box_label = wx.StaticBox(
            self, -1, "Other options:"
        )
        self.other_sizer = wx.StaticBoxSizer(self.other_box_label, wx.HORIZONTAL)

        self.add_peaks_button = wx.ToggleButton(self.row2_label, label="Add Peaks (a)")
        self.add_peaks_button.Bind(wx.EVT_TOGGLEBUTTON, self.OnAddPeaks)
        ID_BUTTON_a = wx.NewIdRef()

        self.select_peak_button = wx.ToggleButton(self.row2_label, label="Select Peak (s)")
        self.select_peak_button.Bind(wx.EVT_TOGGLEBUTTON, self.OnSelectPeak)
        ID_BUTTON_s = wx.NewIdRef()

        self.select_peaks_button = wx.ToggleButton(self.row2_label, label="Select Peak Group (g)")
        self.select_peaks_button.Bind(wx.EVT_TOGGLEBUTTON, self.OnSelectPeaks)
        ID_BUTTON_g = wx.NewIdRef()

        self.remove_peaks_button = wx.Button(self.row2_label, label="Remove Peaks (r)")
        self.remove_peaks_button.Bind(wx.EVT_BUTTON, self.OnRemovePeaks)
        ID_BUTTON_r = wx.NewIdRef()

        self.find_peak_button = wx.Button(self.row2_label, label="Find Peak (f)")
        self.find_peak_button.Bind(wx.EVT_BUTTON, self.OnFindPeaks)
        ID_BUTTON_f = wx.NewIdRef()

        self.move_peaks_button = wx.ToggleButton(self.row2_label, label="Move Peaks (m)")
        self.move_peaks_button.Bind(wx.EVT_TOGGLEBUTTON, self.OnMovePeaks)
        ID_BUTTON_m = wx.NewIdRef()

        self.include_helper_box = wx.CheckBox(self.other_box_label, -1, 'Show helper dialogs')
        self.include_helper_box.SetValue(True)


        self.hide_peaklist = wx.CheckBox(self.other_box_label, -1, 'Hide Peaklist')
        self.hide_peaklist.SetValue(False)

        self.add_at_local_max_box = wx.CheckBox(
            self.other_box_label, -1, 'Add peaks at local maximum'
        )
        self.add_at_local_max_box.SetValue(False)
        self.hide_peaklist.Bind(wx.EVT_CHECKBOX, self.OnHidePeaklist)

        self.undo_button = wx.Button(self.other_box_label, label='Undo (u)')
        self.undo_button.Bind(wx.EVT_BUTTON, self.OnUndo)
        ID_BUTTON_u = wx.NewIdRef()


        self.move_to_local_max = wx.Button(self.row2_label, label='Move to local max (k)')
        self.move_to_local_max.Bind(wx.EVT_BUTTON, self.OnFindLocalMaximum)
        ID_BUTTON_k = wx.NewIdRef()


        self.fit_selected_peaks_button = wx.Button(self.row2_label, label='Fit Selected Peaks (p)')
        self.fit_selected_peaks_button.Bind(wx.EVT_BUTTON, self.OnFitSelectedPeaks)
        ID_BUTTON_p = wx.NewIdRef()

        

        # Creating an accelerator table for keyboard shortcuts for the buttons
        accelerator_table = wx.AcceleratorTable(
            [
                (wx.ACCEL_CTRL, ord("a"), ID_BUTTON_a),
                (wx.ACCEL_CTRL, ord("s"), ID_BUTTON_s),
                (wx.ACCEL_CTRL, ord("g"), ID_BUTTON_g),
                (wx.ACCEL_CTRL, ord("r"), ID_BUTTON_r),
                (wx.ACCEL_CTRL, ord("f"), ID_BUTTON_f),
                (wx.ACCEL_CTRL, ord("m"), ID_BUTTON_m),
                (wx.ACCEL_CTRL, ord("k"), ID_BUTTON_k),
                (wx.ACCEL_CTRL, ord("u"), ID_BUTTON_u),
                (wx.ACCEL_CTRL, ord("p"), ID_BUTTON_p),
            ]
        )

        self.SetAcceleratorTable(accelerator_table)
        self.main_frame.SetAcceleratorTable(accelerator_table)
        self.Bind(wx.EVT_MENU, self.OnAddPeaks, id=ID_BUTTON_a)
        self.Bind(wx.EVT_MENU, self.OnSelectPeak, id=ID_BUTTON_s)
        self.Bind(wx.EVT_MENU, self.OnSelectPeaks, id=ID_BUTTON_g)
        self.Bind(wx.EVT_MENU, self.OnRemovePeaks, id=ID_BUTTON_r)
        self.Bind(wx.EVT_MENU, self.OnMovePeaks, id=ID_BUTTON_m)
        self.Bind(wx.EVT_MENU, self.OnFindPeaks, id=ID_BUTTON_f)
        self.Bind(wx.EVT_MENU, self.OnFindLocalMaximum, id=ID_BUTTON_k)
        self.Bind(wx.EVT_MENU, self.OnUndo, id=ID_BUTTON_u)
        self.Bind(wx.EVT_MENU, self.OnFitSelectedPeaks, id=ID_BUTTON_p)

        self.main_frame.Bind(wx.EVT_MENU, self.OnAddPeaks, id=ID_BUTTON_a)
        self.main_frame.Bind(wx.EVT_MENU, self.OnSelectPeak, id=ID_BUTTON_s)
        self.main_frame.Bind(wx.EVT_MENU, self.OnSelectPeaks, id=ID_BUTTON_g)
        self.main_frame.Bind(wx.EVT_MENU, self.OnRemovePeaks, id=ID_BUTTON_r)
        self.main_frame.Bind(wx.EVT_MENU, self.OnMovePeaks, id=ID_BUTTON_m)
        self.main_frame.Bind(wx.EVT_MENU, self.OnFindPeaks, id=ID_BUTTON_f)
        self.main_frame.Bind(wx.EVT_MENU, self.OnFindLocalMaximum, id=ID_BUTTON_k)
        self.main_frame.Bind(wx.EVT_MENU, self.OnUndo, id=ID_BUTTON_u)
        self.main_frame.Bind(wx.EVT_MENU, self.OnFitSelectedPeaks, id=ID_BUTTON_p)

        self.save_peaks_button = wx.Button(self.other_box_label, label="Save")
        self.save_peaks_button.Bind(wx.EVT_BUTTON, self.OnSave)
        self.save_peaks_button.SetToolTip(
            "Save the peaklist into its own file, which is the one named in the "
            "peaklist box."
        )

        self.refresh_intensities_button = wx.Button(
            self.other_box_label, label="Refresh Intensities"
        )
        self.refresh_intensities_button.Bind(
            wx.EVT_BUTTON, self.OnRefreshIntensities
        )
        self.refresh_intensities_button.SetToolTip(
            "Read the intensity of every peak from the spectrum again, leaving "
            "the positions of the peaks alone. This is used when a peaklist "
            "picked on one spectrum is loaded onto another one."
        )

        self.save_peaks_as_button = wx.Button(self.other_box_label, label="Save As")
        self.save_peaks_as_button.Bind(wx.EVT_BUTTON, self.OnSaveAs)
        self.save_peaks_as_button.SetToolTip(
            "Save the peaklist into a different file, which it is then saved "
            "into from then on."
        )

        self.duplicate_peaklist_button = wx.Button(self.other_box_label, label="Duplicate Peaklist")
        self.duplicate_peaklist_button.Bind(wx.EVT_BUTTON, self.OnDuplicatePeaklist)


        self.row_pickpeaks_label = wx.StaticBox(self, -1, "Peak Picking (nmrglue) - performed on the current selected dataset:")
        self.row_pickpeaks = wx.StaticBoxSizer(self.row_pickpeaks_label, wx.VERTICAL)

        self.peak_picking_threshold_text = wx.StaticText(self.row_pickpeaks_label,-1,"Threshold (% of maximum):")
        self.peak_picking_threshold_box = wx.TextCtrl(self.row_pickpeaks_label,value='10.0',
                size=(50, 20))
        
        self.peak_picking_type_text = wx.StaticText(self.row_pickpeaks_label,-1,"Option:")
        types = ['Positive Peaks', 'Negative Peaks', 'Positive + Negative Peaks']
        self.peak_picking_type = wx.ComboBox(self.row_pickpeaks_label, choices = types, style=wx.CB_READONLY)
        
        self.peak_picking_algorithm_text = wx.StaticText(self.row_pickpeaks_label,-1,"Algorithm:")
        algorithms = ['thres', 'thres-fast', 'downward', 'connected']
        self.peak_picking_algorithm_box = wx.ComboBox(self.row_pickpeaks_label, choices=algorithms, style=wx.CB_READONLY)

        self.peaklist_name_text = wx.StaticText(self.row_pickpeaks_label,-1,"Peaklist name:")
        self.peaklist_name_box = wx.TextCtrl(self.row_pickpeaks_label,value='peaks_nmrglue.list',
                size=(200, 20))

        self.peak_pick_button = wx.Button(self.row_pickpeaks_label, label='Peak Pick')
        self.peak_pick_button.Bind(wx.EVT_BUTTON, self.OnPickPeaks)

        self.row_pickpeaks1 = wx.BoxSizer(wx.HORIZONTAL)
        self.row_pickpeaks2 = wx.BoxSizer(wx.HORIZONTAL)
        
        self.row_pickpeaks1.Add(self.peak_picking_threshold_text)
        self.row_pickpeaks1.AddSpacer(5)
        self.row_pickpeaks1.Add(self.peak_picking_threshold_box)
        self.row_pickpeaks1.AddSpacer(10)
        self.row_pickpeaks1.Add(self.peaklist_name_text)
        self.row_pickpeaks1.AddSpacer(10)
        self.row_pickpeaks1.Add(self.peaklist_name_box)
        self.row_pickpeaks1.AddSpacer(5)
        self.row_pickpeaks1.Add(self.peak_pick_button)

        self.row_pickpeaks2.Add(self.peak_picking_type_text)
        self.row_pickpeaks2.AddSpacer(5)
        self.row_pickpeaks2.Add(self.peak_picking_type)
        self.row_pickpeaks2.AddSpacer(10)
        self.row_pickpeaks2.Add(self.peak_picking_algorithm_text)
        self.row_pickpeaks2.AddSpacer(5)
        self.row_pickpeaks2.Add(self.peak_picking_algorithm_box)

        self.row_pickpeaks.Add(self.row_pickpeaks1, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.row_pickpeaks.AddSpacer(10)
        self.row_pickpeaks.Add(self.row_pickpeaks2, 0, wx.ALIGN_CENTER_HORIZONTAL)

        
        self.row1.AddSpacer(5)
        self.row1.Add(self.add_peaklist_button)
        self.row1.AddSpacer(10)
        self.row1.Add(self.peaklist_selection_text)
        self.row1.AddSpacer(5)
        self.row1.Add(self.current_peaklist_box)
        self.row1.AddSpacer(10)
        self.row1.Add(self.remove_peaklist_button)
        self.row1.AddSpacer(15)
        self.row1.Add(self.dataset_text, 0, wx.ALIGN_CENTER_VERTICAL)
        self.row1.AddSpacer(5)
        self.row1.Add(self.dataset_box)

        self.row2_1 = wx.BoxSizer(wx.HORIZONTAL)

        self.row2_1.AddSpacer(5)
        self.row2_1.Add(self.add_peaks_button)
        self.row2_1.AddSpacer(5)
        self.row2_1.Add(self.select_peak_button)
        self.row2_1.AddSpacer(5)
        self.row2_1.Add(self.select_peaks_button)
        self.row2_1.AddSpacer(5)
        self.row2_1.Add(self.move_peaks_button)
        self.row2_1.AddSpacer(5)
        self.row2_1.Add(self.remove_peaks_button)
        self.row2_1.AddSpacer(5)

        self.row2.Add(self.row2_1, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 5)
        self.row2.AddSpacer(5)

        self.row2_2 = wx.BoxSizer(wx.HORIZONTAL)
        
        self.row2_2.AddSpacer(10)
        self.row2_2.Add(self.move_to_local_max)
        self.row2_2.AddSpacer(10)
        self.row2_2.Add(self.fit_selected_peaks_button)
        self.row2_2.AddSpacer(10)
        self.row2_2.Add(self.find_peak_button)

        self.row2.Add(self.row2_2, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 5)



        self.other_sizer.Add(self.include_helper_box)
        self.other_sizer.AddSpacer(10)
        self.other_sizer.Add(self.hide_peaklist)
        self.other_sizer.AddSpacer(10)
        self.other_sizer.Add(self.undo_button)
        self.other_sizer.AddSpacer(10)
        self.other_sizer.Add(self.duplicate_peaklist_button)
        self.other_sizer.AddSpacer(10)
        self.other_sizer.Add(self.save_peaks_button)
        self.other_sizer.AddSpacer(5)
        self.other_sizer.Add(self.save_peaks_as_button)
        self.other_sizer.AddSpacer(10)
        self.other_sizer.Add(self.refresh_intensities_button)
        self.other_sizer.AddSpacer(10)
        self.other_sizer.Add(self.add_at_local_max_box)

        self.analysis_sizer_label = wx.StaticBox(
            self, -1, "Analysis options:"
        )
        self.analysis_sizer = wx.StaticBoxSizer(self.analysis_sizer_label, wx.HORIZONTAL)

        self.peaklist1_text = wx.StaticText(self.analysis_sizer_label, -1, label = 'Peaklist 1:')
        self.select_peaklist1 = wx.ComboBox(self.analysis_sizer_label, choices=self.peak_list_choices, size=(250, 20), style=wx.CB_READONLY)
        self.peaklist2_text = wx.StaticText(self.analysis_sizer_label, -1, label = 'Peaklist 2:')
        self.select_peaklist2 = wx.ComboBox(self.analysis_sizer_label, choices=self.peak_list_choices, size=(250, 20), style=wx.CB_READONLY)

        self.analyse_button = wx.Button(self.analysis_sizer_label, label="Plot CSPs + Intensities")
        self.analyse_button.Bind(wx.EVT_BUTTON, self.OnAnalyse)

        self.analysis_sizer.AddSpacer(5)
        self.analysis_sizer.Add(self.peaklist1_text)
        self.analysis_sizer.AddSpacer(5)
        self.analysis_sizer.Add(self.select_peaklist1)
        self.analysis_sizer.AddSpacer(5)
        self.analysis_sizer.Add(self.peaklist2_text)
        self.analysis_sizer.AddSpacer(5)
        self.analysis_sizer.Add(self.select_peaklist2)
        self.analysis_sizer.AddSpacer(5)
        self.analysis_sizer.Add(self.analyse_button)



        self.main_peaklist_sizer.AddSpacer(5)
        self.main_peaklist_sizer.Add(
            self.row1, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 5
        )
        self.main_peaklist_sizer.AddSpacer(5)
        self.main_peaklist_sizer.Add(
            self.row_pickpeaks, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 5
        )
        self.main_peaklist_sizer.AddSpacer(5)
        self.main_peaklist_sizer.Add(
            self.row2, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 5
        )
        self.main_peaklist_sizer.AddSpacer(5)
        self.main_peaklist_sizer.Add(
            self.other_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 5
        )
        self.main_peaklist_sizer.AddSpacer(5)
        self.main_peaklist_sizer.Add(
            self.analysis_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 5
        )

        # Then have a table of the currently loaded peaklist (originally blank)

        self.row3_label = wx.StaticBox(self, -1, "Peaklist Table:")
        self.row3 = wx.StaticBoxSizer(self.row3_label, wx.HORIZONTAL)

        self.grid = gridlib.Grid(self.row3_label)
        self.grid.CreateGrid(5, 8)

        self.grid.SetColLabelValue(0, "Peak name")
        self.grid.SetColLabelValue(1, "Shift 1 (ppm)")
        self.grid.SetColLabelValue(2, "Shift 2 (ppm)")
        self.grid.SetColLabelValue(3, "Intensity")
        self.grid.SetColLabelValue(4, "Linewidth 1 (Hz)")
        self.grid.SetColLabelValue(5, "Linewidth 1 (ppm)")
        self.grid.SetColLabelValue(6, "Linewidth 2 (Hz)")
        self.grid.SetColLabelValue(7, "Linewidth 2 (ppm)")

        # Bind event when cell value changes
        self.grid.Bind(gridlib.EVT_GRID_EDITOR_SHOWN, self.on_begin_edit)
        self.grid.Bind(gridlib.EVT_GRID_CELL_CHANGED, self.on_cell_changed)

        
        self.row3.Add(self.grid, proportion=1, flag=wx.EXPAND | wx.ALL, border=5)
        self.main_peaklist_sizer.AddSpacer(10)
        self.main_peaklist_sizer.Add(self.row3, 1, wx.EXPAND | wx.ALL, 5)

        self.Layout()
        self.Refresh()
        total_width = int(self.grid.GetClientSize().width * 0.8)
        col_count = self.grid.GetNumberCols()
        if col_count > 0:
            col_width = int(total_width // col_count)
            for c in range(col_count):
                self.grid.SetColSize(c, col_width)


    def OnClose(self, event):
        """
        Closing the peaklist window puts it away without losing anything: the
        peaklists stay loaded and their peaks stay on the spectrum, and opening
        the window again brings it back as it was. A peaklist is only taken off
        the spectrum by removing it.
        """
        # The peak modes take the mouse clicks on the spectrum, so they are
        # turned off while the window is away
        self.turn_off_togglebuttons()

        # Check if a Fit peaks result window is open
        continue_closing = self.check_fit_window()
        if(continue_closing==False):
            return

        # The peaks stay drawn on the spectrum
        self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)

        self.hide_window(event)

    def hide_window(self, event):
        """
        Put the window away rather than destroying it, so that everything it
        holds is still there when it is opened again.
        """
        self.Hide()

        try:
            if event.CanVeto() == True:
                event.Veto()
                return
        except AttributeError:
            pass

        # The window cannot stay, so it is destroyed as it used to be
        self.Destroy()


    
    def OnHidePeaklist(self, event):
        """
        Update the hidden_peaklist list and then spawn a redraw of
        the main application canvas
        """
        if(self.hide_peaklist.GetValue()==True):
            self.hidden_peaklists.append(self.current_peaklist_box.GetValue())
        else:
            self.hidden_peaklists.remove(self.current_peaklist_box.GetValue())
            if(self.hidden_peaklists==None):
                self.hidden_peaklists=[]

        self.turn_off_togglebuttons()
        self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)

        


    def OnUndo(self, event):
        """
        Go back to the previous peaklist dictionary and then update the main plot
        """
        if(len(self.previous_peaklists)==0):
            return
        self.peak_list_dictionary = self.previous_peaklists[-1]
        self.previous_peaklists.pop(-1)
        if(self.previous_peaklists==None):
            self.previous_peaklists=[]
        self.AddToTable()
        self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)


    def OnFitSelectedPeaks(self, event):
        """
        If there are peaks selected, fit these peaks using the data in
        the selected area dragged to select the peaks. If multiple peaks
        are selected, these peaks will be fitted together.

        The fit output is given as a popout with new peak positions and
        old peak positions shown. The popout has 2D contours and 3D surface
        plots of the overlaid data and fitted data points.
        
        The user can choose to accept these updates
        which will update the peak positions and intensities for the selected
        peaks in the peaklist window.
        """

        fit_selected = fit_peaks(self)
        if(fit_selected.fit_success!=None):
            # The fit was successful so opening up a window to show the results
            fit_result = fit_peaks_2D_window('Fit peaks result', self, fit_selected)


    
    def OnAnalyse(self, event):
        """
        1. Checking the peaklists
        2. Performing the analysis
        3. Opening a window to show the results
        """
        check = self.check_before_analysis()
        if(check==False):
            return
        
        self.build_analysis_dataframe()
        
        

    def build_analysis_dataframe(self):
        """
        Put the values for each peaklist into a common dataframe
        so that analysis of the results within this dataframe is
        simpler.
        """

        dataframe_layout = {'residue number':[], 'peak name':[], 'shift1 1 (ppm)':[], 'shift2 1 (ppm)':[], 'intensity 1':[], 'shift1 2 (ppm)':[], 'shift2 2 (ppm)':[], 'intensity 2':[], 'CSP (ppm)':[], 'Intensity Ratio':[]}

        df = pd.DataFrame(dataframe_layout)

        # Find out the atom types for each shift (e.g. H, N or C)
        atom_types = self.find_atom_types() # shift1 atom type, shift2 atom type

        for i, name in enumerate(self.peak_list_dictionary[self.selected_peaklist1]['peak_name']):
            new_row = {}
            residue_number = self.extract_number(name)
            try:
                index2 = self.peak_list_dictionary[self.selected_peaklist2]['peak_name'].index(name)
            except:
                # The given peak name is not present in both peaklists
                continue
            new_row['residue number'] = [residue_number]
            new_row['peak name'] = [name]

            new_row['shift1 1 (ppm)'] = [self.peak_list_dictionary[self.selected_peaklist1]['shift1'][i]]
            new_row['shift2 1 (ppm)'] = [self.peak_list_dictionary[self.selected_peaklist1]['shift2'][i]]
            new_row['intensity 1'] = [self.peak_list_dictionary[self.selected_peaklist1]['intensity'][i]]
            
            

            new_row['shift1 2 (ppm)'] = [self.peak_list_dictionary[self.selected_peaklist2]['shift1'][index2]]
            new_row['shift2 2 (ppm)'] = [self.peak_list_dictionary[self.selected_peaklist2]['shift2'][index2]]
            new_row['intensity 2'] = [self.peak_list_dictionary[self.selected_peaklist2]['intensity'][index2]]

            delta1 = self.peak_list_dictionary[self.selected_peaklist1]['shift1'][i] - self.peak_list_dictionary[self.selected_peaklist2]['shift1'][index2]
            delta2 = self.peak_list_dictionary[self.selected_peaklist1]['shift2'][i] - self.peak_list_dictionary[self.selected_peaklist2]['shift2'][index2]

            if(atom_types[0]=='H'):
                self.factor1 = 1
                self.factor2 = 0.14
                csp = np.sqrt(delta1**2 + 0.14*(delta2**2))
            else:
                self.factor1 = 0.14
                self.factor2 = 1
                csp = np.sqrt(0.14*(delta1**2) + (delta2**2))

            try:
                intensity_ratio = self.peak_list_dictionary[self.selected_peaklist1]['intensity'][i]/self.peak_list_dictionary[self.selected_peaklist2]['intensity'][index2]
            except:
                intensity_ratio = 0

            new_row['CSP (ppm)'] = [csp]
            new_row['Intensity Ratio'] = [intensity_ratio]

            df = pd.concat([df, pd.DataFrame(new_row)], ignore_index=True)

        df = df.sort_values('residue number')

        analysis = analysis_frame(self, df)

        

    def find_atom_types(self):
        label1 = 'H'
        label2 = 'H'
        axis_label1 = self.main_frame.ax.get_xlabel()
        if('H' in axis_label1):
            label1 = 'H'
        elif('N'in axis_label1):
            label1 = 'N'
        else:
            label1 = 'C'
        axis_label2 = self.main_frame.ax.get_ylabel()
        if('H' in axis_label2):
            label2 = 'H'
        elif('N'in axis_label2):
            label2 = 'N'
        else:
            label2 = 'C'

        return label1, label2


    def extract_number(self, s):
        match = re.match(r"(\d+)", s)
        return int(match.group(1)) if match else float("inf")


    def check_before_analysis(self):
        """
        1. Check that the two peaklists are different
        2. Check that the two peaklists have the same peak names
        """

        self.selected_peaklist1 = self.select_peaklist1.GetValue()
        self.selected_peaklist2 = self.select_peaklist2.GetValue()

        if(self.selected_peaklist1 == self.selected_peaklist2):
            dlg = wx.MessageDialog(
            None,
            "The two peaklists selected are the same. Please select different peaklists for Peaklist 1 and Peaklist 2 and try again.",
            "Peak Analysis",
            wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return False
        

        # Peaklist names
        names1 = sorted(self.peak_list_dictionary[self.selected_peaklist1]['peak_name'])
        names2 = sorted(self.peak_list_dictionary[self.selected_peaklist1]['peak_name'])

        if(names1 != names2):
            dlg = wx.MessageDialog(
            None,
            "The two peaklists do not have the same peak names. Please ensure peaks in the two selected peaklists have the same number of peaks with the same names and try again.",
            "Peak Analysis",
            wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return False
        

        # # Peaklist names
        # intensities1 = self.peak_list_dictionary[selected_peaklist1]['intensities']
        # intensities2 = self.peak_list_dictionary[selected_peaklist1]['intensities']

        # if(0 in intensities1 or 0 in intensities2):
        #     dlg = wx.MessageDialog(
        #     None,
        #     "At least one intensity value in one of the peaklists is currently set to 0. Please set the intensity by moving the peak(s) or by fitting the peak(s) with 0 intensity and try again.",
        #     "Peak Analysis",
        #     wx.OK,
        #     )
        #     dlg.ShowModal()
        #     dlg.Destroy()
        #     return False
        

        # Checking that the peaklist names do not have degenerate numbers

        numbers1 = []
        numbers2 = []

        degenerate_numbers = False
        for name in names1:
            number = self.extract_number(name)
            if(number in numbers1):
                degenerate_numbers = True
            numbers1.append(number)
        for name in names2:
            number = self.extract_number(name)
            if(number in numbers2):
                degenerate_numbers = True
            numbers2.append(number)


        if(degenerate_numbers == True):
            dlg = wx.MessageDialog(
            None,
            "At least one peaklist has more than one peak name containing the same number. Please ensure only one number (the residue number) is written for each peak name and try again.",
            "Peak Analysis",
            wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return False
        
        return True

    
    def OnAddPeakList(self, event, file=''):
        """
        1 - Open a file explorer window (opening at the current directory)
        2 - Try to read the peaklist file (might be necessary to transpose)
        3 - Plot the peaklist file (and when open mincontour 2D need to also
            plot the peaklists too)
        """

        if(file==''):
            # Opening up a file window asking the user to select the 1D peak list - must be in the format of 1st column = peak_name, 2nd column = peak_position
            dlg = wx.FileDialog(self, "Select the peak list", wildcard="", style=wx.FD_OPEN)
            dlg.SetDirectory(os.getcwd())
            if dlg.ShowModal() == wx.ID_OK:
                peaklist_file = dlg.GetPath()
            else:
                dlg.Destroy()
                return
        else:
            peaklist_file = file

        self.AddPeaklist(peaklist_file)

    def AddPeaklist(self, peaklist_file, new_peaklist=False):
        p = pathlib.Path(peaklist_file)
        dirs = p.parts[-3:]
        file_name = p.parts[-1]
        last_directories_path = str(pathlib.Path(*dirs))
        if ".xlsx" in file_name:
            # peaklist = self.ReadCCPNList(peaklist_file)
            pass
        else:
            peaklist = self.ReadPeakList(peaklist_file, new_peaklist, last_directories_path=last_directories_path)
        if type(peaklist) != dict:
            return
        self.peak_list_dictionary[last_directories_path] = peaklist
        if self.peak_list_choices == [""]:
            self.peak_list_choices = [last_directories_path]
        else:
            self.peak_list_choices.append(last_directories_path)

        self.current_peaklist_box.SetItems(self.peak_list_choices)
        self.current_peaklist_box.SetSelection(len(self.peak_list_choices) - 1)

        # The peaklist belongs to the spectrum which is currently selected
        self.set_peaklist_dataset(
            last_directories_path, getattr(self.main_frame, "active_plot_index", 0)
        )
        self.update_dataset_box()

        self.turn_off_togglebuttons()

        self.AddToTable()

        self.peaklist_paths.append(p)

        # so that the peaklist is loaded again when this window is opened next
        # time, along with the spectrum it belongs to
        self.remember_peaklist(
            peaklist_file, self.find_peaklist_dataset(last_directories_path)
        )

        self.update_comparison_boxes()

        self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)

    def update_comparison_boxes(self):
        """
        Show every loaded peaklist in the boxes used to compare two peaklists.
        A box keeps the peaklist it is showing when that peaklist is still
        loaded, so that loading another peaklist does not undo a comparison
        which has been set up.
        """
        for box in [self.select_peaklist1, self.select_peaklist2]:
            try:
                peaklist = box.GetValue()
                box.SetItems(self.peak_list_choices)
                if peaklist in self.peak_list_choices:
                    box.SetSelection(self.peak_list_choices.index(peaklist))
                else:
                    box.SetSelection(0)
            except (RuntimeError, AttributeError):
                pass

    def find_column_labels(self, peaklist=None, linewidths=True) -> list:
        """
        The names of the peaklist columns, used both for the table and for the
        header of a saved peaklist so that the two say the same thing.
        """
        if peaklist == None:
            peaklist = self.current_peaklist_box.GetValue()

        names = list(self.names.get(peaklist, ["", "", ""]))
        while len(names) < 3:
            names = names + [""]

        # A peaklist which did not come from a file has no names of its own, so
        # the axes of the spectrum are used
        spectrum = self.find_spectrum_axis_names()
        for dimension in [0, 1]:
            if str(names[dimension]).strip() == "":
                names[dimension] = spectrum[dimension]

        intensity = str(names[2]).strip()
        if intensity == "" or intensity.lower() == "intensity":
            intensity = "Intensity"

        labels = [
            "Peak name",
            find_axis_column_label(names[0], 1),
            find_axis_column_label(names[1], 2),
            intensity,
        ]

        if linewidths == True:
            labels += [
                "Linewidth 1 (Hz)",
                "Linewidth 1 (ppm)",
                "Linewidth 2 (Hz)",
                "Linewidth 2 (ppm)",
            ]

        return labels

    def find_spectrum_axis_names(self) -> list:
        """
        The names of the two axes of the spectrum as they are shown on the plot,
        used to name the columns of a peaklist which has no names of its own.
        """
        names = ["", ""]

        try:
            names[0] = self.main_frame.ax.get_xlabel()
            names[1] = self.main_frame.ax.get_ylabel()
        except (AttributeError, RuntimeError):
            pass

        return names

    def update_column_labels(self):
        """
        Show which dimension each shift column holds, using the names of the
        axes of the peaklist or of the spectrum.
        """
        try:
            for column, label in enumerate(self.find_column_labels()):
                if column < self.grid.GetNumberCols():
                    self.grid.SetColLabelValue(column, label)
        except (RuntimeError, AttributeError):
            pass

    def find_table_order(self, peaklist):
        """
        The peaks of a peaklist in the order the table shows them, as indexes
        into the lists of the peaklist. The table is sorted by the number in
        the peak name, so this is not always the order the peaks were added in.
        """
        if peaklist not in self.peak_list_dictionary:
            return []

        def extract_number(name):
            match = re.match(r"(\d+)", name)
            return int(match.group(1)) if match else float("inf")

        names = self.peak_list_dictionary[peaklist]["peak_name"]

        return [
            index
            for index, name in sorted(enumerate(names), key=lambda pair: extract_number(pair[1]))
        ]

    def find_linewidths_written(self, peaklist) -> bool:
        """
        Whether any peak of a peaklist has been fitted, so that the linewidth
        columns are worth writing into the peaklist file.
        """
        for key in PEAK_LINEWIDTH_KEYS:
            for index in range(len(self.peak_list_dictionary[peaklist]["peak_name"])):
                if self.find_linewidth(peaklist, key, index) != None:
                    return True

        return False

    def find_linewidth_columns(self, peaklist, index) -> str:
        """
        The fitted linewidths of one peak as they are written into a peaklist
        file. A peak which has not been fitted is written as nan.
        """
        values = []
        for key in PEAK_LINEWIDTH_KEYS:
            value = self.find_linewidth(peaklist, key, index)
            if value == None:
                values.append("nan")
            else:
                values.append("{:.5f}".format(value))

        return " \t ".join(values)

    def pad_peaklist(self, peaklist):
        """
        Make sure that every peak of a peaklist has an entry in every one of
        its lists. A peak which has not been fitted has no linewidths, and a
        peaklist read from a file has none at all.
        """
        if peaklist not in self.peak_list_dictionary:
            return

        dictionary = self.peak_list_dictionary[peaklist]
        number_of_peaks = len(dictionary.get("peak_name", []))

        for key in PEAK_ENTRY_KEYS:
            if key not in dictionary:
                dictionary[key] = []
            while len(dictionary[key]) < number_of_peaks:
                dictionary[key].append(None)

    def set_peak_linewidths(self, peaklist, index, linewidths1, linewidths2):
        """
        Record the fitted linewidths of one peak. Each pair is the linewidth in
        that dimension in Hz and in ppm, either of which can be None when it
        could not be found.
        """
        if peaklist not in self.peak_list_dictionary:
            return

        self.pad_peaklist(peaklist)

        dictionary = self.peak_list_dictionary[peaklist]
        for key, linewidth in zip(
            PEAK_LINEWIDTH_KEYS, list(linewidths1) + list(linewidths2)
        ):
            try:
                dictionary[key][index] = linewidth
            except (KeyError, IndexError):
                pass

    def find_linewidth(self, peaklist, key: str, index: int):
        """
        One fitted linewidth of one peak, or None when the peak has not been
        fitted.
        """
        try:
            value = self.peak_list_dictionary[peaklist][key][index]
        except (KeyError, IndexError, TypeError):
            return None

        if value == None:
            return None

        try:
            value = float(value)
        except (ValueError, TypeError):
            return None

        if value != value:
            # A linewidth which was saved as nan, so the peak was not fitted
            return None

        return value

    def find_linewidth_text(self, value, key: str) -> str:
        """
        How a linewidth is shown in the table. Peaks which have not been fitted
        are left blank.
        """
        if value == None:
            return ""

        if key.endswith("_ppm") == True:
            return "{:.4f}".format(value)

        return "{:.2f}".format(value)

    def find_hz_per_ppm(self, peaklist=None):
        """
        How many Hz there are in one ppm in each dimension of the spectrum a
        peaklist belongs to, used to give fitted linewidths in both units.
        None is returned for a dimension which cannot be read.
        """
        main_frame = self.main_frame

        try:
            if main_frame.multiplot_mode == False:
                unit_conversions = [main_frame.uc0, main_frame.uc1]
            else:
                values = main_frame.values_dictionary[
                    self.find_peaklist_dataset(peaklist)
                ]
                unit_conversions = [values["uc0"], values["uc1"]]
        except (KeyError, IndexError, AttributeError, TypeError):
            return [None, None]

        hz_per_ppm = []
        for unit_conversion in unit_conversions:
            try:
                ppm = unit_conversion.ppm(1) - unit_conversion.ppm(0)
                hz = unit_conversion.hz(1) - unit_conversion.hz(0)
                if ppm == 0:
                    hz_per_ppm.append(None)
                else:
                    hz_per_ppm.append(abs(hz / ppm))
            except (AttributeError, TypeError, ZeroDivisionError):
                hz_per_ppm.append(None)

        return hz_per_ppm

    def AddToTable(self):
        """
        Adding the peaklist just entered into the peaklist table
        """
        row_count = self.grid.GetNumberRows()
        if row_count > 0:
            self.grid.DeleteRows(0, row_count)
        peaklist = self.current_peaklist_box.GetValue()
        data = []

        if peaklist not in self.peak_list_dictionary:
            # No peaklist is loaded, the table is left empty
            return

        # The columns say which dimension they hold
        self.update_column_labels()

        # Every peak has an entry in every list, whether or not it has been fitted
        self.pad_peaklist(peaklist)

        for index in self.find_table_order(peaklist):
            peak = self.peak_list_dictionary[peaklist]["peak_name"][index]
            shift1 = self.peak_list_dictionary[peaklist]["shift1"][index]
            shift2 = self.peak_list_dictionary[peaklist]["shift2"][index]
            intensity = self.peak_list_dictionary[peaklist]["intensity"][index]
            row = [
                peak,
                "{:.5f}".format(shift1),
                "{:.5f}".format(shift2),
                "{:.5e}".format(intensity),
            ]
            for key in PEAK_LINEWIDTH_KEYS:
                row.append(
                    self.find_linewidth_text(
                        self.find_linewidth(peaklist, key, index), key
                    )
                )
            data.append(row)

        num_rows = self.grid.GetNumberRows()
        self.grid.AppendRows(len(data) - num_rows)
        for row, rowData in enumerate(data):
            for col, value in enumerate(rowData):
                self.grid.SetCellValue(row, col, str(value))
                if col >= 4:
                    # The linewidths come from fitting the peaks
                    self.grid.SetReadOnly(row, col, True)

    def on_begin_edit(self, event):
        """
        If the user is editing the peak_name column, store the original value
        """
        row = event.GetRow()
        col = event.GetCol()
        if col == 0:
            self.old_key = self.grid.GetCellValue(row, col)
        else:
            self.old_num = self.grid.GetCellValue(row, col)
        event.Skip()

    def on_cell_changed(self, event):
        """
        When a cell is changed, see if the types are correct
        e.g. the shifts are numbers.
        Can then update the dictionary and re-perform OnMinContour2D.
        """
        row = event.GetRow()
        col = event.GetCol()
        if self.old_key != None:
            peak_name = self.grid.GetCellValue(row, col)
            if (
                peak_name
                in self.peak_list_dictionary[self.current_peaklist_box.GetValue()][
                    "peak_name"
                ]
            ):
                # Give an error saying that this peak name is already taken, changing back to the original value
                self.grid.SetCellValue(row, col, self.old_key)
                dlg = wx.MessageDialog(
                    self,
                    "The peak name entered (row:{}, coloum:{})is already taken, this value has been reset to its previous value".format(
                        str(row), str(col)
                    ),
                    "Warning",
                    wx.OK,
                )
                dlg.ShowModal()
                dlg.Destroy()
            else:
                index = self.peak_list_dictionary[self.current_peaklist_box.GetValue()][
                    "peak_name"
                ].index(self.old_key)
                self.peak_list_dictionary[self.current_peaklist_box.GetValue()][
                    "peak_name"
                ][index] = peak_name

        else:
            peak_name = self.grid.GetCellValue(row, 0)
            index = self.peak_list_dictionary[self.current_peaklist_box.GetValue()][
                "peak_name"
            ].index(peak_name)
            try:
                new_value = float(self.grid.GetCellValue(row, col))
                if col == 1:
                    self.peak_list_dictionary[self.current_peaklist_box.GetValue()][
                        "shift1"
                    ][index] = new_value
                if col == 2:
                    self.peak_list_dictionary[self.current_peaklist_box.GetValue()][
                        "shift2"
                    ][index] = new_value
            except:
                dlg = wx.MessageDialog(
                    self,
                    "The value entered (row:{}, coloum:{})is not a number, this value has been reset to its previous value".format(
                        str(row), str(col)
                    ),
                    "Warning",
                    wx.OK,
                )
                dlg.ShowModal()
                dlg.Destroy()

        self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)

        self.old_key = None
        self.old_num = None

        self.AddToTable()

    def read_linewidths(self, dictionary, line):
        """
        Read the fitted linewidths of a peak from a line of a peaklist file.
        Peaklists written before the linewidths were fitted do not have these
        columns, and a peak which has not been fitted has nan in them.
        """
        for i, key in enumerate(PEAK_LINEWIDTH_KEYS):
            value = None
            try:
                value = float(line[4 + i])
                if value != value:
                    # nan, the peak was not fitted
                    value = None
            except (IndexError, ValueError, TypeError):
                value = None

            dictionary[key].append(value)

    def ReadPeakList(self, peaklist_file, new_peaklist, last_directories_path=''):
        """
        Read the selected peaklist to obtain the chemical shifts in each dimension
        Add a list of peak names, chemical shifts (dim1) and chemical shifts (dim2)
        to the dictionary
        """
        dictionary = {}
        dictionary["peak_name"] = []
        dictionary["shift1"] = []
        dictionary["shift2"] = []
        dictionary['intensity'] = []
        for key in PEAK_LINEWIDTH_KEYS:
            dictionary[key] = []
        
        # Placeholder names for each axis
        name1 = ''
        name2 = ''
        name3 = ''

        # Try to read the peaklist, otherwise give an error saying it could not be read correctly
        try:
            with open(peaklist_file) as file:
                lines = file.readlines()
                if len(lines) != 0:
                    for i, line in enumerate(lines):
                        line = line.split("\n")[0].split()
                        if len(line) >= 3:
                            if(i==0):
                                try:
                                    float(line[1])
                                    dictionary["peak_name"].append(line[0])
                                    dictionary["shift1"].append(float(line[1]))
                                    dictionary["shift2"].append(float(line[2]))
                                    try:
                                        dictionary["intensity"].append(float(line[3]))
                                    except:
                                        dictionary["intensity"].append(0)
                                    self.read_linewidths(dictionary, line)
                                except:
                                    # A header naming the columns, which is
                                    # separated by tabs so that names made of
                                    # more than one word stay in one piece
                                    header = split_peaklist_header(
                                        lines[i].split("\n")[0]
                                    )
                                    name1 = find_axis_name(header[1])
                                    name2 = find_axis_name(header[2])
                                    try:
                                        name3 = header[3]
                                    except IndexError:
                                        name3 = 'intensity'
                            else:
                                try:
                                    float(line[1])
                                    dictionary["peak_name"].append(line[0])
                                    dictionary["shift1"].append(float(line[1]))
                                    dictionary["shift2"].append(float(line[2]))
                                    try:
                                        dictionary["intensity"].append(float(line[3]))
                                    except:
                                        dictionary["intensity"].append(0)
                                    self.read_linewidths(dictionary, line)
                                except:
                                    pass

        except:
            self.peaklist_error_message()
            return None

        if len(dictionary["peak_name"]) == 0:
            if new_peaklist == False and peaklist_file_is_empty(peaklist_file) == False:
                # The file holds something which could not be read as peaks
                self.peaklist_error_message()
                return None

            # An empty peaklist, ready for peaks to be added to it
            self.names[last_directories_path] = [name1, name2, name3]
            return dictionary

        self.names[last_directories_path] = [name1, name2, name3]
    

        # Try to see if the chemical shifts of the peaks are within the 2D spectral range
        dictionary = self.check_peaklist(dictionary, last_directories_path)
        

        return dictionary

    # def ReadCCPNList(self, peaklist_file):
    #     """
    #     Read peaklist that has been exported from a CCPN peaklist table.
    #     """

    #     df = pd.read_excel(peaklist_file, dtype=str)

    #     peak_names = df.iloc[:, 0].tolist()
    #     shift1 = df.iloc[:, 7].to_numpy()
    #     shift2 = df.iloc[:, 8].to_numpy()

    #     shift1_1 = []
    #     shift2_1 = []

    #     for i in range(len(shift1)):
    #         shift1_1.append(float(shift1[i]))
    #         shift2_1.append(float(shift2[i]))

    #     dictionary = {}
    #     dictionary["peak_name"] = peak_names
    #     dictionary["shift1"] = shift1_1
    #     dictionary["shift2"] = shift2_1

    #     # Try to see if the chemical shifts of the peaks are within the 2D spectral range
    #     dictionary = self.check_peaklist(dictionary)

    #     return dictionary

    def check_peaklist(self, dictionary: dict, last_directories_path = ''):
        """
        Try to see if the chemical shifts of the peaks are within the 2D spectral range
        """
        ppms_0 = dictionary["shift1"]
        ppms_1 = dictionary["shift2"]

        axis_labels = self.main_frame.ax.get_xlabel(), self.main_frame.ax.get_ylabel()


        if(self.names[last_directories_path][0]!='' and self.names[last_directories_path][1]!=''):
            if(self.names[last_directories_path][0]==axis_labels[0] and self.names[last_directories_path][1]==axis_labels[1]):
                return dictionary
            elif(self.names[last_directories_path][0]==axis_labels[1] and self.names[last_directories_path][1]==axis_labels[0]):
                dictionary["shift1"] = ppms_1
                dictionary["shift2"] = ppms_0
                self.names[last_directories_path] = [copy.deepcopy(self.names[last_directories_path][1]), copy.deepcopy(self.names[last_directories_path][0])]
                return dictionary



        match_0 = []
        for ppm in ppms_0:
            if ppm > np.min(self.main_frame.ppms_0) and ppm < np.max(
                self.main_frame.ppms_0
            ):
                match_0.append(1)
            else:
                match_0.append(0)

        mean0 = np.mean(np.array(match_0))

        match_1 = []
        for ppm in ppms_1:
            if ppm > np.min(self.main_frame.ppms_1) and ppm < np.max(
                self.main_frame.ppms_1
            ):
                match_1.append(1)
            else:
                match_1.append(0)

        mean1 = np.mean(np.array(match_1))

        if mean0 == 0 and mean1 == 0:
            # No peaks are within the spectrum, trying transposing
            match_0 = []
            for ppm in ppms_0:
                if ppm > np.min(self.main_frame.ppms_1) and ppm < np.max(
                    self.main_frame.ppms_1
                ):
                    match_0.append(1)
                else:
                    match_0.append(0)

            mean0 = np.mean(np.array(match_0))

            match_1 = []
            for ppm in ppms_1:
                if ppm > np.min(self.main_frame.ppms_0) and ppm < np.max(
                    self.main_frame.ppms_0
                ):
                    match_1.append(1)
                else:
                    match_1.append(0)

            mean1 = np.mean(np.array(match_1))

            if mean0 > 0.5 and mean1 > 0.5:
                # More than 50 percent of the peaks are within the spectrum
                dictionary["shift1"] = ppms_1
                dictionary["shift2"] = ppms_0
                old_names = self.names[last_directories_path]
                self.names[last_directories_path] = [old_names[1], old_names[0]]
                if self.main_frame.transposed2D == True:
                    dictionary["shift1"] = ppms_0
                    dictionary["shift2"] = ppms_1
                    old_names = self.names[last_directories_path]
                    self.names[last_directories_path] = [old_names[1], old_names[0]]
                return dictionary

            else:
                return None

        else:
            if self.main_frame.transposed2D == True:
                dictionary["shift1"] = ppms_1
                dictionary["shift2"] = ppms_0
                old_names = self.names[last_directories_path]
                self.names[last_directories_path] = [old_names[1], old_names[0]]
            return dictionary

    def peaklist_error_message(self):
        """
        Gives the user an error when the peaklist was not read correctly
        """

        dlg = wx.MessageDialog(
            self,
            "The selected peaklist was not read correctly. Please select another peak list.",
            "Error",
        )
        dlg.ShowModal()

    def OnPeakListSelection(self, event):
        if self.selected_peaklist != "":
            self.selected_peaklist = self.current_peaklist_box.GetValue()

        self.turn_off_togglebuttons()

        self.AddToTable()

        # Update the hide peaklist selection to match the last stored value
        if(self.current_peaklist_box.GetValue() in self.hidden_peaklists):
            self.hide_peaklist.SetValue(True)
        else:
            self.hide_peaklist.SetValue(False)

        self.update_dataset_box()

        self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)

    def OnRemovePeakList(self, event):
        """
        Remove the selected peaklist from the peak window. The peaklist file
        itself is left alone, it is only removed from this session.
        """
        peaklist = self.find_current_peaklist()
        if peaklist == None:
            self.no_peaklist_message("Removing a peaklist")
            return

        dlg = wx.MessageDialog(
            self,
            "Remove the peaklist {} from this window? The peaklist file itself will not be deleted.".format(
                peaklist
            ),
            "Removing a peaklist",
            wx.YES_NO | wx.ICON_QUESTION,
        )
        result = dlg.ShowModal()
        dlg.Destroy()
        if result != wx.ID_YES:
            return

        self.turn_off_togglebuttons()

        index = self.peak_list_choices.index(peaklist)

        # Stop it coming back when this window is opened next time
        if index < len(self.peaklist_paths):
            self.forget_peaklist(self.peaklist_paths[index])
        self.forget_peaklist(peaklist)

        # Forget everything which was held for this peaklist
        self.peak_list_dictionary.pop(peaklist, None)
        self.peaklist_datasets.pop(peaklist, None)
        self.names.pop(peaklist, None)
        self.peak_list_choices.pop(index)
        if peaklist in self.hidden_peaklists:
            self.hidden_peaklists.remove(peaklist)
        if index < len(self.peaklist_paths):
            self.peaklist_paths.pop(index)

        self.selected_peaklist = ""
        self.selected_peakname = ""
        self.selected_peak_indexes = []

        if len(self.peak_list_choices) == 0:
            self.peak_list_choices = [""]

        for box in [self.current_peaklist_box, self.select_peaklist1, self.select_peaklist2]:
            box.SetItems(self.peak_list_choices)
            box.SetSelection(min(index, len(self.peak_list_choices) - 1))

        if self.current_peaklist_box.GetValue() != "":
            self.selected_peaklist = self.current_peaklist_box.GetValue()

        self.update_dataset_box()
        self.AddToTable()
        self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)

    def turn_off_togglebuttons(self):
        # If any toggle buttons are on, turn them off
        if self.active_add == True:
            self.active_add = False
            self.add_peaks_button.SetValue(False)
            self.main_frame.fig.canvas.mpl_disconnect(self.add_peak_connect)
        if self.active_move:
            self.active_move = False
            self.move_peaks_button.SetValue(False)
            if self.active_select_peak:
                self.main_frame.fig.canvas.mpl_disconnect(self.move_peak_connect)
            if self.active_select_peaks:
                self.main_frame.fig.canvas.mpl_disconnect(self.move_peak_press)
                self.main_frame.fig.canvas.mpl_disconnect(self.move_peak_motion)
                self.main_frame.fig.canvas.mpl_disconnect(self.move_peak_release)
        if self.active_select_peak == True:
            self.select_peak_button.SetValue(False)
            self.active_select_peak = False
            self.selected_peakname = ""
            self.main_frame.fig.canvas.mpl_disconnect(self.select_peak_connect)
            self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)
        if self.active_select_peaks:
            self.active_select_peaks = False
            self.rect = None
            self.start_point = None
            self.select_peaks_button.SetValue(False)
            self.main_frame.fig.canvas.mpl_disconnect(self.select_press)
            self.main_frame.fig.canvas.mpl_disconnect(self.select_release)
            self.main_frame.fig.canvas.mpl_disconnect(self.select_motion)
            self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)

        self.update_mode_buttons()

    def OnAddPeaks(self, event):
        """
        This will allow a user to add a peak to the currently selected peaklist
        A popout will come up saying that the user needs to use the cursor to
        add a peak. De-select the add button once complete.

        The code will also disable all the other buttons which have been
        selected
        """

        continue_function = self.check_hidden_peaklist()
        if(continue_function==False):
            self.add_peaks_button.SetValue(False)
            return

        if self.active_add == True:
            self.active_add = False
            self.add_peaks_button.SetValue(False)
            self.main_frame.fig.canvas.mpl_disconnect(self.add_peak_connect)
            return

        if self.active_select_peaks:
            self.active_select_peaks = False
            self.rect = None
            self.start_point = None
            self.select_peaks_button.SetValue(False)
            self.main_frame.fig.canvas.mpl_disconnect(self.select_press)
            self.main_frame.fig.canvas.mpl_disconnect(self.select_release)
            self.main_frame.fig.canvas.mpl_disconnect(self.select_motion)
            self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)

        if self.active_select_peak:
            self.select_peak_button.SetValue(False)
            self.active_select_peak = False
            self.selected_peakname = ""
            self.main_frame.fig.canvas.mpl_disconnect(self.select_peak_connect)
            self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)

        if self.peak_list_choices == [""]:
            dlg = wx.MessageDialog(
                None,
                "No peaklists are loaded, would you like to create a new peaklist?",
                "Adding peaks",
                wx.YES_NO,
            )
            result = dlg.ShowModal()
            if result == wx.ID_NO:
                dlg.Destroy()
                return

            dlg.Destroy()
            # Making a new peaklist, ask the user to create and save a new file in a file dialog
            dlg = wx.FileDialog(
                None,
                "Creating new peaklist",
                wildcard="Peaklists (*.list;*.txt)|*.list;*.txt|All files (*.*)|*.*",
                style=wx.FD_SAVE,
            )
            dlg.SetDirectory(os.getcwd())
            if dlg.ShowModal() == wx.ID_OK:
                peaklist_file = dlg.GetPath()
                try:
                    with open(peaklist_file, "w") as file:
                        pass
                except OSError as error:
                    dlg.Destroy()
                    dlg = wx.MessageDialog(
                        None,
                        "The new peaklist could not be created ({}). Please try "
                        "again in a folder which can be written to.".format(error),
                        "Adding peaks",
                        wx.OK,
                    )
                    dlg.ShowModal()
                    dlg.Destroy()
                    return

                self.AddPeaklist(peaklist_file, new_peaklist=True)

                if self.peak_list_choices == [""]:
                    # The peaklist was not loaded, so there is nothing to add
                    # peaks to and the add peaks mode is not turned on
                    dlg = wx.MessageDialog(
                        None,
                        "The new peaklist {} could not be loaded, so peaks "
                        "cannot be added to it.".format(peaklist_file),
                        "Adding peaks",
                        wx.OK,
                    )
                    dlg.ShowModal()
                    dlg.Destroy()
                    return

            else:
                dlg.Destroy()
                return

        # Updating the current active values
        self.active_add = True
        self.add_peaks_button.SetValue(True)
        self.add_peaks_button.SetForegroundColour(wx.Colour(60, 60, 60))


        # Connect the canvas click event to an add peak function
        self.add_peak_connect = self.main_frame.fig.canvas.mpl_connect(
            "button_press_event", self.on_click_addpeak
        )

        self.selected_peaklist = self.current_peaklist_box.GetValue()

        if(self.include_helper_box.GetValue()==True):
            dlg = wx.MessageDialog(
                None,
                "Peaks can now be added to the peaklist {} by clicking the cursor. Please de-select the add button when complete.".format(
                    self.current_peaklist_box.GetValue()
                ),
                "Adding Peaks",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()

    def on_click_addpeak(self, event):

        if(self.current_peaklist_box.GetValue() in self.hidden_peaklists):
            dlg = wx.MessageDialog(
                None,
                "Peaks cannot be added to a peaklist whilst it is hidden. Please untick the hide peaklist box and try again"
                "Adding Peaks",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return

        x, y = self.main_frame.ax.transData.inverted().transform((event.x, event.y))

        if x != None and y != None:

            if(len(self.previous_peaklists)>10):
                self.previous_peaklists.pop(0)
            self.previous_peaklists.append(copy.deepcopy(self.peak_list_dictionary))

            # Current peaklist
            current_peaklist = self.find_current_peaklist()
            if current_peaklist == None:
                self.no_peaklist_message("Adding Peaks")
                return

            part = ""
            number = 1
            order = [0, 1]

            if len(self.peak_list_dictionary[current_peaklist]["peak_name"]) > 0:

                peakname = self.peak_list_dictionary[current_peaklist]["peak_name"][-1]
                parts = re.findall(r"[A-Za-z_-]+|\d+", peakname)
                for i, v in enumerate(parts):
                    try:
                        v = int(v)
                        number = v + 1
                    except:
                        part = v
                        if i == 0:
                            order = [1, 0]
                if order == [0, 1]:
                    peakname = str(number) + part
                else:
                    peakname = part + str(number)

                if peakname in self.peak_list_dictionary[current_peaklist]["peak_name"]:
                    peakname = peakname + "_1"

            else:
                peakname = str(number) + part

            intensity = None
            if self.add_at_local_max_box.GetValue() == True:
                # Put the peak on the nearest maximum of the spectrum rather
                # than exactly where it was clicked
                x, y, intensity = self.find_local_maximum(x, y, current_peaklist)

            if intensity == None:
                intensity = self.find_new_intensity(x, y, current_peaklist)

            if self.check_duplicate_peak(x, y, current_peaklist) == False:
                # The user has chosen not to have two peaks on top of one
                # another, so the peaklist is left as it was
                self.previous_peaklists.pop()
                return

            self.peak_list_dictionary[current_peaklist]["peak_name"].append(peakname)
            self.peak_list_dictionary[current_peaklist]["shift1"].append(x)
            self.peak_list_dictionary[current_peaklist]["shift2"].append(y)
            self.peak_list_dictionary[current_peaklist]["intensity"].append(intensity)

            self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)
            self.AddToTable()


    def OnSelectPeak(self, event):
        """
        This will select a peak so that it can be moved etc
        """

        continue_function = self.check_hidden_peaklist()
        if(continue_function==False):
            self.select_peak_button.SetValue(False)
            return
        
        

        if self.active_move:
            if self.active_select_peak:
                self.select_peak_button.SetValue(True)
            return
        self.selected_peaklist = self.current_peaklist_box.GetValue()
        if self.selected_peaklist == "":
            dlg = wx.MessageDialog(
                None,
                "No peaklists are loaded, please load a peaklist and try again.",
                "Warning",
                wx.OK,
            )
            result = dlg.ShowModal()
            dlg.Destroy()
            return

        if self.active_select_peak == True:
            self.active_select_peak = False
            self.selected_peakname = ""
            self.select_peak_button.SetValue(False)
            self.main_frame.fig.canvas.mpl_disconnect(self.select_peak_connect)
            self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)
            return

        # First need to disable other toggle buttons that are selected
        if self.active_add == True:
            self.active_add = False
            self.add_peaks_button.SetValue(False)
            self.main_frame.fig.canvas.mpl_disconnect(self.add_peak_connect)
        if self.active_select_peaks == True:
            self.active_select_peaks = False
            self.rect = None
            self.start_point = None
            self.select_peaks_button.SetValue(False)
            self.main_frame.fig.canvas.mpl_disconnect(self.select_press)
            self.main_frame.fig.canvas.mpl_disconnect(self.select_release)
            self.main_frame.fig.canvas.mpl_disconnect(self.select_motion)
            self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)

        self.active_select_peak = True
        self.select_peak_button.SetValue(True)

        self.select_peak_connect = self.main_frame.fig.canvas.mpl_connect(
            "button_press_event", self.on_click_selectpeak
        )

    def on_click_selectpeak(self, event):
        """
        If the peak is within a tolerence select the peak
        If multiple peaks are within the tolerence, select the closest
        in terms of pixels on the screen.
        """

        if(self.current_peaklist_box.GetValue() in self.hidden_peaklists):
            dlg = wx.MessageDialog(
                None,
                "Peaks cannot be added to a peaklist whilst it is hidden. Please untick the hide peaklist box and try again"
                "Select Peaks",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return

        # Find the markers drawn for the currently selected peaklist
        points = self.find_peak_markers()
        if points == None:
            return

        cont, ind = points.contains(event)
        if cont:
            mouse_coordinates = [event.x, event.y]  # in pixels
            distances = []
            for index in ind["ind"]:
                x = self.peak_list_dictionary[self.current_peaklist_box.GetValue()][
                    "shift1"
                ][index]
                y = self.peak_list_dictionary[self.current_peaklist_box.GetValue()][
                    "shift2"
                ][index]
                x, y = self.main_frame.ax.transData.transform((x, y))
                distance = np.sqrt(
                    (mouse_coordinates[0] - x) ** 2 + (mouse_coordinates[1] - y) ** 2
                )
                distances.append(distance)

            min_index = ind["ind"][np.argmin(np.array(distances))]

            self.selected_peak_indexes = [min_index]
            self.selected_peakname = self.peak_list_dictionary[
                self.current_peaklist_box.GetValue()
            ]["peak_name"][min_index]
            self.remove_peak = True

        else:
            self.selected_peak_indexes = ["N/A"]

        self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)

    def OnSelectPeaks(self, event):
        """
        Giving a popout telling the user to drag a box over
        a region of the plot to select peaks in a given area
        """

        continue_function = self.check_hidden_peaklist()
        if(continue_function==False):
            self.select_peaks_button.SetValue(False)
            return

        if self.active_move:
            if self.active_select_peaks:
                self.select_peaks_button.SetValue(True)
            return

        self.selected_peaklist = self.current_peaklist_box.GetValue()

        if self.selected_peaklist == "":
            dlg = wx.MessageDialog(
                None,
                "No peaklists are loaded, please load a peaklist and try again.",
                "Warning",
                wx.OK,
            )
            result = dlg.ShowModal()
            dlg.Destroy()
            return

        if self.active_select_peaks == True:
            self.active_select_peaks = False
            self.rect = None
            self.start_point = None
            self.select_peaks_button.SetValue(False)
            self.main_frame.fig.canvas.mpl_disconnect(self.select_press)
            self.main_frame.fig.canvas.mpl_disconnect(self.select_release)
            self.main_frame.fig.canvas.mpl_disconnect(self.select_motion)
            self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)
            return

        # First de-select all activated toggles
        if self.active_add == True:
            self.active_add = False
            self.add_peaks_button.SetValue(False)
            self.main_frame.fig.canvas.mpl_disconnect(self.add_peak_connect)
        if self.active_select_peak == True:
            self.active_select_peak = False
            self.select_peak_button.SetValue(False)
            self.selected_peakname = ""
            self.main_frame.fig.canvas.mpl_disconnect(self.select_peak_connect)
            self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)
            return

        self.active_select_peaks = True
        self.select_peaks_button.SetValue(True)

        # self.selected_peaklist = self.current_peaklist_box.GetValue()
        if(self.include_helper_box.GetValue()==True):
            dlg = wx.MessageDialog(
                None,
                "Drag over multiple peaks to select a group. Multiple groups can be selected sequentially by repeating and holding down the shift key.",
                "Select Peaks",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()

        # If drag, finds new peaks
        self.select_press = self.main_frame.fig.canvas.mpl_connect(
            "button_press_event", self.on_press_select
        )
        self.select_release = self.main_frame.fig.canvas.mpl_connect(
            "button_release_event", self.on_release_select
        )
        self.select_motion = self.main_frame.fig.canvas.mpl_connect(
            "motion_notify_event", self.on_motion_select
        )

    def on_press_select(self, event):
        """
        This is activated when the mouse is clicked when select peaks
        is toggled
        """
        x, y = self.main_frame.ax.transData.inverted().transform((event.x, event.y))

        if x != None and y != None:
            self.start_point = (x, y)

            # Create the rectangle
            self.rect = patches.Rectangle(
                self.start_point, 0, 0, linewidth=1, edgecolor="red", facecolor="none"
            )
            self.main_frame.ax.add_patch(self.rect)
            self.main_frame.fig.canvas.draw()
            self.main_frame.UpdateFrame()

    def on_motion_select(self, event):
        """
        This is activated when the mouse is moved when select peaks
        is toggled after it has been clicked
        """
        if not self.start_point:
            return

        x, y = self.main_frame.ax.transData.inverted().transform((event.x, event.y))

        if x != None and y != None:

            # Update rectangle size
            x0, y0 = self.start_point
            x1, y1 = x, y
            width = x1 - x0
            height = y1 - y0

            self.rect.set_width(width)
            self.rect.set_height(height)
            self.rect.set_xy((x0, y0))
            self.main_frame.canvas.draw_idle()
            self.main_frame.UpdateFrame()

    def on_release_select(self, event):
        """
        This is activated when the mouse is released when select peaks
        is toggled after it has been clicked
        """
        if not self.start_point:
            return
        
        if(self.current_peaklist_box.GetValue() in self.hidden_peaklists):
            dlg = wx.MessageDialog(
                None,
                "Peaks cannot be added to a peaklist whilst it is hidden. Please untick the hide peaklist box and try again",
                "Select Peaks",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return
        
        x, y = self.main_frame.ax.transData.inverted().transform((event.x, event.y))

        if x != None and y != None:
            x0, y0 = self.start_point
            x1, y1 = x, y
            xmin, xmax = sorted([x0, x1])
            ymin, ymax = sorted([y0, y1])
            self.find_selected_peaks([xmin, xmax], [ymin, ymax], event)

        # Save the xmin,xmax and ymin,ymax positions so that these can be potentially used
        # to define an area to perform peak fitting in
        self.selected_area = [xmin, xmax, ymin, ymax]

        # Cleanup
        self.start_point = None
        self.rect.set_visible(False)
        self.rect = None
        self.main_frame.canvas.draw()
        self.main_frame.UpdateFrame()

    def find_selected_peaks(self, xcoords: list, ycoords: list, event):
        """
        Find any peaks in the current selected peaklist that are within
        the area just selected by the user.
        """

        if event.key and "shift" in event.key.lower():
            pass
        else:
            self.selected_peak_indexes = []

        for i, peak_name in enumerate(
            self.peak_list_dictionary[self.current_peaklist_box.GetValue()]["peak_name"]
        ):
            x = self.peak_list_dictionary[self.current_peaklist_box.GetValue()][
                "shift1"
            ][i]
            y = self.peak_list_dictionary[self.current_peaklist_box.GetValue()][
                "shift2"
            ][i]
            if x > xcoords[0] and x < xcoords[1]:
                if y > ycoords[0] and y < ycoords[1]:
                    self.selected_peak_indexes.append(i)

        if len(self.selected_peak_indexes) == 0:
            self.selected_peak_indexes = ["N/A"]
        else:
            # If have multiple peaks, add the ability to remove peaks
            self.remove_peak = True

        self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)

    def OnRemovePeaks(self, event):
        """
        If there is a current peak or peaks selected, then remove these peaks
        from the dictionary.

        If a peak or peaks are selected in the table of the Peak List window ask
        if the user if they want to remove these peaks.
        """

        continue_function = self.check_fit_window()
        if(continue_function==False):
            return

        if self.active_select_peak == True or self.active_select_peaks == True:
            if "N/A" not in self.selected_peak_indexes:
                if self.remove_peak == True:
                    count = 0
                    if(len(self.previous_peaklists)>10):
                        self.previous_peaklists.pop(0)
                    self.previous_peaklists.append(copy.deepcopy(self.peak_list_dictionary))
                    peaklist = self.current_peaklist_box.GetValue()
                    self.pad_peaklist(peaklist)
                    for peak_index in self.selected_peak_indexes:
                        # Everything held for the peak is removed, so that the
                        # lists stay in step with one another
                        for key in PEAK_ENTRY_KEYS:
                            try:
                                del self.peak_list_dictionary[peaklist][key][
                                    peak_index - count
                                ]
                            except (KeyError, IndexError):
                                pass

                        count += 1

                    self.remove_peak = False
                    self.selected_peak_indexes = ["N/A"]
                    self.selected_peakname = ""

                    self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)
                    self.AddToTable()

    def OnMovePeaks(self, event):
        """
        This function is activated when the user clicks on the move peaks button.
        The function first deactivates the select peak matplotlib connect functions.
        The code then checks to see if there are peaks selected.
        If peaks are selected then the user is able to click a new peak position (if
        one peak is selected) or drag the selected peaks to new positions (if multiple
        peaks are selected).
        """

        continue_function = self.check_hidden_peaklist()
        if(continue_function==False):
            self.move_peaks_button.SetValue(False)
            return
    

        if self.active_move == True:
            self.active_move = False
            self.move_peaks_button.SetValue(False)

            if self.active_select_peaks == True:
                self.main_frame.fig.canvas.mpl_disconnect(self.move_peak_press)
                self.main_frame.fig.canvas.mpl_disconnect(self.move_peak_motion)
                self.main_frame.fig.canvas.mpl_disconnect(self.move_peak_release)
                self.select_press = self.main_frame.fig.canvas.mpl_connect(
                    "button_press_event", self.on_press_select
                )
                self.select_release = self.main_frame.fig.canvas.mpl_connect(
                    "button_release_event", self.on_release_select
                )
                self.select_motion = self.main_frame.fig.canvas.mpl_connect(
                    "motion_notify_event", self.on_motion_select
                )
            if self.active_select_peak == True:
                self.main_frame.fig.canvas.mpl_disconnect(self.move_peak_connect)
                self.select_peak_connect = self.main_frame.fig.canvas.mpl_connect(
                    "button_press_event", self.on_click_selectpeak
                )

            return

        # Temporarily deactivate the ability to select peak or select group
        if len(self.selected_peak_indexes) == 0 or "N/A" in self.selected_peak_indexes:
            # return as there are no selected peaks
            self.active_move = False
            self.move_peaks_button.SetValue(False)
            dlg = wx.MessageDialog(
                self,
                "There are no peaks selected. Please select a peak or a group of peaks and try again.",
                "Warning",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return

        self.eventDict = {}
        for name in dir(wx):
            if name.startswith("EVT_"):
                evt = getattr(wx, name)
                if isinstance(evt, wx.PyEventBinder):
                    self.eventDict[evt.typeId] = name

        if self.active_select_peaks == True:
            evt_id = event.GetEventType()
            self.main_frame.fig.canvas.mpl_disconnect(self.select_press)
            self.main_frame.fig.canvas.mpl_disconnect(self.select_release)
            self.main_frame.fig.canvas.mpl_disconnect(self.select_motion)
            if self.eventDict[evt_id] != wx.EVT_TOOL_RANGE.typeId:
                if(self.include_helper_box.GetValue()==True):
                    dlg = wx.MessageDialog(
                        self,
                        "Please drag to move the selected peaks to a new location. This can be repeated. Un-toggle the move peaks button when completed. (Note: ensure that zoom/pan in the matplotlib toolbar is not selected. Zoom/pan before entering move peaks mode)",
                        "Move Peaks",
                        wx.OK,
                    )
                    dlg.ShowModal()
                    dlg.Destroy()
        if self.active_select_peak == True:
            self.main_frame.fig.canvas.mpl_disconnect(self.select_peak_connect)
            evt_id = event.GetEventType()
            if evt_id != wx.EVT_TOOL_RANGE.typeId:
                if(self.include_helper_box.GetValue()==True):
                    dlg = wx.MessageDialog(
                        self,
                        "Please click a new location to move the selected peak. This can be repeated. Un-toggle the move peaks button when completed. (Note: ensure zoom/pan in the matplotlib toolbar is not selected. Zoom/pan before entering move peaks mode)",
                        "Move Peaks",
                        wx.OK,
                    )
                    dlg.ShowModal()
                    dlg.Destroy()

        self.active_move = True
        self.move_peaks_button.SetValue(True)

        # If select peak, and there is a peak selected give a popout telling
        # the user to click where they want the peak to go
        if self.active_select_peak == True:
            self.move_peak_connect = self.main_frame.fig.canvas.mpl_connect(
                "button_press_event", self.on_click_movepeak
            )

        if self.active_select_peaks == True:
            self.move_peak_press = self.main_frame.fig.canvas.mpl_connect(
                "button_press_event", self.on_press_movepeak
            )
            self.move_peak_motion = self.main_frame.fig.canvas.mpl_connect(
                "motion_notify_event", self.on_motion_movepeak
            )
            self.move_peak_release = self.main_frame.fig.canvas.mpl_connect(
                "button_release_event", self.on_release_movepeak
            )

        # when dragging, update the x/y coordinates of self.main_frame.points and then
        # redraw the canvas.




    def on_click_movepeak(self, event):
        """
        This function will update the peak position of the selected peak
        depending on where the user clicked.
        """
        x, y = self.main_frame.ax.transData.inverted().transform((event.x, event.y))

        if x != None and y != None:

            self.get_current_data()

            if(len(self.previous_peaklists)>10):
                self.previous_peaklists.pop(0)
            self.previous_peaklists.append(copy.deepcopy(self.peak_list_dictionary))

            self.peak_list_dictionary[self.selected_peaklist]["shift1"][
                self.selected_peak_indexes[0]
            ] = x
            self.peak_list_dictionary[self.selected_peaklist]["shift2"][
                self.selected_peak_indexes[0]
            ] = y

            intensity = self.find_new_intensity(x, y, self.selected_peaklist)

            self.peak_list_dictionary[self.selected_peaklist]["intensity"][
                self.selected_peak_indexes[0]
            ] = intensity

            

            self.update_peak_markers(self.selected_peaklist)
            self.AddToTable()

            self.active_move = False
            self.move_peaks_button.SetValue(False)
            self.main_frame.fig.canvas.mpl_disconnect(self.move_peak_connect)
            self.select_peak_connect = self.main_frame.fig.canvas.mpl_connect(
                "button_press_event", self.on_click_selectpeak
            )

    def on_press_movepeak(self, event):

        x, y = self.main_frame.ax.transData.inverted().transform((event.x, event.y))
        if x != None and y != None:
            if(len(self.previous_peaklists)>10):
                self.previous_peaklists.pop(0)
            self.previous_peaklists.append(copy.deepcopy(self.peak_list_dictionary))
            self.start_point_move = (x, y)
            self.x_init = copy.deepcopy(
                self.peak_list_dictionary[self.selected_peaklist]["shift1"]
            )
            self.y_init = copy.deepcopy(
                self.peak_list_dictionary[self.selected_peaklist]["shift2"]
            )

    def on_motion_movepeak(self, event):
        if self.start_point_move == None:
            return

        x, y = self.main_frame.ax.transData.inverted().transform((event.x, event.y))

        if x != None and y != None:

            # Update rectangle size
            x0, y0 = self.start_point_move
            x1, y1 = x, y
            x_change = x1 - x0
            y_change = y1 - y0

            for index in self.selected_peak_indexes:

                self.peak_list_dictionary[self.selected_peaklist]["shift1"][index] = (
                    self.x_init[index] + x_change
                )
                self.peak_list_dictionary[self.selected_peaklist]["shift2"][index] = (
                    self.y_init[index] + y_change
                )

                intensity = self.find_new_intensity(
                    self.x_init[index] + x_change,
                    self.y_init[index] + y_change,
                    self.selected_peaklist,
                )

                self.peak_list_dictionary[self.selected_peaklist]["intensity"][
                    index
                ] = intensity

            self.update_peak_markers(self.selected_peaklist)

    def on_release_movepeak(self, event):
        # self.on_motion_movepeak(event)
        self.AddToTable()
        self.start_point_move = None

    

    def get_current_data(self):
        if(self.main_frame.multiplot_mode==False):
            self.current_data = self.main_frame.nmrdata.data * self.main_frame.multiply_factor
            self.current_x_values = self.main_frame.new_x_ppms
            self.current_y_values = self.main_frame.new_y_ppms
        
        else:
            self.current_data = self.main_frame.values_dictionary[self.main_frame.active_plot_index]["z_data"] * self.main_frame.values_dictionary[self.main_frame.active_plot_index]["multiply factor"]
            self.current_x_values = self.main_frame.values_dictionary[self.main_frame.active_plot_index]["new_x_ppms"]
            self.current_y_values = self.main_frame.values_dictionary[self.main_frame.active_plot_index]["new_y_ppms"]

    def find_new_intensity(self, x, y, peaklist=None):

        (
            self.current_data,
            self.current_x_values,
            self.current_y_values,
        ) = self.find_plot_data(peaklist)

        x_index = np.argmin(np.abs(self.current_x_values - x))
        y_index = np.argmin(np.abs(self.current_y_values - y))
        c, r = y_index, x_index

        intensity = self.current_data[r][c]
        return intensity
    
    def OnFindPeaks(self, event):
        """
        If one peak is currently selected in the table, then zoom in to this
        peak and select it.
        Before doing this, the code will turn off all active toggled buttons from
        the Peak List frame.
        """
        if self.active_move == True:
            self.active_move = False
            self.move_peaks_button.SetValue(False)

            if self.active_select_peaks == True:
                self.main_frame.fig.canvas.mpl_disconnect(self.move_peak_press)
                self.main_frame.fig.canvas.mpl_disconnect(self.move_peak_motion)
                self.main_frame.fig.canvas.mpl_disconnect(self.move_peak_release)
            if self.active_select_peak == True:
                self.main_frame.fig.canvas.mpl_disconnect(self.move_peak_connect)
        elif self.active_add == True:
            self.active_add = False
            self.add_peaks_button.SetValue(False)
            self.main_frame.fig.canvas.mpl_disconnect(self.add_peak_connect)
            return

        elif self.active_select_peaks:
            self.active_select_peaks = False
            self.rect = None
            self.start_point = None
            self.select_peaks_button.SetValue(False)
            self.main_frame.fig.canvas.mpl_disconnect(self.select_press)
            self.main_frame.fig.canvas.mpl_disconnect(self.select_release)
            self.main_frame.fig.canvas.mpl_disconnect(self.select_motion)
            self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)

        elif self.active_select_peak:
            self.select_peak_button.SetValue(False)
            self.active_select_peak = False
            self.selected_peakname = ""
            self.main_frame.fig.canvas.mpl_disconnect(self.select_peak_connect)
            self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)

        row = self.grid.GetGridCursorRow()
        peak_name = self.grid.GetCellValue(row, 0)
        shift1 = float(self.grid.GetCellValue(row, 1))
        shift2 = float(self.grid.GetCellValue(row, 2))

        # Zoom in on grid selected peak and then select it in the plot.
        width = 0.1  # ppm
        height = 0.1  # ppm

        xvalues = self.main_frame.new_x_ppms


        if(shift1 > np.min(xvalues) and shift1 < np.max(xvalues)):
            xmin = shift1 - width
            xmax = shift1 + width
            ymin = shift2 - height
            ymax = shift2 + height
        else:
            xmin = shift2 - width
            xmax = shift2 + width
            ymin = shift1 - height
            ymax = shift1 + height

        self.main_frame.toolbar.push_current()

        self.zoom_to_region(self.main_frame.ax, [xmin, xmax], [ymin, ymax])
        self.main_frame.UpdateFrame()

        self.main_frame.toolbar.push_current()

    def OnDuplicatePeaklist(self, event):
        """
        Provide a FileDialog where the user can chose the name for
        the peaklist.
        The peaklist will then be saved and loaded
        """
        dlg = wx.FileDialog(self, "Select the peak list", wildcard="", style=wx.FD_SAVE)
        dlg.SetDirectory(os.getcwd())
        if dlg.ShowModal() == wx.ID_OK:
            peaklist_file = dlg.GetPath()
        else:
            dlg.Destroy()
            return


        self.write_peaklist(peaklist_file)

        self.AddPeaklist(peaklist_file, new_peaklist=True)

    def write_peaklist(self, peaklist_file):
        """
        Write the current peaklist into a file, in the order the table shows
        the peaks. The fitted linewidths are written as extra columns when any
        of the peaks have been fitted.
        """
        current_peaklist = self.current_peaklist_box.GetValue()
        if current_peaklist not in self.peak_list_dictionary:
            self.no_peaklist_message("Saving a peaklist")
            return 0

        self.pad_peaklist(current_peaklist)

        dictionary = self.peak_list_dictionary[current_peaklist]
        linewidths = self.find_linewidths_written(current_peaklist)
        written = 0

        with open(peaklist_file, "w") as file:
            # The header names the columns in the same way as the table
            file.write(
                " \t ".join(self.find_column_labels(current_peaklist, linewidths))
                + "\n"
            )

            for index in self.find_table_order(current_peaklist):
                line = "{} \t {} \t {} \t{}".format(
                    dictionary["peak_name"][index],
                    dictionary["shift1"][index],
                    dictionary["shift2"][index],
                    dictionary["intensity"][index],
                )
                if linewidths == True:
                    line += " \t " + self.find_linewidth_columns(
                        current_peaklist, index
                    )
                file.write(line + "\n")
                written += 1

        return written

    def find_save_location(self):
        """
        The folder and file name shown when saving the current peaklist, taken
        from the file the peaklist was read from. The current directory and an
        untitled file are used when the peaklist did not come from a file.
        """
        directory = pathlib.Path(os.getcwd())
        file_name = "Untitled.tab"

        selection = self.current_peaklist_box.GetSelection()
        if selection < 0 or selection >= len(self.peaklist_paths):
            return directory, file_name

        try:
            peaklist_path = pathlib.Path(self.peaklist_paths[selection]).expanduser()
        except TypeError:
            return directory, file_name

        if peaklist_path.name != "":
            file_name = peaklist_path.name

        # absolute rather than parent so that a peaklist given by name alone
        # is saved in the current directory
        parent = peaklist_path.absolute().parent
        if parent.is_dir() == True:
            directory = parent

        return directory, file_name

    def find_peaklist_file(self) -> str:
        """
        The file the current peaklist was read from, created as, or last saved
        as. An empty string is returned when the peaklist has no file yet.
        """
        selection = self.current_peaklist_box.GetSelection()
        if selection < 0 or selection >= len(self.peaklist_paths):
            return ""

        try:
            return str(pathlib.Path(self.peaklist_paths[selection]).expanduser())
        except TypeError:
            return ""

    def set_peaklist_file(self, peaklist_file):
        """
        Record the file the current peaklist is saved into, so that saving it
        again goes back to the same file.
        """
        selection = self.current_peaklist_box.GetSelection()
        if selection < 0 or selection >= len(self.peaklist_paths):
            return

        self.peaklist_paths[selection] = pathlib.Path(peaklist_file).absolute()

    def OnSave(self, event, save_after_picking=False):
        """
        Save the peaklist back into its own file, which is the file it was read
        from or created as, or the one it was last saved as. The user is asked
        for a file only when the peaklist does not have one.
        """
        if save_after_picking == True:
            self.save_peaklist(self.peaklist_name_box.GetValue(), False)
            return

        peaklist_file = self.find_peaklist_file()
        if peaklist_file == "":
            return self.OnSaveAs(event)

        self.save_peaklist(peaklist_file, True)

    def OnSaveAs(self, event):
        """
        Ask for a file to save the peaklist as, which becomes the file it is
        saved into from then on.
        """
        directory, file_name = self.find_save_location()

        dlg = wx.FileDialog(
            self,
            "Save the peaklist as",
            wildcard="Peaklists (*.list;*.tab;*.txt)|*.list;*.tab;*.txt|All files (*.*)|*.*",
            style=wx.FD_SAVE,
        )
        dlg.SetDirectory(str(directory))
        dlg.SetFilename(str(file_name))
        if dlg.ShowModal() == wx.ID_OK:
            peaklist_file = dlg.GetPath()
        else:
            dlg.Destroy()
            return
        dlg.Destroy()

        if peaklist_file == "" or os.path.isdir(peaklist_file) == True:
            dlg = wx.MessageDialog(
                None,
                "No name was given for the peaklist, so it has not been saved. "
                "Please try again and give the peaklist a name.",
                "Save peaklist",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return

        self.save_peaklist(peaklist_file, True)

    def save_peaklist(self, peaklist_file, tell_the_user):
        """
        Write the peaklist into a file, saying where it has gone and what went
        wrong if it could not be written.
        """
        try:
            number_of_peaks = self.write_peaklist(peaklist_file)
        except OSError as error:
            dlg = wx.MessageDialog(
                None,
                "The peaklist could not be saved as {} ({}).".format(
                    peaklist_file, error
                ),
                "Save peaklist",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return

        # Saving it again goes back to the same file, and it is that file which
        # is loaded when this window is opened next time
        previous = self.find_peaklist_file()
        self.set_peaklist_file(peaklist_file)
        if previous != "" and previous != peaklist_file:
            self.forget_peaklist(previous)
        self.remember_peaklist(peaklist_file, self.find_peaklist_dataset())

        if tell_the_user == False:
            return

        if number_of_peaks == 0:
            message = (
                "There are no peaks in the peaklist, so the file {} has been "
                "written empty.".format(peaklist_file)
            )
        else:
            message = "{} peaks have been saved as {}.".format(
                number_of_peaks, peaklist_file
            )

        dlg = wx.MessageDialog(None, message, "Save peaklist", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()

    def OnPickPeaks(self, event):
        """
        Pick peaks using nmrglue peak picking routines and then load this peaklist.
        """
        
        # See if the peaklist name is already in the current directory, and ask the user
        # if they wish to overwrite this.
        peaklist_name = self.peaklist_name_box.GetValue()
        if(peaklist_name in os.listdir()):
            message = 'The peaklist ({}) is already in the current directory. Would you like to overwrite this?'.format(peaklist_name)
            dlg = wx.MessageDialog(None, message, "Pick Peaks", wx.YES_NO)
            result=dlg.ShowModal()
            if(result == wx.ID_NO):
                dlg.Destroy()
                return
            dlg.Destroy()
        
        # Check to see the validity of the value in the threshold box.
        threshold_box_value = self.peak_picking_threshold_box.GetValue()
        try:
            threshold = float(threshold_box_value)
            if(threshold < 0 or threshold > 100):
                message = 'The value in the threshold box ({}) is not a number between 0 and 100, please correct this and try again.'.format(threshold_box_value)
                dlg = wx.MessageBox(message, "Pick Peaks", wx.OK)
                return

        except:
            message = 'The value in the threshold box ({}) is not a number, please correct this and try again.'.format(threshold_box_value)
            dlg = wx.MessageBox(message, "Pick Peaks", wx.OK)
            return
        

        self.xlabel = self.main_frame.ax.get_xlabel()
        self.ylabel = self.main_frame.ax.get_ylabel()

        if(self.main_frame.multiplot_mode==False):
            data = self.main_frame.nmrdata.data * self.main_frame.multiply_factor
            x = self.main_frame.new_x_ppms
            y = self.main_frame.new_y_ppms

        else:
            data = self.main_frame.values_dictionary[self.main_frame.active_plot_index]["z_data"] * self.main_frame.values_dictionary[self.main_frame.active_plot_index]["multiply factor"]
            x = self.main_frame.values_dictionary[self.main_frame.active_plot_index]["new_x_ppms"]
            y = self.main_frame.values_dictionary[self.main_frame.active_plot_index]["new_y_ppms"]

            # Check to see if the current selected plot is hidden or not
            continue_function = self.check_hidden_plot()
            if(continue_function==False):
                return


        threshold = float(threshold_box_value)/100 *np.max(data)

        algorithm_selection = self.peak_picking_algorithm_box.GetValue()
        sign_option = self.peak_picking_type.GetValue()
        if(algorithm_selection == 'thres' or algorithm_selection == 'thres-fast'):
            if(sign_option == 'Positive Peaks'):
                peaks = ng.peakpick.pick(data, pthres=threshold, algorithm=algorithm_selection, msep=[1,1])
            elif(sign_option == 'Negative Peaks'):
                peaks = ng.peakpick.pick(data, nthres=threshold, algorithm=algorithm_selection, msep=[1,1])
            else:
                peaks = ng.peakpick.pick(data, pthres=threshold, nthresh=threshold, algorithm=algorithm_selection, msep=[1,1])
        else:
            if(sign_option == 'Positive Peaks'):
                peaks = ng.peakpick.pick(data, pthres=threshold, algorithm=algorithm_selection)
            elif(sign_option == 'Negative Peaks'):
                peaks = ng.peakpick.pick(data, nthres=threshold, algorithm=algorithm_selection)
            else:
                peaks = ng.peakpick.pick(data, pthres=threshold, nthresh=threshold, algorithm=algorithm_selection)
        

        if(self.main_frame.multiplot_mode==False):
            x = self.main_frame.uc0.ppm(peaks["Y_AXIS"]) + self.main_frame.x_movement
            y = self.main_frame.uc1.ppm(peaks["X_AXIS"]) + self.main_frame.y_movement
        else:
            values = self.main_frame.values_dictionary[self.main_frame.active_plot_index]
            x = values['uc0'].ppm(peaks["Y_AXIS"]) + values['move x']
            y = values['uc1'].ppm(peaks["X_AXIS"]) + values['move y']

        picked_peak_array = []
        for i, xval in enumerate(x):
            picked_peak_array.append([xval, y[i]])

        dictionary = {}
        dictionary["peak_name"] = []
        dictionary["shift1"] = []
        dictionary["shift2"] = []
        dictionary['intensity'] = []
        for i, peak in enumerate(picked_peak_array):
            dictionary["peak_name"].append(str(i+1))
            dictionary["shift1"].append(peak[0])
            dictionary["shift2"].append(peak[1])
            dictionary['intensity'].append(peaks[i][-1])
        
        with open(peaklist_name, 'w') as file:
            # empty the current peaklist file with this name or create an empty peaklist file
            file.write('')
        
        p = pathlib.Path(peaklist_name)
        dirs = p.parts[-3:]
        file_name = p.parts[-1]
        last_directories_path = str(pathlib.Path(*dirs))
        peaklist = dictionary
        self.peak_list_dictionary[last_directories_path] = peaklist
        if self.peak_list_choices == [""]:
            self.peak_list_choices = [last_directories_path]
        else:
            self.peak_list_choices.append(last_directories_path)

        self.current_peaklist_box.SetItems(self.peak_list_choices)
        self.current_peaklist_box.SetSelection(len(self.peak_list_choices) - 1)

        self.names[last_directories_path] = [self.xlabel, self.ylabel]

        # The picked peaklist can be compared with the other peaklists
        self.update_comparison_boxes()

        # Remember where the picked peaklist is saved so that it stays in step
        # with the list of peaklists
        self.peaklist_paths.append(p)

        # The peaklist belongs to the spectrum which is currently selected
        self.set_peaklist_dataset(
            last_directories_path, getattr(self.main_frame, "active_plot_index", 0)
        )

        # and it is loaded again when this window is opened next time
        self.remember_peaklist(
            peaklist_name, self.find_peaklist_dataset(last_directories_path)
        )
        self.update_dataset_box()


        # self.include_2d_fit = True

        # if(self.include_2d_fit==True):
        #     # Perform a 2D fit of the spectra, then subtract the fit from the real data
        #     # Where there is a difference greater than the minimum threshold, then add an
        #     # extra peak to this with no restriction and see if adding another peak improves
        #     # things
        #     self.add_2d_fit()


        self.turn_off_togglebuttons()

        self.AddToTable()

        # Save the new peaklist
        self.OnSave(wx.EVT_BUTTON, save_after_picking=True)

        self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)


    def check_hidden_plot(self):
        current_plot = self.main_frame.plot_combobox.GetValue()
        if(current_plot in self.main_frame.hidden_list):
            dlg = wx.MessageDialog(
                    self,
                    "The selected spectrum is hidden and not visible. Uncheck the hide button in the select plot option box to see the selected spectrum. Would you like to continue anyway?"
                    ,
                    "Warning",
                    wx.OK,
            )
            res = dlg.ShowModal()
            if(res == wx.ID_NO):
                dlg.Destroy()
                return False
            
        return True
    


    def check_hidden_peaklist(self):
        current_peaklist = self.current_peaklist_box.GetValue()
        if(current_peaklist in self.hidden_peaklists):
            dlg = wx.MessageDialog(
                    self,
                    "The selected peaklist is hidden and not visible. Uncheck the hide button for this peaklist before continuing."
                    ,
                    "Warning",
                    wx.OK,
            )
            res = dlg.ShowModal()
            dlg.Destroy()
            return False
            
        return True
    

    def check_fit_window(self):
        """
        Check to see if a fit window plot is currently present. If it is
        tell the user to except the changes from the fit (if desired) and 
        close the window to continue.
        """
        for window in wx.GetTopLevelWindows():
            if isinstance(window, wx.Frame) and window.GetTitle() == "Fit peaks result":
            
                dlg = wx.MessageDialog(
                        self,
                        "There is a fit window \"Fit peaks result\" open. Please accept the changes from the fit (if desired) and close the \"Fit peaks result\" window and try again."
                        ,
                        "Warning",
                        wx.OK,
                )
                res = dlg.ShowModal()
                dlg.Destroy()
                return False
            
        return True

        

    def find_plot_data(self, peaklist=None):
        """
        The data of the spectrum which a peaklist belongs to, along with the
        chemical shifts of its two axes.
        """
        if self.main_frame.multiplot_mode == False:
            data = self.main_frame.nmrdata.data * self.main_frame.multiply_factor

            return data, self.main_frame.new_x_ppms, self.main_frame.new_y_ppms

        index = self.find_peaklist_dataset(peaklist)
        values = self.main_frame.values_dictionary[index]

        return (
            values["z_data"] * values["multiply factor"],
            values["new_x_ppms"],
            values["new_y_ppms"],
        )

    def find_peak_markers(self, peaklist=None):
        """
        The markers drawn on the spectrum for a peaklist. Hidden peaklists are
        not drawn, so the markers are found by the name of the peaklist rather
        than by where the peaklist comes in the list of peaklists. None is
        returned when the peaklist is not currently drawn.
        """
        if peaklist == None:
            peaklist = self.current_peaklist_box.GetValue()

        names = getattr(self.main_frame, "point_names", [])
        points = getattr(self.main_frame, "points", [])

        if peaklist not in names:
            return None

        index = names.index(peaklist)
        if index >= len(points):
            return None

        return points[index]

    def update_peak_markers(self, peaklist):
        """
        Move the markers drawn on the spectrum onto the positions the peaks of
        a peaklist now have, so that peaks follow the cursor as they are moved.
        """
        markers = self.find_peak_markers(peaklist)
        if markers == None:
            return

        markers.set_offsets(
            np.c_[
                self.peak_list_dictionary[peaklist]["shift1"],
                self.peak_list_dictionary[peaklist]["shift2"],
            ]
        )
        self.main_frame.UpdateFrame()

    def find_peak_index(self, x, y, peaklist=None):
        """
        The point of the spectrum (row, column) which a chemical shift
        position falls on. None is returned when the spectrum cannot be read.
        """
        try:
            data, x_values, y_values = self.find_plot_data(peaklist)
        except (KeyError, IndexError, AttributeError, TypeError):
            return None

        try:
            return (
                int(np.argmin(np.abs(x_values - x))),
                int(np.argmin(np.abs(y_values - y))),
            )
        except (ValueError, TypeError):
            return None

    def find_duplicate_peak(self, x, y, peaklist):
        """
        The name of a peak which is already in the peaklist at the position
        given, or None when there is no peak there. Two peaks are at the same
        position when they fall on the same point of the spectrum. Moving a
        new peak to its local maximum can put it on top of a peak which has
        already been picked.
        """
        if peaklist not in self.peak_list_dictionary:
            return None

        peaks = self.peak_list_dictionary[peaklist]
        new_index = self.find_peak_index(x, y, peaklist)

        for i, peakname in enumerate(peaks["peak_name"]):
            try:
                shift1 = peaks["shift1"][i]
                shift2 = peaks["shift2"][i]
            except IndexError:
                continue

            if new_index == None:
                # The spectrum cannot be read so only peaks at exactly the
                # same chemical shifts are treated as being on top of one another
                if shift1 == x and shift2 == y:
                    return peakname
                continue

            if self.find_peak_index(shift1, shift2, peaklist) == new_index:
                return peakname

        return None

    def check_duplicate_peak(self, x, y, peaklist) -> bool:
        """
        Warn the user when a new peak would be added on top of a peak which is
        already in the peaklist, returning whether the peak should be added.
        """
        duplicate = self.find_duplicate_peak(x, y, peaklist)
        if duplicate == None:
            return True

        dlg = wx.MessageDialog(
            None,
            "The peak {} in the peaklist {} is already at this position. Adding this peak will give two peaks on top of each other. Do you want to add it anyway?".format(
                duplicate, peaklist
            ),
            "Peaks at the Same Position",
            wx.YES_NO | wx.NO_DEFAULT | wx.ICON_QUESTION,
        )
        result = dlg.ShowModal()
        dlg.Destroy()

        return result == wx.ID_YES

    def find_local_maximum(self, x, y, peaklist=None):
        """
        Walk uphill from a position to the nearest local maximum of the
        spectrum, returning the chemical shifts of the maximum and the
        intensity there. The starting position is returned if the spectrum
        cannot be read.
        """
        try:
            data, x_values, y_values = self.find_plot_data(peaklist)
            rows, cols = data.shape
        except (KeyError, IndexError, AttributeError, TypeError):
            return x, y, None

        r = int(np.argmin(np.abs(x_values - x)))
        c = int(np.argmin(np.abs(y_values - y)))

        while True:
            # Get all 8 neighbors (including diagonals)
            neighbors = [
                (nr, nc)
                for nr in range(r - 1, r + 2)
                for nc in range(c - 1, c + 2)
                if (0 <= nr < rows and 0 <= nc < cols and (nr, nc) != (r, c))
            ]

            if len(neighbors) == 0:
                break

            # Find the neighbor with the highest value
            best_neighbor = max(neighbors, key=lambda pos: np.abs(data[pos[0], pos[1]]))

            # If the best neighbor is higher, move there
            if np.abs(data[best_neighbor[0], best_neighbor[1]]) > np.abs(data[r, c]):
                r, c = best_neighbor
            else:
                # No neighbor is higher local maximum reached
                break

        return x_values[r], y_values[c], data[r][c]

    def OnFindLocalMaximum(self, event):
        """
        Moves the selected peaks to their nearest local maximum in 2D data.
        """

        if(len(self.previous_peaklists)>10):
            self.previous_peaklists.pop(0)
        self.previous_peaklists.append(copy.deepcopy(self.peak_list_dictionary))

        if(self.main_frame.multiplot_mode==False):
            # Check to see if the current selected plot is hidden or not
            continue_function = self.check_hidden_plot()
            if(continue_function==False):
                return

        for k, peak_index in enumerate(self.selected_peak_indexes):

            x = self.peak_list_dictionary[self.selected_peaklist]["shift1"][
                self.selected_peak_indexes[k]
            ]
            y = self.peak_list_dictionary[self.selected_peaklist]["shift2"][
                self.selected_peak_indexes[k]
            ]

            xvalue, yvalue, intensity = self.find_local_maximum(
                x, y, self.selected_peaklist
            )
            if intensity == None:
                continue

            self.peak_list_dictionary[self.selected_peaklist]["shift1"][
                self.selected_peak_indexes[k]
            ] = xvalue
            self.peak_list_dictionary[self.selected_peaklist]["shift2"][
                self.selected_peak_indexes[k]
            ] = yvalue
            self.peak_list_dictionary[self.selected_peaklist]["intensity"][
                self.selected_peak_indexes[k]
            ] = intensity

        self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)
        self.AddToTable()






class PeakListWindow3D(PeakModeButtons, wx.Frame):
    def __init__(self, title, parent):
        """
        This class contains all the information relating to loading in
        3D peaklists when in the SpinBore window
        """
        self.main_frame = parent
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        width = min(1280, self.monitorWidth)
        height = min(460, self.monitorHeight)
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.panel_peaklist = wx.Panel(self, -1)
        self.main_peaklist_sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_peaklist_sizer)

        self.set_initial_values()
        self.follow_peak_modes()
        self.make_peaklist_window()
        self.watch_navigation_toolbar()
        self.update_mode_buttons()

        # Clicking along the bore dimension moves the selected peak there
        self.bore_click_connect = self.main_frame.fig_bore.canvas.mpl_connect(
            "button_press_event", self.on_click_bore_dimension
        )

        self.Bind(wx.EVT_CLOSE, self.OnClose)

        # The distances at which peaks can be confused follow the data
        self.update_closeness_boxes()

        # The peaklists which were shown last time are loaded again
        self.load_remembered_peaklists()

        self.Show()
        # self.AddPeakListBrowser()

    def set_initial_values(self):
        """
        Setting initial values such as the peak list colour choices
        """
        self.peak_list_choices = [""]
        self.initial_peak_list_colours = ["black"]
        self.selected_colour = "darkviolet"
        self.peak_list_dictionary = {}
        self.selected_peakname = ""
        self.selected_peaklist = ""
        self.selected_peak_indexes = ["N/A"]
        self.bore_xdim = 'shift1'

        # Set when a peak has been picked out to be removed
        self.remove_peak = False

        # The peaklist as it was when it was last read or saved, so that
        # changes which have not been saved can be noticed
        self.saved_peaklist = {}

        # The file the peaklist was read from or created as. Only the last few
        # parts of the path are shown in the window, so the whole path is kept
        # here for saving the peaklist back to the file it came from
        self.peaklist_path = ""

        # The peaklist being used as the reference plane, along with the file
        # it was read from, so that the window can say which one is in use
        self.reference_peaklist = {}
        self.reference_peaklist_name = ""
        self.reference_peaklist_path = ""

        # Flags showing whether a given button is active or not
        self.reference_plane = False
        self.active_add = False
        self.active_select_peak = False
        self.active_select_peaks = False
        self.active_remove = False
        self.active_move = False
        self.active_find = False
        self.active_movez = False

        self.rect = None
        self.start_point = None
        self.start_point_move = None

        self.old_key = None
        self.old_num = None

    def make_peaklist_window(self):
        """
        This window will have the following:
        - a button to add peaklists
        a selection of buttons associated with picking peaks using nmrglue
        - buttons to toggle add peak(s), select peak, select region, remove peak(s), move peak(s), find peak
        """

        self.add_peaklist_button = wx.Button(self, label="Add peaklist")
        self.add_peaklist_button.Bind(wx.EVT_BUTTON, self.AddPeakListBrowser)

        self.remove_peaklist_button = wx.Button(self, label="Remove Peaklist")
        self.remove_peaklist_button.Bind(wx.EVT_BUTTON, self.OnRemovePeakList)
        self.remove_peaklist_button.SetToolTip(
            "Remove the peaklist from this window. The peaklist file itself is "
            "not deleted, and the peaklist is not loaded again when this window "
            "is opened next time."
        )

        self.peaklist_selection_text = wx.StaticText(self, -1, "Selected Peaklist:")

        self.current_peaklist_box = wx.TextCtrl(
            self, -1, value='', size=(250, 20), style = wx.TE_READONLY
        )


        self.row2_label = wx.StaticBox(
            self,
            -1,
            "Manipulate Peaklists: (shorcuts for Mac - cmd+key)",
        )
        self.row2 = wx.StaticBoxSizer(self.row2_label, wx.VERTICAL)
        self.row2_1 = wx.BoxSizer(wx.HORIZONTAL)
        self.row2_2 = wx.BoxSizer(wx.HORIZONTAL)

        self.add_peaks_button = wx.ToggleButton(self.row2_label, label="Add Peaks (a)")
        self.add_peaks_button.Bind(wx.EVT_TOGGLEBUTTON, self.OnAddPeaks)
        ID_BUTTON_a = wx.NewIdRef()

        self.select_peak_button = wx.ToggleButton(self.row2_label, label="Select Peak (s)")
        self.select_peak_button.Bind(wx.EVT_TOGGLEBUTTON, self.OnSelectPeak)
        ID_BUTTON_s = wx.NewIdRef()

        self.select_peaks_button = wx.ToggleButton(
            self.row2_label, label="Select Peak Group (g)"
        )
        self.select_peaks_button.Bind(wx.EVT_TOGGLEBUTTON, self.OnSelectPeaks)
        ID_BUTTON_g = wx.NewIdRef()

        self.remove_peaks_button = wx.Button(self.row2_label, label="Remove Peaks (r)")
        self.remove_peaks_button.Bind(wx.EVT_BUTTON, self.OnRemovePeaks)
        ID_BUTTON_r = wx.NewIdRef()

        self.find_peak_button = wx.Button(self.row2_label, label="Find Peak (f)")
        self.find_peak_button.Bind(wx.EVT_BUTTON, self.OnFindPeaks)
        ID_BUTTON_f = wx.NewIdRef()

        self.move_peaks_button = wx.ToggleButton(self.row2_label, label="Move Peak x/y (m)")
        self.move_peaks_button.Bind(wx.EVT_TOGGLEBUTTON, self.OnMovePeak)
        ID_BUTTON_m = wx.NewIdRef()

        self.move_peaks_bore_button = wx.ToggleButton(self.row2_label, label="Move Peak z (z)")
        self.move_peaks_bore_button.Bind(wx.EVT_TOGGLEBUTTON, self.OnMovePeakz)
        ID_BUTTON_z = wx.NewIdRef()

        ID_BUTTON_k = wx.NewIdRef()
        ID_BUTTON_j = wx.NewIdRef()

        # Creating an accelerator table for keyboard shortcuts for the buttons
        accelerator_table = wx.AcceleratorTable(
            [
                (wx.ACCEL_CTRL, ord("a"), ID_BUTTON_a),
                (wx.ACCEL_CTRL, ord("r"), ID_BUTTON_r),
                (wx.ACCEL_CTRL, ord("g"), ID_BUTTON_g),
                (wx.ACCEL_CTRL, ord("f"), ID_BUTTON_f),
                (wx.ACCEL_CTRL, ord("m"), ID_BUTTON_m),
                (wx.ACCEL_CTRL, ord("z"), ID_BUTTON_z),
                (wx.ACCEL_CTRL, ord("s"), ID_BUTTON_s),
                (wx.ACCEL_CTRL, ord("k"), ID_BUTTON_k),
                (wx.ACCEL_CTRL, ord("j"), ID_BUTTON_j),
            ]
        )

        self.SetAcceleratorTable(accelerator_table)
        self.main_frame.SetAcceleratorTable(accelerator_table)
        self.Bind(wx.EVT_MENU, self.OnAddPeaks, id=ID_BUTTON_a)
        self.Bind(wx.EVT_MENU, self.OnRemovePeaks, id=ID_BUTTON_r)
        self.Bind(wx.EVT_MENU, self.OnMovePeak, id=ID_BUTTON_m)
        self.Bind(wx.EVT_MENU, self.OnMovePeakz, id=ID_BUTTON_z)
        self.Bind(wx.EVT_MENU, self.OnFindPeaks, id=ID_BUTTON_f)
        self.Bind(wx.EVT_MENU, self.OnSelectPeak, id=ID_BUTTON_s)
        self.Bind(wx.EVT_MENU, self.OnFindLocalMaximum, id=ID_BUTTON_k)
        self.Bind(wx.EVT_MENU, self.OnFindLocalMaximumBore, id=ID_BUTTON_j)
        self.Bind(wx.EVT_MENU, self.OnSelectPeaks, id=ID_BUTTON_g)

        self.main_frame.Bind(wx.EVT_MENU, self.OnAddPeaks, id=ID_BUTTON_a)
        self.main_frame.Bind(wx.EVT_MENU, self.OnSelectPeaks, id=ID_BUTTON_g)
        self.main_frame.Bind(wx.EVT_MENU, self.OnRemovePeaks, id=ID_BUTTON_r)
        self.main_frame.Bind(wx.EVT_MENU, self.OnMovePeak, id=ID_BUTTON_m)
        self.main_frame.Bind(wx.EVT_MENU, self.OnMovePeakz, id=ID_BUTTON_z)
        self.main_frame.Bind(wx.EVT_MENU, self.OnFindPeaks, id=ID_BUTTON_f)
        self.main_frame.Bind(wx.EVT_MENU, self.OnSelectPeak, id=ID_BUTTON_s)
        self.main_frame.Bind(wx.EVT_MENU, self.OnFindLocalMaximum, id=ID_BUTTON_k)
        self.main_frame.Bind(wx.EVT_MENU, self.OnFindLocalMaximumBore, id=ID_BUTTON_j)

        self.save_peaks_button = wx.Button(self.row2_label, label="Save")
        self.save_peaks_button.Bind(wx.EVT_BUTTON, self.OnSave)
        self.save_peaks_button.SetToolTip(
            "Save the peaklist into its own file, which is named in the "
            "selected peaklist box above."
        )

        self.resolve_button = wx.Button(self.row2_label, label="Resolve Bore Peaks")
        self.resolve_button.Bind(wx.EVT_BUTTON, self.OnResolveAmbiguities)
        self.resolve_button.SetToolTip(
            "Show the peaks whose position down the bore dimension is not "
            "certain, so that one of the maxima which were found can be chosen "
            "or the peak left as it is."
        )

        self.reanalyse_button = wx.Button(self.row2_label, label="Re-analyse Bore")
        self.reanalyse_button.Bind(wx.EVT_BUTTON, self.OnReanalyseBore)
        self.reanalyse_button.SetToolTip(
            "Work out again which maxima down the bore dimension belong to which "
            "peak, using the expected number and the distance which are set now. "
            "The positions of the peaks in the plane are not changed."
        )

        self.refresh_intensities_button = wx.Button(
            self.row2_label, label="Refresh Intensities"
        )
        self.refresh_intensities_button.Bind(
            wx.EVT_BUTTON, self.OnRefreshIntensities
        )
        self.refresh_intensities_button.SetToolTip(
            "Read the intensity of every peak from the spectrum again, leaving "
            "the positions of the peaks alone. This is used when a peaklist "
            "picked on one spectrum is loaded onto another one."
        )

        self.save_peaks_as_button = wx.Button(self.row2_label, label="Save As")
        self.save_peaks_as_button.Bind(wx.EVT_BUTTON, self.OnSaveAs)
        self.save_peaks_as_button.SetToolTip(
            "Save the peaklist into a different file, which it is then saved "
            "into from then on."
        )

        # Moving the selected peaks onto the nearest maximum, in the plane
        # which is shown and along the bore dimension
        self.move_to_local_max = wx.Button(
            self.row2_label, label="Move to local max x/y (k)"
        )
        self.move_to_local_max.Bind(wx.EVT_BUTTON, self.OnFindLocalMaximum)

        self.move_to_local_max_bore = wx.Button(
            self.row2_label, label="Move to local max z (j)"
        )
        self.move_to_local_max_bore.Bind(wx.EVT_BUTTON, self.OnFindLocalMaximumBore)

        self.add_at_local_max_box = wx.CheckBox(
            self.row2_label, -1, "Add peaks at local maximum (x/y and z)"
        )
        self.add_at_local_max_box.SetValue(False)
        self.add_at_local_max_box.SetToolTip(
            "A peak which is added goes onto the nearest maximum of the plane "
            "which is shown, and is given the chemical shift of the largest "
            "point down the bore dimension at that position rather than zero."
        )

        self.snap_bore_box = wx.CheckBox(
            self.row2_label, -1, "Move to local maximum in z when clicking the bore"
        )
        self.snap_bore_box.SetValue(False)
        self.snap_bore_box.SetToolTip(
            "When a peak is selected, clicking along the bore dimension moves "
            "it there. With this ticked the peak goes onto the maximum nearest "
            "to the click instead of exactly where it was clicked."
        )

        self.row2_1.AddSpacer(5)
        self.row2_1.Add(self.add_peaks_button)
        self.row2_1.AddSpacer(5)
        self.row2_1.Add(self.select_peak_button)
        self.row2_1.AddSpacer(5)
        self.row2_1.Add(self.select_peaks_button)
        self.row2_1.AddSpacer(5)
        self.row2_1.Add(self.move_peaks_button)
        self.row2_1.AddSpacer(5)
        self.row2_1.Add(self.move_peaks_bore_button)
        self.row2_1.AddSpacer(5)
        self.row2_1.Add(self.remove_peaks_button)
        self.row2_1.AddSpacer(5)
        self.row2_1.Add(self.find_peak_button)
        self.row2_1.AddSpacer(5)
        self.row2_1.Add(self.save_peaks_button)
        self.row2_1.AddSpacer(5)
        self.row2_1.Add(self.save_peaks_as_button)

        self.row2_2.AddSpacer(5)
        self.row2_2.Add(self.move_to_local_max)
        self.row2_2.AddSpacer(5)
        self.row2_2.Add(self.move_to_local_max_bore)
        self.row2_2.AddSpacer(5)
        self.row2_2.Add(self.refresh_intensities_button)
        self.row2_2.AddSpacer(5)
        self.row2_2.Add(self.resolve_button)
        self.row2_2.AddSpacer(5)
        self.row2_2.Add(self.reanalyse_button)
        self.row2_2.AddSpacer(10)
        self.row2_2.Add(self.add_at_local_max_box, 0, wx.ALIGN_CENTER_VERTICAL)
        self.row2_2.AddSpacer(10)
        self.row2_2.Add(self.snap_bore_box, 0, wx.ALIGN_CENTER_VERTICAL)

        self.row2.Add(self.row2_1)
        self.row2.AddSpacer(5)
        self.row2.Add(self.row2_2)


        self.row_pickpeaks_label = wx.StaticBox(self, -1, "Peak Picking (nmrglue):")
        self.row_pickpeaks = wx.StaticBoxSizer(self.row_pickpeaks_label, wx.VERTICAL)

        self.peak_picking_threshold_text = wx.StaticText(self.row_pickpeaks_label,-1,"Threshold (% of maximum):")
        self.peak_picking_threshold_box = wx.TextCtrl(self.row_pickpeaks_label,value='10.0',
                size=(30, 20))
        
        self.peak_picking_type_text = wx.StaticText(self.row_pickpeaks_label,-1,"Option:")
        types = ['Positive Peaks', 'Negative Peaks', 'Positive + Negative Peaks']
        self.peak_picking_type = wx.ComboBox(self.row_pickpeaks_label, choices = types, style=wx.CB_READONLY)
        
        self.peak_picking_algorithm_text = wx.StaticText(self.row_pickpeaks_label,-1,"Algorithm:")
        algorithms = ['thres', 'thres-fast', 'downward', 'connected']
        self.peak_picking_algorithm_box = wx.ComboBox(self.row_pickpeaks_label, choices=algorithms, style=wx.CB_READONLY)

        self.reference_plane_button = wx.ToggleButton(self.row_pickpeaks_label,-1,"Load reference plane (optional)")
        self.reference_plane_button.Bind(wx.EVT_TOGGLEBUTTON, self.OnLoadReferencePlane)

        self.reference_method_box = wx.ComboBox(
            self.row_pickpeaks_label,
            choices=[TRACE_METHOD, TEMPLATE_METHOD],
            style=wx.CB_READONLY,
            size=(210, 24),
        )
        self.reference_method_box.SetValue(TRACE_METHOD)
        self.reference_method_box.SetToolTip(
            "How the reference plane is used. Picking down the bore of each "
            "reference peak looks along the bore dimension at each of them in "
            "turn. Picking in 3D and matching finds the peaks of the whole "
            "spectrum first and then relates them back to the reference peaks "
            "by the two dimensions they share, which also says which peaks "
            "belong to nothing in the reference plane."
        )

        self.remove_reference_button = wx.Button(
            self.row_pickpeaks_label, label="Remove reference plane"
        )
        self.remove_reference_button.Bind(
            wx.EVT_BUTTON, self.OnRemoveReferencePlane
        )
        self.remove_reference_button.SetToolTip(
            "Stop using the reference plane. Peaks are then picked from the "
            "spectrum alone, in all three dimensions, and the peaklist which "
            "was loaded as the reference plane is left alone."
        )

        # Says whether a reference plane is being used and which peaklist it is,
        # as the peaks which come out depend on it
        self.reference_plane_text = wx.StaticText(self.row_pickpeaks_label, -1, "")

        self.experiment_text = wx.StaticText(
            self.row_pickpeaks_label, -1, "Experiment:"
        )
        self.experiment_box = wx.ComboBox(
            self.row_pickpeaks_label,
            choices=[name for name, expected in BORE_EXPERIMENTS],
            style=wx.CB_READONLY,
        )
        self.experiment_box.SetValue(BORE_EXPERIMENTS[0][0])
        self.experiment_box.Bind(wx.EVT_COMBOBOX, self.OnExperimentSelection)
        self.experiment_box.SetToolTip(
            "Choosing the experiment fills in how many peaks each peak of the "
            "reference plane is expected to have down the bore dimension. With "
            "no experiment chosen the number can be typed in."
        )

        self.expected_peaks_text = wx.StaticText(
            self.row_pickpeaks_label, -1, "Expected peaks down the bore:"
        )
        self.expected_peaks_box = wx.TextCtrl(
            self.row_pickpeaks_label, value="0", size=(30, 20)
        )
        self.expected_peaks_box.SetToolTip(
            "How many maxima down the bore dimension each peak of the reference "
            "plane is expected to have, which depends on the experiment. A peak "
            "with fewer than this is noted rather than made up, as a peak can be "
            "missing for good reasons such as a proline or the end of the chain. "
            "Zero keeps every maximum which is found."
        )

        self.closeness_text = wx.StaticText(
            self.row_pickpeaks_label, -1, "Peaks closer than:"
        )
        self.closeness_box_1 = wx.TextCtrl(
            self.row_pickpeaks_label, value="", size=(50, 20)
        )
        self.closeness_box_2 = wx.TextCtrl(
            self.row_pickpeaks_label, value="", size=(50, 20)
        )
        self.closeness_units_text = wx.StaticText(
            self.row_pickpeaks_label, -1, "ppm may be confused"
        )
        for box in [self.closeness_box_1, self.closeness_box_2]:
            box.SetToolTip(
                "How close two peaks have to be in the plane before a maximum "
                "down the bore of one of them could belong to the other. The "
                "values are worked out from the spacing of the data and can be "
                "changed."
            )

        self.peaklist_name_text = wx.StaticText(self.row_pickpeaks_label,-1,"Peaklist name:")
        self.peaklist_name_box = wx.TextCtrl(self.row_pickpeaks_label,value='peaks_nmrglue.list',
                size=(200, 20))

        self.peak_pick_button = wx.Button(self.row_pickpeaks_label, label='Peak Pick')
        self.peak_pick_button.Bind(wx.EVT_BUTTON, self.OnPickPeaks)


        self.row_pickpeaks1 = wx.BoxSizer(wx.HORIZONTAL)
        self.row_pickpeaks2 = wx.BoxSizer(wx.HORIZONTAL)
        
        self.row_pickpeaks1.Add(self.peak_picking_threshold_text)
        self.row_pickpeaks1.AddSpacer(5)
        self.row_pickpeaks1.Add(self.peak_picking_threshold_box)
        self.row_pickpeaks1.AddSpacer(10)


        self.row_pickpeaks1.Add(self.reference_plane_button)
        self.row_pickpeaks1.AddSpacer(5)
        self.row_pickpeaks1.Add(self.reference_method_box, 0, wx.ALIGN_CENTER_VERTICAL)
        self.row_pickpeaks1.AddSpacer(10)
        self.row_pickpeaks1.Add(self.peaklist_name_text)
        self.row_pickpeaks1.AddSpacer(10)
        self.row_pickpeaks1.Add(self.peaklist_name_box)
        self.row_pickpeaks1.AddSpacer(5)
        self.row_pickpeaks1.Add(self.peak_pick_button)


        self.row_pickpeaks2.Add(self.peak_picking_type_text)
        self.row_pickpeaks2.AddSpacer(5)
        self.row_pickpeaks2.Add(self.peak_picking_type)
        self.row_pickpeaks2.AddSpacer(10)
        self.row_pickpeaks2.Add(self.peak_picking_algorithm_text)
        self.row_pickpeaks2.AddSpacer(5)
        self.row_pickpeaks2.Add(self.peak_picking_algorithm_box)
        self.row_pickpeaks2.AddSpacer(10)
        self.row_pickpeaks2.Add(self.experiment_text, 0, wx.ALIGN_CENTER_VERTICAL)
        self.row_pickpeaks2.AddSpacer(5)
        self.row_pickpeaks2.Add(self.experiment_box)
        self.row_pickpeaks2.AddSpacer(10)
        self.row_pickpeaks2.Add(self.expected_peaks_text, 0, wx.ALIGN_CENTER_VERTICAL)
        self.row_pickpeaks2.AddSpacer(5)
        self.row_pickpeaks2.Add(self.expected_peaks_box)
        self.row_pickpeaks2.AddSpacer(10)
        self.row_pickpeaks2.Add(self.closeness_text, 0, wx.ALIGN_CENTER_VERTICAL)
        self.row_pickpeaks2.AddSpacer(5)
        self.row_pickpeaks2.Add(self.closeness_box_1)
        self.row_pickpeaks2.AddSpacer(5)
        self.row_pickpeaks2.Add(self.closeness_box_2)
        self.row_pickpeaks2.AddSpacer(5)
        self.row_pickpeaks2.Add(self.closeness_units_text, 0, wx.ALIGN_CENTER_VERTICAL)

        self.row_pickpeaks3 = wx.BoxSizer(wx.HORIZONTAL)
        self.row_pickpeaks3.Add(self.reference_plane_text, 0, wx.ALIGN_CENTER_VERTICAL)
        self.row_pickpeaks3.AddSpacer(10)
        self.row_pickpeaks3.Add(self.remove_reference_button, 0, wx.ALIGN_CENTER_VERTICAL)

        self.row_pickpeaks.Add(self.row_pickpeaks1, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.row_pickpeaks.AddSpacer(10)
        self.row_pickpeaks.Add(self.row_pickpeaks2, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.row_pickpeaks.AddSpacer(5)
        self.row_pickpeaks.Add(self.row_pickpeaks3, 0, wx.ALIGN_CENTER_HORIZONTAL)

        self.update_reference_plane_text()


        self.row1 = wx.BoxSizer(wx.HORIZONTAL)
        self.row1.Add(self.add_peaklist_button)
        self.row1.AddSpacer(5)
        self.row1.Add(self.remove_peaklist_button)
        self.row1.AddSpacer(10)
        self.row1.Add(self.peaklist_selection_text)
        self.row1.AddSpacer(5)
        self.row1.Add(self.current_peaklist_box)

        self.main_peaklist_sizer.AddSpacer(10)
        self.main_peaklist_sizer.Add(
            self.row1, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 5
        )

        self.main_peaklist_sizer.AddSpacer(10)
        self.main_peaklist_sizer.Add(
            self.row_pickpeaks, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 5
        )
        self.main_peaklist_sizer.AddSpacer(10)
        self.main_peaklist_sizer.Add(
            self.row2, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 5
        )

        # Then have a table of the currently loaded peaklist (originally blank)

        self.row3_label = wx.StaticBox(self, -1, "Peaklist Table:")
        self.row3 = wx.StaticBoxSizer(self.row3_label, wx.HORIZONTAL)

        self.grid = gridlib.Grid(self.row3_label)
        self.grid.CreateGrid(5, 7)

        self.grid.SetColLabelValue(0, "Peak name")
        self.grid.SetColLabelValue(1, "Shift 1 (ppm)")
        self.grid.SetColLabelValue(2, "Shift 2 (ppm)")
        self.grid.SetColLabelValue(3, "Shift 3 (ppm)")
        self.grid.SetColLabelValue(4, "Intensity")
        self.grid.SetColLabelValue(5, "Ambiguity")
        self.grid.SetColLabelValue(6, "Alternatives (ppm)")

        # Bind event when cell value changes
        self.grid.Bind(gridlib.EVT_GRID_EDITOR_SHOWN, self.on_begin_edit)
        self.grid.Bind(gridlib.EVT_GRID_CELL_CHANGED, self.on_cell_changed)

        
        self.row3.Add(self.grid, proportion=1, flag=wx.EXPAND | wx.ALL, border=5)
        self.main_peaklist_sizer.AddSpacer(10)
        self.main_peaklist_sizer.Add(self.row3, 1, wx.EXPAND | wx.ALL, 5)

        self.Layout()
        self.Refresh()
        total_width = int(self.grid.GetClientSize().width * 0.8)
        col_count = self.grid.GetNumberCols()
        if col_count > 0:
            col_width = int(total_width // col_count)
            for c in range(col_count):
                self.grid.SetColSize(c, col_width)


    def OnRemovePeakList(self, event):
        """
        Remove the loaded peaklist from this window, leaving the file it came
        from on disk. The peaklist is not loaded again when the window is
        opened next time.
        """
        peaklist = self.find_current_peaklist()
        if peaklist == None:
            self.no_peaklist_message("Removing a peaklist")
            return

        dlg = wx.MessageDialog(
            self,
            "Remove the peaklist {} from this window? The peaklist file itself "
            "will not be deleted.".format(peaklist),
            "Removing a peaklist",
            wx.YES_NO | wx.ICON_QUESTION,
        )
        result = dlg.ShowModal()
        dlg.Destroy()
        if result != wx.ID_YES:
            return

        self.turn_off_togglebuttons()

        # Stop it coming back when this window is opened next time
        self.forget_peaklist(getattr(self, "peaklist_path", peaklist))
        self.forget_peaklist(peaklist)

        # Forget everything which was held for this peaklist
        self.peak_list_dictionary.pop(peaklist, None)

        # Only a peaklist which has been picked has names of its own
        names = getattr(self, "names", None)
        if names != None:
            names.pop(peaklist, None)

        if peaklist in self.peak_list_choices:
            self.peak_list_choices.remove(peaklist)
        if len(self.peak_list_choices) == 0:
            self.peak_list_choices = [""]

        self.peak_list = self.peak_list_choices[-1]
        self.current_peaklist_box.SetValue("")
        self.peaklist_path = ""
        self.saved_peaklist = {}

        self.selected_peaklist = ""
        self.selected_peakname = ""
        self.selected_peak_indexes = ["N/A"]
        self.main_frame.selected_bore_peaks = []

        self.AddToTable()
        self.main_frame.OnBoreSlider(wx.EVT_BUTTON)

    def load_remembered_peaklists(self):
        """
        Load the peaklists which were shown the last time the peaks window was
        open. The bore window draws the peaks of this window, which is not the
        one it holds until this window has been made.
        """
        self.main_frame.peak_lists3D = self

        PeakModeButtons.load_remembered_peaklists(self)

    def analyse_bore_positions(self, positions):
        """
        Pick the maxima down the bore dimension at each position in the plane
        and share them out between the positions which are close together.
        """
        viewer = self.find_viewer()
        data = viewer.nmrdata.data

        try:
            threshold = float(self.peak_picking_threshold_box.GetValue()) / 100 * np.max(
                data
            )
        except (ValueError, TypeError):
            dlg = wx.MessageDialog(
                self,
                "The value in the threshold box is not a number, please correct "
                "it and try again.",
                "Analysing the bore dimension",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return None

        algorithm = self.peak_picking_algorithm_box.GetValue()
        sign_option = self.peak_picking_type.GetValue()

        bore_shifts = np.array(viewer.uc2.ppm(np.arange(data.shape[0])))

        trace_peaks = []
        for shift1, shift2 in positions:
            index3, index1, index2 = self.find_peak_indexes(shift1, shift2, 0)
            plane = (index1, index2)
            trace = data[:, index1, index2]

            found = None
            try:
                if algorithm in ["thres", "thres-fast"]:
                    if sign_option == "Negative Peaks":
                        found = ng.peakpick.pick(
                            trace, nthres=threshold, algorithm=algorithm, msep=1
                        )
                    elif sign_option == "Positive Peaks":
                        found = ng.peakpick.pick(
                            trace, pthres=threshold, algorithm=algorithm, msep=1
                        )
                    else:
                        found = ng.peakpick.pick(
                            trace,
                            pthres=threshold,
                            nthresh=threshold,
                            algorithm=algorithm,
                            msep=1,
                        )
                else:
                    if sign_option == "Negative Peaks":
                        found = ng.peakpick.pick(
                            trace, nthres=threshold, algorithm=algorithm
                        )
                    elif sign_option == "Positive Peaks":
                        found = ng.peakpick.pick(
                            trace, pthres=threshold, algorithm=algorithm
                        )
                    else:
                        found = ng.peakpick.pick(
                            trace,
                            pthres=threshold,
                            nthresh=threshold,
                            algorithm=algorithm,
                        )
            except Exception:
                # The peak picker found nothing down the bore here
                found = None

            candidates = self.find_trace_candidates(
                trace, bore_shifts, found, threshold
            )

            trace_peaks.append(
                {
                    "plane": plane,
                    "candidates": candidates,
                    "projection": self.find_projection_intensity(shift1, shift2),
                }
            )

        return self.find_bore_assignments(data, trace_peaks, positions)

    def apply_bore_assignments(self, peaklist, positions, groups, assignments):
        """
        Put the positions down the bore which the analysis has worked out into
        the peaklist, keeping the peaks where they are in the plane.
        """
        dictionary = self.peak_list_dictionary[peaklist]
        self.pad_notes(peaklist)
        self.peak_candidates = getattr(self, "peak_candidates", {})

        uncertain = 0

        for i, assignment in enumerate(assignments):
            kept = assignment["kept"]
            others = bore_candidates.find_alternatives_text(assignment["others"])
            notes = [
                note for note in assignment["notes"]
                if note != bore_candidates.AMBIGUOUS
            ]

            for position, index in enumerate(groups[i]):
                name = dictionary["peak_name"][index]
                self.peak_candidates[name] = assignment["candidates"]

                if position < len(kept):
                    candidate = kept[position]
                    dictionary["shift3"][index] = candidate["shift"]
                    dictionary["intensity"][index] = candidate["intensity"]
                    peak_notes = list(notes)
                    if candidate["ambiguous"] == True:
                        peak_notes = [bore_candidates.AMBIGUOUS] + peak_notes
                else:
                    # This peak has no maximum of its own down the bore
                    dictionary["shift3"][index] = 0
                    dictionary["intensity"][index] = 0
                    peak_notes = list(notes)

                dictionary["ambiguity"][index] = bore_candidates.find_ambiguity_text(
                    peak_notes
                )
                dictionary["alternatives"][index] = others

                if dictionary["ambiguity"][index] != "":
                    uncertain += 1

        self.AddToTable()
        self.main_frame.OnBoreSlider(wx.EVT_BUTTON)

        return uncertain

    def pad_notes(self, peaklist):
        """
        Make sure every peak has an entry in the two columns which say how sure
        its position down the bore is.
        """
        dictionary = self.peak_list_dictionary[peaklist]
        number_of_peaks = len(dictionary.get("peak_name", []))

        for key in ["ambiguity", "alternatives"]:
            if key not in dictionary:
                dictionary[key] = []
            while len(dictionary[key]) < number_of_peaks:
                dictionary[key].append("")

    def find_uncertain_peaks(self) -> list:
        """
        The peaks whose position down the bore dimension is not certain: the
        ones where the choice between neighbouring peaks was too close to call,
        and the ones which found fewer maxima than expected.
        """
        peaklist = self.find_current_peaklist()
        if peaklist == None:
            return []

        dictionary = self.peak_list_dictionary[peaklist]
        uncertain = []

        for index in self.find_table_order(peaklist):
            note = self.find_note(peaklist, "ambiguity", index)
            if note == "":
                continue

            name = dictionary["peak_name"][index]
            candidates = list(getattr(self, "peak_candidates", {}).get(name, []))

            if len(candidates) == 0:
                # A peaklist which has been read back from a file holds the
                # other positions but not how sure each of them was
                candidates = [
                    {"shift": shift, "intensity": None, "confidence": None}
                    for shift in bore_candidates.read_alternatives(
                        self.find_note(peaklist, "alternatives", index)
                    )
                ]

            if len(candidates) == 0:
                continue

            uncertain.append(
                {
                    "name": name,
                    "note": note,
                    "index": index,
                    "shift3": dictionary["shift3"][index],
                    "candidates": candidates,
                }
            )

        return uncertain

    def offer_to_resolve(self):
        """
        Say how many peaks have an uncertain position down the bore dimension
        and offer to look at them now. They can also be looked at later with
        the resolve button.
        """
        uncertain = self.find_uncertain_peaks()
        if len(uncertain) == 0:
            return

        dlg = wx.MessageDialog(
            self,
            "{} of the peaks have a position down the bore dimension which is "
            "not certain, either because a close neighbour has a similar claim "
            "on a maximum or because fewer maxima were found than expected. "
            "Would you like to look at them now? They can also be looked at "
            "later using the resolve bore peaks button.".format(len(uncertain)),
            "Resolve peaks",
            wx.YES_NO | wx.ICON_QUESTION,
        )
        result = dlg.ShowModal()
        dlg.Destroy()

        if result == wx.ID_YES:
            self.OnResolveAmbiguities(wx.EVT_BUTTON)

    def OnResolveAmbiguities(self, event):
        """
        Show the peaks whose position down the bore is not certain so that each
        of them can be looked at in the spectrum and given a position, kept as
        it is, or accepted where it now sits.
        """
        uncertain = self.find_uncertain_peaks()

        if len(uncertain) == 0:
            dlg = wx.MessageDialog(
                self,
                "There are no peaks with an uncertain position down the bore "
                "dimension in this peaklist.",
                "Resolve peaks",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return

        # Only one of these windows at a time, so that it always shows the
        # peaklist as it is now
        self.close_resolve_dialog()

        self.resolve_dialog = ResolveAmbiguityDialog(self, uncertain)
        self.resolve_dialog.Show()
        self.resolve_dialog.Raise()

    def close_resolve_dialog(self):
        """
        Put away the window for resolving peaks, if it is open.
        """
        dialog = getattr(self, "resolve_dialog", None)
        if dialog == None:
            return

        try:
            dialog.Destroy()
        except RuntimeError:
            pass

        self.resolve_dialog = None

    def reopen_resolve_dialog(self):
        """
        Show the peaks which are still uncertain, after some of them have been
        looked at or moved in the spectrum.
        """
        self.close_resolve_dialog()
        self.OnResolveAmbiguities(wx.EVT_BUTTON)

    def finish_resolving(self):
        """
        Stop showing the positions a peak could have, now that the window for
        resolving peaks has gone.
        """
        self.resolve_dialog = None
        self.main_frame.candidate_shifts = []
        self.main_frame.OnBoreSlider(wx.EVT_BUTTON)

    def find_peak_for_resolving(self, peak):
        """
        A peak of the resolving window as it is in the peaklist now, which is
        not always where it was when the window was opened.
        """
        peaklist = self.find_current_peaklist()
        if peaklist == None:
            return peak

        dictionary = self.peak_list_dictionary[peaklist]
        index = peak.get("index")

        held = dict(peak)
        try:
            if dictionary["peak_name"][index] != peak["name"]:
                # The peaklist has changed since the window was opened
                index = dictionary["peak_name"].index(peak["name"])
                held["index"] = index

            held["shift3"] = dictionary["shift3"][index]
        except (KeyError, IndexError, TypeError, ValueError):
            pass

        return held

    def show_peak_for_resolving(self, peak):
        """
        Show one peak of the resolving window in the spectrum: select it, so
        that the select and move tools work on it, move the position marker onto
        it, zoom the plane onto it and show the positions it could have as lines
        down the bore.
        """
        peaklist = self.find_current_peaklist()
        if peaklist == None:
            return

        held = self.find_peak_for_resolving(peak)
        index = held.get("index")

        dictionary = self.peak_list_dictionary[peaklist]
        try:
            shift1 = dictionary["shift1"][index]
            shift2 = dictionary["shift2"][index]
        except (KeyError, IndexError, TypeError):
            return

        # The peak becomes the selected one, so that the peak tools work on it
        self.selected_peaklist = peaklist
        self.selected_peak_indexes = [index]
        self.selected_peakname = held["name"]
        self.main_frame.selected_bore_peaks = [index]

        # The positions it could have are drawn down the bore alongside its own
        self.main_frame.candidate_shifts = [
            candidate["shift"] for candidate in peak.get("candidates", [])
        ]

        # Zoom the plane onto the peak and show the bore dimension there
        width = 0.1
        self.zoom_to_region(
            self.main_frame.ax_bore,
            [shift1 - width, shift1 + width],
            [shift2 - width, shift2 + width],
        )

        show_bore_position = getattr(self.main_frame, "show_bore_position", None)
        if show_bore_position != None:
            show_bore_position(shift1, shift2)
        else:
            self.main_frame.OnBoreSlider(wx.EVT_BUTTON)

        self.AddToTable()

    def apply_resolutions(self, uncertain, answers) -> int:
        """
        Give the peaks the positions down the bore which were chosen for them.
        A peak which was left alone is not touched, and one which was accepted
        where it sits keeps its position and is noted as resolved.
        """
        peaklist = self.find_current_peaklist()
        if peaklist == None:
            return 0

        dictionary = self.peak_list_dictionary[peaklist]
        self.pad_notes(peaklist)
        changed = 0

        for i, answer in enumerate(answers):
            if answer == None:
                continue

            held = self.find_peak_for_resolving(uncertain[i])
            index = held.get("index")

            try:
                previous = float(dictionary["shift3"][index])
            except (KeyError, IndexError, TypeError, ValueError):
                continue

            if answer == "current":
                # The peak stays where it is, which the user has accepted
                chosen = previous
            else:
                chosen = float(answer["shift"])

            # The position which was given up becomes one of the other options
            others = bore_candidates.read_alternatives(
                self.find_note(peaklist, "alternatives", index)
            )
            others = [shift for shift in others if abs(shift - chosen) > 1e-6]
            if previous != 0 and abs(previous - chosen) > 1e-6:
                others.append(previous)

            dictionary["shift3"][index] = chosen
            dictionary["intensity"][index] = self.find_intensity_3d(
                dictionary["shift1"][index], dictionary["shift2"][index], chosen
            )
            dictionary["ambiguity"][index] = bore_candidates.RESOLVED.capitalize()
            dictionary["alternatives"][index] = ";".join(
                ["{:.5f}".format(shift) for shift in others]
            )
            changed += 1

        self.AddToTable()
        self.main_frame.OnBoreSlider(wx.EVT_BUTTON)

        return changed

    def OnReanalyseBore(self, event):
        """
        Work out again which maxima down the bore belong to which peak, using
        the expected number and the distance which are set now. The positions
        of the peaks in the plane are not changed.
        """
        peaklist = self.find_current_peaklist()
        if peaklist == None:
            self.no_peaklist_message("Analysing the bore dimension")
            return

        dictionary = self.peak_list_dictionary[peaklist]

        # Peaks which were picked from the same position in the plane belong to
        # one peak of the reference plane, so they are gathered back together
        positions = []
        groups = []
        for index in range(len(dictionary["peak_name"])):
            position = (dictionary["shift1"][index], dictionary["shift2"][index])
            for i, held in enumerate(positions):
                if abs(held[0] - position[0]) < 1e-6 and abs(held[1] - position[1]) < 1e-6:
                    groups[i].append(index)
                    break
            else:
                positions.append(position)
                groups.append([index])

        assignments = self.analyse_bore_positions(positions)
        if assignments == None:
            return

        self.apply_bore_assignments(peaklist, positions, groups, assignments)

    def OnExperimentSelection(self, event):
        """
        The user has chosen the experiment, which says how many peaks each peak
        of the reference plane is expected to have down the bore dimension.
        Choosing no experiment leaves the number alone, so that it can be typed
        in for an experiment which is not on the list.
        """
        chosen = self.experiment_box.GetValue()

        for name, expected in BORE_EXPERIMENTS:
            if name != chosen:
                continue

            if expected == None:
                return

            self.expected_peaks_box.SetValue(str(expected))
            return

    def find_expected_peaks(self) -> int:
        """
        How many maxima down the bore dimension are expected for each peak of
        the reference plane, which depends on the experiment. Zero means that
        every maximum which is found is kept.
        """
        try:
            expected = int(float(self.expected_peaks_box.GetValue()))
        except (ValueError, AttributeError, TypeError):
            return 0

        if expected < 0:
            return 0

        return expected

    def find_closeness(self) -> list:
        """
        How close two peaks have to be in each dimension of the plane before a
        maximum down the bore of one of them could belong to the other. The
        boxes are filled in from the spacing of the data and can be changed.
        """
        viewer = self.find_viewer()
        automatic = bore_candidates.find_default_closeness(
            viewer.ppms_0, viewer.ppms_1
        )

        distances = []
        for dimension, box in enumerate([self.closeness_box_1, self.closeness_box_2]):
            try:
                distances.append(abs(float(box.GetValue())))
            except (ValueError, AttributeError, TypeError):
                distances.append(automatic[dimension])

        return distances

    def update_closeness_boxes(self):
        """
        Fill in the distances at which peaks can be confused, worked out from
        the spacing of the data, leaving anything the user has typed alone.
        """
        viewer = self.find_viewer()
        try:
            automatic = bore_candidates.find_default_closeness(
                viewer.ppms_0, viewer.ppms_1
            )
        except (AttributeError, TypeError, IndexError):
            return

        for dimension, box in enumerate([self.closeness_box_1, self.closeness_box_2]):
            try:
                if box.GetValue().strip() == "":
                    box.SetValue("{:.4f}".format(automatic[dimension]))
            except (RuntimeError, AttributeError):
                pass

    def find_template_peaklist(self, picked_peak_array):
        """
        Relate the peaks picked in the whole 3D back to the reference peaklist by
        the two dimensions they share, and make a peaklist of the answer.

        Each reference peak takes its name, and the peaks which sit on it in the
        plane. A peak which more than one reference peak could claim says so and
        names the others. The peaks which belong to nothing in the reference
        plane are kept at the end of the peaklist, so that they can be looked at
        rather than lost.
        """
        template = []
        reference = self.reference_peaklist

        # The reference plane may be held the other way round from the plots
        ppms0 = reference["shift1"]
        ppms1 = reference["shift2"]
        swapped = self.check_reference_peaklists(
            self.main_frame.main_frame.ppms_0,
            self.main_frame.main_frame.ppms_1,
            ppms0,
            ppms1,
        )

        if swapped == None:
            message = (
                "The majority of reference plane peaks are not located within "
                "the current 2D plane of the spin bore. Try loading a different "
                "reference plane peaklist or peak pick in all 3 dimensions. The "
                "reference plane peaklist button will be turned off."
            )
            dlg = wx.MessageDialog(None, message, "Pick Peaks", wx.OK)
            dlg.ShowModal()
            dlg.Destroy()
            self.remove_reference_plane()
            return None

        if swapped == True:
            ppms0, ppms1 = ppms1, ppms0

        for i, name in enumerate(reference["peak_name"]):
            template.append(
                {"name": name, "shift1": ppms0[i], "shift2": ppms1[i]}
            )

        picked = [
            {
                "shift1": peak[0],
                "shift2": peak[1],
                "shift3": peak[2],
                "intensity": peak[3],
            }
            for peak in picked_peak_array
        ]

        assignments, unmatched = bore_candidates.find_template_assignments(
            picked, template, self.find_closeness(), self.find_expected_peaks()
        )

        dictionary = {
            "peak_name": [],
            "shift1": [],
            "shift2": [],
            "shift3": [],
            "intensity": [],
            "ambiguity": [],
            "alternatives": [],
        }

        self.peak_candidates = {}

        for i, assignment in enumerate(assignments):
            name = template[i]["name"]
            others = bore_candidates.find_alternatives_text(assignment["others"])
            notes = [
                note for note in assignment["notes"]
                if note != bore_candidates.AMBIGUOUS
            ]

            if len(assignment["kept"]) == 0:
                # The reference peak has nothing on it, which is worth keeping
                # so that it can be looked at
                dictionary["peak_name"].append(name)
                dictionary["shift1"].append(template[i]["shift1"])
                dictionary["shift2"].append(template[i]["shift2"])
                dictionary["shift3"].append(0)
                dictionary["intensity"].append(0)
                dictionary["ambiguity"].append(
                    bore_candidates.find_ambiguity_text(notes)
                )
                dictionary["alternatives"].append(others)
                self.peak_candidates[name] = assignment["candidates"] if (
                    "candidates" in assignment
                ) else assignment["others"]
                continue

            for j, held in enumerate(assignment["kept"]):
                peakname = name + "_" + str(j + 1)
                dictionary["peak_name"].append(peakname)
                dictionary["shift1"].append(held["shift1"])
                dictionary["shift2"].append(held["shift2"])
                dictionary["shift3"].append(held["shift"])
                dictionary["intensity"].append(held["intensity"])

                peak_notes = list(notes)
                if held["ambiguous"] == True:
                    peak_notes = [
                        "ambiguous with " + ", ".join(held["shared"])
                    ] + peak_notes

                dictionary["ambiguity"].append(
                    bore_candidates.find_ambiguity_text(peak_notes)
                )
                dictionary["alternatives"].append(others)
                self.peak_candidates[peakname] = (
                    assignment["kept"] + assignment["others"]
                )

        # The peaks which belong to nothing in the reference plane
        for count, index in enumerate(unmatched):
            peakname = "unassigned_" + str(count + 1)
            dictionary["peak_name"].append(peakname)
            dictionary["shift1"].append(picked[index]["shift1"])
            dictionary["shift2"].append(picked[index]["shift2"])
            dictionary["shift3"].append(picked[index]["shift3"])
            dictionary["intensity"].append(picked[index]["intensity"])
            dictionary["ambiguity"].append("Not in the reference plane")
            dictionary["alternatives"].append("")

        return dictionary

    def find_projection_intensity(self, shift1, shift2):
        """
        How strong a peak is in the projection which is shown in the plane plot.
        The projection holds the largest intensity down the bore at each
        position, so it says how strong each peak should be down the bore.
        """
        try:
            data, x_values, y_values = self.find_projection_data()
            row = int(np.argmin(np.abs(x_values - shift1)))
            column = int(np.argmin(np.abs(y_values - shift2)))
            return abs(float(data[row][column]))
        except (AttributeError, IndexError, TypeError, ValueError):
            return None

    def find_trace_candidates(self, trace, shifts, found, threshold):
        """
        The maxima down a bore trace: the ones the peak picker reported, and any
        others the trace holds above half of the threshold. A bore can hold more
        than one resonance, and a maximum which did not clear the threshold can
        still be the one a peak belongs to, so it is worth offering.
        """
        candidates = []

        # What the peak picker reported, which it gives back as a table of its
        # own rather than as a list
        reported = []
        try:
            if len(found) > 0:
                reported = found["X_AXIS"]
        except TypeError:
            reported = []

        for j, bore_index in enumerate(reported):
            index = int(round(float(bore_index)))
            try:
                shift = float(shifts[index])
            except (IndexError, TypeError, ValueError):
                continue

            candidates.append(
                {
                    "bore_index": index,
                    "shift": shift,
                    "intensity": float(found[j][-1]),
                }
            )

        candidates += bore_candidates.find_extra_candidates(
            trace,
            shifts,
            abs(threshold) / 2,
            [candidate["bore_index"] for candidate in candidates],
        )

        return candidates

    def find_bore_assignments(self, data, trace_peaks, positions):
        """
        Share the maxima found down the bore between the peaks of the reference
        plane, so that a maximum which belongs to a close neighbour is not
        given to the wrong peak.
        """
        neighbours = bore_candidates.find_close_peaks(
            positions, self.find_closeness()
        )

        return bore_candidates.find_assignments(
            data, trace_peaks, neighbours, self.find_expected_peaks()
        )

    def find_viewer(self):
        """
        The 3D window the bore plot belongs to, which holds the 3D data.
        """
        return self.main_frame.main_frame

    def find_plane_axes_swapped(self) -> bool:
        """
        Whether the axis shown across the bore plot is the second of the two
        plane dimensions of the 3D data rather than the first. The 3D data is
        held as [bore][first plane axis][second plane axis].
        """
        try:
            return len(self.main_frame.new_x_ppms) != len(self.find_viewer().ppms_0)
        except (AttributeError, TypeError):
            return False

    def find_peak_indexes(self, shift1, shift2, shift3):
        """
        Where a peak falls in the 3D data, as indexes into the bore dimension
        and the two plane dimensions, in the order the data is held.
        """
        viewer = self.find_viewer()

        first, second = shift1, shift2
        if self.find_plane_axes_swapped() == True:
            first, second = shift2, shift1

        return (
            int(np.argmin(np.abs(viewer.ppms_2 - shift3))),
            int(np.argmin(np.abs(viewer.ppms_0 - first))),
            int(np.argmin(np.abs(viewer.ppms_1 - second))),
        )

    def find_intensity_3d(self, shift1, shift2, shift3):
        """
        The intensity of the 3D data at the position of a peak.
        """
        try:
            index3, index1, index2 = self.find_peak_indexes(shift1, shift2, shift3)
            return self.find_viewer().nmrdata.data[index3][index1][index2]
        except (AttributeError, IndexError, TypeError, ValueError):
            return 0

    def find_bore_trace(self, shift1, shift2):
        """
        The 1D trace down the bore dimension at the position of a peak, along
        with the chemical shifts of the bore dimension.
        """
        viewer = self.find_viewer()
        index3, index1, index2 = self.find_peak_indexes(shift1, shift2, 0)

        return viewer.ppms_2, viewer.nmrdata.data[:, index1, index2]

    def find_local_maximum_bore(self, shift1, shift2, shift3):
        """
        Walk along the bore dimension from a position to the nearest maximum,
        returning the chemical shift of the maximum and the intensity there.
        The starting position is returned when the data cannot be read.
        """
        try:
            ppms, trace = self.find_bore_trace(shift1, shift2)
        except (AttributeError, IndexError, TypeError, ValueError):
            return shift3, None

        if len(trace) == 0:
            return shift3, None

        index = int(np.argmin(np.abs(ppms - shift3)))

        while True:
            neighbours = [
                position
                for position in [index - 1, index + 1]
                if position >= 0 and position < len(trace)
            ]
            if len(neighbours) == 0:
                break

            best = max(neighbours, key=lambda position: np.abs(trace[position]))
            if np.abs(trace[best]) > np.abs(trace[index]):
                index = best
            else:
                # Neither neighbour is higher, the maximum has been reached
                break

        return ppms[index], trace[index]

    def find_projection_data(self):
        """
        The 2D data shown in the bore plot along with the chemical shifts of
        its two axes, in the order the peaks are held (across the plot, then
        up it).
        """
        bore = self.main_frame

        return bore.nmrdata.data, bore.new_x_ppms, bore.new_y_ppms

    def find_maximum_bore(self, shift1, shift2):
        """
        The chemical shift of the largest point down the bore dimension at the
        position of a peak, used to place a peak which has not been given a
        bore chemical shift yet.
        """
        try:
            ppms, trace = self.find_bore_trace(shift1, shift2)
            if len(trace) == 0:
                return None, None
            index = int(np.argmax(np.abs(trace)))
        except (AttributeError, IndexError, TypeError, ValueError):
            return None, None

        return ppms[index], trace[index]

    def find_local_maximum_plane(self, shift1, shift2, shift3=0):
        """
        Walk uphill from a position to the nearest maximum of the 2D data shown
        in the bore plot, returning the chemical shifts of the maximum and the
        intensity of the 3D data there. The bore chemical shift is used to read
        the intensity and is not changed.
        """
        try:
            data, x_values, y_values = self.find_projection_data()
            rows, columns = data.shape
        except (AttributeError, IndexError, TypeError, ValueError):
            return shift1, shift2, None

        r = int(np.argmin(np.abs(x_values - shift1)))
        c = int(np.argmin(np.abs(y_values - shift2)))

        while True:
            neighbours = [
                (nr, nc)
                for nr in range(r - 1, r + 2)
                for nc in range(c - 1, c + 2)
                if (0 <= nr < rows and 0 <= nc < columns and (nr, nc) != (r, c))
            ]
            if len(neighbours) == 0:
                break

            best = max(neighbours, key=lambda position: np.abs(data[position[0], position[1]]))
            if np.abs(data[best[0], best[1]]) > np.abs(data[r, c]):
                r, c = best
            else:
                # No neighbour is higher, the local maximum has been reached
                break

        shift1, shift2 = x_values[r], y_values[c]

        return shift1, shift2, self.find_intensity_3d(shift1, shift2, shift3)

    def OnSelectPeaks(self, event):
        """
        Drag a box over the plane which is shown to select a group of peaks,
        which can then be moved or moved onto their local maxima together.
        """
        if self.active_move == True or self.active_movez == True:
            if self.active_select_peaks == True:
                self.select_peaks_button.SetValue(True)
            return

        self.selected_peaklist = self.current_peaklist_box.GetValue()

        if self.selected_peaklist not in self.peak_list_dictionary:
            self.select_peaks_button.SetValue(False)
            dlg = wx.MessageDialog(
                None,
                "No peaklists are loaded, please load a peaklist or perform peak picking and try again.",
                "Warning",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return

        if self.active_select_peaks == True:
            # Turning the mode off again
            self.active_select_peaks = False
            self.rect = None
            self.start_point = None
            self.select_peaks_button.SetValue(False)
            self.disconnect_bore("select_press", "select_motion", "select_release")
            self.main_frame.OnBoreSlider(wx.EVT_BUTTON)
            self.update_mode_buttons()
            return

        # The other modes take the same mouse clicks, so they are turned off
        self.turn_off_picking_modes("group")

        self.active_select_peaks = True
        self.select_peaks_button.SetValue(True)

        # The peaks are drawn again so that the selected ones stand out
        self.main_frame.OnBoreSlider(wx.EVT_BUTTON)

        self.select_press = self.main_frame.fig_bore.canvas.mpl_connect(
            "button_press_event", self.on_press_select
        )
        self.select_motion = self.main_frame.fig_bore.canvas.mpl_connect(
            "motion_notify_event", self.on_motion_select
        )
        self.select_release = self.main_frame.fig_bore.canvas.mpl_connect(
            "button_release_event", self.on_release_select
        )

        self.update_mode_buttons()

    def on_press_select(self, event):
        """
        The start of the box which is dragged over a group of peaks.
        """
        if event.inaxes is not self.main_frame.ax_bore:
            return

        self.start_point = (event.xdata, event.ydata)

        self.rect = patches.Rectangle(
            self.start_point, 0, 0, linewidth=1, edgecolor="red", facecolor="none"
        )
        self.main_frame.ax_bore.add_patch(self.rect)
        self.main_frame.canvas_bore.draw_idle()

    def on_motion_select(self, event):
        """
        Showing the box as it is dragged over a group of peaks.
        """
        if self.start_point == None or self.rect == None:
            return

        if event.inaxes is not self.main_frame.ax_bore:
            return

        x0, y0 = self.start_point
        self.rect.set_width(event.xdata - x0)
        self.rect.set_height(event.ydata - y0)
        self.rect.set_xy((x0, y0))
        self.main_frame.canvas_bore.draw_idle()

    def on_release_select(self, event):
        """
        Select every peak of the current peaklist inside the box which was
        dragged. Holding shift adds to the peaks which are already selected.
        """
        if self.start_point == None:
            return

        if event.inaxes is not self.main_frame.ax_bore:
            return

        x0, y0 = self.start_point
        xmin, xmax = sorted([x0, event.xdata])
        ymin, ymax = sorted([y0, event.ydata])

        self.find_peaks_in_area([xmin, xmax], [ymin, ymax], event)

        self.start_point = None
        if self.rect != None:
            self.rect.set_visible(False)
            self.rect = None

        self.main_frame.OnBoreSlider(wx.EVT_BUTTON)
        self.AddToTable()

    def find_peaks_in_area(self, xcoords: list, ycoords: list, event):
        """
        The peaks of the current peaklist which are inside an area of the plane
        which is shown.
        """
        peaklist = self.current_peaklist_box.GetValue()
        if peaklist not in self.peak_list_dictionary:
            return

        if event.key != None and "shift" in str(event.key).lower():
            # Adding to the peaks which are already selected
            selected = [
                index for index in self.selected_peak_indexes if index != "N/A"
            ]
        else:
            selected = []

        dictionary = self.peak_list_dictionary[peaklist]

        for index, peakname in enumerate(dictionary["peak_name"]):
            shift1 = dictionary["shift1"][index]
            shift2 = dictionary["shift2"][index]
            if shift1 > xcoords[0] and shift1 < xcoords[1]:
                if shift2 > ycoords[0] and shift2 < ycoords[1]:
                    if index not in selected:
                        selected.append(index)

        self.selected_peaklist = peaklist
        self.selected_peak_indexes = selected
        if len(selected) == 1:
            self.selected_peakname = dictionary["peak_name"][selected[0]]
        else:
            self.selected_peakname = ""

        # Show the selected peaks down the bore dimension as well
        self.main_frame.selected_bore_peaks = selected

    def update_peaklist_intensities(self, peaklist=None) -> bool:
        """
        Read the intensity of every peak of a peaklist from the 3D data at the
        position of the peak, returning whether the intensities were updated.
        """
        if peaklist == None:
            peaklist = self.find_current_peaklist()

        if peaklist not in self.peak_list_dictionary:
            return False

        dictionary = self.peak_list_dictionary[peaklist]
        intensities = []

        try:
            for index, peakname in enumerate(dictionary["peak_name"]):
                intensities.append(
                    self.find_intensity_3d(
                        dictionary["shift1"][index],
                        dictionary["shift2"][index],
                        dictionary["shift3"][index],
                    )
                )
        except (KeyError, IndexError, AttributeError, TypeError):
            # The spectrum cannot be read, the intensities are left as they were
            return False

        dictionary["intensity"] = intensities

        return True

    def find_peaks_to_move(self):
        """
        The peaks which are currently selected, telling the user when there
        are none. An empty list is returned when nothing is selected.
        """
        if self.selected_peaklist not in self.peak_list_dictionary:
            self.selected_peaklist = self.current_peaklist_box.GetValue()

        indexes = [
            index for index in self.selected_peak_indexes if index != "N/A"
        ]

        if self.selected_peaklist not in self.peak_list_dictionary or len(indexes) == 0:
            dlg = wx.MessageDialog(
                self,
                "There are no peaks selected. Please select a peak using the select peak button and try again.",
                "Warning",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return []

        return indexes

    def OnFindLocalMaximum(self, event):
        """
        Move the selected peaks onto the nearest maximum of the plane which is
        shown in the bore plot, leaving their bore chemical shift alone.
        """
        indexes = self.find_peaks_to_move()
        if len(indexes) == 0:
            return

        dictionary = self.peak_list_dictionary[self.selected_peaklist]

        for index in indexes:
            shift1, shift2, intensity = self.find_local_maximum_plane(
                dictionary["shift1"][index],
                dictionary["shift2"][index],
                dictionary["shift3"][index],
            )
            if intensity == None:
                continue

            dictionary["shift1"][index] = shift1
            dictionary["shift2"][index] = shift2
            dictionary["intensity"][index] = intensity

        self.main_frame.OnBoreSlider(wx.EVT_BUTTON)
        self.AddToTable()

    def OnFindLocalMaximumBore(self, event):
        """
        Move the selected peaks onto the nearest maximum along the bore
        dimension, leaving their position in the plane alone.
        """
        indexes = self.find_peaks_to_move()
        if len(indexes) == 0:
            return

        dictionary = self.peak_list_dictionary[self.selected_peaklist]

        for index in indexes:
            shift3, intensity = self.find_local_maximum_bore(
                dictionary["shift1"][index],
                dictionary["shift2"][index],
                dictionary["shift3"][index],
            )
            if intensity == None:
                continue

            dictionary["shift3"][index] = shift3
            dictionary["intensity"][index] = intensity

        self.main_frame.OnBoreSlider(wx.EVT_BUTTON)
        self.AddToTable()

    def on_click_bore_dimension(self, event):
        """
        When a peak is selected, clicking along the bore dimension moves that
        peak in the bore dimension alone. The peak is put where it was clicked,
        or on the nearest maximum along the bore when that option is chosen.
        """
        if event.button != 1:
            return

        if self.active_movez == True:
            # The move peak (z) mode is dealing with the click
            return

        try:
            bore_axes = [self.main_frame.ax_bore_2, self.main_frame.ax_bore_3]
        except AttributeError:
            return

        if event.inaxes not in bore_axes:
            return

        if self.selected_peaklist not in self.peak_list_dictionary:
            return

        indexes = [index for index in self.selected_peak_indexes if index != "N/A"]
        if len(indexes) == 0:
            return

        shift3 = event.ydata
        if shift3 == None:
            return

        dictionary = self.peak_list_dictionary[self.selected_peaklist]

        for index in indexes:
            shift1 = dictionary["shift1"][index]
            shift2 = dictionary["shift2"][index]

            intensity = None
            if self.snap_bore_box.GetValue() == True:
                # The peak goes onto the maximum nearest to the click
                new_shift3, intensity = self.find_local_maximum_bore(
                    shift1, shift2, shift3
                )
            else:
                new_shift3 = shift3

            if intensity == None:
                new_shift3 = shift3
                intensity = self.find_intensity_3d(shift1, shift2, new_shift3)

            dictionary["shift3"][index] = new_shift3
            dictionary["intensity"][index] = intensity

        self.main_frame.OnBoreSlider(wx.EVT_BUTTON)
        self.AddToTable()

    def find_duplicate_peak(self, shift1, shift2, shift3, peaklist):
        """
        The name of a peak which is already in the peaklist at the position
        given, or None when there is no peak there. Two peaks are at the same
        position when they fall on the same point of the 3D data in all three
        dimensions.
        """
        if peaklist not in self.peak_list_dictionary:
            return None

        dictionary = self.peak_list_dictionary[peaklist]

        try:
            new_index = self.find_peak_indexes(shift1, shift2, shift3)
        except (AttributeError, IndexError, TypeError, ValueError):
            new_index = None

        for i, peakname in enumerate(dictionary["peak_name"]):
            try:
                shifts = [
                    dictionary["shift1"][i],
                    dictionary["shift2"][i],
                    dictionary["shift3"][i],
                ]
            except IndexError:
                continue

            if new_index == None:
                if shifts == [shift1, shift2, shift3]:
                    return peakname
                continue

            if self.find_peak_indexes(*shifts) == new_index:
                return peakname

        return None

    def check_duplicate_peak(self, shift1, shift2, shift3, peaklist) -> bool:
        """
        Warn the user when a new peak would be added on top of a peak which is
        already in the peaklist, returning whether the peak should be added.
        """
        duplicate = self.find_duplicate_peak(shift1, shift2, shift3, peaklist)
        if duplicate == None:
            return True

        dlg = wx.MessageDialog(
            None,
            "The peak {} in the peaklist {} is already at this position in all "
            "three dimensions. Adding this peak will give two peaks on top of "
            "each other. Do you want to add it anyway?".format(duplicate, peaklist),
            "Peaks at the Same Position",
            wx.YES_NO | wx.NO_DEFAULT | wx.ICON_QUESTION,
        )
        result = dlg.ShowModal()
        dlg.Destroy()

        return result == wx.ID_YES

    def find_picked_shifts(self, peaks):
        """
        The chemical shifts of the peaks which nmrglue has picked, in the order
        the peaklist holds them: the axis shown across the bore plot, the axis
        shown up it, and the bore dimension.

        nmrglue names the last axis of the data X_AXIS, then Y_AXIS and
        Z_AXIS. The data is held with the bore dimension first, then the axis
        shown across the plot, then the axis shown up it, so the axis across
        the plot is the Y_AXIS of the picked peaks and the axis up it is the
        X_AXIS. Transposing the bore plot swaps the two unit conversions over
        with the axes they are drawn on.
        """
        if self.main_frame.transposed2D == False:
            x = self.main_frame.uc0.ppm(peaks["Y_AXIS"])
            y = self.main_frame.uc1.ppm(peaks["X_AXIS"])
        else:
            x = self.main_frame.uc0.ppm(peaks["X_AXIS"])
            y = self.main_frame.uc1.ppm(peaks["Y_AXIS"])

        z = self.main_frame.main_frame.uc2.ppm(peaks["Z_AXIS"])

        return x, y, z

    def OnPickPeaks(self, event):
        """
        Pick peaks using nmrglue peak picking routines and then load this peaklist.
        """

        # See if a peaklist is already loaded and then infor: Only one peaklist can be loaded at once for 3D spectra. Continuing will
        # load the new peaklist and the previous peaklist will no longer be loaded

        if(self.current_peaklist_box.GetValue()!=''):
            message = 'Only one peaklist can be loaded at once for 3D spectra. Continuing will load the new peaklist in place of the previous peaklist. Would you like to continue?'
            dlg = wx.MessageDialog(None, message, "Pick Peaks", wx.YES_NO)
            result=dlg.ShowModal()
            if(result == wx.ID_NO):
                dlg.Destroy()
                return
            dlg.Destroy()


        # See if the peaklist name is already in the current directory, and ask the user
        # if they wish to overwrite this.
        peaklist_name = self.peaklist_name_box.GetValue()
        if(peaklist_name in os.listdir()):
            message = 'The peaklist ({}) is already in the current directory. Would you like to overwrite this?'.format(peaklist_name)
            dlg = wx.MessageDialog(None, message, "Pick Peaks", wx.YES_NO)
            result=dlg.ShowModal()
            if(result == wx.ID_NO):
                dlg.Destroy()
                return
            dlg.Destroy()
        
        # Check to see the validity of the value in the threshold box.
        threshold_box_value = self.peak_picking_threshold_box.GetValue()
        try:
            threshold = float(threshold_box_value)
            if(threshold < 0 or threshold > 100):
                message = 'The value in the threshold box ({}) is not a number between 0 and 100, please correct this and try again.'.format(threshold_box_value)
                dlg = wx.MessageBox(message, "Pick Peaks", wx.OK)
                return

        except:
            message = 'The value in the threshold box ({}) is not a number, please correct this and try again.'.format(threshold_box_value)
            dlg = wx.MessageBox(message, "Pick Peaks", wx.OK)
            return
        

        # A reference plane can either be looked along the bore of, or used as a
        # template which the peaks of the whole 3D are related back to
        as_template = (
            self.reference_plane == True
            and self.reference_method_box.GetValue() == TEMPLATE_METHOD
        )

        if self.reference_plane == False or as_template == True:
        
            data = self.main_frame.main_frame.nmrdata.data
            x = self.main_frame.ppms_0
            y = self.main_frame.ppms_1
            z = self.main_frame.ppms_2

            threshold = float(threshold_box_value)/100 *np.max(data)

            algorithm_selection =self.peak_picking_algorithm_box.GetValue()
            sign_option = self.peak_picking_type.GetValue()
            if(algorithm_selection == 'thres' or algorithm_selection == 'thres-fast'):
                if(sign_option=='Positive Peaks'):
                    peaks = ng.peakpick.pick(data, pthres=threshold, algorithm=algorithm_selection, msep=[1,1,1])
                elif(sign_option=='Negative Peaks'):
                    peaks = ng.peakpick.pick(data, nthres=threshold, algorithm=algorithm_selection, msep=[1,1,1])
                else:
                    peaks = ng.peakpick.pick(data, pthres=threshold, nthresh=threshold, algorithm=algorithm_selection, msep=[1,1,1])
            else:
                if(sign_option=='Positive Peaks'):
                    peaks = ng.peakpick.pick(data, pthres=threshold, algorithm=algorithm_selection)
                elif(sign_option=='Negative Peaks'):
                    peaks = ng.peakpick.pick(data, nthres=threshold, algorithm=algorithm_selection)
                else:
                    peaks = ng.peakpick.pick(data, pthres=threshold, nthresh=threshold, algorithm=algorithm_selection)
            
            x, y, z = self.find_picked_shifts(peaks)

            intensities = []
            for p in peaks:
                intensities.append(p[-1])


            picked_peak_array = []
            for i, xval in enumerate(x):
                picked_peak_array.append([xval, y[i], z[i], intensities[i]])

            if as_template == True:
                # The peaks of the whole 3D are related back to the reference
                # peaks by the two dimensions they share
                dictionary = self.find_template_peaklist(picked_peak_array)
                if dictionary == None:
                    return
            else:
                dictionary = {}
                dictionary["peak_name"] = []
                dictionary["shift1"] = []
                dictionary["shift2"] = []
                dictionary["shift3"] = []
                dictionary["intensity"] = []
                dictionary["ambiguity"] = []
                dictionary["alternatives"] = []
                for i, peak in enumerate(picked_peak_array):
                    dictionary["peak_name"].append(str(i+1))
                    dictionary["shift1"].append(peak[0])
                    dictionary["shift2"].append(peak[1])
                    dictionary["shift3"].append(peak[2])
                    dictionary["intensity"].append(peak[3])
                    dictionary["ambiguity"].append("")
                    dictionary["alternatives"].append("")

        else:
            # Picking 3D peaks using a 2D reference plane
            # Effectively peak picking down the bore dimension
            data = self.main_frame.main_frame.nmrdata.data
            x = self.main_frame.main_frame.ppms_0
            y = self.main_frame.main_frame.ppms_1
            z = self.main_frame.main_frame.ppms_2
            
            names = self.reference_peaklist['peak_name']
            ppms0 = self.reference_peaklist['shift1']
            ppms1 = self.reference_peaklist['shift2']

            # Check that ppms0 are in x and ppms1 are in y, otherwise
            # ppms0 and ppms1 might need to be swapped

            check_reference_plane = self.check_reference_peaklists(x, y, ppms0, ppms1)


            if(check_reference_plane == None):
                # The reference plane does not fit the data. Un-toggling the loading reference frame
                # button and informing the user to try a different reference peaklist or to peak pick
                # in all 3 dimensions instead.
                message = 'The majority of reference plane peaks are not located within the current 2D plane of the spin bore. Try loading a different reference plane peaklist or peak pick in all 3 dimensions. The reference plane peaklist button will be turned off.'
                dlg = wx.MessageDialog(None, message, "Pick Peaks", wx.OK)
                result=dlg.ShowModal()
                dlg.Destroy()
                self.remove_reference_plane()
            elif(check_reference_plane == True):
                ppms0_old = ppms0
                ppms1_old = ppms1
                ppms0 = ppms1_old
                ppms1 = ppms0_old


            # For each of the reference plane peaks, the index of the chemical
            # shifts closest to this value needs to be known so a list of 1D
            # data points can be individually fitted by the nmrglue peak picker

            # The 3D data is held as [bore][first plane axis][second plane
            # axis], so each reference peak is found on the two plane axes and
            # the trace down the bore is taken there. The axes are found by
            # their chemical shifts rather than by their size, so that a plane
            # whose two dimensions have the same number of points is not
            # turned round
            planes = []
            for i, ppm0 in enumerate(ppms0):
                planes.append(
                    (
                        int(np.argmin(np.abs(x - ppm0))),
                        int(np.argmin(np.abs(y - ppms1[i]))),
                    )
                )

            data_1D_slices = np.array(
                [data[:, plane[0], plane[1]] for plane in planes]
            )


            threshold = float(threshold_box_value)/100 *np.max(data)


            algorithm_selection =self.peak_picking_algorithm_box.GetValue()
            sign_option = self.peak_picking_type.GetValue()
            if(algorithm_selection == 'thres' or algorithm_selection == 'thres-fast'):
                peaks = []
                for k, data_slice in enumerate(data_1D_slices):
                    try:
                        if(sign_option=='Positive Peaks'):
                            peaks.append(ng.peakpick.pick(data_slice, pthres=threshold, algorithm=algorithm_selection, msep=1))
                        elif(sign_option=='Negative Peaks'):
                            peaks.append(ng.peakpick.pick(data_slice, nthres=threshold, algorithm=algorithm_selection, msep=1))
                        else:
                            peaks.append(ng.peakpick.pick(data_slice, pthres=threshold, nthresh=threshold, algorithm=algorithm_selection, msep=1))

                    except:
                        peaks.append(0)
            else:
                peaks = []
                for k, data_slice in enumerate(data_1D_slices):
                    try:
                        if(sign_option=='Positive Peaks'):
                            peaks.append(ng.peakpick.pick(data_slice, pthres=threshold, algorithm=algorithm_selection))
                        elif(sign_option=='Negative Peaks'):
                            peaks.append(ng.peakpick.pick(data_slice, nthres=threshold, algorithm=algorithm_selection))
                        else:
                            peaks.append(ng.peakpick.pick(data_slice, pthres=threshold, nthresh=threshold, algorithm=algorithm_selection))
                    except:
                        peaks.append(0)
            
            # The maxima found down the bore of each reference peak, which are
            # then shared out between peaks which are close to one another in
            # the plane
            bore_shifts = np.array(
                self.main_frame.main_frame.uc2.ppm(np.arange(data.shape[0]))
            )

            trace_peaks = []
            for i, peak in enumerate(peaks):
                found = peak
                try:
                    if peak == 0:
                        # The peak picker found nothing down the bore here
                        found = None
                except (ValueError, TypeError):
                    pass

                candidates = self.find_trace_candidates(
                    data_1D_slices[i], bore_shifts, found, threshold
                )

                trace_peaks.append(
                    {
                        "plane": planes[i],
                        "candidates": candidates,
                        "projection": self.find_projection_intensity(
                            ppms0[i], ppms1[i]
                        ),
                    }
                )

            assignments = self.find_bore_assignments(
                data, trace_peaks, list(zip(ppms0, ppms1))
            )

            names1 = []
            x = []
            y = []
            z = []
            intensities = []
            ambiguities = []
            alternatives = []
            self.peak_candidates = {}

            for i, assignment in enumerate(assignments):
                if self.main_frame.transposed2D == False:
                    xval = ppms0[i]
                    yval = ppms1[i]
                else:
                    xval = ppms1[i]
                    yval = ppms0[i]

                kept = assignment["kept"]
                notes = bore_candidates.find_ambiguity_text(assignment["notes"])
                others = bore_candidates.find_alternatives_text(assignment["others"])

                if len(kept) == 0:
                    # Nothing was found down the bore of this peak, which is
                    # worth keeping so that it can be looked at
                    peakname = names[i]
                    names1.append(peakname)
                    x.append(xval)
                    y.append(yval)
                    z.append(0)
                    intensities.append(0)
                    ambiguities.append(notes)
                    alternatives.append(others)
                    self.peak_candidates[peakname] = assignment["candidates"]
                    continue

                for j, candidate in enumerate(kept):
                    peakname = names[i] + '_' + str(j + 1)
                    names1.append(peakname)
                    x.append(xval)
                    y.append(yval)
                    z.append(candidate["shift"])
                    intensities.append(candidate["intensity"])
                    if candidate["ambiguous"] == True:
                        ambiguities.append(
                            bore_candidates.find_ambiguity_text(
                                [bore_candidates.AMBIGUOUS]
                                + [
                                    note
                                    for note in assignment["notes"]
                                    if note != bore_candidates.AMBIGUOUS
                                ]
                            )
                        )
                    else:
                        ambiguities.append(
                            bore_candidates.find_ambiguity_text(
                                [
                                    note
                                    for note in assignment["notes"]
                                    if note != bore_candidates.AMBIGUOUS
                                ]
                            )
                        )
                    alternatives.append(others)
                    self.peak_candidates[peakname] = assignment["candidates"]


            


            picked_peak_array = []
            for i, xval in enumerate(x):
                picked_peak_array.append([xval, y[i], z[i], intensities[i]])

            dictionary = {}
            dictionary["peak_name"] = []
            dictionary["shift1"] = []
            dictionary["shift2"] = []
            dictionary["shift3"] = []
            dictionary["intensity"] = []
            dictionary["ambiguity"] = []
            dictionary["alternatives"] = []
            for i, peak in enumerate(picked_peak_array):
                dictionary["ambiguity"].append(ambiguities[i])
                dictionary["alternatives"].append(alternatives[i])
                if(self.reference_plane==True):
                    # Keeping the naming consistent with the reference plane
                    dictionary["peak_name"].append(names1[i])
                else:
                    # Creating a new name for every peak
                    dictionary["peak_name"].append(str(i+1))
                dictionary["shift1"].append(peak[0])
                dictionary["shift2"].append(peak[1])
                dictionary["shift3"].append(peak[2])
                dictionary["intensity"].append(peak[3])
        
        
        # Create a file to store the 3D peaklist
        with open(peaklist_name, 'w') as file:
            file.write('')
        
        p = pathlib.Path(peaklist_name)
        dirs = p.parts[-3:]
        file_name = p.parts[-1]
        last_directories_path = str(pathlib.Path(*dirs))
        peaklist = dictionary
        self.peak_list_dictionary[last_directories_path] = peaklist
        if self.peak_list_choices == [""]:
            self.peak_list_choices = [last_directories_path]
        else:
            self.peak_list_choices.append(last_directories_path)

        self.peak_list = self.peak_list_choices[-1]

        self.turn_off_togglebuttons()

        self.AddToTable()

        # The box shows the peaklist by the same name as it is held under, which
        # is the last few parts of the path rather than the whole of it
        self.current_peaklist_box.SetValue(last_directories_path)

        # Where the peaklist really is, so that saving goes back to it
        self.peaklist_path = str(pathlib.Path(peaklist_name).absolute())

        # and so that it is loaded again when this window is opened next time
        self.remember_peaklist(peaklist_name)

        # Update the plot with the new peaklist
        self.main_frame.OnBoreSlider(wx.EVT_BUTTON)

        # Save the 3D peaklist
        self.OnSave(wx.EVT_BUTTON, peaklist_file=peaklist_name)

        # Offer to look at any peak whose position down the bore is not certain
        self.offer_to_resolve()

    def AddPeakListBrowser(self, event):
        """
        1 - Open a file explorer window (opening at the current directory)
        2 - Try to read the peaklist file (might be necessary to transpose)
        3 - Plot the peaklist file (and when open mincontour 2D need to also
            plot the peaklists too)
        """
        # Opening up a file window asking the user to select the 1D peak list - must be in the format of 1st column = peak_name, 2nd column = peak_position
        dlg = wx.FileDialog(self, "Select the peak list", wildcard="", style=wx.FD_OPEN)
        dlg.SetDirectory(os.getcwd())
        if dlg.ShowModal() == wx.ID_OK:
            peaklist_file = dlg.GetPath()
        else:
            dlg.Destroy()
            return

        self.AddPeaklist(peaklist_file)

    def AddPeaklist(self, peaklist_file, new_peaklist=False):
        p = pathlib.Path(peaklist_file)
        dirs = p.parts[-3:]
        file_name = p.parts[-1]
        last_directories_path = str(pathlib.Path(*dirs))
        if ".xlsx" in file_name:
            peaklist = self.ReadCCPNList(peaklist_file)
        else:
            peaklist = self.ReadPeakList(peaklist_file, new_peaklist)
        if type(peaklist) != dict:
            return
        self.peak_list_dictionary[last_directories_path] = peaklist
        if self.peak_list_choices == [""]:
            self.peak_list_choices = [last_directories_path]
        else:
            self.peak_list_choices.append(last_directories_path)

        self.peak_list = self.peak_list_choices[-1]
        self.turn_off_togglebuttons()

        # The peaks which were selected belonged to the peaklist which was
        # loaded before this one
        self.selected_peak_indexes = ["N/A"]
        self.selected_peakname = ""
        self.selected_peaklist = last_directories_path
        self.main_frame.selected_bore_peaks = []

        self.current_peaklist_box.SetValue(last_directories_path)

        # Where the peaklist really is, so that saving goes back to it
        self.peaklist_path = str(pathlib.Path(peaklist_file).absolute())

        # and so that it is loaded again when this window is opened next time
        self.remember_peaklist(peaklist_file)

        # The peaklist matches the file it has just been read from
        self.mark_saved()

        self.AddToTable()

        self.main_frame.OnBoreSlider(wx.EVT_BUTTON)

    def AddToTable(self):
        """
        Adding the peaklist just entered into the peaklist table
        """
        # Every change to the peaklist ends up here, so the save button is
        # marked or unmarked and the columns are named from this one place
        self.update_saved_button()
        self.update_column_labels()

        peaklist = self.find_current_peaklist()
        if peaklist != None:
            self.pad_notes(peaklist)

        row_count = self.grid.GetNumberRows()
        if row_count > 0:
            self.grid.DeleteRows(0, row_count)

        data = self.find_peak_rows()

        num_rows = self.grid.GetNumberRows()
        if len(data) > num_rows:
            self.grid.AppendRows(len(data) - num_rows)
        for row, rowData in enumerate(data):
            for col, value in enumerate(rowData):
                self.grid.SetCellValue(row, col, str(value))

    def on_begin_edit(self, event):
        """
        If the user is editing the peak_name column, store the original value
        """
        row = event.GetRow()
        col = event.GetCol()
        if col == 0:
            self.old_key = self.grid.GetCellValue(row, col)
        else:
            self.old_num = self.grid.GetCellValue(row, col)
        event.Skip()

    def on_cell_changed(self, event):
        """
        When a cell is changed, see if the types are correct
        e.g. the shifts are numbers.
        Can then update the dictionary and re-perform OnMinContour2D.
        """
        row = event.GetRow()
        col = event.GetCol()
        if self.old_key != None:
            peak_name = self.grid.GetCellValue(row, col)
            if (
                peak_name
                in self.peak_list_dictionary[self.current_peaklist_box.GetValue()][
                    "peak_name"
                ]
            ):
                # Give an error saying that this peak name is already taken, changing back to the original value
                self.grid.SetCellValue(row, col, self.old_key)
                dlg = wx.MessageDialog(
                    self,
                    "The peak name entered (row:{}, coloum:{})is already taken, this value has been reset to its previous value".format(
                        str(row), str(col)
                    ),
                    "Warning",
                    wx.OK,
                )
                dlg.ShowModal()
                dlg.Destroy()
            else:
                index = self.peak_list_dictionary[self.current_peaklist_box.GetValue()][
                    "peak_name"
                ].index(self.old_key)
                self.peak_list_dictionary[self.current_peaklist_box.GetValue()][
                    "peak_name"
                ][index] = peak_name

        else:
            peak_name = self.grid.GetCellValue(row, 0)
            index = self.peak_list_dictionary[self.current_peaklist_box.GetValue()][
                "peak_name"
            ].index(peak_name)
            try:
                new_value = float(self.grid.GetCellValue(row, col))
                if col == 1:
                    self.peak_list_dictionary[self.current_peaklist_box.GetValue()][
                        "shift1"
                    ][index] = new_value
                if col == 2:
                    self.peak_list_dictionary[self.current_peaklist_box.GetValue()][
                        "shift2"
                    ][index] = new_value
            except:
                dlg = wx.MessageDialog(
                    self,
                    "The value entered (row:{}, coloum:{})is not a number, this value has been reset to its previous value".format(
                        str(row), str(col)
                    ),
                    "Warning",
                    wx.OK,
                )
                dlg.ShowModal()
                dlg.Destroy()

        self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)

        self.old_key = None
        self.old_num = None

        self.AddToTable()

    def ReadPeakList(self, peaklist_file, new_peaklist = False, reference_plane = False):
        """
        Read the selected peaklist to obtain the chemical shifts in each dimension
        Add a list of peak names, chemical shifts (dim1) and chemical shifts (dim2)
        to the dictionary
        """

        if(reference_plane==False):
            dictionary = {}
            dictionary["peak_name"] = []
            dictionary["shift1"] = []
            dictionary["shift2"] = []
            dictionary["shift3"] = []
            dictionary["intensity"] = []
            dictionary["ambiguity"] = []
            dictionary["alternatives"] = []
            name1 = ''
            name2 = ''
            name3 = ''
            # Try to read the peaklist, otherwise give an error saying it could not be read correctly
            try:
                with open(peaklist_file) as file:
                    lines = file.readlines()
                    if len(lines) != 0:
                        for i, line in enumerate(lines):
                            line = line.split("\n")[0].split()
                            if len(line) >= 3:
                                if(i==0):
                                    try:
                                        float(line[1])
                                        dictionary["peak_name"].append(line[0])
                                        dictionary["shift1"].append(float(line[1]))
                                        dictionary["shift2"].append(float(line[2]))
                                        dictionary['shift3'].append(float(line[3]))
                                        try:
                                            dictionary["intensity"].append(float(line[4]))
                                        except:
                                            dictionary["intensity"].append(float(0.0))
                                        self.read_notes(dictionary, lines[i])
                                    except:
                                        name1 = line[1]
                                        name2 = line[2]
                                        name3 = line[3]
                                else:
                                    try:
                                        dictionary["peak_name"].append(line[0])
                                        dictionary['shift1'].append(float(line[1]))
                                        dictionary['shift2'].append(float(line[2]))
                                        dictionary['shift3'].append(float(line[3]))
                                        try:
                                            dictionary["intensity"].append(float(line[4]))
                                        except:
                                            dictionary["intensity"].append(float(0.0))
                                        self.read_notes(dictionary, lines[i])
                                    except:
                                        pass

            except:
                self.peaklist_error_message()
                return None

            if len(dictionary["peak_name"]) == 0:
                if (
                    new_peaklist == False
                    and peaklist_file_is_empty(peaklist_file) == False
                ):
                    # The file holds something which could not be read as peaks
                    self.peaklist_error_message()
                    return None

                # An empty peaklist, ready for peaks to be added to it. There
                # are no chemical shifts to compare with the spectrum
                return dictionary

            # Try to see if the chemical shifts of the peaks are within the 2D spectral range
            if(new_peaklist==False):
                dictionary = self.check_peaklist(dictionary)

        else:
            dictionary = {}
            dictionary["peak_name"] = []
            dictionary["shift1"] = []
            dictionary["shift2"] = []

            name1 = ''
            name2 = ''

            # Try to read the peaklist, otherwise give an error saying it could not be read correctly
            try:
                with open(peaklist_file) as file:
                    lines = file.readlines()
                    if len(lines) != 0:
                        for i, line in enumerate(lines):
                            line = line.split("\n")[0].split()
                            if len(line) >= 3:
                                if(i==0):
                                    try:
                                        float(line[1])
                                        dictionary["peak_name"].append(line[0])
                                        dictionary["shift1"].append(float(line[1]))
                                        dictionary["shift2"].append(float(line[2]))
                                    except:
                                        name1 = line[1]
                                        name2 = line[2]

                                else:
                                    try:
                                        dictionary["peak_name"].append(line[0])
                                        dictionary["shift1"].append(float(line[1]))
                                        dictionary["shift2"].append(float(line[2]))
                                    except:
                                        pass

            except:
                self.peaklist_error_message()
                return None

            if len(dictionary["peak_name"]) == 0 and new_peaklist == False:
                self.peaklist_error_message()
                return None

            # Try to see if the chemical shifts of the peaks are within the 2D spectral range
            dictionary = self.check_reference_peaklist(dictionary)


        return dictionary

    def ReadCCPNList(self, peaklist_file, reference_plane=False):
        """
        Read peaklist that has been exported from a CCPN peaklist table.
        """

        message = 'Reading in a peaklist from an excel file is not currently implemented. Please use tabular (.tab format).'
        dlg = wx.MessageDialog(
                self,
                message,
                "Warning",
                wx.OK,
            )
        dlg.ShowModal()
        dlg.Destroy()
        return None


        # df = pd.read_excel(peaklist_file, dtype=str)

        # try:
        #     if(reference_plane==False):
        #         peak_names = df.iloc[:, 0].tolist()
        #         shift1 = df.iloc[:, 8].to_numpy()
        #         shift2 = df.iloc[:, 9].to_numpy()
        #         shift3 = df.iloc[:, 10].to_numpy()
        #         intensity = df.iloc[:, 14].to_numpy()

        #         shift1_1 = []
        #         shift2_1 = []
        #         shift3_1 = []
        #         intensity_1 = []

        #         for i in range(len(shift1)):
        #             shift1_1.append(float(shift1[i]))
        #             shift2_1.append(float(shift2[i]))
        #             shift3_1.append(float(shift3[i]))
        #             intensity_1.append(float(intensity[i]))

        #         dictionary = {}
        #         dictionary["peak_name"] = peak_names
        #         dictionary["shift1"] = shift1_1
        #         dictionary["shift2"] = shift2_1
        #         dictionary["shift3"] = shift3_1
        #         dictionary["intensity"] = intensity_1

        #         # Try to see if the chemical shifts of the peaks are within the 2D spectral range
        #         dictionary = self.check_peaklist(dictionary)

        #         return dictionary
        #     else:
        #         peak_names = df.iloc[:, 0].tolist()
        #         shift1 = df.iloc[:, 8].to_numpy()
        #         shift2 = df.iloc[:, 9].to_numpy()

        #         shift1_1 = []
        #         shift2_1 = []

        #         for i in range(len(shift1)):
        #             shift1_1.append(float(shift1[i]))
        #             shift2_1.append(float(shift2[i]))

        #         dictionary = {}
        #         dictionary["peak_name"] = peak_names
        #         dictionary["shift1"] = shift1_1
        #         dictionary["shift2"] = shift2_1

        #         # Try to see if the chemical shifts of the peaks are within the 2D spectral range
        #         dictionary = self.check_reference_peaklist(dictionary)

        #         return dictionary

        # except:
        #     self.peaklist_error_message()
        #     return None
        

    def check_reference_peaklists(self, x, y, xpeaks, ypeaks):
        """
        Check to see that the xpeaks are in the x chemical shift range and that
        the ypeaks are in the y chemical shift range. If they are not, apply a 
        transpose of the peaklist through swap_peaks = True. If this still doesn't
        work, return False and then inform the user that the 2D peaklist does not 
        seem to fit to the 2D plane.
        """
        ppms_0 = xpeaks
        ppms_1 = ypeaks

        swap_peaks = False


        match_0 = []
        for ppm in ppms_0:
            if ppm > np.min(x) and ppm < np.max(x):
                match_0.append(1)
            else:
                match_0.append(0)

        mean0 = np.mean(np.array(match_0))

        match_1 = []
        for ppm in ppms_1:
            if ppm > np.min(y) and ppm < np.max(y):
                match_1.append(1)
            else:
                match_1.append(0)

        mean1 = np.mean(np.array(match_1))

        if mean0 == 0 and mean1 == 0:
            # No peaks are within the spectrum, trying transposing
            match_0 = []
            for ppm in ppms_0:
                if ppm > np.min(y) and ppm < np.max(y):
                    match_0.append(1)
                else:
                    match_0.append(0)

            mean0 = np.mean(np.array(match_0))

            match_1 = []
            for ppm in ppms_1:
                if ppm > np.min(x) and ppm < np.max(x):
                    match_1.append(1)
                else:
                    match_1.append(0)

            mean1 = np.mean(np.array(match_1))

            if mean0 > 0.5 and mean1 > 0.5:
                # More than 50 percent of the peaks are within the spectrum
                swap_peaks = True
                return swap_peaks

            else:
                return None
            
        return swap_peaks

    def check_peaklist(self, dictionary: dict):
        """
        Match a peaklist to the spectrum which is being shown.

        A peaklist can have been picked from a different projection of the same
        3D data, so its three columns are not always in the order this window
        shows them. Each column is compared with the range of each axis of the
        plots and the columns are put in the order of what is shown: the axis
        across the plane, the axis up it, and the bore dimension. A peaklist
        which does not fit the spectrum at all is refused.
        """
        columns = [
            copy.deepcopy(dictionary["shift1"]),
            copy.deepcopy(dictionary["shift2"]),
            copy.deepcopy(dictionary["shift3"]),
        ]

        ranges = self.find_axis_ranges()

        order = bore_candidates.find_axis_assignment(
            columns, ranges, score=self.find_peaklist_data_fit(columns)
        )

        if order == None:
            dlg = wx.MessageDialog(
                self,
                "The chemical shifts in this peaklist do not fit the spectrum "
                "which is shown, whichever way round its columns are read. Try "
                "a different peaklist, or one picked from this spectrum.",
                "Warning",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return None

        dictionary["shift1"] = columns[order[0]]
        dictionary["shift2"] = columns[order[1]]
        dictionary["shift3"] = columns[order[2]]

        # The columns are now in the order the plots show them, so the first one
        # is the axis across the plane
        self.bore_xdim = 'shift1'

        return dictionary

    def realign_peaklists(self, bore_changed=False) -> dict:
        """
        Put every loaded peaklist into the order of the plane which is shown.

        This is used when the plane is changed, which happens when another data
        orientation is chosen. The peaks themselves do not move: each of them
        holds a shift for each of the three dimensions, and it is only which
        dimension is shown where that has changed, so the columns are matched to
        the axes again in the same way as a peaklist which is read from a file.

        A peaklist which cannot be matched to the new plane is left as it was
        and named in the answer, rather than being scrambled.
        """

        ranges = self.find_axis_ranges()

        answer = {"moved": [], "kept": [], "refused": []}

        for name, dictionary in self.peak_list_dictionary.items():
            columns = [
                copy.deepcopy(dictionary["shift1"]),
                copy.deepcopy(dictionary["shift2"]),
                copy.deepcopy(dictionary["shift3"]),
            ]

            if len(columns[0]) == 0:
                # An empty peaklist has nothing to put in order
                answer["kept"].append(name)
                continue

            order = bore_candidates.find_axis_assignment(
                columns, ranges, score=self.find_peaklist_data_fit(columns)
            )

            if order == None:
                answer["refused"].append(name)
                continue

            if list(order) == [0, 1, 2]:
                answer["kept"].append(name)
                continue

            dictionary["shift1"] = columns[order[0]]
            dictionary["shift2"] = columns[order[1]]
            dictionary["shift3"] = columns[order[2]]
            answer["moved"].append(name)

        # The first column is the axis across the plane again
        self.bore_xdim = "shift1"

        if bore_changed == True:
            # The other positions which were found down the old bore dimension
            # say nothing about the new one
            self.peak_candidates = {}
            for name, dictionary in self.peak_list_dictionary.items():
                if "alternatives" in dictionary:
                    dictionary["alternatives"] = [
                        "" for note in dictionary["alternatives"]
                    ]

        self.AddToTable()

        return answer

    def find_peaklist_data_fit(self, columns):
        """
        A way of telling how well a peaklist fits the spectrum with its columns
        one way round rather than another, used when two dimensions cover the
        same range of chemical shifts and the ranges cannot tell them apart.

        The peaks of the right arrangement sit on the data, so the intensity of
        the spectrum at the peaks says which way round the columns belong.
        """

        def find_fit(order):
            total = 0.0
            counted = 0

            for i in range(len(columns[0])):
                try:
                    total += abs(
                        float(
                            self.find_intensity_3d(
                                columns[order[0]][i],
                                columns[order[1]][i],
                                columns[order[2]][i],
                            )
                        )
                    )
                    counted += 1
                except (IndexError, KeyError, TypeError, ValueError):
                    continue

            if counted == 0:
                return 0.0

            return total / counted

        return find_fit

    def find_axis_ranges(self) -> list:
        """
        The lowest and highest chemical shift of each axis of the plots, in the
        order a peaklist holds them: the axis across the plane, the axis up it,
        and the bore dimension. The plane axes are taken from the plot as it is
        shown, so a projection which has been transposed is followed.
        """
        bore = self.main_frame
        viewer = self.find_viewer()

        axes = []
        for values in [bore.new_x_ppms, bore.new_y_ppms, viewer.ppms_2]:
            try:
                axes.append((float(np.min(values)), float(np.max(values))))
            except (TypeError, ValueError):
                axes.append((0.0, 0.0))

        return axes

    def check_reference_peaklist(self, dictionary: dict):
        """
        Try to see if the chemical shifts of the peaks are within the 2D spectral range
        """
        ppms_0 = copy.deepcopy(dictionary["shift1"])
        ppms_1 = copy.deepcopy(dictionary["shift2"])


        shifts = [ppms_0, ppms_1]

        ppms_projection = []
        for i in range(2):
            ppms_projection.append(shifts[i])


        ppms_0 = ppms_projection[0]
        ppms_1 = ppms_projection[1]


        mean_0 = np.mean(ppms_0)
        mean_1 = np.mean(ppms_1)

        match_0 = []
        for ppm in ppms_0:
            if ppm > np.min(self.main_frame.ppms_0) and ppm < np.max(
                self.main_frame.ppms_0
            ):
                match_0.append(1)
            else:
                match_0.append(0)

        mean0 = np.mean(np.array(match_0))

        match_1 = []
        for ppm in ppms_1:
            if ppm > np.min(self.main_frame.ppms_1) and ppm < np.max(
                self.main_frame.ppms_1
            ):
                match_1.append(1)
            else:
                match_1.append(0)

        mean1 = np.mean(np.array(match_1))

        if mean0 == 0 and mean1 == 0:
            # No peaks are within the spectrum, trying transposing
            match_0 = []
            for ppm in ppms_0:
                if ppm > np.min(self.main_frame.ppms_1) and ppm < np.max(
                    self.main_frame.ppms_1
                ):
                    match_0.append(1)
                else:
                    match_0.append(0)

            mean0 = np.mean(np.array(match_0))

            match_1 = []
            for ppm in ppms_1:
                if ppm > np.min(self.main_frame.ppms_0) and ppm < np.max(
                    self.main_frame.ppms_0
                ):
                    match_1.append(1)
                else:
                    match_1.append(0)

            mean1 = np.mean(np.array(match_1))

            if mean0 > 0.5 and mean1 > 0.5:
                # More than 50 percent of the peaks are within the spectrum
                dictionary["shift1"] = ppms_1
                dictionary["shift2"] = ppms_0
                self.bore_xdim = 'shift1'
                if self.main_frame.transposed2D == True:
                    dictionary["shift1"] = ppms_0
                    dictionary["shift2"] = ppms_1
                    self.bore_xdim = 'shift2'
                return dictionary

            else:
                return None

        else:
            if self.main_frame.transposed2D == True:
                dictionary["shift1"] = ppms_1
                dictionary["shift2"] = ppms_0
                self.bore_xdim = 'shift1'

        return dictionary

    def peaklist_error_message(self):
        """
        Gives the user an error when the peaklist was not read correctly
        """

        dlg = wx.MessageDialog(
            self,
            "The selected peaklist was not read correctly. Please select another peak list.",
            "Error",
        )
        dlg.ShowModal()

    def OnPeakListSelection(self, event):
        if self.selected_peaklist != "":
            self.selected_peaklist = self.current_peaklist_box.GetValue()

        self.turn_off_togglebuttons()

        self.AddToTable()

    def disconnect_bore(self, *names):
        """
        Disconnect mouse handlers from the bore plots by name, ignoring any
        which were never connected.
        """
        for name in names:
            connection = getattr(self, name, None)
            if connection == None:
                continue
            try:
                self.main_frame.fig_bore.canvas.mpl_disconnect(connection)
            except (AttributeError, TypeError):
                pass
            setattr(self, name, None)

    def turn_off_picking_modes(self, keep=""):
        """
        Adding peaks, selecting a peak and selecting a group of peaks all take
        the mouse clicks on the plane, so only one of them can be on at a time.
        Any of them which is on is turned off, apart from the one named.
        """
        if keep != "add" and self.active_add == True:
            self.active_add = False
            self.add_peaks_button.SetValue(False)
            self.disconnect_bore("add_peak_connect")

        if keep != "select" and self.active_select_peak == True:
            self.active_select_peak = False
            self.select_peak_button.SetValue(False)
            self.selected_peakname = ""
            self.main_frame.plot_cross = True
            self.disconnect_bore("select_peak_connect")

        if keep != "group" and self.active_select_peaks == True:
            self.active_select_peaks = False
            self.select_peaks_button.SetValue(False)
            self.rect = None
            self.start_point = None
            self.disconnect_bore("select_press", "select_motion", "select_release")

        self.update_mode_buttons()

    def turn_off_togglebuttons(self):
        # If any toggle buttons are on, turn them off
        if self.active_add == True:
            self.active_add = False
            self.add_peaks_button.SetValue(False)
            self.disconnect_bore("add_peak_connect")
        if self.active_move:
            self.active_move = False
            self.move_peaks_button.SetValue(False)
            if self.active_select_peak:
                self.disconnect_bore("move_peak_connect")
            if self.active_select_peaks:
                self.disconnect_bore(
                    "move_peak_press", "move_peak_motion", "move_peak_release"
                )
        if self.active_movez == True:
            self.active_movez = False
            self.move_peaks_bore_button.SetValue(False)
            self.disconnect_bore("move_peak_connectz")
        if self.active_select_peak == True:
            self.select_peak_button.SetValue(False)
            self.active_select_peak = False
            self.selected_peakname = ""
            self.disconnect_bore("select_peak_connect")
            self.main_frame.OnBoreSlider(wx.EVT_BUTTON)
        if self.active_select_peaks == True:
            self.select_peaks_button.SetValue(False)
            self.active_select_peaks = False
            self.rect = None
            self.start_point = None
            self.disconnect_bore("select_press", "select_motion", "select_release")
            self.main_frame.OnBoreSlider(wx.EVT_BUTTON)

        self.update_mode_buttons()

    def OnLoadReferencePlane(self, event):
        """
        If a reference plane is already loaded, then it is removed.
        Otherwise, setting reference_plane to true and asking the user to select the reference plane
        peaklist to be loaded from a file dialog.
        """

        if self.reference_plane == True:
            self.remove_reference_plane()
            return

        self.reference_plane = True

        # Opening up a file window asking the user to select the 1D peak list - must be in the format of 1st column = peak_name, 2nd column = peak_position
        dlg = wx.FileDialog(self, "Select the 2D reference plane peak list", wildcard="", style=wx.FD_OPEN)
        dlg.SetDirectory(os.getcwd())
        if dlg.ShowModal() == wx.ID_OK:
            peaklist_file = dlg.GetPath()
        else:
            self.remove_reference_plane()
            dlg.Destroy()
            return
        

        p = pathlib.Path(peaklist_file)
        dirs = p.parts[-3:]
        file_name = p.parts[-1]
        if ".xlsx" in file_name:
            peaklist = self.ReadCCPNList(peaklist_file, reference_plane=True)
        else:
            peaklist = self.ReadPeakList(peaklist_file, reference_plane=True)
        if type(peaklist) != dict:
            self.remove_reference_plane()
            return
        
        self.reference_peaklist = peaklist
        self.reference_peaklist_name = file_name
        self.reference_peaklist_path = str(p)

        self.reference_plane_button.SetValue(True)
        self.update_reference_plane_text()

    def OnRemoveReferencePlane(self, event):
        """
        Stop using the reference plane, which was only ever optional. Peaks are
        then picked from the spectrum alone in all three dimensions. The
        peaklist which was loaded as the reference plane is left on disk.
        """

        if self.reference_plane == False:
            return

        self.remove_reference_plane()

    def remove_reference_plane(self):
        """
        Forget the reference plane and say so in the window. This is used both
        when the reference plane is removed on purpose and when loading one did
        not work, so that the button is never left looking as though a reference
        plane is in use when it is not.
        """

        self.reference_plane = False
        self.reference_peaklist = {}
        self.reference_peaklist_name = ""
        self.reference_peaklist_path = ""

        try:
            self.reference_plane_button.SetValue(False)
        except (AttributeError, RuntimeError):
            pass

        self.update_reference_plane_text()

    def update_reference_plane_text(self):
        """
        Say whether a reference plane is being used and, if it is, which
        peaklist it is, as the peaks which come out of picking depend on it.
        The remove button and the choice of how the reference plane is used are
        only of any use while one is loaded.
        """

        try:
            label = self.reference_plane_text
        except AttributeError:
            return

        if self.reference_plane == True:
            name = self.reference_peaklist_name
            if name == "":
                name = "a peaklist which has not been named"

            count = len(self.reference_peaklist.get("peak_name", []))
            peaks = "{} peak{}".format(count, "" if count == 1 else "s")

            label.SetLabel(
                "Reference plane IN USE:  {}  ({})".format(name, peaks)
            )
            label.SetForegroundColour(wx.Colour(0, 100, 0))
            label.SetToolTip(
                "Peaks are picked using this peaklist as the reference plane. "
                "It was read from " + str(self.reference_peaklist_path)
            )
            self.reference_plane_button.SetLabel("Reference plane loaded")
        else:
            label.SetLabel(
                "No reference plane: peaks are picked from the spectrum alone, "
                "in all 3 dimensions"
            )
            label.SetForegroundColour(wx.Colour(100, 100, 100))
            label.SetToolTip(
                "Loading a reference plane is optional. Without one, peak "
                "picking finds the peaks of the whole 3D and names them in the "
                "order they are found."
            )
            self.reference_plane_button.SetLabel("Load reference plane (optional)")

        font = label.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD if self.reference_plane else wx.FONTWEIGHT_NORMAL)
        label.SetFont(font)

        for widget in [self.remove_reference_button, self.reference_method_box]:
            try:
                widget.Enable(self.reference_plane)
            except (AttributeError, RuntimeError):
                pass

        try:
            self.row_pickpeaks.Layout()
            self.Layout()
        except (AttributeError, RuntimeError):
            pass
    
        

    
    

        


    def OnAddPeaks(self, event):
        """
        This will allow a user to add a peak to the currently selected peaklist
        A popout will come up saying that the user needs to use the cursor to
        add a peak. De-select the add button once complete.

        The code will also disable all the other buttons which have been
        selected
        """

        if self.active_add == True:
            self.active_add = False
            self.add_peaks_button.SetValue(False)
            self.disconnect_bore("add_peak_connect")
            self.main_frame.OnBoreSlider(wx.EVT_BUTTON)
            return



        # Only one of the picking modes can be on at a time
        self.turn_off_picking_modes("add")

        if self.peak_list_choices == [""]:
            dlg = wx.MessageDialog(
                None,
                "No peaklists are loaded, would you like to create a new peaklist?",
                "Adding peaks",
                wx.YES_NO,
            )
            result = dlg.ShowModal()
            if result == wx.ID_NO:
                dlg.Destroy()
                return

            dlg.Destroy()
            # Making a new peaklist, ask the user to create and save a new file in a file dialog
            dlg = wx.FileDialog(
                None,
                "Creating new peaklist",
                wildcard="Peaklists (*.list;*.txt)|*.list;*.txt|All files (*.*)|*.*",
                style=wx.FD_SAVE,
            )
            dlg.SetDirectory(os.getcwd())
            if dlg.ShowModal() == wx.ID_OK:
                peaklist_file = dlg.GetPath()
                try:
                    with open(peaklist_file, "w") as file:
                        pass
                except OSError as error:
                    dlg.Destroy()
                    dlg = wx.MessageDialog(
                        None,
                        "The new peaklist could not be created ({}). Please try "
                        "again in a folder which can be written to.".format(error),
                        "Adding peaks",
                        wx.OK,
                    )
                    dlg.ShowModal()
                    dlg.Destroy()
                    return

                self.AddPeaklist(peaklist_file, new_peaklist=True)

                if self.peak_list_choices == [""]:
                    # The peaklist was not loaded, so there is nothing to add
                    # peaks to and the add peaks mode is not turned on
                    dlg = wx.MessageDialog(
                        None,
                        "The new peaklist {} could not be loaded, so peaks "
                        "cannot be added to it.".format(peaklist_file),
                        "Adding peaks",
                        wx.OK,
                    )
                    dlg.ShowModal()
                    dlg.Destroy()
                    return

            else:
                dlg.Destroy()
                return

        # Updating the current active values
        self.active_add = True
        self.add_peaks_button.SetValue(True)
        self.add_peaks_button.SetForegroundColour(wx.Colour(60, 60, 60))

        # Connect the canvas click event to an add peak function
        self.add_peak_connect = self.main_frame.fig_bore.canvas.mpl_connect(
            "button_press_event", self.on_click_addpeak
        )

        self.selected_peaklist = self.current_peaklist_box.GetValue()

        dlg = wx.MessageDialog(
            None,
            "Peaks can now be added to the peaklist {} by left-clicking the cursor. Please de-select the add button when complete. These peaks will intially be given a chemical shift in the 3rd (bore) dimension equal to 0. This chemical shift can be updated to the correct value using the Move Peaks (z) button. The initial intensity will be set to 0, but this will be updated once the 3rd chemical shift is adjusted to the correct value.".format(
                self.current_peaklist_box.GetValue()
            ),
            "Adding Peaks",
            wx.OK,
        )
        dlg.ShowModal()
        dlg.Destroy()

    def on_click_addpeak(self, event):

        if(event.button != 1):
            return
        
        if event.inaxes is not self.main_frame.ax_bore:
            return

        x, y = self.main_frame.ax_bore.transData.inverted().transform((event.x, event.y))

        if x != None and y != None:

            # Current peaklist
            current_peaklist = self.find_current_peaklist()
            if current_peaklist == None:
                # There is no peaklist to add the peak to
                return

            part = ""
            number = 1
            order = [0, 1]

            if len(self.peak_list_dictionary[current_peaklist]["peak_name"]) > 0:

                peakname = self.peak_list_dictionary[current_peaklist]["peak_name"][-1]
                parts = re.findall(r"[A-Za-z_-]+|\d+", peakname)
                for i, v in enumerate(parts):
                    try:
                        v = int(v)
                        number = v + 1
                    except:
                        part = v
                        if i == 0:
                            order = [1, 0]
                if order == [0, 1]:
                    peakname = str(number) + part
                else:
                    peakname = part + str(number)

                if peakname in self.peak_list_dictionary[current_peaklist]["peak_name"]:
                    peakname = peakname + "_1"

            else:
                peakname = str(number) + part

            shift3 = 0
            intensity = 0

            if self.add_at_local_max_box.GetValue() == True:
                # Put the peak on the nearest maximum of the plane which is
                # shown rather than exactly where it was clicked
                x, y, plane_intensity = self.find_local_maximum_plane(x, y)

                # and give it the chemical shift of the largest point down the
                # bore dimension at that position rather than zero
                bore_shift, bore_intensity = self.find_maximum_bore(x, y)
                if bore_shift != None:
                    shift3 = bore_shift
                    intensity = bore_intensity
                elif plane_intensity != None:
                    intensity = plane_intensity

            if self.check_duplicate_peak(x, y, shift3, current_peaklist) == False:
                # The user has chosen not to have two peaks on top of one
                # another, so the peaklist is left as it was
                return

            self.peak_list_dictionary[current_peaklist]["peak_name"].append(peakname)
            self.peak_list_dictionary[current_peaklist]["shift1"].append(x)
            self.peak_list_dictionary[current_peaklist]["shift2"].append(y)
            self.peak_list_dictionary[current_peaklist]["shift3"].append(shift3)
            self.peak_list_dictionary[current_peaklist]["intensity"].append(intensity)

            # The peak which has just been added becomes the selected peak, so
            # that its bore dimension is shown and it can be moved along the
            # bore straight away
            self.selected_peaklist = current_peaklist
            self.selected_peak_indexes = [
                len(self.peak_list_dictionary[current_peaklist]["peak_name"]) - 1
            ]
            self.selected_peakname = peakname
            self.main_frame.selected_bore_peaks = list(self.selected_peak_indexes)

            # The table is filled in before the plots are redrawn, so that the
            # peak is in the table whatever the plots do
            self.AddToTable()

            # Redraw the plane so that the new peak is shown along with the
            # peaks which are already there
            self.main_frame.OnBoreSlider(wx.EVT_BUTTON)

            # and move the position marker onto the new peak so that its bore
            # dimension is shown
            show_bore_position = getattr(
                self.main_frame, "show_bore_position", None
            )
            if show_bore_position != None:
                show_bore_position(x, y)

            self.AddToTable()

    def OnSelectPeak(self, event):
        """
        This will select a peak so that it can be moved etc
        """
        if self.active_move:
            if self.active_select_peak:
                self.select_peak_button.SetValue(True)
            return
        self.selected_peaklist = self.current_peaklist_box.GetValue()
        if self.selected_peaklist == "":
            dlg = wx.MessageDialog(
                None,
                "No peaklists are loaded, please load a peaklist or perform peak picking and try again.",
                "Warning",
                wx.OK,
            )
            result = dlg.ShowModal()
            dlg.Destroy()
            return

        if self.active_select_peak == True:
            self.active_select_peak = False
            self.selected_peakname = ""
            self.selected_peak_indexes = ['N/A']
            self.select_peak_button.SetValue(False)
            self.main_frame.plot_cross = True
            self.main_frame.fig_bore.canvas.mpl_disconnect(self.select_peak_connect)
            self.main_frame.OnBoreSlider(wx.EVT_BUTTON)
            return

        # Only one of the picking modes can be on at a time
        self.turn_off_picking_modes("select")

        self.active_select_peak = True
        self.select_peak_button.SetValue(True)
        self.main_frame.plot_cross = False

        self.select_peak_connect = self.main_frame.fig_bore.canvas.mpl_connect(
            "button_press_event", self.on_click_selectpeak
        )

    def on_click_selectpeak(self, event):
        """
        If the peak is within a tolerence select the peak
        If multiple peaks are within the tolerence, select the closest
        in terms of pixels on the screen.
        """

        if(event.button != 1):
            return
        if event.inaxes is not self.main_frame.ax_bore:
            return

        # Find the index of the currently selected peaklist
        points = self.main_frame.points[0]

        cont, ind = points.contains(event)
        if cont:
            mouse_coordinates = [event.x, event.y]  # in pixels
            distances = []
            for index in ind["ind"]:
                x = self.peak_list_dictionary[self.current_peaklist_box.GetValue()][
                    "shift1"
                ][index]
                y = self.peak_list_dictionary[self.current_peaklist_box.GetValue()][
                    "shift2"
                ][index]
                x, y = self.main_frame.ax_bore.transData.transform((x, y))
                distance = np.sqrt(
                    (mouse_coordinates[0] - x) ** 2 + (mouse_coordinates[1] - y) ** 2
                )
                distances.append(distance)

            min_index = ind["ind"][np.argmin(np.array(distances))]

            self.selected_peak_indexes = [min_index]
            self.selected_peakname = self.peak_list_dictionary[
                self.current_peaklist_box.GetValue()
            ]["peak_name"][min_index]


            # def strip_name(name):
            #     name1 = name.split('_')[:-2]

            # # Find all peaks which contain the same name after the last _ (these are from the same bore)
            # indexes = []
            # for p in self.peak_list_dictionary[self.current_peaklist_box.GetValue()]["peak_name"]:
            #     initial_name
                


            self.remove_peak = True

        else:
            self.selected_peak_indexes = ["N/A"]

        self.main_frame.OnBoreSlider(wx.EVT_BUTTON)



    def OnMovePeak(self, event):
        """
        This function allows the peak to be moved in the x/y plane of the 2D projection plot
        """

        if self.active_move == True:
            self.active_move = False
            self.move_peaks_button.SetValue(False)
            self.disconnect_bore("move_peak_connect")
            if self.active_select_peak == True:
                # Clicking the plane selects a peak again
                self.select_peak_connect = self.main_frame.fig_bore.canvas.mpl_connect(
                    "button_press_event", self.on_click_selectpeak
                )

            return


        # Temporarily deactivate the ability to select peak or select group
        if len(self.selected_peak_indexes) == 0 or "N/A" in self.selected_peak_indexes:
            # return as there are no selected peaks
            self.active_move = False
            self.move_peaks_button.SetValue(False)
            dlg = wx.MessageDialog(
                self,
                "There are no peaks selected. Please select a peak or a group of peaks and try again.",
                "Warning",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return

        self.eventDict = {}
        for name in dir(wx):
            if name.startswith("EVT_"):
                evt = getattr(wx, name)
                if isinstance(evt, wx.PyEventBinder):
                    self.eventDict[evt.typeId] = name

        if self.active_select_peak == True:
            # While the peak is being moved, clicking the plane moves it rather
            # than selecting another peak
            self.disconnect_bore("select_peak_connect")
            try:
                evt_id = event.GetEventType()
            except AttributeError:
                evt_id = None
            if evt_id != wx.EVT_TOOL_RANGE.typeId:
                dlg = wx.MessageDialog(
                    self,
                    "Please double-left-click a new location to move the selected peak. This can be repeated. Un-toggle the move peaks (x/y) button when completed. (Note: ensure zoom in the matplotlib toolbar is not selected. Zoom before entering move peaks mode if this is required).",
                    "Move Peaks (x/y)",
                    wx.OK,
                )
                dlg.ShowModal()
                dlg.Destroy()

        self.active_move = True
        self.move_peaks_button.SetValue(True)

        # Clicking the plane moves the peak, however the peak was selected
        self.move_peak_connect = self.main_frame.fig_bore.canvas.mpl_connect(
            "button_press_event", self.on_click_movepeak3d
        )



    def on_click_movepeak3d(self, event):
        """
        This function will update the peak position of the selected peak
        depending on where the user clicked.
        """
        if(event.button != 1):
            return
        if event.inaxes is not self.main_frame.ax_bore:
            return
        x, y = self.main_frame.ax_bore.transData.inverted().transform((event.x, event.y))

        if x != None and y != None:
            self.peak_list_dictionary[self.selected_peaklist]["shift1"][
                self.selected_peak_indexes[0]
            ] = x
            self.peak_list_dictionary[self.selected_peaklist]["shift2"][
                self.selected_peak_indexes[0]
            ] = y

            index = 0
            for i, [peaklist, dictionary] in enumerate(
                self.peak_list_dictionary.items()
            ):
                if peaklist == self.selected_peaklist:
                    index = i

            self.main_frame.points[index].set_offsets(
                np.c_[
                    self.peak_list_dictionary[self.selected_peaklist]["shift1"],
                    self.peak_list_dictionary[self.selected_peaklist]["shift2"],
                ]
            )

            # For the x,y,z dimensions, find out the nearest point of the 3D data and then update the intensity of this peak in the peaklist
            shift1 = self.peak_list_dictionary[self.selected_peaklist]["shift1"][self.selected_peak_indexes[0]]
            shift2 = self.peak_list_dictionary[self.selected_peaklist]["shift2"][self.selected_peak_indexes[0]]
            shift3 = self.peak_list_dictionary[self.selected_peaklist]["shift3"][self.selected_peak_indexes[0]]

            index1 = np.argmin(np.abs(self.main_frame.ppms_0 - shift1))
            index2 = np.argmin(np.abs(self.main_frame.ppms_1 - shift2))
            if(self.main_frame.transposed2D == True):
                index1 = np.argmin(np.abs(self.main_frame.ppms_0 - shift2))
                index2 = np.argmin(np.abs(self.main_frame.ppms_1 - shift1))
            index3 = np.argmin(np.abs(self.main_frame.main_frame.ppms_2 - shift3))

            intensity = self.main_frame.main_frame.nmrdata.data[index3][index2][index1]

            self.peak_list_dictionary[self.selected_peaklist]["intensity"][self.selected_peak_indexes[0]] = intensity
            


            # self.main_frame.points[index].set_ydata(self.peak_list_dictionary[self.selected_peaklist]['shift2'])
            self.AddToTable()
        self.main_frame.OnBoreSlider(event)
            # self.main_frame.OnBoreSliderStripPlot(wx.EVT_BUTTON)
            # self.main_frame.UpdateBoreFrame()
            

            

    def OnMovePeakz(self, event):
        """
        This function allows the peak to be moved in the z-dimension of the bore plot
        """
        """
        This function allows the peak to be moved in the x/y plane of the 2D projection plot
        """

        if self.active_movez == True:
            self.active_movez = False
            self.move_peaks_bore_button.SetValue(False)
            self.disconnect_bore("move_peak_connectz")
            if self.active_select_peak == True:
                # Clicking the plane selects a peak again
                self.select_peak_connect = self.main_frame.fig_bore.canvas.mpl_connect(
                    "button_press_event", self.on_click_selectpeak
                )

            return


        # Temporarily deactivate the ability to select peak or select group
        if len(self.selected_peak_indexes) == 0 or "N/A" in self.selected_peak_indexes:
            # return as there are no selected peaks
            self.active_move = False
            self.move_peaks_button.SetValue(False)
            dlg = wx.MessageDialog(
                self,
                "There are no peaks selected. Please select a peak or a group of peaks and try again.",
                "Warning",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return

        self.eventDict = {}
        for name in dir(wx):
            if name.startswith("EVT_"):
                evt = getattr(wx, name)
                if isinstance(evt, wx.PyEventBinder):
                    self.eventDict[evt.typeId] = name

        if self.active_select_peak == True:
            # While the peak is being moved, clicking the plane moves it rather
            # than selecting another peak
            self.disconnect_bore("select_peak_connect")
            try:
                evt_id = event.GetEventType()
            except AttributeError:
                evt_id = None
            if evt_id != wx.EVT_TOOL_RANGE.typeId:
                dlg = wx.MessageDialog(
                    self,
                    "Please double-left-click a new location along the bore dimension to move the selected peak. This can be repeated. Un-toggle the move peaks (z) button when completed. (Note: ensure zoom/pan in the matplotlib toolbar is not selected. Zoom/pan before entering move peaks mode if this is required).",
                    "Move Peaks (z)",
                    wx.OK,
                )
                dlg.ShowModal()
                dlg.Destroy()

        self.active_movez = True
        self.move_peaks_bore_button.SetValue(True)

        # Clicking along the bore moves the peak, however the peak was selected
        self.move_peak_connectz = self.main_frame.fig_bore.canvas.mpl_connect(
            "button_press_event", self.on_click_movepeakz
        )

    def on_click_movepeakz(self, event):
        """
        This function will update the peak position of the selected peak
        depending on where the user clicked (along the bore dimension)
        """
        if(event.button!=1):
            return
        if event.inaxes is self.main_frame.ax_bore:
            return

        x2, y2 = self.main_frame.ax_bore_2.transData.inverted().transform((event.x, event.y))
        x3, y3 = self.main_frame.ax_bore_3.transData.inverted().transform((event.x, event.y))

        if x2 != None and y2 != None:

            self.peak_list_dictionary[self.selected_peaklist]["shift3"][
                self.selected_peak_indexes[0]
            ] = y2

            index = 0
            for i, [peaklist, dictionary] in enumerate(
                self.peak_list_dictionary.items()
            ):
                if peaklist == self.selected_peaklist:
                    index = i

            # For the x,y,z dimensions, find out the nearest point of the 3D data and then update the intensity of this peak in the peaklist
            shift1 = self.peak_list_dictionary[self.selected_peaklist]["shift1"][self.selected_peak_indexes[0]]
            shift2 = self.peak_list_dictionary[self.selected_peaklist]["shift2"][self.selected_peak_indexes[0]]
            shift3 = self.peak_list_dictionary[self.selected_peaklist]["shift3"][self.selected_peak_indexes[0]]

            index1 = np.argmin(np.abs(self.main_frame.ppms_0 - shift1))
            index2 = np.argmin(np.abs(self.main_frame.ppms_1 - shift2))
            if(self.main_frame.transposed2D == True):
                index1 = np.argmin(np.abs(self.main_frame.ppms_0 - shift2))
                index2 = np.argmin(np.abs(self.main_frame.ppms_1 - shift1))
            index3 = np.argmin(np.abs(self.main_frame.main_frame.ppms_2 - shift3))

            intensity = self.main_frame.main_frame.nmrdata.data[index3][index2][index1]

            self.peak_list_dictionary[self.selected_peaklist]["intensity"][self.selected_peak_indexes[0]] = intensity


            # self.main_frame.points[index].set_ydata(self.peak_list_dictionary[self.selected_peaklist]['shift2'])
            self.AddToTable()
        elif x3 != None and y3 != None:

            self.peak_list_dictionary[self.selected_peaklist]["shift3"][
                self.selected_peak_indexes[0]
            ] = y3

            index = 0
            for i, [peaklist, dictionary] in enumerate(
                self.peak_list_dictionary.items()
            ):
                if peaklist == self.selected_peaklist:
                    index = i

            # For the x,y,z dimensions, find out the nearest point of the 3D data and then update the intensity of this peak in the peaklist
            shift1 = self.peak_list_dictionary[self.selected_peaklist]["shift1"][self.selected_peak_indexes[0]]
            shift2 = self.peak_list_dictionary[self.selected_peaklist]["shift2"][self.selected_peak_indexes[0]]
            shift3 = self.peak_list_dictionary[self.selected_peaklist]["shift3"][self.selected_peak_indexes[0]]

            index1 = np.argmin(np.abs(self.main_frame.ppms_0 - shift1))
            index2 = np.argmin(np.abs(self.main_frame.ppms_1 - shift2))
            if(self.main_frame.transposed2D == True):
                index1 = np.argmin(np.abs(self.main_frame.ppms_0 - shift2))
                index2 = np.argmin(np.abs(self.main_frame.ppms_1 - shift1))
            index3 = np.argmin(np.abs(self.main_frame.main_frame.ppms_2 - shift3))

            intensity = self.main_frame.main_frame.nmrdata.data[index3][index2][index1]

            self.peak_list_dictionary[self.selected_peaklist]["intensity"][self.selected_peak_indexes[0]] = intensity


            # self.main_frame.points[index].set_ydata(self.peak_list_dictionary[self.selected_peaklist]['shift2'])
            self.AddToTable()
        self.main_frame.OnBoreSlider(event)
            # self.main_frame.OnBoreSliderStripPlot(wx.EVT_BUTTON)
            # self.main_frame.UpdateBoreFrame()


    def OnRemovePeaks(self, event):
        """
        If there is a current peak or peaks selected, then remove these peaks
        from the dictionary.

        If a peak or peaks are selected in the table of the Peak List window ask
        if the user if they want to remove these peaks.
        """

        indexes = self.find_peaks_to_move()
        if len(indexes) == 0:
            return

        dictionary = self.peak_list_dictionary[self.selected_peaklist]

        # Removing from the end so that the indexes of the peaks which are
        # still to be removed do not move, and taking everything which is held
        # for the peak so that the lists stay in step with one another
        for index in sorted(indexes, reverse=True):
            for key in ["peak_name", "shift1", "shift2", "shift3", "intensity"]:
                try:
                    del dictionary[key][index]
                except (KeyError, IndexError, TypeError):
                    pass

        self.remove_peak = False
        self.selected_peak_indexes = ["N/A"]
        self.selected_peakname = ""
        self.main_frame.selected_bore_peaks = []

        self.main_frame.OnBoreSlider(wx.EVT_BUTTON)
        self.AddToTable()


    def on_click_movepeak(self, event):
        """
        This function will update the peak position of the selected peak
        depending on where the user clicked.
        """

        if event.inaxes is not self.main_frame.ax_bore:
            return
        x, y = self.main_frame.ax_bore.transData.inverted().transform((event.x, event.y))

        if x != None and y != None:
            self.peak_list_dictionary[self.selected_peaklist]["shift1"][
                self.selected_peak_indexes[0]
            ] = x
            self.peak_list_dictionary[self.selected_peaklist]["shift2"][
                self.selected_peak_indexes[0]
            ] = y

            index = 0
            for i, [peaklist, dictionary] in enumerate(
                self.peak_list_dictionary.items()
            ):
                if peaklist == self.selected_peaklist:
                    index = i

            self.main_frame.points[index].set_offsets(
                np.c_[
                    self.peak_list_dictionary[self.selected_peaklist]["shift1"],
                    self.peak_list_dictionary[self.selected_peaklist]["shift2"],
                ]
            )
            # self.main_frame.points[index].set_ydata(self.peak_list_dictionary[self.selected_peaklist]['shift2'])
            self.main_frame.UpdateBoreFrame()
            self.AddToTable()

    def on_press_movepeak(self, event):

        x, y = self.main_frame.ax.transData.inverted().transform((event.x, event.y))
        if x != None and y != None:
            self.start_point_move = (x, y)
            self.x_init = copy.deepcopy(
                self.peak_list_dictionary[self.selected_peaklist]["shift1"]
            )
            self.y_init = copy.deepcopy(
                self.peak_list_dictionary[self.selected_peaklist]["shift2"]
            )

    def on_motion_movepeak(self, event):
        if self.start_point_move == None:
            return

        x, y = self.main_frame.ax.transData.inverted().transform((event.x, event.y))

        if x != None and y != None:

            # Update rectangle size
            x0, y0 = self.start_point_move
            x1, y1 = x, y
            x_change = x1 - x0
            y_change = y1 - y0

            for index in self.selected_peak_indexes:

                self.peak_list_dictionary[self.selected_peaklist]["shift1"][index] = (
                    self.x_init[index] + x_change
                )
                self.peak_list_dictionary[self.selected_peaklist]["shift2"][index] = (
                    self.y_init[index] + y_change
                )

            ind = 0
            for i, [peaklist, dictionary] in enumerate(
                self.peak_list_dictionary.items()
            ):
                if peaklist == self.selected_peaklist:
                    ind = i

            self.main_frame.points[ind].set_offsets(
                np.c_[
                    self.peak_list_dictionary[self.selected_peaklist]["shift1"],
                    self.peak_list_dictionary[self.selected_peaklist]["shift2"],
                ]
            )
            # self.main_frame.points[index].set_ydata(self.peak_list_dictionary[self.selected_peaklist]['shift2'])
            self.main_frame.UpdateFrame()

    def on_release_movepeak(self, event):
        # self.on_motion_movepeak(event)
        self.AddToTable()
        self.start_point_move = None

    def OnFindPeaks(self, event):
        """
        If one peak is currently selected in the table, then zoom in to this
        peak and select it.
        Before doing this, the code will turn off all active toggled buttons from
        the Peak List frame.
        """

        if self.active_add == True:
            self.active_add = False
            self.add_peaks_button.SetValue(False)
            self.main_frame.fig_bore.canvas.mpl_disconnect(self.add_peak_connect)
            return

        row = self.grid.GetGridCursorRow()
        peak_name = self.grid.GetCellValue(row, 0)
        shift1 = self.grid.GetCellValue(row, 1)
        shift2 = self.grid.GetCellValue(row, 2)

        # Zoom in on grid selected peak and then select it in the plot.
        width = 0.05  # ppm
        height = 0.05  # ppm

        xmin = float(shift1) - width
        xmax = float(shift1) + width

        ymin = float(shift2) - height
        ymax = float(shift2) + height

        self.main_frame.toolbar_bore.push_current()

        self.zoom_to_region(self.main_frame.ax_bore, [xmin, xmax], [ymin, ymax])

        self.simulate_peak_selection_click(float(shift1), float(shift2))

        self.main_frame.UpdateBoreFrame()

        self.main_frame.toolbar_bore.push_current()

    def simulate_peak_selection_click(self, shift1, shift2):
        # evt_down = wx.MouseEvent(wx.wxEVT_LEFT_DOWN)
        # evt_down.SetX(shift1)
        # evt_down.SetY(shift2)
        # evt_down.SetEventObject(self.main_frame.canvas_bore)
        # self.main_frame.on_click_bore(evt_down)
        disp_x, disp_y = self.main_frame.ax_bore.transData.transform((shift1, shift2))
        event = MPLMouseEvent(
            name="button_press_event",
            canvas=self.main_frame.canvas_bore,
            x=disp_x,
            y=disp_y,
            button=1,
            key=None,
            step=0,
            dblclick=False,
            guiEvent=None,
        )
        event.inaxes = self.main_frame.ax_bore
        event.xdata = shift1
        event.ydata = shift2
        self.main_frame.on_pick(event)
        self.main_frame.on_click_bore(event)

    def find_current_peaklist(self):
        """
        The peaklist which is loaded. The box showing it is what everything
        else works from, so it is used in preference to the name recorded when
        the peaklist was read.
        """
        peaklist = self.current_peaklist_box.GetValue()
        if peaklist in self.peak_list_dictionary:
            return peaklist

        peaklist = getattr(self, "peak_list", "")
        if peaklist in self.peak_list_dictionary:
            return peaklist

        return None

    def find_axis_labels(self):
        """
        The names of the three dimensions of the peaklist, taken from the axes
        of the plots so that the peaklist says the same as what is shown: the
        axis across the plane, the axis up it, and the bore dimension.
        """
        labels = ["", "", ""]

        try:
            labels[0] = self.main_frame.ax_bore.get_xlabel()
            labels[1] = self.main_frame.ax_bore.get_ylabel()
            labels[2] = self.main_frame.ax_bore_2.get_ylabel()
        except (AttributeError, RuntimeError):
            pass

        return labels

    def find_column_labels(self) -> list:
        """
        The names of the peaklist columns, used both for the table and for the
        header of a saved peaklist so that the two say the same thing.
        """
        labels = self.find_axis_labels()

        return [
            "Peak name",
            find_axis_column_label(labels[0], 1),
            find_axis_column_label(labels[1], 2),
            find_axis_column_label(labels[2], 3),
            "Intensity",
            "Ambiguity",
            "Alternatives (ppm)",
        ]

    def update_column_labels(self):
        """
        Show which dimension each shift column holds, using the names of the
        axes of the plots.
        """
        try:
            for column, label in enumerate(self.find_column_labels()):
                if column < self.grid.GetNumberCols():
                    self.grid.SetColLabelValue(column, label)
        except (RuntimeError, AttributeError):
            pass

    def find_table_order(self, peaklist):
        """
        The peaks of a peaklist in the order the table shows them, as indexes
        into its lists. The table is sorted by the number in the peak name.
        """
        def extract_number(name):
            match = re.match(r"(\d+)", name)
            return int(match.group(1)) if match else float("inf")

        names = self.peak_list_dictionary[peaklist]["peak_name"]

        return [
            index
            for index, name in sorted(
                enumerate(names), key=lambda pair: extract_number(pair[1])
            )
        ]

    def read_notes(self, dictionary, line):
        """
        Read what a peaklist file holds about how sure a peak's position down
        the bore is. A peaklist written before this was worked out has neither
        column, so the peak simply has nothing in them.
        """
        columns = split_peaklist_header(str(line).split("\n")[0])

        for i, key in enumerate(["ambiguity", "alternatives"]):
            value = ""
            try:
                value = str(columns[5 + i]).strip()
            except IndexError:
                value = ""

            if value.lower() == "nan":
                value = ""

            dictionary[key].append(value)

    def find_note(self, peaklist, key, index) -> str:
        """
        What is held for a peak in one of the columns which say how sure its
        position down the bore is. A peaklist which has not been through that
        analysis has nothing in them.
        """
        try:
            value = self.peak_list_dictionary[peaklist][key][index]
        except (KeyError, IndexError, TypeError):
            return ""

        if value == None:
            return ""

        return str(value)

    def find_peak_rows(self):
        """
        The peaks of the loaded peaklist as rows of text, in the order the table
        shows them. The table and the file which is saved are both filled from
        here, so what is saved is what the peaklist holds rather than whatever
        the table happens to be showing.
        """
        peaklist = self.find_current_peaklist()
        if peaklist == None:
            return []

        dictionary = self.peak_list_dictionary[peaklist]
        rows = []

        for index in self.find_table_order(peaklist):
            try:
                rows.append([
                    dictionary["peak_name"][index],
                    "{:.5f}".format(dictionary["shift1"][index]),
                    "{:.5f}".format(dictionary["shift2"][index]),
                    "{:.5f}".format(dictionary["shift3"][index]),
                    "{:.5e}".format(dictionary["intensity"][index]),
                    self.find_note(peaklist, "ambiguity", index),
                    self.find_note(peaklist, "alternatives", index),
                ])
            except (IndexError, KeyError, TypeError, ValueError):
                continue

        return rows

    def find_save_location(self):
        """
        The folder and file name shown when saving, taken from the file the
        peaklist was read from or created as. The current directory and an
        untitled file are used when there is no peaklist.
        """
        directory = pathlib.Path(os.getcwd())
        file_name = "Untitled.list"

        # The whole path of the file the peaklist came from, rather than the
        # shortened name which is shown in the window
        peaklist = getattr(self, "peaklist_path", "")
        if peaklist == "" or peaklist == None:
            peaklist = self.current_peaklist_box.GetValue()

        if peaklist == "":
            return directory, file_name

        try:
            peaklist_path = pathlib.Path(peaklist).expanduser()
        except TypeError:
            return directory, file_name

        if peaklist_path.name != "":
            file_name = peaklist_path.name

        parent = peaklist_path.absolute().parent
        if parent.is_dir() == True:
            directory = parent

        return directory, file_name

    def mark_saved(self):
        """
        Remember the peaklist as it is now, so that changes made after this are
        noticed as changes which have not been saved.
        """
        peaklist = self.current_peaklist_box.GetValue()
        self.saved_peaklist = copy.deepcopy(
            self.peak_list_dictionary.get(peaklist, {})
        )
        self.update_saved_button()

    def find_unsaved_changes(self) -> bool:
        """
        Whether the peaklist has changed since it was read or last saved.
        """
        peaklist = self.current_peaklist_box.GetValue()
        if peaklist not in self.peak_list_dictionary:
            return False

        return self.peak_list_dictionary[peaklist] != getattr(
            self, "saved_peaklist", {}
        )

    def update_saved_button(self):
        """
        Mark the save button while the peaklist has changes which have not been
        saved, so that it is clear whether the file on disk is up to date.
        """
        try:
            if self.find_unsaved_changes() == True:
                self.save_peaks_button.SetLabel("Save *")
            else:
                self.save_peaks_button.SetLabel("Save")
        except (RuntimeError, AttributeError):
            pass

    def OnClose(self, event):
        """
        Offer to save the peaklist when the window is closed while it has
        changes which have not been saved.
        """
        self.turn_off_togglebuttons()

        if self.find_unsaved_changes() == True:
            dlg = wx.MessageDialog(
                self,
                "The peaklist {} has changes which have not been saved. Would "
                "you like to save it before closing?".format(
                    self.current_peaklist_box.GetValue()
                ),
                "Save peaklist",
                wx.YES_NO | wx.CANCEL | wx.ICON_QUESTION,
            )
            result = dlg.ShowModal()
            dlg.Destroy()

            if result == wx.ID_CANCEL:
                # Staying open so that the peaks are not lost
                try:
                    event.Veto()
                except AttributeError:
                    pass
                return

            if result == wx.ID_YES:
                self.OnSave(wx.EVT_BUTTON)

        # The peaks stay on the plots, and the window holds its peaklist until
        # the peaklist is removed
        self.hide_window(event)

    def hide_window(self, event):
        """
        Put the window away rather than destroying it, so that everything it
        holds is still there when it is opened again.
        """
        self.Hide()

        try:
            if event.CanVeto() == True:
                event.Veto()
                return
        except AttributeError:
            pass

        # The window cannot stay, so it is destroyed as it used to be
        self.Destroy()

    def OnSave(self, event, peaklist_file=''):
        """
        Save the peaklist back into its own file, which is the file it was read
        from or created as, or the one it was last saved as. The user is asked
        for a file only when the peaklist does not have one.
        """
        tell_the_user = peaklist_file == ''

        if peaklist_file == '':
            peaklist_file = getattr(self, "peaklist_path", "")

        if peaklist_file == '' or peaklist_file == None:
            # The peaklist has no file of its own to be saved into
            return self.OnSaveAs(event)

        self.save_peaklist(peaklist_file, False, tell_the_user)

    def OnSaveAs(self, event):
        """
        Ask for a file to save the peaklist as, which becomes the file the
        peaklist is saved into from then on.
        """
        save_2d_plane = False

        directory, file_name = self.find_save_location()
        dlg = wx.FileDialog(
            self,
            "Save the peaklist as",
            wildcard="Peaklists (*.list;*.tab;*.txt)|*.list;*.tab;*.txt|All files (*.*)|*.*",
            style=wx.FD_SAVE,
        )
        dlg.SetDirectory(str(directory))
        dlg.SetFilename(str(file_name))
        if dlg.ShowModal() == wx.ID_OK:
            peaklist_file = dlg.GetPath()
        else:
            dlg.Destroy()
            return
        dlg.Destroy()

        if peaklist_file == "" or os.path.isdir(peaklist_file) == True:
            # No name was given, so there is nowhere to save the peaklist
            dlg = wx.MessageDialog(
                None,
                "No name was given for the peaklist, so it has not been "
                "saved. Please try again and give the peaklist a name.",
                "Save peaklist",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return

        message = 'Would you like to save the full 3D peaklist (click yes), or would you like to save the 2D reference plane (click no)?'
        dlg = wx.MessageDialog(None, message, "Save peaklist", wx.YES_NO)
        result = dlg.ShowModal()
        dlg.Destroy()
        if(result == wx.ID_NO):
            save_2d_plane = True

        self.save_peaklist(peaklist_file, save_2d_plane, True)

    def save_peaklist(self, peaklist_file, save_2d_plane, tell_the_user):
        """
        Write the peaklist into a file, saying where it has gone and what went
        wrong if it could not be written.
        """
        try:
            number_of_peaks = self.write_peaklist(peaklist_file, save_2d_plane)
        except OSError as error:
            dlg = wx.MessageDialog(
                None,
                "The peaklist could not be saved as {} ({}).".format(
                    peaklist_file, error
                ),
                "Save peaklist",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return

        if number_of_peaks == 0:
            dlg = wx.MessageDialog(
                None,
                "There are no peaks in the peaklist, so the file {} has been "
                "written empty.".format(peaklist_file),
                "Save peaklist",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
        elif tell_the_user == True:
            # Say where the peaklist has gone, as it is not always the file it
            # was read from
            dlg = wx.MessageDialog(
                None,
                "{} peaks have been saved as {}.".format(
                    number_of_peaks, peaklist_file
                ),
                "Save peaklist",
                wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()

        if save_2d_plane == False:
            # The peaklist on disk now holds what the window holds, and later
            # saves go back to the file it has just been written to
            self.peaklist_path = str(pathlib.Path(peaklist_file).absolute())
            self.mark_saved()

    def write_peaklist(self, peaklist_file, save_2d_plane=False):
        """
        Write the peaklist into a file, in the order the table shows it.
        """
        rows = self.find_peak_rows()
        labels = self.find_column_labels()

        with open(peaklist_file, "w") as file:
            # The header names the columns in the same way as the table. Only
            # the two shifts of the plane are saved for a reference plane
            if save_2d_plane == False:
                file.write(" \t ".join(labels) + "\n")
            else:
                file.write(" \t ".join(labels[:3]) + "\n")

            for row in rows:
                if save_2d_plane == False:
                    file.write(" \t ".join(row) + "\n")
                else:
                    file.write(" \t ".join(row[:3]) + "\n")

        return len(rows)


class ResolveAmbiguityDialog(wx.Dialog):
    """
    A window listing the peaks whose position down the bore dimension is not
    certain. Each peak can be shown in the spectrum, so that the peak tools can
    be used to look at it, and then given one of the positions which were found,
    kept as it is, or accepted where it now sits.

    The window does not take over the application, so the spectrum and the peak
    window can be used while it is open.
    """

    KEEP = "Leave as it is"
    CURRENT = "Accept the position it has now"

    def __init__(self, parent, peaks):
        """
        peaks is a list of dictionaries holding the name of the peak, the text
        saying what is uncertain about it, where it is in the peaklist and the
        positions it could have, each as {"shift", "intensity", "confidence"}.
        """
        wx.Dialog.__init__(
            self,
            parent,
            title="Resolve peaks down the bore",
            size=(860, 480),
            style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER,
        )

        self.peaks_window = parent
        self.peaks = peaks
        self.choices = []
        self.positions = []
        self.show_buttons = []

        sizer = wx.BoxSizer(wx.VERTICAL)

        sizer.Add(
            wx.StaticText(
                self,
                -1,
                "These peaks could have more than one position down the bore "
                "dimension, or found\nfewer maxima than expected. Show a peak to "
                "look at it in the spectrum, where the\nselect and move tools can "
                "be used, then give it a position or leave it as it is.",
            ),
            0,
            wx.ALL,
            10,
        )

        scrolled = wx.ScrolledWindow(self, -1, style=wx.VSCROLL)
        scrolled.SetScrollRate(0, 10)
        rows = wx.FlexGridSizer(cols=5, hgap=10, vgap=6)
        rows.AddGrowableCol(4)

        for name in ["", "Peak", "What is uncertain", "Position now", "Give it"]:
            label = wx.StaticText(scrolled, -1, name)
            label.SetFont(label.GetFont().Bold())
            rows.Add(label, 0, wx.ALIGN_CENTER_VERTICAL)

        for i, peak in enumerate(peaks):
            show = wx.Button(scrolled, -1, "Show", size=(60, 24))
            show.Bind(wx.EVT_BUTTON, self.on_show)
            show.peak_number = i
            rows.Add(show, 0, wx.ALIGN_CENTER_VERTICAL)
            self.show_buttons.append(show)

            rows.Add(
                wx.StaticText(scrolled, -1, str(peak["name"])),
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )
            rows.Add(
                wx.StaticText(scrolled, -1, str(peak["note"])),
                0,
                wx.ALIGN_CENTER_VERTICAL,
            )

            position = wx.StaticText(scrolled, -1, self.find_position_text(peak))
            rows.Add(position, 0, wx.ALIGN_CENTER_VERTICAL)
            self.positions.append(position)

            choice = wx.Choice(scrolled, -1, choices=self.find_choices(peak))
            choice.SetSelection(0)
            rows.Add(choice, 0, wx.EXPAND)
            self.choices.append(choice)

        scrolled.SetSizer(rows)
        sizer.Add(scrolled, 1, wx.EXPAND | wx.ALL, 10)

        buttons = wx.BoxSizer(wx.HORIZONTAL)
        self.apply_button = wx.Button(self, -1, "Apply")
        self.apply_button.Bind(wx.EVT_BUTTON, self.on_apply)
        self.apply_button.SetToolTip(
            "Give each peak the position chosen for it. The peaks left alone are "
            "not changed."
        )
        self.refresh_button = wx.Button(self, -1, "Refresh")
        self.refresh_button.Bind(wx.EVT_BUTTON, self.on_refresh)
        self.refresh_button.SetToolTip(
            "Read the peaks again, to show where they are now after being moved "
            "in the spectrum."
        )
        self.close_button = wx.Button(self, wx.ID_CANCEL, "Close")
        buttons.Add(self.apply_button)
        buttons.AddSpacer(10)
        buttons.Add(self.refresh_button)
        buttons.AddSpacer(10)
        buttons.Add(self.close_button)
        sizer.Add(buttons, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 10)

        self.SetSizer(sizer)

        self.Bind(wx.EVT_CLOSE, self.on_close)

    def find_position_text(self, peak) -> str:
        """
        Where the peak sits down the bore dimension at the moment.
        """
        shift = peak.get("shift3")

        if shift == None or shift == 0:
            return "not placed"

        return "{:.4f} ppm".format(float(shift))

    def find_choices(self, peak) -> list:
        """
        How the positions a peak could have are offered. Leaving it alone comes
        first, so that a window which is not touched changes nothing.
        """
        choices = [self.KEEP, self.CURRENT]

        for candidate in peak["candidates"]:
            text = "{:.4f} ppm".format(candidate["shift"])
            if candidate.get("intensity") != None:
                text += "   intensity {:.3e}".format(candidate["intensity"])
            if candidate.get("confidence") != None:
                text += "   confidence {:.2f}".format(candidate["confidence"])
            if candidate.get("extra") == True:
                text += "   (below the threshold)"
            choices.append(text)

        return choices

    def find_answers(self) -> list:
        """
        What was chosen for each peak: None to leave it alone, "current" to
        accept where it now sits, or the position to give it.
        """
        answers = []

        for i, choice in enumerate(self.choices):
            selection = choice.GetSelection()

            if selection <= 0:
                answers.append(None)
                continue

            if selection == 1:
                answers.append("current")
                continue

            try:
                answers.append(self.peaks[i]["candidates"][selection - 2])
            except IndexError:
                answers.append(None)

        return answers

    def on_show(self, event):
        """
        Show the peak of this row in the spectrum, so that it can be looked at
        with the peak tools before a decision is made.
        """
        number = getattr(event.GetEventObject(), "peak_number", None)
        if number == None:
            return

        self.show_peak(number)

    def show_peak(self, number):
        """
        Select the peak in the peak window, move the position marker onto it and
        show the positions it could have down the bore.
        """
        try:
            peak = self.peaks[number]
        except IndexError:
            return

        self.peaks_window.show_peak_for_resolving(peak)

        try:
            self.positions[number].SetLabel(
                self.find_position_text(
                    self.peaks_window.find_peak_for_resolving(peak)
                )
            )
        except (RuntimeError, AttributeError):
            pass

    def on_apply(self, event):
        """
        Give each peak the position which was chosen for it.
        """
        self.peaks_window.apply_resolutions(self.peaks, self.find_answers())

        for i, peak in enumerate(self.peaks):
            try:
                self.positions[i].SetLabel(
                    self.find_position_text(
                        self.peaks_window.find_peak_for_resolving(peak)
                    )
                )
                self.choices[i].SetSelection(0)
            except (RuntimeError, AttributeError):
                pass

    def on_refresh(self, event):
        """
        Read the peaks again, so that peaks which have been moved in the
        spectrum are shown where they now are.
        """
        self.peaks_window.reopen_resolve_dialog()

    def on_close(self, event):
        """
        The positions a peak could have are no longer shown once the window has
        gone.
        """
        self.peaks_window.finish_resolving()

        self.Destroy()
