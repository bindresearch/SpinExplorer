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

# Everything which is held for each peak of a 2D peaklist, so that the lists
# stay the same length as peaks are added and removed
PEAK_ENTRY_KEYS = [
    "peak_name",
    "shift1",
    "shift2",
    "intensity",
] + PEAK_LINEWIDTH_KEYS


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
        The pan/zoom toolbar of the spectrum this peak window belongs to.
        """
        return getattr(getattr(self, "main_frame", None), "toolbar", None)

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
            "OnMovePeakBore",
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

        # Turn off all togglebuttons before closing (unlinks the matplotlib canvas of peaklist specific tasks)
        self.turn_off_togglebuttons()

        # Check if a Fit peaks result window is open
        continue_closing = self.check_fit_window()
        if(continue_closing==False):
            return
        
        # Telling the user that closing will lead to all unsaved changes to each peaklist being lost, asking if they wish to continue
        # Asking the user if they would like to save a session containing these peaklists
        dlg = wx.MessageDialog(
                    self,
                    "Closing the peaklist window will lead to all unsaved changes to each peaklist being lost. Would you like to continue?"
                    ,
                    "Warning",
                    wx.YES_NO,
        )
        res = dlg.ShowModal()
        if(res == wx.ID_NO):    
            dlg.Destroy()
            return
        dlg.Destroy()

        # Asking the user if they would like to save a session containing these peaklists
        dlg = wx.MessageDialog(
                    self,
                    "Would you like to save the a session containing the current peaklists before closing?"
                    ,
                    "Save",
                    wx.YES_NO,
        )
        res = dlg.ShowModal()
        if(res == wx.ID_YES):
            self.main_frame.OnSaveSessionButton2D(wx.EVT_BUTTON)
        dlg.Destroy()


        self.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True, showpeaks=False)
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
                                    name1 = line[1]
                                    name2 = line[2]
                                    try:
                                        name3 = line[3]
                                    except:
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

        if len(dictionary["peak_name"]) == 0 and new_peaklist == False:
            self.peaklist_error_message()
            return None
        
        
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
                None, "Creating new peaklist", wildcard="*.list|*.txt", style=wx.FD_SAVE
            )
            dlg.SetDirectory(os.getcwd())
            if dlg.ShowModal() == wx.ID_OK:
                peaklist_file = dlg.GetPath()
                with open(peaklist_file, "w") as file:
                    pass

                self.AddPeaklist(peaklist_file, new_peaklist=True)

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

        self.main_frame.ax.set_xlim([xmax, xmin])
        self.main_frame.ax.set_ylim([ymax, ymin])
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
            return

        self.pad_peaklist(current_peaklist)

        dictionary = self.peak_list_dictionary[current_peaklist]
        linewidths = self.find_linewidths_written(current_peaklist)

        with open(peaklist_file, "w") as file:
            if self.names[current_peaklist] != ["", ""]:
                header = "Peak \t {} \t {} \t {}".format(
                    self.names[current_peaklist][0],
                    self.names[current_peaklist][1],
                    "Intensity",
                )
                if linewidths == True:
                    header += " \t Linewidth1(Hz) \t Linewidth1(ppm)"
                    header += " \t Linewidth2(Hz) \t Linewidth2(ppm)"
                file.write(header + "\n")

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

    def OnSave(self, event, save_after_picking=False):
        """
        Provide a FileDialog where the user can chose the name for
        the peaklist.
        The peaklist will then be saved.
        """

        directory, file_name = self.find_save_location()

        if(save_after_picking==False):
            dlg = wx.FileDialog(self, "Select the folder and name to save the peaklist as.", wildcard="", style=wx.FD_SAVE)
            dlg.SetDirectory(str(directory))
            dlg.SetFilename(str(file_name))
            if dlg.ShowModal() == wx.ID_OK:
                peaklist_file = dlg.GetPath()
            else:
                dlg.Destroy()
                return
        else:
            peaklist_file = self.peaklist_name_box.GetValue()

        self.write_peaklist(peaklist_file)

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
        width = 1000
        height = 400
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.panel_peaklist = wx.Panel(self, -1)
        self.main_peaklist_sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_peaklist_sizer)

        self.set_initial_values()
        self.follow_peak_modes()
        self.make_peaklist_window()
        self.watch_navigation_toolbar()
        self.update_mode_buttons()
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

        self.peaklist_selection_text = wx.StaticText(self, -1, "Selected Peaklist:")

        self.current_peaklist_box = wx.TextCtrl(
            self, -1, value='', size=(250, 20), style = wx.TE_READONLY
        )


        self.row2_label = wx.StaticBox(
            self,
            -1,
            "Manipulate Peaklists: (shorcuts for Mac - cmd+key)",
        )
        self.row2 = wx.StaticBoxSizer(self.row2_label, wx.HORIZONTAL)

        self.add_peaks_button = wx.ToggleButton(self.row2_label, label="Add Peaks (a)")
        self.add_peaks_button.Bind(wx.EVT_TOGGLEBUTTON, self.OnAddPeaks)
        ID_BUTTON_a = wx.NewIdRef()

        self.select_peak_button = wx.ToggleButton(self.row2_label, label="Select Peak (s)")
        self.select_peak_button.Bind(wx.EVT_TOGGLEBUTTON, self.OnSelectPeak)
        ID_BUTTON_s = wx.NewIdRef()

        self.add_borepeak_button = wx.ToggleButton(self.row2_label, label="Add Bore Peak (b)")
        self.add_borepeak_button.Bind(wx.EVT_TOGGLEBUTTON, self.OnAddBorePeak)
        ID_BUTTON_b = wx.NewIdRef()

        # self.select_peaks_button = wx.ToggleButton(self, label="Select Peak Group (g)")
        # self.select_peaks_button.Bind(wx.EVT_TOGGLEBUTTON, self.OnSelectPeaks)
        # ID_BUTTON_g = wx.NewIdRef()

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

        # Creating an accelerator table for keyboard shortcuts for the buttons
        accelerator_table = wx.AcceleratorTable(
            [
                (wx.ACCEL_CTRL, ord("a"), ID_BUTTON_a),
                (wx.ACCEL_CTRL, ord("r"), ID_BUTTON_r),
                (wx.ACCEL_CTRL, ord("b"), ID_BUTTON_b),
                (wx.ACCEL_CTRL, ord("f"), ID_BUTTON_f),
                (wx.ACCEL_CTRL, ord("m"), ID_BUTTON_m),
                (wx.ACCEL_CTRL, ord("z"), ID_BUTTON_z),
                (wx.ACCEL_CTRL, ord("s"), ID_BUTTON_s),
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

        self.main_frame.Bind(wx.EVT_MENU, self.OnAddPeaks, id=ID_BUTTON_a)
        self.main_frame.Bind(wx.EVT_MENU, self.OnAddBorePeak, id=ID_BUTTON_b)
        self.main_frame.Bind(wx.EVT_MENU, self.OnRemovePeaks, id=ID_BUTTON_r)
        self.main_frame.Bind(wx.EVT_MENU, self.OnMovePeak, id=ID_BUTTON_m)
        self.main_frame.Bind(wx.EVT_MENU, self.OnMovePeakz, id=ID_BUTTON_z)
        self.main_frame.Bind(wx.EVT_MENU, self.OnFindPeaks, id=ID_BUTTON_f)
        self.main_frame.Bind(wx.EVT_MENU, self.OnSelectPeak, id=ID_BUTTON_s)

        self.save_peaks_button = wx.Button(self.row2_label, label="Save")
        self.save_peaks_button.Bind(wx.EVT_BUTTON, self.OnSave)


        self.row2.AddSpacer(5)
        self.row2.Add(self.add_peaks_button)
        self.row2.AddSpacer(5)
        self.row2.Add(self.select_peak_button)
        self.row2.AddSpacer(5)
        self.row2.Add(self.add_borepeak_button)
        self.row2.AddSpacer(5)
        self.row2.Add(self.move_peaks_button)
        self.row2.AddSpacer(5)
        self.row2.Add(self.move_peaks_bore_button)
        self.row2.AddSpacer(5)
        self.row2.Add(self.remove_peaks_button)
        self.row2.AddSpacer(5)
        self.row2.Add(self.find_peak_button)
        self.row2.AddSpacer(5)
        self.row2.Add(self.save_peaks_button)


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


        self.row1 = wx.BoxSizer(wx.HORIZONTAL)
        self.row1.Add(self.add_peaklist_button)
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
        self.grid.CreateGrid(5, 5)

        self.grid.SetColLabelValue(0, "Peak name")
        self.grid.SetColLabelValue(1, "Shift 1 (ppm)")
        self.grid.SetColLabelValue(2, "Shift 2 (ppm)")
        self.grid.SetColLabelValue(3, "Shift 3 (ppm)")
        self.grid.SetColLabelValue(4, "Intensity")

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
        

        if(self.reference_plane==False):
        
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
            
            if(self.main_frame.transposed2D == False):
                x = self.main_frame.uc0.ppm(peaks["X_AXIS"])
                y = self.main_frame.uc1.ppm(peaks["Y_AXIS"])
            else:
                x = self.main_frame.uc0.ppm(peaks["Y_AXIS"])
                y = self.main_frame.uc1.ppm(peaks["X_AXIS"])
            z = self.main_frame.main_frame.uc2.ppm(peaks["Z_AXIS"])

            intensities = []
            for p in peaks:
                intensities.append(p[-1])


            picked_peak_array = []
            for i, xval in enumerate(x):
                picked_peak_array.append([xval, y[i], z[i], intensities[i]])

            dictionary = {}
            dictionary["peak_name"] = []
            dictionary["shift1"] = []
            dictionary["shift2"] = []
            dictionary["shift3"] = []
            dictionary["intensity"] = []
            for i, peak in enumerate(picked_peak_array):
                dictionary["peak_name"].append(str(i+1))
                dictionary["shift1"].append(peak[0])
                dictionary["shift2"].append(peak[1])
                dictionary["shift3"].append(peak[2])
                dictionary["intensity"].append(peak[3])

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
                self.OnLoadReferencePlane(wx.EVT_TOGGLEBUTTON)
            elif(check_reference_plane == True):
                ppms0_old = ppms0
                ppms1_old = ppms1
                ppms0 = ppms1_old
                ppms1 = ppms0_old


            # For each of the reference plane peaks, the index of the chemical
            # shifts closest to this value needs to be known so a list of 1D
            # data points can be individually fitted by the nmrglue peak picker

            indexes0 = []
            indexes1 = []

            for i, ppm0 in enumerate(ppms0):
                if(len(x)==data.shape[-1] and len(y)==data.shape[-2]):
                    index0 = np.argmin(np.abs(x-ppm0))
                    index1 = np.argmin(np.abs(y-ppms1[i]))
                else:
                    index1 = np.argmin(np.abs(x-ppm0))
                    index0 = np.argmin(np.abs(y-ppms1[i]))
                indexes0.append(index0)
                indexes1.append(index1)
                
            

            data_1D_slices = []
            for i, index in enumerate(indexes0):
                data_1D_slices.append(data[:, indexes1[i], index])

            data_1D_slices = np.array(data_1D_slices)


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
            
            names1 = []
            x = []
            y = []
            z = []
            intensities = []
            for i, peak in enumerate(peaks):
                try:
                    if(peak==0):
                        # No peak picked for this slice
                        continue
                except:
                    pass
                if(self.main_frame.transposed2D == False):
                    xval = ppms0[i]
                    yval = ppms1[i]
                else:
                    xval = ppms1[i]
                    yval = ppms0[i]
                z_list = self.main_frame.main_frame.uc2.ppm(peak["X_AXIS"])
                intensity_list = []
                for p in peak:
                    intensity_list.append(p[-1])
                for j,zval in enumerate(z_list):
                    names1.append(names[i]+'_'+str(j+1))
                    x.append(xval)
                    y.append(yval)
                    z.append(zval)
                    intensities.append(intensity_list[j])


            


            picked_peak_array = []
            for i, xval in enumerate(x):
                picked_peak_array.append([xval, y[i], z[i], intensities[i]])

            dictionary = {}
            dictionary["peak_name"] = []
            dictionary["shift1"] = []
            dictionary["shift2"] = []
            dictionary["shift3"] = []
            dictionary["intensity"] = []
            for i, peak in enumerate(picked_peak_array):
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

        self.current_peaklist_box.SetValue(peaklist_name)

        # Update the plot with the new peaklist
        self.main_frame.OnBoreSlider(wx.EVT_BUTTON)

        # Save the 3D peaklist
        self.OnSave(wx.EVT_BUTTON, peaklist_file=peaklist_name)

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

        self.AddToTable()

        self.current_peaklist_box.SetValue(last_directories_path)

        self.main_frame.OnBoreSlider(wx.EVT_BUTTON)

    def AddToTable(self):
        """
        Adding the peaklist just entered into the peaklist table
        """
        row_count = self.grid.GetNumberRows()
        if row_count > 0:
            self.grid.DeleteRows(0, row_count)
        peaklist = self.peak_list
        data = []

        def extract_number(s):
            match = re.match(r"(\d+)", s)
            return int(match.group(1)) if match else float("inf")

        # Pair each item with its original index
        indexed_arr = list(enumerate(self.peak_list_dictionary[peaklist]["peak_name"]))

        # Sort by number while keeping track of original indices
        sorted_indexed = sorted(indexed_arr, key=lambda x: extract_number(x[1]))

        # Extract sorted values and index mapping
        sorted_values = [val for _, val in sorted_indexed]
        index_mapping = {
            new_idx: old_idx for new_idx, (old_idx, _) in enumerate(sorted_indexed)
        }

        for i, peak_name in enumerate(self.peak_list_dictionary[peaklist]["peak_name"]):
            index = index_mapping[i]
            peak = self.peak_list_dictionary[peaklist]["peak_name"][index]
            shift1 = self.peak_list_dictionary[peaklist]["shift1"][index]
            shift2 = self.peak_list_dictionary[peaklist]["shift2"][index]
            shift3 = self.peak_list_dictionary[peaklist]["shift3"][index]
            intensity = self.peak_list_dictionary[peaklist]["intensity"][index]
            data.append([peak, "{:.5f}".format(shift1), "{:.5f}".format(shift2), "{:.5f}".format(shift3), "{:.5e}".format(intensity)])

        num_rows = self.grid.GetNumberRows()
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
                                    except:
                                        pass

            except:
                self.peaklist_error_message()
                return None

            if len(dictionary["peak_name"]) == 0 and new_peaklist == False:
                self.peaklist_error_message()
                return None

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
        Try to see if the chemical shifts of the peaks are within the 2D spectral range
        """
        ppms_0 = copy.deepcopy(dictionary["shift1"])
        ppms_1 = copy.deepcopy(dictionary["shift2"])
        ppms_2 = copy.deepcopy(dictionary["shift3"])


        shifts = [ppms_0, ppms_1, ppms_2]

        mean_0 = np.mean(ppms_0)
        mean_1 = np.mean(ppms_1)
        mean_2 = np.mean(ppms_2)


        # If there are two dimensions with the same ppm axis (e.g. (H)N(CA)NH), then the bore should be dimension 3

        # find out which chemical shift is the bore dimension
        if mean_0 > np.min(self.main_frame.main_frame.ppms_2) and mean_0 < np.max(
            self.main_frame.main_frame.ppms_2
        ):
            bore_shifts = 0
        elif mean_1 > np.min(self.main_frame.main_frame.ppms_2) and mean_1 < np.max(
            self.main_frame.main_frame.ppms_2
        ):
            bore_shifts = 1
        elif mean_2 > np.min(self.main_frame.main_frame.ppms_2) and mean_2 < np.max(
            self.main_frame.main_frame.ppms_2
        ):
            bore_shifts = 2
        else:
            dlg = wx.MessageDialog(
                self,
                "Chemical shifts in the peaklist for the bore dimension do not match any chemical shift axis. Try using a different peaklist",
                "Warning",
                wx.OK,
            )
            return None
        

        ppms_projection = []
        for i in range(3):
            if i == bore_shifts:
                continue
            else:
                ppms_projection.append(shifts[i])

        dictionary["shift3"] = shifts[bore_shifts]

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
                self.bore_xdim = 'shift2'

        return dictionary
    
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

    def turn_off_togglebuttons(self):
        # If any toggle buttons are on, turn them off
        if self.active_add == True:
            self.active_add = False
            self.add_peaks_button.SetValue(False)
            self.main_frame.fig_bore.canvas.mpl_disconnect(self.add_peak_connect)
        if self.active_move:
            self.active_move = False
            if self.active_select_peak:
                self.main_frame.fig_bore.canvas.mpl_disconnect(self.move_peak_connect)
            if self.active_select_peaks:
                self.main_frame.fig_bore.canvas.mpl_disconnect(self.move_peak_press)
                self.main_frame.fig_bore.canvas.mpl_disconnect(self.move_peak_motion)
                self.main_frame.fig_bore.canvas.mpl_disconnect(self.move_peak_release)
        if self.active_select_peak == True:
            self.select_peak_button.SetValue(False)
            self.active_select_peak = False
            self.selected_peakname = ""
            self.main_frame.fig_bore.canvas.mpl_disconnect(self.select_peak_connect)
            self.main_frame.OnBoreSlider(wx.EVT_BUTTON)

        self.update_mode_buttons()

    def OnLoadReferencePlane(self, event):
        """
        If a reference plane is already loaded, then setting reference_plane to False
        Otherwise, setting reference_plane to true and asking the user to select the reference plane
        peaklist to be loaded from a file dialog.
        """

        if self.reference_plane == True:
            self.reference_plane = False
            return
        
        self.reference_plane = True

        # Opening up a file window asking the user to select the 1D peak list - must be in the format of 1st column = peak_name, 2nd column = peak_position
        dlg = wx.FileDialog(self, "Select the 2D reference plane peak list", wildcard="", style=wx.FD_OPEN)
        dlg.SetDirectory(os.getcwd())
        if dlg.ShowModal() == wx.ID_OK:
            peaklist_file = dlg.GetPath()
        else:
            self.OnLoadReferencePlane(wx.EVT_TOGGLEBUTTON)
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
            self.OnLoadReferencePlane(wx.EVT_TOGGLEBUTTON)
            return
        
        self.reference_peaklist = peaklist
    
        

    
    

        


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
            self.main_frame.fig_bore.canvas.mpl_disconnect(self.add_peak_connect)
            return



        if self.active_select_peak:
            self.select_peak_button.SetValue(False)
            self.active_select_peak = False
            self.selected_peakname = ""
            self.main_frame.fig_bore.canvas.mpl_disconnect(self.select_peak_connect)
            self.main_frame.OnBoreSlider(wx.EVT_BUTTON)

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
                None, "Creating new peaklist", wildcard="*.list|*.txt", style=wx.FD_SAVE
            )
            dlg.SetDirectory(os.getcwd())
            if dlg.ShowModal() == wx.ID_OK:
                peaklist_file = dlg.GetPath()
                with open(peaklist_file, "w") as file:
                    pass

                self.AddPeaklist(peaklist_file, new_peaklist=True)

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
            current_peaklist = self.current_peaklist_box.GetValue()

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

            self.peak_list_dictionary[current_peaklist]["peak_name"].append(peakname)
            self.peak_list_dictionary[current_peaklist]["shift1"].append(x)
            self.peak_list_dictionary[current_peaklist]["shift2"].append(y)
            self.peak_list_dictionary[current_peaklist]["shift3"].append(0)
            self.peak_list_dictionary[current_peaklist]["intensity"].append(0)

            self.main_frame.OnBoreSlider(wx.EVT_BUTTON)
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

        # First need to disable other toggle buttons that are selected
        if self.active_add == True:
            self.active_add = False
            self.add_peaks_button.SetValue(False)
            self.main_frame.fig_bore.canvas.mpl_disconnect(self.add_peak_connect)

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
            if self.active_select_peak == True:
                self.main_frame.fig_bore.canvas.mpl_disconnect(self.move_peak_connect)
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
            self.main_frame.fig_bore.canvas.mpl_disconnect(self.select_peak_connect)
            evt_id = event.GetEventType()
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


        # If select peak, and there is a peak selected give a popout telling
        # the user to click where they want the peak to go
        if self.active_select_peak == True:
            self.move_peak_connect = self.main_frame.fig_bore.canvas.mpl_connect(
                "button_press_event", self.on_click_movepeak3d
            )


    def OnAddBorePeak(self, event):
        # This functionality will be added shortly
        dlg = wx.MessageDialog(
                    self,
                    "This feature is not yet implemented, but is planned to be added to a future release.",
                    "Not yet implemented",
                    wx.OK)
        dlg.ShowModal()
        dlg.Destroy()


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
            if self.active_select_peak == True:
                self.main_frame.fig_bore.canvas.mpl_disconnect(self.move_peak_connectz)
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
            self.main_frame.fig_bore.canvas.mpl_disconnect(self.select_peak_connect)
            evt_id = event.GetEventType()
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


        # If select peak, and there is a peak selected give a popout telling
        # the user to click where they want the peak to go
        if self.active_select_peak == True:
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

        if self.active_select_peak == True:
            if "N/A" not in self.selected_peak_indexes:
                if self.remove_peak == True:
                    count = 0
                    for peak_index in self.selected_peak_indexes:
                        del self.peak_list_dictionary[
                            self.current_peaklist_box.GetValue()
                        ]["peak_name"][peak_index]
                        del self.peak_list_dictionary[
                            self.current_peaklist_box.GetValue()
                        ]["shift1"][peak_index]
                        del self.peak_list_dictionary[
                            self.current_peaklist_box.GetValue()
                        ]["shift2"][peak_index]
                        del self.peak_list_dictionary[
                            self.current_peaklist_box.GetValue()
                        ]["shift3"][peak_index]

                        count += 1

                    self.remove_peak = False
                    self.selected_peak_indexes = ["N/A"]
                    self.selected_peakname = ""

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

        self.main_frame.ax_bore.set_xlim([xmin, xmax])
        self.main_frame.ax_bore.set_ylim([ymin, ymax])

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

    def OnSave(self, event, peaklist_file=''):
        """
        Provide a FileDialog where the user can chose the name for
        the peaklist.
        The peaklist will then be saved.
        """

        save_2d_plane = False

        if(peaklist_file == ''):
            dlg = wx.FileDialog(self, "Select the peak list", wildcard="", style=wx.FD_SAVE)
            dlg.SetDirectory(os.getcwd())
            if dlg.ShowModal() == wx.ID_OK:
                peaklist_file = dlg.GetPath()
            else:
                dlg.Destroy()
                return
            
            
            message = 'Would you like to save the full 3D peaklist (click yes), or would you like to save the 2D reference plane (click no)?'
            dlg = wx.MessageDialog(None, message, "Pick Peaks", wx.YES_NO)
            result=dlg.ShowModal()
            if(result == wx.ID_NO):
                save_2d_plane = True


        with open(peaklist_file, "w") as file:

            # Save all elements in the grid
            num_rows = self.grid.GetNumberRows()
            for i in range(num_rows):
                
                peak = self.grid.GetCellValue(i, 0)
                shift1 = self.grid.GetCellValue(i, 1)
                shift2 = self.grid.GetCellValue(i, 2)
                shift3 = self.grid.GetCellValue(i, 3)
                intensity = self.grid.GetCellValue(i, 4)
                if(save_2d_plane == False):
                    file.write(
                        "{} \t {} \t {} \t {} \t {}\n".format(peak, shift1, shift2, shift3, intensity)
                    )
                else:
                    file.write(
                        "{} \t {} \t {}\n".format(peak, shift1, shift2)
                    )



