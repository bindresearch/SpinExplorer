#!/usr/bin/env python3

"""MIT License

Copyright (c) 2025 James Eaton, Andrew Baldwin (University of Oxford)
              2025, Bind Research

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE."""

import wx
import os

from SpinExplorer.SpinProcess.FormattingGUI.ProcessingComponents.info_buttons import (
    InfoButtons,
)
from SpinExplorer.SpinProcess.Processing.process_nmrglue import (
    ProcessNMRGlue,
)
from SpinExplorer.SpinProcess.Processing.process_pipe import (
    ProcessNMRPipe,
)
from SpinExplorer.SpinProcess.FormattingGUI.frames import (
    DirectDimensionFrame,
    IndirectDimensionFrame,
)
from SpinExplorer.SpinProcess.Processing.checking_parameters import CheckingParameters

from SpinExplorer.SpinProcess.StoringParameters.save_parameters import Save_json
from SpinExplorer.SpinProcess.StoringParameters.read_parameters import Read_json


class NotebookProcess(wx.Notebook):
    def __init__(self, parent, nmr_data):
        """
        Initialising a notebook which will hold the graphical interface
        for each complex dimension of the NMR spectrum. The user can
        toggle between each dimension. There is also a save processing
        parameters button and process (nmrpipe or nmrglue) buttons.
        """
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

        self.nmr_data = nmr_data

        # Information buttons class
        info_buttons = InfoButtons(self)

        # Adding a processing frame for each non-real dimension
        self.add_frames(info_buttons)

    def add_frames(self, info_buttons):
        """
        Adding the processing frames to the notebook for each complex dimension.
        The direct dimension has a DirectDimensionFrame class structure
        whereas the indirect dimension has a
        IndirectDimensionFrame structure.

        Currently this is supported for 3 complex dimensions. However,
        additional IndirectDimensionFrame classes could be added to support
        additional dimensions in the future (though the processing classes
        would also need to be updated for this).
        """
        self.tabs = []
        self.tabDim1 = DirectDimensionFrame(self.parent, self, info_buttons)
        self.tabs.append(self.tabDim1)
        self.AddPage(self.tabDim1, "Dimension 1 (" + self.nmr_data.axislabels[0] + ")")
        if self.nmr_data.dim == 2 and self.nmr_data.pseudo_axis == False:
            self.tabDim2 = IndirectDimensionFrame(
                self.parent, self, info_buttons, self.tabDim1
            )
            self.tabs.append(self.tabDim2)
            self.AddPage(
                self.tabDim2, "Dimension 2 (" + self.nmr_data.axislabels[1] + ")"
            )
        if self.nmr_data.dim == 3 and self.nmr_data.pseudo_axis == True:
            if self.nmr_data.index == 2:
                self.tabDim2 = IndirectDimensionFrame(
                    self.parent, self, info_buttons, self.tabDim1
                )
                self.tabs.append(self.tabDim2)
                self.AddPage(
                    self.tabDim2, "Dimension 2 (" + self.nmr_data.axislabels[1] + ")"
                )
            else:
                self.tabDim2 = IndirectDimensionFrame(
                    self.parent, self, info_buttons, self.tabDim1
                )
                self.tabs.append(self.tabDim2)
                self.AddPage(
                    self.tabDim2, "Dimension 2 (" + self.nmr_data.axislabels[2] + ")"
                )
        if self.nmr_data.dim == 3 and self.nmr_data.pseudo_axis == False:
            self.tabDim2 = IndirectDimensionFrame(
                self.parent, self, info_buttons, self.tabDim1
            )
            self.tabs.append(self.tabDim2)
            self.AddPage(
                self.tabDim2, "Dimension 2 (" + self.nmr_data.axislabels[1] + ")"
            )
            self.tabDim3 = IndirectDimensionFrame(
                self.parent, self, info_buttons, self.tabDim1
            )
            self.tabs.append(self.tabDim3)
            self.AddPage(
                self.tabDim3, "Dimension 3 (" + self.nmr_data.axislabels[2] + ")"
            )

        # Setting the fourier transform modes to the guessed values
        self.add_ft_mode_guess()

        # read previously saved parameters
        self.on_read_processing()

    def create_buttons(self, parent):
        """
        Creating a button to save the current processing parameters
        as well as buttons for processing using either nmrPipe or
        nmrglue.
        """
        # Have a button for make nmrproc.com file, show nmrproc.com file and run processing
        self.button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.saveprocessing = wx.Button(parent, -1, "Save Processing Parameters")
        self.saveprocessing.Bind(wx.EVT_BUTTON, self.on_save_processing)
        self.button_sizer.Add(self.saveprocessing, 0, wx.ALIGN_CENTER_VERTICAL)
        self.button_sizer.AddSpacer(10)
        self.processing_button_pipe = wx.Button(parent, -1, "Process (nmrPipe)")
        self.processing_button_pipe.Bind(wx.EVT_BUTTON, self.on_process_pipe)
        self.button_sizer.Add(self.processing_button_pipe, 0, wx.ALIGN_CENTER_VERTICAL)
        self.button_sizer.AddSpacer(10)
        self.processing_button_glue = wx.Button(parent, -1, "Process (nmrglue)")
        self.processing_button_glue.Bind(wx.EVT_BUTTON, self.on_process_nmrglue)
        self.button_sizer.Add(self.processing_button_glue, 0, wx.ALIGN_CENTER_VERTICAL)

        self.parent.main_sizer.AddSpacer(20)
        self.parent.main_sizer.Add(self.button_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.parent.main_sizer.AddSpacer(10)

    def change_to_path(self):
        if self.parent.path != "":
            os.chdir(self.parent.path)
        else:
            if self.parent.original_frame != None:
                if self.parent.original_frame.parent.path != None:
                    try:
                        os.chdir(self.parent.original_frame.parent.path)
                    except:
                        pass
            if self.parent.file_parser == True:
                os.chdir(self.parent.path)

    def change_to_path_run(self):
        if self.parent.path != "":
            os.chdir(self.parent.path)
        else:
            if self.parent.original_frame != None:
                try:
                    self.parent.original_frame.Disable()
                except:
                    pass
                if self.parent.original_frame.parent.path != "":
                    os.chdir(self.parent.original_frame.parent.path)
            if self.parent.file_parser == True:
                os.chdir(self.parent.path)

    def change_to_cwd(self):
        if self.parent.cwd != "":
            os.chdir(self.parent.cwd)
        else:
            # Change path if using unidecFile parser
            if self.parent.original_frame != None:
                try:
                    if self.parent.original_frame.parent.cwd != None:
                        os.chdir(self.parent.original_frame.parent.cwd)
                except:
                    pass

            if self.parent.file_parser == True:
                # Change the path to the path of the original file
                os.chdir(self.parent.cwd)

    def on_process_pipe(self, event):
        """
        Process using nmrPipe. The code will first check if the user
        is selecting NUS processing. Various popouts will warn users
        about potential issues such as ensuring SMILE is installed
        or potential long reconstruction times if direct dimension
        data extraction is not selected.
        """
        try:
            if self.tabDim2.linear_prediction_radio_box_dim2.GetSelection() == 2:
                # SMILE processing is selected, asking the user to confirm SMILE is installed as part of nmrPipe
                dlg = wx.MessageDialog(
                    self,
                    "SMILE processing is selected. Ensure that SMILE is installed as part of nmrPipe",
                    "Warning",
                    wx.OK | wx.CANCEL | wx.ICON_WARNING,
                )
                self.Raise()
                self.SetFocus()
                result = dlg.ShowModal()
                if result == wx.ID_CANCEL:
                    self.change_to_cwd()
                    return

                if self.tabDim1.extraction_checkbox.GetValue() == False:
                    dlg = wx.MessageDialog(
                        self,
                        "No direct dimension data extraction is selected, SMILE reconstruction may take a while. Consider extracting a region of the direct dimension before reconstruction. Do you want to continue or cancel?",
                        "Warning",
                        wx.OK | wx.CANCEL | wx.ICON_WARNING,
                    )
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    if result == wx.ID_CANCEL:
                        self.change_to_cwd()
                        return
        except:
            pass

        checking = CheckingParameters(self, self.tabs)
        continue_processing = checking.check_parameter_validity()
        if continue_processing == True:
            if(checking.check_nmrpipe_fid(self.nmr_data)==True):
                self.on_save_processing(wx.EVT_BUTTON)
                processing = ProcessNMRPipe(self, self.tabs, self.nmr_data)

    def on_process_nmrglue(self, event):
        """
        Processing using nmrglue. First the data will be checked to see
        if SMILE non-uniform sampling reconstruction is present. If it
        is present then the processing will stop as SMILE reconstruction
        is part of nmrPipe and is not possible in nmrglue. In the future,
        NUS reconstruction using FID-Net will be added as a possibility
        to remove the requirement of SMILE from nmrPipe.
        """
        try:
            if self.tabDim2.linear_prediction_radio_box_dim2.GetSelection() == 2:
                # SMILE processing is selected, asking the user to confirm SMILE is installed as part of nmrPipe
                dlg = wx.MessageDialog(
                    self,
                    "NUS reconstruction using SMILE is selected. This is not possible using nmrglue processing. Stopping processing.",
                    "Warning",
                    wx.OK | wx.ICON_WARNING,
                )
                self.Raise()
                self.SetFocus()
                result = dlg.ShowModal()
                self.change_to_cwd()
                return
        except:
            pass

        checking = CheckingParameters(self, self.tabs)
        if checking.check_parameter_validity() == True:
            if(checking.check_nmrglue_fid(self.nmr_data)==True):
                self.on_save_processing(wx.EVT_BUTTON)
                processing = ProcessNMRGlue(
                    self, self.tabs, self.nmr_data, interactive_phasing=False
                )


    

    def on_save_processing(self, event):
        """
        Saving the current parameters in the SpinProcess graphical interface
        into parameters.json
        """
        save = Save_json(self, self.nmr_data, self.tabs)

    def on_read_processing(self):
        """
        Read a previous parameters.json file and load these values into
        the graphical interface for SpinProcess.
        """

        self.read = Read_json(self, self.nmr_data, self.tabs)

    def add_ft_mode_guess(self):
        """
        Setting the fourier transform modes for the indirect dimensions
        to the values guessed from the conversion acquisition modes if
        there are more than 2 complex dimensions.
        """
        dimensions = len(self.nmr_data.data.shape)
        if dimensions > 1:
            if dimensions == 2 and self.nmr_data.pseudo_axis == True:
                return
            else:
                indirect_modes = self.nmr_data.ft_options[1:]
                for i, mode in enumerate(indirect_modes):
                    try:
                        self.tabs[i + 1].fourier_transform.ft_method_selection = mode
                    except:
                        pass
