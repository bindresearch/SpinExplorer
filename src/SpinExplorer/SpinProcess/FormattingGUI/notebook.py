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
from pathlib import Path
import shutil
import threading

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

from SpinExplorer.SpinConverter.SpinConverter import SpinConverter

from SpinExplorer.SpinProcess.StoringParameters.save_parameters import Save_json
from SpinExplorer.SpinProcess.StoringParameters.read_parameters import Read_json

from SpinExplorer.SpinView.SpinView import SpinView


from SpinExplorer.SpinExplorer_CL_tools.pulse_sequence_parsing import (
    PulseSequenceParser,
)

from appdirs import user_data_dir

appname = 'SpinExplorer'
appauthor = "James Eaton"
data_dir = user_data_dir(appname, appauthor)
os.makedirs(data_dir, exist_ok=True)


class NotebookProcess(wx.Notebook):
    def __init__(self, parent, nmr_data, popout=True):
        """
        Initialising a notebook which will hold the graphical interface
        for each complex dimension of the NMR spectrum. The user can
        toggle between each dimension. There is also a save processing
        parameters button and process (nmrpipe or nmrglue) buttons.
        """
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 0.7 * self.monitorWidth
        self.height = 0.75 * self.monitorHeight
        self.popout = popout
        self.parent = parent
        wx.Notebook.__init__(
            self,
            parent,
            id=wx.ID_ANY,
            style=wx.BK_DEFAULT,
            size=(self.width, self.height),
        )

        self.nmr_data = nmr_data

        # In developer mode any processing can be used whatever the data was
        # converted with. It can also be turned on before starting SpinExplorer
        # by setting the SPINEXPLORER_DEVELOPER environment variable
        self.developer_mode = os.environ.get("SPINEXPLORER_DEVELOPER", "") not in [
            "",
            "0",
            "false",
            "False",
        ]

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
        self.tabDim1.dimension_size.set_data_dimension(0)
        if self.nmr_data.dim == 2 and self.nmr_data.pseudo_axis == False:
            self.tabDim2 = IndirectDimensionFrame(
                self.parent, self, info_buttons, self.tabDim1, dimension=1
            )
            self.tabs.append(self.tabDim2)
            self.AddPage(
                self.tabDim2, "Dimension 2 (" + self.nmr_data.axislabels[1] + ")"
            )
            self.tabDim2.dimension_size.set_data_dimension(1)
        if self.nmr_data.dim == 3 and self.nmr_data.pseudo_axis == True:
            if self.nmr_data.index == 2:
                self.tabDim2 = IndirectDimensionFrame(
                    self.parent, self, info_buttons, self.tabDim1, 1
                )
                self.tabs.append(self.tabDim2)
                self.AddPage(
                    self.tabDim2, "Dimension 2 (" + self.nmr_data.axislabels[1] + ")"
                )
                self.tabDim2.dimension_size.set_data_dimension(1)
            else:
                self.tabDim2 = IndirectDimensionFrame(
                    self.parent, self, info_buttons, self.tabDim1, 2
                )
                self.tabs.append(self.tabDim2)
                self.AddPage(
                    self.tabDim2, "Dimension 2 (" + self.nmr_data.axislabels[2] + ")"
                )
                self.tabDim2.dimension_size.set_data_dimension(2)
        if self.nmr_data.dim == 3 and self.nmr_data.pseudo_axis == False:
            self.tabDim2 = IndirectDimensionFrame(
                self.parent, self, info_buttons, self.tabDim1, 1
            )
            self.tabs.append(self.tabDim2)
            self.AddPage(
                self.tabDim2, "Dimension 2 (" + self.nmr_data.axislabels[1] + ")"
            )
            self.tabDim2.dimension_size.set_data_dimension(1)
            self.tabDim3 = IndirectDimensionFrame(
                self.parent, self, info_buttons, self.tabDim1, 2
            )
            self.tabs.append(self.tabDim3)
            self.AddPage(
                self.tabDim3, "Dimension 3 (" + self.nmr_data.axislabels[2] + ")"
            )
            self.tabDim3.dimension_size.set_data_dimension(2)

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

        self.show_fid = wx.Button(parent, -1, "Show FID")
        self.show_fid.Bind(wx.EVT_BUTTON, self.on_show_fid)
        self.button_sizer.Add(self.show_fid, 0, wx.ALIGN_CENTER_VERTICAL)
        self.button_sizer.AddSpacer(10)


        self.saveprocessing = wx.Button(parent, -1, "Save Processing Parameters")
        self.saveprocessing.Bind(wx.EVT_BUTTON, self.on_save_processing)
        self.button_sizer.Add(self.saveprocessing, 0, wx.ALIGN_CENTER_VERTICAL)
        self.button_sizer.AddSpacer(10)
        self.saveprocessing_default = wx.Button(parent, -1, "Save Default Parameters")
        self.saveprocessing_default.Bind(wx.EVT_BUTTON, self.on_save_default_processing)
        self.button_sizer.Add(self.saveprocessing_default, 0, wx.ALIGN_CENTER_VERTICAL)
        self.button_sizer.AddSpacer(10)
        self.processing_button_pipe = wx.Button(parent, -1, "Process (nmrPipe)")
        self.processing_button_pipe.Bind(wx.EVT_BUTTON, self.on_process_pipe)
        self.button_sizer.Add(self.processing_button_pipe, 0, wx.ALIGN_CENTER_VERTICAL)
        self.button_sizer.AddSpacer(10)
        self.processing_button_glue = wx.Button(parent, -1, "Process (nmrglue)")
        self.processing_button_glue.Bind(wx.EVT_BUTTON, self.on_process_nmrglue)
        self.button_sizer.Add(self.processing_button_glue, 0, wx.ALIGN_CENTER_VERTICAL)
        self.button_sizer.AddSpacer(10)
        self.batch_processing_button_glue = wx.Button(parent, -1, "Batch Process (nmrglue)")
        self.batch_processing_button_glue.Bind(wx.EVT_BUTTON, self.on_batch_process_nmrglue)
        self.button_sizer.Add(self.batch_processing_button_glue, 0, wx.ALIGN_CENTER_VERTICAL)

        self.button_sizer.AddSpacer(10)
        self.developer_mode_checkbox = wx.CheckBox(parent, -1, "Developer mode")
        self.developer_mode_checkbox.SetValue(self.developer_mode)
        self.developer_mode_checkbox.SetToolTip(
            "Allow any processing and reconstruction whatever the data was "
            "converted with. Processing which does not match the conversion "
            "can fail or give incorrect spectra."
        )
        self.developer_mode_checkbox.Bind(wx.EVT_CHECKBOX, self.on_developer_mode)
        self.button_sizer.Add(
            self.developer_mode_checkbox, 0, wx.ALIGN_CENTER_VERTICAL
        )

        # Only the processing which matches the conversion can be used
        self.update_processing_buttons()

        self.parent.main_sizer.AddSpacer(20)
        self.parent.main_sizer.Add(self.button_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.parent.main_sizer.AddSpacer(10)

    def find_conversion_method(self) -> str:
        """
        Whether the fid being processed was converted using nmrglue or using
        nmrPipe. The nmrglue conversion records this in the comment of the
        NMRPipe header, so anything else has come from the nmrPipe conversion.
        """
        try:
            comment = str(self.nmr_data.dic["FDCOMMENT"])
        except (KeyError, AttributeError, TypeError):
            return "nmrPipe"

        if "nmrglue" in comment:
            return "nmrglue"

        return "nmrPipe"

    def find_developer_mode(self) -> bool:
        """
        Whether developer mode is on, where any processing can be used whatever
        the data was converted with.
        """
        return getattr(self, "developer_mode", False)

    def on_developer_mode(self, event):
        """
        The user has turned developer mode on or off. The processing which is
        available is updated to match.
        """
        self.developer_mode = self.developer_mode_checkbox.GetValue()

        self.update_processing_buttons()
        self.update_nus_options()

    def update_nus_options(self):
        """
        Update the NUS reconstruction options of every indirect dimension, as
        SMILE reconstruction is part of nmrPipe and so cannot be used with data
        converted using nmrglue.
        """
        for tab in self.tabs:
            linear_prediction = getattr(tab, "linear_prediction", None)
            update_smile = getattr(linear_prediction, "update_smile_option", None)
            if update_smile != None:
                update_smile()

    def find_smile_allowed(self) -> bool:
        """
        SMILE NUS reconstruction is part of nmrPipe, so it can only be used
        with data which was converted using nmrPipe.
        """
        if self.find_developer_mode() == True:
            return True

        return self.find_conversion_method() == "nmrPipe"

    def find_smile_message(self) -> str:
        """
        The message telling the user why SMILE reconstruction is not possible
        for this data.
        """
        return (
            "SMILE NUS reconstruction is part of nmrPipe, and this data was "
            "converted using nmrglue. To use SMILE, convert the data again "
            "using nmrPipe in SpinConverter. SpinExplorer IST NUS "
            "reconstruction can be used with this data."
        )

    def find_conversion_message(self, method: str) -> str:
        """
        The message telling the user why processing using a method is not
        possible for this data.
        """
        return (
            "This data was converted using {0}, so it can only be processed "
            "using {0}. To process it using {1}, convert the data again using "
            "{1} in SpinConverter.".format(self.find_conversion_method(), method)
        )

    def update_processing_buttons(self):
        """
        Data converted using nmrPipe has to be processed using nmrPipe and
        data converted using nmrglue has to be processed using nmrglue, so
        the button for the processing which does not match the conversion is
        greyed out.
        """
        conversion = self.find_conversion_method()

        buttons = [
            (self.processing_button_pipe, "nmrPipe"),
            (self.processing_button_glue, "nmrglue"),
        ]

        for button, method in buttons:
            if method == conversion:
                button.Enable(True)
                button.SetToolTip(
                    "Process the data using {}, which it was converted with.".format(
                        method
                    )
                )
            elif self.find_developer_mode() == True:
                button.Enable(True)
                button.SetToolTip(
                    "Developer mode: process the data using {}, which it was "
                    "not converted with.".format(method)
                )
            else:
                button.Enable(False)
                button.SetToolTip(self.find_conversion_message(method))

    def check_conversion_method(self, method: str) -> bool:
        """
        Check that the processing which has been asked for matches the way the
        data was converted, telling the user when it does not.
        """
        if self.find_conversion_method() == method:
            return True

        if self.find_developer_mode() == True:
            return True

        dlg = wx.MessageDialog(
            self,
            self.find_conversion_message(method),
            "Processing not possible",
            wx.OK | wx.ICON_INFORMATION,
        )
        self.Raise()
        self.SetFocus()
        dlg.ShowModal()
        dlg.Destroy()

        return False

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

    def on_show_fid(self, event):
        """
        Create a SpinView window containing the FID
        """
        fid_frame = SpinView(fid_data=[self.nmr_data])
    

    def on_process_pipe(self, event):
        """
        Process using nmrPipe. The code will first check if the user
        is selecting NUS processing. Various popouts will warn users
        about potential issues such as ensuring SMILE is installed
        or potential long reconstruction times if direct dimension
        data extraction is not selected.
        """
        if self.check_conversion_method("nmrPipe") == False:
            return

        try:
            if (self.tabDim2.linear_prediction.linear_prediction_radio_box_indirect.GetSelection()==2):
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

                if self.tabDim1.extraction.extraction_checkbox.GetValue() == False:
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


        self.sync_gui_values()

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
        if(self.popout==True):
            if self.check_conversion_method("nmrglue") == False:
                return

        try:
            if (self.tabDim2.linear_prediction.linear_prediction_radio_box_indirect.GetSelection()==2):
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

            if (self.tabDim2.linear_prediction.linear_prediction_radio_box_indirect.GetSelection()==3):
                if self.tabDim1.extraction.extraction_checkbox.GetValue() == False:
                    dlg = wx.MessageDialog(
                        self,
                        "No direct dimension data extraction is selected, SpinExplorer IST reconstruction may take a while. Consider extracting a region of the direct dimension before reconstruction. Do you want to continue or cancel?",
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

        self.sync_gui_values()

        checking = CheckingParameters(self, self.tabs)
        if checking.check_parameter_validity() == True:
            if(checking.check_nmrglue_fid(self.nmr_data)==True):
                self.on_save_processing(wx.EVT_BUTTON, popout_messages=self.popout)
                self.processing = ProcessNMRGlue(
                    self, self.tabs, self.nmr_data, interactive_phasing=False, popout_messages=self.popout
                )


    def on_batch_process_nmrglue(self, event):
            """
            Batch processing using nmrglue. First the data will be checked to see
            if SMILE non-uniform sampling reconstruction is present. If it
            is present then the processing will stop as SMILE reconstruction
            is part of nmrPipe and is not possible in nmrglue. In the future,
            NUS reconstruction using FID-Net will be added as a possibility
            to remove the requirement of SMILE from nmrPipe.


            This will try to apply the same conversion and processing for the current
            spectrum to all spectra in a user selected folder.

            """
            if self.check_conversion_method("nmrglue") == False:
                return
    
            try:
                if (self.tabDim2.linear_prediction.linear_prediction_radio_box_indirect.GetSelection()==2):
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
    
                if (self.tabDim2.linear_prediction.linear_prediction_radio_box_indirect.GetSelection()==3):
                    if self.tabDim1.extraction.extraction_checkbox.GetValue() == False:
                        dlg = wx.MessageDialog(
                            self,
                            "No direct dimension data extraction is selected, SpinExplorer IST reconstruction may take a while. Consider extracting a region of the direct dimension before reconstruction. Do you want to continue or cancel?",
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
    
            self.sync_gui_values()
    
            checking = CheckingParameters(self, self.tabs)
            if checking.check_parameter_validity() == True:
                if(checking.check_nmrglue_fid(self.nmr_data)==True):
                    # Asking the user which folder of data they wish to apply batch conversion and processing too
                    dlg = wx.DirDialog(None, "Choose folder containing data folders for batch processing. Processing will be applied to the current dataset and all data in the selected folder.", "", wx.DD_DIR_MUST_EXIST)

                    if dlg.ShowModal() != wx.ID_OK:
                        dlg2 = wx.MessageDialog(
                            self,
                            "No folder was selected. No processing will be applied. ",
                            "Warning",
                            wx.OK | wx.ICON_WARNING,
                        )
                        self.Raise()
                        self.SetFocus()
                        result = dlg2.ShowModal()
                        dlg2.Destroy()
                        dlg.Destroy()
                        return
                    selected_folder = dlg.GetPath()
                    # See how many directories are within the selected folder and if they have
                    folder = Path(selected_folder)
                    folders_of_data = [sub for sub in folder.iterdir() if sub.is_dir() and any((sub / n).is_file() for n in ['fid', 'ser'])]

                    if(len(folders_of_data)==0):
                        dlg = wx.MessageDialog(
                            self,
                            "There are no folders within the selected folder that contain fid and ser files to perform batch processing on. Make sure that within each data folder within the selected folder that there is an fid or ser file. Exiting.",
                            "Warning",
                            wx.OK | wx.ICON_WARNING,
                        )
                        self.Raise()
                        self.SetFocus()
                        result = dlg.ShowModal()

                    
                    self.on_save_processing(wx.EVT_BUTTON, popout_messages = False)
                    processing = ProcessNMRGlue(
                        self, self.tabs, self.nmr_data, interactive_phasing=False, batch_main = True, popout_messages=False
                    )

                    # If successful then can move onto the batch processing
                    continue_batch = False
                    if(processing.success == True):
                        dlg = wx.MessageDialog(
                            self,
                            "Processing worked successfully on the current open dataset. Proceeding to batch processing will copy the conversion/processing parameters from the current dataset to all folders in the selected data folder and perform batch conversion and processing on them. This will overwrite any previously saved parameters or processed data in each data folder. Would you like to continue?",
                            "Continue batch processing",
                            wx.YES_NO,
                        )
                        self.Raise()
                        self.SetFocus()
                        result = dlg.ShowModal()
                        if(result == wx.ID_YES):
                            continue_batch = True
                    else:
                        dlg = wx.MessageDialog(
                            self,
                            "Processing was not completed successfully on the current dataset, batch processing will not continue. Please ensure the current dataset can be processed successfully before trying batch processing.",
                            "Warning",
                            wx.OK | wx.ICON_WARNING
                        )
                        self.Raise()
                        self.SetFocus()
                        result = dlg.ShowModal()


                    original_folder = os.getcwd()

                    list_of_folders = []
                    for folder in folders_of_data:
                        if(str(folder)!=str(original_folder)):
                            list_of_folders.append(str(folder))

                    number_of_folders_to_process = len(list_of_folders)

                    
                        

                    if(continue_batch==True):
                        count = 0
                        success_count = 0
                        failed_message = ""
                        dlg = LogProgressDialog(self, maximum=number_of_folders_to_process)
                        dlg.Show()
                        dlg.Raise()
                                               
                        for data_folder in folders_of_data:
                            if(str(data_folder)==str(original_folder)):
                                continue
                            else:
                                dlg.update(count, f"Data folder {count} of {number_of_folders_to_process}", f"Beginning processing of {data_folder}")
                                if dlg.cancel_event.is_set():
                                    os.chdir(original_folder)
                                    break
                                if not dlg:
                                    break
                                count += 1
                                shutil.copy('parameters.json', os.path.join(data_folder, 'parameters.json'))

                                # Open a hidden window of SpinConverter with the parameters from the json file loaded
                                os.chdir(data_folder)

                                

                                converter = SpinConverter(explorer=True, popout=False)
                                converter.Hide()

                                converter.on_convert_glue(wx.EVT_BUTTON)

                                if dlg.cancel_event.is_set():
                                    os.chdir(original_folder)
                                    break
                                if not dlg:
                                    break

                                converter_success = converter.glue_conversion.success
                                converter.Destroy()

                                from SpinExplorer.SpinProcess.SpinProcess import SpinProcess
                                processer = SpinProcess(explorer=True, popout=False)
                                processer.Hide()
                                processer.notebook.on_process_nmrglue(wx.EVT_BUTTON)

                                if dlg.cancel_event.is_set():
                                    os.chdir(original_folder)
                                    break
                                if not dlg:
                                    break

                                success = processer.notebook.processing.success
                                processer.Destroy()
                                if(success):
                                    success_count+=1
                                    dlg.update(count, f"Data folder {count} of {number_of_folders_to_process}", f"Processed spectrum in folder {data_folder}")
                        
                                else:
                                    dlg.update(count, f"Data folder {count} of {number_of_folders_to_process}", f"Failed to process data in {data_folder}")

                                if dlg.cancel_event.is_set():
                                    os.chdir(original_folder)
                                    break
                                if not dlg:
                                    break
                                        

                                os.chdir(original_folder)


                        if(count==success_count):
                            dlg.update(count, "Processed all spectra successfully. Close window to continue.")
                        else:
                            dlg.update(count, "Processed " + str(success_count) + " of " + str(number_of_folders_to_process) + " spectra successfully. See traceback for details. Close window to continue.")


    

    def sync_gui_values(self):
        """
        Making sure that the values currently shown in the graphical interface
        have been stored in the processing parameter variables. The apodization
        and NUS textcontrols only update their stored values when the user
        presses enter, so this is called before the parameters are saved or
        used for processing.
        """
        for tab in self.tabs:
            for component in [tab.apodization, tab.linear_prediction]:
                update_values = getattr(
                    component, "update_stored_values_from_gui", None
                )
                if update_values != None:
                    update_values()

    def on_save_default_processing(self, event):
        """
        Saving the current parameters in the SpinProcess graphical interface
        into parameters.json and saving these as the default processing for 
        this pulse sequence.
        """
        self.sync_gui_values()
        # save the current parameters to parameters.json
        save = Save_json(self, self.nmr_data, self.tabs, popout_messages=True)

        # Get the current pulse sequence
        self.auto_sequence = None
        try:
            self.auto_sequence = PulseSequenceParser().parse()
    
            # copy parameters.json to the location of saved default processing
            self.package_directory = user_data_dir(appname, appauthor)
            os.makedirs(self.package_directory, exist_ok=True)
            self.saved_processing_directory = os.path.join(self.package_directory, "saved_processing")
            os.makedirs(self.saved_processing_directory, exist_ok=True)
            # search through all saved files to see if the pulseprogram exists
            saved_default_path = os.path.join(self.saved_processing_directory, self.auto_sequence)
            shutil.copy(os.path.join(str(os.getcwd()), 'parameters.json'), saved_default_path)

            dlg = wx.MessageDialog(
                self,
                "Default processing parameters have been successfully saved for the current pulse sequence " + str(self.auto_sequence),
                "Save default parameters complete",
                wx.OK
            )
            self.Raise()
            self.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()

        except:
            dlg = wx.MessageDialog(
                self,
                "Unable to save the current parameters as default parameters.",
                "Error",
                wx.OK | wx.ICON_WARNING
            )
            self.Raise()
            self.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
        

    def on_save_processing(self, event, popout_messages=True):
        """
        Saving the current parameters in the SpinProcess graphical interface
        into parameters.json
        """
        self.sync_gui_values()
        save = Save_json(self, self.nmr_data, self.tabs, popout_messages)

    def on_read_processing(self):
        """
        Read a previous parameters.json file and load these values into
        the graphical interface for SpinProcess.
        """

        self.read = Read_json(self, self.nmr_data, self.tabs, popout=self.popout)

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





class LogProgressDialog(wx.Dialog):
    def __init__(self, parent, title="Batch Process", maximum=100):
        super().__init__(parent, title=title, size=(600, 350),
                         style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)

        self.status = wx.StaticText(self, label="Starting batch process...")
        self.gauge = wx.Gauge(self, range=maximum)
        self.log = wx.TextCtrl(
            self, style=wx.TE_MULTILINE | wx.TE_READONLY | wx.TE_RICH2
        )
        self.btn = wx.Button(self, wx.ID_CANCEL, "Cancel")

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.status, 0, wx.ALL | wx.EXPAND, 10)
        sizer.Add(self.gauge, 0, wx.LEFT | wx.RIGHT | wx.EXPAND, 10)
        sizer.Add(self.log, 1, wx.ALL | wx.EXPAND, 10)
        sizer.Add(self.btn, 0, wx.ALL | wx.ALIGN_RIGHT, 10)
        self.SetSizer(sizer)

        self.cancel_event = threading.Event()
        self.btn.Bind(wx.EVT_BUTTON, self.on_cancel)
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self.Centre()

    def on_cancel(self, event):
        if self.cancel_event.is_set() or self.btn.GetLabel() == "Close":
            self.EndModal(wx.ID_OK)
        else:
            self.cancel_event.set()
            self.btn.Disable()
            self.append("Cancelled (close the window to continue)")

    # GUI-thread methods (call via wx.CallAfter from workers)
    def update(self, value, status=None, line=None):
        self.gauge.SetValue(value)
        if status is not None:
            self.status.SetLabel(status)
        if line is not None:
            self.append(line)

    def append(self, line):
        self.log.AppendText(line + "\n")   # auto-scrolls to the end

    def done(self, message="Done"):
        self.status.SetLabel(message)
        self.btn.Enable()
        self.btn.SetLabel("Close")

    def OnClose(self, event):
        self.Destroy()