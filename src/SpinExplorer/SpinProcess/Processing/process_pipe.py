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
import subprocess
import nmrglue as ng
import shutil

# Import relevant SpinExplorer modules
from SpinExplorer.SpinProcess.Processing.PipeProcessing.write_nmrpipe_processing import (
    WriteNMRPipeProcessing,
)


class ProcessNMRPipe:
    def __init__(self, notebook, dimension_tabs, nmr_data):
        """
        This class will first create an nmrproc.com script based
        on the current parameters in SpinProcess. It also checks
        that each parameter is valid.
        The nmrproc.com script is ran using subprocess.
        """
        self.notebook = notebook
        self.dimension_tabs = dimension_tabs
        self.nmr_data = nmr_data

        make_nmrproc = self.on_make_nmrproc_com()
        if make_nmrproc == False:
            """
            If the user says they don't want to overwrite
            the nmrproc.com file, return without processing
            """
            return
        else:
            self.on_run_processing_nmrpipe()

    def on_make_nmrproc_com(self):
        """
        Making the nmrproc.com file based on the current
        SpinProcess parameters.
        """
        # Checking if an nmrproc.com file exits
        # return if it exists and the user doesn't want to overwrite it
        if self.checking_for_nmrproc(on_make=True) == False:
            return False

        # Initialise class which writes the nmrPipe processing lines
        write_nmrpipe = WriteNMRPipeProcessing(
            self.notebook, self.dimension_tabs, self.nmr_data
        )

        # Checking whether processing is required for second/third dimensions
        include_dim2, include_dim3 = self.checking_dimensions()


        # Create the nmrproc.com file
        nmrproc_com = open("nmrproc.com", "w")
        nmrproc_com = self.write_initial_lines(nmrproc_com)

        # Adding the processing lines for the direct dimension
        self.add_direct_dimension(nmrproc_com, write_nmrpipe)

        # Adding processing lines for the first complex indirect dimension
        if include_dim2:
            check_nus = self.check_nus(1)
            if check_nus != [0]:
                write_nmrpipe.apply_NUS_SMILE(nmrproc_com, check_nus)
                nmrproc_com = write_nmrpipe.transpose_line(nmrproc_com)

            else:
                if self.nmr_data.pseudo_axis == False:
                    nmrproc_com = write_nmrpipe.transpose_line(nmrproc_com)
                elif self.nmr_data.pseudo_axis == True:
                    if self.nmr_data.index == 2:
                        nmrproc_com = write_nmrpipe.transpose_line(nmrproc_com)
                    elif self.nmr_data.index == 1:
                        # If the pseudo axis is the central axis then need to move the third axis
                        nmrproc_com = write_nmrpipe.transpose_line(nmrproc_com)
                        nmrproc_com = write_nmrpipe.zero_transpose_line(nmrproc_com)
            
            if(self.check_lp(1)==1):
                write_nmrpipe.linear_prediction(1, nmrproc_com, indirect=True)
            self.add_indirect_dimension(nmrproc_com, 1, write_nmrpipe)

        # Adding processing lines for the second complex indirect dimension
        if include_dim3:
            nmrproc_com = write_nmrpipe.zero_transpose_line(nmrproc_com)
            if check_nus == [2]:
                write_nmrpipe.apply_NUS_SMILE(nmrproc_com, check_nus)
            if(self.check_lp(2)==1):
                write_nmrpipe.linear_prediction(2, nmrproc_com, indirect=True)
            self.add_indirect_dimension(nmrproc_com, 2, write_nmrpipe)

        # Adding the output line to the nmrproc.com script
        self.write_output_line(nmrproc_com)

    def on_run_processing_nmrpipe(self):
        """
        Checking that the nmrproc.com script exists and that
        it has excecutable permissions.
        """

        # Check to see if the nmrproc.com file exists
        checking = self.checking_for_nmrproc(on_make=False)
        if checking == False:
            return

        # Finding the correct name output for the processed data
        nmrfile = self.find_nmrfile_output_name(
            self.nmr_data.dim, self.nmr_data.pseudo_axis
        )

        # Add execute permissions to the nmrproc.com file
        os.system("chmod +x nmrproc.com")

        checking = self.checking_processed_data(nmrfile)
        if checking == False:
            return

        # Run the nmrproc.com file
        command = "csh nmrproc.com"

        self.nmrfile = nmrfile

        self.subprocess_frame = nmrPipeSubprocess(self)


    def update_comment(self):
        # Set the comment to pipe so that the fact nmrpipe processing was used is noted in the processed spectrum header
        if self.nmr_data.pseudo_axis == True:
            # Read the processed data and update the FDCOMMENT
            dic, data = ng.pipe.read(self.nmrfile)
            # Adding the fact that there is a pseudo axis to the FDCOMMENT
            dic['FDCOMMENT'] += 'pseudo'
            ng.pipe.write(self.nmrfile, dic, data, overwrite=True)
        self.final_changes(self.nmrfile)
        self.completion_notification(self.nmrfile)
        if self.notebook.parent.original_frame != None:
            self.notebook.parent.Destroy()

    def completion_notification(self, nmrfile: str):
        """
        This function provides the user with a message saying
        that the NMR processing completed successfully and
        has the name nmrfile.
        """
        dlg = wx.MessageDialog(
            self.notebook,
            "NMR processing is complete, data was saved as {}".format(nmrfile),
            "Completion",
            wx.OK,
        )
        self.notebook.Raise()
        self.notebook.SetFocus()
        result = dlg.ShowModal()

    def checking_processed_data(self, nmrfile: str) -> bool:
        """
        Function to test if there is currently processed data
        in the current directory and asking if the user wants
        to overwrite this.
        """
        # Check to see if nmrfile file already exists, if it does ask the user if they want to overwrite it
        if os.path.exists(nmrfile):
            dlg = wx.MessageDialog(
                self.notebook,
                "The {} file already exists. Do you want to overwrite it?".format(
                    nmrfile
                ),
                "Warning",
                wx.YES_NO | wx.ICON_WARNING,
            )
            self.notebook.Raise()
            self.notebook.SetFocus()
            result = dlg.ShowModal()
            if result == wx.ID_NO:
                if self.notebook.parent.original_frame != None:
                    self.notebook.parent.original_frame.Disable()
                    if self.notebook.parent.original_frame.parent.cwd != "":
                        os.chdir(self.notebook.parent.original_frame.parent.cwd)
                return False

        return True

    def find_nmrfile_output_name(self, dim: int, pseudo_axis: bool) -> str:
        """
        Determining the correct name for the nmrfile based on
        the dimensions and if the spectrum has a pseudo axis
        or not.
        """

        if dim == 1:
            nmrfile = "test.ft"
        elif dim == 2 and pseudo_axis == False:
            nmrfile = "test.ft2"
        elif dim == 2 and pseudo_axis == True:
            nmrfile = "test.ft"
        elif dim == 3 and pseudo_axis == True:
            nmrfile = "test.ft2"
        elif dim == 3 and pseudo_axis == False:
            nmrfile = "test.ft3"
        else:
            nmrfile = "test.ft"

        return nmrfile

    def final_changes(self, nmrfile: str):
        """
        Update the current path and also spawn a new SpinExplorer
        if are using the Reprocess SpinExplorer button.
        Also check for the processed nmrfile.
        """

        if self.notebook.parent.original_frame != None:
            if self.notebook.parent.original_frame.parent.cwd != "":
                os.chdir(self.notebook.parent.original_frame.parent.cwd)

        if self.notebook.parent.file_parser == True:
            os.chdir(self.notebook.parent.cwd)

        original_frame = []
        if self.notebook.parent.original_frame != None:
            original_frame = True
            self.notebook.parent.original_frame.Enable()
            path = self.notebook.parent.original_frame.parent.path
            cwd = self.notebook.parent.original_frame.parent.cwd
            self.notebook.parent.original_frame.parent.reprocess = True
            # self.notebook.parent.original_frame.parent.DestroyChildren()
            # self.subprocess_frame.Destroy()
            self.notebook.parent.original_frame.parent.Destroy()
            # self.notebook.parent.original_frame.parent.Hide()
            if self.notebook.parent.original_frame.parent.path != "":
                os.chdir(self.notebook.parent.original_frame.parent.path)
            from SpinExplorer.SpinView.SpinView import SpinView

            wx.CallAfter(self.update_spinview, (path, cwd))


    def update_spinview(self, a):
        from SpinExplorer.SpinView.SpinView import SpinView
        app = SpinView()
        if self.notebook.parent.original_frame.parent.cwd != "":
            path = a[0]
            cwd = a[1]
            app.path = path
            app.cwd = cwd

        # self.check_for_processed_file(nmrfile, app, original_frame)

    def check_for_processed_file(self, nmrfile: str, app, original_frame: bool):
        """
        Check to see if the output nmrfile file exists
        """
        if os.path.exists(nmrfile) == False:
            if original_frame == True:
                if app.cwd != "":
                    os.chdir(app.cwd)
            elif self.notebook.cwd != "":
                os.chdir(self.notebook.cwd)
            message = "The processing spectrum file ({}) file cannot be found in the current directory. Processing unsuccessful. Ensure that nmrPipe has been downloaded and added to the path.".format(
                nmrfile
            )
            dlg = wx.MessageDialog(
                self.notebook, message, "Warning", wx.OK | wx.ICON_WARNING
            )
            result = dlg.ShowModal()

            return
        else:
            if original_frame == True:
                if app.cwd != "":
                    os.chdir(app.cwd)
            elif self.notebook.cwd != "":
                os.chdir(self.notebook.cwd)
            message = "Processing successful. The processed spectrum file ({}) has been created in the current directory.".format(
                nmrfile
            )
            dlg = wx.MessageDialog(
                self.notebook, message, "Success", wx.OK | wx.ICON_INFORMATION
            )
            result = dlg.ShowModal()
            return

    def checking_for_nmrproc(self, on_make=False) -> bool:
        """
        Checking for nmrproc.com file.
        on_make = True (Checking before making nmrproc.com)
        on_make = False (Checking before running nmrproc.com)
        """
        if on_make == True:
            # Check to see if the nmrproc.com file already exists, if it does ask the user if they want to overwrite it
            if os.path.exists("./nmrproc.com"):
                dlg = wx.MessageDialog(
                    self.notebook,
                    "The nmrproc.com file already exists. Do you want to overwrite it?",
                    "Warning",
                    wx.YES_NO | wx.ICON_WARNING,
                )
                self.notebook.Raise()
                self.notebook.SetFocus()
                result = dlg.ShowModal()
                if result == wx.ID_NO:
                    return False
            return True
        else:
            if os.path.exists("./nmrproc.com") == False:
                dlg = wx.MessageDialog(
                    self.notebook,
                    "The nmrproc.com file cannot be found in the current directory",
                    "Warning",
                    wx.OK | wx.ICON_WARNING,
                )
                self.notebook.Raise()
                self.notebook.SetFocus()
                result = dlg.ShowModal()
                if self.notebook.parent.original_frame != None:
                    try:
                        self.notebook.parent.original_frame.Disable()
                    except:
                        pass
                    if self.notebook.parent.original_frame.parent.cwd != "":
                        os.chdir(self.notebook.parent.original_frame.parent.cwd)
                return False
            return True

    def checking_dimensions(self) -> list[bool]:
        """
        Based on the number of SpinProcess tabs, determine whether
        second/third dimension processing is required.
        """
        include_dim2 = False
        include_dim3 = False
        if self.nmr_data.dim == 2 and self.nmr_data.pseudo_axis == False:
            include_dim2 = True
        elif self.nmr_data.dim == 3:
            include_dim2 = True
        else:
            include_dim2 = False
        if self.nmr_data.dim == 3 and self.nmr_data.pseudo_axis == False:
            include_dim3 = True
        else:
            include_dim3 = False

        return include_dim2, include_dim3

    def write_initial_lines(self, nmrproc_com):
        """
        Writing the initial nmrproc.com lines
        """
        nmrproc_com.write("#!/bin/csh\n\n")
        if "/" not in self.nmr_data.fid_file:
            nmrproc_com.write("nmrPipe -in " + self.nmr_data.fid_file + " -verb \\\n")
        else:
            nmrproc_com.write("xyz2pipe -in ./fids/test%03d.fid -verb \\\n")

        return nmrproc_com

    def add_direct_dimension(self, nmrproc_com, write_nmrpipe):
        """
        Adding the relevent processing lines for the direct
        dimension to the nmrproc.com file.
        """
        dimension_tab = self.dimension_tabs[0]

        # Check to see if the solvent suppression checkbox is checked
        if (
            dimension_tab.solvent_suppression.solvent_suppression_checkbox.GetValue()
            == True
        ):
            nmrproc_com = write_nmrpipe.solvent_suppression(0, nmrproc_com)

        # Check to see if the linear prediction checkbox is checked
        if (
            dimension_tab.linear_prediction.linear_prediction_checkbox.GetValue()
            == True
        ):
            nmrproc_com = write_nmrpipe.linear_prediction(0, nmrproc_com)

        # Check to see if the apodization checkbox is checked
        if dimension_tab.apodization.apodization_checkbox.GetValue() == True:
            nmrproc_com = write_nmrpipe.apodization(0, nmrproc_com)

        # Check to see if the zero filling checkbox is checked
        if dimension_tab.zero_filling.zero_filling_checkbox.GetValue() == True:
            nmrproc_com = write_nmrpipe.zero_filling(0, nmrproc_com)

        # Check to see if the fourier transform checkbox is checked
        if (
            dimension_tab.fourier_transform.fourier_transform_checkbox.GetValue()
            == True
        ):
            nmrproc_com = write_nmrpipe.fourier_transform(0, nmrproc_com)

        # Check to see if the phase correction checkbox is checked
        if dimension_tab.phasing.magnitude_mode_checkbox.GetValue() == True:
            nmrproc_com = write_nmrpipe.magnitude_mode(nmrproc_com)
        elif dimension_tab.phasing.phase_correction_checkbox.GetValue() == True:
            nmrproc_com = write_nmrpipe.phasing(0, nmrproc_com)

        # Check to see if the extraction checkbox is checked
        if dimension_tab.extraction.extraction_checkbox.GetValue() == True:
            nmrproc_com = write_nmrpipe.extraction(0, nmrproc_com)

        # Check to see if the baseline correction checkbox is checked
        if (
            dimension_tab.baseline_correction.baseline_correction_checkbox.GetValue()
            == True
        ):
            nmrproc_com = write_nmrpipe.baseline_correction(0, nmrproc_com)

        return nmrproc_com

    def add_indirect_dimension(self, nmrproc_com, dimension, write_nmrpipe):
        """
        Adding the relevent processing lines for the indirect
        dimension to the nmrproc.com file.
        """
        dimension_tab = self.dimension_tabs[dimension]

        # Check to see if the apodization checkbox is checked
        if dimension_tab.apodization.apodization_checkbox.GetValue() == True:
            nmrproc_com = write_nmrpipe.apodization(dimension, nmrproc_com)

        # Check to see if the zero filling checkbox is checked
        if dimension_tab.zero_filling.zero_filling_checkbox.GetValue() == True:
            nmrproc_com = write_nmrpipe.zero_filling(dimension, nmrproc_com)

        # Check to see if the fourier transform checkbox is checked
        if (
            dimension_tab.fourier_transform.fourier_transform_checkbox.GetValue()
            == True
        ):
            nmrproc_com = write_nmrpipe.fourier_transform(dimension, nmrproc_com)

        # Check to see if the phase correction checkbox is checked
        if dimension_tab.phasing.phase_correction_checkbox_indirect.GetValue() == True:
            nmrproc_com = write_nmrpipe.phasing(dimension, nmrproc_com)

        # Check to see if the extraction checkbox is checked
        if dimension_tab.extraction.extraction_checkbox.GetValue() == True:
            nmrproc_com = write_nmrpipe.extraction(dimension, nmrproc_com)

        # Check to see if the baseline correction checkbox is checked
        if (
            dimension_tab.baseline_correction.baseline_correction_checkbox.GetValue()
            == True
        ):
            nmrproc_com = write_nmrpipe.baseline_correction(dimension, nmrproc_com)

        return nmrproc_com

    def check_nus(self, dimension) -> list[int]:
        """
        Checking if NUS reconstruction has been selected.
        Output:
        [0] - no NUS reconstruction for any indirect dimension
        [1] - NUS reconstruction for the first indirect dimension
        [2] - NUS reconstruction for the second indirect dimension
              (if present)
        [1,2] - NUS reconstruction for the first and second indirect
                dimensions.

        The code will also check the nusfile to ensure that the length
        of the rows in nusfile are constistent with the length of the
        output array. This prevents errors where  if there are 2 NUS
        dimensions, but the user has only selected to reconstruct one
        dimension. The user will get a warning to either turn off
        NUS reconstruction or to select both dimensions for NUS
        reconstruction in the graphical interface.
        """
        if (
            self.dimension_tabs[
                dimension
            ].linear_prediction.linear_prediction_radio_box_indirect.GetSelection()
            == 2
        ):
            try:
                if (
                    self.dimension_tabs[
                        dimension + 1
                    ].linear_prediction.linear_prediction_radio_box_indirect.GetSelection()
                    == 2
                ):
                    nus = [1, 2]
                else:
                    nus = [1]
            except:
                nus = [1]
        else:
            nus = [0]

        return nus
    



    def check_lp(self, dimension) -> list[int]:
        """
        Checking if linear prediction has been selected.
        Output:
        0 - linear prediction has not been selected
        1 - linear prediction has been selected
        """
        if (
            self.dimension_tabs[
                dimension
            ].linear_prediction.linear_prediction_radio_box_indirect.GetSelection()
            == 1
        ):
            return 1
        else:
            return 0



    def write_output_line(self, nmrproc_com):
        """
        Writing the output line for nmrPipe
        For 3D data, projections are also created using the
        proj3D.tcl command.
        """

        # Move all existing .dat files to a folder called OldProjections

        current_dir = os.getcwd()
        old_dir = os.path.join(current_dir, "OldProjections")
        
        # Create 'Old' directory if it doesn't exist
        os.makedirs(old_dir, exist_ok=True)
        
        # Loop through files in the current directory
        for filename in os.listdir(current_dir):
            if filename.endswith(".dat") and os.path.isfile(filename):
                source = os.path.join(current_dir, filename)
                destination = os.path.join(old_dir, filename)
                shutil.move(source, destination)


        # Finding the correct name output for the processed data
        nmrfile = self.find_nmrfile_output_name(
            self.nmr_data.dim, self.nmr_data.pseudo_axis
        )
        nmrproc_com.write(" -ov -out {}\n".format(nmrfile))
        if self.nmr_data.dim == 3:
            nmrproc_com.write("proj3D.tcl -in {}".format(nmrfile))





import threading
import sys
import signal
import time
import tempfile
import subprocess
import threading
import tempfile
import os
import signal
import sys

class nmrPipeSubprocess(wx.Frame):
    def __init__(self, parent):
        super().__init__(None, title="nmrPipe Processing", size=(700,500))
        panel = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        self.parent = parent

        # Output display
        self.text_ctrl = wx.TextCtrl(panel, style=wx.TE_MULTILINE | wx.TE_READONLY)
        vbox.Add(self.text_ctrl, 1, wx.EXPAND | wx.ALL, 5)

        # Start/Stop buttons
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.stop_btn = wx.Button(panel, label="Stop")
        self.stop_btn.Disable()
        hbox.Add(self.stop_btn, 0)
        vbox.Add(hbox, 0, wx.ALL, 5)
        panel.SetSizer(vbox)
        self.stop_btn.Bind(wx.EVT_BUTTON, self.on_stop)

        # Timer for polling the file
        self.timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.on_timer, self.timer)

        self.process = None
        self.temp_file = None
        self.last_pos = 0

        self.Show()

        self.on_start()

    def append_text(self, text):
        wx.CallAfter(self.text_ctrl.AppendText, text)

    def on_timer(self, event):
        if not self.temp_file or not self.process:
            return

        # Read new content
        self.temp_file.seek(self.last_pos)
        new_data = self.temp_file.read()
        if new_data:
            try:
                text = new_data.decode(errors="replace")
            except Exception:
                text = str(new_data)
            self.append_text(text)
            self.last_pos = self.temp_file.tell()

            self.last_line_time = time.time()

            lines = text.strip().splitlines()
            for line in lines:
                if hasattr(self, 'last_line'):
                    if line == self.last_line:
                        self.repeat_count += 1
                    else:
                        self.repeat_count = 1
                else:
                    self.repeat_count = 1
                self.last_line = line

                if self.repeat_count >= 100:
                    self.append_text("\nDetected repeated line 100 times, terminating process.\n")
                    self.on_stop(None, user_termination=False)
                    return

        # If process has ended, stop the timer
        if self.process.poll() is not None:
            self.temp_file.seek(self.last_pos)
            remaining = self.temp_file.read()
            if remaining:
                try:
                    self.append_text(remaining.decode(errors="replace"))
                except Exception:
                    self.append_text(str(remaining))
            self.timer.Stop()
            self.append_text(f"\nProcess finished with return code {self.process.returncode}\n")
            wx.CallAfter(self.stop_btn.Disable)
            self.temp_file.close()
            self.temp_file = None
            self.process = None

    def start_subprocess(self, command):
        """Start the subprocess writing stdout/stderr to temp file."""
        # Open temp file in binary mode to avoid text buffering issues
        self.temp_file = tempfile.TemporaryFile(mode="w+b")
        self.last_pos = 0
        

        # Start subprocess in a new process group
        if sys.platform == "win32":
            self.process = subprocess.Popen(
                command,
                shell=True,
                stdout=self.temp_file,
                stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NEW_PROCESS_GROUP
            )
        else:
            self.process = subprocess.Popen(
                command,
                shell=True,
                stdout=self.temp_file,
                stderr=subprocess.STDOUT,
                start_new_session=True
                # preexec_fn=os.setsid
            )

        # Start timer from main thread
        self.timer.Start(100)  # poll every 100ms
        

        # Wait in background thread to stop timer when done
        def wait_process():
            self.process.wait()
            self.timer.Stop()
            self.temp_file.seek(self.last_pos)
            remaining = self.temp_file.read()
            if remaining:
                try:
                    self.append_text(remaining.decode(errors="replace"))
                except Exception:
                    self.append_text(str(remaining))
            if(self.process.returncode==0):
                self.append_text(f"\nProcess finished with no detectable error.\n")
                self.temp_file.flush()
                self.temp_file.close()
                self.temp_file = None
                self.process = None
                wx.CallAfter(self.parent.update_comment)
            else:
                self.temp_file.flush()
                self.temp_file.close()
                self.temp_file = None
                self.process = None
                self.append_text(f"\nProcess finished with errors, check the traceback for sources of error.\n")
            wx.CallAfter(self.stop_btn.Disable)
            

        self.thread = threading.Thread(target=wait_process)
        self.thread.start()

    # def start_subprocess(self, command):
    #     """Start the subprocess writing stdout/stderr to temp file."""
    #     # Open temp file in binary mode to avoid text buffering issues
    #     self.temp_file = tempfile.TemporaryFile(mode="w+b")
    #     self.last_pos = 0
        

    #     # Start subprocess in a new process group
    #     if sys.platform == "win32":
    #         self.process = subprocess.Popen(
    #             command,
    #             shell=True,
    #             stdout=self.temp_file,
    #             stderr=subprocess.STDOUT,
    #             creationflags=subprocess.CREATE_NEW_PROCESS_GROUP
    #         )
    #     else:
    #         self.process = subprocess.Popen(
    #             command,
    #             shell=True,
    #             stdout=self.temp_file,
    #             stderr=subprocess.STDOUT,
    #             start_new_session=True
    #             # preexec_fn=os.setsid
    #         )

    #     # Start timer from main thread
    #     self.timer.Start(100)  # poll every 100ms
        

    #     # Wait in background thread to stop timer when done
    #     def wait_process():
    #         self.process.wait()
    #         self.timer.Stop()
    #         self.temp_file.seek(self.last_pos)
    #         remaining = self.temp_file.read()
    #         if remaining:
    #             try:
    #                 self.append_text(remaining.decode(errors="replace"))
    #             except Exception:
    #                 self.append_text(str(remaining))
    #         if(self.process.returncode==0):
    #             self.append_text(f"\nProcess finished with no detectable error.\n")
    #             self.temp_file.flush()
    #             self.temp_file.close()
    #             self.temp_file = None
    #             self.process = None
    #             wx.CallAfter(lambda: self.parent.update_comment())
    #         else:
    #             self.temp_file.flush()
    #             self.temp_file.close()
    #             self.temp_file = None
    #             self.process = None
    #             self.append_text(f"\nProcess finished with errors, check the traceback for sources of error.\n")
    #         wx.CallAfter(self.stop_btn.Disable)
            

    #     self.thread = threading.Thread(target=wait_process)
    #     self.thread.start()
        # threading.Thread(target=wait_process, daemon=True).start()

    def on_start(self):
        self.stop_btn.Enable()
        command = "./nmrproc.com" 
        self.start_subprocess(command)

    def on_stop(self, event, user_termination=True):
        if self.process and self.process.poll() is None:
            try:
                if sys.platform == "win32":
                    self.process.send_signal(signal.CTRL_BREAK_EVENT)
                else:
                    os.killpg(os.getpgid(self.process.pid), signal.SIGTERM)
                
                if(user_termination):
                    self.append_text("\nProcess terminated by user.\n")
                else:
                    self.append_text("\nProcess terminated by timeout. Scroll through the traceback to find potential sources of error \n")
            except Exception as e:
                self.append_text(f"\nFailed to terminate process: {e}\n")
            self.timer.Stop()
            self.stop_btn.Disable()
            if self.temp_file:
                self.temp_file.close()
                self.temp_file = None
            self.process = None
   
