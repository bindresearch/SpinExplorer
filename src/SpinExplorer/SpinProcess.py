#!/usr/bin/env python3

"""MIT License

Copyright (c) 2025 James Eaton, Andrew Baldwin

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


print('-------------------------------------------------------------')
print('                         SpinProcess                         ')
print('-------------------------------------------------------------')
print('                (version 1.2) 20th June 2025                 ')
print(' (c) 2025 James Eaton, Andrew Baldwin (University of Oxford) ')
print('                        MIT License                          ')
print('-------------------------------------------------------------')
print('                     Processing NMR Data                     ')
print('-------------------------------------------------------------')
print(' Documentation at:')
print(' https://github.com/james-eaton-1/SpinExplorer')
print('-------------------------------------------------------------')
print('')




import sys

import wx
import wx.lib.agw.hyperlink as hl

# Import relevant modules
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import nmrglue as ng
import subprocess
import os
import darkdetect

# Find out the version of operating system being used (Mac, Linux, Windows)
if(sys.platform=='linux'):
    platform = 'linux'
    height = 30
elif(sys.platform=='darwin'):
    platform = 'mac'
    height = 16
else:
    platform = 'windows'
    height = 30


# See if the nmrPipe command works, if not set the platform to windows
if(platform == 'mac' or platform == 'linux'):
    try:
        p = subprocess.Popen('nmrPipe', stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        out, err = p.communicate()
        if('NMRPipe System Version' in str(err)):
            platform=platform
        else:
            platform = 'windows'

    except:
        platform = 'windows'



# James Eaton, 10/06/2025, University of Oxford
# This program is designed to allow the user to process NMR FID data that has been converted to nmrPipe format.


# Suppress complex warning from numpy 
import warnings
# warnings.simplefilter("ignore", np.ComplexWarning)  # For old numpy versions
warnings.simplefilter("ignore", np.exceptions.ComplexWarning)   # For new numpy versions


# Read the FID data from the nmrPipe file 
class ReadFID():
    def __init__(self):

        self.tempframe = wx.Frame(None, -1, 'temp',size=(1,1))
        self.tempframe.Hide()



        self.find_fid()
        self.read_fid()
        self.get_dimensions()

    def find_fid(self):
        # Find the nmrPipe file in the current directory
        fid_files = [file for file in os.listdir() if file.endswith('.fid')]
        if('fids' in os.listdir()):
            fid_files.append('fids')
        if len(fid_files) == 0:
            dlg = wx.MessageDialog(self.tempframe, 'No nmrPipe file found in current directory', 'Error', wx.OK | wx.ICON_ERROR)
            self.tempframe.Raise()
            self.tempframe.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            exit()
        elif len(fid_files) > 1:
            res = ChooseFile(fid_files, self)
            res.ShowModal()
            res.Destroy()
        else:
            self.fid_file = fid_files[0]

    def read_fid(self):
        # Read the nmrPipe file
        self.dic,self.data = ng.pipe.read(self.fid_file)

    def get_dimensions(self):
        try:
            if(platform != 'windows'):
                # Try the command showhdr to get the axis labels
                command = 'showhdr ' + self.fid_file
                # output = subprocess.check_output(command)
                p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
                
                ## Wait for command to terminate. Get return returncode ##
                p_status = p.wait()

                (output, err) = p.communicate()

                
                output = output.decode()
                output = output.split('\n')
                if(output!=['']):
                    self.dim = []
                    for line in output:
                        if('NAME:' in line):
                            line = line.split()[1:]
                            self.dim = len(line)
                            self.axis_labels = line
                            break
                else:
                    fidcom = open('fid.com', 'r')
                    lines = fidcom.readlines()
                    fidcom.close()
                    # Find the line of fid.com with -xLAB in it
                    for line in lines:
                        if('-xLAB' in line):
                            line = line.split('\n')[0].split()
                            if(len(line)==3):
                                self.dim = 1
                                self.axis_labels = [line[1]]
                            elif(len(line)==5):
                                self.dim = 2
                                self.axis_labels = [line[1],line[3]]
                            else:
                                self.dim = 3
                                self.axis_labels = [line[1],line[3],line[5]]
                            break

            else:
                fidcom = open('fid.com', 'r')
                lines = fidcom.readlines()
                fidcom.close()
                # Find the line of fid.com with -xLAB in it
                for line in lines:
                    if('-xLAB' in line):
                        line = line.split('\n')[0].split()
                        if(len(line)==3):
                            self.dim = 1
                            self.axis_labels = [line[1]]
                        elif(len(line)==5):
                            self.dim = 2
                            self.axis_labels = [line[1],line[3]]
                        else:
                            self.dim = 3
                            self.axis_labels = [line[1],line[3],line[5]]
                        break
        
        except:
            # If showhdr fails, use the dimension of the data
            self.dim = len(self.fid.shape)-1
            self.axis_labels = ['X', 'Y', 'Z']

        # Read fid.com file, if any axes are real, flag that as a pseudo axis
        try:
            fidcom = open('fid.com', 'r')
            lines = fidcom.readlines()
            fidcom.close()
            for line in lines:
                if('-xN' in line):
                    line = line.split('\n')[0].split()
                    if(len(line)==3):
                        self.number_of_points = [int(line[1])]
                    elif(len(line)==5):
                        self.number_of_points = [int(line[1]),int(line[3])]
                    else:
                        self.number_of_points = [int(line[1]),int(line[3]),int(line[5])]
                elif('-xSW' in line):
                    line = line.split('\n')[0].split()
                    if(len(line)==3):
                        self.spectral_width = [float(line[1])]
                    elif(len(line)==5):
                        self.spectral_width = [float(line[1]),float(line[3])]
                    else:
                        self.spectral_width = [float(line[1]),float(line[3]),float(line[5])]

                elif('-xMODE' in line):        
                    if('Real' in line):
                        self.pseudo_axis = True
                        if(self.dim == 2):
                            self.index = 0
                        elif(self.dim == 3):
                            line = line.split('\n')[0].split()
                            # Find the index of the real axis in the fid.com file
                            if(line[1]=='Real'):
                                self.index = 0
                            elif(line[3]=='Real'):
                                self.index = 1
                            else:
                                self.index = 2     
                    else:
                        self.pseudo_axis = False
                    continue
                    

        except:
            self.pseudo_axis = False
            self.number_of_points = [1]
            self.spectral_width = [1]





class ChooseFile(wx.Dialog):
    def __init__(self, spectrum_file, parent):
        wx.Dialog.__init__(self, None, wx.ID_ANY, 'Select FID Data', wx.DefaultPosition, size=(300, 200))
        self.spectrum_file = spectrum_file
        self.parent = parent
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.main_sizer.AddSpacer(10)
        self.message = wx.StaticText(self, label='Multiple FID files in current directory. Please select an an option to continue.\n')
        self.main_sizer.Add(self.message, 0, wx.ALL, 5)
        self.file_combobox = wx.ComboBox(self, choices=spectrum_file, style=wx.CB_READONLY)
        self.main_sizer.Add(self.file_combobox, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5)
        self.ok_button = wx.Button(self, label='OK')
        self.ok_button.Bind(wx.EVT_BUTTON, self.OnOK)
        self.main_sizer.Add(self.ok_button, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5)
        self.SetSizer(self.main_sizer)
        self.Centre()

    def OnOK(self, event):
        file_selection = self.file_combobox.GetSelection()
        self.parent.fid_file = self.spectrum_file[file_selection]
        self.Close()
    
        

        


class SpinProcess(wx.Frame):
    def __init__(self, original_frame = None, file_parser = False, path='', cwd='',reprocess=False):
        # Get the monitor size and set the window size to 85% of the monitor size
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 0.7*self.monitorWidth
        self.height = 0.75*self.monitorHeight

        # Read the NMR data in the current directory
        self.nmr_data = ReadFID()

        # Initially set the reprocessing flag to False
        self.reprocess = reprocess
        self.original_frame = original_frame
        self.file_parser = file_parser
        self.path = path
        self.cwd = cwd

        # Create the main window
        self.main_window = wx.Frame.__init__(self, None, title='SpinProcess', size=(self.width, self.height))

        notebook = NotebookProcess(self, self.nmr_data)
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.main_sizer.AddSpacer(10)
        self.main_sizer.Add(notebook, 1, wx.EXPAND)
        notebook.create_buttons(parent=self)

        self.SetSizerAndFit(self.main_sizer)
        # self.SetWindowStyle(wx.STAY_ON_TOP)
        self.Show()
        self.Centre()


        self.Bind(wx.EVT_CLOSE, self.OnClose)

    
    def OnClose(self, event):
        self.Destroy()
        sys.exit()




    def change_frame_size(self, width, height):
        self.SetSize(width, height)

        # Centre the window on the screen
        self.Centre()







    

class NotebookProcess(wx.Notebook):
    def __init__(self, parent, nmr_data):
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 0.7*self.monitorWidth
        self.height = 0.75*self.monitorHeight
        self.parent = parent
        wx.Notebook.__init__(self, parent, id=wx.ID_ANY, style=wx.BK_DEFAULT, size=(self.width, self.height))
        if(darkdetect.isDark() == False):
            self.SetBackgroundColour('#edeeef')

        self.nmr_data = nmr_data
        self.tabDim1 = OneDFrame(self)
        self.AddPage(self.tabDim1, 'Dimension 1 (' + self.nmr_data.axis_labels[0] + ')')
        if(darkdetect.isDark() == True and platform != 'windows'):
            self.tabDim1.SetBackgroundColour((53, 53, 53, 255))
        else:
            self.tabDim1.SetBackgroundColour('#edeeef')
        if(self.nmr_data.dim == 2 and self.nmr_data.pseudo_axis==False):
            self.tabDim2 = TwoDFrame(self, self.tabDim1)
            self.AddPage(self.tabDim2, 'Dimension 2 (' + self.nmr_data.axis_labels[1] + ')')
            if(darkdetect.isDark() == True):
                self.tabDim2.SetBackgroundColour((53, 53, 53, 255))
            else:
                self.tabDim2.SetBackgroundColour('#edeeef')
        if(self.nmr_data.dim == 3 and self.nmr_data.pseudo_axis==True):
            if(self.nmr_data.index == 2):
                self.tabDim2 = TwoDFrame(self, self.tabDim1)
                self.AddPage(self.tabDim2, 'Dimension 2 (' + self.nmr_data.axis_labels[1] + ')')
                if(darkdetect.isDark() == True and platform != 'windows'):
                    self.tabDim2.SetBackgroundColour((53, 53, 53, 255))
                else:
                    self.tabDim2.SetBackgroundColour('#edeeef')
            else:
                self.tabDim2 = TwoDFrame(self, self.tabDim1)
                self.AddPage(self.tabDim2, 'Dimension 2 (' + self.nmr_data.axis_labels[2] + ')')
                if(darkdetect.isDark() == True and platform != 'windows'):
                    self.tabDim2.SetBackgroundColour((53, 53, 53, 255))
                else:
                    self.tabDim2.SetBackgroundColour('#edeeef')
        if(self.nmr_data.dim == 3 and self.nmr_data.pseudo_axis==False):
            self.tabDim2 = TwoDFrame(self, self.tabDim1)
            self.AddPage(self.tabDim2, 'Dimension 2 (' + self.nmr_data.axis_labels[1] + ')')
            if(darkdetect.isDark() == True and platform != 'windows'):
                self.tabDim2.SetBackgroundColour((53, 53, 53, 255))
            else:
                self.tabDim2.SetBackgroundColour('#edeeef')
            self.tabDim3 = ThreeDFrame(self, self.tabDim1)
            self.AddPage(self.tabDim3, 'Dimension 3 (' + self.nmr_data.axis_labels[2] + ')')
            if(darkdetect.isDark() == True and platform != 'windows'):
                self.tabDim3.SetBackgroundColour((53, 53, 53, 255))
            else:
                self.tabDim3.SetBackgroundColour('#edeeef')
            
        
        


    def create_buttons(self, parent):
        # Have a button for make nmrproc.com file, show nmrproc.com file and run processing
        self.button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.make_nmrproc_com_button = wx.Button(parent, -1, 'Make Processing File')
        self.make_nmrproc_com_button.Bind(wx.EVT_BUTTON, self.on_make_nmrproc_com)
        self.button_sizer.Add(self.make_nmrproc_com_button, 0, wx.ALIGN_CENTER_VERTICAL)
        self.button_sizer.AddSpacer(10)
        self.show_nmrproc_com_button = wx.Button(parent, -1, 'Show Processing File')
        self.show_nmrproc_com_button.Bind(wx.EVT_BUTTON, self.on_show_nmrproc_com)
        self.button_sizer.Add(self.show_nmrproc_com_button, 0, wx.ALIGN_CENTER_VERTICAL)
        self.button_sizer.AddSpacer(10)
        self.run_processing_button = wx.Button(parent, -1, 'Run Processing')
        self.run_processing_button.Bind(wx.EVT_BUTTON, self.on_run_processing)
        self.button_sizer.Add(self.run_processing_button, 0, wx.ALIGN_CENTER_VERTICAL)
        self.parent.main_sizer.AddSpacer(20)
        self.parent.main_sizer.Add(self.button_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.parent.main_sizer.AddSpacer(10)



    def on_make_nmrproc_com(self, event):
        # Change path if using unidecFile parser
        self.change_to_path()
        

        # Check to see that the extraction start and end values are valid
        if(self.tabDim1.extraction_checkbox.GetValue() == True):
            # Check that the extraction values are both numbers
            try:
                float(self.tabDim1.extraction_ppm_start_textcontrol.GetValue())
                float(self.tabDim1.extraction_ppm_end_textcontrol.GetValue())
            except:
                dlg = wx.MessageDialog(self, 'Extraction error: The extraction values must be numbers', 'Warning', wx.OK | wx.ICON_WARNING)
                self.Raise()
                self.SetFocus()
                result = dlg.ShowModal()
                self.change_to_cwd()
                return
            if(float(self.tabDim1.extraction_ppm_start_textcontrol.GetValue()) >= float(self.tabDim1.extraction_ppm_end_textcontrol.GetValue())):
                dlg = wx.MessageDialog(self, 'Extraction error: The extraction start value (ppm) must be less than the extraction end value (ppm)', 'Warning', wx.OK | wx.ICON_WARNING)
                self.Raise()
                self.SetFocus()
                result = dlg.ShowModal()
                self.change_to_cwd()
                return
        # Check to see that the phasing values are valid
        if(self.tabDim1.phase_correction_checkbox.GetValue() == True):
            try:
                float(self.tabDim1.phase_correction_p0_textcontrol.GetValue())
                float(self.tabDim1.phase_correction_p1_textcontrol.GetValue())
            except:
                dlg = wx.MessageDialog(self, 'Phasing error: The phase correction values must be numbers', 'Warning', wx.OK | wx.ICON_WARNING)
                self.Raise()
                self.SetFocus()
                result = dlg.ShowModal()
                self.change_to_cwd()
                return
        # Check to see that the node width 
        if(self.tabDim1.baseline_correction_checkbox.GetValue() == True):
            try:
                int(self.tabDim1.baseline_correction_nodes_textcontrol.GetValue())
            except:
                dlg = wx.MessageDialog(self, 'Baseline error: The node width must be an integer number of points', 'Warning', wx.OK | wx.ICON_WARNING)
                self.Raise()
                self.SetFocus()
                result = dlg.ShowModal()
                self.change_to_cwd()
                return
            try:
                node_list = self.tabDim1.baseline_correction_node_list_textcontrol.GetValue().split(',')
                node_list_final = []
                for node in node_list:
                    node_list_final.append(float(node))

                if(len(node_list_final)==0):
                    dlg = wx.MessageDialog(self, 'Baseline error: The node list must contain at least one value', 'Warning', wx.OK | wx.ICON_WARNING)
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    self.change_to_cwd()
                    return
                
            except:
                dlg = wx.MessageDialog(self, 'Baseline error: The node list must be a list of comma separated numbers', 'Warning', wx.OK | wx.ICON_WARNING)
                self.Raise()
                self.SetFocus()
                result = dlg.ShowModal()
                self.change_to_cwd()
                return
            # If polynomial is selected, check to see that the polynomial order is valid
            if(self.tabDim1.baseline_correction_radio_box_selection == 1):
                try:
                    int(self.tabDim1.baseline_correction_polynomial_order_textcontrol.GetValue())
                except:
                    dlg = wx.MessageDialog(self, 'Baseline error: The polynomial order must be an integer for polynomial baselining', 'Warning', wx.OK | wx.ICON_WARNING)
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    self.change_to_cwd()
                    return
                

                
        if(self.nmr_data.dim == 2 and self.nmr_data.pseudo_axis==False):
            include_dim2 = True
        elif(self.nmr_data.dim == 3):
            include_dim2 = True
        else:
            include_dim2 = False
        
        if(include_dim2 == True):
            if(self.tabDim2.extraction_checkbox_dim2.GetValue() == True):
                # Check that the extraction values are both numbers
                try:
                    float(self.tabDim2.extraction_ppm_start_textcontrol_dim2.GetValue())
                    float(self.tabDim2.extraction_ppm_end_textcontrol_dim2.GetValue())
                except:
                    dlg = wx.MessageDialog(self, 'Extraction error (dimension 2): The extraction values must be numbers', 'Warning', wx.OK | wx.ICON_WARNING)
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    self.change_to_cwd()
                    return
                if(float(self.tabDim2.extraction_ppm_start_textcontrol_dim2.GetValue()) >= float(self.tabDim2.extraction_ppm_end_textcontrol_dim2.GetValue())):
                    dlg = wx.MessageDialog(self, 'Extraction error (dimension 2): The extraction start value (ppm) must be less than the extraction end value (ppm)', 'Warning', wx.OK | wx.ICON_WARNING)
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    self.change_to_cwd()
                    return
            if(self.tabDim2.phase_correction_checkbox_dim2.GetValue() == True):
                try:
                    float(self.tabDim2.phase_correction_p0_textcontrol_dim2.GetValue())
                    float(self.tabDim2.phase_correction_p1_textcontrol_dim2.GetValue())
                except:
                    dlg = wx.MessageDialog(self, 'Phasing error (dimension 2): The phase correction values must be numbers', 'Warning', wx.OK | wx.ICON_WARNING)
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    self.change_to_cwd()
                    return
            if(self.tabDim2.baseline_correction_checkbox_dim2.GetValue() == True):
                try:
                    int(self.tabDim2.baseline_correction_nodes_textcontrol_dim2.GetValue())
                except:
                    dlg = wx.MessageDialog(self, 'Baseline error (dimension 2): The node width must be an integer number of points', 'Warning', wx.OK | wx.ICON_WARNING)
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    self.change_to_cwd()
                    return
                try:
                    node_list = self.tabDim2.baseline_correction_node_list_textcontrol_dim2.GetValue().split(',')
                    node_list_final_dim2 = []
                    for node in node_list:
                        node_list_final_dim2.append(float(node))
                    
                    if(len(node_list_final_dim2)==0):
                        dlg = wx.MessageDialog(self, 'Baseline error (dimension 2): The node list must contain at least one value', 'Warning', wx.OK | wx.ICON_WARNING)
                        self.Raise()
                        self.SetFocus()
                        result = dlg.ShowModal()
                        self.change_to_cwd()
                        return
                    
                except:
                    dlg = wx.MessageDialog(self, 'Baseline error (dimension 2): The node list must be a list of comma separated numbers', 'Warning', wx.OK | wx.ICON_WARNING)
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    self.change_to_cwd()
                    return
                # If polynomial is selected, check to see that the polynomial order is valid
                if(self.tabDim2.baseline_correction_radio_box_selection_dim2 == 1):
                    try:
                        int(self.tabDim2.baseline_correction_polynomial_order_textcontrol_dim2.GetValue())
                    except:
                        dlg = wx.MessageDialog(self, 'Baseline error (dimension 2): The polynomial order must be an integer for polynomial baselining', 'Warning', wx.OK | wx.ICON_WARNING)
                        self.Raise()
                        self.SetFocus()
                        result = dlg.ShowModal()
                        self.change_to_cwd()
                        return
            
            # If SMILE processing is selected, check to see that the SMILE file exists
            if(self.tabDim2.linear_prediction_radio_box_dim2.GetSelection() == 2):
                # SMILE processing is selected, ask the user to confirm SMILE is installed as part of nmrPipe
                # Have 2 options for the user to select, either to continue or to cancel
                dlg = wx.MessageDialog(self, 'SMILE processing is selected. Ensure that SMILE is installed as part of nmrPipe', 'Warning', wx.OK | wx.CANCEL | wx.ICON_WARNING)
                self.Raise()
                self.SetFocus()
                result = dlg.ShowModal()
                if(result == wx.ID_CANCEL):
                    self.change_to_cwd()
                    return
                
                if(self.tabDim1.extraction_checkbox.GetValue()==False):
                    dlg = wx.MessageDialog(self, 'No direct dimension data extraction is selected, SMILE reconstruction may take a while. Consider extracting a region of the direct dimension before reconstruction. Do you want to continue or cancel?', 'Warning', wx.OK | wx.CANCEL | wx.ICON_WARNING)
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    if(result == wx.ID_CANCEL):
                        return

                

                # List the files in the current directory
                files = os.listdir()
                # Check to see if the SMILE file exists
                if(self.tabDim2.smile_nus_file_textcontrol_dim2.GetValue() not in files):
                    message = 'SMILE NUS reconstruction error (dimension 2): The NUS file ' + self.tabDim2.smile_nus_file_textcontrol_dim2.GetValue() + ' cannot be found in the current directory'
                    dlg = wx.MessageDialog(self, message, 'Warning', wx.OK | wx.ICON_WARNING)
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    self.change_to_cwd()
                    return

                # Check that the number of CPUs is an integer
                try:
                    int(self.tabDim2.smile_nus_cpu_textcontrol_dim2.GetValue())
                except:
                    dlg = wx.MessageDialog(self, 'SMILE NUS reconstruction error (dimension 2): The number of CPUs must be an integer', 'Warning', wx.OK | wx.ICON_WARNING)
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    self.change_to_cwd()
                    return
                
                # Check that the number of iterations is an integer
                try:
                    int(self.tabDim2.smile_nus_iterations_textcontrol_dim2.GetValue())
                except:
                    dlg = wx.MessageDialog(self, 'SMILE NUS reconstruction error (dimension 2): The number of iterations must be an integer', 'Warning', wx.OK | wx.ICON_WARNING)
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    self.change_to_cwd()
                    return

        if(self.nmr_data.dim == 3 and self.nmr_data.pseudo_axis==False):
            include_dim3 = True
            if(self.tabDim3.extraction_checkbox_dim3.GetValue() == True):
                # Check that the extraction values are both numbers
                try:
                    float(self.tabDim3.extraction_ppm_start_textcontrol_dim3.GetValue())
                    float(self.tabDim3.extraction_ppm_end_textcontrol_dim3.GetValue())
                except:
                    dlg = wx.MessageDialog(self, 'Extraction error (dimension 3): The extraction values must be numbers', 'Warning', wx.OK | wx.ICON_WARNING)
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    self.change_to_cwd()
                    return
                if(float(self.tabDim3.extraction_ppm_start_textcontrol_dim3.GetValue()) >= float(self.tabDim3.extraction_ppm_end_textcontrol_dim3.GetValue())):
                    dlg = wx.MessageDialog(self, 'Extraction error (dimension 3): The extraction start value (ppm) must be less than the extraction end value (ppm)', 'Warning', wx.OK | wx.ICON_WARNING)
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    self.change_to_cwd()
                    return
            if(self.tabDim3.phase_correction_checkbox_dim3.GetValue() == True):
                try:
                    float(self.tabDim3.phase_correction_p0_textcontrol_dim3.GetValue())
                    float(self.tabDim3.phase_correction_p1_textcontrol_dim3.GetValue())
                except:
                    dlg = wx.MessageDialog(self, 'Phasing error (dimension 3): The phase correction values must be numbers', 'Warning', wx.OK | wx.ICON_WARNING)
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    self.change_to_cwd()
                    return
            if(self.tabDim3.baseline_correction_checkbox_dim3.GetValue() == True):
                try:
                    int(self.tabDim3.baseline_correction_nodes_textcontrol_dim3.GetValue())
                except:
                    dlg = wx.MessageDialog(self, 'Baseline error (dimension 3): The node width must be an integer number of points', 'Warning', wx.OK | wx.ICON_WARNING)
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    self.change_to_cwd()
                    return
                try:
                    node_list = self.tabDim3.baseline_correction_node_list_textcontrol_dim3.GetValue().split(',')
                    node_list_final_dim3 = []
                    for node in node_list:
                        node_list_final_dim3.append(float(node))
                    
                    if(len(node_list_final_dim3)==0):
                        dlg = wx.MessageDialog(self, 'Baseline error (dimension 3): The node list must contain at least one value', 'Warning', wx.OK | wx.ICON_WARNING)
                        self.Raise()
                        self.SetFocus()
                        result = dlg.ShowModal()
                        self.change_to_cwd()

                except:
                    dlg = wx.MessageDialog(self, 'Baseline error (dimension 3): The node list must be a list of comma separated numbers', 'Warning', wx.OK | wx.ICON_WARNING)
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    self.change_to_cwd()
                    return
        else:
            include_dim3 = False
            

        # Check to see if the nmrproc.com file already exists, if it does ask the user if they want to overwrite it
        if(os.path.exists('./nmrproc.com')):
            dlg = wx.MessageDialog(self, 'The nmrproc.com file already exists. Do you want to overwrite it?', 'Warning', wx.YES_NO | wx.ICON_WARNING)
            self.Raise()
            self.SetFocus()
            result = dlg.ShowModal()
            if(result == wx.ID_NO):
                self.change_to_cwd()
                return

        # Create the nmrproc.com file
        nmrproc_com = open('nmrproc.com', 'w')
        nmrproc_com.write('#!/bin/csh\n\n')
        if('/' not in self.nmr_data.fid_file):
            nmrproc_com.write('nmrPipe -in ' + self.nmr_data.fid_file + '\\\n')
        else:
            nmrproc_com.write('xyz2pipe -in ./fids/test%03d.fid \\\n')

        # Check to see if the solvent suppression checkbox is checked
        if(self.tabDim1.solvent_suppression_checkbox.GetValue() == True):
            solvent_suppression_line = '| nmrPipe -fn SOL'
            if(self.tabDim1.solvent_suppression_filter_selection == 0):
                solvent_suppression_line += ' -mode 1'
                if(self.tabDim1.solvent_suppression_lowpass_shape_selection == 0):
                    solvent_suppression_line += ' -fs 1'
                elif(self.tabDim1.solvent_suppression_lowpass_shape_selection == 1):
                    solvent_suppression_line += ' -fs 2'
                elif(self.tabDim1.solvent_suppression_lowpass_shape_selection == 2):
                    solvent_suppression_line += ' -fs 3'
            elif(self.tabDim1.solvent_suppression_filter_selection == 1):
                solvent_suppression_line += ' -mode 2'
            elif(self.tabDim1.solvent_suppression_filter_selection == 2):
                solvent_suppression_line += ' -mode 3'

            nmrproc_com.write(solvent_suppression_line + ' \\\n')

        # Check to see if the linear prediction checkbox is checked
        if(self.tabDim1.linear_prediction_checkbox.GetValue() == True):
            linear_prediction_line = '| nmrPipe -fn LP'
            if(self.tabDim1.linear_prediction_options_selection == 0):
                linear_prediction_line += ' -after'
            elif(self.tabDim1.linear_prediction_options_selection == 1):
                linear_prediction_line += ' -before'
            if(self.tabDim1.linear_prediction_coefficients_selection == 0):
                linear_prediction_line += ' -f'
            elif(self.tabDim1.linear_prediction_coefficients_selection == 1):
                linear_prediction_line += ' -b'
            elif(self.tabDim1.linear_prediction_coefficients_selection == 2):
                linear_prediction_line += ' -fb'
            
            nmrproc_com.write(linear_prediction_line + ' \\\n')
        
        # Check to see if the apodization checkbox is checked
        if(self.tabDim1.apodization_checkbox.GetValue() == True):
            apodization_line = '| nmrPipe -fn'
            if(self.tabDim1.apodization_combobox_selection == 0):
                apodization_line += ' EM'
                # Input a line broadening of 0 and first point scaling
                apodization_line += ' -lb 0 -c 0.5'
            elif(self.tabDim1.apodization_combobox_selection == 1):
                apodization_line += ' EM'
                apodization_line += ' -lb ' + str(self.tabDim1.exponential_line_broadening) + ' -c ' + str(self.tabDim1.apodization_first_point_scaling)
            elif(self.tabDim1.apodization_combobox_selection == 2):
                apodization_line += ' GM'
                apodization_line += ' -g1 ' + str(self.tabDim1.g1) + ' -g2 ' + str(self.tabDim1.g2) + ' -g3 ' + str(self.tabDim1.g3) + ' -c ' + str(self.tabDim1.apodization_first_point_scaling)
            elif(self.tabDim1.apodization_combobox_selection == 3):
                apodization_line += ' SP'
                apodization_line += ' -off ' + str(self.tabDim1.offset) + ' -end ' + str(self.tabDim1.end) + ' -pow ' + str(self.tabDim1.power) + ' -c ' + str(self.tabDim1.apodization_first_point_scaling)
            elif(self.tabDim1.apodization_combobox_selection == 4):
                apodization_line += ' GMB'
                apodization_line += ' -lb ' + str(self.tabDim1.a) + ' -gb ' + str(self.tabDim1.b) + ' -c ' + str(self.tabDim1.apodization_first_point_scaling)
            elif(self.tabDim1.apodization_combobox_selection == 5):
                apodization_line += ' TM'
                apodization_line += ' -t1 ' + str(self.tabDim1.t1) + ' -t2 ' + str(self.tabDim1.t2) + ' -c ' + str(self.tabDim1.apodization_first_point_scaling)
            elif(self.tabDim1.apodization_combobox_selection == 6):
                apodization_line += ' TRI'
                apodization_line += ' -loc ' + str(self.tabDim1.loc) + ' -c ' + str(self.tabDim1.apodization_first_point_scaling)

            nmrproc_com.write(apodization_line + ' \\\n')

        # Check to see if the zero filling checkbox is checked
        if(self.tabDim1.zero_filling_checkbox.GetValue() == True):
            zero_filling_line = '| nmrPipe -fn ZF'
            if(self.tabDim1.zero_filling_combobox_selection == 0):
                zero_filling_line += ' -zf ' + str(self.tabDim1.zero_filling_value_doubling_times)
            if(self.tabDim1.zero_filling_combobox_selection == 1):
                zero_filling_line += ' -pad ' + str(self.tabDim1.zero_filling_value_zeros_to_add)
            elif(self.tabDim1.zero_filling_combobox_selection == 2):
                zero_filling_line += ' -size ' + str(self.tabDim1.zero_filling_value_final_data_size)
            if(self.tabDim1.zero_filling_round_checkbox.GetValue() == True):
                zero_filling_line += ' -auto'

            nmrproc_com.write(zero_filling_line + ' \\\n')

        # Check to see if the fourier transform checkbox is checked
        if(self.tabDim1.fourier_transform_checkbox.GetValue() == True):
            fourier_transform_line = '| nmrPipe -fn FT'
            if(self.tabDim1.ft_method_selection == 0):
                fourier_transform_line += ' -auto'
            elif(self.tabDim1.ft_method_selection == 1):
                fourier_transform_line += ' -real'
            elif(self.tabDim1.ft_method_selection == 2):
                fourier_transform_line += ' -inv'
            elif(self.tabDim1.ft_method_selection == 3):
                fourier_transform_line += ' -alt'
            elif(self.tabDim1.ft_method_selection == 4):
                fourier_transform_line += ' -neg'

            nmrproc_com.write(fourier_transform_line + ' \\\n')

        # Check to see if the phase correction checkbox is checked
        if(self.tabDim1.phase_correction_checkbox.GetValue() == True):
            phase_correction_line = '| nmrPipe -fn PS'
            phase_correction_line += ' -p0 ' + str(self.tabDim1.phase_correction_p0_textcontrol.GetValue()) + ' -p1 ' + str(self.tabDim1.phase_correction_p1_textcontrol.GetValue()) + ' -di '

            nmrproc_com.write(phase_correction_line + '\\\n')
        else:
            # Check to see if magnitude mode is toggles
            if(self.tabDim1.magnitude_mode_checkbox.GetValue()==True):
                mc_line = '| nmrPipe -fn MC'
                nmrproc_com.write(mc_line + '\\\n')

        # Check to see if the extraction checkbox is checked
        if(self.tabDim1.extraction_checkbox.GetValue() == True):
            extraction_line = '| nmrPipe -fn EXT'
            extraction_line += ' -x1 ' + str(self.tabDim1.extraction_ppm_start_textcontrol.GetValue()) + 'ppm -xn ' + str(self.tabDim1.extraction_ppm_end_textcontrol.GetValue()) + 'ppm -sw '

            nmrproc_com.write(extraction_line + ' \\\n')

        # Check to see if the baseline correction checkbox is checked
        if(self.tabDim1.baseline_correction_checkbox.GetValue() == True):
            if(self.tabDim1.baseline_correction_radio_box_selection == 0):
                baseline_correction_line = '| nmrPipe -fn BASE'
                # add the node width and node list to the baseline correction line
                baseline_correction_line += ' -nw ' + str(self.tabDim1.baseline_correction_nodes_textcontrol.GetValue()) + ' -nl '
                for node in node_list_final:
                    baseline_correction_line += str(node) + '% '
                nmrproc_com.write(baseline_correction_line + '\\\n')
            elif(self.tabDim1.baseline_correction_radio_box_selection == 1):
                baseline_correction_line = '| nmrPipe -fn POLY'
                # add the node width and node list to the baseline correction line
                baseline_correction_line += ' -nw ' + str(self.tabDim1.baseline_correction_node_list_textcontrol.GetValue()) + ' -nl '
                for node in node_list_final:
                    baseline_correction_line += str(node) + '% '
                # add the polynomial order to the baseline correction line
                baseline_correction_line += ' -ord ' + str(self.tabDim1.baseline_correction_polynomial_order_textcontrol.GetValue())
                nmrproc_com.write(baseline_correction_line + ' \\\n')


        if(include_dim2 == True):
            if(self.tabDim2.linear_prediction_radio_box_dim2.GetSelection() == 2):
                if(include_dim3==True):
                    if(self.tabDim3.linear_prediction_radio_box_dim3.GetSelection()==2):
                        nmrproc_com.write('| nmrPipe -fn ZTP \\\n')
                        nmrproc_com.write('| nmrPipe -fn TP \\\n')
                        smile_line = '| nmrPipe -fn SMILE -nDim ' + str(int(self.nmr_data.dim)) + ' -nThread ' + str(self.tabDim2.smile_nus_cpu_textcontrol_dim2.GetValue()) + ' -report 1 -sample ' + str(self.tabDim2.nuslist_name_dim2) + '\\\n' + '              -maxIter ' + str(int(self.tabDim2.smile_nus_iterations_textcontrol_dim2.GetValue())) + '        \\\n' + '       -xP0 {} -xP1 {} -xT {} '.format(self.tabDim2.phase_correction_p0_textcontrol_dim2.GetValue(), self.tabDim2.phase_correction_p1_textcontrol_dim2.GetValue(), self.tabDim2.smile_nus_extension_textcontrol_dim2.GetValue())
                        nmrproc_com.write(smile_line + ' \\\n')
                        if(self.tabDim2.apodization_dim2_checkbox_value==True):
                            if(self.tabDim2.apodization_dim2_combobox_selection==0):
                                smile_line_apod = ' -xApod EM -xQ1 0.0 -xQ2 0.0 -xQ3 0.0 \\\n'
                            elif(self.tabDim2.apodization_dim2_combobox_selection==1):
                                smile_line_apod = ' -xApod EM -xQ1 {} -xQ2 0.0 -xQ3 0.0 \\\n'.format(str(self.tabDim2.exponential_line_broadening_dim2))
                            elif(self.tabDim2.apodization_dim2_combobox_selection==2):
                                smile_line_apod = ' -xApod GM -xQ1 {} -xQ2 {} -xQ3 {} \\\n'.format(str(self.tabDim2.g1_dim2), str(self.tabDim2.g2_dim2),str(self.tabDim2.g3_dim2))
                            elif(self.tabDim2.apodization_dim2_combobox_selection==3):
                                smile_line_apod = ' -xApod SP -xQ1 {} -xQ2 {} -xQ3 {} \\\n'.format(str(self.tabDim2.offset_dim2), str(self.tabDim2.end_dim2),str(self.tabDim2.power_dim2))
                            elif(self.tabDim2.apodization_dim2_combobox_selection==4):
                                smile_line_apod = ' -xApod GMB -xQ1 {} -xQ2 {} -xQ3 0.0 \\\n'.format(str(self.tabDim2.a_dim2), str(self.tabDim2.b_dim2))
                            elif(self.tabDim2.apodization_dim2_combobox_selection==5):
                                smile_line_apod = ' -xApod TM -xQ1 {} -xQ2 {} -xQ3 0.0 \\\n'.format(str(self.tabDim2.t1_dim2), str(self.tabDim2.t2_dim2))
                            elif(self.tabDim2.apodization_dim2_combobox_selection==5):
                                smile_line_apod = ' -xApod TRI -xQ1 {} -xQ2 0.0 -xQ3 0.0 \\\n'.format(str(self.tabDim2.loc_dim2))
                        nmrproc_com.write(smile_line_apod)

                        smile_line2 = '       -yP0 {} -yP1 {} -yT {} '.format(self.tabDim3.phase_correction_p0_textcontrol_dim3.GetValue(), self.tabDim3.phase_correction_p1_textcontrol_dim3.GetValue(), self.tabDim3.smile_nus_data_extension_textcontrol_dim3.GetValue())
                        nmrproc_com.write(smile_line2 + ' \\\n')
                        if(self.tabDim3.apodization_dim3_checkbox_value==True):
                            if(self.tabDim3.apodization_dim3_combobox_selection==0):
                                smile_line_apod = ' -yApod EM -yQ1 0.0 -yQ2 0.0 -yQ3 0.0 \\\n'
                            elif(self.tabDim3.apodization_dim3_combobox_selection==1):
                                smile_line_apod = ' -yApod EM -yQ1 {} -yQ2 0.0 -yQ3 0.0 \\\n'.format(str(self.tabDim3.exponential_line_broadening_dim3))
                            elif(self.tabDim3.apodization_dim3_combobox_selection==2):
                                smile_line_apod = ' -yApod GM -yQ1 {} -yQ2 {} -yQ3 {} \\\n'.format(str(self.tabDim3.g1_dim3), str(self.tabDim3.g2_dim3),str(self.tabDim3.g3_dim3))
                            elif(self.tabDim3.apodization_dim3_combobox_selection==3):
                                smile_line_apod = ' -yApod SP -yQ1 {} -yQ2 {} -yQ3 {} \\\n'.format(str(self.tabDim3.offset_dim3), str(self.tabDim3.end_dim3),str(self.tabDim3.power_dim3))
                            elif(self.tabDim3.apodization_dim3_combobox_selection==4):
                                smile_line_apod = ' -yApod GMB -yQ1 {} -yQ2 {} -yQ3 0.0 \\\n'.format(str(self.tabDim3.a_dim3), str(self.tabDim3.b_dim3))
                            elif(self.tabDim3.apodization_dim3_combobox_selection==5):
                                smile_line_apod = ' -yApod TM -yQ1 {} -yQ2 {} -yQ3 0.0 \\\n'.format(str(self.tabDim3.t1_dim3), str(self.tabDim3.t2_dim3))
                            elif(self.tabDim3.apodization_dim3_combobox_selection==5):
                                smile_line_apod = ' -yApod TRI -yQ1 {} -yQ2 0.0 -yQ3 0.0 \\\n'.format(str(self.tabDim3.loc_dim3))
                        nmrproc_com.write(smile_line_apod)
                        nmrproc_com.write('| nmrPipe -fn TP \\\n')
                        nmrproc_com.write('| nmrPipe -fn ZTP \\\n')
                else:
                    nmrproc_com.write('| nmrPipe -fn TP \\\n')
                    smile_line = '| nmrPipe -fn SMILE -nDim ' + str(int(self.nmr_data.dim)) + ' -nThread ' + str(self.tabDim2.smile_nus_cpu_textcontrol_dim2.GetValue()) + ' -report 1 -sample ' + str(self.tabDim2.nuslist_name_dim2) + '\\\n' + '              -maxIter ' + str(int(self.tabDim2.smile_nus_iterations_textcontrol_dim2.GetValue())) + '        \\\n' + '       -xP0 {} -xP1 {} -xT {} '.format(self.tabDim2.phase_correction_p0_textcontrol_dim2.GetValue(), self.tabDim2.phase_correction_p1_textcontrol_dim2.GetValue(), self.tabDim2.smile_nus_extension_textcontrol_dim2.GetValue())
                    nmrproc_com.write(smile_line + ' \\\n')
                    if(self.tabDim2.apodization_dim2_checkbox_value==True):
                        if(self.tabDim2.apodization_dim2_combobox_selection==0):
                            smile_line_apod = ' -xApod EM -xQ1 0.0 -xQ2 0.0 -xQ3 0.0 \\\n'
                        elif(self.tabDim2.apodization_dim2_combobox_selection==1):
                            smile_line_apod = ' -xApod EM -xQ1 {} -xQ2 0.0 -xQ3 0.0 \\\n'.format(str(self.tabDim2.exponential_line_broadening_dim2))
                        elif(self.tabDim2.apodization_dim2_combobox_selection==2):
                            smile_line_apod = ' -xApod GM -xQ1 {} -xQ2 {} -xQ3 {} \\\n'.format(str(self.tabDim2.g1_dim2), str(self.tabDim2.g2_dim2),str(self.tabDim2.g3_dim2))
                        elif(self.tabDim2.apodization_dim2_combobox_selection==3):
                            smile_line_apod = ' -xApod SP -xQ1 {} -xQ2 {} -xQ3 {} \\\n'.format(str(self.tabDim2.offset_dim2), str(self.tabDim2.end_dim2),str(self.tabDim2.power_dim2))
                        elif(self.tabDim2.apodization_dim2_combobox_selection==4):
                            smile_line_apod = ' -xApod GMB -xQ1 {} -xQ2 {} -xQ3 0.0 \\\n'.format(str(self.tabDim2.a_dim2), str(self.tabDim2.b_dim2))
                        elif(self.tabDim2.apodization_dim2_combobox_selection==5):
                            smile_line_apod = ' -xApod TM -xQ1 {} -xQ2 {} -xQ3 0.0 \\\n'.format(str(self.tabDim2.t1_dim2), str(self.tabDim2.t2_dim2))
                        elif(self.tabDim2.apodization_dim2_combobox_selection==5):
                            smile_line_apod = ' -xApod TRI -xQ1 {} -xQ2 0.0 -xQ3 0.0 \\\n'.format(str(self.tabDim2.loc_dim2))
                        nmrproc_com.write(smile_line_apod)

            if(self.nmr_data.pseudo_axis == False):
                # Transpose the data
                nmrproc_com.write('| nmrPipe -fn TP \\\n')
            elif(self.nmr_data.pseudo_axis == True):
                if(self.nmr_data.index == 2):
                    nmrproc_com.write('| nmrPipe -fn TP \\\n')
                elif(self.nmr_data.index == 1):
                    nmrproc_com.write('| nmrPipe -fn TP \\\n')
                    nmrproc_com.write('| nmrPipe -fn ZTP \\\n')

            # Check to see if LP or SMILE processing is selected
            if(self.tabDim2.linear_prediction_radio_box_dim2.GetSelection() == 0):
                # Linear prediction is not selected
                pass
            elif(self.tabDim2.linear_prediction_radio_box_dim2.GetSelection() == 1):
                # Linear prediction is selected
                # Get the linear prediction options and coefficients
                linear_prediction_line = '| nmrPipe -fn LP'
                if(self.tabDim2.linear_prediction_dim2_options_selection == 0):
                    linear_prediction_line += ' -after'
                elif(self.tabDim2.linear_prediction_dim2_options_selection == 1):
                    linear_prediction_line += ' -before'
                if(self.tabDim2.linear_prediction_dim2_coefficients_selection == 0):
                    linear_prediction_line += ' -f'
                elif(self.tabDim2.linear_prediction_dim2_coefficients_selection == 1):
                    linear_prediction_line += ' -b'
                elif(self.tabDim2.linear_prediction_dim2_coefficients_selection == 2):
                    linear_prediction_line += ' -fb'
                # Add the linear prediction line to the nmrproc.com file
                nmrproc_com.write(linear_prediction_line + ' \\\n')

            # Check to see if apodization is checked 
            if(self.tabDim2.linear_prediction_radio_box_dim2.GetSelection() != 2):
                if(self.tabDim2.apodization_dim2_checkbox_value==True):
                    apodization_line = '| nmrPipe -fn '
                    if(self.tabDim2.apodization_dim2_combobox_selection==0):
                        apodization_line += 'EM -lb 0.0 -c ' + str(self.tabDim2.apodization_first_point_scaling_dim2)
                    elif(self.tabDim2.apodization_dim2_combobox_selection==1):
                        apodization_line += 'EM -lb ' + str(self.tabDim2.exponential_line_broadening_dim2) + ' -c ' + str(self.tabDim2.apodization_first_point_scaling_dim2)
                    elif(self.tabDim2.apodization_dim2_combobox_selection==2):
                        apodization_line += 'GM -g1 ' + str(self.tabDim2.g1_dim2) + ' -g2 ' + str(self.tabDim2.g2_dim2) + ' -g3 ' + str(self.tabDim2.g3_dim2) + ' -c ' + str(self.tabDim2.apodization_first_point_scaling_dim2)
                    elif(self.tabDim2.apodization_dim2_combobox_selection==3):
                        apodization_line += ' SP'
                        apodization_line += ' -off ' + str(self.tabDim2.offset_dim2) + ' -end ' + str(self.tabDim2.end_dim2) + ' -pow ' + str(self.tabDim2.power_dim2) + ' -c ' + str(self.tabDim2.apodization_first_point_scaling_dim2)
                    elif(self.tabDim2.apodization_dim2_combobox_selection == 4):
                        apodization_line += ' GMB'
                        apodization_line += ' -a ' + str(self.tabDim2.a_dim2) + ' -b ' + str(self.tabDim2.b_dim2) + ' -c ' + str(self.tabDim2.apodization_first_point_scaling_dim2)
                    elif(self.tabDim1.apodization_dim2_combobox_selection == 5):
                        apodization_line += ' TP'
                        apodization_line += ' -t1 ' + str(self.tabDim2.t1_dim2) + ' -t2 ' + str(self.tabDim2.t2_dim2) + ' -c ' + str(self.tabDim2.apodization_first_point_scaling_dim2)
                    elif(self.tabDim1.apodization_dim2_combobox_selection == 6):
                        apodization_line += ' TRI'
                        apodization_line += ' -loc ' + str(self.tabDim2.loc_dim2) + ' -c ' + str(self.tabDim2.apodization_first_point_scaling_dim2)
                    nmrproc_com.write(apodization_line + ' \\\n')


            # Check to see if zero filling is selected
            if(self.tabDim2.zero_filling_checkbox_dim2.GetValue() == True):
                zero_filling_line = '| nmrPipe -fn ZF'
                if(self.tabDim2.zero_filling_dim2_combobox_selection == 0):
                    zero_filling_line += ' -zf ' + str(self.tabDim2.zero_filling_dim2_value_doubling_times)
                if(self.tabDim2.zero_filling_dim2_combobox_selection == 1):
                    zero_filling_line += ' -pad ' + str(self.tabDim2.zero_filling_dim2_value_zeros_to_add)
                elif(self.tabDim2.zero_filling_dim2_combobox_selection == 2):
                    zero_filling_line += ' -size ' + str(self.tabDim2.zero_filling_dim2_value_final_data_size)
                if(self.tabDim2.zero_filling_round_checkbox_dim2.GetValue() == True):
                    zero_filling_line += ' -auto'

                nmrproc_com.write(zero_filling_line + ' \\\n')
            
            # Check to see if the fourier transform checkbox is checked
            if(self.tabDim2.fourier_transform_checkbox_dim2.GetValue() == True):
                fourier_transform_line = '| nmrPipe -fn FT'
                if(self.tabDim2.ft_method_selection_dim2 == 0):
                    fourier_transform_line += ' -auto'
                elif(self.tabDim2.ft_method_selection_dim2 == 1):
                    fourier_transform_line += ' -real'
                elif(self.tabDim2.ft_method_selection_dim2 == 2):
                    fourier_transform_line += ' -inv'
                elif(self.tabDim2.ft_method_selection_dim2 == 3):
                    fourier_transform_line += ' -alt'
                elif(self.tabDim2.ft_method_selection_dim2 == 4):
                    fourier_transform_line += ' -neg'

                nmrproc_com.write(fourier_transform_line + ' \\\n')
            
            # Check to see if the phase correction checkbox is checked
            if(self.tabDim2.phase_correction_checkbox_dim2.GetValue() == True):
                phase_correction_line = '| nmrPipe -fn PS'
                phase_correction_line += ' -p0 ' + str(self.tabDim2.phase_correction_p0_textcontrol_dim2.GetValue()) + ' -p1 ' + str(self.tabDim2.phase_correction_p1_textcontrol_dim2.GetValue()) + ' -di '

                nmrproc_com.write(phase_correction_line + '\\\n')
            
            # Check to see if the extraction checkbox is checked
            if(self.tabDim2.extraction_checkbox_dim2.GetValue() == True):
                extraction_line = '| nmrPipe -fn EXT'
                extraction_line += ' -x1 ' + str(self.tabDim2.extraction_ppm_start_textcontrol_dim2.GetValue()) + 'ppm'+ ' -xn ' + str(self.tabDim2.extraction_ppm_end_textcontrol_dim2.GetValue()) + 'ppm -sw '

                nmrproc_com.write(extraction_line + ' \\\n')
            
            # Check to see if the baseline correction checkbox is checked
            if(self.tabDim2.baseline_correction_checkbox_dim2.GetValue() == True):
                if(self.tabDim2.baseline_correction_radio_box_selection_dim2 == 0):
                    baseline_correction_line = '| nmrPipe -fn BASE'
                    # add the node width and node list to the baseline correction line
                    baseline_correction_line += ' -nw ' + str(self.tabDim2.baseline_correction_nodes_textcontrol_dim2.GetValue()) + ' -nl '
                    for node in node_list_final_dim2:
                        baseline_correction_line += str(node) + '% '
                    nmrproc_com.write(baseline_correction_line + '\\\n')
                elif(self.tabDim2.baseline_correction_radio_box_selection_dim2 == 1):
                    baseline_correction_line = '| nmrPipe -fn POLY'
                    # add the node width and node list to the baseline correction line
                    baseline_correction_line += ' -nw ' + str(self.tabDim2.baseline_correction_nodes_textcontrol_dim2.GetValue()) + ' -nl '
                    for node in node_list_final_dim2:
                        baseline_correction_line += str(node) + '% '
                    # add the polynomial order to the baseline correction line
                    baseline_correction_line += ' -ord ' + str(self.tabDim2.baseline_correction_polynomial_order_textcontrol_dim2.GetValue())
                    nmrproc_com.write(baseline_correction_line + ' \\\n')
            
            # if('/' in self.nmr_data.fid_file):
            #     # Then processing 3D with fids
            #     nmrproc_com.write('| pipe2xyz -out ./ft/t1%03d.ft3 -x -verb -ov \n')
            #     nmrproc_com.write('xyz2pipe -in ./ft/t1%03d.ft3 -x -verb \\')



        if(include_dim3 == True):
            # Zero transpose the data
            nmrproc_com.write('| nmrPipe -fn ZTP \\\n')
            # Check to see if LP or SMILE processing is selected
            if(self.tabDim3.linear_prediction_radio_box_dim3.GetSelection() == 0):
                # Linear prediction is not selected
                pass
            elif(self.tabDim3.linear_prediction_radio_box_dim3.GetSelection() == 1):
                # Linear prediction is selected
                # Get the linear prediction options and coefficients
                linear_prediction_line = '| nmrPipe -fn LP'
                if(self.tabDim3.linear_prediction_dim3_options_selection == 0):
                    linear_prediction_line += ' -after'
                elif(self.tabDim3.linear_prediction_dim3_options_selection == 1):
                    linear_prediction_line += ' -before'
                if(self.tabDim3.linear_prediction_dim3_coefficients_selection == 0):
                    linear_prediction_line += ' -f'
                elif(self.tabDim3.linear_prediction_dim3_coefficients_selection == 1):
                    linear_prediction_line += ' -b'
                elif(self.tabDim3.linear_prediction_dim3_coefficients_selection == 2):
                    linear_prediction_line += ' -fb'
                # Add the linear prediction line to the nmrproc.com file
                nmrproc_com.write(linear_prediction_line + ' \\\n')


            # Check to see if apodization is checked 
            if(self.tabDim3.linear_prediction_radio_box_dim3.GetSelection() != 2):
                if(self.tabDim3.apodization_dim3_checkbox_value==True):
                    apodization_line = '| nmrPipe -fn '
                    if(self.tabDim3.apodization_dim3_combobox_selection==0):
                        apodization_line += 'EM -lb 0.0 -c ' + str(self.tabDim3.apodization_first_point_scaling_dim3)
                    elif(self.tabDim3.apodization_dim3_combobox_selection==1):
                        apodization_line += 'EM -lb ' + str(self.tabDim3.exponential_line_broadening_dim3) + ' -c ' + str(self.tabDim3.apodization_first_point_scaling_dim3)
                    elif(self.tabDim3.apodization_dim3_combobox_selection==2):
                        apodization_line += 'GM -g1 ' + str(self.tabDim3.g1_dim3) + ' -g2 ' + str(self.tabDim3.g2_dim3) + ' -g3 ' + str(self.tabDim3.g3_dim3) + ' -c ' + str(self.tabDim3.apodization_first_point_scaling_dim3)
                    elif(self.tabDim3.apodization_dim3_combobox_selection==3):
                        apodization_line += ' SP'
                        apodization_line += ' -off ' + str(self.tabDim3.offset_dim3) + ' -end ' + str(self.tabDim3.end_dim3) + ' -pow ' + str(self.tabDim3.power_dim3) + ' -c ' + str(self.tabDim3.apodization_first_point_scaling_dim3)
                    elif(self.tabDim3.apodization_dim3_combobox_selection == 4):
                        apodization_line += ' GMB'
                        apodization_line += ' -lb ' + str(self.tabDim3.a_dim3) + ' -gb ' + str(self.tabDim3.b_dim3) + ' -c ' + str(self.tabDim3.apodization_first_point_scaling_dim3)
                    elif(self.tabDim3.apodization_dim3_combobox_selection == 5):
                        apodization_line += ' TM'
                        apodization_line += ' -t1 ' + str(self.tabDim3.t1_dim3) + ' -t2 ' + str(self.tabDim3.t2_dim3) + ' -c ' + str(self.tabDim3.apodization_first_point_scaling_dim3)
                    elif(self.tabDim3.apodization_dim3_combobox_selection == 6):
                        apodization_line += ' TRI'
                        apodization_line += ' -loc ' + str(self.tabDim3.loc_dim3) + ' -c ' + str(self.tabDim3.apodization_first_point_scaling_dim3)
                    nmrproc_com.write(apodization_line + ' \\\n')

            # Check to see if zero filling is selected
            if(self.tabDim3.zero_filling_checkbox_dim3.GetValue() == True):
                zero_filling_line = '| nmrPipe -fn ZF'
                if(self.tabDim3.zero_filling_dim3_combobox_selection == 0):
                    zero_filling_line += ' -zf ' + str(self.tabDim3.zero_filling_dim3_value_doubling_times)
                if(self.tabDim3.zero_filling_dim3_combobox_selection == 1):
                    zero_filling_line += ' -pad ' + str(self.tabDim3.zero_filling_dim3_value_zeros_to_add)
                elif(self.tabDim3.zero_filling_dim3_combobox_selection == 2):
                    zero_filling_line += ' -size ' + str(self.tabDim3.zero_filling_dim3_value_final_data_size)
                if(self.tabDim3.zero_filling_round_checkbox_dim3.GetValue() == True):
                    zero_filling_line += ' -auto'

                nmrproc_com.write(zero_filling_line + ' \\\n')

            # Check to see if the fourier transform checkbox is checked
            if(self.tabDim3.fourier_transform_checkbox_dim3.GetValue() == True):
                fourier_transform_line = '| nmrPipe -fn FT'
                if(self.tabDim3.ft_method_selection_dim3 == 0):
                    fourier_transform_line += ' -auto'
                elif(self.tabDim3.ft_method_selection_dim3 == 1):
                    fourier_transform_line += ' -real'
                elif(self.tabDim3.ft_method_selection_dim3 == 2):
                    fourier_transform_line += ' -inv'
                elif(self.tabDim3.ft_method_selection_dim3 == 3):
                    fourier_transform_line += ' -alt'
                elif(self.tabDim3.ft_method_selection_dim3 == 4):
                    fourier_transform_line += ' -neg'

                nmrproc_com.write(fourier_transform_line + ' \\\n')
            
            # Check to see if the phase correction checkbox is checked
            if(self.tabDim3.phase_correction_checkbox_dim3.GetValue() == True):
                phase_correction_line = '| nmrPipe -fn PS'
                phase_correction_line += ' -p0 ' + str(self.tabDim3.phase_correction_p0_textcontrol_dim3.GetValue()) + ' -p1 ' + str(self.tabDim3.phase_correction_p1_textcontrol_dim3.GetValue()) + ' -di '

                nmrproc_com.write(phase_correction_line + '\\\n')
            
            # Check to see if the extraction checkbox is checked
            if(self.tabDim3.extraction_checkbox_dim3.GetValue() == True):
                extraction_line = '| nmrPipe -fn EXT'
                extraction_line += ' -x1 ' + str(self.tabDim3.extraction_ppm_start_textcontrol_dim3.GetValue()) + 'ppm'+ ' -xn ' + str(self.tabDim3.extraction_ppm_end_textcontrol_dim3.GetValue()) + 'ppm -sw '

                nmrproc_com.write(extraction_line + ' \\\n')
            
            # Check to see if the baseline correction checkbox is checked
            if(self.tabDim3.baseline_correction_checkbox_dim3.GetValue() == True):
                if(self.tabDim3.baseline_correction_radio_box_selection_dim3 == 0):
                    baseline_correction_line = '| nmrPipe -fn BASE'
                    # add the node width and node list to the baseline correction line
                    baseline_correction_line += ' -nw ' + str(self.tabDim3.baseline_correction_nodes_textcontrol_dim3.GetValue()) + ' -nl '
                    for node in node_list_final_dim3:
                        baseline_correction_line += str(node) + '% '
                    nmrproc_com.write(baseline_correction_line + '\\\n')
                elif(self.tabDim3.baseline_correction_radio_box_selection_dim3 == 1):
                    baseline_correction_line = '| nmrPipe -fn POLY'
                    # add the node width and node list to the baseline correction line
                    baseline_correction_line += ' -nw ' + str(self.tabDim3.baseline_correction_nodes_textcontrol_dim3.GetValue()) + ' -nl '
                    for node in node_list_final_dim3:
                        baseline_correction_line += str(node) + '% '
                    # add the polynomial order to the baseline correction line
                    baseline_correction_line += ' -ord ' + str(self.tabDim3.baseline_correction_polynomial_order_textcontrol_dim3.GetValue())
                    nmrproc_com.write(baseline_correction_line + ' \\\n')
        

        
        # if('/' not in self.nmr_data.fid_file):
        if(self.nmr_data.dim == 2 and self.nmr_data.pseudo_axis==False):
            nmrproc_com.write(' -ov -out test.ft2\n')
        elif(self.nmr_data.dim == 2 and self.nmr_data.pseudo_axis==True):
            nmrproc_com.write(' -ov -out test.ft\n')
        elif(self.nmr_data.dim == 3 and self.nmr_data.pseudo_axis==True):
            nmrproc_com.write(' -ov -out test.ft2\n')
            nmrproc_com.write('proj3D.tcl -in test.ft2')
        elif(self.nmr_data.dim == 1):
            nmrproc_com.write(' -ov -out test.ft\n')
        elif(self.nmr_data.dim == 3 and self.nmr_data.pseudo_axis==False):
            nmrproc_com.write(' -ov -out test.ft3\n')
            nmrproc_com.write('proj3D.tcl -in test.ft3')
        # else:
        #     nmrproc_com.write('| nmrPipe  -fn TP -auto \\\n')
        #     nmrproc_com.write('> ./XYZA/t2.ft3 \n')
        #     nmrproc_com.write('nmrPipe -in ./XYZA/t2.ft3  \\\n')
        #     nmrproc_com.write('| pipe2xyz -out ft/test%03d.ft3 -x  \\\n')
        #     nmrproc_com.write('| proj3D.tcl -in ft/test%03d.ft3  \\\n')
        #     nmrproc_com.write('mv ./XYZA/t2.ft3 ./test.ft3  \\\n')






        self.save_processing_parameters()

        self.change_to_cwd()


    def change_to_path(self):
        if(self.parent.path != ''):
            os.chdir(self.parent.path)
        else:
            if(self.parent.original_frame != None):
                if(self.parent.original_frame.parent.path != None):
                    try:
                        os.chdir(self.parent.original_frame.parent.path)
                    except:
                        pass
            if(self.parent.file_parser == True):
                os.chdir(self.parent.path)

    def change_to_path_run(self):
        if(self.parent.path != ''): 
            os.chdir(self.parent.path)
        else:
            if(self.parent.original_frame != None):
                try:
                    self.parent.original_frame.Disable()
                except:
                    pass
                if(self.parent.original_frame.parent.path !=''):
                    os.chdir(self.parent.original_frame.parent.path)
            if(self.parent.file_parser == True):
                os.chdir(self.parent.path)

    def change_to_cwd(self):
        if(self.parent.cwd != ''):
            os.chdir(self.parent.cwd)
        else:
            # Change path if using unidecFile parser
            if(self.parent.original_frame != None):
                try:
                    if(self.parent.original_frame.parent.cwd != None):
                        os.chdir(self.parent.original_frame.parent.cwd)
                except:
                    pass

            if(self.parent.file_parser == True):
                # Change the path to the path of the original file
                os.chdir(self.parent.cwd)

    
    def save_processing_parameters(self):
        # Save all the processing parameters to a processing file
        processing_file = open('processing_parameters.txt', 'w')
        processing_file.write('Processing Parameters\n\n')
        processing_file.write('Dimension 1\n')
        processing_file.write('Solvent Suppression: {}\n'.format(self.tabDim1.solvent_suppression_checkbox.GetValue()))
        processing_file.write('Filter Selection: {}\n'.format(self.tabDim1.solvent_suppression_filter_selection))
        processing_file.write('Lowpass Shape Selection: {}\n'.format(self.tabDim1.solvent_suppression_lowpass_shape_selection))
        processing_file.write('Linear Prediction: {}\n'.format(self.tabDim1.linear_prediction_checkbox.GetValue()))
        processing_file.write('Linear Prediction Options Selection: {}\n'.format(self.tabDim1.linear_prediction_options_selection))
        processing_file.write('Linear Prediction Coefficients Selection: {}\n'.format(self.tabDim1.linear_prediction_coefficients_selection))
        processing_file.write('Apodization: {}\n'.format(self.tabDim1.apodization_checkbox.GetValue()))
        processing_file.write('Apodization Combobox Selection: {}\n'.format(self.tabDim1.apodization_combobox_selection))
        processing_file.write('Exponential Line Broadening: {}\n'.format(self.tabDim1.exponential_line_broadening))
        processing_file.write('Apodization First Point Scaling: {}\n'.format(self.tabDim1.apodization_first_point_scaling))
        processing_file.write('G1: {}\n'.format(self.tabDim1.g1))
        processing_file.write('G2: {}\n'.format(self.tabDim1.g2))
        processing_file.write('G3: {}\n'.format(self.tabDim1.g3))
        processing_file.write('Offset: {}\n'.format(self.tabDim1.offset))
        processing_file.write('End: {}\n'.format(self.tabDim1.end))
        processing_file.write('Power: {}\n'.format(self.tabDim1.power))
        processing_file.write('A: {}\n'.format(self.tabDim1.a))
        processing_file.write('B: {}\n'.format(self.tabDim1.b))
        processing_file.write('T1: {}\n'.format(self.tabDim1.t1))
        processing_file.write('T2: {}\n'.format(self.tabDim1.t2))
        processing_file.write('Loc: {}\n'.format(self.tabDim1.loc))
        processing_file.write('Zero Filling: {}\n'.format(self.tabDim1.zero_filling_checkbox.GetValue()))
        processing_file.write('Zero Filling Combobox Selection: {}\n'.format(self.tabDim1.zero_filling_combobox_selection))
        processing_file.write('Zero Filling Value Doubling Times: {}\n'.format(self.tabDim1.zero_filling_value_doubling_times))
        processing_file.write('Zero Filling Value Zeros to Add: {}\n'.format(self.tabDim1.zero_filling_value_zeros_to_add))
        processing_file.write('Zero Filling Value Final Data Size: {}\n'.format(self.tabDim1.zero_filling_value_final_data_size))
        processing_file.write('Zero Filling Round Checkbox: {}\n'.format(self.tabDim1.zero_filling_round_checkbox.GetValue()))
        processing_file.write('Fourier Transform: {}\n'.format(self.tabDim1.fourier_transform_checkbox.GetValue()))
        processing_file.write('Fourier Transform Method Selection: {}\n'.format(self.tabDim1.ft_method_selection))
        processing_file.write('Phase Correction: {}\n'.format(self.tabDim1.phase_correction_checkbox.GetValue()))
        processing_file.write('Phase Correction P0: {}\n'.format(self.tabDim1.phase_correction_p0_textcontrol.GetValue()))
        processing_file.write('Phase Correction P1: {}\n'.format(self.tabDim1.phase_correction_p1_textcontrol.GetValue()))
        processing_file.write('Magnitude Mode: {}\n'.format(self.tabDim1.magnitude_mode_checkbox.GetValue()))
        processing_file.write('Extraction: {}\n'.format(self.tabDim1.extraction_checkbox.GetValue()))
        processing_file.write('Extraction PPM Start: {}\n'.format(self.tabDim1.extraction_ppm_start_textcontrol.GetValue()))
        processing_file.write('Extraction PPM End: {}\n'.format(self.tabDim1.extraction_ppm_end_textcontrol.GetValue()))
        processing_file.write('Baseline Correction: {}\n'.format(self.tabDim1.baseline_correction_checkbox.GetValue()))
        processing_file.write('Baseline Correction Radio Box Selection: {}\n'.format(self.tabDim1.baseline_correction_radio_box_selection))
        processing_file.write('Baseline Correction Nodes: {}\n'.format(self.tabDim1.baseline_correction_nodes_textcontrol.GetValue()))
        processing_file.write('Baseline Correction Node List: {}\n'.format(self.tabDim1.baseline_correction_node_list_textcontrol.GetValue()))
        processing_file.write('Baseline Correction Polynomial Order: {}\n'.format(self.tabDim1.baseline_correction_polynomial_order_textcontrol.GetValue()))
        processing_file.write('\n\n')

        include_dim2 = False
        include_dim3 = False
        if(self.nmr_data.dim == 2 and self.nmr_data.pseudo_axis==False):
            include_dim2 = True
        if(self.nmr_data.dim == 3 and self.nmr_data.pseudo_axis==True):
            include_dim2 = True
            include_dim3 = False
        if(self.nmr_data.dim == 3 and self.nmr_data.pseudo_axis==False):
            include_dim2 = True
            include_dim3 = True
        
        if(include_dim2 == True):
            processing_file.write('Dimension 2\n')
            processing_file.write('Linear Prediction: {}\n'.format(self.tabDim2.linear_prediction_radio_box_dim2.GetSelection()))
            processing_file.write('Linear Prediction Options Selection: {}\n'.format(self.tabDim2.linear_prediction_dim2_options_selection))
            processing_file.write('Linear Prediction Coefficients Selection: {}\n'.format(self.tabDim2.linear_prediction_dim2_coefficients_selection))
            processing_file.write('NUS file: {}\n'.format(self.tabDim2.nuslist_name_dim2))
            processing_file.write('NUS CPU: {}\n'.format(self.tabDim2.number_of_nus_CPU_dim2))
            processing_file.write('NUS Iterations: {}\n'.format(self.tabDim2.nus_iterations_dim2))
            processing_file.write('Apodization: {}\n'.format(self.tabDim2.apodization_checkbox_dim2.GetValue()))
            processing_file.write('Apodization Combobox Selection: {}\n'.format(self.tabDim2.apodization_dim2_combobox_selection))
            processing_file.write('Exponential Line Broadening: {}\n'.format(self.tabDim2.exponential_line_broadening_dim2))
            processing_file.write('Apodization First Point Scaling: {}\n'.format(self.tabDim2.apodization_first_point_scaling_dim2))
            processing_file.write('G1: {}\n'.format(self.tabDim2.g1_dim2))
            processing_file.write('G2: {}\n'.format(self.tabDim2.g2_dim2))
            processing_file.write('G3: {}\n'.format(self.tabDim2.g3_dim2))
            processing_file.write('Offset: {}\n'.format(self.tabDim2.offset_dim2))
            processing_file.write('End: {}\n'.format(self.tabDim2.end_dim2))
            processing_file.write('Power: {}\n'.format(self.tabDim2.power_dim2))
            processing_file.write('A: {}\n'.format(self.tabDim2.a_dim2))
            processing_file.write('B: {}\n'.format(self.tabDim2.b_dim2))
            processing_file.write('T1: {}\n'.format(self.tabDim2.t1_dim2))
            processing_file.write('T2: {}\n'.format(self.tabDim2.t2_dim2))
            processing_file.write('Loc: {}\n'.format(self.tabDim2.loc_dim2))
            processing_file.write('Zero Filling: {}\n'.format(self.tabDim2.zero_filling_checkbox_dim2.GetValue()))
            processing_file.write('Zero Filling Combobox Selection: {}\n'.format(self.tabDim2.zero_filling_dim2_combobox_selection))
            processing_file.write('Zero Filling Value Doubling Times: {}\n'.format(self.tabDim2.zero_filling_dim2_value_doubling_times))
            processing_file.write('Zero Filling Value Zeros to Add: {}\n'.format(self.tabDim2.zero_filling_dim2_value_zeros_to_add))
            processing_file.write('Zero Filling Value Final Data Size: {}\n'.format(self.tabDim2.zero_filling_dim2_value_final_data_size))
            processing_file.write('Zero Filling Round Checkbox: {}\n'.format(self.tabDim2.zero_filling_round_checkbox_dim2.GetValue()))
            processing_file.write('Fourier Transform: {}\n'.format(self.tabDim2.fourier_transform_checkbox_dim2.GetValue()))
            processing_file.write('Fourier Transform Method Selection: {}\n'.format(self.tabDim2.ft_method_selection_dim2))
            processing_file.write('Phase Correction: {}\n'.format(self.tabDim2.phase_correction_checkbox_dim2.GetValue()))
            processing_file.write('Phase Correction P0: {}\n'.format(self.tabDim2.phase_correction_p0_textcontrol_dim2.GetValue()))
            processing_file.write('Phase Correction P1: {}\n'.format(self.tabDim2.phase_correction_p1_textcontrol_dim2.GetValue()))
            processing_file.write('F1180: {}\n'.format(self.tabDim2.phase_correction_f1180_button_dim2.GetValue()))
            processing_file.write('Extraction: {}\n'.format(self.tabDim2.extraction_checkbox_dim2.GetValue()))
            processing_file.write('Extraction PPM Start: {}\n'.format(self.tabDim2.extraction_ppm_start_textcontrol_dim2.GetValue()))
            processing_file.write('Extraction PPM End: {}\n'.format(self.tabDim2.extraction_ppm_end_textcontrol_dim2.GetValue()))
            processing_file.write('Baseline Correction: {}\n'.format(self.tabDim2.baseline_correction_checkbox_dim2.GetValue()))
            processing_file.write('Baseline Correction Radio Box Selection: {}\n'.format(self.tabDim2.baseline_correction_radio_box_selection_dim2))
            processing_file.write('Baseline Correction Nodes: {}\n'.format(self.tabDim2.baseline_correction_nodes_textcontrol_dim2.GetValue()))
            processing_file.write('Baseline Correction Node List: {}\n'.format(self.tabDim2.baseline_correction_node_list_textcontrol_dim2.GetValue()))
            processing_file.write('Baseline Correction Polynomial Order: {}\n'.format(self.tabDim2.baseline_correction_polynomial_order_textcontrol_dim2.GetValue()))
            processing_file.write('\n\n')


        if(include_dim3 == True):
            processing_file.write('Dimension 3\n')
            processing_file.write('Linear Prediction: {}\n'.format(self.tabDim3.linear_prediction_radio_box_dim3.GetSelection()))
            processing_file.write('Linear Prediction Options Selection: {}\n'.format(self.tabDim3.linear_prediction_dim3_options_selection))
            processing_file.write('Linear Prediction Coefficients Selection: {}\n'.format(self.tabDim3.linear_prediction_dim3_coefficients_selection))
            processing_file.write('NUS file: {}\n'.format(self.tabDim2.nuslist_name_dim2))
            processing_file.write('NUS CPU: {}\n'.format(self.tabDim2.number_of_nus_CPU_dim2))
            processing_file.write('NUS Iterations: {}\n'.format(self.tabDim2.nus_iterations_dim2))
            processing_file.write('Apodization: {}\n'.format(self.tabDim3.apodization_checkbox_dim3.GetValue()))
            processing_file.write('Apodization Combobox Selection: {}\n'.format(self.tabDim3.apodization_dim3_combobox_selection))
            processing_file.write('Exponential Line Broadening: {}\n'.format(self.tabDim3.exponential_line_broadening_dim3))
            processing_file.write('Apodization First Point Scaling: {}\n'.format(self.tabDim3.apodization_first_point_scaling_dim3))
            processing_file.write('G1: {}\n'.format(self.tabDim3.g1_dim3))
            processing_file.write('G2: {}\n'.format(self.tabDim3.g2_dim3))
            processing_file.write('G3: {}\n'.format(self.tabDim3.g3_dim3))
            processing_file.write('Offset: {}\n'.format(self.tabDim3.offset_dim3))
            processing_file.write('End: {}\n'.format(self.tabDim3.end_dim3))
            processing_file.write('Power: {}\n'.format(self.tabDim3.power_dim3))
            processing_file.write('A: {}\n'.format(self.tabDim3.a_dim3))
            processing_file.write('B: {}\n'.format(self.tabDim3.b_dim3))
            processing_file.write('T1: {}\n'.format(self.tabDim3.t1_dim3))
            processing_file.write('T2: {}\n'.format(self.tabDim3.t2_dim3))
            processing_file.write('Loc: {}\n'.format(self.tabDim3.loc_dim3))
            processing_file.write('Zero Filling: {}\n'.format(self.tabDim3.zero_filling_checkbox_dim3.GetValue()))
            processing_file.write('Zero Filling Combobox Selection: {}\n'.format(self.tabDim3.zero_filling_dim3_combobox_selection))
            processing_file.write('Zero Filling Value Doubling Times: {}\n'.format(self.tabDim3.zero_filling_dim3_value_doubling_times))
            processing_file.write('Zero Filling Value Zeros to Add: {}\n'.format(self.tabDim3.zero_filling_dim3_value_zeros_to_add))
            processing_file.write('Zero Filling Value Final Data Size: {}\n'.format(self.tabDim3.zero_filling_dim3_value_final_data_size))
            processing_file.write('Zero Filling Round Checkbox: {}\n'.format(self.tabDim3.zero_filling_round_checkbox_dim3.GetValue()))
            processing_file.write('Fourier Transform: {}\n'.format(self.tabDim3.fourier_transform_checkbox_dim3.GetValue()))
            processing_file.write('Fourier Transform Method Selection: {}\n'.format(self.tabDim3.ft_method_selection_dim3))
            processing_file.write('Phase Correction: {}\n'.format(self.tabDim3.phase_correction_checkbox_dim3.GetValue()))
            processing_file.write('Phase Correction P0: {}\n'.format(self.tabDim3.phase_correction_p0_textcontrol_dim3.GetValue()))
            processing_file.write('Phase Correction P1: {}\n'.format(self.tabDim3.phase_correction_p1_textcontrol_dim3.GetValue()))
            processing_file.write('F1180: {}\n'.format(self.tabDim3.phase_correction_f1180_button_dim3.GetValue()))
            processing_file.write('Extraction: {}\n'.format(self.tabDim3.extraction_checkbox_dim3.GetValue()))
            processing_file.write('Extraction PPM Start: {}\n'.format(self.tabDim3.extraction_ppm_start_textcontrol_dim3.GetValue()))
            processing_file.write('Extraction PPM End: {}\n'.format(self.tabDim3.extraction_ppm_end_textcontrol_dim3.GetValue()))
            processing_file.write('Baseline Correction: {}\n'.format(self.tabDim3.baseline_correction_checkbox_dim3.GetValue()))
            processing_file.write('Baseline Correction Radio Box Selection: {}\n'.format(self.tabDim3.baseline_correction_radio_box_selection_dim3))
            processing_file.write('Baseline Correction Nodes: {}\n'.format(self.tabDim3.baseline_correction_nodes_textcontrol_dim3.GetValue()))
            processing_file.write('Baseline Correction Node List: {}\n'.format(self.tabDim3.baseline_correction_node_list_textcontrol_dim3.GetValue()))
            processing_file.write('Baseline Correction Polynomial Order: {}\n'.format(self.tabDim3.baseline_correction_polynomial_order_textcontrol_dim3.GetValue()))
            processing_file.write('\n\n')





    def on_show_nmrproc_com(self, event):
        if(self.parent.path != ''):
            os.chdir(self.parent.path)
        elif(self.parent.file_parser == True):
            os.chdir(self.parent.path)
        # Create a popout window showing the nmrproc.com file
        try:
            nmrproc_com = open('nmrproc.com', 'r')
            lines = nmrproc_com.readlines()
            nmrproc_com.close()
        except:
            dlg = wx.MessageDialog(self, 'The nmrproc.com file cannot be found in the current directory', 'Warning', wx.OK | wx.ICON_WARNING)
            self.Raise()
            self.SetFocus()
            result = dlg.ShowModal()
            return

        nmrproc_com_text = ''
        for line in lines:
            nmrproc_com_text += line

        nmrproc_com_window = wx.Frame(self, -1, 'NMRPipe Processing Script', size=(500, 500))
        if(darkdetect.isDark() == True and platform!='windows'):
            nmrproc_com_window.SetBackgroundColour((53, 53, 53, 255))
        else:
            nmrproc_com_window.SetBackgroundColour('White')

        nmrproc_com_textcontrol = wx.TextCtrl(nmrproc_com_window, -1, nmrproc_com_text, style=wx.TE_MULTILINE | wx.TE_READONLY)
        nmrproc_com_sizer = wx.BoxSizer(wx.VERTICAL)
        nmrproc_com_sizer.Add(nmrproc_com_textcontrol, 1, wx.EXPAND)
        nmrproc_com_window.SetSizer(nmrproc_com_sizer)
        nmrproc_com_window.Show()

        if(self.parent.cwd != ''):
            os.chdir(self.parent.cwd)
        elif(self.parent.file_parser == True):
            os.chdir(self.parent.cwd)


    def on_run_processing(self, event):
        try:
            if(self.tabDim2.linear_prediction_radio_box_dim2.GetSelection() == 2):
                # SMILE processing is selected, asking the user to confirm SMILE is installed as part of nmrPipe
                dlg = wx.MessageDialog(self, 'SMILE processing is selected. Ensure that SMILE is installed as part of nmrPipe', 'Warning', wx.OK | wx.CANCEL | wx.ICON_WARNING)
                self.Raise()
                self.SetFocus()
                result = dlg.ShowModal()
                if(result == wx.ID_CANCEL):
                    self.change_to_cwd()
                    return
                
                if(self.tabDim1.extraction_checkbox.GetValue()==False):
                    dlg = wx.MessageDialog(self, 'No direct dimension data extraction is selected, SMILE reconstruction may take a while. Consider extracting a region of the direct dimension before reconstruction. Do you want to continue or cancel?', 'Warning', wx.OK | wx.CANCEL | wx.ICON_WARNING)
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    if(result == wx.ID_CANCEL):
                        self.change_to_cwd()
                        return
        except:
            pass

        if(platform=='windows'):
            self.on_run_processing_windows(event)
        else:
            self.on_run_processing_nmrproc(event)

    def on_run_processing_windows(self, event):
        # Disable spinview while reprocessing the data
        self.change_to_path_run()
        
        # Apply the processing parameters to the data
        self.apply_processing_parameters()
        


    def apply_processing_parameters(self):
        # Process the data according to the user inputted processing parameters
        data = self.nmr_data.data
        dic = self.nmr_data.dic
        if(self.nmr_data.dim == 1):
            nmrfile = 'test.ft'
        elif(self.nmr_data.dim == 2):
            nmrfile = 'test.ft2'
            # If NUS reconstruction is selected, say that this is not possible in windows, please use a machine containing nmrPipe
            if(self.nmr_data.pseudo_axis == False):
                if(self.tabDim2.linear_prediction_radio_box_dim2_selection == 2):
                    dlg = wx.MessageDialog(self, 'NUS reconstruction is not possible on Windows machines. Please use a machine containing nmrPipe.', 'Warning', wx.OK | wx.ICON_WARNING)
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    dlg.Destroy()
                    return
        else:
            if(self.nmr_data.pseudo_axis==False):
                # 3D processing is not supported for windows processing, please use a machine containing nmrPipe
                dlg = wx.MessageDialog(self, '3D processing is not supported for windows processing. Please use a machine containing nmrPipe.', 'Warning', wx.OK | wx.ICON_WARNING)
                self.Raise()
                self.SetFocus()
                result = dlg.ShowModal()
                dlg.Destroy()
                return

        
        if(self.nmr_data.dim == 1):
            # Process the first dimension
            dic,data = self.process_dimension_1(dic,data,dim=0)

        elif(self.nmr_data.dim == 2 and self.nmr_data.pseudo_axis==True):
            # Process the first dimension of the pseudo2D spectrum
            dic,data = self.process_dimension_1(dic,data,dim=1)
        elif(self.nmr_data.dim == 2 and self.nmr_data.pseudo_axis==False):
            # Process the first dimension
            dic,data = self.process_dimension_1(dic,data)
            # Process the second dimension
            dic,data = self.process_dimension_2(dic,data)
        elif(self.nmr_data.dim == 3 and self.nmr_data.pseudo_axis==True):
            dic,data = self.process_dimension_1(dic,data,dim=2)
            # Find out which one of the axes is chemical shift and which one is pseudo


        
        ng.pipe.write(nmrfile, dic, data, overwrite=True)
        original_frame = []
        if(self.parent.original_frame != None):
            original_frame = True
            self.parent.original_frame.Enable()
            path = self.parent.original_frame.parent.path
            cwd = self.parent.original_frame.parent.cwd
            self.parent.original_frame.parent.reprocess=True
            self.parent.original_frame.parent.Close()
            if(self.parent.original_frame.parent.path !=''):
                os.chdir(self.parent.original_frame.parent.path)
            from SpinExplorer.SpinView import MyApp
            app = MyApp()
            if(self.parent.original_frame.parent.cwd !=''):
                app.path = path
                app.cwd = cwd
                
        # Check to see if the output file exists
        if(os.path.exists(nmrfile) == False):
            if(original_frame==True):
                if(app.cwd !=''):
                    os.chdir(app.cwd)
            elif(self.parent.cwd !=''):
                os.chdir(self.parent.cwd)
            message = 'The processing spectrum file ({}) file cannot be found in the current directory. Processing unsuccessful. Ensure that nmrPipe has been downloaded and added to the path.'.format(nmrfile)
            dlg = wx.MessageDialog(self, message, 'Warning', wx.OK | wx.ICON_WARNING)
            self.Raise()
            self.SetFocus()
            result = dlg.ShowModal()
            
                
            return
        else:
            if(original_frame==True):
                if(app.cwd !=''):
                    os.chdir(app.cwd)
            elif(self.parent.cwd !=''):
                os.chdir(self.parent.cwd)
            message = 'Processing successful. The processed spectrum file ({}) has been created in the current directory.'.format(nmrfile)
            dlg = wx.MessageDialog(self, message, 'Success', wx.OK | wx.ICON_INFORMATION)
            self.Raise()
            self.SetFocus()
            result = dlg.ShowModal()
            return
        
    def process_dimension_1(self, dic, data,dim):
        # Process the first dimension
        if(self.tabDim1.solvent_suppression_checkbox.GetValue() == True):
            # Apply solvent suppression
            data_orgiginal = data
            if(self.tabDim1.solvent_suppression_filter_selection == 0):
                try:
                    dic,data = ng.pipe_proc.sol(dic,data_orgiginal, mode = 'low', fs=int(self.tabDim1.solvent_suppression_lowpass_shape_selection)+1)
                except:
                    try:
                        fl = 16
                        fs = int(self.tabDim1.solvent_suppression_lowpass_shape_selection)+1
                        import scipy.signal.windows
                        if(fs == 1):
                            filter = scipy.signal.windows.boxcar(33)
                            data =  ng.proc_bl.sol_general(data, filter, w=33, mode='same')
                        else:
                            # Give a message saying the solvent suppression did not work correctly, the spectrum was processed without solvent suppression
                            dlg = wx.MessageDialog(self, 'The solvent suppression filter did not work correctly. Continuing without digital solvent suppression.', 'Warning', wx.OK | wx.ICON_WARNING)
                            self.Raise()
                            self.SetFocus()
                            result = dlg.ShowModal()
                            data = data_orgiginal
                    except:
                        # Give a message saying the solvent suppression did not work correctly, the spectrum was processed without solvent suppression
                        dlg = wx.MessageDialog(self, 'The solvent suppression filter did not work correctly. Continuing without digital solvent suppression.', 'Warning', wx.OK | wx.ICON_WARNING)
                        self.Raise()
                        self.SetFocus()
                        result = dlg.ShowModal()
                        data = data_orgiginal

            else:
                # Give an error saying that the selected solvent suppression filter is not supported for windows processing, please use a machine containing nmrPipe
                dlg = wx.MessageDialog(self, 'The selected solvent suppression filter is not supported for windows processing. Please change to low bandpass filter or use a machine containing nmrPipe.', 'Warning', wx.OK | wx.ICON_WARNING)
                self.Raise()
                self.SetFocus()
                result = dlg.ShowModal()
                return
        if(self.tabDim1.linear_prediction_checkbox.GetValue() == True):
            # Apply linear prediction
            if(self.tabDim1.linear_prediction_options_selection == 0):
                if(self.tabDim1.linear_prediction_options_selection == 0):
                    append = 'after'
                else:
                    append = 'before'
                if(self.tabDim1.linear_prediction_coefficients_selection == 0):
                    mode = 'f'
                elif(self.tabDim1.linear_prediction_coefficients_selection == 1):
                    mode = 'b'
                else:
                    mode = 'fb'
                dic,data = ng.pipe_proc.lp(dic,data, pred='default',mode=mode, append=append)
            
        if(self.tabDim1.apodization_checkbox.GetValue() == True):
            # Apply apodization
            if(self.tabDim1.apodization_combobox_selection == 0):
                dic,data = ng.pipe_proc.em(dic,data, lb=0.0, c=float(self.tabDim1.apodization_first_point_scaling))
            elif(self.tabDim1.apodization_combobox_selection == 1):
                dic,data = ng.pipe_proc.em(dic,data, lb=float(self.tabDim1.exponential_line_broadening), c=float(self.tabDim1.apodization_first_point_scaling))
            elif(self.tabDim1.apodization_combobox_selection == 2):
                dic,data = ng.pipe_proc.gm(dic,data, g1=float(self.tabDim1.g1), g2=float(self.tabDim1.g2), g3=float(self.tabDim1.g3), c=float(self.tabDim1.apodization_first_point_scaling))
            elif(self.tabDim1.apodization_combobox_selection == 3):
                dic,data = ng.pipe_proc.sp(dic,data, off=float(self.tabDim1.offset), end=float(self.tabDim1.end), pow=int(self.tabDim1.power), c=float(self.tabDim1.apodization_first_point_scaling))
            elif(self.tabDim1.apodization_combobox_selection == 4):
                dic,data = ng.pipe_proc.gmb(dic,data, a=float(self.tabDim1.a), b=float(self.tabDim1.b), c=float(self.tabDim1.apodization_first_point_scaling))
            elif(self.tabDim1.apodization_combobox_selection == 5):
                dic,data = ng.pipe_proc.tp(dic,data, t1=float(self.tabDim1.t1), t2=float(self.tabDim1.t2), c=float(self.tabDim1.apodization_first_point_scaling))
            elif(self.tabDim1.apodization_combobox_selection == 6):
                dic,data = ng.pipe_proc.tri(dic,data, loc=float(self.tabDim1.loc), c=float(self.tabDim1.apodization_first_point_scaling))


        if(self.tabDim1.zero_filling_checkbox.GetValue() == True):
            if(self.tabDim1.zero_filling_round_checkbox.GetValue() == True):
                round = True
            else:
                round = False

            if(self.tabDim1.zero_filling_combobox_selection == 0):
                dic,data = ng.pipe_proc.zf(dic,data, zf=int(self.tabDim1.zero_filling_value_doubling_times), auto=round)
            elif(self.tabDim1.zero_filling_combobox_selection == 1):
                dic,data = ng.pipe_proc.zf(dic,data, pad=int(self.tabDim1.zero_filling_value_zeros_to_add), auto=round)
            elif(self.tabDim1.zero_filling_combobox_selection == 2):
                dic,data = ng.pipe_proc.zf(dic,data, size=int(self.tabDim1.zero_filling_value_final_data_size), auto=round)

            dic['FDF2TDSIZE'] = data.T.shape[0]


        if(self.tabDim1.fourier_transform_checkbox.GetValue() == True):
            if(self.tabDim1.ft_method_selection == 0):
                dic,data = ng.pipe_proc.ft(dic,data, auto=True)
            elif(self.tabDim1.ft_method_selection == 1):
                dic,data = ng.pipe_proc.ft(dic,data, real=True)
            elif(self.tabDim1.ft_method_selection == 2):
                dic,data = ng.pipe_proc.ft(dic,data, inv=True)
            elif(self.tabDim1.ft_method_selection == 3):
                dic,data = ng.pipe_proc.ft(dic,data, alt=True)

        if(self.tabDim1.phase_correction_checkbox.GetValue() == True):
            dic,data = ng.pipe_proc.ps(dic,data, p0=float(self.tabDim1.phase_correction_p0_textcontrol.GetValue()), p1=float(self.tabDim1.phase_correction_p1_textcontrol.GetValue()))

        if(self.tabDim1.magnitude_mode_checkbox.GetValue()==True):
            dic,data = ng.pipe_proc.mc(dic,data)

        if(self.tabDim1.extraction_checkbox.GetValue() == True):
            # Find the indexes of the ppm values selected
            # Get the ppm values from the data
            ppm_values = ng.pipe.make_uc(dic, data, dim=dim)
            ppm_values = ppm_values.ppm_scale()
            x_initial = np.abs(ppm_values - float(self.tabDim1.extraction_ppm_start_textcontrol.GetValue())).argmin()
            x_final = np.abs(ppm_values - float(self.tabDim1.extraction_ppm_end_textcontrol.GetValue())).argmin()
            if(x_initial > x_final):
                x_initial, x_final = x_final, x_initial
            # Change x_initial and x_final so that the difference is an even number
            if((x_final - x_initial + 1) % 2 != 0):
                x_final += 1
            dic,data = ng.pipe_proc.ext(dic,data, x1=x_initial, xn=x_final, sw=True)


        if(self.tabDim1.baseline_correction_checkbox.GetValue() == True):
            if(self.tabDim1.baseline_correction_radio_box_selection == 1):
                # If POLY baseline correction is selected, this is not currently supported on windows without nmrPipe
                message = 'The selected baseline correction method is not supported for windows processing. Please use a machine containing nmrPipe or use a linear baselining method.'
                dlg = wx.MessageDialog(self, message, 'Warning', wx.OK | wx.ICON_WARNING)
                self.Raise()
                self.SetFocus()
                result = dlg.ShowModal()
                return
            
            # Split the node list
            node_list = self.tabDim1.baseline_correction_node_list_textcontrol.GetValue()
            node_list = node_list.split(',')
            node_list_final = []
            for node in node_list:
                node_list_final.append(float(node))


            # Convert nodes into points
            node_list_final = np.array(node_list_final)
            node_list_final = (node_list_final/100)*len(data)
            node_list_final = node_list_final.astype(int)
            # Replace any zeros with a number greater than 1 to allow the nmrglue baselining routines to work correctly
            node_list_final[node_list_final == 0] = int(self.tabDim1.baseline_correction_nodes_textcontrol.GetValue()) + 1
            
            dic, data = ng.pipe_proc.base(dic, data, nl = node_list_final, nw = int(self.tabDim1.baseline_correction_nodes_textcontrol.GetValue()))

        dic['FDF2TDSIZE'] = data.T.shape[0]
        
        return dic, data
    

    def process_dimension_2(self, dic, data):
        # Transpose to the second dimension
        dic, data = ng.pipe_proc.tp(dic, data)

        # Process the second dimension
        if(self.tabDim2.linear_prediction_radio_box_dim2_selection == 1):
            # Apply linear prediction
            if(self.tabDim2.linear_prediction_dim2_options_selection == 0):
                append = 'after'
            else:
                append = 'before'
            if(self.tabDim2.linear_prediction_dim2_coefficients_selection == 0):
                mode = 'f'
            elif(self.tabDim2.linear_prediction_dim2_coefficients_selection == 1):
                mode = 'b'
            else:
                mode = 'fb'
            dic,data = ng.pipe_proc.lp(dic,data, pred='default',mode=mode, append=append)
            
        if(self.tabDim2.apodization_checkbox_dim2.GetValue() == True):
            # Apply apodization
            if(self.tabDim2.apodization_dim2_combobox_selection == 0):
                dic,data = ng.pipe_proc.em(dic,data, lb=0.0, c=float(self.tabDim2.apodization_first_point_scaling_dim2))
            elif(self.tabDim2.apodization_dim2_combobox_selection == 1):
                dic,data = ng.pipe_proc.em(dic,data, lb=float(self.tabDim2.exponential_line_broadening_dim2), c=float(self.tabDim2.apodization_first_point_scaling_dim2))
            elif(self.tabDim2.apodization_dim2_combobox_selection == 2):
                dic,data = ng.pipe_proc.gm(dic,data, g1=float(self.tabDim2.g1_dim2), g2=float(self.tabDim2.g2_dim2), g3=float(self.tabDim2.g3_dim2), c=float(self.tabDim2.apodization_first_point_scaling_dim2))
            elif(self.tabDim2.apodization_dim2_combobox_selection == 3):
                dic,data = ng.pipe_proc.sp(dic,data, off=float(self.tabDim2.offset_dim2), end=float(self.tabDim2.end_dim2), pow=int(self.tabDim2.power_dim2), c=float(self.tabDim2.apodization_first_point_scaling_dim2))
            elif(self.tabDim2.apodization_dim2_combobox_selection == 4):
                dic,data = ng.pipe_proc.gmb(dic,data, a=float(self.tabDim2.a_dim2), b=float(self.tabDim2.b_dim2), c=float(self.tabDim2.apodization_first_point_scaling_dim2))
            elif(self.tabDim2.apodization_dim2_combobox_selection == 5):
                dic,data = ng.pipe_proc.tp(dic,data, t1=float(self.tabDim2.t1_dim2), t2=float(self.tabDim2.t2_dim2), c=float(self.tabDim2.apodization_first_point_scaling_dim2))
            elif(self.tabDim2.apodization_dim2_combobox_selection == 6):
                dic,data = ng.pipe_proc.tri(dic,data, loc=float(self.tabDim2.loc_dim2), c=float(self.tabDim2.apodization_first_point_scaling_dim2))


        if(self.tabDim2.zero_filling_checkbox_dim2.GetValue() == True):
            if(self.tabDim2.zero_filling_round_checkbox_dim2.GetValue() == True):
                round = True
            else:
                round = False

            if(self.tabDim2.zero_filling_dim2_combobox_selection == 0):
                dic,data = ng.pipe_proc.zf(dic,data, zf=int(self.tabDim2.zero_filling_dim2_value_doubling_times), auto=round)
            elif(self.tabDim2.zero_filling_dim2_combobox_selection == 1):
                dic,data = ng.pipe_proc.zf(dic,data, pad=int(self.tabDim2.zero_filling_dim2_value_zeros_to_add), auto=round)
            elif(self.tabDim2.zero_filling_dim2_combobox_selection == 2):
                dic,data = ng.pipe_proc.zf(dic,data, size=int(self.tabDim2.zero_filling_dim2_value_final_data_size), auto=round)

        if(self.tabDim2.fourier_transform_checkbox_dim2.GetValue() == True):
            if(self.tabDim2.ft_method_selection_dim2 == 0):
                dic,data = ng.pipe_proc.ft(dic,data, auto=True)
            elif(self.tabDim2.ft_method_selection_dim2 == 1):
                dic,data = ng.pipe_proc.ft(dic,data, real=True)
            elif(self.tabDim2.ft_method_selection_dim2 == 2):
                dic,data = ng.pipe_proc.ft(dic,data, inv=True)
            elif(self.tabDim2.ft_method_selection_dim2 == 3):
                dic,data = ng.pipe_proc.ft(dic,data, alt=True)

        if(self.tabDim2.phase_correction_checkbox_dim2.GetValue() == True):
            dic,data = ng.pipe_proc.ps(dic,data, p0=float(self.tabDim2.phase_correction_p0_textcontrol_dim2.GetValue()), p1=float(self.tabDim2.phase_correction_p1_textcontrol_dim2.GetValue()))

        if(self.tabDim2.extraction_checkbox_dim2.GetValue() == True):
            # Find the indexes of the ppm values selected
            # Get the ppm values from the data
            ppm_values = ng.pipe.make_uc(dic, data,dim=1)
            ppm_values = ppm_values.ppm_scale()
            x_initial = np.abs(ppm_values - float(self.tabDim2.extraction_ppm_start_textcontrol_dim2.GetValue())).argmin()
            x_final = np.abs(ppm_values - float(self.tabDim2.extraction_ppm_end_textcontrol_dim2.GetValue())).argmin()
            if(x_initial > x_final):
                x_initial, x_final = x_final, x_initial

            if((x_final - x_initial + 1) % 2 != 0):
                x_final += 1
            dic,data = ng.pipe_proc.ext(dic,data, x1=x_initial, xn=x_final, sw=True)
            

        if(self.tabDim2.baseline_correction_checkbox_dim2.GetValue() == True):
            if(self.tabDim2.baseline_correction_radio_box_selection_dim2 == 1):
                # If POLY baseline correction is selected, this is not currently supported on windows without nmrPipe
                message = 'The selected baseline correction method is not supported for windows processing. Please use a machine containing nmrPipe or use a linear baselining method.'
                dlg = wx.MessageDialog(self, message, 'Warning', wx.OK | wx.ICON_WARNING)
                self.Raise()
                self.SetFocus()
                result = dlg.ShowModal()
                return
            
            # Split the node list
            node_list = self.tabDim2.baseline_correction_node_list_textcontrol_dim2.GetValue()
            node_list = node_list.split(',')
            node_list_final = []
            for node in node_list:
                node_list_final.append(float(node))


            # Convert nodes into points
            node_list_final = np.array(node_list_final)
            node_list_final = (node_list_final/100)*len(data)
            node_list_final = node_list_final.astype(int)
            # Replace any zeros with a number greater than 1 to allow the nmrglue baselining routines to work correctly
            node_list_final[node_list_final == 0] = int(self.tabDim2.baseline_correction_nodes_textcontrol_dim2.GetValue()) + 1
            
            dic, data = ng.pipe_proc.base(dic, data, nl = node_list_final, nw = int(self.tabDim2.baseline_correction_nodes_textcontrol_dim2.GetValue()))

        
        dic['FDF1TDSIZE'] = data.T.shape[0]
        dic['FDF1FTSIZE'] = data.T.shape[0]
        dic['FDF1QUADFLAG'] = 1.0

        dic, data = ng.pipe_proc.tp(dic, data)

        return dic, data
        

    def on_run_processing_nmrproc(self, event):
        # Disable spinview while reprocessing the data
        self.change_to_path_run()
        
        # Check to see if the nmrproc.com file exists
        if(os.path.exists('nmrproc.com') == False):
            dlg = wx.MessageDialog(self, 'The nmrproc.com file cannot be found in the current directory', 'Warning', wx.OK | wx.ICON_WARNING)
            self.Raise()
            self.SetFocus()
            result = dlg.ShowModal()
            if(self.parent.original_frame != None):
                try:
                    self.parent.original_frame.Disable()
                except:
                    pass
                if(self.parent.original_frame.parent.cwd !=''):
                    os.chdir(self.parent.original_frame.parent.cwd)
            return
        
        # Add execute permissions to the nmrproc.com file
        os.system('chmod +x nmrproc.com')
        
        # Check to see if test.ft file already exists, if it does ask the user if they want to overwrite it
        if(os.path.exists('test.ft')):
            dlg = wx.MessageDialog(self, 'The test.ft file already exists. Do you want to overwrite it?', 'Warning', wx.YES_NO | wx.ICON_WARNING)
            self.Raise()
            self.SetFocus()
            result = dlg.ShowModal()
            if(result == wx.ID_NO):
                if(self.parent.original_frame != None):
                    self.parent.original_frame.Disable()
                    if(self.parent.original_frame.parent.cwd !=''):
                        os.chdir(self.parent.original_frame.parent.cwd)
                return

        
        # Run the nmrproc.com file
        command = 'csh nmrproc.com'

        # Check to see if the output file is not empty
        p = subprocess.Popen(command, stdout=subprocess.PIPE,shell=True)
        p.wait()

        if(self.parent.original_frame != None):        
            if(self.parent.original_frame.parent.cwd !=''):
                os.chdir(self.parent.original_frame.parent.cwd)
        
        if(self.parent.file_parser == True):
            os.chdir(self.parent.cwd)

            
            



        if(self.nmr_data.dim == 1):
            nmrfile = 'test.ft'
        elif(self.nmr_data.dim == 2 and self.nmr_data.pseudo_axis==False):
            nmrfile = 'test.ft2'
        elif(self.nmr_data.dim == 2 and self.nmr_data.pseudo_axis==True):
            nmrfile = 'test.ft'
        elif(self.nmr_data.dim == 3 and self.nmr_data.pseudo_axis==True):
            nmrfile = 'test.ft2'
        elif(self.nmr_data.dim == 3 and self.nmr_data.pseudo_axis==False):
            nmrfile = 'test.ft3'
        else:
            nmrfile = 'test.ft'

        original_frame = []
        if(self.parent.original_frame != None):
            original_frame = True
            self.parent.original_frame.Enable()
            path = self.parent.original_frame.parent.path
            cwd = self.parent.original_frame.parent.cwd
            self.parent.original_frame.parent.reprocess = True
            self.parent.original_frame.parent.Close()
            if(self.parent.original_frame.parent.path !=''):
                os.chdir(self.parent.original_frame.parent.path)
            from SpinExplorer.SpinView import MyApp
            app = MyApp()
            if(self.parent.original_frame.parent.cwd !=''):
                app.path = path
                app.cwd = cwd

                
        # Check to see if the output file exists
        if(os.path.exists(nmrfile) == False):
            if(original_frame==True):
                if(app.cwd !=''):
                    os.chdir(app.cwd)
            elif(self.parent.cwd !=''):
                os.chdir(self.parent.cwd)
            message = 'The processing spectrum file ({}) file cannot be found in the current directory. Processing unsuccessful. Ensure that nmrPipe has been downloaded and added to the path.'.format(nmrfile)
            dlg = wx.MessageDialog(self, message, 'Warning', wx.OK | wx.ICON_WARNING)
            result = dlg.ShowModal()
            
                
            return
        else:
            if(original_frame==True):
                if(app.cwd !=''):
                    os.chdir(app.cwd)
            elif(self.parent.cwd !=''):
                os.chdir(self.parent.cwd)
            message = 'Processing successful. The processed spectrum file ({}) has been created in the current directory.'.format(nmrfile)
            dlg = wx.MessageDialog(self, message, 'Success', wx.OK | wx.ICON_INFORMATION)
            result = dlg.ShowModal()
            return
        
        

        

        
        



class OneDFrame(wx.Panel):

    def __init__(self, parent):
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 0.7*self.monitorWidth
        self.height = 0.75*self.monitorHeight
        self.parent = parent
        wx.Panel.__init__(self, parent, id=wx.ID_ANY, size=(self.width, self.height))
        if(darkdetect.isDark() == False or platform=='windows'):
            self.SetBackgroundColour('#edeeef')
        else:
            self.SetBackgroundColour('#282A36')
        # Create panel for processing dimension 1 of the data
        self.nmr_data = parent.nmr_data
        self.set_variables()
        self.create_canvas()
        self.create_menu_bar()

    def set_variables(self):

        # See if NMR processing file (nmrproc.com) can be found, if it can try to load the variables from it
        if(os.path.exists('processing_parameters.txt')):
            found_nmrproc_com = True
        else:
            found_nmrproc_com = False

        self.set_initial_solvent_suppression_variables()
        self.set_initial_linear_prediction_variables()
        self.set_initial_apodization_variables()
        self.set_initial_zero_filling_variables()
        self.set_initial_fourier_transform_variables()
        self.set_initial_phasing_variables()
        self.set_initial_extraction_variables()
        self.set_initial_baseline_correction_variables()
        
        self.parent.load_variables = False
        if(found_nmrproc_com == False):
            pass
        else:
            # Ask the user if they want to load the variables from the nmrproc.com file
            dlg = wx.MessageDialog(self, 'A file containing NMR processing parameters has been found (processing_parameters.txt). Do you want to load the variables from it?', 'Warning', wx.YES_NO | wx.ICON_WARNING)
            self.Raise()
            self.SetFocus()
            result = dlg.ShowModal()
            if(result == wx.ID_YES):
                try:  
                    self.parent.load_variables = True
                    self.load_variables_from_nmrproc_com_1D()
                except:
                    pass

            else:
                pass

    def load_variables_from_nmrproc_com_1D(self):
        # Open processing_parameters.txt file and load the variables from it
        file = open('processing_parameters.txt', 'r')
        lines = file.readlines()
        file.close()

        self.direct_solvent_suppression = False
        self.linear_prediction_checkbox_value = False
        self.apodization_checkbox_value = False
        self.zero_filling_checkbox_value = False
        self.fourier_transform_checkbox_value = False
        self.phase_correction_checkbox_value = False
        self.extraction_checkbox_value = False
        self.baseline_correction_checkbox_value = False

        include_line = False
        for line in lines:
            if('Dimension 1' in line):
                include_line = True
                continue
            if(include_line == False):
                continue
            if(include_line == True and 'Dimension 2' in line):
                include_line = False
                break
            if(include_line == True):
                line = line.split('\n')[0]
                if(line.split(':')[0]=='Solvent Suppression'):
                    if(line.split(': ')[1].strip()=='True' in line):
                        self.direct_solvent_suppression = True
                    else:
                        self.direct_solvent_suppression = False
                elif(line.split(':')[0]=='Filter Selection'):
                    self.solvent_suppression_filter_selection = int(line.split(': ')[1])
                elif(line.split(':')[0]=='Lowpass Shape Selection'):
                    self.solvent_suppression_lowpass_shape_selection = int(line.split(': ')[1])
                elif(line.split(':')[0]=='Linear Prediction'):
                    if(line.split(': ')[1].strip()=='True'):
                        self.linear_prediction_checkbox_value = True
                    else:
                        self.linear_prediction_checkbox_value = False
                elif(line.split(':')[0]=='Linear Prediction Options Selection'):
                    self.linear_prediction_options_selection = int(line.split(': ')[1])
                elif(line.split(':')[0]=='Linear Prediction Coefficients Selection'):
                    self.linear_prediction_coefficients_selection = int(line.split(': ')[1])
                elif(line.split(':')[0]=='Apodization'):
                    if(line.split(': ')[1].strip()=='True' in line):
                        self.apodization_checkbox_value = True
                    else:
                        self.apodization_checkbox_value = False
                elif(line.split(':')[0]=='Apodization Combobox Selection'):
                    self.apodization_combobox_selection = int(line.split(': ')[1])
                    self.apodization_combobox_selection_old = int(line.split(': ')[1])
                elif(line.split(':')[0]=='Exponential Line Broadening'):
                    self.exponential_line_broadening = float(line.split(': ')[1])
                elif(line.split(':')[0]=='Apodization First Point Scaling'):
                    self.apodization_first_point_scaling = float(line.split(': ')[1])
                elif(line.split(':')[0]=='G1'):
                    self.g1 = float(line.split(': ')[1])
                elif(line.split(':')[0]=='G2'):
                    self.g2 = float(line.split(': ')[1])
                elif(line.split(':')[0]=='G3'):
                    self.g3 = float(line.split(': ')[1])
                elif(line.split(':')[0]=='Offset'):
                    self.offset = float(line.split(': ')[1])
                elif(line.split(':')[0]=='End'):
                    self.end = float(line.split(': ')[1])
                elif(line.split(':')[0]=='Power'):
                    self.power = float(line.split(': ')[1])
                elif('A:' in line):
                    self.a = float(line.split(': ')[1])
                elif('B:' in line):
                    self.b = float(line.split(': ')[1])
                elif(line.split(':')[0]=='T1'):
                    self.t1 = float(line.split(': ')[1])
                elif(line.split(':')[0]=='T2'):
                    self.t2 = float(line.split(': ')[1])
                elif(line.split(':')[0]=='Loc'):
                    self.loc = float(line.split(': ')[1])
                elif(line.split(':')[0]=='Zero Filling'):
                    if('True' in line):
                        self.zero_filling_checkbox_value = True
                    else:
                        self.zero_filling_checkbox_value = False
                elif(line.split(':')[0]=='Zero Filling Combobox Selection'):
                    self.zero_filling_combobox_selection = int(line.split(': ')[1])
                elif(line.split(':')[0]=='Zero Filling Value Doubling Times'):
                    self.zero_filling_value_doubling_times = int(line.split(': ')[1])
                elif(line.split(':')[0]=='Zero Filling Value Zeros to Add'):
                    self.zero_filling_value_zeros_to_add = int(line.split(': ')[1])
                elif(line.split(':')[0]=='Zero Filling Value Final Data Size'):
                    self.zero_filling_value_final_data_size = int(line.split(': ')[1])
                elif(line.split(':')[0]=='Zero Filling Round Checkbox'):
                    if('True' in line):
                        self.zero_filling_round_checkbox_value = True
                    else:
                        self.zero_filling_round_checkbox_value = False
                elif(line.split(':')[0]=='Fourier Transform'):
                    if('True' in line):
                        self.fourier_transform_checkbox_value = True
                    else:
                        self.fourier_transform_checkbox_value = False
                elif(line.split(':')[0]=='Fourier Transform Method Selection'):
                    self.ft_method_selection = int(line.split(': ')[1])
                elif(line.split(':')[0]=='Phase Correction'):
                    if('True' in line):
                        self.phase_correction_checkbox_value = True
                    else:
                        self.phase_correction_checkbox_value = False
                elif(line.split(':')[0]=='Phase Correction P0'):
                    self.p0_total = float(line.split(': ')[1])
                elif(line.split(':')[0]=='Phase Correction P1'):
                    self.p1_total = float(line.split(': ')[1])
                elif(line.split(':')[0]=='Magnitude Mode'):
                    if('True' in line):
                        self.magnitude_mode_toggle = True
                    else:
                        self.magnitude_mode_toggle = False
                elif(line.split(':')[0]=='Extraction'):
                    if('True' in line):
                        self.extraction_checkbox_value = True
                    else:
                        self.extraction_checkbox_value = False
                elif(line.split(':')[0]=='Extraction PPM Start'):
                    self.extraction_ppm_start = float(line.split(': ')[1])
                elif(line.split(':')[0]=='Extraction PPM End'):
                    self.extraction_ppm_end = float(line.split(': ')[1])
                elif(line.split(':')[0]=='Baseline Correction'):
                    if('True' in line):
                        self.baseline_correction_checkbox_value = True
                    else:
                        self.baseline_correction_checkbox_value = False
                elif(line.split(':')[0]=='Baseline Correction Radio Box Selection'):
                    self.baseline_correction_radio_box_selection = int(line.split(': ')[1])
                elif(line.split(':')[0]=='Baseline Correction Nodes'):
                    self.baseline_correction_nodes = int(line.split(': ')[1])
                elif(line.split(':')[0]=='Baseline Correction Node List'):
                    self.baseline_correction_node_list = line.split(': ')[1]
                elif(line.split(':')[0]=='Baseline Correction Polynomial Order'):
                    self.baseline_correction_polynomial_order = int(line.split(': ')[1])





    def set_initial_solvent_suppression_variables(self):
        if(self.nmr_data.axis_labels[0] == 'H1' or self.nmr_data.axis_labels[0] == '1H' or self.nmr_data.axis_labels[0] == 'H'):
            self.direct_solvent_suppression = True
        else:
            self.direct_solvent_suppression = False

        self.include_direct_linear_prediction = False

        self.solvent_suppression_filter_selection = 0
        self.solvent_suppression_lowpass_shape_selection = 0
        self.solvent_suppression_filter_length = 16
        self.solvent_suppression_polynomial_order = 2
        self.solvent_suppression_spline_noise = 1.0
        self.solvent_suppression_spline_smoothfactor = 1.1
        
    
    def set_initial_linear_prediction_variables(self):
        self.linear_prediction_checkbox_value = False
        self.linear_prediction_options_selection = 0
        self.linear_prediction_coefficients_selection = 0


    def set_initial_apodization_variables(self):
        self.apodization_checkbox_value = True
        self.apodization_combobox_selection = 1
        self.apodization_combobox_selection_old = 1
        
        # Initial values for exponential apodization
        self.exponential_line_broadening = 0.5
        self.apodization_first_point_scaling = 0.5

        # Initial values for Lorentz to Gauss apodization
        self.g1 = 0.33
        self.g2 = 1
        self.g3 = 0.0

        # Initial values for Sinebell apodization
        self.offset = 0.5
        self.end = 0.98
        self.power = 1.0

        # Initial values for Gauss Broadening apodization
        self.a = 1.0
        self.b = 1.0

        # Initial values for Trapezoid apodization
        self.t1 = int((self.nmr_data.number_of_points[0]/2)/4)
        self.t2 = int((self.nmr_data.number_of_points[0]/2)/4)

        # Initial values for Triangle apodization
        self.loc = 0.5

    def set_initial_zero_filling_variables(self):
        self.zero_filling_checkbox_value = True
        self.zero_filling_combobox_selection = 0
        self.zero_filling_combobox_selection_old = 0
        self.zero_filling_value_doubling_times = 1
        self.zero_filling_value_zeros_to_add = 0
        self.zero_filling_value_final_data_size = self.nmr_data.number_of_points[0]*2
        self.zero_filling_round_checkbox_value = True

    def set_initial_fourier_transform_variables(self):
        self.fourier_transform_checkbox_value = True
        self.ft_method_selection = 0 # Initially use the 'auto' method of FT as default


    def set_initial_phasing_variables(self):
        self.phase_correction_checkbox_value = True
        self.p0_total = 0.0
        self.p1_total = 0.0
        self.magnitude_mode_toggle = False



    def set_initial_extraction_variables(self):
        self.extraction_checkbox_value = False
        self.extraction_ppm_start = 0.0
        self.extraction_ppm_end = 0.0




    def set_initial_baseline_correction_variables(self):
        self.baseline_correction_radio_box_selection = 0
        self.node_width = '2'
        self.node_list = '0,5,95,100'
        self.polynomial_order = '4'



        

    def create_canvas(self):

        if(darkdetect.isDark() == True and platform != 'windows'):
            self.SetBackgroundColour((53, 53, 53, 255))
            matplotlib.rc('axes',edgecolor='white')
            matplotlib.rc('xtick',color='white')
            matplotlib.rc('ytick',color='white')
            matplotlib.rc('axes',labelcolor='white')
            matplotlib.rc('axes',facecolor='#282A36')
            matplotlib.rc('legend', labelcolor='white')


        else:
            self.SetBackgroundColour('White')
            matplotlib.rc('axes',edgecolor='black')
            matplotlib.rc('xtick',color='black')
            matplotlib.rc('ytick',color='black')
            matplotlib.rc('axes',labelcolor='black')
            matplotlib.rc('axes',facecolor='white')
            matplotlib.rc('legend', labelcolor='black')


        

        


    def create_menu_bar(self):
        # Create the main sizer
        self.main_sizer=wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_sizer)

        self.sizer_1 = wx.BoxSizer(wx.VERTICAL)
        self.sizer_1.AddSpacer(10)

    
        # Create all the sizers
        self.create_solvent_suppression_sizer(parent=self)
        self.create_linear_prediction_sizer(parent=self)
        self.create_apodization_sizer(parent=self)
        self.create_zero_filling_sizer(parent=self)
        self.create_fourier_transform_sizer(parent=self)
        self.create_phase_correction_sizer(parent=self)
        self.create_extraction_sizer(parent=self)
        self.create_baseline_correction_sizer(parent=self)
        
        self.main_sizer.Add(self.sizer_1, 0, wx.EXPAND)


        self.SetSizerAndFit(self.main_sizer)
        self.Layout()

        # Get the size of the main sizer and set the window size to 1.05 times the size of the main sizer
        self.width, self.height = self.main_sizer.GetSize()
        self.parent.parent.change_frame_size(int(self.width*1.05), int(self.height*1.25))




        


    def create_solvent_suppression_sizer(self, parent):
         # Create a box for solvent suppression options
        self.solvent_suppression_box = wx.StaticBox(parent, -1, 'Solvent Suppression')
        self.solvent_suppression_sizer = wx.StaticBoxSizer(self.solvent_suppression_box, wx.HORIZONTAL)
        self.solvent_suppression_checkbox = wx.CheckBox(parent, -1, 'Apply solvent suppression')
        self.solvent_suppression_checkbox.SetValue(self.direct_solvent_suppression)
        self.solvent_suppression_sizer.Add(self.solvent_suppression_checkbox, 0, wx.ALIGN_CENTER_VERTICAL)
        self.solvent_suppression_sizer.AddSpacer(10)
        self.solvent_suppression_extra_options = wx.Button(parent, -1, 'Advanced Options')
        self.solvent_suppression_sizer.Add(self.solvent_suppression_extra_options, 0, wx.ALIGN_CENTER_VERTICAL)
        self.solvent_suppression_sizer.AddSpacer(10)
        self.solvent_suppression_extra_options.Bind(wx.EVT_BUTTON, self.solvent_suppression_extra_options_click)
        # Have a button showing information on solvent suppression
        self.solvent_suppression_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))
        self.solvent_suppression_info.Bind(wx.EVT_BUTTON, self.on_solvent_suppression_info)
        self.solvent_suppression_sizer.Add(self.solvent_suppression_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer_1.Add(self.solvent_suppression_sizer)
        self.sizer_1.AddSpacer(10)


    def create_linear_prediction_sizer(self, parent):
        # Create a box for linear prediction options
        self.linear_prediction_box = wx.StaticBox(parent, -1, 'Linear Prediction')
        self.linear_prediction_sizer = wx.StaticBoxSizer(self.linear_prediction_box, wx.HORIZONTAL)
        self.linear_prediction_checkbox = wx.CheckBox(parent, -1, 'Apply linear prediction')
        self.linear_prediction_checkbox.SetValue(self.linear_prediction_checkbox_value)
        self.linear_prediction_checkbox.Bind(wx.EVT_CHECKBOX, self.on_linear_prediction_checkbox)
        self.linear_prediction_sizer.Add(self.linear_prediction_checkbox, 0, wx.ALIGN_CENTER_VERTICAL)
        self.linear_prediction_sizer.AddSpacer(10)
        # Have a combobox for linear prediction options
        self.linear_prediction_options_text = wx.StaticText(parent, -1, 'Add Predicted Points:')
        self.linear_prediction_sizer.Add(self.linear_prediction_options_text, 0, wx.ALIGN_CENTER_VERTICAL)
        self.linear_prediction_sizer.AddSpacer(5)
        self.linear_prediction_options = ['After FID', 'Before FID']
        self.linear_prediction_combobox = wx.ComboBox(parent, -1, choices=self.linear_prediction_options, style=wx.CB_READONLY)
        self.linear_prediction_combobox.Bind(wx.EVT_COMBOBOX, self.on_linear_prediction_combobox_options)
        self.linear_prediction_combobox.SetSelection(self.linear_prediction_options_selection)
        self.linear_prediction_sizer.Add(self.linear_prediction_combobox, 0, wx.ALIGN_CENTER_VERTICAL)
        self.linear_prediction_sizer.AddSpacer(10)
        # Have a combobox of predicted coefficient options
        self.linear_prediction_coefficients_text = wx.StaticText(parent, -1, 'Predicted Coefficients:')
        self.linear_prediction_sizer.Add(self.linear_prediction_coefficients_text, 0, wx.ALIGN_CENTER_VERTICAL)
        self.linear_prediction_sizer.AddSpacer(5)
        self.linear_prediction_coefficients_options = ['Forward', 'Backward', 'Both']
        self.linear_prediction_coefficients_combobox = wx.ComboBox(parent, -1, choices=self.linear_prediction_coefficients_options, style=wx.CB_READONLY)
        self.linear_prediction_coefficients_combobox.Bind(wx.EVT_COMBOBOX, self.on_linear_prediction_coefficients_combobox)
        self.linear_prediction_coefficients_combobox.SetSelection(self.linear_prediction_coefficients_selection)
        self.linear_prediction_sizer.Add(self.linear_prediction_coefficients_combobox, 0, wx.ALIGN_CENTER_VERTICAL)
        self.linear_prediction_sizer.AddSpacer(10)

        # Have a button showing information on linear prediction
        self.linear_prediction_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))
        self.linear_prediction_info.Bind(wx.EVT_BUTTON, self.on_linear_prediction_info)
        self.linear_prediction_sizer.Add(self.linear_prediction_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer_1.Add(self.linear_prediction_sizer)
        self.sizer_1.AddSpacer(10)   

    def on_linear_prediction_checkbox(self, event):
        if(self.linear_prediction_checkbox.GetValue() == True):
            self.linear_prediction_checkbox_value = True
        else:
            self.linear_prediction_checkbox_value = False 

    def on_linear_prediction_combobox_options(self, event):
        self.linear_prediction_options_selection = self.linear_prediction_combobox.GetSelection()
    
    def on_linear_prediction_coefficients_combobox(self, event):
        self.linear_prediction_coefficients_selection = self.linear_prediction_coefficients_combobox.GetSelection()
    

    def create_apodization_sizer(self, parent):
        # Create a box for apodization options
        self.apodization_box = wx.StaticBox(parent, -1, 'Apodization')
        self.apodization_sizer = wx.StaticBoxSizer(self.apodization_box, wx.HORIZONTAL)
        self.apodization_checkbox = wx.CheckBox(parent, -1, 'Apply apodization')
        self.apodization_checkbox.SetValue(self.apodization_checkbox_value)
        self.apodization_checkbox.Bind(wx.EVT_CHECKBOX, self.on_apodization_checkbox)
        self.apodization_sizer.Add(self.apodization_checkbox, 0, wx.ALIGN_CENTER_VERTICAL)
        self.apodization_sizer.AddSpacer(10)
        # Have a combobox for apodization options
        self.apodization_options = ['None', 'Exponential', 'Lorentz to Gauss', 'Sinebell', 'Gauss Broadening', 'Trapazoid', 'Triangle']
        self.apodization_combobox = wx.ComboBox(parent, -1, choices=self.apodization_options, style=wx.CB_READONLY)
        self.apodization_combobox.SetSelection(self.apodization_combobox_selection)
        self.apodization_combobox.Bind(wx.EVT_COMBOBOX, self.on_apodization_combobox)
        self.apodization_sizer.Add(self.apodization_combobox, 0, wx.ALIGN_CENTER_VERTICAL)
        self.apodization_sizer.AddSpacer(10)
        if(self.apodization_combobox_selection == 1):
            # Have a textcontrol for the line broadening
            self.apodization_line_broadening_label = wx.StaticText(parent, -1, 'Line Broadening (Hz):')
            self.apodization_sizer.Add(self.apodization_line_broadening_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_line_broadening_textcontrol = wx.TextCtrl(parent, -1, str(self.exponential_line_broadening), size=(30, 20))
            self.apodization_line_broadening_textcontrol.Bind(wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol)
            self.apodization_sizer.Add(self.apodization_line_broadening_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(parent, -1, 'First Point Scaling:')
            self.apodization_sizer.Add(self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_first_point_textcontrol = wx.TextCtrl(parent, -1, str(self.apodization_first_point_scaling), size=(30, 20))
            self.apodization_sizer.Add(self.apodization_first_point_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer.AddSpacer(10)
        elif(self.apodization_combobox_selection == 2):
            # Have a textcontrol for the g1 value
            self.apodization_g1_label = wx.StaticText(parent, -1, 'Inverse Lorentzian (Hz):')
            self.apodization_sizer.Add(self.apodization_g1_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_g1_textcontrol = wx.TextCtrl(parent, -1, str(self.g1), size=(40, 20))
            self.apodization_g1_textcontrol.Bind(wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol)
            self.apodization_sizer.Add(self.apodization_g1_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the g2 value
            self.apodization_g2_label = wx.StaticText(parent, -1, 'Gaussian Broadening (Hz):')
            self.apodization_sizer.Add(self.apodization_g2_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_g2_textcontrol = wx.TextCtrl(parent, -1, str(self.g2), size=(40, 20))
            self.apodization_g2_textcontrol.Bind(wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol)
            self.apodization_sizer.Add(self.apodization_g2_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the g3 value
            self.apodization_g3_label = wx.StaticText(parent, -1, 'Gaussian Shift:')
            self.apodization_sizer.Add(self.apodization_g3_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_g3_textcontrol = wx.TextCtrl(parent, -1, str(self.g3), size=(40, 20))
            self.apodization_g3_textcontrol.Bind(wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol)
            self.apodization_sizer.Add(self.apodization_g3_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(parent, -1, 'First Point Scaling:')
            self.apodization_sizer.Add(self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_first_point_textcontrol = wx.TextCtrl(self, -1, str(self.apodization_first_point_scaling), size=(30, 20))
            self.apodization_sizer.Add(self.apodization_first_point_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer.AddSpacer(10)
        elif(self.apodization_combobox_selection == 3):
            # Have a textcontrol for the offset value
            self.apodization_offset_label = wx.StaticText(parent, -1, 'Offset (\u03c0):')
            self.apodization_sizer.Add(self.apodization_offset_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_offset_textcontrol = wx.TextCtrl(parent, -1, str(self.offset), size=(40, 20))
            self.apodization_offset_textcontrol.Bind(wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol)
            self.apodization_sizer.Add(self.apodization_offset_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the end value
            self.apodization_end_label = wx.StaticText(parent, -1, 'End (\u03c0):')
            self.apodization_sizer.Add(self.apodization_end_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_end_textcontrol = wx.TextCtrl(parent, -1, str(self.end), size=(40, 20))
            self.apodization_end_textcontrol.Bind(wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol)
            self.apodization_sizer.Add(self.apodization_end_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the power value
            self.apodization_power_label = wx.StaticText(parent, -1, 'Power:')
            self.apodization_sizer.Add(self.apodization_power_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_power_textcontrol = wx.TextCtrl(parent, -1, str(self.power), size=(30, 20))
            self.apodization_power_textcontrol.Bind(wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol)
            self.apodization_sizer.Add(self.apodization_power_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(parent, -1, 'First Point Scaling:')
            self.apodization_sizer.Add(self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_first_point_textcontrol = wx.TextCtrl(parent, -1, str(self.apodization_first_point_scaling), size=(30, 20))
            self.apodization_sizer.Add(self.apodization_first_point_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer.AddSpacer(10)
        elif(self.apodization_combobox_selection == 4):
            # Have a textcontrol for the a value
            self.apodization_a_label = wx.StaticText(parent, -1, 'Line Broadening (Hz):')
            self.apodization_sizer.Add(self.apodization_a_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_a_textcontrol = wx.TextCtrl(parent, -1, str(self.a), size=(40, 20))
            self.apodization_a_textcontrol.Bind(wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol)
            self.apodization_sizer.Add(self.apodization_a_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the b value
            self.apodization_b_label = wx.StaticText(parent, -1, 'Gaussian Broadening (Hz):')
            self.apodization_sizer.Add(self.apodization_b_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_b_textcontrol = wx.TextCtrl(parent, -1, str(self.b), size=(40, 20))
            self.apodization_b_textcontrol.Bind(wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol)
            self.apodization_sizer.Add(self.apodization_b_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(parent, -1, 'First Point Scaling:')
            self.apodization_sizer.Add(self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_first_point_textcontrol = wx.TextCtrl(self, -1, str(self.apodization_first_point_scaling), size=(30, 20))
            self.apodization_sizer.Add(self.apodization_first_point_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer.AddSpacer(10)
        elif(self.apodization_combobox_selection == 5):
            # Have a textcontrol for the t1 value
            self.apodization_t1_label = wx.StaticText(parent, -1, 'Ramp up points:')
            self.apodization_sizer.Add(self.apodization_t1_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_t1_textcontrol = wx.TextCtrl(parent, -1, str(self.t1), size=(50, 20))
            self.apodization_t1_textcontrol.Bind(wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol)
            self.apodization_sizer.Add(self.apodization_t1_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the t2 value
            self.apodization_t2_label = wx.StaticText(parent, -1, 'Ramp down points:')
            self.apodization_sizer.Add(self.apodization_t2_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_t2_textcontrol = wx.TextCtrl(parent, -1, str(self.t2), size=(50, 20))
            self.apodization_t2_textcontrol.Bind(wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol)
            self.apodization_sizer.Add(self.apodization_t2_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(parent, -1, 'First Point Scaling:')
            self.apodization_sizer.Add(self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_first_point_textcontrol = wx.TextCtrl(parent, -1, str(self.apodization_first_point_scaling), size=(30, 20))
            self.apodization_sizer.Add(self.apodization_first_point_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer.AddSpacer(10)
        elif(self.apodization_combobox_selection == 6):
            # Have a textcontrol for the loc value
            self.apodization_loc_label = wx.StaticText(parent, -1, 'Location of maximum:')
            self.apodization_sizer.Add(self.apodization_loc_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_loc_textcontrol = wx.TextCtrl(parent, -1, str(self.loc), size=(40, 20))
            self.apodization_loc_textcontrol.Bind(wx.EVT_CHAR_HOOK, self.on_apodization_textcontrol)
            self.apodization_sizer.Add(self.apodization_loc_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(parent, -1, 'First Point Scaling:')
            self.apodization_sizer.Add(self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_first_point_textcontrol = wx.TextCtrl(parent, -1, '0.5', size=(30, 20))
            self.apodization_sizer.Add(self.apodization_first_point_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer.AddSpacer(10)

            
        # Have a button for information on currently selected apodization containing unicode i in a circle
        self.apodization_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))
        self.apodization_info.Bind(wx.EVT_BUTTON, self.on_apodization_info)
        self.apodization_sizer.Add(self.apodization_info, 0, wx.ALIGN_CENTER_VERTICAL)

        # Have a mini plot of the apodization function along with the FID first slice
        self.plot_window_function()




        self.sizer_1.Add(self.apodization_sizer)
        self.sizer_1.AddSpacer(10)

    def on_apodization_checkbox(self, event):
        if(self.apodization_checkbox.GetValue() == True):
            self.apodization_checkbox_value = True
        else:
            self.apodization_checkbox_value = False


    def on_apodization_textcontrol(self, event):
        # If the user presses enter, update the plot
        keycode = event.GetKeyCode()
        if(keycode == wx.WXK_RETURN):
            self.update_window_function_plot()
        event.Skip()


    

        
    def create_zero_filling_sizer(self, parent):
     # Create a box for zero filling options
        self.zero_filling_box = wx.StaticBox(parent, -1, 'Zero Filling')
        self.zero_filling_sizer = wx.StaticBoxSizer(self.zero_filling_box, wx.HORIZONTAL)
        self.zero_filling_checkbox = wx.CheckBox(parent, -1, 'Apply zero filling')
        self.zero_filling_checkbox.SetValue(self.zero_filling_checkbox_value)
        self.zero_filling_checkbox.Bind(wx.EVT_CHECKBOX, self.on_zero_filling_checkbox)
        self.zero_filling_sizer.Add(self.zero_filling_checkbox, 0, wx.ALIGN_CENTER_VERTICAL)
        self.zero_filling_sizer.AddSpacer(10)
        # Have a combobox for zero filling options
        self.zf_options_label = wx.StaticText(parent, -1, 'Options:')
        self.zero_filling_sizer.Add(self.zf_options_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.zero_filling_sizer.AddSpacer(5)
        self.zero_filling_options = ['Doubling spectrum size', 'Adding Zeros', 'Final data size']
        self.zero_filling_combobox = wx.ComboBox(parent, -1, choices=self.zero_filling_options, style=wx.CB_READONLY)
        self.zero_filling_combobox.Bind(wx.EVT_COMBOBOX, self.on_zero_filling_combobox)
        self.zero_filling_combobox.SetSelection(self.zero_filling_combobox_selection)
        self.zero_filling_sizer.Add(self.zero_filling_combobox, 0, wx.ALIGN_CENTER_VERTICAL)
        self.zero_filling_sizer.AddSpacer(10)
        if(self.zero_filling_combobox_selection == 0):
            # Have a textcontrol for the doubling number/number of zeros/final data size
            self.zf_value_label = wx.StaticText(parent, -1, 'Doubling number:')
            self.zero_filling_sizer.Add(self.zf_value_label, 0, wx.ALIGN_CENTER_VERTICAL)

            self.zero_filling_textcontrol = wx.TextCtrl(parent, -1, str(self.zero_filling_value_doubling_times), size=(40, 20))
            self.zero_filling_textcontrol.Bind(wx.EVT_TEXT, self.on_zero_filling_textcontrol_doubling_times)
            self.zero_filling_sizer.AddSpacer(5)
            self.zero_filling_sizer.Add(self.zero_filling_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.zero_filling_sizer.AddSpacer(20)
        elif(self.zero_filling_combobox_selection == 1):
            # Have a textcontrol for the doubling number/number of zeros/final data size
            self.zf_value_label = wx.StaticText(parent, -1, 'Number of zeros to add:')
            self.zero_filling_sizer.Add(self.zf_value_label, 0, wx.ALIGN_CENTER_VERTICAL)

            self.zero_filling_textcontrol = wx.TextCtrl(parent, -1, str(self.zero_filling_value_zeros_to_add), size=(40, 20))
            self.zero_filling_textcontrol.Bind(wx.EVT_TEXT, self.on_zero_filling_textcontrol_zeros_to_add)
            self.zero_filling_sizer.AddSpacer(5)
            self.zero_filling_sizer.Add(self.zero_filling_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.zero_filling_sizer.AddSpacer(20)
        elif(self.zero_filling_combobox_selection == 2):
            # Have a textcontrol for the doubling number/number of zeros/final data size
            self.zf_value_label = wx.StaticText(parent, -1, 'Final data size:')
            self.zero_filling_sizer.Add(self.zf_value_label, 0, wx.ALIGN_CENTER_VERTICAL)

            self.zero_filling_textcontrol = wx.TextCtrl(parent, -1, str(self.zero_filling_value_final_data_size), size=(40, 20))
            self.zero_filling_textcontrol.Bind(wx.EVT_TEXT, self.on_zero_filling_textcontrol_final_size)
            self.zero_filling_sizer.AddSpacer(5)
            self.zero_filling_sizer.Add(self.zero_filling_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
            self.zero_filling_sizer.AddSpacer(20)

        # Have a checkbox for rounding to the nearest power of 2
        self.zero_filling_round_checkbox = wx.CheckBox(parent, -1, 'Round to nearest power of 2')
        self.zero_filling_round_checkbox.SetValue(True)
        self.zero_filling_sizer.Add(self.zero_filling_round_checkbox, 0, wx.ALIGN_CENTER_VERTICAL)
        self.zero_filling_sizer.AddSpacer(10)
        

        # Have a button showing information on zero filling
        self.zero_filling_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))
        self.zero_filling_info.Bind(wx.EVT_BUTTON, self.on_zero_fill_info)
        self.zero_filling_sizer.Add(self.zero_filling_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.zero_filling_sizer.AddSpacer(10)

        self.sizer_1.Add(self.zero_filling_sizer)
        self.sizer_1.AddSpacer(10)


    def on_zero_filling_checkbox(self, event):
        if(self.zero_filling_checkbox.GetValue() == True):
            self.zero_filling_checkbox_value = True
        else:
            self.zero_filling_checkbox_value = False

    def on_zero_filling_textcontrol_doubling_times(self, event):
        try:
            self.zero_filling_value_doubling_times = int(self.zero_filling_textcontrol.GetValue())
        except:
            if(self.zero_filling_textcontrol.GetValue() == ''):
                self.zero_filling_value_doubling_times = 0
                self.zero_filling_textcontrol.SetValue('')
            else:
                # Give an popout error message saying that the value must be an integer
                self.zero_filling_value_doubling_times = 1
                self.zero_filling_textcontrol.SetValue(str(self.zero_filling_value_doubling_times))

                message = 'The value for the zero filling doubling number must be an integer. Resetting value to 1.'
                title = 'Invalid value'
                style = wx.OK | wx.ICON_ERROR
                wx.MessageBox(message, title, style)


    def on_zero_filling_textcontrol_zeros_to_add(self, event):
        try:
            self.zero_filling_value_zeros_to_add = int(self.zero_filling_textcontrol.GetValue())
        except:
            if(self.zero_filling_textcontrol.GetValue() == ''):
                self.zero_filling_value_zeros_to_add = 0
                self.zero_filling_textcontrol.SetValue('')
            else:
                # Give an popout error message saying that the value must be an integer
                self.zero_filling_value_zeros_to_add = 0
                self.zero_filling_textcontrol.SetValue(str(self.zero_filling_value_zeros_to_add))

                message = 'The value for the zero filling (zeros to add) number must be an integer. Resetting value to 0.'
                title = 'Invalid value'
                style = wx.OK | wx.ICON_ERROR
                wx.MessageBox(message, title, style)


    def on_zero_filling_textcontrol_final_size(self, event):
        try:
            self.zero_filling_value_final_data_size = int(self.zero_filling_textcontrol.GetValue())
        except:
            if(self.zero_filling_textcontrol.GetValue() == ''):
                self.zero_filling_value_final_data_size = self.nmr_data.number_of_points[0]*2
                self.zero_filling_textcontrol.SetValue('')
            else:
                # Give an popout error message saying that the value must be an integer
                self.zero_filling_value_final_data_size = self.nmr_data.number_of_points[0]*2
                self.zero_filling_textcontrol.SetValue(str(self.zero_filling_value_final_data_size))

                message = 'The value for the zero filling final data size must be an integer. Resetting value to {}.'.format(self.nmr_data.number_of_points[0]*2)
                title = 'Invalid value'
                style = wx.OK | wx.ICON_ERROR
                wx.MessageBox(message, title, style)
        


                                                   




    def create_fourier_transform_sizer(self, parent):
        # Create a box for fourier transform options
        self.fourier_transform_box = wx.StaticBox(parent, -1, 'Fourier Transform')
        self.fourier_transform_sizer = wx.StaticBoxSizer(self.fourier_transform_box, wx.HORIZONTAL)
        self.fourier_transform_checkbox = wx.CheckBox(parent, -1, 'Apply fourier transform')
        self.fourier_transform_checkbox.SetValue(True)
        self.fourier_transform_sizer.Add(self.fourier_transform_checkbox, 0, wx.ALIGN_CENTER_VERTICAL)
        self.fourier_transform_sizer.AddSpacer(10)
        # Have a button for advanced options for fourier transform
        self.fourier_transform_advanced_options = wx.Button(parent, -1, 'Advanced Options')
        self.fourier_transform_advanced_options.Bind(wx.EVT_BUTTON, self.on_fourier_transform_advanced_options)
        self.fourier_transform_sizer.Add(self.fourier_transform_advanced_options, 0, wx.ALIGN_CENTER_VERTICAL)
        self.fourier_transform_sizer.AddSpacer(10)

        # Have a button showing information on fourier transform
        self.fourier_transform_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))
        self.fourier_transform_info.Bind(wx.EVT_BUTTON, self.on_fourier_transform_info)
        self.fourier_transform_sizer.Add(self.fourier_transform_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer_1.Add(self.fourier_transform_sizer)
        self.sizer_1.AddSpacer(10)


    def create_phase_correction_sizer(self, parent):
        # Create a box for phase correction options
        self.phase_correction_box = wx.StaticBox(parent, -1, 'Phase Correction')
        self.phase_correction_sizer = wx.StaticBoxSizer(self.phase_correction_box, wx.HORIZONTAL)
        self.phase_correction_checkbox = wx.CheckBox(parent, -1, 'Apply phase correction')
        self.phase_correction_checkbox.SetValue(self.phase_correction_checkbox_value)
        self.phase_correction_checkbox.Bind(wx.EVT_CHECKBOX, self.on_phase_correction_checkbox)
        self.phase_correction_sizer.Add(self.phase_correction_checkbox, 0, wx.ALIGN_CENTER_VERTICAL)
        self.phase_correction_sizer.AddSpacer(10)
        # Have a textcontrol for p0 and p1 values
        self.phase_correction_p0_label = wx.StaticText(parent, -1, 'Zero order correction (p0):')
        self.phase_correction_sizer.Add(self.phase_correction_p0_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.phase_correction_p0_textcontrol = wx.TextCtrl(parent, -1, str(self.p0_total), size=(50, 20))
        self.phase_correction_p0_textcontrol.Bind(wx.EVT_TEXT, self.on_phase_correction_textcontrol)
        self.phase_correction_sizer.Add(self.phase_correction_p0_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
        self.phase_correction_sizer.AddSpacer(10)
        self.phase_correction_p1_label = wx.StaticText(parent, -1, 'First order correction (p1):')
        self.phase_correction_sizer.Add(self.phase_correction_p1_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.phase_correction_p1_textcontrol = wx.TextCtrl(parent, -1, str(self.p1_total), size=(50, 20))
        self.phase_correction_p1_textcontrol.Bind(wx.EVT_TEXT, self.on_phase_correction_textcontrol)
        self.phase_correction_sizer.Add(self.phase_correction_p1_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
        self.phase_correction_sizer.AddSpacer(10)

        # A button to toggle magnitude mode
        self.magnitude_mode_label = wx.StaticText(parent,-1,'Magnitude Mode:')
        self.phase_correction_sizer.Add(self.magnitude_mode_label,0,wx.ALIGN_CENTER_VERTICAL)
        self.magnitude_mode_checkbox = wx.CheckBox(parent,-1)
        self.magnitude_mode_checkbox.SetValue(self.magnitude_mode_toggle)
        self.phase_correction_sizer.Add(self.magnitude_mode_checkbox,0,wx.ALIGN_CENTER_VERTICAL)
        self.phase_correction_sizer.AddSpacer(10)

        # Have a button for automatic phase correction
        self.phase_correction_auto_button = wx.Button(parent, -1, 'Interactive Phase Correction')
        self.phase_correction_auto_button.Bind(wx.EVT_BUTTON, self.on_phase_correction_interactive)
        self.phase_correction_sizer.Add(self.phase_correction_auto_button, 0, wx.ALIGN_CENTER_VERTICAL)
        self.phase_correction_sizer.AddSpacer(10)

        # Have a button showing information on phase correction
        self.phase_correction_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))
        self.phase_correction_info.Bind(wx.EVT_BUTTON, self.on_phase_correction_info)
        self.phase_correction_sizer.Add(self.phase_correction_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer_1.Add(self.phase_correction_sizer)
        self.sizer_1.AddSpacer(10)


    def on_phase_correction_checkbox(self, event):
        if(self.phase_correction_checkbox.GetValue() == True):
            self.phase_correction_checkbox_value = True
        else:
            self.phase_correction_checkbox_value = False


    def on_phase_correction_textcontrol(self, event):
        self.p0_total = self.phase_correction_p0_textcontrol.GetValue()
        self.p1_total = self.phase_correction_p1_textcontrol.GetValue()
        try:
            if(np.abs(float(self.p1_total))>45):
                self.apodization_first_point_scaling = 1.0
                self.apodization_first_point_textcontrol.SetValue('1.0')
            else:
                self.apodization_first_point_scaling = 0.5
                self.apodization_first_point_textcontrol.SetValue('0.5')
        except:
            pass


    def create_extraction_sizer(self, parent):
        # A box for extraction of data between two ppm values
        self.extraction_box = wx.StaticBox(parent, -1, 'Extraction')
        self.extraction_sizer = wx.StaticBoxSizer(self.extraction_box, wx.HORIZONTAL)
        self.extraction_checkbox = wx.CheckBox(parent, -1, 'Include data extraction')
        self.extraction_checkbox.Bind(wx.EVT_CHECKBOX, self.on_extraction_checkbox)
        self.extraction_checkbox.SetValue(self.extraction_checkbox_value)
        self.extraction_sizer.Add(self.extraction_checkbox, 0, wx.ALIGN_CENTER_VERTICAL)
        self.extraction_sizer.AddSpacer(10)
        # Have a textcontrol for the ppm start value
        self.extraction_ppm_start_label = wx.StaticText(parent, -1, 'Start chemical shift (ppm):')
        self.extraction_sizer.Add(self.extraction_ppm_start_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.extraction_ppm_start_textcontrol = wx.TextCtrl(parent, -1, str(self.extraction_ppm_start), size=(40, 20))
        self.extraction_ppm_start_textcontrol.Bind(wx.EVT_TEXT, self.on_extraction_textcontrol)
        self.extraction_sizer.Add(self.extraction_ppm_start_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
        self.extraction_sizer.AddSpacer(10)
        # Have a textcontrol for the ppm end value
        self.extraction_ppm_end_label = wx.StaticText(parent, -1, 'End chemical shift (ppm):')
        self.extraction_sizer.Add(self.extraction_ppm_end_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.extraction_ppm_end_textcontrol = wx.TextCtrl(parent, -1, str(self.extraction_ppm_end), size=(40, 20))
        self.extraction_ppm_end_textcontrol.Bind(wx.EVT_TEXT, self.on_extraction_textcontrol)
        self.extraction_sizer.Add(self.extraction_ppm_end_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
        self.extraction_sizer.AddSpacer(10)
        # Have a button showing information on extraction
        self.extraction_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))
        self.extraction_info.Bind(wx.EVT_BUTTON, self.on_extraction_info)
        self.extraction_sizer.Add(self.extraction_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer_1.Add(self.extraction_sizer)
        self.sizer_1.AddSpacer(10)


    def on_extraction_checkbox(self, event):
        self.extraction_checkbox_value = self.extraction_checkbox.GetValue()

    def on_extraction_textcontrol(self, event):
        self.extraction_ppm_start = self.extraction_ppm_start_textcontrol.GetValue()
        self.extraction_ppm_end = self.extraction_ppm_end_textcontrol.GetValue()




    def create_baseline_correction_sizer(self, parent):
        # Create a box for baseline correction options (linear/polynomial)
        self.baseline_correction_box = wx.StaticBox(parent, -1, 'Baseline Correction')
        self.baseline_correction_sizer = wx.StaticBoxSizer(self.baseline_correction_box, wx.HORIZONTAL)
        self.baseline_correction_checkbox = wx.CheckBox(parent, -1, 'Apply baseline correction')
        self.baseline_correction_checkbox.SetValue(False)
        self.baseline_correction_sizer.Add(self.baseline_correction_checkbox, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_sizer.AddSpacer(10)
        # Have a radio box for linear or polynomial baseline correction
        self.baseline_correction_radio_box = wx.RadioBox(parent, -1, 'Baseline Correction Method', choices=['Linear', 'Polynomial'])
        # Bind the radio box to a function that will update the baseline correction options
        self.baseline_correction_radio_box.Bind(wx.EVT_RADIOBOX, self.on_baseline_correction_radio_box)
        self.baseline_correction_radio_box.SetSelection(self.baseline_correction_radio_box_selection)
        self.baseline_correction_sizer.Add(self.baseline_correction_radio_box, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_sizer.AddSpacer(10)
        
        # If linear baseline correction is selected, have a textcontrol for the node values to use
        self.baseline_correction_nodes_label = wx.StaticText(parent, -1, 'Node width (pts):')
        self.baseline_correction_sizer.Add(self.baseline_correction_nodes_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_nodes_textcontrol = wx.TextCtrl(parent, -1, str(self.node_width), size=(30, 20))
        self.baseline_correction_nodes_textcontrol.Bind(wx.EVT_TEXT, self.on_baseline_correction_textcontrol)
        self.baseline_correction_sizer.Add(self.baseline_correction_nodes_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_sizer.AddSpacer(10)
        # Have a textcontrol for the node list (percentages)
        self.baseline_correction_node_list_label = wx.StaticText(parent, -1, 'Node list (%):')
        self.baseline_correction_sizer.Add(self.baseline_correction_node_list_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_node_list_textcontrol = wx.TextCtrl(parent, -1, str(self.node_list), size=(100, 20))
        self.baseline_correction_node_list_textcontrol.Bind(wx.EVT_TEXT, self.on_baseline_correction_textcontrol)
        self.baseline_correction_sizer.Add(self.baseline_correction_node_list_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_sizer.AddSpacer(10)
        # If polynomial baseline correction is selected, have a textcontrol for the polynomial order

        self.baseline_correction_polynomial_order_label = wx.StaticText(parent, -1, 'Polynomial order:')
        self.baseline_correction_sizer.Add(self.baseline_correction_polynomial_order_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_polynomial_order_textcontrol = wx.TextCtrl(parent, -1, str(self.polynomial_order), size=(30, 20))
        self.baseline_correction_polynomial_order_textcontrol.Bind(wx.EVT_TEXT, self.on_baseline_correction_textcontrol)
        self.baseline_correction_sizer.Add(self.baseline_correction_polynomial_order_textcontrol, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_sizer.AddSpacer(10)

        if(self.baseline_correction_radio_box_selection == 0):
            self.baseline_correction_polynomial_order_label.Hide()
            self.baseline_correction_polynomial_order_textcontrol.Hide()


        # Have a button showing information on baseline correction
        self.baseline_correction_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))


        self.baseline_correction_info.Bind(wx.EVT_BUTTON, self.on_baseline_correction_info)
        self.baseline_correction_sizer.Add(self.baseline_correction_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer_1.Add(self.baseline_correction_sizer)
        self.sizer_1.AddSpacer(10)


    def on_baseline_correction_textcontrol(self, event):
        self.node_width = self.baseline_correction_nodes_textcontrol.GetValue()
        self.node_list = self.baseline_correction_node_list_textcontrol.GetValue()
        self.polynomial_order = self.baseline_correction_polynomial_order_textcontrol.GetValue()




    def find_nmr_data_for_phasing(self):
        # Check to see if nmrpipe fid file exists
        # Check to see what the path of the original frame is
        if(self.parent.parent.original_frame != None):
            if(self.parent.parent.original_frame.parent.path != ''):
                path = self.parent.parent.original_frame.parent.path
                os.chdir(path)
        if(self.parent.parent.file_parser == True):
            os.chdir(self.parent.parent.path)

        if(os.path.exists('./test.fid') == False):
            # Error and exit
            self.error_message = wx.MessageDialog(self, 'No NMRPipe FID file found in the current directory. Please ensure that the file is named test.fid and is in the current directory.', 'Error', wx.OK | wx.ICON_ERROR)
            self.error_message.ShowModal()
            return
        
        if(platform!='windows'):
            # Create a temporary nmrproc.com file for phasing
            self.nmrproc_com_file = open('nmrproc_phasing.com', 'w')

            # Write the nmrproc.com file for phasing
            self.nmrproc_com_file.write('#!/bin/csh\n\n')
            self.nmrproc_com_file.write('nmrPipe -in test.fid \\\n')
            self.nmrproc_com_file.write('| nmrPipe -fn EM -lb 0.3 -c 0.5 \\\n')
            self.nmrproc_com_file.write('| nmrPipe -fn ZF -zf 1 -auto \\\n')
            self.nmrproc_com_file.write('| nmrPipe -fn FT -auto \\\n')
            self.nmrproc_com_file.write('| nmrPipe -fn PS -p0 0.0 -p1 0.0 \\\n')
            self.nmrproc_com_file.write(' -ov -out test_phasing.ft\n')

            self.nmrproc_com_file.close()

            # Run the nmrproc.com file
            os.system('csh nmrproc_phasing.com')

            # Read in the nmr data
            self.nmr_d, self.nmr_spectrum = ng.pipe.read('test_phasing.ft')

            # Remove the temporary nmrproc.com file
            os.system('rm nmrproc_phasing.com')

            # Remove the temporary test_phasing.ft file
            os.system('rm test_phasing.ft')

        else:
            dic,data = ng.pipe.read('test.fid')
            # Perform at fourier transform in the direct dimension
            dic,data = ng.pipe_proc.ft(dic,data,auto=True)
            self.nmr_d,self.nmr_spectrum = dic, data

        # Check to see what the path of the original frame is
        if(self.parent.parent.original_frame != None):
            if(self.parent.parent.original_frame.parent.cwd != ''):
                cwd = self.parent.parent.original_frame.parent.cwd
                os.chdir(cwd)

        if(self.parent.parent.file_parser == True):
            os.chdir(self.parent.parent.cwd)



        # Get the ppm scale
        self.uc = ng.pipe.make_uc(self.nmr_d, self.nmr_spectrum, dim=-1)

        # Get the ppm scale
        self.ppm_scale = self.uc.ppm_scale()






    def on_phase_correction_interactive(self, event):

        self.find_nmr_data_for_phasing()
        
        # Make a new window with the interactive phase correction 
        self.interactive_phase_correction_window = InteractivePhasingFrame(self, self.nmr_spectrum, self.ppm_scale, self.nmr_d)







        
        
        









    def on_fourier_transform_advanced_options(self, event):
        # Create a frame with a set of advanced options for the fourier transform implementation
        self.fourier_transform_advanced_options_window = wx.Frame(self, -1, 'Fourier Transform Advanced Options', size=(700, 300))
        if(darkdetect.isDark() == True and platform != 'windows'):
            self.fourier_transform_advanced_options_window.SetBackgroundColour((53, 53, 53, 255))
        else:
            self.fourier_transform_advanced_options_window.SetBackgroundColour('White')
        
        self.fourier_transform_advanced_options_window_sizer = wx.BoxSizer(wx.VERTICAL)
        self.fourier_transform_advanced_options_window.SetSizer(self.fourier_transform_advanced_options_window_sizer)

        # Create a sizer for the fourier transform advanced options
        self.ft_label = wx.StaticBox(self.fourier_transform_advanced_options_window, -1, 'Fourier Transform Method:')
        self.fourier_transform_advanced_options_sizer = wx.StaticBoxSizer(self.ft_label,wx.VERTICAL)

        # Have a radiobox for auto, real, inverse, sign alternation
        self.fourier_transform_advanced_options_sizer.AddSpacer(10)
        self.fourier_transform_auto_real_inverse_sign_alternation_radio_box = wx.RadioBox(self.fourier_transform_advanced_options_window, -1, choices=['Auto', 'Real', 'Inverse', 'Sign Alternation', 'Negative'], style=wx.RA_SPECIFY_COLS)
        self.fourier_transform_auto_real_inverse_sign_alternation_radio_box.SetSelection(self.ft_method_selection)
        self.fourier_transform_advanced_options_sizer.Add(self.fourier_transform_auto_real_inverse_sign_alternation_radio_box, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.fourier_transform_advanced_options_sizer.AddSpacer(10)


        self.ft_method_text = 'Auto: The auto method will automatically select the best method for the fourier transform of the FID. \n\nReal: The Fourier Transform will be applied to the real part of the FID only. \n\nInverse: The inverse Fourier Transform will be applied to the FID. \n\nSign Alternation: The sign alternation method will be applied to the FID. \n\n'

        self.ft_method_info = wx.StaticText(self.fourier_transform_advanced_options_window, -1, self.ft_method_text)
        self.fourier_transform_advanced_options_sizer.Add(self.ft_method_info, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.fourier_transform_advanced_options_sizer.AddSpacer(10)

        # Have a save and close button
        self.fourier_transform_advanced_options_save_button = wx.Button(self.fourier_transform_advanced_options_window, -1, 'Save and Close')
        self.fourier_transform_advanced_options_save_button.Bind(wx.EVT_BUTTON, self.on_fourier_transform_advanced_options_save)
        self.fourier_transform_advanced_options_sizer.Add(self.fourier_transform_advanced_options_save_button, 0, wx.ALIGN_CENTER_HORIZONTAL)



        self.fourier_transform_advanced_options_window_sizer.Add(self.fourier_transform_advanced_options_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)

        self.fourier_transform_advanced_options_window.Show()

    def on_fourier_transform_advanced_options_save(self, event):
        # Save the current selection and close the window
        self.ft_method_selection = self.fourier_transform_auto_real_inverse_sign_alternation_radio_box.GetSelection()
        self.fourier_transform_advanced_options_window.Close()


    def on_fourier_transform_info(self, event):
        ft_text = 'The fourier transform applies a complex fourier transform to the FID to convert it to a frequency domain spectrum.\n Further information can be found using the link below.'

        # Create a popup window with the information
        self.fourier_transform_info_window = wx.Frame(self, -1, 'Fourier Transform Information', size=(450, 150))
        if(darkdetect.isDark() == True and platform != 'windows'):
            self.fourier_transform_info_window.SetBackgroundColour((53, 53, 53, 255))
            colour = "RED"
        else:
            self.fourier_transform_info_window.SetBackgroundColour('White')
            colour = "BLUE"
        
        self.fourier_transform_info_window_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.fourier_transform_info_window.SetSizer(self.fourier_transform_info_window_sizer)

        self.fourier_transform_info_window_sizer.AddSpacer(10)

        # Create a sizer for the fourier transform information
        self.fourier_transform_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.fourier_transform_info_sizer.AddSpacer(10)
        self.fourier_transform_info_sizer.Add(wx.StaticText(self.fourier_transform_info_window, -1, ft_text, size=(400,50)), 0, wx.ALIGN_CENTER)
        self.fourier_transform_info_sizer.AddSpacer(10)

        # Have a hyperlink to the fourier transform information
        self.fourier_transform_info_hyperlink = hl.HyperLinkCtrl(self.fourier_transform_info_window, -1, 'NMRPipe Help Page for Fourier Transform', URL='http://www.nmrscience.com/ref/nmrpipe/ft.html')
        self.fourier_transform_info_hyperlink.SetColours(colour, colour, colour)
        self.fourier_transform_info_hyperlink.SetUnderlines(False, False, False)
        self.fourier_transform_info_hyperlink.SetBold(False)
        self.fourier_transform_info_hyperlink.UpdateLink()
        self.fourier_transform_info_sizer.Add(self.fourier_transform_info_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.fourier_transform_info_sizer.AddSpacer(10)

        self.fourier_transform_info_window_sizer.Add(self.fourier_transform_info_sizer, 0, wx.ALIGN_CENTER)

        self.fourier_transform_info_window.Show()


    def on_phase_correction_info(self, event):
        phase_correction_text = 'Phase correction is a method to correct for phase errors in the FID. Zero order phase correction (p0) is used to correct a phase offset that is applied equally across the spectrum. However, a first order phase correction (p1) is used to correct the phasing in a spectrum where peaks in different locations of the spectrum require a different phasing value. \n Further information can be found using the link below.'

        # Create a popup window with the information
        self.phase_correction_info_window = wx.Frame(self, -1, 'Phase Correction Information', size=(450, 200))
        if(darkdetect.isDark() == True and platform != 'windows'):
            self.phase_correction_info_window.SetBackgroundColour((53, 53, 53, 255))
            colour = "RED"
        else:
            self.phase_correction_info_window.SetBackgroundColour('White')
            colour = "BLUE"
        
        self.phase_correction_info_window_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.phase_correction_info_window.SetSizer(self.phase_correction_info_window_sizer)

        self.phase_correction_info_window_sizer.AddSpacer(10)

        # Create a sizer for the phase correction information
        self.phase_correction_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.phase_correction_info_sizer.AddSpacer(10)
        self.phase_correction_info_sizer.Add(wx.StaticText(self.phase_correction_info_window, -1, phase_correction_text, size=(400,100)), 0, wx.ALIGN_CENTER)
        self.phase_correction_info_sizer.AddSpacer(10)

        # Have a hyperlink to the phase correction information
        self.phase_correction_info_hyperlink = hl.HyperLinkCtrl(self.phase_correction_info_window, -1, 'NMRPipe Help Page for Phase Correction', URL='http://www.nmrscience.com/ref/nmrpipe/ps.html')
        self.phase_correction_info_hyperlink.SetColours(colour, colour, colour)
        self.phase_correction_info_hyperlink.SetUnderlines(False, False, False)
        self.phase_correction_info_hyperlink.SetBold(False)
        self.phase_correction_info_hyperlink.UpdateLink()
        self.phase_correction_info_sizer.Add(self.phase_correction_info_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.phase_correction_info_sizer.AddSpacer(10)

        self.phase_correction_info_window_sizer.Add(self.phase_correction_info_sizer, 0, wx.ALIGN_CENTER)

        self.phase_correction_info_window.Show()
        



    def on_extraction_info(self, event):
        extraction_text = 'Extraction of data between two chemical shift values can be used to extract a region of interest from the spectrum. This can be useful for removing solvent signals or other unwanted peaks. \n Further information can be found using the link below.'

        # Create a popup window with the information
        self.extraction_info_window = wx.Frame(self, -1, 'Extraction Information', size=(450, 200))
        if(darkdetect.isDark() == True and platform != 'windows'):
            self.extraction_info_window.SetBackgroundColour((53, 53, 53, 255))
            colour = "RED"
        else:
            self.extraction_info_window.SetBackgroundColour('White')
            colour = "BLUE"
        
        self.extraction_info_window_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.extraction_info_window.SetSizer(self.extraction_info_window_sizer)

        self.extraction_info_window_sizer.AddSpacer(10)

        # Create a sizer for the extraction information
        self.extraction_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.extraction_info_sizer.AddSpacer(10)
        self.extraction_info_sizer.Add(wx.StaticText(self.extraction_info_window, -1, extraction_text, size=(400,100)), 0, wx.ALIGN_CENTER)
        self.extraction_info_sizer.AddSpacer(10)

        # Have a hyperlink to the extraction information
        self.extraction_info_hyperlink = hl.HyperLinkCtrl(self.extraction_info_window, -1, 'NMRPipe Help Page for Extraction', URL='http://www.nmrscience.com/ref/nmrpipe/ext.html')
        self.extraction_info_hyperlink.SetColours(colour, colour, colour)
        self.extraction_info_hyperlink.SetUnderlines(False, False, False)
        self.extraction_info_hyperlink.SetBold(False)
        self.extraction_info_hyperlink.UpdateLink()

        self.extraction_info_sizer.Add(self.extraction_info_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.extraction_info_sizer.AddSpacer(10)

        self.extraction_info_window_sizer.Add(self.extraction_info_sizer, 0, wx.ALIGN_CENTER)

        self.extraction_info_window.Show()




    def on_baseline_correction_info(self, event):
        # Include information on linear and polynomial baseline correction
        linear_information = 'Linear baseline correction: \n Linear baseline correction is a method to correct for a linear baseline issue in the spectrum. The linear baseline is removed by fitting a straight line to the spectrum and subtracting it from the spectrum. \n\n'
        polynomial_information = 'Polynomial baseline correction: \n Polynomial baseline correction is a method to correct for a polynomial baseline issue in the spectrum. The polynomial baseline is removed by fitting a polynomial to the spectrum and subtracting it from the spectrum. \n\n'
        extra_information = 'The node list is a list of percentages that are used to define the nodes (points which are expected to have 0 intensity) for the baseline correction. The node width is the number of points used to define the nodes. Further advanced options can be added to the processing file nmrproc.com file manually \n\n'
        # Create a popup window with the information
        self.baseline_correction_info_window = wx.Frame(self, -1, 'Baseline Correction Information', size=(450, 450))
        if(darkdetect.isDark() == True and platform != 'windows'):
            self.baseline_correction_info_window.SetBackgroundColour((53, 53, 53, 255))
            colour = "RED"
        else:
            self.baseline_correction_info_window.SetBackgroundColour('White')
            colour = "BLUE"
        
        self.baseline_correction_info_window_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.baseline_correction_info_window.SetSizer(self.baseline_correction_info_window_sizer)

        self.baseline_correction_info_window_sizer.AddSpacer(10)

        # Create a sizer for the baseline correction information
        self.baseline_correction_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.baseline_correction_info_sizer.AddSpacer(10)
        self.baseline_correction_info_sizer.Add(wx.StaticText(self.baseline_correction_info_window, -1, linear_information + polynomial_information + extra_information, size=(400,300)), 0, wx.ALIGN_CENTER)
        self.baseline_correction_info_sizer.AddSpacer(10)

        # Have a hyperlink to the linear baseline correction information
        self.baseline_correction_info_hyperlink = hl.HyperLinkCtrl(self.baseline_correction_info_window, -1, 'NMRPipe Help Page for Linear Baseline Correction', URL='http://www.nmrscience.com/ref/nmrpipe/base.html')
        self.baseline_correction_info_hyperlink.SetColours(colour, colour, colour)
        self.baseline_correction_info_hyperlink.SetUnderlines(False, False, False)
        self.baseline_correction_info_hyperlink.SetBold(False)
        self.baseline_correction_info_hyperlink.UpdateLink()

        # Have a hyperlink to the polynomial baseline correction information
        self.baseline_correction_info_hyperlink_2 = hl.HyperLinkCtrl(self.baseline_correction_info_window, -1, 'NMRPipe Help Page for Polynomial Baseline Correction', URL='http://www.nmrscience.com/ref/nmrpipe/poly.html')
        self.baseline_correction_info_hyperlink_2.SetColours(colour, colour, colour)
        self.baseline_correction_info_hyperlink_2.SetUnderlines(False, False, False)
        self.baseline_correction_info_hyperlink_2.SetBold(False)
        self.baseline_correction_info_hyperlink_2.UpdateLink()

        self.baseline_correction_info_sizer.Add(self.baseline_correction_info_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.baseline_correction_info_sizer.AddSpacer(10)
        self.baseline_correction_info_sizer.Add(self.baseline_correction_info_hyperlink_2, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.baseline_correction_info_sizer.AddSpacer(10)

        self.baseline_correction_info_window_sizer.Add(self.baseline_correction_info_sizer, 0, wx.ALIGN_CENTER)

        self.baseline_correction_info_window.Show()




    def on_baseline_correction_radio_box(self, event):
        # If the user selects linear or polynomial baseline correction, update the options
        self.baseline_correction_radio_box_selection = self.baseline_correction_radio_box.GetSelection()

        if(self.baseline_correction_radio_box_selection == 0):
            # Remove the polynomial order textcontrol
            self.baseline_correction_sizer.Hide(self.baseline_correction_polynomial_order_label)
            self.baseline_correction_sizer.Hide(self.baseline_correction_polynomial_order_textcontrol)
            self.baseline_correction_sizer.Layout()
        elif(self.baseline_correction_radio_box_selection == 1):
            # Add the polynomial order textcontrol
            self.baseline_correction_sizer.Show(self.baseline_correction_polynomial_order_label)
            self.baseline_correction_sizer.Show(self.baseline_correction_polynomial_order_textcontrol)
            self.baseline_correction_sizer.Layout()


        

        



    def solvent_suppression_extra_options_click(self, event):
        # Give a popup window with options for solvent suppression
        self.solvent_suppression_extra_options_window = wx.Frame(self, -1, 'Solvent Suppression Advanced Options', size=(400, 300))
        if(darkdetect.isDark() == True and platform != 'windows'):
            self.solvent_suppression_extra_options_window.SetBackgroundColour((53, 53, 53, 255))
        else:
            self.solvent_suppression_extra_options_window.SetBackgroundColour('White')
        
        self.solvent_suppression_extra_options_window_sizer = wx.BoxSizer(wx.VERTICAL)
        self.solvent_suppression_extra_options_window.SetSizer(self.solvent_suppression_extra_options_window_sizer)

        # Create a sizer for the solvent suppression options
        self.solvent_suppression_extra_options_sizer = wx.BoxSizer(wx.VERTICAL)
        self.solvent_suppression_extra_options_window_sizer.Add(self.solvent_suppression_extra_options_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.solvent_suppression_extra_options_window_sizer.AddSpacer(10)

        # Have a radio button for low-pass filter, spline, or polynomial
        self.radio_box = wx.RadioBox(self.solvent_suppression_extra_options_window, -1, 'Solvent Suppression Method', choices=['Low-pass filter', 'Spline', 'Polynomial'], style=wx.RA_SPECIFY_ROWS)
        self.radio_box.SetSelection(self.solvent_suppression_filter_selection)
        self.solvent_suppression_extra_options_sizer.Add(self.radio_box, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.solvent_suppression_extra_options_sizer.AddSpacer(10)


        # Have a radiobox for the low-pass filter shape
        self.lowpass_shape_radio_box = wx.RadioBox(self.solvent_suppression_extra_options_window, -1, 'Low-pass filter shape', choices=['Boxcar','Sine', 'Sine Squared'], style=wx.RA_SPECIFY_ROWS)
        self.lowpass_shape_radio_box.SetSelection(self.solvent_suppression_lowpass_shape_selection)
        self.solvent_suppression_extra_options_sizer.Add(self.lowpass_shape_radio_box, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.solvent_suppression_extra_options_sizer.AddSpacer(10)



        self.solvent_suppression_extra_options_window.SetSizer(self.solvent_suppression_extra_options_window_sizer)
        self.solvent_suppression_extra_options_window.Show()




    def on_solvent_suppression_info(self, event):
        # Include information on how the low-pass, splines and polynomial solvent suppression works

        general_information = 'Aqueous NMR spectra often suffer from substantial solvent signals reducing spectral quality. \n If NMR spectra have been run with the solvent set at the carrier frequency there will low frequency oscillations in the FID\n (in the rotating frame) corresponding to the solvent signal. These oscillations can be removed by applying a low-pass filter to the FID. \n\n'

        lowpass_filter = 'Low-pass filter: \n A low-pass filter is a filter that passes signals with a frequency lower than a certain cutoff frequency and attenuates signals with frequencies higher than the cutoff frequency. \n\n'

        lowpass_filter_shape = 'Low-pass filter shape: \n The shape of the low-pass filter can be set to a boxcar, sine, or sine squared shape. \n\n'

        additional_advanced_options = 'Additional advanced options using splines or polynomial digital solvent suppression can also be selected. Further advanced options such as filter length (-fl) can be added manually to the nmrproc.com file to further refine the solvent suppression. \n\n'

        citation_information = 'Citation: \n Dominique Marion, Mitsuhiko Ikura, Ad Bax, Improved solvent suppression in one- and two-dimensional NMR spectra by convolution of time-domain data, Journal of Magnetic Resonance (1969), Volume 84, Issue 2, 1989, Pages 425-430, ISSN 0022-2364, https://doi.org/10.1016/0022-2364(89)90391-0. (https://www.sciencedirect.com/science/article/pii/0022236489903910)'

        # Create a popup window with the information
        self.solvent_suppression_info_window = wx.Frame(self, -1, 'Solvent Suppression Information', size=(450, 600))
        if(darkdetect.isDark() == True and platform != 'windows'):
            self.solvent_suppression_info_window.SetBackgroundColour((53, 53, 53, 255))
            colour = "RED"
        else:
            self.solvent_suppression_info_window.SetBackgroundColour('White')
            colour = "BLUE"

        self.solvent_suppression_info_window_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.solvent_suppression_info_window.SetSizer(self.solvent_suppression_info_window_sizer)

        self.solvent_suppression_info_window_sizer.AddSpacer(10)

        # Create a sizer for the solvent suppression information
        self.solvent_suppression_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.solvent_suppression_info_sizer.AddSpacer(10)
        
        self.solvent_suppression_info_window_sizer.AddSpacer(10)

        # Have a text control for the information
        self.solvent_suppression_info_textcontrol = wx.StaticText(self.solvent_suppression_info_window, label = general_information + lowpass_filter + lowpass_filter_shape + additional_advanced_options + citation_information, size=(400, 500))

        self.solvent_suppression_info_sizer.Add(self.solvent_suppression_info_textcontrol, 0, wx.ALIGN_CENTER_HORIZONTAL)

        # Add a hyperlink to the sizer for the NMRPipe SOL help page
        self.solvent_suppression_info_sizer.AddSpacer(10)
        
        self.sol_hyperlink = hl.HyperLinkCtrl(self.solvent_suppression_info_window, -1, "NMRPipe Solvent Suppression Help Page", URL="http://www.nmrscience.com/ref/nmrpipe/sol.html")
        
        self.sol_hyperlink.SetColours(colour, colour, colour)
        self.sol_hyperlink.SetUnderlines(False, False, False)
        self.sol_hyperlink.UpdateLink()

        self.solvent_suppression_info_sizer.Add(self.sol_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL)

        self.solvent_suppression_info_sizer.AddSpacer(10)

        
        self.solvent_suppression_info_window_sizer.Add(self.solvent_suppression_info_sizer, 0, wx.ALIGN_CENTER_VERTICAL)
        self.solvent_suppression_info_sizer.AddSpacer(10)
        self.solvent_suppression_info_window.SetSizer(self.solvent_suppression_info_window_sizer)
        self.solvent_suppression_info_window.Show()


    def on_zero_filling_combobox(self, event):
        self.zero_filling_combobox_selection = self.zero_filling_combobox.GetSelection()
        # self.zero_filling_sizer.Clear(delete_windows=True)

        self.zero_filling_sizer.Clear()
        self.zero_filling_sizer.Detach(self.zero_filling_checkbox)
        self.zero_filling_checkbox.Destroy()
        self.zero_filling_sizer.Detach(self.zf_options_label)
        self.zf_options_label.Destroy()
        self.zero_filling_sizer.Detach(self.zero_filling_info)
        self.zero_filling_info.Destroy()
        self.zero_filling_sizer.Detach(self.zf_value_label)
        self.zf_value_label.Destroy()
        self.zero_filling_sizer.Detach(self.zero_filling_round_checkbox)
        self.zero_filling_round_checkbox.Destroy()
        self.zero_filling_sizer.Detach(self.zero_filling_textcontrol)
        self.zero_filling_textcontrol.Destroy()

        self.zero_filling_sizer.Detach(self.zero_filling_combobox)
        self.zero_filling_combobox.Hide()



        if(self.apodization_combobox_selection_old == 1):
            # Remove the previous textcontrols
            
            self.apodization_sizer.Detach(self.apodization_line_broadening_label)
            self.apodization_line_broadening_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_line_broadening_textcontrol)
            self.apodization_line_broadening_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_textcontrol)
            self.apodization_first_point_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)

        elif(self.apodization_combobox_selection_old == 2):
            self.apodization_sizer.Detach(self.apodization_g1_label)
            self.apodization_g1_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_g1_textcontrol)
            self.apodization_g1_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_g2_label)
            self.apodization_g2_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_g2_textcontrol)
            self.apodization_g2_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_g3_label)
            self.apodization_g3_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_g3_textcontrol)
            self.apodization_g3_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_textcontrol)
            self.apodization_first_point_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)

        elif(self.apodization_combobox_selection_old == 3):
            self.apodization_sizer.Detach(self.apodization_offset_label)
            self.apodization_offset_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_offset_textcontrol)
            self.apodization_offset_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_end_label)
            self.apodization_end_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_end_textcontrol)
            self.apodization_end_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_power_label)
            self.apodization_power_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_power_textcontrol)
            self.apodization_power_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_textcontrol)
            self.apodization_first_point_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)
        elif(self.apodization_combobox_selection_old == 4):
            self.apodization_sizer.Detach(self.apodization_a_label)
            self.apodization_a_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_a_textcontrol)
            self.apodization_a_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_b_label)
            self.apodization_b_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_b_textcontrol)
            self.apodization_b_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_textcontrol)
            self.apodization_first_point_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)
        elif(self.apodization_combobox_selection_old == 5):
            self.apodization_sizer.Detach(self.apodization_t1_label)
            self.apodization_t1_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_t1_textcontrol)
            self.apodization_t1_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_t2_label)
            self.apodization_t2_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_t2_textcontrol)
            self.apodization_t2_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_textcontrol)
            self.apodization_first_point_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)
        elif(self.apodization_combobox_selection_old == 6):
            self.apodization_sizer.Detach(self.apodization_loc_label)
            self.apodization_loc_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_loc_textcontrol)
            self.apodization_loc_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_textcontrol)
            self.apodization_first_point_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)
        elif(self.apodization_combobox_selection_old == 0):
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)
        


        self.apodization_sizer.Detach(self.apodization_checkbox)
        self.apodization_checkbox.Destroy()
        self.apodization_sizer.Detach(self.apodization_combobox)
        self.apodization_combobox.Hide()
       


        

        # Delete the current apodization sizer and then create a new one
        
        self.apodization_sizer.Detach(self.apodization_plot_sizer)
        self.apodization_plot_ax.clear()
        self.apodization_plot.clear()
        self.apodization_plot_sizer.Clear(True)

        self.sizer_1.Remove(self.apodization_sizer)
        # self.apodization_sizer.Clear(delete_windows=True)


        # Remove the linear prediction sizers
        self.linear_prediction_sizer.Clear(delete_windows=True)
        # self.sizer_1.Remove(self.linear_prediction_sizer)
        

        # Remove the solvent suppression sizers
        self.solvent_suppression_sizer.Clear(delete_windows=True)
        # self.sizer_1.Remove(self.solvent_suppression_sizer)



        
        self.sizer_1.Clear(delete_windows=True)


        self.create_menu_bar()
        self.Refresh()
        self.Update()
        self.Layout()

        self.zero_filling_combobox_selection = self.zero_filling_combobox_selection_old




    def on_apodization_combobox(self, event):
        self.apodization_combobox_selection = self.apodization_combobox.GetSelection()

        # Destroy the combobox and textcontrols for the previous apodization function
        # self.apodization_sizer.Detach(self.apodization_combobox)
        # self.apodization_combobox.Destroy()

        # Remove the zf sizer
        self.zero_filling_sizer.Clear(delete_windows=True)

        if(self.apodization_combobox_selection_old == 1):
            # Remove the previous textcontrols
            
            self.apodization_sizer.Detach(self.apodization_line_broadening_label)
            self.apodization_line_broadening_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_line_broadening_textcontrol)
            self.apodization_line_broadening_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_textcontrol)
            self.apodization_first_point_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)

        elif(self.apodization_combobox_selection_old == 2):
            self.apodization_sizer.Detach(self.apodization_g1_label)
            self.apodization_g1_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_g1_textcontrol)
            self.apodization_g1_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_g2_label)
            self.apodization_g2_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_g2_textcontrol)
            self.apodization_g2_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_g3_label)
            self.apodization_g3_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_g3_textcontrol)
            self.apodization_g3_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_textcontrol)
            self.apodization_first_point_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)

        elif(self.apodization_combobox_selection_old == 3):
            self.apodization_sizer.Detach(self.apodization_offset_label)
            self.apodization_offset_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_offset_textcontrol)
            self.apodization_offset_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_end_label)
            self.apodization_end_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_end_textcontrol)
            self.apodization_end_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_power_label)
            self.apodization_power_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_power_textcontrol)
            self.apodization_power_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_textcontrol)
            self.apodization_first_point_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)
        elif(self.apodization_combobox_selection_old == 4):
            self.apodization_sizer.Detach(self.apodization_a_label)
            self.apodization_a_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_a_textcontrol)
            self.apodization_a_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_b_label)
            self.apodization_b_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_b_textcontrol)
            self.apodization_b_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_textcontrol)
            self.apodization_first_point_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)
        elif(self.apodization_combobox_selection_old == 5):
            self.apodization_sizer.Detach(self.apodization_t1_label)
            self.apodization_t1_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_t1_textcontrol)
            self.apodization_t1_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_t2_label)
            self.apodization_t2_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_t2_textcontrol)
            self.apodization_t2_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_textcontrol)
            self.apodization_first_point_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)
        elif(self.apodization_combobox_selection_old == 6):
            self.apodization_sizer.Detach(self.apodization_loc_label)
            self.apodization_loc_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_loc_textcontrol)
            self.apodization_loc_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer.Detach(self.apodization_first_point_textcontrol)
            self.apodization_first_point_textcontrol.Destroy()
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)
        elif(self.apodization_combobox_selection_old == 0):
            self.apodization_sizer.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer.Detach(self.apodization_plot_sizer)
            self.apodization_plot_sizer.Clear(True)
            self.apodization_plot_ax.clear()
            self.apodization_plot.clear()
            self.apodization_plot_sizer.Clear(True)
        


        self.apodization_sizer.Detach(self.apodization_checkbox)
        self.apodization_checkbox.Destroy()
        self.apodization_sizer.Detach(self.apodization_combobox)
        self.apodization_combobox.Hide()
       


        

        # Delete the current apodization sizer and then create a new one
        
        self.apodization_sizer.Detach(self.apodization_plot_sizer)
        self.apodization_plot_ax.clear()
        self.apodization_plot.clear()
        self.apodization_plot_sizer.Clear(True)

        self.sizer_1.Remove(self.apodization_sizer)
        # self.apodization_sizer.Clear(delete_windows=True)

 


        



        # Remove the linear prediction sizers
        self.linear_prediction_sizer.Clear(delete_windows=True)
        # self.sizer_1.Remove(self.linear_prediction_sizer)
        

        # Remove the solvent suppression sizers
        self.solvent_suppression_sizer.Clear(delete_windows=True)
        # self.sizer_1.Remove(self.solvent_suppression_sizer)

        
        self.sizer_1.Clear(delete_windows=True)

           
        self.create_menu_bar()
        self.Refresh()
        self.Update()
        self.Layout()




        self.apodization_combobox_selection_old = self.apodization_combobox_selection
        


    def on_apodization_info(self, event):
        # Include information on how the apodization functions work and when to use each one

        
        # Make a pop out window with all of the apodization information
        # Create a new frame
        self.apodization_info_frame = wx.Frame(self, -1, 'Apodization Information', size=(500, 600))
        if(darkdetect.isDark() == True and platform != 'windows'):
            self.apodization_info_frame.SetBackgroundColour('#414141')
            colour = "RED"
        else:
            self.apodization_info_frame.SetBackgroundColour('White')
            colour = "BLUE"

        
        # Create a sizer to hold the box
        self.apodization_info_sizer_window = wx.BoxSizer(wx.VERTICAL)
        self.apodization_info_sizer_window.AddSpacer(10)

        general_information = 'Apodization (window functions) can be multiplied to the FID to increase signal to noise/increase resolution. Furthermore, window functions can be applied to data whose FID has been truncated (often in the indirect dimension) in order to reduce the presence of Sinc wiggles in the NMR spectrum. Below is an explanation of common apodization functions and when to use them.\n\n'

        first_point_information = 'First point scaling: In data where there is no first order (p1) phase correction required a value of 0.5 should be used, otherwise a value of 1.0 should be used.\n\n'
        
        exponential_information = 'Exponential: \n The exponential window function is used to apply an exponential decay to the FID. This suppresses the noise at the end of an FID enhancing signal to noise, but leads to a reduction in resolution. The line broadening term dictates the speed of the decay of the exponential function applied to the FID. \n '
        exponential_info_hyperlink = hl.HyperLinkCtrl(self.apodization_info_frame, -1, 'Further information on the nmrPipe exponential window function', URL='http://www.nmrscience.com/ref/nmrpipe/em.html')
        exponential_info_hyperlink.SetColours(colour, colour, colour)
        exponential_info_hyperlink.SetUnderlines(False, False, False)
        exponential_info_hyperlink.UpdateLink()


        lorentz_to_gauss_information = 'Lorentz to Gauss: \n The Lorentz to Gauss window function is used to apply a Lorentzian to Gaussian transformation to the FID. The inverse Lorentzian parameter can be tuned to remove all Lorentzian character from peaks in the NMR spectrum leaving them with a pure Gaussian peak shape. This can enhance signal resolution and is advantageous for peak picking routines which are often more accurate with Gaussian peak shapes. The Gaussian broadening value can be increased further to enhance the signal to noise at the cost of reduced resolution. As a rule of thumb, the Gaussian broadening parameter should be 3x larger than the inverse Lorentzian value. \n'
        lorentz_to_gauss_info_hyperlink = hl.HyperLinkCtrl(self.apodization_info_frame, -1, 'Further information on the nmrPipe Lorentz to Gauss window function', URL='http://www.nmrscience.com/ref/nmrpipe/gm.html')
        lorentz_to_gauss_info_hyperlink.SetColours(colour, colour, colour)
        lorentz_to_gauss_info_hyperlink.SetUnderlines(False, False, False)
        lorentz_to_gauss_info_hyperlink.UpdateLink()

        sinebell_information = 'Sinebell: \n The Sinebell window function is used to apply a sinebell function to the FID. The offset value adjusts the phase of the sinebell function, the end value adjusts the end of the sinebell function, and the power value adjusts the power of the sinebell function.'
        sinebell_info_hyperlink = hl.HyperLinkCtrl(self.apodization_info_frame, -1, 'Further information on the nmrPipe Sinebell window function', URL='http://www.nmrscience.com/ref/nmrpipe/sp.html')
        sinebell_info_hyperlink.SetColours(colour, colour, colour)
        sinebell_info_hyperlink.SetUnderlines(False, False, False)
        sinebell_info_hyperlink.UpdateLink()

        gauss_broadening_information = 'Gauss broadening: \n The Gauss broadening window function is used to apply a Decaying Exponential / Gaussian function to the FID. The line broadening value adjusts the decay of the exponential function, whilst the gaussian broadening term adjusts the decay of the gaussian function. \n'
        gauss_broadening_hyperlink = hl.HyperLinkCtrl(self.apodization_info_frame, -1, 'Further information on the nmrPipe Gaussian broadening window function', URL='http://www.nmrscience.com/ref/nmrpipe/gmb.html')
        gauss_broadening_hyperlink.SetColours(colour, colour, colour)
        gauss_broadening_hyperlink.SetUnderlines(False, False, False)
        gauss_broadening_hyperlink.UpdateLink()

        trapeziod_broadening_information = 'Trapezoid: A trapezoid function is applied to the FID. The number of ramp up and ramp down points can be adjusted to change the shape of the trapezoid\n'
        trapezoid_hyperlink = hl.HyperLinkCtrl(self.apodization_info_frame, -1, 'Further information on the nmrPipe trapezoid window function', URL='http://www.nmrscience.com/ref/nmrpipe/tm.html')
        trapezoid_hyperlink.SetColours(colour, colour, colour)
        trapezoid_hyperlink.SetUnderlines(False, False, False)
        trapezoid_hyperlink.UpdateLink()

        triangle_information = 'Triangle: A triangle function is applied to the FID. The location of the maximum of the triangle can be adjusted between 0 (first point) and 1 (last point).\n'
        triangle_hyperlink = hl.HyperLinkCtrl(self.apodization_info_frame, -1, 'Further information on the nmrPipe triangle window function', URL='http://www.nmrscience.com/ref/nmrpipe/tri.html')
        triangle_hyperlink.SetColours(colour, colour, colour)
        triangle_hyperlink.SetUnderlines(False, False, False)
        triangle_hyperlink.UpdateLink()


        # Create a sizer to hold the text
        self.apodization_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.apodization_info_sizer.AddSpacer(10)


        self.apodization_info_text1 = wx.StaticText(self.apodization_info_frame, -1, general_information, style=wx.ALIGN_CENTER_HORIZONTAL)
        self.apodization_info_text2 = wx.StaticText(self.apodization_info_frame, -1, first_point_information, style=wx.ALIGN_CENTER_HORIZONTAL)
        self.apodization_info_text3 = wx.StaticText(self.apodization_info_frame, -1, exponential_information, style=wx.ALIGN_CENTER_HORIZONTAL)
        self.apodization_info_sizer.Add(self.apodization_info_text1, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.apodization_info_sizer.AddSpacer(10)
        self.apodization_info_sizer.Add(self.apodization_info_text2, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.apodization_info_sizer.AddSpacer(10)
        self.apodization_info_sizer.Add(self.apodization_info_text3, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.apodization_info_sizer.Add(exponential_info_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.apodization_info_sizer.AddSpacer(10)
        self.apodization_info_sizer.Add(wx.StaticText(self.apodization_info_frame, -1, lorentz_to_gauss_information, style=wx.ALIGN_CENTER_HORIZONTAL), 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.apodization_info_sizer.Add(lorentz_to_gauss_info_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.apodization_info_sizer.AddSpacer(10)
        self.apodization_info_sizer.Add(wx.StaticText(self.apodization_info_frame, -1, sinebell_information, style=wx.ALIGN_CENTER_HORIZONTAL), 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.apodization_info_sizer.Add(sinebell_info_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.apodization_info_sizer.AddSpacer(10)
        self.apodization_info_sizer.Add(wx.StaticText(self.apodization_info_frame, -1, gauss_broadening_information, style=wx.ALIGN_CENTER_HORIZONTAL), 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.apodization_info_sizer.Add(gauss_broadening_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.apodization_info_sizer.AddSpacer(10)
        self.apodization_info_sizer.Add(wx.StaticText(self.apodization_info_frame, -1, trapeziod_broadening_information, style=wx.ALIGN_CENTER_HORIZONTAL), 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.apodization_info_sizer.Add(trapezoid_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.apodization_info_sizer.AddSpacer(10)
        self.apodization_info_sizer.Add(wx.StaticText(self.apodization_info_frame, -1, triangle_information, style=wx.ALIGN_CENTER_HORIZONTAL), 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.apodization_info_sizer.Add(triangle_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.apodization_info_sizer.AddSpacer(10)



        # Add the sizer to the window sizer
        self.apodization_info_sizer_window.Add(self.apodization_info_sizer, 0, wx.ALIGN_CENTER)
        self.apodization_info_sizer_window.AddSpacer(10)

        # Add the window sizer to the frame
        self.apodization_info_frame.SetSizer(self.apodization_info_sizer_window)

        # Show the frame
        self.apodization_info_frame.Show()


    def OnPressWindow(self, event):
        # Create a matplotlib popout of the plot
        fig = plt.figure(facecolor='white')
        ax = fig.add_subplot(111)
        ax.set_facecolor('white')
        ax.spines['bottom'].set_color('k')
        ax.spines['top'].set_color('k') 
        ax.spines['right'].set_color('k')
        ax.spines['left'].set_color('k')
        ax.tick_params(axis='x', colors='k')
        ax.tick_params(axis='y', colors='k')
        ax.yaxis.label.set_color('k')
        ax.xaxis.label.set_color('k')
        line1_x, line1_y = self.line1.get_data()
        line2_x, line2_y = self.line2.get_data()
        ax.plot(line1_x, line1_y, color='#1f77b4')
        ax.plot(line2_x, line2_y, color='k')
        ax.set_xlim(-(self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0]/20, (self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0]*21/20)
        ax.set_ylim(-1.5, 1.5)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Normalised Intensity (a.u.)')
        plt.show()






    def plot_window_function(self):
        self.apodization_plot_sizer = wx.BoxSizer(wx.VERTICAL)
        if(darkdetect.isDark() == True and platform != 'windows'):
            self.apodization_plot = Figure(figsize=(1,0.5),facecolor='#3b3b3b')
        else:
            self.apodization_plot = Figure(figsize=(1,0.5),facecolor='#e5e6e7')
        self.apodization_plot_ax = self.apodization_plot.add_subplot(111)
        self.apodization_plot_canvas = FigCanvas(self, -1, self.apodization_plot)
        self.apodization_plot_canvas.mpl_connect('button_press_event', self.OnPressWindow)
        
        if(darkdetect.isDark() == True and platform != 'windows'):
            self.apodization_plot_ax.set_facecolor('#3b3b3b')
        else:
            self.apodization_plot_ax.set_facecolor('#e5e6e7')

        self.apodization_plot_ax.set_xticks([])
        self.apodization_plot_ax.set_yticks([])

        # If the apodization function is None, make remove the axes of the plot
        if(self.apodization_combobox_selection == 0):
            self.apodization_plot_ax.spines['top'].set_visible(False)
            self.apodization_plot_ax.spines['right'].set_visible(False)
            self.apodization_plot_ax.spines['bottom'].set_visible(False)
            self.apodization_plot_ax.spines['left'].set_visible(False)

        self.plot_window_function_input()

    def plot_window_function_input(self):
        
        # Is the digital filter removed before or after processing
        self.before_processing = False
        try:
            with open('fid.com', 'r') as file:
                data = file.readlines()
                for line in data:
                    if('-AMX' in line.split()):
                        self.before_processing = False
                        break
                    elif('-DMX' in line.split()):
                        self.before_processing = True
                        break
        except:
            self.before_processing = True

        # Plot the NMR FID
        self.data = self.nmr_data.data

        if(self.nmr_data.dim == 1):
            self.data = self.data
        elif(self.nmr_data.dim == 2):
            self.data = self.data[0]
        else:
            self.data = self.data[0][0]

        data = self.data

        x=np.linspace(0, (self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0], int(self.nmr_data.number_of_points[0]/2))
        x_data = np.linspace(0, (len(data))/self.nmr_data.spectral_width[0], int(len(data)))

        # if(self.before_processing == False):
        #     x=np.linspace(0, (self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0], int(self.nmr_data.number_of_points[0]/2))
        # else:
        #     x=np.linspace(0, (len(data))/self.nmr_data.spectral_width[0], int(len(data)))
        self.line2, = self.apodization_plot_ax.plot(x_data,data/max(data),color='k')
        if(self.apodization_combobox_selection == 1):
            # Exponential window function
            self.line1, = self.apodization_plot_ax.plot(x, np.exp(-(np.pi*x*self.exponential_line_broadening)), color='#1f77b4')
            self.apodization_plot_ax.set_ylim(-1.5, 1.5)
            self.apodization_plot_ax.set_xlim(-(self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0]/20, (self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0]*21/20)
            
        elif(self.apodization_combobox_selection == 2):
            # Lorentz to Gauss window function
            e = np.pi * self.nmr_data.number_of_points[0]/self.nmr_data.spectral_width[0] * self.g1
            g = 0.6 * np.pi * self.g2 * (self.g3 * ((self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0] - 1) - x)
            func = np.exp(e - g * g)
            self.line1, = self.apodization_plot_ax.plot(x, func, color='#1f77b4')
            self.apodization_plot_ax.set_ylim(-1.5, 1.5)
            self.apodization_plot_ax.set_xlim(-(self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0]/20, (self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0]*21/20)
        elif(self.apodization_combobox_selection == 3):
            # Sinebell window function
            func = np.sin((np.pi*self.offset + np.pi*(self.end-self.offset)*x)/((((self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0]))))**self.power
            self.line1, = self.apodization_plot_ax.plot(x, func, color='#1f77b4')
            self.apodization_plot_ax.set_ylim(-1.5, 1.5)
            self.apodization_plot_ax.set_xlim(-(self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0]/20, (self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0]*21/20)
        elif(self.apodization_combobox_selection == 4):
            # Gauss broadening window function
            func = np.exp(-self.a*(x**2) - self.b*x)
            self.line1, = self.apodization_plot_ax.plot(x, func, color='#1f77b4')
            self.apodization_plot_ax.set_ylim(-1.5, 1.5)
            self.apodization_plot_ax.set_xlim(-(self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0]/20, (self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0]*21/20)
        elif(self.apodization_combobox_selection == 5):
            # Trapazoid window function
            func = np.concatenate((np.linspace(0, 1, int(self.t1)), np.ones(int(self.nmr_data.number_of_points[0]/2) - int(self.t1) - int(self.t2)),np.linspace(1, 0, int(self.t2))))
            self.line1, = self.apodization_plot_ax.plot(x, func, color='#1f77b4')
            self.apodization_plot_ax.set_ylim(-1.5, 1.5)
            self.apodization_plot_ax.set_xlim(-(self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0]/20, (self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0]*21/20)
        elif(self.apodization_combobox_selection == 6):
            # Triangle window function
            func = np.concatenate((np.linspace(0, 1, int(self.loc*(self.nmr_data.number_of_points[0]/2))), np.linspace(1, 0, int((1-self.loc)*(self.nmr_data.number_of_points[0]/2)))))
            self.line1, = self.apodization_plot_ax.plot(x, func, color='#1f77b4')
            
            self.apodization_plot_ax.set_ylim(-1.5, 1.5)
            self.apodization_plot_ax.set_xlim(-(self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0]/20, (self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0]*21/20)



            
        self.apodization_plot_ax.set_xlim(-(self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0]/20, (self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0]*21/20)
            

        self.apodization_plot_sizer.Add(self.apodization_plot_canvas, 0, wx.EXPAND)


        self.apodization_sizer.Add(self.apodization_plot_sizer, 0, wx.ALIGN_CENTER_VERTICAL)
        self.apodization_sizer.AddSpacer(10)


    def update_window_function_plot(self):
        data = self.data
        if(self.before_processing == False):
            x=np.linspace(0, (self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0], int(self.nmr_data.number_of_points[0]/2))
        else:
            x=np.linspace(0, (len(data))/self.nmr_data.spectral_width[0], int(len(data)))
        try:
            c = float(self.apodization_first_point_textcontrol.GetValue())
            self.apodization_first_point_scaling = c
        except:
            # Give a popout window saying that the values are not valid
            msg = wx.MessageDialog(self, 'The value entered for apodization first point scaling is not valid (use 0.5 or 1.0)', 'Error', wx.OK | wx.ICON_ERROR)
            msg.ShowModal()
            msg.Destroy()
            self.apodization_first_point_textcontrol.SetValue(str(self.apodization_first_point_scaling))
            return
        if(c != 0.5 and c != 1.0):
            msg = wx.MessageDialog(self, 'The value entered for apodization first point scaling is not valid (use 0.5 or 1.0)', 'Error', wx.OK | wx.ICON_ERROR)
            msg.ShowModal()
            msg.Destroy()
            self.apodization_first_point_textcontrol.SetValue(str(self.apodization_first_point_scaling))
            return
        self.apodization_first_point_scaling = c
        if(self.apodization_combobox_selection==1):
            try:
                em = float(self.apodization_line_broadening_textcontrol.GetValue())
            except:
                # Give a popout window saying that the values are not valid
                msg = wx.MessageDialog(self, 'The values entered are not valid', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_line_broadening_textcontrol.SetValue(str(self.exponential_line_broadening))
                return
            self.exponential_line_broadening = em
            
            self.line1.set_ydata(np.exp(-(np.pi*x*self.exponential_line_broadening)))
        elif(self.apodization_combobox_selection==2):
            try:
                g1 = float(self.apodization_g1_textcontrol.GetValue())
                g2 = float(self.apodization_g2_textcontrol.GetValue())
                g3 = float(self.apodization_g3_textcontrol.GetValue())
            except:
                # Give a popout window saying that the values are not valid
                msg = wx.MessageDialog(self, 'The values entered are not valid', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_g1_textcontrol.SetValue(str(self.g1))
                self.apodization_g2_textcontrol.SetValue(str(self.g2))
                self.apodization_g3_textcontrol.SetValue(str(self.g3))
                return
            # Check to see if g3 is between 0 and 1
            if(g3 < 0 or g3 > 1):
                msg = wx.MessageDialog(self, 'Gaussian shift must be between 0 and 1', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_g3_textcontrol.SetValue(str(self.g3))
                return
            self.g1 = g1
            self.g2 = g2
            self.g3 = g3
            e = np.pi * (self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0] * self.g1
            g = 0.6 * np.pi * self.g2 * (self.g3 * ((self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0] - 1) - x)
            func = np.exp(e - g * g)
            self.line1.set_ydata(func)

            self.apodization_plot_ax.set_xlim(0, (self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0])


        elif(self.apodization_combobox_selection==3):
            try:
                offset = float(self.apodization_offset_textcontrol.GetValue())
                end = float(self.apodization_end_textcontrol.GetValue())
                power = float(self.apodization_power_textcontrol.GetValue())
                power = int(power)
            except:
                msg = wx.MessageDialog(self, 'The values entered are not valid', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_offset_textcontrol.SetValue(str(self.offset))
                self.apodization_end_textcontrol.SetValue(str(self.end))
                self.apodization_power_textcontrol.SetValue(str(self.power))
                return
            # Check that offset and end are between 0 and 1
            if(offset < 0 or offset > 1):
                msg = wx.MessageDialog(self, 'Offset values must be between 0 and 1', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_offset_textcontrol.SetValue(str(self.offset))
                return
            if(end < 0 or end > 1):
                msg = wx.MessageDialog(self, 'End values must be between 0 and 1', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_end_textcontrol.SetValue(str(self.end))
                return
            # Check that power is greater than 0
            if(power < 0):
                msg = wx.MessageDialog(self, 'Power must be greater than 0', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_power_textcontrol.SetValue(str(self.power))
                return
            self.offset = offset
            self.end = end
            self.power = power
            func = np.sin((np.pi*self.offset + np.pi*(self.end-self.offset)*x)/((((self.nmr_data.number_of_points[0]/2)/self.nmr_data.spectral_width[0]))))**self.power
            self.line1.set_ydata(func)
        elif(self.apodization_combobox_selection==4):
            try:
                a = float(self.apodization_a_textcontrol.GetValue())
                b = float(self.apodization_b_textcontrol.GetValue())
            except:
                msg = wx.MessageDialog(self, 'The values entered are not valid', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_a_textcontrol.SetValue(str(self.a))
                self.apodization_b_textcontrol.SetValue(str(self.b))
                return
            self.a = a
            self.b = b
            func = np.exp(-self.a*(x**2) - self.b*x)
            self.line1.set_ydata(func)
        elif(self.apodization_combobox_selection==5):
            try:
                t1 = float(self.apodization_t1_textcontrol.GetValue())
                t2 = float(self.apodization_t2_textcontrol.GetValue())
            except:
                msg = wx.MessageDialog(self, 'The values entered are not valid', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_t1_textcontrol.SetValue(str(self.t1))
                self.apodization_t2_textcontrol.SetValue(str(self.t2))
                return
            # Ensure that t1 and t2 are greater than 0
            if(t1 < 0 or t2 < 0):
                msg = wx.MessageDialog(self, 'Ramp up and ramp down points must be greater than 0', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_t1_textcontrol.SetValue(str(self.t1))
                self.apodization_t2_textcontrol.SetValue(str(self.t2))
                return
            # Ensure that t1 + t2 is less than the number of points
            if(t1 + t2 > (self.nmr_data.number_of_points[0]/2)):
                message = 'Ramp up and ramp down points must be less than the number of points (' + str(self.nmr_data.number_of_points[0]) + ')'
                msg = wx.MessageDialog(self, message, 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_t1_textcontrol.SetValue(str(self.t1))
                self.apodization_t2_textcontrol.SetValue(str(self.t2))
                return
            self.t1 = t1
            self.t2 = t2
            func = np.concatenate((np.linspace(0, 1, int(self.t1)), np.ones(int(self.nmr_data.number_of_points[0]/2) - int(self.t1) - int(self.t2)),np.linspace(1, 0, int(self.t2))))
            self.line1.set_ydata(func)
        elif(self.apodization_combobox_selection==6):
            try:
                loc = float(self.apodization_loc_textcontrol.GetValue())
            except:
                msg = wx.MessageDialog(self, 'The values entered are not valid', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_loc_textcontrol.SetValue(str(self.loc))
                return
            # Ensure that loc is between 0 and 1
            if(self.loc < 0 or self.loc > 1):
                msg = wx.MessageDialog(self, 'Location of maximum must be between 0 and 1', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_loc_textcontrol.SetValue(str(self.loc))
                return
            self.loc = loc
            func = np.concatenate((np.linspace(0, 1, int(self.loc*(self.nmr_data.number_of_points[0]/2))), np.linspace(1, 0, int(self.nmr_data.number_of_points[0]/2) - int(self.loc*int(self.nmr_data.number_of_points[0]/2)))))
            self.line1.set_ydata(func)

        self.apodization_plot_canvas.draw()
            


    def on_zero_fill_info(self, event):
        # Create a popout window with information about zero filling

        # Create a new frame
        self.zero_fill_info_frame = wx.Frame(self, -1, 'Zero Filling Information', size=(500, 300))
        if(darkdetect.isDark() == True and platform != 'windows'):
            self.zero_fill_info_frame.SetBackgroundColour('#414141')
            colour = "RED"
        else:
            self.zero_fill_info_frame.SetBackgroundColour('White')
            colour = "BLUE"
        
        # Create a sizer to hold the box
        self.zero_fill_info_sizer_window = wx.BoxSizer(wx.VERTICAL)
        self.zero_fill_info_sizer_window.AddSpacer(10)

        # Create a sizer to hold the text
        self.zero_fill_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.zero_fill_info_sizer.AddSpacer(10)

        # Create a text box with the information
        # Zero filling information
        zero_fill_information = 'Zero filling is a method used to increase the resolution of NMR spectra. It is used to add zeros to the end of the FID to increase the number of points.\n\n It is important that the size of the data is at least doubled to prevent loss of resolution when the imaginary component of the complex data is deleted. In addition it is advised that the data is rounded to the nearest power of 2 in order to speed up the Fast Fourier Transform process.\n\nFurther advanced zero filling options can be added manually to the nmrproc.com file.'

        self.zero_fill_info_text = wx.StaticText(self.zero_fill_info_frame, -1, zero_fill_information, size=(450, 150), style=wx.ALIGN_CENTER)

        # Add the text to the sizer
        self.zero_fill_info_sizer.Add(self.zero_fill_info_text, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.zero_fill_info_sizer.AddSpacer(10)

        # Have a url to the nmrPipe help page for zero filling
        url = 'http://www.nmrscience.com/ref/nmrpipe/zf.html'
        self.zero_fill_info_url = hl.HyperLinkCtrl(self.zero_fill_info_frame, -1, 'NMRPipe Help Page for Zero Filling', URL=url)
        self.zero_fill_info_url.SetColours(colour, colour, colour)
        self.zero_fill_info_url.SetUnderlines(False, False, False)
        self.zero_fill_info_url.UpdateLink()

        # Add the url to the sizer
        self.zero_fill_info_sizer.Add(self.zero_fill_info_url, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.zero_fill_info_sizer.AddSpacer(10)

        # Add the sizer to the window sizer
        self.zero_fill_info_sizer_window.Add(self.zero_fill_info_sizer, 0, wx.ALIGN_CENTER)
        self.zero_fill_info_sizer_window.AddSpacer(10)

        # Add the window sizer to the frame
        self.zero_fill_info_frame.SetSizer(self.zero_fill_info_sizer_window)

        # Show the frame
        self.zero_fill_info_frame.Show()

    def on_linear_prediction_info(self, event):
        # Create a popout window with information about linear prediction

        # Create a new frame
        self.linear_prediction_info_frame = wx.Frame(self, -1, 'Linear Prediction Information', size=(500, 300))
        if(darkdetect.isDark() == True and platform != 'windows'):
            self.linear_prediction_info_frame.SetBackgroundColour('#414141')
            colour = "RED"
        else:
            self.linear_prediction_info_frame.SetBackgroundColour('White')
            colour = "BLUE"

        # Create a sizer to hold the box
        self.linear_prediction_info_sizer_window = wx.BoxSizer(wx.VERTICAL)
        self.linear_prediction_info_sizer_window.AddSpacer(10)

        # Create a sizer to hold the text
        self.linear_prediction_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.linear_prediction_info_sizer.AddSpacer(10)

        # Create a text box with the information
        # Linear prediction information
        linear_prediction_information = 'Linear prediction is a method used to increase the resolution of NMR spectra. It is used to predict the points of truncated FIDs (especially in indirect dimensions) and increase signal resolution.\n\n The linear prediction coefficients can be predicted using the forward FID data, backward data or an average of both directions. Then these can be used to add predicted points either before or after the current FID.\n\n Note that advanced options such as  -pred (number of predicted points) and -ord (number of predicted coefficients) can be implemented by manually added them to the nmrproc.com file.'

        self.linear_prediction_info_text = wx.StaticText(self.linear_prediction_info_frame, -1, linear_prediction_information, size=(450, 200), style=wx.ALIGN_CENTER_HORIZONTAL)
        
        # Add the text to the sizer
        self.linear_prediction_info_sizer.Add(self.linear_prediction_info_text, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.linear_prediction_info_sizer.AddSpacer(10)

        # Add a url to the nmrPipe help page
        url = 'http://www.nmrscience.com/ref/nmrpipe/lp.html'  
        self.linear_prediction_info_url = hl.HyperLinkCtrl(self.linear_prediction_info_frame, -1, 'NMRPipe Help Page for Linear Prediction', URL=url)
        self.linear_prediction_info_url.SetColours(colour, colour, colour)
        self.linear_prediction_info_url.SetUnderlines(False, False, False)
        self.linear_prediction_info_url.UpdateLink()


        # Add url to the sizer
        self.linear_prediction_info_sizer.Add(self.linear_prediction_info_url, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.linear_prediction_info_sizer.AddSpacer(10)

        # Add the sizer to the window sizer
        self.linear_prediction_info_sizer_window.Add(self.linear_prediction_info_sizer, 0, wx.ALIGN_CENTER)
        self.linear_prediction_info_sizer_window.AddSpacer(10)

        # Add the window sizer to the frame
        self.linear_prediction_info_frame.SetSizer(self.linear_prediction_info_sizer_window)

        # Show the frame
        self.linear_prediction_info_frame.Show()

        

class TwoDFrame(wx.Panel):
    def __init__(self,parent, oneDFrame):
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 0.7*self.monitorWidth
        self.height = 0.75*self.monitorHeight
        self.parent = parent
        wx.Panel.__init__(self, parent, id=wx.ID_ANY, size=(self.width, self.height))

        self.oneDFrame = oneDFrame

        # Create panel for processing dimension 1 of the data
        self.nmr_data = parent.nmr_data
        self.set_variables_dim2()
        self.create_canvas_dim2()
        self.create_menu_bar_dim2()


    def set_variables_dim2(self):
        self.set_initial_linear_prediction_variables_dim2()
        self.set_initial_apodization_variables_dim2()
        self.set_initial_zero_filling_variables_dim2()
        self.set_initial_fourier_transform_variables_dim2()
        self.set_initial_phasing_variables_dim2()
        self.set_initial_extraction_variables_dim2()
        self.set_initial_baseline_correction_variables_dim2()

        if(self.parent.load_variables == True):
            try:
                self.load_variables_from_nmrproc_com_2D()
            except:
                pass

    def load_variables_from_nmrproc_com_2D(self):
        # Open processing_parameters.txt file and load the variables from it
        file = open('processing_parameters.txt', 'r')
        lines = file.readlines()
        file.close()

        include_line = False
        for line in lines:
            if('Dimension 2' in line):
                include_line = True
                continue
            if(include_line == False):
                continue
            if(include_line == True and 'Dimension 3' in line):
                include_line = False
                break
            if(include_line == True):
                line = line.split('\n')[0]
                if(line.split(':')[0] == 'Linear Prediction'):
                    self.linear_prediction_radio_box_dim2_selection = int(line.split(': ')[1])
                if(line.split(':')[0] =='Linear Prediction Options Selection'):
                    self.linear_prediction_dim2_options_selection = int(line.split(': ')[1])
                if(line.split(':')[0] =='Linear Prediction Coefficients Selection'):
                    self.linear_prediction_dim2_coefficients_selection = int(line.split(': ')[1])
                if(line.split(':')[0] =='NUS file'):
                    self.nuslist_name_dim2 = line.split(': ')[1]
                if(line.split(':')[0] =='NUS CPU'):
                    self.smile_nus_cpu_textcontrol_dim2 = int(line.split(': ')[1])
                if(line.split(':')[0] =='NUS Iterations'):
                    self.smile_nus_iterations_textcontrol_dim2 = int(line.split(': ')[1])
                if(line.split(':')[0] =='Apodization'):
                    if('True' in line):
                        self.apodization_dim2_checkbox_value = True
                    else:
                        self.apodization_dim2_checkbox_value = False
                if(line.split(':')[0] =='Apodization Combobox Selection'):
                    self.apodization_dim2_combobox_selection = int(line.split(': ')[1])
                    self.apodization_dim2_combobox_selection_old = int(line.split(': ')[1])
                if(line.split(':')[0] =='Exponential Line Broadening'):
                    self.exponential_line_broadening_dim2 = float(line.split(': ')[1])
                if(line.split(':')[0] =='Apodization First Point Scaling'):
                    self.apodization_first_point_scaling_dim2 = float(line.split(': ')[1])
                if(line.split(':')[0] =='G1'):
                    self.g1_dim2 = float(line.split(': ')[1])
                if(line.split(':')[0] =='G2'):
                    self.g2_dim2 = float(line.split(': ')[1])
                if(line.split(':')[0] =='G3'):
                    self.g3_dim2 = float(line.split(': ')[1])
                if(line.split(':')[0] =='Offset'):
                    self.offset_dim2 = float(line.split(': ')[1])
                if(line.split(':')[0] =='End'):
                    self.end_dim2 = float(line.split(': ')[1])
                if(line.split(':')[0] =='Power'):
                    self.power_dim2 = float(line.split(': ')[1])
                if(line.split(':')[0] =='A'):
                    self.a_dim2 = float(line.split(': ')[1])
                if(line.split(':')[0] =='B'):
                    self.b_dim2 = float(line.split(': ')[1])
                if(line.split(':')[0] =='T1'):
                    self.t1_dim2 = float(line.split(': ')[1])
                if(line.split(':')[0] =='T2'):
                    self.t2_dim2 = float(line.split(': ')[1])
                if(line.split(':')[0] =='Loc'):
                    self.loc_dim2 = float(line.split(': ')[1])

                if(line.split(':')[0] =='Zero Filling'):
                    if('True' in line):
                        self.zero_filling_checkbox_dim2_value = True
                    else:
                        self.zero_filling_checkbox_dim2_value = False
                if(line.split(':')[0] =='Zero Filling Combobox Selection'):
                    self.zero_filling_dim2_combobox_selection = int(line.split(': ')[1])
                if(line.split(':')[0] =='Zero Filling Value Doubling Times'):
                    self.zero_filling_dim2_value_doubling_times = int(line.split(': ')[1])
                if(line.split(':')[0] =='Zero Filling Value Zeros to Add'):
                    self.zero_filling_dim2_value_zeros_to_add = int(line.split(': ')[1])
                if(line.split(':')[0] =='Zero Filling Value Final Data Size'):
                    self.zero_filling_dim2_value_final_data_size = int(line.split(': ')[1])
                if(line.split(':')[0] =='Zero Filling Round Checkbox'):
                    if('True' in line):
                        self.zero_filling_round_checkbox_dim2_value = True
                    else:
                        self.zero_filling_round_checkbox_dim2_value = False
                if(line.split(':')[0] =='Fourier Transform'):
                    if('True' in line):
                        self.fourier_transform_checkbox_dim2_value = True
                    else:
                        self.fourier_transform_checkbox_dim2_value = False
                if(line.split(':')[0] =='Fourier Transform Method Selection'):
                    self.ft_method_selection_dim2 = int(line.split(': ')[1])
                if(line.split(':')[0] =='Phase Correction'):
                    if('True' in line):
                        self.phase_correction_checkbox_dim2_value = True
                    else:
                        self.phase_correction_checkbox_dim2_value = False
                if(line.split(':')[0] =='Phase Correction P0'):
                    self.p0_total_dim2 = float(line.split(': ')[1])
                if(line.split(':')[0] =='Phase Correction P1'):
                    self.p1_total_dim2 = float(line.split(': ')[1])
                if(line.split(':')[0] =='F1180'):
                    if('True' in line):
                        self.f1180_dim2 = True
                    else:
                        self.f1180_dim2 = False
                if(line.split(':')[0] =='Extraction'):
                    if('True' in line):
                        self.extraction_checkbox_dim2_value = True
                    else:
                        self.extraction_checkbox_dim2_value = False
                if(line.split(':')[0] =='Extraction PPM Start'):
                    self.extraction_ppm_start_dim2 = float(line.split(': ')[1])
                if(line.split(':')[0] =='Extraction PPM End'):
                    self.extraction_ppm_end_dim2 = float(line.split(': ')[1])
                if(line.split(':')[0] =='Baseline Correction'):
                    if('True' in line):
                        self.baseline_correction_checkbox_dim2_value = True
                    else:
                        self.baseline_correction_checkbox_dim2_value = False
                if(line.split(':')[0] =='Baseline Correction Radio Box Selection'):
                    self.baseline_correction_radio_box_selection_dim2 = int(line.split(': ')[1])
                if(line.split(':')[0] =='Baseline Correction Nodes'):
                    self.baseline_correction_nodes_dim2 = int(line.split(': ')[1])
                if(line.split(':')[0] =='Baseline Correction Node List'):
                    self.baseline_correction_node_list_dim2 = line.split(': ')[1]
                if(line.split(':')[0] =='Baseline Correction Polynomial Order'):
                    self.baseline_correction_polynomial_order_dim2 = int(line.split(': ')[1])

        
    
    def set_initial_linear_prediction_variables_dim2(self):
        self.linear_prediction_radio_box_dim2_selection = 1
        self.linear_prediction_dim2_checkbox_value = False
        self.linear_prediction_dim2_options_selection = 0
        self.linear_prediction_dim2_coefficients_selection = 0
        self.linear_prediction_selection = 0

        # Check to see if the nuslist file exists in the current directory using os.path.isfile('nuslist')
        if(os.path.isfile('nuslist')):
            self.nuslist_name_dim2 = 'nuslist'
        else:
            self.nuslist_name_dim2 = ''

        self.number_of_nus_CPU_dim2 = 1
        self.nus_iterations_dim2 = 800
        self.smile_data_extension_number_dim2 = 0 # int(self.nmr_data.number_of_points[1]*1.5)



    def set_initial_apodization_variables_dim2(self):
        self.apodization_dim2_checkbox_value = True
        self.apodization_dim2_combobox_selection = 1
        self.apodization_dim2_combobox_selection_old = 1
        
        # Initial values for exponential apodization
        self.exponential_line_broadening_dim2 = 0.5
        self.apodization_first_point_scaling_dim2 = 0.5

        # Initial values for Lorentz to Gauss apodization
        self.g1_dim2 = 0.33
        self.g2_dim2 = 1
        self.g3_dim2 = 0.0

        # Initial values for Sinebell apodization
        self.offset_dim2 = 0.5
        self.end_dim2 = 0.98
        self.power_dim2 = 1.0

        # Initial values for Gauss Broadening apodization
        self.a_dim2 = 1.0
        self.b_dim2 = 1.0

        # Initial values for Trapezoid apodization
        self.t1_dim2 = int((self.nmr_data.number_of_points[1]/2)/4)
        self.t2_dim2 = int((self.nmr_data.number_of_points[1]/2)/4)

        # Initial values for Triangle apodization
        self.loc_dim2 = 0.5

    def set_initial_zero_filling_variables_dim2(self):
        self.zero_filling_dim2_checkbox_value = True
        self.zero_filling_dim2_combobox_selection = 0
        self.zero_filling_dim2_combobox_selection_old = 0
        self.zero_filling_dim2_value_doubling_times = 1
        self.zero_filling_dim2_value_zeros_to_add = 0
        self.zero_filling_dim2_value_final_data_size = 0
        self.zero_filling_dim2_round_checkbox_value = True

    def set_initial_fourier_transform_variables_dim2(self):
        self.fourier_transform_dim2_checkbox_value = True
        self.ft_method_selection_dim2 = 0 # Initially use the 'auto' method of FT as default


    def set_initial_phasing_variables_dim2(self):
        self.phasing_dim2_checkbox_value = True
        self.p0_total_dim2 = 0.0
        self.p1_total_dim2 = 0.0
        self.p0_total_dim2_old = 0.0
        self.p1_total_dim2_old = 0.0
        self.phasing_from_smile = False
        self.f1180_dim2 = False


    def set_initial_extraction_variables_dim2(self):
        self.extraction_checkbox_value_dim2 = False
        self.extraction_start_dim2 = '0.0'
        self.extraction_end_dim2 = '0.0'




    def set_initial_baseline_correction_variables_dim2(self):
        self.baseline_correction_checkbox_value_dim2 = False
        self.baseline_correction_radio_box_selection_dim2 = 0
        self.node_list_dim2 = '0,5,95,100'
        self.node_width_dim2 = '2'
        self.polynomial_order_dim2 = '4'



        

    def create_canvas_dim2(self):

        if(darkdetect.isDark() == True and platform != 'windows'):
            self.SetBackgroundColour((53, 53, 53, 255))
            matplotlib.rc('axes',edgecolor='white')
            matplotlib.rc('xtick',color='white')
            matplotlib.rc('ytick',color='white')
            matplotlib.rc('axes',labelcolor='white')
            matplotlib.rc('axes',facecolor='#282A36')
            matplotlib.rc('legend', labelcolor='white')


        else:
            self.SetBackgroundColour('White')
            matplotlib.rc('axes',edgecolor='black')
            matplotlib.rc('xtick',color='black')
            matplotlib.rc('ytick',color='black')
            matplotlib.rc('axes',labelcolor='black')
            matplotlib.rc('axes',facecolor='white')
            matplotlib.rc('legend', labelcolor='black')


        

        


    def create_menu_bar_dim2(self):
        # Create the main sizer
        self.main_sizer=wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_sizer)
        
        # Create a sizer for the processing options for the first dimension
        self.sizer_2 = wx.BoxSizer(wx.VERTICAL)
        self.sizer_2.AddSpacer(10)
    
        # Create all the sizers (allow a checkbox at the top for SMILE NUS reconstruction which will change the possible options)
        # For NUS reconstruction using SMILE need to have exact phasing paramaters (first process without NUS and then calculate phase in indirect dimension, then process again using SMILE containing exact phasing parameters)
        self.create_linear_prediction_sizer_dim2(parent=self)
        self.create_apodization_sizer_dim2(parent=self)
        self.create_zero_filling_sizer_dim2(parent=self)
        self.create_fourier_transform_sizer_dim2(parent=self)
        self.create_phase_correction_sizer_dim2(parent=self)
        self.create_extraction_sizer_dim2(parent=self)
        self.create_baseline_correction_sizer_dim2(parent=self)
        
        self.main_sizer.Add(self.sizer_2, 0, wx.EXPAND)


        self.SetSizerAndFit(self.main_sizer)
        self.Layout()


        # Get the size of the main sizer and set the window size to 1.05 times the size of the main sizer
        self.width, self.height = self.main_sizer.GetSize()
        self.parent.parent.change_frame_size(int(self.width*1.05), int(self.height*1.25))


    def create_linear_prediction_sizer_dim2(self, parent):
        # Create a sizer for the linear prediction options
        self.linear_prediction_sizer_dim2_label = wx.StaticBox(self, -1, 'Linear Prediction/SMILE NUS Reconstruction')
        self.linear_prediction_sizer_dim2 = wx.StaticBoxSizer(self.linear_prediction_sizer_dim2_label, wx.HORIZONTAL)
        self.linear_prediction_sizer_dim2.AddSpacer(10)

        # Have a radiobox for None, Linear Prediction and SMILE NUS Reconstruction
        self.linear_prediction_radio_box_dim2 = wx.RadioBox(parent, -1, '', choices=['None', 'Linear Prediction', 'SMILE NUS Reconstruction'], style=wx.RA_SPECIFY_ROWS)
        self.linear_prediction_radio_box_dim2.Bind(wx.EVT_RADIOBOX, self.on_linear_prediction_radio_box_dim2)
        self.linear_prediction_radio_box_dim2.SetSelection(self.linear_prediction_radio_box_dim2_selection)

        self.linear_prediction_sizer_dim2.Add(self.linear_prediction_radio_box_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.linear_prediction_sizer_dim2.AddSpacer(10)
        

        if(self.linear_prediction_radio_box_dim2.GetSelection() == 1):
            # Have a combobox for linear prediction options
            self.linear_prediction_options_text = wx.StaticText(parent, -1, 'Add Predicted Points:')
            self.linear_prediction_sizer_dim2.Add(self.linear_prediction_options_text, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim2.AddSpacer(5)
            self.linear_prediction_options = ['After FID', 'Before FID']
            self.linear_prediction_combobox_dim2 = wx.ComboBox(parent, -1, choices=self.linear_prediction_options, style=wx.CB_READONLY)
            self.linear_prediction_combobox_dim2.SetSelection(self.linear_prediction_dim2_options_selection)
            self.linear_prediction_combobox_dim2.Bind(wx.EVT_COMBOBOX, self.on_linear_prediction_combobox_dim2)
            self.linear_prediction_sizer_dim2.Add(self.linear_prediction_combobox_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim2.AddSpacer(10)
            # Have a combobox of predicted coefficient options
            self.linear_prediction_coefficients_text = wx.StaticText(parent, -1, 'Predicted Coefficients:')
            self.linear_prediction_sizer_dim2.Add(self.linear_prediction_coefficients_text, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim2.AddSpacer(5)
            self.linear_prediction_coefficients_options = ['Forward', 'Backward', 'Both']
            self.linear_prediction_coefficients_combobox_dim2 = wx.ComboBox(parent, -1, choices=self.linear_prediction_coefficients_options, style=wx.CB_READONLY)
            self.linear_prediction_coefficients_combobox_dim2.SetSelection(self.linear_prediction_dim2_coefficients_selection)
            self.linear_prediction_coefficients_combobox_dim2.Bind(wx.EVT_COMBOBOX, self.on_linear_prediction_combobox_coefficients_dim2)
            self.linear_prediction_sizer_dim2.Add(self.linear_prediction_coefficients_combobox_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim2.AddSpacer(10)
        elif(self.linear_prediction_radio_box_dim2.GetSelection()==2):
            # Have a set of options for SMILE NUS processing

            # NUS file
            self.smile_nus_file_text = wx.StaticText(parent, -1, 'NUS File:')
            self.linear_prediction_sizer_dim2.Add(self.smile_nus_file_text, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim2.AddSpacer(5)

            self.smile_nus_file_textcontrol_dim2 = wx.TextCtrl(parent, -1, self.nuslist_name_dim2, size=(100, 20))
            self.smile_nus_file_textcontrol_dim2.Bind(wx.EVT_TEXT, self.on_smile_nus_file_textcontrol_dim2)

            self.linear_prediction_sizer_dim2.Add(self.smile_nus_file_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim2.AddSpacer(10)

            # # Zero order phase correction
            # self.smile_nus_p0_text = wx.StaticText(parent, -1, 'Zero Order Phase Correction (p0):')
            # self.linear_prediction_sizer_dim2.Add(self.smile_nus_p0_text, 0, wx.ALIGN_CENTER_VERTICAL)
            # self.linear_prediction_sizer_dim2.AddSpacer(5)
            # self.smile_nus_p0_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.p0_total_dim2), size=(50, 20))
            # self.smile_nus_p0_textcontrol_dim2.Bind(wx.EVT_TEXT, self.on_smile_nus_p0_textcontrol_dim2)
            # self.linear_prediction_sizer_dim2.Add(self.smile_nus_p0_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            # self.linear_prediction_sizer_dim2.AddSpacer(10)

            # # First order phase correction
            # self.smile_nus_p1_text = wx.StaticText(parent, -1, 'First Order Phase Correction (p1):')
            # self.linear_prediction_sizer_dim2.Add(self.smile_nus_p1_text, 0, wx.ALIGN_CENTER_VERTICAL)
            # self.linear_prediction_sizer_dim2.AddSpacer(5)
            # self.smile_nus_p1_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.p1_total_dim2), size=(50, 20))
            # self.smile_nus_p1_textcontrol_dim2.Bind(wx.EVT_TEXT, self.on_smile_nus_p1_textcontrol_dim2)
            # self.linear_prediction_sizer_dim2.Add(self.smile_nus_p1_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            # self.linear_prediction_sizer_dim2.AddSpacer(10)

            # Number of points to add to the data
            self.smile_nus_extension_text = wx.StaticText(parent, -1, 'Data extension:')
            self.linear_prediction_sizer_dim2.Add(self.smile_nus_extension_text, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim2.AddSpacer(5)
            self.smile_nus_extension_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.smile_data_extension_number_dim2), size=(50, 20))
            self.smile_nus_extension_textcontrol_dim2.Bind(wx.EVT_TEXT, self.on_smile_nus_extension_textcontrol_dim2)
            self.linear_prediction_sizer_dim2.Add(self.smile_nus_extension_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim2.AddSpacer(10)



            # Number of CPU's
            self.smile_nus_cpu_text = wx.StaticText(parent, -1, 'Number of CPU\'s:')
            self.linear_prediction_sizer_dim2.Add(self.smile_nus_cpu_text, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim2.AddSpacer(5)
            self.smile_nus_cpu_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.number_of_nus_CPU_dim2), size=(30, 20))
            self.smile_nus_cpu_textcontrol_dim2.Bind(wx.EVT_TEXT, self.on_smile_nus_cpu_textcontrol_dim2)
            self.linear_prediction_sizer_dim2.Add(self.smile_nus_cpu_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim2.AddSpacer(10)

            # Number of iterations
            self.smile_nus_iterations_text = wx.StaticText(parent, -1, 'Number of Iterations:')
            self.linear_prediction_sizer_dim2.Add(self.smile_nus_iterations_text, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim2.AddSpacer(5)
            self.smile_nus_iterations_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.nus_iterations_dim2), size=(30, 20))
            self.smile_nus_iterations_textcontrol_dim2.Bind(wx.EVT_TEXT, self.on_smile_nus_iterations_textcontrol_dim2)
            self.linear_prediction_sizer_dim2.Add(self.smile_nus_iterations_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)

            

        # Have a button showing information on linear prediction
        self.linear_prediction_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))
        self.linear_prediction_info.Bind(wx.EVT_BUTTON, self.on_linear_prediction_info_dim2)
        self.linear_prediction_sizer_dim2.AddSpacer(10)
        self.linear_prediction_sizer_dim2.Add(self.linear_prediction_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer_2.Add(self.linear_prediction_sizer_dim2)
        self.sizer_2.AddSpacer(10)   

    def on_linear_prediction_combobox_dim2(self, event):
        # Get the selection from the combobox and update the linear prediction options
        self.linear_prediction_dim2_options_selection = self.linear_prediction_combobox_dim2.GetSelection()

    def on_linear_prediction_combobox_coefficients_dim2(self, event):
        # Get the selection from the combobox and update the linear prediction options
        self.linear_prediction_dim2_coefficients_selection = self.linear_prediction_coefficients_combobox_dim2.GetSelection()

    def on_smile_nus_file_textcontrol_dim2(self, event):
        # Get the value from the textcontrol
        self.nuslist_name_dim2 = self.smile_nus_file_textcontrol_dim2.GetValue()

    def on_smile_nus_p0_textcontrol_dim2(self, event):
        self.p0_total_dim2 = self.smile_nus_p0_textcontrol_dim2.GetValue()
        self.phasing_from_smile = True
        # Update the phasing values in the phasing section too
        self.phase_correction_p0_textcontrol_dim2.SetValue(str(self.p0_total_dim2))
        self.phasing_from_smile = False

    def on_smile_nus_p1_textcontrol_dim2(self, event):
        self.p1_total_dim2 = self.smile_nus_p1_textcontrol_dim2.GetValue()
        self.phasing_from_smile = True
        # Update the phasing values in the phasing section too
        self.phase_correction_p1_textcontrol_dim2.SetValue(str(self.p1_total_dim2))
        self.phasing_from_smile = False


    def on_smile_nus_extension_textcontrol_dim2(self, event):
        if(self.smile_nus_extension_textcontrol_dim2.GetValue() != ''):
            try:
                self.smile_data_extension_number_dim2 = int(self.smile_nus_extension_textcontrol_dim2.GetValue())
            except:
                msg = wx.MessageDialog(self, 'The value entered for NUS data extension not a valid integer', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.smile_nus_extension_textcontrol_dim2.SetValue(str(self.smile_data_extension_number_dim2))
                return
        else:
            self.smile_data_extension_number_dim2 = self.smile_nus_extension_textcontrol_dim2.GetValue()
    
    def on_smile_nus_cpu_textcontrol_dim2(self, event):
        if(self.smile_nus_cpu_textcontrol_dim2.GetValue() != ''):
            try:
                self.number_of_nus_CPU_dim2 = int(self.smile_nus_cpu_textcontrol_dim2.GetValue())
            except:
                msg = wx.MessageDialog(self, 'The value entered for number of CPU\'s is not a valid integer', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.smile_nus_cpu_textcontrol_dim2.SetValue(str(self.number_of_nus_CPU_dim2))
                return

    def on_smile_nus_iterations_textcontrol_dim2(self, event):
        if(self.smile_nus_iterations_textcontrol_dim2.GetValue() != ''):
            try:
                self.nus_iterations_dim2 = int(self.smile_nus_iterations_textcontrol_dim2.GetValue())
            except:
                msg = wx.MessageDialog(self, 'The value entered for number of iterations is not a valid integer', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.smile_nus_iterations_textcontrol_dim2.SetValue(str(self.nus_iterations_dim2))
                return

    



    def on_linear_prediction_info_dim2(self,event):
        # Create a popout window with information about linear prediction

        # Create a new frame
        self.linear_prediction_info_frame = wx.Frame(self, -1, 'Linear Prediction / SMILE Information', size=(500, 500))
        if(darkdetect.isDark() == True and platform != 'windows'):
            self.linear_prediction_info_frame.SetBackgroundColour('#414141')
            colour = "RED"
        else:
            self.linear_prediction_info_frame.SetBackgroundColour('White')
            colour = "BLUE"

        # Create a sizer to hold the box
        self.linear_prediction_info_sizer_window = wx.BoxSizer(wx.VERTICAL)
        self.linear_prediction_info_sizer_window.AddSpacer(10)

        # Create a sizer to hold the text
        self.linear_prediction_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.linear_prediction_info_sizer.AddSpacer(10)

        # Create a text box with the information
        # Linear prediction information
        linear_prediction_information = 'Linear prediction is a method used to increase the resolution of NMR spectra. It is used to predict the points of truncated FIDs (especially in indirect dimensions) and increase signal resolution.\n\n The linear prediction coefficients can be predicted using the forward FID data, backward data or an average of both directions. Then these can be used to add predicted points either before or after the current FID.\n\n Note that advanced options such as  -pred (number of predicted points) and -ord (number of predicted coefficients) can be implemented by manually added them to the nmrproc.com file.'

        self.linear_prediction_info_text = wx.StaticText(self.linear_prediction_info_frame, -1, linear_prediction_information, size=(450, 200), style=wx.ALIGN_CENTER_HORIZONTAL)
        
        # Add the text to the sizer
        self.linear_prediction_info_sizer.Add(self.linear_prediction_info_text, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.linear_prediction_info_sizer.AddSpacer(10)

        # Add a url to the nmrPipe help page
        url = 'http://www.nmrscience.com/ref/nmrpipe/lp.html'  
        self.linear_prediction_info_url = hl.HyperLinkCtrl(self.linear_prediction_info_frame, -1, 'NMRPipe Help Page for Linear Prediction', URL=url)
        self.linear_prediction_info_url.SetColours(colour, colour, colour)
        self.linear_prediction_info_url.SetUnderlines(False, False, False)
        self.linear_prediction_info_url.UpdateLink()


        # Add url to the sizer
        self.linear_prediction_info_sizer.Add(self.linear_prediction_info_url, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.linear_prediction_info_sizer.AddSpacer(10)

        # Add the sizer to the window sizer
        self.linear_prediction_info_sizer_window.Add(self.linear_prediction_info_sizer, 0, wx.ALIGN_CENTER)
        self.linear_prediction_info_sizer_window.AddSpacer(10)


        # Have text to explain SMILE NUS reconstruction
        smile_nus_text = 'SMILE NUS reconstruction is a method used to reconstruct non-uniformly sampled data. The NUS file is a list of points that have been sampled in the FID.\nThe number of CPU\'s (default=1) is the number of cores that will be used to perform the reconstruction and the number of iterations can be changed to improve the accuracy (default=800).\n Furthermore, in order for accurate SMILE reconstruction, the correct zero (p0) and first (p1) order phase correction values need to be inputted.'

        self.smile_nus_text = wx.StaticText(self.linear_prediction_info_frame, -1, smile_nus_text, size=(450, 200), style=wx.ALIGN_CENTER_HORIZONTAL)
        self.linear_prediction_info_sizer_window.Add(self.smile_nus_text, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.linear_prediction_info_sizer_window.AddSpacer(10)




        # Add the window sizer to the frame
        self.linear_prediction_info_frame.SetSizer(self.linear_prediction_info_sizer_window)

        # Show the frame
        self.linear_prediction_info_frame.Show()

    
    def on_linear_prediction_radio_box_dim2(self, event):
        # Get the selection from the radio box and update the linear prediction options
        self.linear_prediction_radio_box_dim2_selection = self.linear_prediction_radio_box_dim2.GetSelection()
        
        # Remove all the old sizers and replot


        if(self.apodization_dim2_combobox_selection_old == 1):
            # Remove the previous textcontrols
            
            self.apodization_sizer_dim2.Detach(self.apodization_line_broadening_label)
            self.apodization_line_broadening_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_line_broadening_textcontrol_dim2)
            self.apodization_line_broadening_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_textcontrol_dim2)
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)

        elif(self.apodization_dim2_combobox_selection_old == 2):
            self.apodization_sizer_dim2.Detach(self.apodization_g1_label)
            self.apodization_g1_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g1_textcontrol_dim2)
            self.apodization_g1_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g2_label)
            self.apodization_g2_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g2_textcontrol_dim2)
            self.apodization_g2_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g3_label)
            self.apodization_g3_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g3_textcontrol_dim2)
            self.apodization_g3_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_textcontrol_dim2)
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)

        elif(self.apodization_dim2_combobox_selection_old == 3):
            self.apodization_sizer_dim2.Detach(self.apodization_offset_label)
            self.apodization_offset_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_offset_textcontrol)
            self.apodization_offset_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_end_label)
            self.apodization_end_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_end_textcontrol_dim2)
            self.apodization_end_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_power_label)
            self.apodization_power_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_power_textcontrol_dim2)
            self.apodization_power_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_textcontrol_dim2)
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        elif(self.apodization_dim2_combobox_selection_old == 4):
            self.apodization_sizer_dim2.Detach(self.apodization_a_label)
            self.apodization_a_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_a_textcontrol_dim2)
            self.apodization_a_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_b_label)
            self.apodization_b_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_b_textcontrol_dim2)
            self.apodization_b_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_textcontrol_dim2)
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        elif(self.apodization_dim2_combobox_selection_old == 5):
            self.apodization_sizer_dim2.Detach(self.apodization_t1_label)
            self.apodization_t1_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_t1_textcontrol_dim2)
            self.apodization_t1_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_t2_label)
            self.apodization_t2_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_t2_textcontrol_dim2)
            self.apodization_t2_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_textcontrol_dim2)
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        elif(self.apodization_dim2_combobox_selection_old == 6):
            self.apodization_sizer_dim2.Detach(self.apodization_loc_label)
            self.apodization_loc_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_loc_textcontrol_dim2)
            self.apodization_loc_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_textcontrol_dim2)
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        elif(self.apodization_dim2_combobox_selection_old == 0):
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        


        self.apodization_sizer_dim2.Detach(self.apodization_checkbox_dim2)
        self.apodization_checkbox_dim2.Destroy()
        self.apodization_sizer_dim2.Detach(self.apodization_combobox_dim2)
        self.apodization_combobox_dim2.Hide()
       


        

        # Delete the current apodization sizer and then create a new one
        
        self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
        self.apodization_plot_ax_dim2.clear()
        self.apodization_plot_ax_dim2.clear()
        self.apodization_plot_sizer_dim2.Clear(True)

        self.sizer_2.Remove(self.apodization_sizer_dim2)
        # self.apodization_sizer.Clear(delete_windows=True)

 


        



        # Remove the linear prediction sizers
        self.linear_prediction_sizer_dim2.Clear(delete_windows=True)
        # self.sizer_1.Remove(self.linear_prediction_sizer)
        

        # self.sizer_1.Remove(self.solvent_suppression_sizer)

        
        self.sizer_2.Clear(delete_windows=True)

           
        self.create_menu_bar_dim2()
        self.Refresh()
        self.Update()
        self.Layout()


    def create_apodization_sizer_dim2(self, parent):

        # Create a box for apodization options
        self.apodization_box_dim2 = wx.StaticBox(parent, -1, 'Apodization')
        self.apodization_sizer_dim2 = wx.StaticBoxSizer(self.apodization_box_dim2, wx.HORIZONTAL)
        self.apodization_checkbox_dim2 = wx.CheckBox(parent, -1, 'Apply apodization')
        self.apodization_checkbox_dim2.Bind(wx.EVT_CHECKBOX, self.on_apodization_checkbox_dim2)
        self.apodization_checkbox_dim2.SetValue(self.apodization_dim2_checkbox_value)
        self.apodization_sizer_dim2.Add(self.apodization_checkbox_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.apodization_sizer_dim2.AddSpacer(10)
        # Have a combobox for apodization options
        self.apodization_options_dim2 = ['None', 'Exponential', 'Lorentz to Gauss', 'Sinebell', 'Gauss Broadening', 'Trapazoid', 'Triangle']
        self.apodization_combobox_dim2 = wx.ComboBox(parent, -1, choices=self.apodization_options_dim2, style=wx.CB_READONLY)
        self.apodization_combobox_dim2.SetSelection(self.apodization_dim2_combobox_selection)
        self.apodization_combobox_dim2.Bind(wx.EVT_COMBOBOX, self.on_apodization_combobox_dim2)
        self.apodization_sizer_dim2.Add(self.apodization_combobox_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.apodization_sizer_dim2.AddSpacer(10)
        if(self.apodization_dim2_combobox_selection == 1):
            # Have a textcontrol for the line broadening
            self.apodization_line_broadening_label = wx.StaticText(parent, -1, 'Line Broadening (Hz):')
            self.apodization_sizer_dim2.Add(self.apodization_line_broadening_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_line_broadening_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.exponential_line_broadening_dim2), size=(30, 20))
            self.apodization_line_broadening_textcontrol_dim2.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2)
            self.apodization_sizer_dim2.Add(self.apodization_line_broadening_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(parent, -1, 'First Point Scaling:')
            self.apodization_sizer_dim2.Add(self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_first_point_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.apodization_first_point_scaling_dim2), size=(30, 20))
            self.apodization_sizer_dim2.Add(self.apodization_first_point_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim2.AddSpacer(10)
        elif(self.apodization_dim2_combobox_selection == 2):
            # Have a textcontrol for the g1 value
            self.apodization_g1_label = wx.StaticText(parent, -1, 'Inverse Lorentzian (Hz):')
            self.apodization_sizer_dim2.Add(self.apodization_g1_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_g1_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.g1_dim2), size=(40, 20))
            self.apodization_g1_textcontrol_dim2.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2)
            self.apodization_sizer_dim2.Add(self.apodization_g1_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the g2 value
            self.apodization_g2_label = wx.StaticText(parent, -1, 'Gaussian Broadening (Hz):')
            self.apodization_sizer_dim2.Add(self.apodization_g2_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_g2_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.g2_dim2), size=(40, 20))
            self.apodization_g2_textcontrol_dim2.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2)
            self.apodization_sizer_dim2.Add(self.apodization_g2_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the g3 value
            self.apodization_g3_label = wx.StaticText(parent, -1, 'Gaussian Shift:')
            self.apodization_sizer_dim2.Add(self.apodization_g3_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_g3_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.g3_dim2), size=(40, 20))
            self.apodization_g3_textcontrol_dim2.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2)
            self.apodization_sizer_dim2.Add(self.apodization_g3_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(parent, -1, 'First Point Scaling:')
            self.apodization_sizer_dim2.Add(self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_first_point_textcontrol_dim2 = wx.TextCtrl(self, -1, str(self.apodization_first_point_scaling_dim2), size=(30, 20))
            self.apodization_sizer_dim2.Add(self.apodization_first_point_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim2.AddSpacer(10)
        elif(self.apodization_dim2_combobox_selection == 3):
            # Have a textcontrol for the offset value
            self.apodization_offset_label = wx.StaticText(parent, -1, 'Offset (\u03c0):')
            self.apodization_sizer_dim2.Add(self.apodization_offset_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_offset_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.offset_dim2), size=(40, 20))
            self.apodization_offset_textcontrol_dim2.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2)
            self.apodization_sizer_dim2.Add(self.apodization_offset_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the end value
            self.apodization_end_label = wx.StaticText(parent, -1, 'End (\u03c0):')
            self.apodization_sizer_dim2.Add(self.apodization_end_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_end_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.end_dim2), size=(40, 20))
            self.apodization_end_textcontrol_dim2.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2)
            self.apodization_sizer_dim2.Add(self.apodization_end_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the power value
            self.apodization_power_label = wx.StaticText(parent, -1, 'Power:')
            self.apodization_sizer_dim2.Add(self.apodization_power_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_power_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.power_dim2), size=(30, 20))
            self.apodization_power_textcontrol_dim2.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2)
            self.apodization_sizer_dim2.Add(self.apodization_power_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(parent, -1, 'First Point Scaling:')
            self.apodization_sizer_dim2.Add(self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_first_point_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.apodization_first_point_scaling_dim2), size=(30, 20))
            self.apodization_sizer_dim2.Add(self.apodization_first_point_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim2.AddSpacer(10)
        elif(self.apodization_dim2_combobox_selection == 4):
            # Have a textcontrol for the a value
            self.apodization_a_label = wx.StaticText(parent, -1, 'Line Broadening (Hz):')
            self.apodization_sizer_dim2.Add(self.apodization_a_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_a_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.a_dim2), size=(40, 20))
            self.apodization_a_textcontrol_dim2.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2)
            self.apodization_sizer_dim2.Add(self.apodization_a_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the b value
            self.apodization_b_label = wx.StaticText(parent, -1, 'Gaussian Broadening (Hz):')
            self.apodization_sizer_dim2.Add(self.apodization_b_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_b_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.b_dim2), size=(40, 20))
            self.apodization_b_textcontrol_dim2.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2)
            self.apodization_sizer_dim2.Add(self.apodization_b_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(parent, -1, 'First Point Scaling:')
            self.apodization_sizer_dim2.Add(self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_first_point_textcontrol_dim2 = wx.TextCtrl(self, -1, str(self.apodization_first_point_scaling_dim2), size=(30, 20))
            self.apodization_sizer_dim2.Add(self.apodization_first_point_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim2.AddSpacer(10)
        elif(self.apodization_dim2_combobox_selection == 5):
            # Have a textcontrol for the t1 value
            self.apodization_t1_label = wx.StaticText(parent, -1, 'Ramp up points:')
            self.apodization_sizer_dim2.Add(self.apodization_t1_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_t1_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.t1_dim2), size=(50, 20))
            self.apodization_t1_textcontrol_dim2.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2)
            self.apodization_sizer_dim2.Add(self.apodization_t1_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the t2 value
            self.apodization_t2_label = wx.StaticText(parent, -1, 'Ramp down points:')
            self.apodization_sizer_dim2.Add(self.apodization_t2_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_t2_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.t2_dim2), size=(50, 20))
            self.apodization_t2_textcontrol_dim2.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2)
            self.apodization_sizer_dim2.Add(self.apodization_t2_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(parent, -1, 'First Point Scaling:')
            self.apodization_sizer_dim2.Add(self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_first_point_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.apodization_first_point_scaling_dim2), size=(30, 20))
            self.apodization_sizer_dim2.Add(self.apodization_first_point_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim2.AddSpacer(10)
        elif(self.apodization_dim2_combobox_selection == 6):
            # Have a textcontrol for the loc value
            self.apodization_loc_label = wx.StaticText(parent, -1, 'Location of maximum:')
            self.apodization_sizer_dim2.Add(self.apodization_loc_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_loc_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.loc_dim2), size=(40, 20))
            self.apodization_loc_textcontrol_dim2.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim2)
            self.apodization_sizer_dim2.Add(self.apodization_loc_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim2.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(parent, -1, 'First Point Scaling:')
            self.apodization_sizer_dim2.Add(self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_first_point_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.apodization_first_point_scaling_dim2), size=(30, 20))
            self.apodization_sizer_dim2.Add(self.apodization_first_point_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim2.AddSpacer(10)

            
        # Have a button for information on currently selected apodization containing unicode i in a circle
        self.apodization_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))
        self.apodization_info.Bind(wx.EVT_BUTTON, self.oneDFrame.on_apodization_info)
        self.apodization_sizer_dim2.Add(self.apodization_info, 0, wx.ALIGN_CENTER_VERTICAL)

        # Have a mini plots of the apodization function along with the FID first slice
        self.plot_window_function_dim2()




        self.sizer_2.Add(self.apodization_sizer_dim2)
        self.sizer_2.AddSpacer(10)

    def on_apodization_checkbox_dim2(self, event):
        # Get the selection from the checkbox
        self.apodization_dim2_checkbox_selection = self.apodization_checkbox_dim2.GetValue()



    def on_apodization_textcontrol_dim2(self, event):
        # If the user presses enter, update the plot
        keycode = event.GetKeyCode()
        if(keycode == wx.WXK_RETURN):
            self.update_window_function_plot_dim2()
        event.Skip()


    def plot_window_function_dim2(self):
        self.apodization_plot_sizer_dim2 = wx.BoxSizer(wx.VERTICAL)
        if(darkdetect.isDark() == True and platform!='windows'):
            self.apodization_plot_dim2 = Figure(figsize=(1, 0.5),facecolor='#3b3b3b')
        else:
            self.apodization_plot_dim2 = Figure(figsize=(1, 0.5),facecolor='#e6e6e7')
        self.apodization_plot_ax_dim2 = self.apodization_plot_dim2.add_subplot(111)
        # self.apodization_plot_ax.set_axis_off()
        
        if(darkdetect.isDark() == True and platform!='windows'):
            self.apodization_plot_ax_dim2.set_facecolor('#3b3b3b')
        else:
            self.apodization_plot_ax_dim2.set_facecolor('#e6e6e7')

        self.apodization_plot_ax_dim2.set_xticks([])
        self.apodization_plot_ax_dim2.set_yticks([])

        # If the apodization function is None, make remove the axes of the plot
        if(self.apodization_dim2_combobox_selection == 0):
            self.apodization_plot_ax_dim2.spines['top'].set_visible(False)
            self.apodization_plot_ax_dim2.spines['right'].set_visible(False)
            self.apodization_plot_ax_dim2.spines['bottom'].set_visible(False)
            self.apodization_plot_ax_dim2.spines['left'].set_visible(False)

        
        # If have a pseudo axis and the pseudo axis is the 2nd dimension
        if(self.nmr_data.pseudo_axis == True and self.nmr_data.index == 1):
            x=np.linspace(0, (self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2], int(self.nmr_data.number_of_points[2]/2))
            if(self.apodization_dim2_combobox_selection == 1):
                # Exponential window function
                self.line1, = self.apodization_plot_ax_dim2.plot(x, np.exp(-(np.pi*x*self.exponential_line_broadening_dim2)), color='#1f77b4')
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(0, (self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2])
                
            elif(self.apodization_dim2_combobox_selection == 2):
                # Lorentz to Gauss window function
                e = np.pi * (self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2] * self.g1_dim2
                g = 0.6 * np.pi * self.g2_dim2 * (self.g3_dim2 * ((self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2] - 1) - x)
                func = np.exp(e - g * g)
                self.line1, = self.apodization_plot_ax_dim2.plot(x, func, color='#1f77b4')
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(0, (self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2])
            elif(self.apodization_dim2_combobox_selection == 3):
                # Sinebell window function
                func = np.sin((np.pi*self.offset_dim2 + np.pi*(self.end_dim2-self.offset_dim2)*x)/((((self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2]))))**self.power_dim2
                self.line1, = self.apodization_plot_ax_dim2.plot(x, func, color='#1f77b4')
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(0, (self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2])
            elif(self.apodization_dim2_combobox_selection == 4):
                # Gauss broadening window function
                func = np.exp(-self.a_dim2*(x**2) - self.b_dim2*x)
                self.line1, = self.apodization_plot_ax_dim2.plot(x, func, color='#1f77b4')
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(0, self.nmr_data.number_of_points[2]/self.nmr_data.spectral_width[2])
            elif(self.apodization_dim2_combobox_selection == 5):
                # Trapazoid window function
                func = np.concatenate((np.linspace(0, 1, int(self.t1_dim2)), np.ones(int(self.nmr_data.number_of_points[2]/2) - int(self.t1_dim2) - int(self.t1_dim2)),np.linspace(1, 0, int(self.t2_dim2))))
                self.line1, = self.apodization_plot_ax_dim2.plot(x, func, color='#1f77b4')
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(0, self.nmr_data.number_of_points[2]/self.nmr_data.spectral_width[2])
            elif(self.apodization_dim2_combobox_selection == 6):
                # Triangle window function
                func = np.concatenate((np.linspace(0, 1, int(self.loc_dim2*(self.nmr_data.number_of_points[2]/2))), np.linspace(1, 0, int((1-self.loc_dim2)*(self.nmr_data.number_of_points[2]/2)))))
                self.line1, = self.apodization_plot_ax_dim2.plot(x, func, color='#1f77b4')
                
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(0, (self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2])



                
            self.apodization_plot_ax_dim2.set_xlim(0, self.nmr_data.number_of_points[2]/self.nmr_data.spectral_width[2])

        else:
            x=np.linspace(0, (self.nmr_data.number_of_points[1]/2)/self.nmr_data.spectral_width[1], int(self.nmr_data.number_of_points[1]/2))
            if(self.apodization_dim2_combobox_selection == 1):
                # Exponential window function
                self.line1, = self.apodization_plot_ax_dim2.plot(x, np.exp(-(np.pi*x*self.exponential_line_broadening_dim2)), color='#1f77b4')
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(0, (self.nmr_data.number_of_points[1]/2)/self.nmr_data.spectral_width[1])
                
            elif(self.apodization_dim2_combobox_selection == 2):
                # Lorentz to Gauss window function
                e = np.pi * (self.nmr_data.number_of_points[1]/2)/self.nmr_data.spectral_width[1] * self.g1_dim2
                g = 0.6 * np.pi * self.g2_dim2 * (self.g3_dim2 * ((self.nmr_data.number_of_points[1]/2)/self.nmr_data.spectral_width[1] - 1) - x)
                func = np.exp(e - g * g)
                self.line1, = self.apodization_plot_ax_dim2.plot(x, func, color='#1f77b4')
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(0, (self.nmr_data.number_of_points[1]/2)/self.nmr_data.spectral_width[1])
            elif(self.apodization_dim2_combobox_selection == 3):
                # Sinebell window function
                func = np.sin((np.pi*self.offset_dim2 + np.pi*(self.end_dim2-self.offset_dim2)*x)/((((self.nmr_data.number_of_points[1]/2)/self.nmr_data.spectral_width[1]))))**self.power_dim2
                self.line1, = self.apodization_plot_ax_dim2.plot(x, func, color='#1f77b4')
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(0, (self.nmr_data.number_of_points[1]/2)/self.nmr_data.spectral_width[1])
            elif(self.apodization_dim2_combobox_selection == 4):
                # Gauss broadening window function
                func = np.exp(-self.a_dim2*(x**2) - self.b_dim2*x)
                self.line1, = self.apodization_plot_ax_dim2.plot(x, func, color='#1f77b4')
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(0, self.nmr_data.number_of_points[1]/self.nmr_data.spectral_width[1])
            elif(self.apodization_dim2_combobox_selection == 5):
                # Trapazoid window function
                func = np.concatenate((np.linspace(0, 1, int(self.t1_dim2)), np.ones(int(self.nmr_data.number_of_points[1]/2) - int(self.t1_dim2) - int(self.t1_dim2)),np.linspace(1, 0, int(self.t2_dim2))))
                self.line1, = self.apodization_plot_ax_dim2.plot(x, func, color='#1f77b4')
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(0, self.nmr_data.number_of_points[1]/self.nmr_data.spectral_width[1])
            elif(self.apodization_dim2_combobox_selection == 6):
                # Triangle window function
                func = np.concatenate((np.linspace(0, 1, int(self.loc_dim2*(self.nmr_data.number_of_points[1]/2))), np.linspace(1, 0, int((1-self.loc_dim2)*(self.nmr_data.number_of_points[1]/2)))))
                self.line1, = self.apodization_plot_ax_dim2.plot(x, func, color='#1f77b4')
                
                self.apodization_plot_ax_dim2.set_ylim(0, 1.5)
                self.apodization_plot_ax_dim2.set_xlim(0, (self.nmr_data.number_of_points[1]/2)/self.nmr_data.spectral_width[1])



                
            self.apodization_plot_ax_dim2.set_xlim(0, self.nmr_data.number_of_points[1]/self.nmr_data.spectral_width[1])

        self.apodization_plot_canvas = FigCanvas(self, -1, self.apodization_plot_dim2)
        self.apodization_plot_sizer_dim2.Add(self.apodization_plot_canvas, 0, wx.EXPAND)


        self.apodization_sizer_dim2.Add(self.apodization_plot_sizer_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.apodization_sizer_dim2.AddSpacer(10)


    def update_window_function_plot_dim2(self):
        if(self.nmr_data.pseudo_axis == True and self.nmr_data.index == 1):
            x=np.linspace(0, (self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2], int(self.nmr_data.number_of_points[2]/2))
        else:
            x=np.linspace(0, (self.nmr_data.number_of_points[1]/2)/self.nmr_data.spectral_width[1], int(self.nmr_data.number_of_points[1]/2))
        try:
            c = float(self.apodization_first_point_textcontrol_dim2.GetValue())
            self.apodization_first_point_scaling_dim2 = c
        except:
            # Give a popout window saying that the values are not valid
            msg = wx.MessageDialog(self, 'The value entered for apodization first point scaling is not valid (use 0.5 or 1.0)', 'Error', wx.OK | wx.ICON_ERROR)
            msg.ShowModal()
            msg.Destroy()
            self.apodization_first_point_textcontrol_dim2.SetValue(str(self.apodization_first_point_scaling_dim2))
            return
        if(c != 0.5 and c != 1.0):
            msg = wx.MessageDialog(self, 'The value entered for apodization first point scaling is not valid (use 0.5 or 1.0)', 'Error', wx.OK | wx.ICON_ERROR)
            msg.ShowModal()
            msg.Destroy()
            self.apodization_first_point_textcontrol_dim2.SetValue(str(self.apodization_first_point_scaling_dim2))
            return
        self.apodization_first_point_scaling = c
        if(self.apodization_dim2_combobox_selection==1):
            try:
                em = float(self.apodization_line_broadening_textcontrol_dim2.GetValue())
            except:
                # Give a popout window saying that the values are not valid
                msg = wx.MessageDialog(self, 'The values entered are not valid', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_line_broadening_textcontrol_dim2.SetValue(str(self.exponential_line_broadening_dim2))
                return
            self.exponential_line_broadening_dim2 = em
            
            self.line1.set_ydata(np.exp(-(np.pi*x*self.exponential_line_broadening_dim2)))
        elif(self.apodization_dim2_combobox_selection==2):
            try:
                g1 = float(self.apodization_g1_textcontrol_dim2.GetValue())
                g2 = float(self.apodization_g2_textcontrol_dim2.GetValue())
                g3 = float(self.apodization_g3_textcontrol_dim2.GetValue())
            except:
                # Give a popout window saying that the values are not valid
                msg = wx.MessageDialog(self, 'The values entered are not valid', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_g1_textcontrol_dim2.SetValue(str(self.g1_dim2))
                self.apodization_g2_textcontrol_dim2.SetValue(str(self.g2_dim2))
                self.apodization_g3_textcontrol_dim2.SetValue(str(self.g3_dim2))
                return
            # Check to see if g3 is between 0 and 1
            if(g3 < 0 or g3 > 1):
                msg = wx.MessageDialog(self, 'Gaussian shift must be between 0 and 1', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_g3_textcontrol_dim2.SetValue(str(self.g3_dim2))
                return
            self.g1_dim2 = g1
            self.g2_dim2 = g2
            self.g3_dim2 = g3
            if(self.nmr_data.pseudo_axis == True and self.nmr_data.index == 1):
                e = np.pi * (self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2] * self.g1_dim2
                g = 0.6 * np.pi * self.g2_dim2 * (self.g3_dim2 * ((self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2] - 1) - x)
                func = np.exp(e - g * g)
                self.line1.set_ydata(func)

                self.apodization_plot_ax_dim2.set_xlim(0, (self.nmr_data.number_of_points[1]/2)/self.nmr_data.spectral_width[1])

            else:
                e = np.pi * (self.nmr_data.number_of_points[1]/2)/self.nmr_data.spectral_width[1] * self.g1_dim2
                g = 0.6 * np.pi * self.g2_dim2 * (self.g3_dim2 * ((self.nmr_data.number_of_points[1]/2)/self.nmr_data.spectral_width[1] - 1) - x)
                func = np.exp(e - g * g)
                self.line1.set_ydata(func)

                self.apodization_plot_ax_dim2.set_xlim(0, (self.nmr_data.number_of_points[1]/2)/self.nmr_data.spectral_width[1])


        elif(self.apodization_dim2_combobox_selection==3):
            try:
                offset = float(self.apodization_offset_textcontrol_dim2.GetValue())
                end = float(self.apodization_end_textcontrol_dim2.GetValue())
                power = float(self.apodization_power_textcontrol_dim2.GetValue())
                power = int(power)
            except:
                msg = wx.MessageDialog(self, 'The values entered are not valid', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_offset_textcontrol_dim2.SetValue(str(self.offset_dim2))
                self.apodization_end_textcontrol_dim2.SetValue(str(self.end_dim2))
                self.apodization_power_textcontrol_dim2.SetValue(str(self.power_dim2))
                return
            # Check that offset and end are between 0 and 1
            if(offset < 0 or offset > 1):
                msg = wx.MessageDialog(self, 'Offset values must be between 0 and 1', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_offset_textcontrol_dim2.SetValue(str(self.offset_dim2))
                return
            if(end < 0 or end > 1):
                msg = wx.MessageDialog(self, 'End values must be between 0 and 1', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_end_textcontrol_dim2.SetValue(str(self.end_dim2))
                return
            # Check that power is greater than 0
            if(power < 0):
                msg = wx.MessageDialog(self, 'Power must be greater than 0', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_power_textcontrol_dim2.SetValue(str(self.power_dim2))
                return
            self.offset_dim2 = offset
            self.end_dim2 = end
            self.power_dim2 = power
            if(self.nmr_data.pseudo_axis == True and self.nmr_data.index == 1):
                func = np.sin((np.pi*self.offset_dim2 + np.pi*(self.end_dim2-self.offset_dim2)*x)/((((self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2]))))**self.power_dim2
            else:
                func = np.sin((np.pi*self.offset_dim2 + np.pi*(self.end_dim2-self.offset_dim2)*x)/((((self.nmr_data.number_of_points[1]/2)/self.nmr_data.spectral_width[1]))))**self.power_dim2
            self.line1.set_ydata(func)
        elif(self.apodization_dim2_combobox_selection==4):
            try:
                a = float(self.apodization_a_textcontrol_dim2.GetValue())
                b = float(self.apodization_b_textcontrol_dim2.GetValue())
            except:
                msg = wx.MessageDialog(self, 'The values entered are not valid', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_a_textcontrol_dim2.SetValue(str(self.a_dim2))
                self.apodization_b_textcontrol_dim2.SetValue(str(self.b_dim2))
                return
            self.a_dim2 = a
            self.b_dim2 = b
            func = np.exp(-self.a_dim2*(x**2) - self.b_dim2*x)
            self.line1.set_ydata(func)
        elif(self.apodization_dim2_combobox_selection==5):
            try:
                t1 = float(self.apodization_t1_textcontrol_dim2.GetValue())
                t2 = float(self.apodization_t2_textcontrol_dim2.GetValue())
            except:
                msg = wx.MessageDialog(self, 'The values entered are not valid', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_t1_textcontrol_dim2.SetValue(str(self.t1_dim2))
                self.apodization_t2_textcontrol_dim2.SetValue(str(self.t2_dim2))
                return
            # Ensure that t1 and t2 are greater than 0
            if(t1 < 0 or t2 < 0):
                msg = wx.MessageDialog(self, 'Ramp up and ramp down points must be greater than 0', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_t1_textcontrol_dim2.SetValue(str(self.t1_dim2))
                self.apodization_t2_textcontrol_dim2.SetValue(str(self.t2_dim2))
                return
            # Ensure that t1 + t2 is less than the number of points
            if(self.nmr_data.pseudo_axis == True and self.nmr_data.index == 1):
                if(t1 + t2 > self.nmr_data.number_of_points[2]):
                    message = 'Ramp up and ramp down points must be less than the number of points (' + str(self.nmr_data.number_of_points[2]) + ')'
                    msg = wx.MessageDialog(self, message, 'Error', wx.OK | wx.ICON_ERROR)
                    msg.ShowModal()
                    msg.Destroy()
                    self.apodization_t1_textcontrol_dim2.SetValue(str(self.t1_dim2))
                    self.apodization_t2_textcontrol_dim2.SetValue(str(self.t2_dim2))
                    return
            
            else:
                if(t1 + t2 > self.nmr_data.number_of_points[1]):
                    message = 'Ramp up and ramp down points must be less than the number of points (' + str(self.nmr_data.number_of_points[1]) + ')'
                    msg = wx.MessageDialog(self, message, 'Error', wx.OK | wx.ICON_ERROR)
                    msg.ShowModal()
                    msg.Destroy()
                    self.apodization_t1_textcontrol_dim2.SetValue(str(self.t1_dim2))
                    self.apodization_t2_textcontrol_dim2.SetValue(str(self.t2_dim2))
                    return
            self.t1_dim2 = t1
            self.t2_dim2 = t2
            if(self.nmr_data.pseudo_axis == True and self.nmr_data.index == 1):
                func = np.concatenate((np.linspace(0, 1, int(self.t1_dim2)), np.ones(int(self.nmr_data.number_of_points[2]/2) - int(self.t1_dim2) - int(self.t2_dim2)),np.linspace(1, 0, int(self.t2_dim2))))
            else:
                func = np.concatenate((np.linspace(0, 1, int(self.t1_dim2)), np.ones(int(self.nmr_data.number_of_points[1]/2) - int(self.t1_dim2) - int(self.t2_dim2)),np.linspace(1, 0, int(self.t2_dim2))))
            self.line1.set_ydata(func)
        elif(self.apodization_dim2_combobox_selection==6):
            try:
                loc = float(self.apodization_loc_textcontrol_dim2.GetValue())
            except:
                msg = wx.MessageDialog(self, 'The values entered are not valid', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_loc_textcontrol_dim2.SetValue(str(self.loc_dim2))
                return
            # Ensure that loc is between 0 and 1
            if(loc < 0 or loc > 1):
                msg = wx.MessageDialog(self, 'Location of maximum must be between 0 and 1', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_loc_textcontrol_dim2.SetValue(str(self.loc_dim2))
                return
            self.loc_dim2 = loc
            if(self.nmr_data.pseudo_axis == True and self.nmr_data.index == 1):
                func = np.concatenate((np.linspace(0, 1, int(self.loc_dim2*(self.nmr_data.number_of_points[2]/2))), np.linspace(1, 0, int(self.nmr_data.number_of_points[2]/2) - int(self.loc_dim2*self.nmr_data.number_of_points[2]))))
            else:
                func = np.concatenate((np.linspace(0, 1, int(self.loc_dim2*(self.nmr_data.number_of_points[1]/2))), np.linspace(1, 0, int(self.nmr_data.number_of_points[1]/2) - int(self.loc_dim2*self.nmr_data.number_of_points[1]))))
            self.line1.set_ydata(func)

        self.apodization_plot_canvas.draw()


    def on_apodization_combobox_dim2(self,event):
        self.apodization_dim2_combobox_selection= self.apodization_combobox_dim2.GetSelection()

        # Destroy the combobox and textcontrols for the previous apodization function
        # self.apodization_sizer.Detach(self.apodization_combobox)
        # self.apodization_combobox.Destroy()

        # # Remove the zf sizer
        self.zero_filling_sizer_dim2.Clear(delete_windows=True)

        if(self.apodization_dim2_combobox_selection_old == 1):
            # Remove the previous textcontrols
            
            self.apodization_sizer_dim2.Detach(self.apodization_line_broadening_label)
            self.apodization_line_broadening_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_line_broadening_textcontrol_dim2)
            self.apodization_line_broadening_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_textcontrol_dim2)
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)

        elif(self.apodization_dim2_combobox_selection_old == 2):
            self.apodization_sizer_dim2.Detach(self.apodization_g1_label)
            self.apodization_g1_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g1_textcontrol_dim2)
            self.apodization_g1_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g2_label)
            self.apodization_g2_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g2_textcontrol_dim2)
            self.apodization_g2_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g3_label)
            self.apodization_g3_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g3_textcontrol_dim2)
            self.apodization_g3_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_textcontrol_dim2)
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)

        elif(self.apodization_dim2_combobox_selection_old == 3):
            self.apodization_sizer_dim2.Detach(self.apodization_offset_label)
            self.apodization_offset_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_offset_textcontrol_dim2)
            self.apodization_offset_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_end_label)
            self.apodization_end_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_end_textcontrol_dim2)
            self.apodization_end_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_power_label)
            self.apodization_power_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_power_textcontrol_dim2)
            self.apodization_power_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_textcontrol_dim2)
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        elif(self.apodization_dim2_combobox_selection_old == 4):
            self.apodization_sizer_dim2.Detach(self.apodization_a_label)
            self.apodization_a_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_a_textcontrol_dim2)
            self.apodization_a_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_b_label)
            self.apodization_b_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_b_textcontrol_dim2)
            self.apodization_b_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_textcontrol_dim2)
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        elif(self.apodization_dim2_combobox_selection_old == 5):
            self.apodization_sizer_dim2.Detach(self.apodization_t1_label)
            self.apodization_t1_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_t1_textcontrol_dim2)
            self.apodization_t1_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_t2_label)
            self.apodization_t2_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_t2_textcontrol_dim2)
            self.apodization_t2_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_textcontrol_dim2)
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        elif(self.apodization_dim2_combobox_selection_old == 6):
            self.apodization_sizer_dim2.Detach(self.apodization_loc_label)
            self.apodization_loc_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_loc_textcontrol_dim2)
            self.apodization_loc_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_textcontrol_dim2)
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        elif(self.apodization_dim2_combobox_selection_old == 0):
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        


        self.apodization_sizer_dim2.Detach(self.apodization_checkbox_dim2)
        self.apodization_checkbox_dim2.Destroy()
        self.apodization_sizer_dim2.Detach(self.apodization_combobox_dim2)
        self.apodization_combobox_dim2.Hide()
       


        

        # Delete the current apodization sizer and then create a new one
        
        self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
        self.apodization_plot_ax_dim2.clear()
        self.apodization_plot_ax_dim2.clear()
        self.apodization_plot_sizer_dim2.Clear(True)

        self.sizer_2.Remove(self.apodization_sizer_dim2)
        # self.apodization_sizer.Clear(delete_windows=True)

 


        



        # Remove the linear prediction sizers
        self.linear_prediction_sizer_dim2.Clear(delete_windows=True)
        # self.sizer_1.Remove(self.linear_prediction_sizer)
        

        # self.sizer_1.Remove(self.solvent_suppression_sizer)

        
        self.sizer_2.Clear(delete_windows=True)

           
        self.create_menu_bar_dim2()
        self.Refresh()
        self.Update()
        self.Layout()




        self.apodization_dim2_combobox_selection_old = self.apodization_dim2_combobox_selection



    def create_zero_filling_sizer_dim2(self, parent):
        # Create a box for zero filling options
        self.zero_filling_box_dim2 = wx.StaticBox(parent, -1, 'Zero Filling')
        self.zero_filling_sizer_dim2 = wx.StaticBoxSizer(self.zero_filling_box_dim2, wx.HORIZONTAL)
        self.zero_filling_checkbox_dim2 = wx.CheckBox(parent, -1, 'Apply zero filling')
        self.zero_filling_checkbox_dim2.SetValue(self.zero_filling_dim2_checkbox_value)
        self.zero_filling_checkbox_dim2.Bind(wx.EVT_CHECKBOX, self.on_zero_filling_checkbox_dim2)
        self.zero_filling_sizer_dim2.Add(self.zero_filling_checkbox_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.zero_filling_sizer_dim2.AddSpacer(10)
        # Have a combobox for zero filling options
        self.zf_options_label = wx.StaticText(parent, -1, 'Options:')
        self.zero_filling_sizer_dim2.Add(self.zf_options_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.zero_filling_sizer_dim2.AddSpacer(5)
        self.zero_filling_options_dim2 = ['Doubling spectrum size', 'Adding Zeros', 'Final data size']
        self.zero_filling_combobox_dim2 = wx.ComboBox(parent, -1, choices=self.zero_filling_options_dim2, style=wx.CB_READONLY)
        self.zero_filling_combobox_dim2.Bind(wx.EVT_COMBOBOX, self.on_zero_filling_combobox_dim2)
        self.zero_filling_combobox_dim2.SetSelection(self.zero_filling_dim2_combobox_selection)
        self.zero_filling_sizer_dim2.Add(self.zero_filling_combobox_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.zero_filling_sizer_dim2.AddSpacer(10)
        if(self.zero_filling_dim2_combobox_selection == 0):
            # Have a textcontrol for the doubling number/number of zeros/final data size
            self.zf_value_label = wx.StaticText(parent, -1, 'Doubling number:')
            self.zero_filling_sizer_dim2.Add(self.zf_value_label, 0, wx.ALIGN_CENTER_VERTICAL)

            self.zero_filling_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.zero_filling_dim2_value_doubling_times), size=(40, 20))
            self.zero_filling_sizer_dim2.AddSpacer(5)
            self.zero_filling_sizer_dim2.Add(self.zero_filling_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.zero_filling_textcontrol_dim2.Bind(wx.EVT_TEXT, self.on_zero_filling_doubling_number_dim2)
            self.zero_filling_sizer_dim2.AddSpacer(20)
        elif(self.zero_filling_dim2_combobox_selection == 1):
            # Have a textcontrol for the doubling number/number of zeros/final data size
            self.zf_value_label = wx.StaticText(parent, -1, 'Number of zeros to add:')
            self.zero_filling_sizer_dim2.Add(self.zf_value_label, 0, wx.ALIGN_CENTER_VERTICAL)

            self.zero_filling_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.zero_filling_dim2_value_zeros_to_add), size=(40, 20))
            self.zero_filling_sizer_dim2.AddSpacer(5)
            self.zero_filling_sizer_dim2.Add(self.zero_filling_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.zero_filling_textcontrol_dim2.Bind(wx.EVT_TEXT, self.on_zero_filling_zeros_to_add_dim2)    
            self.zero_filling_sizer_dim2.AddSpacer(20)
        elif(self.zero_filling_dim2_combobox_selection == 2):
            # Have a textcontrol for the doubling number/number of zeros/final data size
            self.zf_value_label = wx.StaticText(parent, -1, 'Final data size:')
            self.zero_filling_sizer_dim2.Add(self.zf_value_label, 0, wx.ALIGN_CENTER_VERTICAL)

            self.zero_filling_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.zero_filling_dim2_value_final_data_size), size=(40, 20))
            self.zero_filling_textcontrol_dim2.Bind(wx.EVT_TEXT, self.on_zero_filling_final_size_dim2)
            self.zero_filling_sizer_dim2.AddSpacer(5)
            self.zero_filling_sizer_dim2.Add(self.zero_filling_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
            self.zero_filling_sizer_dim2.AddSpacer(20)

        # Have a checkbox for rounding to the nearest power of 2
        self.zero_filling_round_checkbox_dim2 = wx.CheckBox(parent, -1, 'Round to nearest power of 2')
        self.zero_filling_round_checkbox_dim2.SetValue(self.zero_filling_dim2_round_checkbox_value)
        self.zero_filling_round_checkbox_dim2.Bind(wx.EVT_CHECKBOX, self.on_zero_filling_round_checkbox_dim2)
        self.zero_filling_sizer_dim2.Add(self.zero_filling_round_checkbox_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.zero_filling_sizer_dim2.AddSpacer(10)
        

        # Have a button showing information on zero filling
        self.zero_filling_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))
        self.zero_filling_info.Bind(wx.EVT_BUTTON, self.oneDFrame.on_zero_fill_info)
        self.zero_filling_sizer_dim2.Add(self.zero_filling_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.zero_filling_sizer_dim2.AddSpacer(10)

        self.sizer_2.Add(self.zero_filling_sizer_dim2)
        self.sizer_2.AddSpacer(10)

    def on_zero_filling_checkbox_dim2(self, event):
        self.zero_filling_dim2_checkbox_value = self.zero_filling_checkbox_dim2.GetValue()


    def on_zero_filling_round_checkbox_dim2(self, event):
        self.zero_filling_dim2_round_checkbox_value = self.zero_filling_round_checkbox_dim2.GetValue()


    def on_zero_filling_final_size_dim2(self, event):
        self.zero_filling_dim2_value_final_data_size = self.zero_filling_textcontrol_dim2.GetValue()

    def on_zero_filling_zeros_to_add_dim2(self, event):
        self.zero_filling_dim2_value_zeros_to_add = self.zero_filling_textcontrol_dim2.GetValue()

    def on_zero_filling_doubling_number_dim2(self, event):
        self.zero_filling_dim2_value_doubling_times = self.zero_filling_textcontrol_dim2.GetValue()



    def on_zero_filling_combobox_dim2(self, event):
        self.zero_filling_dim2_combobox_selection = self.zero_filling_combobox_dim2.GetSelection()
        # # # Remove the zf sizer
        self.zero_filling_sizer_dim2.Clear()
        self.zero_filling_sizer_dim2.Detach(self.zero_filling_checkbox_dim2)
        self.zero_filling_checkbox_dim2.Destroy()
        self.zero_filling_sizer_dim2.Detach(self.zf_options_label)
        self.zf_options_label.Destroy()
        self.zero_filling_sizer_dim2.Detach(self.zero_filling_info)
        self.zero_filling_info.Destroy()
        self.zero_filling_sizer_dim2.Detach(self.zf_value_label)
        self.zf_value_label.Destroy()
        self.zero_filling_sizer_dim2.Detach(self.zero_filling_round_checkbox_dim2)
        self.zero_filling_round_checkbox_dim2.Destroy()
        self.zero_filling_sizer_dim2.Detach(self.zero_filling_textcontrol_dim2)
        self.zero_filling_textcontrol_dim2.Destroy()

        self.zero_filling_sizer_dim2.Detach(self.zero_filling_combobox_dim2)
        self.zero_filling_combobox_dim2.Hide()

        if(self.apodization_dim2_combobox_selection_old == 1):
            # Remove the previous textcontrols
            
            self.apodization_sizer_dim2.Detach(self.apodization_line_broadening_label)
            self.apodization_line_broadening_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_line_broadening_textcontrol_dim2)
            self.apodization_line_broadening_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_textcontrol_dim2)
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)

        elif(self.apodization_dim2_combobox_selection_old == 2):
            self.apodization_sizer_dim2.Detach(self.apodization_g1_label)
            self.apodization_g1_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g1_textcontrol_dim2)
            self.apodization_g1_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g2_label)
            self.apodization_g2_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g2_textcontrol_dim2)
            self.apodization_g2_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g3_label)
            self.apodization_g3_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_g3_textcontrol_dim2)
            self.apodization_g3_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_textcontrol_dim2)
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)

        elif(self.apodization_dim2_combobox_selection_old == 3):
            self.apodization_sizer_dim2.Detach(self.apodization_offset_label)
            self.apodization_offset_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_offset_textcontrol_dim2)
            self.apodization_offset_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_end_label)
            self.apodization_end_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_end_textcontrol_dim2)
            self.apodization_end_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_power_label)
            self.apodization_power_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_power_textcontrol_dim2)
            self.apodization_power_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_textcontrol_dim2)
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        elif(self.apodization_dim2_combobox_selection_old == 4):
            self.apodization_sizer_dim2.Detach(self.apodization_a_label)
            self.apodization_a_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_a_textcontrol_dim2)
            self.apodization_a_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_b_label)
            self.apodization_b_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_b_textcontrol_dim2)
            self.apodization_b_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_textcontrol_dim2)
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        elif(self.apodization_dim2_combobox_selection_old == 5):
            self.apodization_sizer_dim2.Detach(self.apodization_t1_label)
            self.apodization_t1_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_t1_textcontrol_dim2)
            self.apodization_t1_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_t2_label)
            self.apodization_t2_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_t2_textcontrol_dim2)
            self.apodization_t2_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_textcontrol_dim2)
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        elif(self.apodization_dim2_combobox_selection_old == 6):
            self.apodization_sizer_dim2.Detach(self.apodization_loc_label)
            self.apodization_loc_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_loc_textcontrol_dim2)
            self.apodization_loc_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_first_point_textcontrol_dim2)
            self.apodization_first_point_textcontrol_dim2.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        elif(self.apodization_dim2_combobox_selection_old == 0):
            self.apodization_sizer_dim2.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
            self.apodization_plot_sizer_dim2.Clear(True)
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_ax_dim2.clear()
            self.apodization_plot_sizer_dim2.Clear(True)
        


        self.apodization_sizer_dim2.Detach(self.apodization_checkbox_dim2)
        self.apodization_checkbox_dim2.Destroy()
        self.apodization_sizer_dim2.Detach(self.apodization_combobox_dim2)
        self.apodization_combobox_dim2.Hide()
       


        

        # Delete the current apodization sizer and then create a new one
        
        self.apodization_sizer_dim2.Detach(self.apodization_plot_sizer_dim2)
        self.apodization_plot_ax_dim2.clear()
        self.apodization_plot_ax_dim2.clear()
        self.apodization_plot_sizer_dim2.Clear(True)

        self.sizer_2.Remove(self.apodization_sizer_dim2)
        # self.apodization_sizer.Clear(delete_windows=True)

 


        



        # Remove the linear prediction sizers
        self.linear_prediction_sizer_dim2.Clear(delete_windows=True)
        # self.sizer_1.Remove(self.linear_prediction_sizer)
        

        # self.sizer_1.Remove(self.solvent_suppression_sizer)

        
        self.sizer_2.Clear(delete_windows=True)

        self.create_menu_bar_dim2()
        self.Refresh()
        self.Update()
        self.Layout()





    def create_fourier_transform_sizer_dim2(self, parent):
        # Create a box for fourier transform options
        self.fourier_transform_box = wx.StaticBox(parent, -1, 'Fourier Transform')
        self.fourier_transform_sizer_dim2 = wx.StaticBoxSizer(self.fourier_transform_box, wx.HORIZONTAL)
        self.fourier_transform_checkbox_dim2 = wx.CheckBox(parent, -1, 'Apply fourier transform')
        self.fourier_transform_checkbox_dim2.Bind(wx.EVT_CHECKBOX, self.on_fourier_transform_checkbox_dim2)
        self.fourier_transform_checkbox_dim2.SetValue(self.fourier_transform_dim2_checkbox_value)
        self.fourier_transform_sizer_dim2.Add(self.fourier_transform_checkbox_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.fourier_transform_sizer_dim2.AddSpacer(10)
        # Have a button for advanced options for fourier transform
        self.fourier_transform_advanced_options_dim2 = wx.Button(parent, -1, 'Advanced Options')
        self.fourier_transform_advanced_options_dim2.Bind(wx.EVT_BUTTON, self.on_fourier_transform_advanced_options_dim2)
        self.fourier_transform_sizer_dim2.Add(self.fourier_transform_advanced_options_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.fourier_transform_sizer_dim2.AddSpacer(10)

        # Have a button showing information on fourier transform
        self.fourier_transform_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))
        self.fourier_transform_info.Bind(wx.EVT_BUTTON, self.oneDFrame.on_fourier_transform_info)
        self.fourier_transform_sizer_dim2.Add(self.fourier_transform_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer_2.Add(self.fourier_transform_sizer_dim2)
        self.sizer_2.AddSpacer(10)

    def on_fourier_transform_checkbox_dim2(self, event):
        self.fourier_transform_dim2_checkbox_value = self.fourier_transform_checkbox_dim2.GetValue()


    def on_fourier_transform_advanced_options_dim2(self, event):
        # Create a frame with a set of advanced options for the fourier transform implementation
        self.fourier_transform_advanced_options_window_dim2 = wx.Frame(self, -1, 'Fourier Transform Advanced Options (Dimension 2)', size=(700, 300))
        if(darkdetect.isDark() == True and platform!='windows'):
            self.fourier_transform_advanced_options_window_dim2.SetBackgroundColour((53, 53, 53, 255))
        else:
            self.fourier_transform_advanced_options_window_dim2.SetBackgroundColour('White')
        
        self.fourier_transform_advanced_options_window_sizer_dim2 = wx.BoxSizer(wx.VERTICAL)
        self.fourier_transform_advanced_options_window_dim2.SetSizer(self.fourier_transform_advanced_options_window_sizer_dim2)

        # Create a sizer for the fourier transform advanced options
        self.ft_label = wx.StaticBox(self.fourier_transform_advanced_options_window_dim2, -1, 'Fourier Transform Method:')
        self.fourier_transform_advanced_options_sizer_dim2 = wx.StaticBoxSizer(self.ft_label,wx.VERTICAL)

        # Have a radiobox for auto, real, inverse, sign alternation
        self.fourier_transform_advanced_options_sizer_dim2.AddSpacer(10)
        self.fourier_transform_auto_real_inverse_sign_alternation_radio_box_dim2 = wx.RadioBox(self.fourier_transform_advanced_options_window_dim2, -1, choices=['Auto', 'Real', 'Inverse', 'Sign Alternation', 'Negative'], style=wx.RA_SPECIFY_COLS)
        self.fourier_transform_auto_real_inverse_sign_alternation_radio_box_dim2.SetSelection(self.ft_method_selection_dim2)
        self.fourier_transform_advanced_options_sizer_dim2.Add(self.fourier_transform_auto_real_inverse_sign_alternation_radio_box_dim2, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.fourier_transform_advanced_options_sizer_dim2.AddSpacer(10)


        self.ft_method_text = 'Auto: The auto method will automatically select the best method for the fourier transform of the FID. \n\nReal: The Fourier Transform will be applied to the real part of the FID only. \n\nInverse: The inverse Fourier Transform will be applied to the FID. \n\nSign Alternation: The sign alternation method will be applied to the FID. \n\n'

        self.ft_method_info = wx.StaticText(self.fourier_transform_advanced_options_window_dim2, -1, self.ft_method_text)
        self.fourier_transform_advanced_options_sizer_dim2.Add(self.ft_method_info, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.fourier_transform_advanced_options_sizer_dim2.AddSpacer(10)

        # Have a save and close button
        self.fourier_transform_advanced_options_save_button_dim2 = wx.Button(self.fourier_transform_advanced_options_window_dim2, -1, 'Save and Close')
        self.fourier_transform_advanced_options_save_button_dim2.Bind(wx.EVT_BUTTON, self.on_fourier_transform_advanced_options_save_dim2)
        self.fourier_transform_advanced_options_sizer_dim2.Add(self.fourier_transform_advanced_options_save_button_dim2, 0, wx.ALIGN_CENTER_HORIZONTAL)



        self.fourier_transform_advanced_options_window_sizer_dim2.Add(self.fourier_transform_advanced_options_sizer_dim2, 0, wx.ALIGN_CENTER_HORIZONTAL)

        self.fourier_transform_advanced_options_window_dim2.Show()

    def on_fourier_transform_advanced_options_save_dim2(self, event):
        # Save the current selection and close the window
        self.ft_method_selection_dim2 = self.fourier_transform_auto_real_inverse_sign_alternation_radio_box_dim2.GetSelection()
        self.fourier_transform_advanced_options_window_dim2.Close()


    def create_phase_correction_sizer_dim2(self, parent):
        # Create a box for phase correction options
        self.phase_correction_box_dim2 = wx.StaticBox(parent, -1, 'Phase Correction')
        self.phase_correction_sizer_dim2 = wx.StaticBoxSizer(self.phase_correction_box_dim2, wx.HORIZONTAL)
        self.phase_correction_checkbox_dim2 = wx.CheckBox(parent, -1, 'Apply phase correction')
        self.phase_correction_checkbox_dim2.Bind(wx.EVT_CHECKBOX, self.on_phase_correction_checkbox_dim2)
        self.phase_correction_checkbox_dim2.SetValue(self.phasing_dim2_checkbox_value)
        self.phase_correction_sizer_dim2.Add(self.phase_correction_checkbox_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.phase_correction_sizer_dim2.AddSpacer(10)
        # Have a textcontrol for p0 and p1 values
        self.phase_correction_p0_label = wx.StaticText(parent, -1, 'Zero order correction (p0):')
        self.phase_correction_sizer_dim2.Add(self.phase_correction_p0_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.phase_correction_p0_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.p0_total_dim2), size=(50, 20))
        self.phase_correction_p0_textcontrol_dim2.Bind(wx.EVT_TEXT, self.on_phase_correction_p0_dim2)
        self.phase_correction_sizer_dim2.Add(self.phase_correction_p0_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.phase_correction_sizer_dim2.AddSpacer(10)
        self.phase_correction_p1_label = wx.StaticText(parent, -1, 'First order correction (p1):')
        self.phase_correction_sizer_dim2.Add(self.phase_correction_p1_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.phase_correction_p1_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.p1_total_dim2), size=(50, 20))
        self.phase_correction_p1_textcontrol_dim2.Bind(wx.EVT_TEXT, self.on_phase_correction_p1_dim2)
        self.phase_correction_sizer_dim2.Add(self.phase_correction_p1_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.phase_correction_sizer_dim2.AddSpacer(10)

        # Have a checkbox for f1180
        self.phase_correction_f1180_button_dim2 = wx.CheckBox(parent, -1, 'F1180')
        self.phase_correction_f1180_button_dim2.Bind(wx.EVT_CHECKBOX, self.on_phase_correction_f1180_dim2)
        self.phase_correction_f1180_button_dim2.SetValue(self.f1180_dim2)
        self.phase_correction_sizer_dim2.Add(self.phase_correction_f1180_button_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.phase_correction_sizer_dim2.AddSpacer(10)

        # Have a button showing information on phase correction
        self.phase_correction_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))
        self.phase_correction_info.Bind(wx.EVT_BUTTON, self.on_phase_correction_info_dim2)
        self.phase_correction_sizer_dim2.Add(self.phase_correction_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer_2.Add(self.phase_correction_sizer_dim2)
        self.sizer_2.AddSpacer(10)


    def on_phase_correction_checkbox_dim2(self, event):
        self.phasing_dim2_checkbox_value = self.phase_correction_checkbox_dim2.GetValue()



    def on_phase_correction_p0_dim2(self, event):
        self.p0_total_dim2 = self.phase_correction_p0_textcontrol_dim2.GetValue()

    def on_phase_correction_p1_dim2(self, event):
        self.p1_total_dim2 = self.phase_correction_p1_textcontrol_dim2.GetValue()
        try:
            if(np.abs(float(self.p1_total_dim2)) > 45):
                self.apodization_first_point_scaling_dim2 = 1.0
                self.apodization_first_point_textcontrol_dim2.SetValue(str(self.apodization_first_point_scaling_dim2))
            else:
                self.apodization_first_point_scaling_dim2 = 0.5
                self.apodization_first_point_textcontrol_dim2.SetValue(str(self.apodization_first_point_scaling_dim2))
        except:
            pass


    def on_phase_correction_f1180_dim2(self, event):
        if(self.phase_correction_f1180_button_dim2.GetValue() == True):
            self.c_old = self.apodization_first_point_scaling_dim2
            self.p0_total_dim2_old = self.p0_total_dim2
            self.p1_total_dim2_old = self.p1_total_dim2
            # Apply -90 p0 and 180 p1 to the phase correction textcontrols
            self.p0_total_dim2 = -90.0
            self.p1_total_dim2 = 180.0
            self.phase_correction_p0_textcontrol_dim2.SetValue(str(self.p0_total_dim2))
            self.phase_correction_p1_textcontrol_dim2.SetValue(str(self.p1_total_dim2))
            # Disable the phase correction textcontrols
            self.phase_correction_p0_textcontrol_dim2.Disable()
            self.phase_correction_p1_textcontrol_dim2.Disable()
            self.apodization_first_point_scaling_dim2 = 1.0
            self.apodization_first_point_textcontrol_dim2.SetValue(str(self.apodization_first_point_scaling_dim2))
        else:
            self.p0_total_dim2 = self.p0_total_dim2_old
            self.p1_total_dim2 = self.p1_total_dim2_old
            self.phase_correction_p0_textcontrol_dim2.SetValue(str(self.p0_total_dim2))
            self.phase_correction_p1_textcontrol_dim2.SetValue(str(self.p1_total_dim2))
            self.phase_correction_p0_textcontrol_dim2.Enable()
            self.phase_correction_p1_textcontrol_dim2.Enable()
            try:
                self.apodization_first_point_scaling_dim2 = self.c_old
            except:
                self.apodization_first_point_scaling_dim2 = 0.5
            self.apodization_first_point_textcontrol_dim2.SetValue(str(self.apodization_first_point_scaling_dim2))


    def on_phase_correction_info_dim2(self, event):
        phase_correction_text = 'Phase correction is a method to correct for phase errors in the FID. Zero order phase correction (p0) is used to correct a phase offset that is applied equally across the spectrum. However, a first order phase correction (p1) is used to correct the phasing in a spectrum where peaks in different locations of the spectrum require a different phasing value. For the indirect dimension, it is often the case that the acquisition is delayed by an exact time so that the resulting spectrum can be phased using the phase values of: p0=-90, p1=180. This is often termed F1180. \n Further information can be found using the link below.'

        # Create a popup window with the information
        self.phase_correction_info_window = wx.Frame(self, -1, 'Phase Correction Information', size=(450, 300))
        if(darkdetect.isDark() == True and platform!='windows'):
            self.phase_correction_info_window.SetBackgroundColour((53, 53, 53, 255))
            colour = "RED"
        else:
            self.phase_correction_info_window.SetBackgroundColour('White')
            colour = "BLUE"
        
        self.phase_correction_info_window_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.phase_correction_info_window.SetSizer(self.phase_correction_info_window_sizer)

        self.phase_correction_info_window_sizer.AddSpacer(10)

        # Create a sizer for the phase correction information
        self.phase_correction_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.phase_correction_info_sizer.AddSpacer(10)
        self.phase_correction_info_sizer.Add(wx.StaticText(self.phase_correction_info_window, -1, phase_correction_text, size=(400,200)), 0, wx.ALIGN_CENTER)
        self.phase_correction_info_sizer.AddSpacer(10)

        # Have a hyperlink to the phase correction information
        self.phase_correction_info_hyperlink = hl.HyperLinkCtrl(self.phase_correction_info_window, -1, 'NMRPipe Help Page for Phase Correction', URL='http://www.nmrscience.com/ref/nmrpipe/ps.html')
        self.phase_correction_info_hyperlink.SetColours(colour, colour, colour)
        self.phase_correction_info_hyperlink.SetUnderlines(False, False, False)
        self.phase_correction_info_hyperlink.SetBold(False)
        self.phase_correction_info_hyperlink.UpdateLink()
        self.phase_correction_info_sizer.Add(self.phase_correction_info_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.phase_correction_info_sizer.AddSpacer(10)

        self.phase_correction_info_window_sizer.Add(self.phase_correction_info_sizer, 0, wx.ALIGN_CENTER)

        self.phase_correction_info_window.Show()


    def create_extraction_sizer_dim2(self, parent):
        # A box for extraction of data between two ppm values
        self.extraction_box_dim2 = wx.StaticBox(parent, -1, 'Extraction')
        self.extraction_sizer_dim2 = wx.StaticBoxSizer(self.extraction_box_dim2, wx.HORIZONTAL)
        self.extraction_checkbox_dim2 = wx.CheckBox(parent, -1, 'Include data extraction')
        self.extraction_checkbox_dim2.Bind(wx.EVT_CHECKBOX, self.on_extraction_checkbox_dim2)
        self.extraction_checkbox_dim2.SetValue(self.extraction_checkbox_value_dim2)
        self.extraction_sizer_dim2.Add(self.extraction_checkbox_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.extraction_sizer_dim2.AddSpacer(10)
        # Have a textcontrol for the ppm start value
        self.extraction_ppm_start_label = wx.StaticText(parent, -1, 'Start chemical shift (ppm):')
        self.extraction_sizer_dim2.Add(self.extraction_ppm_start_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.extraction_ppm_start_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.extraction_start_dim2), size=(40, 20))
        self.extraction_ppm_start_textcontrol_dim2.Bind(wx.EVT_TEXT, self.on_extraction_dim2)
        self.extraction_sizer_dim2.Add(self.extraction_ppm_start_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.extraction_sizer_dim2.AddSpacer(10)
        # Have a textcontrol for the ppm end value
        self.extraction_ppm_end_label = wx.StaticText(parent, -1, 'End chemical shift (ppm):')
        self.extraction_sizer_dim2.Add(self.extraction_ppm_end_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.extraction_ppm_end_textcontrol_dim2 = wx.TextCtrl(parent, -1, str(self.extraction_end_dim2), size=(40, 20))
        self.extraction_ppm_end_textcontrol_dim2.Bind(wx.EVT_TEXT, self.on_extraction_dim2)
        self.extraction_sizer_dim2.Add(self.extraction_ppm_end_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.extraction_sizer_dim2.AddSpacer(10)
        # Have a button showing information on extraction
        self.extraction_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))
        self.extraction_info.Bind(wx.EVT_BUTTON, self.oneDFrame.on_extraction_info)
        self.extraction_sizer_dim2.Add(self.extraction_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer_2.Add(self.extraction_sizer_dim2)
        self.sizer_2.AddSpacer(10)

    def on_extraction_checkbox_dim2(self, event):
        self.extraction_checkbox_value_dim2 = self.extraction_checkbox_dim2.GetValue()



    def on_extraction_dim2(self, event):
        self.extraction_start_dim2 = self.extraction_ppm_start_textcontrol_dim2.GetValue()
        self.extraction_end_dim2 = self.extraction_ppm_end_textcontrol_dim2.GetValue()




    def create_baseline_correction_sizer_dim2(self, parent):
        # Create a box for baseline correction options (linear/polynomial)
        self.baseline_correction_box_dim2 = wx.StaticBox(parent, -1, 'Baseline Correction')
        self.baseline_correction_sizer_dim2 = wx.StaticBoxSizer(self.baseline_correction_box_dim2, wx.HORIZONTAL)
        self.baseline_correction_checkbox_dim2 = wx.CheckBox(parent, -1, 'Apply baseline correction')
        self.baseline_correction_checkbox_dim2.Bind(wx.EVT_CHECKBOX, self.on_baseline_correction_checkbox_dim2)
        self.baseline_correction_checkbox_dim2.SetValue(self.baseline_correction_checkbox_value_dim2)
        self.baseline_correction_sizer_dim2.Add(self.baseline_correction_checkbox_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_sizer_dim2.AddSpacer(10)
        # Have a radio box for linear or polynomial baseline correction
        self.baseline_correction_radio_box_dim2 = wx.RadioBox(parent, -1, 'Baseline Correction Method', choices=['Linear', 'Polynomial'])
        # Bind the radio box to a function that will update the baseline correction options
        self.baseline_correction_radio_box_dim2.Bind(wx.EVT_RADIOBOX, self.on_baseline_correction_radio_box_dim2)
        self.baseline_correction_radio_box_dim2.SetSelection(self.baseline_correction_radio_box_selection_dim2)
        self.baseline_correction_sizer_dim2.Add(self.baseline_correction_radio_box_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_sizer_dim2.AddSpacer(10)
        
        # If linear baseline correction is selected, have a textcontrol for the node values to use
        self.baseline_correction_nodes_label = wx.StaticText(parent, -1, 'Node width (pts):')
        self.baseline_correction_sizer_dim2.Add(self.baseline_correction_nodes_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_nodes_textcontrol_dim2 = wx.TextCtrl(parent, -1, self.node_width_dim2, size=(30, 20))
        self.baseline_correction_nodes_textcontrol_dim2.Bind(wx.EVT_TEXT, self.on_baseline_correction_textcontrol_dim2)
        self.baseline_correction_sizer_dim2.Add(self.baseline_correction_nodes_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_sizer_dim2.AddSpacer(10)
        # Have a textcontrol for the node list (percentages)
        self.baseline_correction_node_list_label = wx.StaticText(parent, -1, 'Node list (%):')
        self.baseline_correction_sizer_dim2.Add(self.baseline_correction_node_list_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_node_list_textcontrol_dim2 = wx.TextCtrl(parent, -1, self.node_list_dim2, size=(100, 20))
        self.baseline_correction_node_list_textcontrol_dim2.Bind(wx.EVT_TEXT, self.on_baseline_correction_textcontrol_dim2)
        self.baseline_correction_sizer_dim2.Add(self.baseline_correction_node_list_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_sizer_dim2.AddSpacer(10)
        # If polynomial baseline correction is selected, have a textcontrol for the polynomial order

        self.baseline_correction_polynomial_order_label = wx.StaticText(parent, -1, 'Polynomial order:')
        self.baseline_correction_sizer_dim2.Add(self.baseline_correction_polynomial_order_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_polynomial_order_textcontrol_dim2 = wx.TextCtrl(parent, -1, self.polynomial_order_dim2, size=(30, 20))
        self.baseline_correction_polynomial_order_textcontrol_dim2.Bind(wx.EVT_TEXT, self.on_baseline_correction_textcontrol_dim2)
        self.baseline_correction_sizer_dim2.Add(self.baseline_correction_polynomial_order_textcontrol_dim2, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_sizer_dim2.AddSpacer(10)

        if(self.baseline_correction_radio_box_selection_dim2 == 0):
            self.baseline_correction_polynomial_order_label.Hide()
            self.baseline_correction_polynomial_order_textcontrol_dim2.Hide()


        # Have a button showing information on baseline correction
        self.baseline_correction_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))


        self.baseline_correction_info.Bind(wx.EVT_BUTTON, self.oneDFrame.on_baseline_correction_info)
        self.baseline_correction_sizer_dim2.Add(self.baseline_correction_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer_2.Add(self.baseline_correction_sizer_dim2)
        self.sizer_2.AddSpacer(10)


    def on_baseline_correction_checkbox_dim2(self, event):
        self.baseline_correction_checkbox_value_dim2 = self.baseline_correction_checkbox_dim2.GetValue()



    def on_baseline_correction_radio_box_dim2(self, event):
        # If the user selects linear or polynomial baseline correction, update the options
        self.baseline_correction_radio_box_selection_dim2 = self.baseline_correction_radio_box_dim2.GetSelection()

        if(self.baseline_correction_radio_box_selection_dim2 == 0):
            # Remove the polynomial order textcontrol
            self.baseline_correction_sizer_dim2.Hide(self.baseline_correction_polynomial_order_label)
            self.baseline_correction_sizer_dim2.Hide(self.baseline_correction_polynomial_order_textcontrol_dim2)
            self.baseline_correction_sizer_dim2.Layout()
        elif(self.baseline_correction_radio_box_selection_dim2 == 1):
            # Add the polynomial order textcontrol
            self.baseline_correction_sizer_dim2.Show(self.baseline_correction_polynomial_order_label)
            self.baseline_correction_sizer_dim2.Show(self.baseline_correction_polynomial_order_textcontrol_dim2)
            self.baseline_correction_sizer_dim2.Layout()

    def on_baseline_correction_textcontrol_dim2(self, event):
        # If the node width or node list textcontrols are changed, update the node width and node list
        self.node_width_dim2 = self.baseline_correction_nodes_textcontrol_dim2.GetValue()
        self.node_list_dim2 = self.baseline_correction_node_list_textcontrol_dim2.GetValue()
        self.polynomial_order_dim2 = self.baseline_correction_polynomial_order_textcontrol_dim2.GetValue()






class ThreeDFrame(wx.Panel):
    def __init__(self,parent, oneDFrame):
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 0.7*self.monitorWidth
        self.height = 0.75*self.monitorHeight
        self.parent = parent
        wx.Panel.__init__(self, parent, id=wx.ID_ANY, size=(self.width, self.height))

        self.oneDFrame = oneDFrame
        # Create panel for processing dimension 1 of the data
        self.nmr_data = parent.nmr_data
        self.set_variables_dim3()
        self.create_canvas_dim3()
        self.create_menu_bar_dim3()


    def set_variables_dim3(self):
        
        self.set_initial_linear_prediction_variables_dim3()
        self.set_initial_apodization_variables_dim3()
        self.set_initial_zero_filling_variables_dim3()
        self.set_initial_fourier_transform_variables_dim3()
        self.set_initial_phasing_variables_dim3()
        self.set_initial_extraction_variables_dim3()
        self.set_initial_baseline_correction_variables_dim3()

        if(self.parent.load_variables == True):
            try:
                self.load_variables_from_nmrproc_com_3D()
            except:
                pass

    def load_variables_from_nmrproc_com_3D(self):
        # Open processing_parameters.txt file and load the variables from it
        file = open('processing_parameters.txt', 'r')
        lines = file.readlines()
        file.close()

        include_line = False
        for line in lines:
            if('Dimension 3' in line):
                include_line = True
                continue
            if(include_line == False):
                continue
            if(include_line == True and 'Dimension 4' in line):
                include_line = False
                break
            if(include_line == True):
                line = line.split('\n')[0]
                if(line.split(':')[0] == 'Linear Prediction'):
                    self.linear_prediction_radio_box_dim3_selection = int(line.split(': ')[1])
                if(line.split(':')[0] =='Linear Prediction Options Selection'):
                    self.linear_prediction_dim3_options_selection = int(line.split(': ')[1])
                if(line.split(':')[0] =='Linear Prediction Coefficients Selection'):
                    self.linear_prediction_dim3_coefficients_selection = int(line.split(': ')[1])
                if(line.split(':')[0] =='NUS file'):
                    self.nuslist_name_dim3 = line.split(': ')[1]
                if(line.split(':')[0] =='NUS CPU'):
                    self.smile_nus_cpu_textcontrol_dim3 = int(line.split(': ')[1])
                if(line.split(':')[0] =='NUS Iterations'):
                    self.smile_nus_iterations_textcontrol_dim3 = int(line.split(': ')[1])
                if(line.split(':')[0] =='Apodization'):
                    if('True' in line):
                        self.apodization_dim3_checkbox_value = True
                    else:
                        self.apodization_dim3_checkbox_value = False
                if(line.split(':')[0] =='Apodization Combobox Selection'):
                    self.apodization_dim3_combobox_selection = int(line.split(': ')[1])
                    self.apodization_dim3_combobox_selection_old = int(line.split(': ')[1])
                if(line.split(':')[0] =='Exponential Line Broadening'):
                    self.exponential_line_broadening_dim3 = float(line.split(': ')[1])
                if(line.split(':')[0] =='Apodization First Point Scaling'):
                    self.apodization_first_point_scaling_dim3 = float(line.split(': ')[1])
                if(line.split(':')[0] =='G1'):
                    self.g1_dim3 = float(line.split(': ')[1])
                if(line.split(':')[0] =='G2'):
                    self.g2_dim3 = float(line.split(': ')[1])
                if(line.split(':')[0] =='G3'):
                    self.g3_dim3 = float(line.split(': ')[1])
                if(line.split(':')[0] =='Offset'):
                    self.offset_dim3 = float(line.split(': ')[1])
                if(line.split(':')[0] =='End'):
                    self.end_dim3 = float(line.split(': ')[1])
                if(line.split(':')[0] =='Power'):
                    self.power_dim3 = float(line.split(': ')[1])
                if(line.split(':')[0] =='A'):
                    self.a_dim3 = float(line.split(': ')[1])
                if(line.split(':')[0] =='B'):
                    self.b_dim3 = float(line.split(': ')[1])
                if(line.split(':')[0] =='T1'):
                    self.t1_dim3 = float(line.split(': ')[1])
                if(line.split(':')[0] =='T2'):
                    self.t2_dim3 = float(line.split(': ')[1])
                if(line.split(':')[0] =='Loc'):
                    self.loc_dim3 = float(line.split(': ')[1])

                if(line.split(':')[0] =='Zero Filling'):
                    if('True' in line):
                        self.zero_filling_checkbox_dim3_value = True
                    else:
                        self.zero_filling_checkbox_dim3_value = False
                if(line.split(':')[0] =='Zero Filling Combobox Selection'):
                    self.zero_filling_dim3_combobox_selection = int(line.split(': ')[1])
                if(line.split(':')[0] =='Zero Filling Value Doubling Times'):
                    self.zero_filling_dim3_value_doubling_times = int(line.split(': ')[1])
                if(line.split(':')[0] =='Zero Filling Value Zeros to Add'):
                    self.zero_filling_dim3_value_zeros_to_add = int(line.split(': ')[1])
                if(line.split(':')[0] =='Zero Filling Value Final Data Size'):
                    self.zero_filling_dim3_value_final_data_size = int(line.split(': ')[1])
                if(line.split(':')[0] =='Zero Filling Round Checkbox'):
                    if('True' in line):
                        self.zero_filling_round_checkbox_dim3_value = True
                    else:
                        self.zero_filling_round_checkbox_dim3_value = False
                if(line.split(':')[0] =='Fourier Transform'):
                    if('True' in line):
                        self.fourier_transform_checkbox_dim3_value = True
                    else:
                        self.fourier_transform_checkbox_dim3_value = False
                if(line.split(':')[0] =='Fourier Transform Method Selection'):
                    self.ft_method_selection_dim3 = int(line.split(': ')[1])
                if(line.split(':')[0] =='Phase Correction'):
                    if('True' in line):
                        self.phase_correction_checkbox_dim3_value = True
                    else:
                        self.phase_correction_checkbox_dim3_value = False
                if(line.split(':')[0] =='Phase Correction P0'):
                    self.p0_total_dim3 = float(line.split(': ')[1])
                if(line.split(':')[0] =='Phase Correction P1'):
                    self.p1_total_dim3 = float(line.split(': ')[1])
                if(line.split(':')[0] =='F1180'):
                    if('True' in line):
                        self.f1180_dim3 = True
                    else:
                        self.f1180_dim3 = False
                if(line.split(':')[0] =='Extraction'):
                    if('True' in line):
                        self.extraction_checkbox_dim3_value = True
                    else:
                        self.extraction_checkbox_dim3_value = False
                if(line.split(':')[0] =='Extraction PPM Start'):
                    self.extraction_ppm_start_dim3 = float(line.split(': ')[1])
                if(line.split(':')[0] =='Extraction PPM End'):
                    self.extraction_ppm_end_dim3 = float(line.split(': ')[1])
                if(line.split(':')[0] =='Baseline Correction'):
                    if('True' in line):
                        self.baseline_correction_checkbox_dim3_value = True
                    else:
                        self.baseline_correction_checkbox_dim3_value = False
                if(line.split(':')[0] =='Baseline Correction Radio Box Selection'):
                    self.baseline_correction_radio_box_selection_dim3 = int(line.split(': ')[1])
                if(line.split(':')[0] =='Baseline Correction Nodes'):
                    self.baseline_correction_nodes_dim3 = int(line.split(': ')[1])
                if(line.split(':')[0] =='Baseline Correction Node List'):
                    self.baseline_correction_node_list_dim3 = line.split(': ')[1]
                if(line.split(':')[0] =='Baseline Correction Polynomial Order'):
                    self.baseline_correction_polynomial_order_dim3 = int(line.split(': ')[1])

        
    
    def set_initial_linear_prediction_variables_dim3(self):
        self.linear_prediction_radio_box_dim3_selection = 1
        self.linear_prediction_dim3_checkbox_value = False
        self.linear_prediction_dim3_options_selection = 0
        self.linear_prediction_dim3_coefficients_selection = 0
        self.linear_prediction_selection = 0

        # Check to see if the nuslist file exists in the current directory using os.path.isfile('nuslist')
        if(os.path.isfile('nuslist')):
            self.nuslist_name_dim3 = 'nuslist'
        else:
            self.nuslist_name_dim3 = ''

        self.number_of_nus_CPU_dim3 = 1
        self.nus_iterations_dim3 = 800
        self.nus_data_extension_dim3 = 0 #int(self.nmr_data.number_of_points[2]*1.5)




    def set_initial_apodization_variables_dim3(self):
        self.apodization_dim3_checkbox_value = True
        self.apodization_dim3_combobox_selection = 1
        self.apodization_dim3_combobox_selection_old = 1
        
        # Initial values for exponential apodization
        self.exponential_line_broadening_dim3 = 0.5
        self.apodization_first_point_scaling_dim3 = 0.5

        # Initial values for Lorentz to Gauss apodization
        self.g1_dim3 = 0.33
        self.g2_dim3 = 1
        self.g3_dim3 = 0.0

        # Initial values for Sinebell apodization
        self.offset_dim3 = 0.5
        self.end_dim3 = 0.98
        self.power_dim3 = 1.0

        # Initial values for Gauss Broadening apodization
        self.a_dim3 = 1.0
        self.b_dim3 = 1.0

        # Initial values for Trapezoid apodization
        self.t1_dim3 = int((self.nmr_data.number_of_points[2]/2)/4)
        self.t2_dim3 = int((self.nmr_data.number_of_points[2]/2)/4)

        # Initial values for Triangle apodization
        self.loc_dim3 = 0.5

    def set_initial_zero_filling_variables_dim3(self):
        self.zero_filling_dim3_checkbox_value = True
        self.zero_filling_dim3_combobox_selection = 0
        self.zero_filling_dim3_combobox_selection_old = 0
        self.zero_filling_dim3_value_doubling_times = 1
        self.zero_filling_dim3_value_zeros_to_add = 0
        self.zero_filling_dim3_value_final_data_size = 0
        self.zero_filling_dim3_round_checkbox_value = True

    def set_initial_fourier_transform_variables_dim3(self):
        self.fourier_transform_dim3_checkbox_value = True
        self.ft_method_selection_dim3 = 0 # Initially use the 'auto' method of FT as default


    def set_initial_phasing_variables_dim3(self):
        self.phasing_dim3_checkbox_value = True
        self.p0_total_dim3 = 0.0
        self.p1_total_dim3 = 0.0
        self.p0_total_dim3_old = 0.0
        self.p1_total_dim3_old = 0.0
        self.phasing_from_smile = False
        self.f1180_dim3 = False


    def set_initial_extraction_variables_dim3(self):
        self.extraction_checkbox_value_dim3 = False
        self.extraction_start_dim3 = '0.0'
        self.extraction_end_dim3 = '0.0'




    def set_initial_baseline_correction_variables_dim3(self):
        self.baseline_correction_checkbox_value_dim3 = False
        self.baseline_correction_radio_box_selection_dim3 = 0
        self.node_list_dim3 = '0,5,95,100'
        self.node_width_dim3 = '2'
        self.polynomial_order_dim3 = '4'



        

    def create_canvas_dim3(self):

        if(darkdetect.isDark() == True and platform!='windows'):
            self.SetBackgroundColour((53, 53, 53, 255))
            matplotlib.rc('axes',edgecolor='white')
            matplotlib.rc('xtick',color='white')
            matplotlib.rc('ytick',color='white')
            matplotlib.rc('axes',labelcolor='white')
            matplotlib.rc('axes',facecolor='#282A36')
            matplotlib.rc('legend', labelcolor='white')


        else:
            self.SetBackgroundColour('White')
            matplotlib.rc('axes',edgecolor='black')
            matplotlib.rc('xtick',color='black')
            matplotlib.rc('ytick',color='black')
            matplotlib.rc('axes',labelcolor='black')
            matplotlib.rc('axes',facecolor='white')
            matplotlib.rc('legend', labelcolor='black')


        

        


    def create_menu_bar_dim3(self):
        # Create the main sizer
        self.main_sizer=wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_sizer)
        
        # Create a sizer for the processing options for the first dimension
        self.sizer_2 = wx.BoxSizer(wx.VERTICAL)
        self.sizer_2.AddSpacer(10)
    
        # Create all the sizers (allow a checkbox at the top for SMILE NUS reconstruction which will change the possible options)
        # For NUS reconstruction using SMILE need to have exact phasing paramaters (first process without NUS and then calculate phase in indirect dimension, then process again using SMILE containing exact phasing parameters)
        self.create_linear_prediction_sizer_dim3(parent=self)
        self.create_apodization_sizer_dim3(parent=self)
        self.create_zero_filling_sizer_dim3(parent=self)
        self.create_fourier_transform_sizer_dim3(parent=self)
        self.create_phase_correction_sizer_dim3(parent=self)
        self.create_extraction_sizer_dim3(parent=self)
        self.create_baseline_correction_sizer_dim3(parent=self)
        
        self.main_sizer.Add(self.sizer_2, 0, wx.EXPAND)


        self.SetSizerAndFit(self.main_sizer)
        self.Layout()


        # Get the size of the main sizer and set the window size to 1.05 times the size of the main sizer
        self.width, self.height = self.main_sizer.GetSize()
        self.parent.parent.change_frame_size(int(self.width*1.05), int(self.height*1.25))


    def create_linear_prediction_sizer_dim3(self, parent):
        # Create a sizer for the linear prediction options
        self.linear_prediction_sizer_dim3_label = wx.StaticBox(self, -1, 'Linear Prediction/SMILE NUS Reconstruction')
        self.linear_prediction_sizer_dim3 = wx.StaticBoxSizer(self.linear_prediction_sizer_dim3_label, wx.HORIZONTAL)
        self.linear_prediction_sizer_dim3.AddSpacer(10)

        # Have a radiobox for None, Linear Prediction and SMILE NUS Reconstruction
        self.linear_prediction_radio_box_dim3 = wx.RadioBox(parent, -1, '', choices=['None', 'Linear Prediction', 'SMILE NUS Reconstruction'], style=wx.RA_SPECIFY_ROWS)
        self.linear_prediction_radio_box_dim3.Bind(wx.EVT_RADIOBOX, self.on_linear_prediction_radio_box_dim3)
        self.linear_prediction_radio_box_dim3.SetSelection(self.linear_prediction_radio_box_dim3_selection)

        self.linear_prediction_sizer_dim3.Add(self.linear_prediction_radio_box_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.linear_prediction_sizer_dim3.AddSpacer(10)
        

        if(self.linear_prediction_radio_box_dim3.GetSelection() == 1):
            # Have a combobox for linear prediction options
            self.linear_prediction_options_text = wx.StaticText(parent, -1, 'Add Predicted Points:')
            self.linear_prediction_sizer_dim3.Add(self.linear_prediction_options_text, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim3.AddSpacer(5)
            self.linear_prediction_options = ['After FID', 'Before FID']
            self.linear_prediction_combobox_dim3 = wx.ComboBox(parent, -1, choices=self.linear_prediction_options, style=wx.CB_READONLY)
            self.linear_prediction_combobox_dim3.SetSelection(self.linear_prediction_dim3_options_selection)
            self.linear_prediction_combobox_dim3.Bind(wx.EVT_COMBOBOX, self.on_linear_prediction_combobox_dim3)
            self.linear_prediction_sizer_dim3.Add(self.linear_prediction_combobox_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim3.AddSpacer(10)
            # Have a combobox of predicted coefficient options
            self.linear_prediction_coefficients_text = wx.StaticText(parent, -1, 'Predicted Coefficients:')
            self.linear_prediction_sizer_dim3.Add(self.linear_prediction_coefficients_text, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim3.AddSpacer(5)
            self.linear_prediction_coefficients_options = ['Forward', 'Backward', 'Both']
            self.linear_prediction_coefficients_combobox_dim3 = wx.ComboBox(parent, -1, choices=self.linear_prediction_coefficients_options, style=wx.CB_READONLY)
            self.linear_prediction_coefficients_combobox_dim3.SetSelection(self.linear_prediction_dim3_coefficients_selection)
            self.linear_prediction_coefficients_combobox_dim3.Bind(wx.EVT_COMBOBOX, self.on_linear_prediction_combobox_coefficients_dim3)
            self.linear_prediction_sizer_dim3.Add(self.linear_prediction_coefficients_combobox_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim3.AddSpacer(10)
        elif(self.linear_prediction_radio_box_dim3.GetSelection()==2):
            # Have a set of options for SMILE NUS processing

            # NUS file
            self.smile_nus_file_text = wx.StaticText(parent, -1, 'NUS File:')
            self.linear_prediction_sizer_dim3.Add(self.smile_nus_file_text, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim3.AddSpacer(5)

            self.smile_nus_file_textcontrol_dim3 = wx.TextCtrl(parent, -1, self.nuslist_name_dim3, size=(100, 20))
            self.smile_nus_file_textcontrol_dim3.Bind(wx.EVT_TEXT, self.on_smile_nus_file_textcontrol_dim3)

            self.linear_prediction_sizer_dim3.Add(self.smile_nus_file_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim3.AddSpacer(10)

            # # Zero order phase correction
            # self.smile_nus_p0_text = wx.StaticText(parent, -1, 'Zero Order Phase Correction (p0):')
            # self.linear_prediction_sizer_dim3.Add(self.smile_nus_p0_text, 0, wx.ALIGN_CENTER_VERTICAL)
            # self.linear_prediction_sizer_dim3.AddSpacer(5)
            # self.smile_nus_p0_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.p0_total_dim3), size=(50, 20))
            # self.smile_nus_p0_textcontrol_dim3.Bind(wx.EVT_TEXT, self.on_smile_nus_p0_textcontrol_dim3)
            # self.linear_prediction_sizer_dim3.Add(self.smile_nus_p0_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            # self.linear_prediction_sizer_dim3.AddSpacer(10)

            # # First order phase correction
            # self.smile_nus_p1_text = wx.StaticText(parent, -1, 'First Order Phase Correction (p1):')
            # self.linear_prediction_sizer_dim3.Add(self.smile_nus_p1_text, 0, wx.ALIGN_CENTER_VERTICAL)
            # self.linear_prediction_sizer_dim3.AddSpacer(5)
            # self.smile_nus_p1_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.p1_total_dim3), size=(50, 20))
            # self.smile_nus_p1_textcontrol_dim3.Bind(wx.EVT_TEXT, self.on_smile_nus_p1_textcontrol_dim3)
            # self.linear_prediction_sizer_dim3.Add(self.smile_nus_p1_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            # self.linear_prediction_sizer_dim3.AddSpacer(10)

            # Have a data extension textcontrol
            self.smile_nus_data_extension_text = wx.StaticText(parent, -1, 'Data Extension:')
            self.linear_prediction_sizer_dim3.Add(self.smile_nus_data_extension_text, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim3.AddSpacer(5)
            self.smile_nus_data_extension_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.nus_data_extension_dim3), size=(50, 20))
            self.smile_nus_data_extension_textcontrol_dim3.Bind(wx.EVT_TEXT, self.on_smile_nus_data_extension_textcontrol_dim3)
            self.linear_prediction_sizer_dim3.Add(self.smile_nus_data_extension_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim3.AddSpacer(10)


            # Number of CPU's
            self.smile_nus_cpu_text = wx.StaticText(parent, -1, 'Number of CPU\'s:')
            self.linear_prediction_sizer_dim3.Add(self.smile_nus_cpu_text, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim3.AddSpacer(5)
            self.smile_nus_cpu_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.number_of_nus_CPU_dim3), size=(30, 20))
            self.smile_nus_cpu_textcontrol_dim3.Bind(wx.EVT_TEXT, self.on_smile_nus_cpu_textcontrol_dim3)
            self.linear_prediction_sizer_dim3.Add(self.smile_nus_cpu_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim3.AddSpacer(10)

            # Number of iterations
            self.smile_nus_iterations_text = wx.StaticText(parent, -1, 'Number of Iterations:')
            self.linear_prediction_sizer_dim3.Add(self.smile_nus_iterations_text, 0, wx.ALIGN_CENTER_VERTICAL)
            self.linear_prediction_sizer_dim3.AddSpacer(5)
            self.smile_nus_iterations_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.nus_iterations_dim3), size=(30, 20))
            self.smile_nus_iterations_textcontrol_dim3.Bind(wx.EVT_TEXT, self.on_smile_nus_iterations_textcontrol_dim3)
            self.linear_prediction_sizer_dim3.Add(self.smile_nus_iterations_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)



        # Have a button showing information on linear prediction
        self.linear_prediction_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))
        self.linear_prediction_info.Bind(wx.EVT_BUTTON, self.on_linear_prediction_info_dim3)
        self.linear_prediction_sizer_dim3.AddSpacer(10)
        self.linear_prediction_sizer_dim3.Add(self.linear_prediction_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer_2.Add(self.linear_prediction_sizer_dim3)
        self.sizer_2.AddSpacer(10)   

    def on_linear_prediction_combobox_dim3(self, event):
        # Get the selection from the combobox and update the linear prediction options
        self.linear_prediction_dim3_options_selection = self.linear_prediction_combobox_dim3.GetSelection()

    def on_linear_prediction_combobox_coefficients_dim3(self, event):
        # Get the selection from the combobox and update the linear prediction options
        self.linear_prediction_dim3_coefficients_selection = self.linear_prediction_coefficients_combobox_dim3.GetSelection()

    def on_smile_nus_file_textcontrol_dim3(self, event):
        # Get the value from the textcontrol
        self.nuslist_name_dim3 = self.smile_nus_file_textcontrol_dim3.GetValue()

    def on_smile_nus_p0_textcontrol_dim3(self, event):
        self.p0_total_dim3 = self.smile_nus_p0_textcontrol_dim3.GetValue()
        self.phasing_from_smile = True
        # Update the phasing values in the phasing section too
        self.phase_correction_p0_textcontrol_dim3.SetValue(str(self.p0_total_dim3))
        self.phasing_from_smile = False

    def on_smile_nus_p1_textcontrol_dim3(self, event):
        self.p1_total_dim3 = self.smile_nus_p1_textcontrol_dim3.GetValue()
        self.phasing_from_smile = True
        # Update the phasing values in the phasing section too
        self.phase_correction_p1_textcontrol_dim3.SetValue(str(self.p1_total_dim3))
        self.phasing_from_smile = False


    def on_smile_nus_data_extension_textcontrol_dim3(self, event):
        if(self.smile_nus_data_extension_textcontrol_dim3.GetValue() != ''):
            try:
                self.nus_data_extension_dim3 = int(self.smile_nus_data_extension_textcontrol_dim3.GetValue())
            except:
                self.nus_data_extension_dim3 = int(self.nmr_data.number_of_points[2]*1.5)
                # Give an error message
                error_message = wx.MessageDialog(self, 'Data extension value needs to be an integer. Resetting to original value.', 'Error', wx.OK | wx.ICON_ERROR)
                error_message.ShowModal()
                error_message.Destroy()
        else:
            self.nus_data_extension_dim3 = ''
        
    def on_smile_nus_cpu_textcontrol_dim3(self, event):
        if(self.smile_nus_cpu_textcontrol_dim3.GetValue() != ''):
            try:
                self.number_of_nus_CPU_dim3 = int(self.smile_nus_cpu_textcontrol_dim3.GetValue())
            except:
                self.number_of_nus_CPU_dim3 = 1
                # Give an error message
                error_message = wx.MessageDialog(self, 'Number of CPU\'s needs to be an integer. Resetting to original value.', 'Error', wx.OK | wx.ICON_ERROR)
                error_message.ShowModal()
                error_message.Destroy()
        else:
            self.number_of_nus_CPU_dim3 = ''

    def on_smile_nus_iterations_textcontrol_dim3(self, event):
        if(self.smile_nus_iterations_textcontrol_dim3.GetValue() != ''):
            try:
                self.nus_iterations_dim3 = int(self.smile_nus_iterations_textcontrol_dim3.GetValue())
            except:
                self.nus_iterations_dim3 = 800
                # Give an error message
                error_message = wx.MessageDialog(self, 'Number of iterations needs to be an integer. Resetting to original value.', 'Error', wx.OK | wx.ICON_ERROR)
                error_message.ShowModal()
                error_message.Destroy()
        else:
            self.nus_iterations_dim3 = ''



    def on_linear_prediction_info_dim3(self,event):
        # Create a popout window with information about linear prediction

        # Create a new frame
        self.linear_prediction_info_frame = wx.Frame(self, -1, 'Linear Prediction / SMILE Information', size=(500, 500))
        if(darkdetect.isDark() == True and platform!='windows'):
            self.linear_prediction_info_frame.SetBackgroundColour('#414141')
            colour = "RED"
        else:
            self.linear_prediction_info_frame.SetBackgroundColour('White')
            colour = "BLUE"

        # Create a sizer to hold the box
        self.linear_prediction_info_sizer_window = wx.BoxSizer(wx.VERTICAL)
        self.linear_prediction_info_sizer_window.AddSpacer(10)

        # Create a sizer to hold the text
        self.linear_prediction_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.linear_prediction_info_sizer.AddSpacer(10)

        # Create a text box with the information
        # Linear prediction information
        linear_prediction_information = 'Linear prediction is a method used to increase the resolution of NMR spectra. It is used to predict the points of truncated FIDs (especially in indirect dimensions) and increase signal resolution.\n\n The linear prediction coefficients can be predicted using the forward FID data, backward data or an average of both directions. Then these can be used to add predicted points either before or after the current FID.\n\n Note that advanced options such as  -pred (number of predicted points) and -ord (number of predicted coefficients) can be implemented by manually added them to the nmrproc.com file.'

        self.linear_prediction_info_text = wx.StaticText(self.linear_prediction_info_frame, -1, linear_prediction_information, size=(450, 200), style=wx.ALIGN_CENTER_HORIZONTAL)
        
        # Add the text to the sizer
        self.linear_prediction_info_sizer.Add(self.linear_prediction_info_text, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.linear_prediction_info_sizer.AddSpacer(10)

        # Add a url to the nmrPipe help page
        url = 'http://www.nmrscience.com/ref/nmrpipe/lp.html'  
        self.linear_prediction_info_url = hl.HyperLinkCtrl(self.linear_prediction_info_frame, -1, 'NMRPipe Help Page for Linear Prediction', URL=url)
        self.linear_prediction_info_url.SetColours(colour, colour, colour)
        self.linear_prediction_info_url.SetUnderlines(False, False, False)
        self.linear_prediction_info_url.UpdateLink()


        # Add url to the sizer
        self.linear_prediction_info_sizer.Add(self.linear_prediction_info_url, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.linear_prediction_info_sizer.AddSpacer(10)

        # Add the sizer to the window sizer
        self.linear_prediction_info_sizer_window.Add(self.linear_prediction_info_sizer, 0, wx.ALIGN_CENTER)
        self.linear_prediction_info_sizer_window.AddSpacer(10)


        # Have text to explain SMILE NUS reconstruction
        smile_nus_text = 'SMILE NUS reconstruction is a method used to reconstruct non-uniformly sampled data. The NUS file is a list of points that have been sampled in the FID.\nThe number of CPU\'s (default=1) is the number of cores that will be used to perform the reconstruction and the number of iterations can be changed to improve the accuracy (default=800).\n Furthermore, in order for accurate SMILE reconstruction, the correct zero (p0) and first (p1) order phase correction values need to be inputted.'

        self.smile_nus_text = wx.StaticText(self.linear_prediction_info_frame, -1, smile_nus_text, size=(450, 200), style=wx.ALIGN_CENTER_HORIZONTAL)
        self.linear_prediction_info_sizer_window.Add(self.smile_nus_text, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.linear_prediction_info_sizer_window.AddSpacer(10)




        # Add the window sizer to the frame
        self.linear_prediction_info_frame.SetSizer(self.linear_prediction_info_sizer_window)

        # Show the frame
        self.linear_prediction_info_frame.Show()

    
    def on_linear_prediction_radio_box_dim3(self, event):
        # Get the selection from the radio box and update the linear prediction options
        self.linear_prediction_radio_box_dim3_selection = self.linear_prediction_radio_box_dim3.GetSelection()
        
        # Remove all the old sizers and replot


        if(self.apodization_dim3_combobox_selection_old == 1):
            # Remove the previous textcontrols
            
            self.apodization_sizer_dim3.Detach(self.apodization_line_broadening_label)
            self.apodization_line_broadening_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_line_broadening_textcontrol_dim3)
            self.apodization_line_broadening_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_textcontrol_dim3)
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)

        elif(self.apodization_dim3_combobox_selection_old == 2):
            self.apodization_sizer_dim3.Detach(self.apodization_g1_label)
            self.apodization_g1_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g1_textcontrol_dim3)
            self.apodization_g1_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g2_label)
            self.apodization_g2_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g2_textcontrol_dim3)
            self.apodization_g2_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g3_label)
            self.apodization_g3_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g3_textcontrol_dim3)
            self.apodization_g3_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_textcontrol_dim3)
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)

        elif(self.apodization_dim3_combobox_selection_old == 3):
            self.apodization_sizer_dim3.Detach(self.apodization_offset_label)
            self.apodization_offset_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_offset_textcontrol)
            self.apodization_offset_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_end_label)
            self.apodization_end_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_end_textcontrol_dim3)
            self.apodization_end_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_power_label)
            self.apodization_power_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_power_textcontrol_dim3)
            self.apodization_power_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_textcontrol_dim3)
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif(self.apodization_dim3_combobox_selection_old == 4):
            self.apodization_sizer_dim3.Detach(self.apodization_a_label)
            self.apodization_a_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_a_textcontrol_dim3)
            self.apodization_a_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_b_label)
            self.apodization_b_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_b_textcontrol_dim3)
            self.apodization_b_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_textcontrol_dim3)
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif(self.apodization_dim3_combobox_selection_old == 5):
            self.apodization_sizer_dim3.Detach(self.apodization_t1_label)
            self.apodization_t1_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_t1_textcontrol_dim3)
            self.apodization_t1_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_t2_label)
            self.apodization_t2_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_t2_textcontrol_dim3)
            self.apodization_t2_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_textcontrol_dim3)
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif(self.apodization_dim3_combobox_selection_old == 6):
            self.apodization_sizer_dim3.Detach(self.apodization_loc_label)
            self.apodization_loc_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_loc_textcontrol_dim3)
            self.apodization_loc_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_textcontrol_dim3)
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif(self.apodization_dim3_combobox_selection_old == 0):
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        


        self.apodization_sizer_dim3.Detach(self.apodization_checkbox_dim3)
        self.apodization_checkbox_dim3.Destroy()
        self.apodization_sizer_dim3.Detach(self.apodization_combobox_dim3)
        self.apodization_combobox_dim3.Hide()
       


        

        # Delete the current apodization sizer and then create a new one
        
        self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
        self.apodization_plot_ax_dim3.clear()
        self.apodization_plot_ax_dim3.clear()
        self.apodization_plot_sizer_dim3.Clear(True)

        self.sizer_2.Remove(self.apodization_sizer_dim3)
        # self.apodization_sizer.Clear(delete_windows=True)

 


        



        # Remove the linear prediction sizers
        self.linear_prediction_sizer_dim3.Clear(delete_windows=True)
        # self.sizer_1.Remove(self.linear_prediction_sizer)
        

        # self.sizer_1.Remove(self.solvent_suppression_sizer)

        
        self.sizer_2.Clear(delete_windows=True)

           
        self.create_menu_bar_dim3()
        self.Refresh()
        self.Update()
        self.Layout()


    def create_apodization_sizer_dim3(self, parent):

        # Create a box for apodization options
        self.apodization_box_dim3 = wx.StaticBox(parent, -1, 'Apodization')
        self.apodization_sizer_dim3 = wx.StaticBoxSizer(self.apodization_box_dim3, wx.HORIZONTAL)
        self.apodization_checkbox_dim3 = wx.CheckBox(parent, -1, 'Apply apodization')
        self.apodization_checkbox_dim3.Bind(wx.EVT_CHECKBOX, self.on_apodization_checkbox_dim3)
        self.apodization_checkbox_dim3.SetValue(self.apodization_dim3_checkbox_value)
        self.apodization_sizer_dim3.Add(self.apodization_checkbox_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.apodization_sizer_dim3.AddSpacer(10)
        # Have a combobox for apodization options
        self.apodization_options_dim3 = ['None', 'Exponential', 'Lorentz to Gauss', 'Sinebell', 'Gauss Broadening', 'Trapazoid', 'Triangle']
        self.apodization_combobox_dim3 = wx.ComboBox(parent, -1, choices=self.apodization_options_dim3, style=wx.CB_READONLY)
        self.apodization_combobox_dim3.SetSelection(self.apodization_dim3_combobox_selection)
        self.apodization_combobox_dim3.Bind(wx.EVT_COMBOBOX, self.on_apodization_combobox_dim3)
        self.apodization_sizer_dim3.Add(self.apodization_combobox_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.apodization_sizer_dim3.AddSpacer(10)
        if(self.apodization_dim3_combobox_selection == 1):
            # Have a textcontrol for the line broadening
            self.apodization_line_broadening_label = wx.StaticText(parent, -1, 'Line Broadening (Hz):')
            self.apodization_sizer_dim3.Add(self.apodization_line_broadening_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_line_broadening_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.exponential_line_broadening_dim3), size=(30, 20))
            self.apodization_line_broadening_textcontrol_dim3.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3)
            self.apodization_sizer_dim3.Add(self.apodization_line_broadening_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(parent, -1, 'First Point Scaling:')
            self.apodization_sizer_dim3.Add(self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_first_point_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.apodization_first_point_scaling_dim3), size=(30, 20))
            self.apodization_sizer_dim3.Add(self.apodization_first_point_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim3.AddSpacer(10)
        elif(self.apodization_dim3_combobox_selection == 2):
            # Have a textcontrol for the g1 value
            self.apodization_g1_label = wx.StaticText(parent, -1, 'Inverse Lorentzian (Hz):')
            self.apodization_sizer_dim3.Add(self.apodization_g1_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_g1_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.g1_dim3), size=(40, 20))
            self.apodization_g1_textcontrol_dim3.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3)
            self.apodization_sizer_dim3.Add(self.apodization_g1_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the g2 value
            self.apodization_g2_label = wx.StaticText(parent, -1, 'Gaussian Broadening (Hz):')
            self.apodization_sizer_dim3.Add(self.apodization_g2_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_g2_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.g2_dim3), size=(40, 20))
            self.apodization_g2_textcontrol_dim3.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3)
            self.apodization_sizer_dim3.Add(self.apodization_g2_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the g3 value
            self.apodization_g3_label = wx.StaticText(parent, -1, 'Gaussian Shift:')
            self.apodization_sizer_dim3.Add(self.apodization_g3_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_g3_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.g3_dim3), size=(40, 20))
            self.apodization_g3_textcontrol_dim3.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3)
            self.apodization_sizer_dim3.Add(self.apodization_g3_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(parent, -1, 'First Point Scaling:')
            self.apodization_sizer_dim3.Add(self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_first_point_textcontrol_dim3 = wx.TextCtrl(self, -1, str(self.apodization_first_point_scaling_dim3), size=(30, 20))
            self.apodization_sizer_dim3.Add(self.apodization_first_point_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim3.AddSpacer(10)
        elif(self.apodization_dim3_combobox_selection == 3):
            # Have a textcontrol for the offset value
            self.apodization_offset_label = wx.StaticText(parent, -1, 'Offset (\u03c0):')
            self.apodization_sizer_dim3.Add(self.apodization_offset_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_offset_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.offset_dim3), size=(40, 20))
            self.apodization_offset_textcontrol_dim3.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3)
            self.apodization_sizer_dim3.Add(self.apodization_offset_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the end value
            self.apodization_end_label = wx.StaticText(parent, -1, 'End (\u03c0):')
            self.apodization_sizer_dim3.Add(self.apodization_end_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_end_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.end_dim3), size=(40, 20))
            self.apodization_end_textcontrol_dim3.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3)
            self.apodization_sizer_dim3.Add(self.apodization_end_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the power value
            self.apodization_power_label = wx.StaticText(parent, -1, 'Power:')
            self.apodization_sizer_dim3.Add(self.apodization_power_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_power_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.power_dim3), size=(30, 20))
            self.apodization_power_textcontrol_dim3.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3)
            self.apodization_sizer_dim3.Add(self.apodization_power_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(parent, -1, 'First Point Scaling:')
            self.apodization_sizer_dim3.Add(self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_first_point_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.apodization_first_point_scaling_dim3), size=(30, 20))
            self.apodization_sizer_dim3.Add(self.apodization_first_point_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim3.AddSpacer(10)
        elif(self.apodization_dim3_combobox_selection == 4):
            # Have a textcontrol for the a value
            self.apodization_a_label = wx.StaticText(parent, -1, 'Line Broadening (Hz):')
            self.apodization_sizer_dim3.Add(self.apodization_a_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_a_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.a_dim3), size=(40, 20))
            self.apodization_a_textcontrol_dim3.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3)
            self.apodization_sizer_dim3.Add(self.apodization_a_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the b value
            self.apodization_b_label = wx.StaticText(parent, -1, 'Gaussian Broadening (Hz):')
            self.apodization_sizer_dim3.Add(self.apodization_b_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_b_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.b_dim3), size=(40, 20))
            self.apodization_b_textcontrol_dim3.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3)
            self.apodization_sizer_dim3.Add(self.apodization_b_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(parent, -1, 'First Point Scaling:')
            self.apodization_sizer_dim3.Add(self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_first_point_textcontrol_dim3 = wx.TextCtrl(self, -1, str(self.apodization_first_point_scaling_dim3), size=(30, 20))
            self.apodization_sizer_dim3.Add(self.apodization_first_point_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim3.AddSpacer(10)
        elif(self.apodization_dim3_combobox_selection == 5):
            # Have a textcontrol for the t1 value
            self.apodization_t1_label = wx.StaticText(parent, -1, 'Ramp up points:')
            self.apodization_sizer_dim3.Add(self.apodization_t1_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_t1_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.t1_dim3), size=(50, 20))
            self.apodization_t1_textcontrol_dim3.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3)
            self.apodization_sizer_dim3.Add(self.apodization_t1_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the t2 value
            self.apodization_t2_label = wx.StaticText(parent, -1, 'Ramp down points:')
            self.apodization_sizer_dim3.Add(self.apodization_t2_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_t2_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.t2_dim3), size=(50, 20))
            self.apodization_t2_textcontrol_dim3.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3)
            self.apodization_sizer_dim3.Add(self.apodization_t2_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(parent, -1, 'First Point Scaling:')
            self.apodization_sizer_dim3.Add(self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_first_point_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.apodization_first_point_scaling_dim3), size=(30, 20))
            self.apodization_sizer_dim3.Add(self.apodization_first_point_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim3.AddSpacer(10)
        elif(self.apodization_dim3_combobox_selection == 6):
            # Have a textcontrol for the loc value
            self.apodization_loc_label = wx.StaticText(parent, -1, 'Location of maximum:')
            self.apodization_sizer_dim3.Add(self.apodization_loc_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_loc_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.loc_dim3), size=(40, 20))
            self.apodization_loc_textcontrol_dim3.Bind(wx.EVT_KEY_DOWN, self.on_apodization_textcontrol_dim3)
            self.apodization_sizer_dim3.Add(self.apodization_loc_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim3.AddSpacer(10)
            # Have a textcontrol for the first point scaling
            self.apodization_first_point_label = wx.StaticText(parent, -1, 'First Point Scaling:')
            self.apodization_sizer_dim3.Add(self.apodization_first_point_label, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_first_point_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.apodization_first_point_scaling_dim3), size=(30, 20))
            self.apodization_sizer_dim3.Add(self.apodization_first_point_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.apodization_sizer_dim3.AddSpacer(10)

            
        # Have a button for information on currently selected apodization containing unicode i in a circle
        self.apodization_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))
        self.apodization_info.Bind(wx.EVT_BUTTON, self.oneDFrame.on_apodization_info)
        self.apodization_sizer_dim3.Add(self.apodization_info, 0, wx.ALIGN_CENTER_VERTICAL)

        # Have a mini plots of the apodization function along with the FID first slice
        self.plot_window_function_dim3()




        self.sizer_2.Add(self.apodization_sizer_dim3)
        self.sizer_2.AddSpacer(10)

    def on_apodization_checkbox_dim3(self, event):
        # Get the selection from the checkbox
        self.apodization_dim3_checkbox_selection = self.apodization_checkbox_dim3.GetValue()



    def on_apodization_textcontrol_dim3(self, event):
        # If the user presses enter, update the plot
        keycode = event.GetKeyCode()
        if(keycode == wx.WXK_RETURN):
            self.update_window_function_plot_dim3()
        event.Skip()


    def plot_window_function_dim3(self):
        self.apodization_plot_sizer_dim3 = wx.BoxSizer(wx.VERTICAL)
        if(darkdetect.isDark() == True and platform!='windows'):
            self.apodization_plot_dim3 = Figure(figsize=(1, 0.5),facecolor='#3b3b3b')
        else:
            self.apodization_plot_dim3 = Figure(figsize=(1, 0.5),facecolor='#e6e6e7')
        self.apodization_plot_ax_dim3 = self.apodization_plot_dim3.add_subplot(111)
        # self.apodization_plot_ax.set_axis_off()
        
        if(darkdetect.isDark() == True and platform!='windows'):
            self.apodization_plot_ax_dim3.set_facecolor('#3b3b3b')
        else:
            self.apodization_plot_ax_dim3.set_facecolor('#e6e6e7')

        self.apodization_plot_ax_dim3.set_xticks([])
        self.apodization_plot_ax_dim3.set_yticks([])

        # If the apodization function is None, make remove the axes of the plot
        if(self.apodization_dim3_combobox_selection == 0):
            self.apodization_plot_ax_dim3.spines['top'].set_visible(False)
            self.apodization_plot_ax_dim3.spines['right'].set_visible(False)
            self.apodization_plot_ax_dim3.spines['bottom'].set_visible(False)
            self.apodization_plot_ax_dim3.spines['left'].set_visible(False)


        x=np.linspace(0, (self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2], int(self.nmr_data.number_of_points[2]/2))
        if(self.apodization_dim3_combobox_selection == 1):
            # Exponential window function
            self.line1, = self.apodization_plot_ax_dim3.plot(x, np.exp(-(np.pi*x*self.exponential_line_broadening_dim3)), color='#1f77b4')
            self.apodization_plot_ax_dim3.set_ylim(0, 1.5)
            self.apodization_plot_ax_dim3.set_xlim(0, (self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2])
            
        elif(self.apodization_dim3_combobox_selection == 2):
            # Lorentz to Gauss window function
            e = np.pi * (self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2] * self.g1_dim3
            g = 0.6 * np.pi * self.g2_dim3 * (self.g3_dim3 * ((self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2] - 1) - x)
            func = np.exp(e - g * g)
            self.line1, = self.apodization_plot_ax_dim3.plot(x, func, color='#1f77b4')
            self.apodization_plot_ax_dim3.set_ylim(0, 1.5)
            self.apodization_plot_ax_dim3.set_xlim(0, (self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2])
        elif(self.apodization_dim3_combobox_selection == 3):
            # Sinebell window function
            func = np.sin((np.pi*self.offset_dim3 + np.pi*(self.end_dim3-self.offset_dim3)*x)/((((self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2]))))**self.power_dim3
            self.line1, = self.apodization_plot_ax_dim3.plot(x, func, color='#1f77b4')
            self.apodization_plot_ax_dim3.set_ylim(0, 1.5)
            self.apodization_plot_ax_dim3.set_xlim(0, (self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2])
        elif(self.apodization_dim3_combobox_selection == 4):
            # Gauss broadening window function
            func = np.exp(-self.a_dim3*(x**2) - self.b_dim3*x)
            self.line1, = self.apodization_plot_ax_dim3.plot(x, func, color='#1f77b4')
            self.apodization_plot_ax_dim3.set_ylim(0, 1.5)
            self.apodization_plot_ax_dim3.set_xlim(0, (self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2])
        elif(self.apodization_dim3_combobox_selection == 5):
            # Trapazoid window function
            func = np.concatenate((np.linspace(0, 1, int(self.t1_dim3)), np.ones(int(self.nmr_data.number_of_points[2]/2) - int(self.t1_dim3) - int(self.t1_dim3)),np.linspace(1, 0, int(self.t2_dim3))))
            self.line1, = self.apodization_plot_ax_dim3.plot(x, func, color='#1f77b4')
            self.apodization_plot_ax_dim3.set_ylim(0, 1.5)
            self.apodization_plot_ax_dim3.set_xlim(0, (self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2])
        elif(self.apodization_dim3_combobox_selection == 6):
            # Triangle window function
            func = np.concatenate((np.linspace(0, 1, int(self.loc_dim3*(self.nmr_data.number_of_points[2]/2))), np.linspace(1, 0, int((1-self.loc_dim3)*(self.nmr_data.number_of_points[2]/2)))))
            self.line1, = self.apodization_plot_ax_dim3.plot(x, func, color='#1f77b4')
            
            self.apodization_plot_ax_dim3.set_ylim(0, 1.5)
            self.apodization_plot_ax_dim3.set_xlim(0, (self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2])



            
        self.apodization_plot_ax_dim3.set_xlim(0, (self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2])

        self.apodization_plot_canvas = FigCanvas(self, -1, self.apodization_plot_dim3)
        self.apodization_plot_sizer_dim3.Add(self.apodization_plot_canvas, 0, wx.EXPAND)


        self.apodization_sizer_dim3.Add(self.apodization_plot_sizer_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.apodization_sizer_dim3.AddSpacer(10)


    def update_window_function_plot_dim3(self):
        x=np.linspace(0, (self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[1], int(self.nmr_data.number_of_points[2]/2))
        try:
            c = float(self.apodization_first_point_textcontrol_dim3.GetValue())
            self.apodization_first_point_scaling_dim3 = c
        except:
            # Give a popout window saying that the values are not valid
            msg = wx.MessageDialog(self, 'The value entered for apodization first point scaling is not valid (use 0.5 or 1.0)', 'Error', wx.OK | wx.ICON_ERROR)
            msg.ShowModal()
            msg.Destroy()
            self.apodization_first_point_textcontrol_dim3.SetValue(str(self.apodization_first_point_scaling_dim3))
            return
        if(c != 0.5 and c != 1.0):
            msg = wx.MessageDialog(self, 'The value entered for apodization first point scaling is not valid (use 0.5 or 1.0)', 'Error', wx.OK | wx.ICON_ERROR)
            msg.ShowModal()
            msg.Destroy()
            self.apodization_first_point_textcontrol_dim3.SetValue(str(self.apodization_first_point_scaling_dim3))
            return
        self.apodization_first_point_scaling = c
        if(self.apodization_dim3_combobox_selection==1):
            try:
                em = float(self.apodization_line_broadening_textcontrol_dim3.GetValue())
            except:
                # Give a popout window saying that the values are not valid
                msg = wx.MessageDialog(self, 'The values entered are not valid', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_line_broadening_textcontrol_dim3.SetValue(str(self.exponential_line_broadening_dim3))
                return
            self.exponential_line_broadening_dim3 = em
            
            self.line1.set_ydata(np.exp(-(np.pi*x*self.exponential_line_broadening_dim3)))
        elif(self.apodization_dim3_combobox_selection==2):
            try:
                g1 = float(self.apodization_g1_textcontrol_dim3.GetValue())
                g2 = float(self.apodization_g2_textcontrol_dim3.GetValue())
                g3 = float(self.apodization_g3_textcontrol_dim3.GetValue())
            except:
                # Give a popout window saying that the values are not valid
                msg = wx.MessageDialog(self, 'The values entered are not valid', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_g1_textcontrol_dim3.SetValue(str(self.g1_dim3))
                self.apodization_g2_textcontrol_dim3.SetValue(str(self.g2_dim3))
                self.apodization_g3_textcontrol_dim3.SetValue(str(self.g3_dim3))
                return
            # Check to see if g3 is between 0 and 1
            if(g3 < 0 or g3 > 1):
                msg = wx.MessageDialog(self, 'Gaussian shift must be between 0 and 1', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_g3_textcontrol_dim3.SetValue(str(self.g3_dim3))
                return
            self.g1_dim3 = g1
            self.g2_dim3 = g2
            self.g3_dim3 = g3
            e = np.pi * (self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2] * self.g1_dim3
            g = 0.6 * np.pi * self.g2_dim3 * (self.g3_dim3 * ((self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2] - 1) - x)
            func = np.exp(e - g * g)
            self.line1.set_ydata(func)

            self.apodization_plot_ax_dim3.set_xlim(0, (self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2])


        elif(self.apodization_dim3_combobox_selection==3):
            try:
                offset = float(self.apodization_offset_textcontrol_dim3.GetValue())
                end = float(self.apodization_end_textcontrol_dim3.GetValue())
                power = float(self.apodization_power_textcontrol_dim3.GetValue())
                power = int(power)
            except:
                msg = wx.MessageDialog(self, 'The values entered are not valid', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_offset_textcontrol_dim3.SetValue(str(self.offset_dim3))
                self.apodization_end_textcontrol_dim3.SetValue(str(self.end_dim3))
                self.apodization_power_textcontrol_dim3.SetValue(str(self.power_dim3))
                return
            # Check that offset and end are between 0 and 1
            if(offset < 0 or offset > 1):
                msg = wx.MessageDialog(self, 'Offset values must be between 0 and 1', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_offset_textcontrol_dim3.SetValue(str(self.offset_dim3))
                return
            if(end < 0 or end > 1):
                msg = wx.MessageDialog(self, 'End values must be between 0 and 1', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_end_textcontrol_dim3.SetValue(str(self.end_dim3))
                return
            # Check that power is greater than 0
            if(power < 0):
                msg = wx.MessageDialog(self, 'Power must be greater than 0', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_power_textcontrol_dim3.SetValue(str(self.power_dim3))
                return
            self.offset_dim3 = offset
            self.end_dim3 = end
            self.power_dim3 = power
            func = np.sin((np.pi*self.offset_dim3 + np.pi*(self.end_dim3-self.offset_dim3)*x)/((((self.nmr_data.number_of_points[2]/2)/self.nmr_data.spectral_width[2]))))**self.power_dim3
            self.line1.set_ydata(func)
        elif(self.apodization_dim3_combobox_selection==4):
            try:
                a = float(self.apodization_a_textcontrol_dim3.GetValue())
                b = float(self.apodization_b_textcontrol_dim3.GetValue())
            except:
                msg = wx.MessageDialog(self, 'The values entered are not valid', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_a_textcontrol_dim3.SetValue(str(self.a_dim3))
                self.apodization_b_textcontrol_dim3.SetValue(str(self.b_dim3))
                return
            self.a_dim3 = a
            self.b_dim3 = b
            func = np.exp(-self.a_dim3*(x**2) - self.b_dim3*x)
            self.line1.set_ydata(func)
        elif(self.apodization_dim3_combobox_selection==5):
            try:
                t1 = float(self.apodization_t1_textcontrol_dim3.GetValue())
                t2 = float(self.apodization_t2_textcontrol_dim3.GetValue())
            except:
                msg = wx.MessageDialog(self, 'The values entered are not valid', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_t1_textcontrol_dim3.SetValue(str(self.t1_dim3))
                self.apodization_t2_textcontrol_dim3.SetValue(str(self.t2_dim3))
                return
            # Ensure that t1 and t2 are greater than 0
            if(t1 < 0 or t2 < 0):
                msg = wx.MessageDialog(self, 'Ramp up and ramp down points must be greater than 0', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_t1_textcontrol_dim3.SetValue(str(self.t1_dim3))
                self.apodization_t2_textcontrol_dim3.SetValue(str(self.t2_dim3))
                return
            # Ensure that t1 + t2 is less than the number of points
            if(t1 + t2 > self.nmr_data.number_of_points[2]):
                message = 'Ramp up and ramp down points must be less than the number of points (' + str(self.nmr_data.number_of_points[2]) + ')'
                msg = wx.MessageDialog(self, message, 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_t1_textcontrol_dim3.SetValue(str(self.t1_dim3))
                self.apodization_t2_textcontrol_dim3.SetValue(str(self.t2_dim3))
                return
            self.t1_dim3 = t1
            self.t2_dim3 = t2
            func = np.concatenate((np.linspace(0, 1, int(self.t1_dim3)), np.ones(int(self.nmr_data.number_of_points[2]/2) - int(self.t1_dim3) - int(self.t2_dim3)),np.linspace(1, 0, int(self.t2_dim3))))
            self.line1.set_ydata(func)
        elif(self.apodization_dim3_combobox_selection==6):
            try:
                loc = float(self.apodization_loc_textcontrol_dim3.GetValue())
            except:
                msg = wx.MessageDialog(self, 'The values entered are not valid', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_loc_textcontrol_dim3.SetValue(str(self.loc_dim3))
                return
            # Ensure that loc is between 0 and 1
            if(loc < 0 or loc > 1):
                msg = wx.MessageDialog(self, 'Location of maximum must be between 0 and 1', 'Error', wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                self.apodization_loc_textcontrol_dim3.SetValue(str(self.loc_dim3))
                return
            self.loc_dim3 = loc
            func = np.concatenate((np.linspace(0, 1, int(self.loc_dim3*(self.nmr_data.number_of_points[2]/2))), np.linspace(1, 0, int(self.nmr_data.number_of_points[2]/2) - int(self.loc_dim3*(self.nmr_data.number_of_points[2])/2))))
            self.line1.set_ydata(func)

        self.apodization_plot_canvas.draw()


    def on_apodization_combobox_dim3(self,event):
        self.apodization_dim3_combobox_selection= self.apodization_combobox_dim3.GetSelection()

        # Destroy the combobox and textcontrols for the previous apodization function
        # self.apodization_sizer.Detach(self.apodization_combobox)
        # self.apodization_combobox.Destroy()

        # # Remove the zf sizer
        self.zero_filling_sizer_dim3.Clear(delete_windows=True)

        if(self.apodization_dim3_combobox_selection_old == 1):
            # Remove the previous textcontrols
            
            self.apodization_sizer_dim3.Detach(self.apodization_line_broadening_label)
            self.apodization_line_broadening_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_line_broadening_textcontrol_dim3)
            self.apodization_line_broadening_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_textcontrol_dim3)
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)

        elif(self.apodization_dim3_combobox_selection_old == 2):
            self.apodization_sizer_dim3.Detach(self.apodization_g1_label)
            self.apodization_g1_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g1_textcontrol_dim3)
            self.apodization_g1_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g2_label)
            self.apodization_g2_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g2_textcontrol_dim3)
            self.apodization_g2_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g3_label)
            self.apodization_g3_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g3_textcontrol_dim3)
            self.apodization_g3_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_textcontrol_dim3)
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)

        elif(self.apodization_dim3_combobox_selection_old == 3):
            self.apodization_sizer_dim3.Detach(self.apodization_offset_label)
            self.apodization_offset_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_offset_textcontrol_dim3)
            self.apodization_offset_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_end_label)
            self.apodization_end_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_end_textcontrol_dim3)
            self.apodization_end_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_power_label)
            self.apodization_power_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_power_textcontrol_dim3)
            self.apodization_power_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_textcontrol_dim3)
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif(self.apodization_dim3_combobox_selection_old == 4):
            self.apodization_sizer_dim3.Detach(self.apodization_a_label)
            self.apodization_a_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_a_textcontrol_dim3)
            self.apodization_a_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_b_label)
            self.apodization_b_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_b_textcontrol_dim3)
            self.apodization_b_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_textcontrol_dim3)
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif(self.apodization_dim3_combobox_selection_old == 5):
            self.apodization_sizer_dim3.Detach(self.apodization_t1_label)
            self.apodization_t1_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_t1_textcontrol_dim3)
            self.apodization_t1_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_t2_label)
            self.apodization_t2_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_t2_textcontrol_dim3)
            self.apodization_t2_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_textcontrol_dim3)
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif(self.apodization_dim3_combobox_selection_old == 6):
            self.apodization_sizer_dim3.Detach(self.apodization_loc_label)
            self.apodization_loc_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_loc_textcontrol_dim3)
            self.apodization_loc_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_textcontrol_dim3)
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif(self.apodization_dim3_combobox_selection_old == 0):
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        


        self.apodization_sizer_dim3.Detach(self.apodization_checkbox_dim3)
        self.apodization_checkbox_dim3.Destroy()
        self.apodization_sizer_dim3.Detach(self.apodization_combobox_dim3)
        self.apodization_combobox_dim3.Hide()
       


        

        # Delete the current apodization sizer and then create a new one
        
        self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
        self.apodization_plot_ax_dim3.clear()
        self.apodization_plot_ax_dim3.clear()
        self.apodization_plot_sizer_dim3.Clear(True)

        self.sizer_2.Remove(self.apodization_sizer_dim3)
        # self.apodization_sizer.Clear(delete_windows=True)

 


        



        # Remove the linear prediction sizers
        self.linear_prediction_sizer_dim3.Clear(delete_windows=True)
        # self.sizer_1.Remove(self.linear_prediction_sizer)
        

        # self.sizer_1.Remove(self.solvent_suppression_sizer)

        
        self.sizer_2.Clear(delete_windows=True)

           
        self.create_menu_bar_dim3()
        self.Refresh()
        self.Update()
        self.Layout()




        self.apodization_dim3_combobox_selection_old = self.apodization_dim3_combobox_selection



    def create_zero_filling_sizer_dim3(self, parent):
        # Create a box for zero filling options
        self.zero_filling_box_dim3 = wx.StaticBox(parent, -1, 'Zero Filling')
        self.zero_filling_sizer_dim3 = wx.StaticBoxSizer(self.zero_filling_box_dim3, wx.HORIZONTAL)
        self.zero_filling_checkbox_dim3 = wx.CheckBox(parent, -1, 'Apply zero filling')
        self.zero_filling_checkbox_dim3.SetValue(self.zero_filling_dim3_checkbox_value)
        self.zero_filling_checkbox_dim3.Bind(wx.EVT_CHECKBOX, self.on_zero_filling_checkbox_dim3)
        self.zero_filling_sizer_dim3.Add(self.zero_filling_checkbox_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.zero_filling_sizer_dim3.AddSpacer(10)
        # Have a combobox for zero filling options
        self.zf_options_label = wx.StaticText(parent, -1, 'Options:')
        self.zero_filling_sizer_dim3.Add(self.zf_options_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.zero_filling_sizer_dim3.AddSpacer(5)
        self.zero_filling_options_dim3 = ['Doubling spectrum size', 'Adding Zeros', 'Final data size']
        self.zero_filling_combobox_dim3 = wx.ComboBox(parent, -1, choices=self.zero_filling_options_dim3, style=wx.CB_READONLY)
        self.zero_filling_combobox_dim3.Bind(wx.EVT_COMBOBOX, self.on_zero_filling_combobox_dim3)
        self.zero_filling_combobox_dim3.SetSelection(self.zero_filling_dim3_combobox_selection)
        self.zero_filling_sizer_dim3.Add(self.zero_filling_combobox_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.zero_filling_sizer_dim3.AddSpacer(10)
        if(self.zero_filling_dim3_combobox_selection == 0):
            # Have a textcontrol for the doubling number/number of zeros/final data size
            self.zf_value_label = wx.StaticText(parent, -1, 'Doubling number:')
            self.zero_filling_sizer_dim3.Add(self.zf_value_label, 0, wx.ALIGN_CENTER_VERTICAL)

            self.zero_filling_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.zero_filling_dim3_value_doubling_times), size=(40, 20))
            self.zero_filling_sizer_dim3.AddSpacer(5)
            self.zero_filling_sizer_dim3.Add(self.zero_filling_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.zero_filling_textcontrol_dim3.Bind(wx.EVT_TEXT, self.on_zero_filling_doubling_number_dim3)
            self.zero_filling_sizer_dim3.AddSpacer(20)
        elif(self.zero_filling_dim3_combobox_selection == 1):
            # Have a textcontrol for the doubling number/number of zeros/final data size
            self.zf_value_label = wx.StaticText(parent, -1, 'Number of zeros to add:')
            self.zero_filling_sizer_dim3.Add(self.zf_value_label, 0, wx.ALIGN_CENTER_VERTICAL)

            self.zero_filling_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.zero_filling_dim3_value_zeros_to_add), size=(40, 20))
            self.zero_filling_sizer_dim3.AddSpacer(5)
            self.zero_filling_sizer_dim3.Add(self.zero_filling_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.zero_filling_textcontrol_dim3.Bind(wx.EVT_TEXT, self.on_zero_filling_zeros_to_add_dim3)    
            self.zero_filling_sizer_dim3.AddSpacer(20)
        elif(self.zero_filling_dim3_combobox_selection == 2):
            # Have a textcontrol for the doubling number/number of zeros/final data size
            self.zf_value_label = wx.StaticText(parent, -1, 'Final data size:')
            self.zero_filling_sizer_dim3.Add(self.zf_value_label, 0, wx.ALIGN_CENTER_VERTICAL)

            self.zero_filling_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.zero_filling_dim3_value_final_data_size), size=(40, 20))
            self.zero_filling_textcontrol_dim3.Bind(wx.EVT_TEXT, self.on_zero_filling_final_size_dim3)
            self.zero_filling_sizer_dim3.AddSpacer(5)
            self.zero_filling_sizer_dim3.Add(self.zero_filling_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
            self.zero_filling_sizer_dim3.AddSpacer(20)

        # Have a checkbox for rounding to the nearest power of 2
        self.zero_filling_round_checkbox_dim3 = wx.CheckBox(parent, -1, 'Round to nearest power of 2')
        self.zero_filling_round_checkbox_dim3.SetValue(self.zero_filling_dim3_round_checkbox_value)
        self.zero_filling_round_checkbox_dim3.Bind(wx.EVT_CHECKBOX, self.on_zero_filling_round_checkbox_dim3)
        self.zero_filling_sizer_dim3.Add(self.zero_filling_round_checkbox_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.zero_filling_sizer_dim3.AddSpacer(10)
        

        # Have a button showing information on zero filling
        self.zero_filling_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))
        self.zero_filling_info.Bind(wx.EVT_BUTTON, self.oneDFrame.on_zero_fill_info)
        self.zero_filling_sizer_dim3.Add(self.zero_filling_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.zero_filling_sizer_dim3.AddSpacer(10)

        self.sizer_2.Add(self.zero_filling_sizer_dim3)
        self.sizer_2.AddSpacer(10)

    def on_zero_filling_checkbox_dim3(self, event):
        self.zero_filling_dim3_checkbox_value = self.zero_filling_checkbox_dim3.GetValue()


    def on_zero_filling_round_checkbox_dim3(self, event):
        self.zero_filling_dim3_round_checkbox_value = self.zero_filling_round_checkbox_dim3.GetValue()


    def on_zero_filling_final_size_dim3(self, event):
        self.zero_filling_dim3_value_final_data_size = self.zero_filling_textcontrol_dim3.GetValue()

    def on_zero_filling_zeros_to_add_dim3(self, event):
        self.zero_filling_dim3_value_zeros_to_add = self.zero_filling_textcontrol_dim3.GetValue()

    def on_zero_filling_doubling_number_dim3(self, event):
        self.zero_filling_dim3_value_doubling_times = self.zero_filling_textcontrol_dim3.GetValue()



    def on_zero_filling_combobox_dim3(self, event):
        self.zero_filling_dim3_combobox_selection = self.zero_filling_combobox_dim3.GetSelection()
        # # # Remove the zf sizer
        self.zero_filling_sizer_dim3.Clear()
        self.zero_filling_sizer_dim3.Detach(self.zero_filling_checkbox_dim3)
        self.zero_filling_checkbox_dim3.Destroy()
        self.zero_filling_sizer_dim3.Detach(self.zf_options_label)
        self.zf_options_label.Destroy()
        self.zero_filling_sizer_dim3.Detach(self.zero_filling_info)
        self.zero_filling_info.Destroy()
        self.zero_filling_sizer_dim3.Detach(self.zf_value_label)
        self.zf_value_label.Destroy()
        self.zero_filling_sizer_dim3.Detach(self.zero_filling_round_checkbox_dim3)
        self.zero_filling_round_checkbox_dim3.Destroy()
        self.zero_filling_sizer_dim3.Detach(self.zero_filling_textcontrol_dim3)
        self.zero_filling_textcontrol_dim3.Destroy()

        self.zero_filling_sizer_dim3.Detach(self.zero_filling_combobox_dim3)
        self.zero_filling_combobox_dim3.Hide()

        if(self.apodization_dim3_combobox_selection_old == 1):
            # Remove the previous textcontrols
            
            self.apodization_sizer_dim3.Detach(self.apodization_line_broadening_label)
            self.apodization_line_broadening_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_line_broadening_textcontrol_dim3)
            self.apodization_line_broadening_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_textcontrol_dim3)
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)

        elif(self.apodization_dim3_combobox_selection_old == 2):
            self.apodization_sizer_dim3.Detach(self.apodization_g1_label)
            self.apodization_g1_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g1_textcontrol_dim3)
            self.apodization_g1_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g2_label)
            self.apodization_g2_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g2_textcontrol_dim3)
            self.apodization_g2_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g3_label)
            self.apodization_g3_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_g3_textcontrol_dim3)
            self.apodization_g3_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_textcontrol_dim3)
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)

        elif(self.apodization_dim3_combobox_selection_old == 3):
            self.apodization_sizer_dim3.Detach(self.apodization_offset_label)
            self.apodization_offset_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_offset_textcontrol_dim3)
            self.apodization_offset_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_end_label)
            self.apodization_end_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_end_textcontrol_dim3)
            self.apodization_end_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_power_label)
            self.apodization_power_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_power_textcontrol_dim3)
            self.apodization_power_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_textcontrol_dim3)
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif(self.apodization_dim3_combobox_selection_old == 4):
            self.apodization_sizer_dim3.Detach(self.apodization_a_label)
            self.apodization_a_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_a_textcontrol_dim3)
            self.apodization_a_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_b_label)
            self.apodization_b_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_b_textcontrol_dim3)
            self.apodization_b_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_textcontrol_dim3)
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif(self.apodization_dim3_combobox_selection_old == 5):
            self.apodization_sizer_dim3.Detach(self.apodization_t1_label)
            self.apodization_t1_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_t1_textcontrol_dim3)
            self.apodization_t1_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_t2_label)
            self.apodization_t2_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_t2_textcontrol_dim3)
            self.apodization_t2_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_textcontrol_dim3)
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif(self.apodization_dim3_combobox_selection_old == 6):
            self.apodization_sizer_dim3.Detach(self.apodization_loc_label)
            self.apodization_loc_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_loc_textcontrol_dim3)
            self.apodization_loc_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_label)
            self.apodization_first_point_label.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_first_point_textcontrol_dim3)
            self.apodization_first_point_textcontrol_dim3.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        elif(self.apodization_dim3_combobox_selection_old == 0):
            self.apodization_sizer_dim3.Detach(self.apodization_info)
            self.apodization_info.Destroy()
            self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
            self.apodization_plot_sizer_dim3.Clear(True)
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_ax_dim3.clear()
            self.apodization_plot_sizer_dim3.Clear(True)
        


        self.apodization_sizer_dim3.Detach(self.apodization_checkbox_dim3)
        self.apodization_checkbox_dim3.Destroy()
        self.apodization_sizer_dim3.Detach(self.apodization_combobox_dim3)
        self.apodization_combobox_dim3.Hide()
       


        

        # Delete the current apodization sizer and then create a new one
        
        self.apodization_sizer_dim3.Detach(self.apodization_plot_sizer_dim3)
        self.apodization_plot_ax_dim3.clear()
        self.apodization_plot_ax_dim3.clear()
        self.apodization_plot_sizer_dim3.Clear(True)

        self.sizer_2.Remove(self.apodization_sizer_dim3)
        # self.apodization_sizer.Clear(delete_windows=True)

 


        



        # Remove the linear prediction sizers
        self.linear_prediction_sizer_dim3.Clear(delete_windows=True)
        # self.sizer_1.Remove(self.linear_prediction_sizer)
        

        # self.sizer_1.Remove(self.solvent_suppression_sizer)

        
        self.sizer_2.Clear(delete_windows=True)

        self.create_menu_bar_dim3()
        self.Refresh()
        self.Update()
        self.Layout()





    def create_fourier_transform_sizer_dim3(self, parent):
        # Create a box for fourier transform options
        self.fourier_transform_box = wx.StaticBox(parent, -1, 'Fourier Transform')
        self.fourier_transform_sizer_dim3 = wx.StaticBoxSizer(self.fourier_transform_box, wx.HORIZONTAL)
        self.fourier_transform_checkbox_dim3 = wx.CheckBox(parent, -1, 'Apply fourier transform')
        self.fourier_transform_checkbox_dim3.Bind(wx.EVT_CHECKBOX, self.on_fourier_transform_checkbox_dim3)
        self.fourier_transform_checkbox_dim3.SetValue(self.fourier_transform_dim3_checkbox_value)
        self.fourier_transform_sizer_dim3.Add(self.fourier_transform_checkbox_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.fourier_transform_sizer_dim3.AddSpacer(10)
        # Have a button for advanced options for fourier transform
        self.fourier_transform_advanced_options_dim3 = wx.Button(parent, -1, 'Advanced Options')
        self.fourier_transform_advanced_options_dim3.Bind(wx.EVT_BUTTON, self.on_fourier_transform_advanced_options_dim3)
        self.fourier_transform_sizer_dim3.Add(self.fourier_transform_advanced_options_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.fourier_transform_sizer_dim3.AddSpacer(10)

        # Have a button showing information on fourier transform
        self.fourier_transform_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))
        self.fourier_transform_info.Bind(wx.EVT_BUTTON, self.oneDFrame.on_fourier_transform_info)
        self.fourier_transform_sizer_dim3.Add(self.fourier_transform_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer_2.Add(self.fourier_transform_sizer_dim3)
        self.sizer_2.AddSpacer(10)

    def on_fourier_transform_checkbox_dim3(self, event):
        self.fourier_transform_dim3_checkbox_value = self.fourier_transform_checkbox_dim3.GetValue()


    def on_fourier_transform_advanced_options_dim3(self, event):
        # Create a frame with a set of advanced options for the fourier transform implementation
        self.fourier_transform_advanced_options_window_dim3 = wx.Frame(self, -1, 'Fourier Transform Advanced Options (Dimension 2)', size=(700, 300))
        if(darkdetect.isDark() == True and platform!='windows'):
            self.fourier_transform_advanced_options_window_dim3.SetBackgroundColour((53, 53, 53, 255))
        else:
            self.fourier_transform_advanced_options_window_dim3.SetBackgroundColour('White')
        
        self.fourier_transform_advanced_options_window_sizer_dim3 = wx.BoxSizer(wx.VERTICAL)
        self.fourier_transform_advanced_options_window_dim3.SetSizer(self.fourier_transform_advanced_options_window_sizer_dim3)

        # Create a sizer for the fourier transform advanced options
        self.ft_label = wx.StaticBox(self.fourier_transform_advanced_options_window_dim3, -1, 'Fourier Transform Method:')
        self.fourier_transform_advanced_options_sizer_dim3 = wx.StaticBoxSizer(self.ft_label,wx.VERTICAL)

        # Have a radiobox for auto, real, inverse, sign alternation
        self.fourier_transform_advanced_options_sizer_dim3.AddSpacer(10)
        self.fourier_transform_auto_real_inverse_sign_alternation_radio_box_dim3 = wx.RadioBox(self.fourier_transform_advanced_options_window_dim3, -1, choices=['Auto', 'Real', 'Inverse', 'Sign Alternation', 'Negative'], style=wx.RA_SPECIFY_COLS)
        self.fourier_transform_auto_real_inverse_sign_alternation_radio_box_dim3.SetSelection(self.ft_method_selection_dim3)
        self.fourier_transform_advanced_options_sizer_dim3.Add(self.fourier_transform_auto_real_inverse_sign_alternation_radio_box_dim3, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.fourier_transform_advanced_options_sizer_dim3.AddSpacer(10)


        self.ft_method_text = 'Auto: The auto method will automatically select the best method for the fourier transform of the FID. \n\nReal: The Fourier Transform will be applied to the real part of the FID only. \n\nInverse: The inverse Fourier Transform will be applied to the FID. \n\nSign Alternation: The sign alternation method will be applied to the FID. \n\n'

        self.ft_method_info = wx.StaticText(self.fourier_transform_advanced_options_window_dim3, -1, self.ft_method_text)
        self.fourier_transform_advanced_options_sizer_dim3.Add(self.ft_method_info, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.fourier_transform_advanced_options_sizer_dim3.AddSpacer(10)

        # Have a save and close button
        self.fourier_transform_advanced_options_save_button_dim3 = wx.Button(self.fourier_transform_advanced_options_window_dim3, -1, 'Save and Close')
        self.fourier_transform_advanced_options_save_button_dim3.Bind(wx.EVT_BUTTON, self.on_fourier_transform_advanced_options_save_dim3)
        self.fourier_transform_advanced_options_sizer_dim3.Add(self.fourier_transform_advanced_options_save_button_dim3, 0, wx.ALIGN_CENTER_HORIZONTAL)



        self.fourier_transform_advanced_options_window_sizer_dim3.Add(self.fourier_transform_advanced_options_sizer_dim3, 0, wx.ALIGN_CENTER_HORIZONTAL)

        self.fourier_transform_advanced_options_window_dim3.Show()

    def on_fourier_transform_advanced_options_save_dim3(self, event):
        # Save the current selection and close the window
        self.ft_method_selection_dim3 = self.fourier_transform_auto_real_inverse_sign_alternation_radio_box_dim3.GetSelection()
        self.fourier_transform_advanced_options_window_dim3.Close()


    def create_phase_correction_sizer_dim3(self, parent):
        # Create a box for phase correction options
        self.phase_correction_box_dim3 = wx.StaticBox(parent, -1, 'Phase Correction')
        self.phase_correction_sizer_dim3 = wx.StaticBoxSizer(self.phase_correction_box_dim3, wx.HORIZONTAL)
        self.phase_correction_checkbox_dim3 = wx.CheckBox(parent, -1, 'Apply phase correction')
        self.phase_correction_checkbox_dim3.Bind(wx.EVT_CHECKBOX, self.on_phase_correction_checkbox_dim3)
        self.phase_correction_checkbox_dim3.SetValue(self.phasing_dim3_checkbox_value)
        self.phase_correction_sizer_dim3.Add(self.phase_correction_checkbox_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.phase_correction_sizer_dim3.AddSpacer(10)
        # Have a textcontrol for p0 and p1 values
        self.phase_correction_p0_label = wx.StaticText(parent, -1, 'Zero order correction (p0):')
        self.phase_correction_sizer_dim3.Add(self.phase_correction_p0_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.phase_correction_p0_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.p0_total_dim3), size=(50, 20))
        self.phase_correction_p0_textcontrol_dim3.Bind(wx.EVT_TEXT, self.on_phase_correction_p0_dim3)
        self.phase_correction_sizer_dim3.Add(self.phase_correction_p0_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.phase_correction_sizer_dim3.AddSpacer(10)
        self.phase_correction_p1_label = wx.StaticText(parent, -1, 'First order correction (p1):')
        self.phase_correction_sizer_dim3.Add(self.phase_correction_p1_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.phase_correction_p1_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.p1_total_dim3), size=(50, 20))
        self.phase_correction_p1_textcontrol_dim3.Bind(wx.EVT_TEXT, self.on_phase_correction_p1_dim3)
        self.phase_correction_sizer_dim3.Add(self.phase_correction_p1_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.phase_correction_sizer_dim3.AddSpacer(10)

        # Have a checkbox for f1180
        self.phase_correction_f1180_button_dim3 = wx.CheckBox(parent, -1, 'F1180')
        self.phase_correction_f1180_button_dim3.Bind(wx.EVT_CHECKBOX, self.on_phase_correction_f1180_dim3)
        self.phase_correction_f1180_button_dim3.SetValue(self.f1180_dim3)
        self.phase_correction_sizer_dim3.Add(self.phase_correction_f1180_button_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.phase_correction_sizer_dim3.AddSpacer(10)

        # Have a button showing information on phase correction
        self.phase_correction_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))
        self.phase_correction_info.Bind(wx.EVT_BUTTON, self.on_phase_correction_info_dim3)
        self.phase_correction_sizer_dim3.Add(self.phase_correction_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer_2.Add(self.phase_correction_sizer_dim3)
        self.sizer_2.AddSpacer(10)


    def on_phase_correction_checkbox_dim3(self, event):
        self.phasing_dim3_checkbox_value = self.phase_correction_checkbox_dim3.GetValue()



    def on_phase_correction_p0_dim3(self, event):
        self.p0_total_dim3 = self.phase_correction_p0_textcontrol_dim3.GetValue()

    def on_phase_correction_p1_dim3(self, event):
        self.p1_total_dim3 = self.phase_correction_p1_textcontrol_dim3.GetValue()
        try:
            if(np.abs(float(self.p1_total_dim3)) > 45):
                self.apodization_first_point_scaling_dim3 = 1.0
                self.apodization_first_point_textcontrol_dim3.SetValue(str(self.apodization_first_point_scaling_dim3))
            else:
                self.apodization_first_point_scaling_dim3 = 0.0
                self.apodization_first_point_textcontrol_dim3.SetValue(str(self.apodization_first_point_scaling_dim3))
        except:
            pass


    def on_phase_correction_f1180_dim3(self, event):
        if(self.phase_correction_f1180_button_dim3.GetValue() == True):
            self.c_old = self.apodization_first_point_scaling_dim3
            self.p0_total_dim3_old = self.p0_total_dim3
            self.p1_total_dim3_old = self.p1_total_dim3
            # Apply -90 p0 and 180 p1 to the phase correction textcontrols
            self.p0_total_dim3 = -90.0
            self.p1_total_dim3 = 180.0
            self.phase_correction_p0_textcontrol_dim3.SetValue(str(self.p0_total_dim3))
            self.phase_correction_p1_textcontrol_dim3.SetValue(str(self.p1_total_dim3))
            # Disable the phase correction textcontrols
            self.phase_correction_p0_textcontrol_dim3.Disable()
            self.phase_correction_p1_textcontrol_dim3.Disable()

            self.apodization_first_point_scaling_dim3 = 1.0
            self.apodization_first_point_textcontrol_dim3.SetValue(str(self.apodization_first_point_scaling_dim3))
        else:
            self.p0_total_dim3 = self.p0_total_dim3_old
            self.p1_total_dim3 = self.p1_total_dim3_old
            self.phase_correction_p0_textcontrol_dim3.SetValue(str(self.p0_total_dim3))
            self.phase_correction_p1_textcontrol_dim3.SetValue(str(self.p1_total_dim3))
            self.phase_correction_p0_textcontrol_dim3.Enable()
            self.phase_correction_p1_textcontrol_dim3.Enable()
            self.apodization_first_point_scaling_dim3 = self.c_old
            self.apodization_first_point_textcontrol_dim3.SetValue(str(self.apodization_first_point_scaling_dim3))


    def on_phase_correction_info_dim3(self, event):
        phase_correction_text = 'Phase correction is a method to correct for phase errors in the FID. Zero order phase correction (p0) is used to correct a phase offset that is applied equally across the spectrum. However, a first order phase correction (p1) is used to correct the phasing in a spectrum where peaks in different locations of the spectrum require a different phasing value. For the indirect dimension, it is often the case that the acquisition is delayed by an exact time so that the resulting spectrum can be phased using the phase values of: p0=-90, p1=180. This is often termed F1180. \n Further information can be found using the link below.'

        # Create a popup window with the information
        self.phase_correction_info_window = wx.Frame(self, -1, 'Phase Correction Information', size=(450, 300))
        if(darkdetect.isDark() == True and platform!='windows'):
            self.phase_correction_info_window.SetBackgroundColour((53, 53, 53, 255))
            colour = "RED"
        else:
            self.phase_correction_info_window.SetBackgroundColour('White')
            colour = "BLUE"
        
        self.phase_correction_info_window_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.phase_correction_info_window.SetSizer(self.phase_correction_info_window_sizer)

        self.phase_correction_info_window_sizer.AddSpacer(10)

        # Create a sizer for the phase correction information
        self.phase_correction_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.phase_correction_info_sizer.AddSpacer(10)
        self.phase_correction_info_sizer.Add(wx.StaticText(self.phase_correction_info_window, -1, phase_correction_text, size=(400,200)), 0, wx.ALIGN_CENTER)
        self.phase_correction_info_sizer.AddSpacer(10)

        # Have a hyperlink to the phase correction information
        self.phase_correction_info_hyperlink = hl.HyperLinkCtrl(self.phase_correction_info_window, -1, 'NMRPipe Help Page for Phase Correction', URL='http://www.nmrscience.com/ref/nmrpipe/ps.html')
        self.phase_correction_info_hyperlink.SetColours(colour, colour, colour)
        self.phase_correction_info_hyperlink.SetUnderlines(False, False, False)
        self.phase_correction_info_hyperlink.SetBold(False)
        self.phase_correction_info_hyperlink.UpdateLink()
        self.phase_correction_info_sizer.Add(self.phase_correction_info_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.phase_correction_info_sizer.AddSpacer(10)

        self.phase_correction_info_window_sizer.Add(self.phase_correction_info_sizer, 0, wx.ALIGN_CENTER)

        self.phase_correction_info_window.Show()


    def create_extraction_sizer_dim3(self, parent):
        # A box for extraction of data between two ppm values
        self.extraction_box_dim3 = wx.StaticBox(parent, -1, 'Extraction')
        self.extraction_sizer_dim3 = wx.StaticBoxSizer(self.extraction_box_dim3, wx.HORIZONTAL)
        self.extraction_checkbox_dim3 = wx.CheckBox(parent, -1, 'Include data extraction')
        self.extraction_checkbox_dim3.Bind(wx.EVT_CHECKBOX, self.on_extraction_checkbox_dim3)
        self.extraction_checkbox_dim3.SetValue(self.extraction_checkbox_value_dim3)
        self.extraction_sizer_dim3.Add(self.extraction_checkbox_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.extraction_sizer_dim3.AddSpacer(10)
        # Have a textcontrol for the ppm start value
        self.extraction_ppm_start_label = wx.StaticText(parent, -1, 'Start chemical shift (ppm):')
        self.extraction_sizer_dim3.Add(self.extraction_ppm_start_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.extraction_ppm_start_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.extraction_start_dim3), size=(40, 20))
        self.extraction_ppm_start_textcontrol_dim3.Bind(wx.EVT_TEXT, self.on_extraction_dim3)
        self.extraction_sizer_dim3.Add(self.extraction_ppm_start_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.extraction_sizer_dim3.AddSpacer(10)
        # Have a textcontrol for the ppm end value
        self.extraction_ppm_end_label = wx.StaticText(parent, -1, 'End chemical shift (ppm):')
        self.extraction_sizer_dim3.Add(self.extraction_ppm_end_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.extraction_ppm_end_textcontrol_dim3 = wx.TextCtrl(parent, -1, str(self.extraction_end_dim3), size=(40, 20))
        self.extraction_ppm_end_textcontrol_dim3.Bind(wx.EVT_TEXT, self.on_extraction_dim3)
        self.extraction_sizer_dim3.Add(self.extraction_ppm_end_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.extraction_sizer_dim3.AddSpacer(10)
        # Have a button showing information on extraction
        self.extraction_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))
        self.extraction_info.Bind(wx.EVT_BUTTON, self.oneDFrame.on_extraction_info)
        self.extraction_sizer_dim3.Add(self.extraction_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer_2.Add(self.extraction_sizer_dim3)
        self.sizer_2.AddSpacer(10)

    def on_extraction_checkbox_dim3(self, event):
        self.extraction_checkbox_value_dim3 = self.extraction_checkbox_dim3.GetValue()



    def on_extraction_dim3(self, event):
        self.extraction_start_dim3 = self.extraction_ppm_start_textcontrol_dim3.GetValue()
        self.extraction_end_dim3 = self.extraction_ppm_end_textcontrol_dim3.GetValue()




    def create_baseline_correction_sizer_dim3(self, parent):
        # Create a box for baseline correction options (linear/polynomial)
        self.baseline_correction_box_dim3 = wx.StaticBox(parent, -1, 'Baseline Correction')
        self.baseline_correction_sizer_dim3 = wx.StaticBoxSizer(self.baseline_correction_box_dim3, wx.HORIZONTAL)
        self.baseline_correction_checkbox_dim3 = wx.CheckBox(parent, -1, 'Apply baseline correction')
        self.baseline_correction_checkbox_dim3.Bind(wx.EVT_CHECKBOX, self.on_baseline_correction_checkbox_dim3)
        self.baseline_correction_checkbox_dim3.SetValue(self.baseline_correction_checkbox_value_dim3)
        self.baseline_correction_sizer_dim3.Add(self.baseline_correction_checkbox_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_sizer_dim3.AddSpacer(10)
        # Have a radio box for linear or polynomial baseline correction
        self.baseline_correction_radio_box_dim3 = wx.RadioBox(parent, -1, 'Baseline Correction Method', choices=['Linear', 'Polynomial'])
        # Bind the radio box to a function that will update the baseline correction options
        self.baseline_correction_radio_box_dim3.Bind(wx.EVT_RADIOBOX, self.on_baseline_correction_radio_box_dim3)
        self.baseline_correction_radio_box_dim3.SetSelection(self.baseline_correction_radio_box_selection_dim3)
        self.baseline_correction_sizer_dim3.Add(self.baseline_correction_radio_box_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_sizer_dim3.AddSpacer(10)
        
        # If linear baseline correction is selected, have a textcontrol for the node values to use
        self.baseline_correction_nodes_label = wx.StaticText(parent, -1, 'Node width (pts):')
        self.baseline_correction_sizer_dim3.Add(self.baseline_correction_nodes_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_nodes_textcontrol_dim3 = wx.TextCtrl(parent, -1, self.node_width_dim3, size=(30, 20))
        self.baseline_correction_nodes_textcontrol_dim3.Bind(wx.EVT_TEXT, self.on_baseline_correction_textcontrol_dim3)
        self.baseline_correction_sizer_dim3.Add(self.baseline_correction_nodes_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_sizer_dim3.AddSpacer(10)
        # Have a textcontrol for the node list (percentages)
        self.baseline_correction_node_list_label = wx.StaticText(parent, -1, 'Node list (%):')
        self.baseline_correction_sizer_dim3.Add(self.baseline_correction_node_list_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_node_list_textcontrol_dim3 = wx.TextCtrl(parent, -1, self.node_list_dim3, size=(100, 20))
        self.baseline_correction_node_list_textcontrol_dim3.Bind(wx.EVT_TEXT, self.on_baseline_correction_textcontrol_dim3)
        self.baseline_correction_sizer_dim3.Add(self.baseline_correction_node_list_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_sizer_dim3.AddSpacer(10)
        # If polynomial baseline correction is selected, have a textcontrol for the polynomial order

        self.baseline_correction_polynomial_order_label = wx.StaticText(parent, -1, 'Polynomial order:')
        self.baseline_correction_sizer_dim3.Add(self.baseline_correction_polynomial_order_label, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_polynomial_order_textcontrol_dim3 = wx.TextCtrl(parent, -1, self.polynomial_order_dim3, size=(30, 20))
        self.baseline_correction_polynomial_order_textcontrol_dim3.Bind(wx.EVT_TEXT, self.on_baseline_correction_textcontrol_dim3)
        self.baseline_correction_sizer_dim3.Add(self.baseline_correction_polynomial_order_textcontrol_dim3, 0, wx.ALIGN_CENTER_VERTICAL)
        self.baseline_correction_sizer_dim3.AddSpacer(10)

        if(self.baseline_correction_radio_box_selection_dim3 == 0):
            self.baseline_correction_polynomial_order_label.Hide()
            self.baseline_correction_polynomial_order_textcontrol_dim3.Hide()


        # Have a button showing information on baseline correction
        self.baseline_correction_info = wx.Button(parent, -1, '\u24D8', size=(25, 32))


        self.baseline_correction_info.Bind(wx.EVT_BUTTON, self.oneDFrame.on_baseline_correction_info)
        self.baseline_correction_sizer_dim3.Add(self.baseline_correction_info, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer_2.Add(self.baseline_correction_sizer_dim3)
        self.sizer_2.AddSpacer(10)


    def on_baseline_correction_checkbox_dim3(self, event):
        self.baseline_correction_checkbox_value_dim3 = self.baseline_correction_checkbox_dim3.GetValue()



    def on_baseline_correction_radio_box_dim3(self, event):
        # If the user selects linear or polynomial baseline correction, update the options
        self.baseline_correction_radio_box_selection_dim3 = self.baseline_correction_radio_box_dim3.GetSelection()

        if(self.baseline_correction_radio_box_selection_dim3 == 0):
            # Remove the polynomial order textcontrol
            self.baseline_correction_sizer_dim3.Hide(self.baseline_correction_polynomial_order_label)
            self.baseline_correction_sizer_dim3.Hide(self.baseline_correction_polynomial_order_textcontrol_dim3)
            self.baseline_correction_sizer_dim3.Layout()
        elif(self.baseline_correction_radio_box_selection_dim3 == 1):
            # Add the polynomial order textcontrol
            self.baseline_correction_sizer_dim3.Show(self.baseline_correction_polynomial_order_label)
            self.baseline_correction_sizer_dim3.Show(self.baseline_correction_polynomial_order_textcontrol_dim3)
            self.baseline_correction_sizer_dim3.Layout()

    def on_baseline_correction_textcontrol_dim3(self, event):
        # If the node width or node list textcontrols are changed, update the node width and node list
        self.node_width_dim3 = self.baseline_correction_nodes_textcontrol_dim3.GetValue()
        self.node_list_dim3 = self.baseline_correction_node_list_textcontrol_dim3.GetValue()
        self.polynomial_order_dim3 = self.baseline_correction_polynomial_order_textcontrol_dim3.GetValue()


        
class InteractivePhasingFrame(wx.Frame):
    def __init__(self, main_frame, nmr_spectrum, ppms, nmr_d):
        # Get the monitor size and set the window size to 85% of the monitor size
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 1.0*self.monitorWidth
        self.height = 0.85*self.monitorHeight
        self.phasing_frame = wx.Frame.__init__(self, None, wx.ID_ANY,'Interactive Phasing',wx.DefaultPosition, size=(int(self.width), int(self.height)))

        self.main_frame = main_frame
        
        self.nmr_spectrum = nmr_spectrum
        self.ppms = ppms


        self.total_P0 = 0.0
        self.total_P1 = 0.0


        try:
            if(len(self.nmr_spectrum[0]) > 1):
                
                try:
                    len(self.nmr_spectrum[0][0])
                    self.nmr_spectrum = self.nmr_spectrum[0][0]
                except:
                    self.nmr_spectrum = self.nmr_spectrum[0]
        except:
            pass
        self.nmr_d = nmr_d


        self.create_canvas()

    
    def create_canvas(self):

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.sizer)
        self.fig = Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.toolbar = NavigationToolbar(self.canvas)

        if(darkdetect.isDark() == True and platform!='windows'):
            self.SetBackgroundColour((53, 53, 53, 255))
            self.toolbar.SetBackgroundColour((53, 53, 53, 255))
            self.fig.set_facecolor("#282A36")
            matplotlib.rc('axes',edgecolor='white')
            matplotlib.rc('xtick',color='white')
            matplotlib.rc('ytick',color='white')
            matplotlib.rc('axes',labelcolor='white')
            matplotlib.rc('axes',facecolor='#282A36')
            matplotlib.rc('legend', labelcolor='white')


        else:
            self.SetBackgroundColour('White')
            if(platform=='windows' and darkdetect.isDark() == True):
                self.toolbar.SetBackgroundColour('grey')
            else:
                self.toolbar.SetBackgroundColour('White')
            self.fig.set_facecolor("White")
            matplotlib.rc('axes',edgecolor='black')
            matplotlib.rc('xtick',color='black')
            matplotlib.rc('ytick',color='black')
            matplotlib.rc('axes',labelcolor='black')
            matplotlib.rc('axes',facecolor='white')
            matplotlib.rc('legend', labelcolor='black')

        
        self.sizer.Add(self.toolbar, 0, wx.EXPAND)
        self.sizer.AddSpacer(10)


        # Suppress complex warning from numpy 
        import warnings
        # warnings.simplefilter("ignore", np.ComplexWarning)  # For old numpy versions
        warnings.simplefilter("ignore", np.exceptions.ComplexWarning)   # For new numpy versions

        self.create_sliders()
        self.draw_figure_1D_phasing()
        self.Layout()
        self.Show()


    def create_sliders(self):
        from SpinExplorer.SpinView import FloatSlider
        # Create the phasing 1D sizer
        self.phasing_label = wx.StaticBox(self, -1, 'Phasing:')
        self.phasing_sizer = wx.StaticBoxSizer(self.phasing_label, wx.VERTICAL)
        self.P0_label = wx.StaticText(self, label="P0 (Coarse):")
        self.P1_label = wx.StaticText(self, label="P1 (Coarse):")
        self.P0_slider = FloatSlider(self, id=-1,value=0.0,minval=-180, maxval=180, res=0.1,size=(300, height))
        self.P1_slider = FloatSlider(self, id=-1,value=0.0,minval=-180, maxval=180, res=0.1,size=(300, height))
        self.P0_slider.Bind(wx.EVT_SLIDER, self.OnSliderScroll1D)
        self.P1_slider.Bind(wx.EVT_SLIDER, self.OnSliderScroll1D)
        self.P0_label_fine = wx.StaticText(self, label="P0 (Fine):     ")
        self.P1_label_fine = wx.StaticText(self, label="P1 (Fine):     ")
        self.P0_slider_fine = FloatSlider(self, id=-1,value=0.0,minval=-10, maxval=10, res=0.01,size=(300, height))
        self.P1_slider_fine = FloatSlider(self, id=-1,value=0.0,minval=-10, maxval=10, res=0.01,size=(300, height))
        self.P0_slider_fine.Bind(wx.EVT_SLIDER, self.OnSliderScroll1D)
        self.P1_slider_fine.Bind(wx.EVT_SLIDER, self.OnSliderScroll1D)
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
        self.phasing_combined.AddSpacer(160)
        self.phasing_combined.Add(self.P0_total_value)
        self.phasing_combined.AddSpacer(170)
        self.phasing_combined.Add(self.P1_total)
        self.phasing_combined.AddSpacer(160)
        self.phasing_combined.Add(self.P1_total_value)
        self.phasing_sizer.AddSpacer(5)
        self.phasing_sizer.Add(self.sizer_coarse)
        self.phasing_sizer.AddSpacer(10)
        self.phasing_sizer.Add(self.sizer_fine)
        self.phasing_sizer.AddSpacer(10)
        self.phasing_sizer.Add(self.phasing_combined)

        # Add a button to set the pivot point for phasing
        self.pivot_button = wx.Button(self, label="Set Pivot Point")
        self.pivot_button.Bind(wx.EVT_BUTTON, self.OnPivotButton)
        self.pivot_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.pivot_sizer.AddSpacer(500)
        self.pivot_sizer.Add(self.pivot_button)


        self.pivot_x_default = 0
        
        

        # Add a button to remove the pivot point
        self.remove_pivot_button = wx.Button(self, label="Remove Pivot Point")
        self.remove_pivot_button.Bind(wx.EVT_BUTTON, self.OnRemovePivotButton)
        self.pivot_sizer.AddSpacer(20)
        self.pivot_sizer.Add(self.remove_pivot_button)

        # Add the pivot point buttons to the phasing sizer
        self.phasing_sizer.AddSpacer(10)
        self.phasing_sizer.Add(self.pivot_sizer)



        # Create a sizer for changing the y axis limits in the spectrum
        self.zoom_label = wx.StaticBox(self, -1, 'Y Axis Zoom (%):')
        self.zoom_sizer = wx.StaticBoxSizer(self.zoom_label, wx.VERTICAL)
        self.intensity_slider = FloatSlider(self, id=-1,value=0,minval=-1, maxval=10, res=0.01,size=(300, height))
        self.intensity_slider.Bind(wx.EVT_SLIDER, self.OnIntensityScroll1D)
        self.zoom_sizer.AddSpacer(5)
        self.zoom_sizer.Add(self.intensity_slider)





        # Have a save and close button
        self.save_button = wx.Button(self, label="Save and Close")
        self.save_button.Bind(wx.EVT_BUTTON, self.OnSavePhasing)
        


        


        # Add all the sizers to the main sizer
        self.sizer1 = wx.BoxSizer(wx.VERTICAL)
        self.sizer2 = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer2.Add(self.phasing_sizer, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer2.AddSpacer(20)
        self.sizer2.Add(self.zoom_sizer, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer1.Add(self.sizer2, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.sizer1.AddSpacer(20)
        self.sizer1.Add(self.save_button, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.sizer1.AddSpacer(20)
        self.sizer.Add(self.sizer1, 0, wx.ALIGN_CENTER_HORIZONTAL)





    def UpdateFrame(self):
        self.canvas.draw()
        self.canvas.Refresh()
        self.canvas.Update()



    def OnSavePhasing(self,event):
        # Function to save the phasing values and close the window
        self.main_frame.phase_correction_p0_textcontrol.SetValue(str(self.total_P0))
        self.main_frame.phase_correction_p1_textcontrol.SetValue(str(self.total_P1))
        self.Close()
            

    def OnIntensityScroll1D(self, event):
        # Function to change the y axis limits
        intensity_percent = 10**float(self.intensity_slider.GetValue())
        self.ax.set_ylim(-(np.max(self.data)/8)/(intensity_percent/100),np.max(self.data)/(intensity_percent/100))
        self.UpdateFrame()



    def draw_figure_1D_phasing(self):
        # Function to plot the 1D spectrum
        self.ax = self.fig.add_subplot(111)

        self.data=self.nmr_spectrum
        self.line1, = self.ax.plot(self.ppms,self.data,linewidth=0.5)
        self.pivot_x = self.pivot_x_default
        self.ax.set_xlabel('Chemical shift (ppm)')
        self.ax.set_ylabel('Intensity')
        self.ax.set_xlim(max(self.ppms),min(self.ppms))
        self.line1.set_color('tab:blue')
        self.pivot_line = self.ax.axvline(self.pivot_x_default, color='black', linestyle='--')
        self.pivot_line.set_visible(False)
        self.UpdateFrame()


    def OnPivotButton(self,event):
        # Get the user to select a pivot point for phasing by clicking on the spectrum
        # Give a message box to tell the user to click on the spectrum where they want the pivot point to be
        wx.MessageBox('Click on the spectrum to set the location of the pivot point for P1 phasing.', 'Pivot Point', wx.OK | wx.ICON_INFORMATION)
        self.pivot_press = self.canvas.mpl_connect('button_press_event', self.OnPivotClick)
        


    
    def OnPivotClick(self,event):
        # Function to get the x value of the pivot point for phasing
        self.pivot_x = event.xdata
        self.pivot_line.set_xdata([self.pivot_x])
        
        # Find the index of the point closest to the pivot point
        self.pivot_index = np.abs(self.ppms-self.pivot_x).argmin()
        self.pivot_x = self.pivot_index
        self.canvas.mpl_disconnect(self.pivot_press)
        self.pivot_line.set_visible(True)
        self.OnSliderScroll1D(wx.EVT_SCROLL)


    def OnRemovePivotButton(self,event):
        if(self.pivot_line.get_visible()!=True):
            # Give a message saying there is no pivot point to remove
            wx.MessageBox('There is no pivot point to remove.', 'Remove Pivot Point', wx.OK | wx.ICON_INFORMATION)
        else:
            # Function to remove the pivot point for phasing
            self.pivot_x = self.pivot_x_default
            self.pivot_line.set_visible(False)
            self.OnSliderScroll1D(wx.EVT_SCROLL)
        

    def OnSliderScroll1D(self, event):
        #Get all the slider values for P0 and P1 (coarse and fine), put the combined coarse and fine values on the screen
        self.total_P0 = self.P0_slider.GetValue() + self.P0_slider_fine.GetValue()
        self.total_P1 = self.P1_slider.GetValue() + self.P1_slider_fine.GetValue()
        self.P0_total_value.SetLabel('{:.2f}'.format(self.total_P0))
        self.P1_total_value.SetLabel('{:.2f}'.format(self.total_P1))
        self.phase1D()

        
    def phase1D(self):
        # Function to phase the data using the combined course/fine phasing values and plot 
        imaginary_data = ng.process.proc_base.ht(self.nmr_spectrum, self.nmr_spectrum.shape[0])
        self.data = imaginary_data * np.exp(1j * (self.total_P0*np.pi/180 + self.total_P1*(np.pi/180) * (np.arange(-self.pivot_x, -self.pivot_x+self.nmr_spectrum.shape[0])/self.nmr_spectrum.shape[0]))) + np.ones(len(self.nmr_spectrum))
        self.line1.set_ydata(self.data + np.ones(len(self.data)))
        self.UpdateFrame()
        


    

def main():
    app = wx.App()
    frame = SpinProcess()
    app.MainLoop()

        



if __name__ == '__main__':
    main()
    