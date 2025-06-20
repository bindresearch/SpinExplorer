#!/usr/bin/env python3

"""MIT License

Copyright (c) 2025 James Eaton, Andrew Baldwin (University of Oxford)

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
print('                        SpinConverter                        ')
print('-------------------------------------------------------------')
print('                (version 1.2) 20th June 2025                 ')
print(' (c) 2025 James Eaton, Andrew Baldwin (University of Oxford) ')
print('                        MIT License                          ')
print('-------------------------------------------------------------')
print('            Converting NMR data to nmrPipe format            ')
print('-------------------------------------------------------------')
print(' Documentation at:')
print(' https://github.com/james-eaton-1/SpinExplorer')
print('-------------------------------------------------------------')
print('')



import sys

import wx

# Import relevant modules
import numpy as np
import matplotlib
matplotlib.use('WXAgg')
import nmrglue as ng
import subprocess
import os
import darkdetect
# Check to see if using mac, linux or windows
import sys
if(sys.platform == 'darwin'):
    platform = 'mac'
elif(sys.platform == 'linux'):
    platform = 'linux'
else:
    platform = 'windows'

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
# This program is designed to allow the user to convert NMR data from Bruker/Varian into NMRPipe format so it can be viewed using 
# SpinView.py. It is designed to be used with the SpinProcess.py program used to process the converted nmrPipe FID to produce 
# an NMR spectrum. These spectra can then be viewed using SpinView.py, a GUI for viewing NMR data.


# This class will be used for the conversion of Bruker/Varian data to NMRPipe format

class Converter():
    def __init__(self, parent):
        # Remove warning where FID is read using nmrglue
        import warnings
        warnings.simplefilter("ignore",UserWarning)

        self.parent = parent

        # Creating a hidden frame to be used as a parent for popout messages
        self.tempframe = wx.Frame(None, title="Temporary Parent", size=(1, 1))
        self.tempframe.Hide()  # Hide the frame since we don't need it to be visible


        # Get the NMR data and parameters
        self.find_nmr_files()
        self.read_nmr_data()
        self.find_parameters()

        


    def find_nmr_files(self):
        # Find the NMR files in the current directory
        self.parameter_file = ''
        if('acqus' in os.listdir('.')):
            self.spectrometer = 'Bruker'
            self.parameter_file = 'acqus'
        elif('procpar' in os.listdir('.')):
            self.spectrometer = 'Varian'
            self.parameter_file = 'procpar'
        else:
            # Try to find acqu file 
            if('acqu' in os.listdir('.')):
                self.spectrometer = 'Bruker'
                self.parameter_file = 'acqu'
                # Inform the user that an acqu parameter file but not an acqus parameter file has been found, would you like to continue?
                dlg = wx.MessageDialog(self.tempframe, 'An acqu parameter file but not an acqus parameter file has been found, would you like to continue?', 'Continue', wx.YES_NO | wx.ICON_ERROR)
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                if(dlg.ShowModal()==wx.ID_YES):
                    dlg.Destroy()
                else:
                    dlg.Destroy()
                    exit()
            else:
                # Give a popout error message saying that there are no bruker NMR files in the current directory, but found an acqus file
                dlg = wx.MessageDialog(self.tempframe, 'No Bruker (acqus) or Varian (procpar) NMR files found in the current directory. Please check the current directory and try again.', 'Error', wx.OK | wx.ICON_ERROR)
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                dlg.ShowModal()
                dlg.Destroy()
                exit()
        if(self.spectrometer == 'Bruker'):
            self.files = []
            for file in os.listdir('.'):
                if(file.endswith('ser') or file=='fid'):
                    self.files.append(file)
            if(len(self.files) == 0):
                # Give a popout error message saying that there are no bruker NMR files in the current directory, but found an acqus file
                dlg = wx.MessageDialog(self.tempframe, 'No Bruker NMR files found in the current directory. Please check the current directory and try again.', 'Error', wx.OK | wx.ICON_ERROR)
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                dlg.ShowModal()
                dlg.Destroy()
                exit()
                
            elif(len(self.files) == 2):
                # Give a popout error message saying that there are two NMR files in the current directory
                dlg = wx.MessageDialog(self.tempframe, 'An acqus file was found and two NMR files (fid and ser) found in the current directory. Please check the current directory and try again.', 'Error', wx.OK | wx.ICON_ERROR)
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                dlg.ShowModal()
                dlg.Destroy()
                exit()
        elif(self.spectrometer == 'Varian'):
            self.files = []
            for file in os.listdir('.'):
                if(file=='fid' or file=='origfid'):
                    self.files.append(file)
            if(len(self.files) == 0):
                # Give a popout error message saying that there are no bruker NMR files in the current directory, but found a procpar file
                dlg = wx.MessageDialog(self.tempframe, 'No Varian NMR files found in the current directory, but a procpar file was found. Please press no to exit, or press yes to manually select the raw FID file.', 'Error', wx.YES_NO | wx.ICON_ERROR)
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                if(dlg.ShowModal() == wx.ID_YES):
                    dlg = wx.FileDialog(self.tempframe, 'Select the raw FID file', wildcard="",style=wx.FD_OPEN)
                    dlg.SetDirectory(os.getcwd())


                    if(dlg.ShowModal() == wx.ID_OK):
                        self.files = [dlg.GetPath()]
                    else:
                        dlg.Destroy()
                        exit()
                else:
                    dlg.Destroy()
                    exit()
                



    def read_nmr_data(self):
        # Read the NMR data
        if(self.spectrometer == 'Bruker'):
            # if pdata directory exists, if it is empty then change its name to pdata_original
            if('pdata' in os.listdir('.')):
                if(os.listdir('pdata')==[]):
                    os.rename('pdata','pdata_original')
            try:
                self.nmr_dic, self.nmr_data = ng.bruker.read(dir='./', bin_file = self.files[0])
            except:
                dlg = wx.MessageDialog(self.tempframe,'Error: Unable to read the Bruker NMR data. Please check the current directory and try again.', 'Error', wx.OK | wx.ICON_ERROR)
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                dlg.ShowModal()
                dlg.Destroy()
                exit()
        if(self.spectrometer == 'Varian'):
            try:
                self.nmr_dic, self.nmr_data = ng.varian.read(dir='./', fid_file=self.files[0])
            except:
                dlg = wx.MessageDialog(self.tempframe,'Error: Unable to read the Varian NMR data. Please check the current directory and try again.', 'Error', wx.OK | wx.ICON_ERROR)
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                dlg.ShowModal()
                dlg.Destroy()
                exit()


        self.data_dimensions = len(self.nmr_data.shape)

                



        
    def find_parameters(self):
        # Find the parameters in the acqus/procpar files
        if(self.spectrometer == 'Bruker'):
            self.acqus_file = open(self.parameter_file,'r')
            self.acqus_file_lines = self.acqus_file.readlines()
            self.acqus_file.close()

            self.udic = ng.bruker.guess_udic(self.nmr_dic, self.nmr_data)
            self.find_size_bruker()
            self.find_pseudo_axis_bruker()
            self.find_sw_bruker()
            self.find_nucleus_frequencies_bruker()
            self.find_labels_bruker()
            self.find_temperature_bruker()
            self.calculate_carrier_frequency_bruker()
        
        elif(self.spectrometer == 'Varian'):
            self.procpar_file = open('procpar','r')
            self.procpar_file_lines = self.procpar_file.readlines()
            self.procpar_file.close()

            self.udic = ng.varian.guess_udic(self.nmr_dic, self.nmr_data) 
            self.find_size_varian()
            self.find_axes_pseudo_varian()
            self.find_sw_varian()
            self.find_nucleus_frequencies_varian()
            self.find_labels_varian()
            self.find_temperature_varian()
            self.calculate_carrier_frequency_varian()
            

            



            
            

    def find_size_bruker(self):
        # Look for TD in the acqus file
        for i in range(len(self.acqus_file_lines)):
            if('##$TD=' in self.acqus_file_lines[i]):
                line = self.acqus_file_lines[i].split()
                self.size_direct = int(line[1])
                self.size_direct_complex = int(self.size_direct)
                break
        
        # Look for the size of the indirect dimensions
        self.size_indirect = []
        found_indirect_sizes = False
        for i in range(len(self.acqus_file_lines)):
            if('##$TD_INDIRECT' in self.acqus_file_lines[i]):
                found_indirect_sizes = True
                continue
            if(found_indirect_sizes == True):
                line = self.acqus_file_lines[i].split()
                for j in range(len(line)):
                    if(line[j] == '0'):
                        pass
                    else:
                        self.size_indirect.append(int(line[j]))
                break


        try:
            self.size_direct
            # Find the max value in nmr_data.shape
            if(len(self.nmr_data.shape)==1):
                self.size_1 = max(self.nmr_data.shape)
                self.size_2 = self.size_1
                for i in range(len(self.size_indirect)):
                    self.size_2 = self.size_2 / self.size_indirect[i]
                self.size_direct_complex = int(self.size_2*2)
            else:
                # Sometimes have the issue where stored complex data size is larger than TD, this ensures that the direct dimension
                # size is altered to the larger size to correctly be read
                self.size_1 = self.nmr_data.shape[-1]
                self.size_2 = self.size_1
                if(self.size_1*2 > self.size_direct):
                    self.size_direct_complex = int(self.size_1*2)
                self.size_direct = self.size_1*2


        except:
            dlg = wx.MessageDialog(self.tempframe,'Error: TD not found in acqus file. Unable to determine size of data for direct dimension. Unable to convert data to NMRPipe format. Please check the acqus file and try again.', 'Error', wx.OK | wx.ICON_ERROR)
            self.tempframe.Raise()
            self.tempframe.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            exit()

        if(len(self.size_indirect) == 0):
            # Look to see if there is an acqu2s/acqu3s file
            if('acqu2s' in os.listdir('./')):
                self.sizes_dim2 = []
                self.sizes_dim2_nus = []
                # Look for TD in the acqu2s file
                with open('acqu2s','r') as file:
                    file_lines = file.readlines()
                    for i in range(len(file_lines)):
                        if('##$NUSTD=' in file_lines[i]):
                            line = file_lines[i].split('\n')[0].split()
                            self.sizes_dim2_nus = int(line[1])
                            break
                        if('##$TD=' in file_lines[i]):
                            line = file_lines[i].split()
                            self.sizes_dim2 = int(line[1])
                            break
                
                if(self.sizes_dim2_nus != []):
                    self.size_indirect.append(self.sizes_dim2_nus)
                else:
                    self.size_indirect.append(self.sizes_dim2)

            if('acqu3s' in os.listdir('./')):
                self.sizes_dim3 = []
                self.sizes_dim3_nus = []
                # Look for TD in the acqu3s file
                with open('acqu3s','r') as file:
                    file_lines = file.readlines()
                    for i in range(len(file_lines)):
                        if('##$NUSTD=' in file_lines[i]):
                            line = file_lines[i].split('\n')[0].split()
                            self.sizes_dim3_nus = int(line[1])
                            break
                        if('##$TD=' in file_lines[i]):
                            line = file_lines[i].split()
                            self.sizes_dim3 = int(line[1])
                            break
                
                if(self.sizes_dim3_nus != []):
                    self.size_indirect.append(self.sizes_dim3_nus)
                else:
                    self.size_indirect.append(self.sizes_dim3)
                    



        try:
            self.size_indirect
        except:
            dlg = wx.MessageDialog(self.tempframe,'Error: TD_INDIRECT not found in acqus file. Unable to determine size of data for indirect dimension. Unable to convert data to NMRPipe format. Please check the acqus file and try again.', 'Error', wx.OK | wx.ICON_ERROR)
            self.tempframe.Raise()
            self.tempframe.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            exit()
        
        if(self.size_indirect!=[] and 'acqu2s' not in os.listdir('./')):
            self.size_indirect = []

        # if(found_indirect_sizes==True):
        #     self.size_indirect.reverse()

        # Remove values from size indirect if they are equal to 1
        sizes_new = []
        for size in self.size_indirect:
            if(size!=1):
                sizes_new.append(size)

        self.size_indirect = sizes_new


        # Try to go through acqu2s and acqu3s and find the nucleus labels and corresponding TD values
        self.indirect_sizes_dict = {}
        if('acqu2s' in os.listdir('./')):
                # Look for TD in the acqu2s file
                nus = False
                with open('acqu2s','r') as file:
                    file_lines = file.readlines()
                    
                    for i in range(len(file_lines)):
                        if('##$NUSTD=' in file_lines[i] or '##$NusTD=' in file_lines[i]):
                            nus=True
                            line = file_lines[i].split('\n')[0].split()
                            self.sizes_dim2_nus = int(line[1])
                        if('##$TD=' in file_lines[i]):
                            line = file_lines[i].split()
                            self.sizes_dim2 = int(line[1])
                        if('##$NUC1=' in file_lines[i]):
                            line = file_lines[i].split('\n')[0]
                            nuc = line.split('<')[1].split('>')[0]
                    

                try:
                    if(nus==False):
                        self.indirect_sizes_dict[nuc] = self.sizes_dim2
                    else:
                        self.indirect_sizes_dict[nuc] = self.sizes_dim2_nus
                except:
                    pass
                

        if('acqu3s' in os.listdir('./')):
            # Look for TD in the acqu2s file
            nus = False
            with open('acqu3s','r') as file:
                file_lines = file.readlines()
                for i in range(len(file_lines)):
                    if('##$NUSTD=' in file_lines[i] or '##$NusTD=' in file_lines[i]):
                        nus=True
                        line = file_lines[i].split('\n')[0].split()
                        self.sizes_dim3_nus = int(line[1])
                    if('##$TD=' in file_lines[i]):
                        line = file_lines[i].split()
                        self.sizes_dim3 = int(line[1])
                    if('##$NUC1=' in file_lines[i]):
                        line = file_lines[i].split('\n')[0]
                        nuc = line.split('<')[1].split('>')[0]
                
            try:
                if(nus==False):
                    if(nuc not in self.indirect_sizes_dict.keys()):
                        self.indirect_sizes_dict[nuc] = self.sizes_dim3
                    else:
                        self.indirect_sizes_dict[nuc+'_1'] = self.sizes_dim3
                else:
                    if(nuc not in self.indirect_sizes_dict.keys()):
                        self.indirect_sizes_dict[nuc] = self.sizes_dim3_nus
                    else:
                        self.indirect_sizes_dict[nuc+'_1'] = self.sizes_dim3_nus
            except:
                pass






    def find_size_varian(self):
        # Look for np in the procpar file
        include_np_value = False
        include_ni_value = False
        include_ni2_value = False
        include_ni3_value = False
        self.size_direct = 0
        self.size_indirect = []
        for i in range(len(self.procpar_file_lines)):
                line = self.procpar_file_lines[i].split()
                if(include_np_value == True):
                    self.size_direct = int(line[1])
                    include_np_value = False
                elif(line[0] == 'np'):
                    include_np_value = True
                elif(include_ni_value == True):
                    self.size_indirect.append(int(line[1]))
                    include_ni_value = False
                elif(line[0] == 'ni'):
                    include_ni_value = True
                elif(include_ni2_value == True):
                    self.size_indirect.append(int(line[1]))
                    include_ni2_value = False
                elif(line[0] == 'ni2'):
                    include_ni2_value = True
                elif(include_ni3_value == True):
                    self.size_indirect.append(int(line[1]))
                    include_ni3_value = False
                elif(line[0] == 'ni3'):
                    include_ni3_value = True
                    

        

    
        if(self.size_direct == 0):
            dlg = wx.MessageDialog(self.tempframe,'Error: np not found in procpar file. Unable to determine size of data for direct dimension. Unable to convert data to NMRPipe format. Please check the procpar file and try again.', 'Error', wx.OK | wx.ICON_ERROR)
            self.tempframe.Raise()
            self.tempframe.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            exit()


        
        

    def find_sw_varian(self):
        # Look for sw, sw1, sw2 in the procpar file
        include_sw_value = False
        include_sw1_value = False
        include_sw2_value = False
        self.sw_direct = 0
        self.sw_indirect = {}
        for i in range(len(self.procpar_file_lines)):
            line = self.procpar_file_lines[i].split()
            if(include_sw_value == True):
                self.sw_direct = float(line[1])
                include_sw_value = False
            elif(line[0] == 'sw'):
                include_sw_value = True
            elif(include_sw1_value == True):
                self.sw_indirect['sw1']=float(line[1])
                include_sw1_value = False
            elif(line[0] == 'sw1'):
                include_sw1_value = True
            elif(include_sw2_value == True):
                self.sw_indirect['sw2']=float(line[1])
                include_sw2_value = False
            elif(line[0] == 'sw2'):
                include_sw2_value = True


        


    def find_nucleus_frequencies_varian(self):
        # Look for sfrq, dfrq2, dfrq3 in the procpar file
        include_sfrq_value = False
        include_dfrq_value = False
        include_dfrq2_value = False
        include_dfrq3_value = False
        self.nucleus_frequency_direct = 0
        self.nucleus_frequencies_indirect = []
        self.nucleus_frequencies_indirect_order = []
        include_zero = False
        for i in range(len(self.procpar_file_lines)):
            line = self.procpar_file_lines[i].split()
            if(include_sfrq_value == True):
                self.nucleus_frequency_direct = float(line[1])
                include_sfrq_value = False
            elif(line[0] == 'sfrq'):
                include_sfrq_value = True
            elif(include_dfrq_value == True):
                if(line[1]!='0'):
                    self.nucleus_frequencies_indirect.append(float(line[1]))
                    self.nucleus_frequencies_indirect_order.append('dfrq')
                else:
                    include_zero = True
                include_dfrq_value = False
            elif(line[0] == 'dfrq'):
                include_dfrq_value = True
            elif(include_dfrq2_value == True):
                if(line[1]!='0'):
                    self.nucleus_frequencies_indirect.append(float(line[1]))
                    self.nucleus_frequencies_indirect_order.append('dfrq2')
                else:
                    include_zero = True
                include_dfrq2_value = False
            elif(line[0] == 'dfrq2'):
                include_dfrq2_value = True
            elif(include_dfrq3_value == True):
                if(line[1]!='0'):
                    self.nucleus_frequencies_indirect.append(float(line[1]))
                    self.nucleus_frequencies_indirect_order.append('dfrq3')
                else:
                    include_zero = True
                include_dfrq3_value = False
            elif(line[0] == 'dfrq3'):
                include_dfrq3_value = True

        
        if(self.nucleus_frequency_direct == 0):
            dlg = wx.MessageDialog(self.tempframe,'Error: sfrq not found in procpar file. Unable to determine nucleus frequency for direct dimension. Unable to convert data to NMRPipe format. Please check the procpar file and try again.', 'Error', wx.OK | wx.ICON_ERROR)
            self.tempframe.Raise()
            self.tempframe.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            exit()
        
        if(self.nucleus_frequencies_indirect == []):
            dlg = wx.MessageDialog(self.tempframe,'Error: dfrq/dfrq2/dfrq3 not found in procpar file. Would you like to continue anyway?', 'Error', wx.YES_NO | wx.ICON_ERROR)
            self.tempframe.Raise()
            self.tempframe.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            if(dlg.GetReturnCode() == wx.ID_NO):
                exit()
            else:
                self.nucleus_frequencies_indirect = [0.0,0.0,0.0]
                self.nucleus_frequencies_indirect_order = ['dfrq','dfrq2','dfrq3']

        
        self.nucleus_frequencies_indirect.reverse()
        if(include_zero == True):
            self.nucleus_frequencies_indirect.append(0.0)
        self.nucleus_frequencies = [self.nucleus_frequency_direct] + self.nucleus_frequencies_indirect + [0.0]
        
    
    def find_labels_varian(self):
        # Find tn, dn, dn2, dn3 in the procpar file
        include_tn_value = False
        include_dn_value = False
        include_dn2_value = False
        include_dn3_value = False
        self.label_direct = ''
        self.labels_indirect = []
        for i in range(len(self.procpar_file_lines)):
            line = self.procpar_file_lines[i].split()
            if(include_tn_value == True):
                self.label_direct = line[1].split('"')[1]
                include_tn_value = False
            elif(line[0] == 'tn'):
                include_tn_value = True
            elif(include_dn_value == True):
                self.labels_indirect.append(line[1].split('"')[1])
                include_dn_value = False
            elif(line[0] == 'dn'):
                include_dn_value = True
            elif(include_dn2_value == True):
                self.labels_indirect.append(line[1].split('"')[1])
                include_dn2_value = False
            elif(line[0] == 'dn2'):
                include_dn2_value = True
            elif(include_dn3_value == True):
                self.labels_indirect.append(line[1].split('"')[1])
                include_dn3_value = False
            elif(line[0] == 'dn3'):
                include_dn3_value = True

        
        if(self.label_direct == ''):
            dlg = wx.MessageDialog(self.tempframe,'Error: Label for direct dimension (tn) not found in procpar file. Setting label as 1.', 'Error', wx.OK | wx.ICON_ERROR)
            self.tempframe.Raise()
            self.tempframe.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            self.label_direct = '1'
        
        if(self.labels_indirect == []):
            dlg = wx.MessageDialog(self.tempframe,'Error: Labels for indirect dimensions (dn/dn2/dn3) not found in procpar file. Setting labels as 2, 3, 4.', 'Error', wx.OK | wx.ICON_ERROR)
            self.tempframe.Raise()
            self.tempframe.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            self.labels_indirect = ['2','3','4']

        if(self.other_params == True):
            self.labels_correct_order = [self.labels_indirect[0]] + [self.arrayed_parameter] + self.labels_indirect[1:]
        else:
            self.labels_correct_order = self.labels_indirect 
            if(self.reverse_acquisition_order == True):
                self.labels_correct_order.reverse()



    def find_temperature_varian(self):
        # Find the temperature the spectrum was recorded at
        self.temperature = 298.15   # Default temperature is 298.15K
        include_temp_value = False
        for i in range(len(self.procpar_file_lines)):
            line = self.procpar_file_lines[i].split()
            if(include_temp_value == True):
                self.temperature = float(line[1]) + 273.15
                include_temp_value = False
            elif(line[0] == 'temp'):
                include_temp_value = True

    
    def calculate_carrier_frequency_varian(self):
        # Calculate the carrier frequency for each dimension
        # For direct dimension, if proton calculate based on water and include water referencing as an option
        if(self.label_direct == '1H' or self.label_direct == 'H1' or self.label_direct == 'H'):
            #return water chemical shift in range 0-100oC
            self.water_ppm = 7.83 - self.temperature/96.9
            
            


            self.references_proton = [self.water_ppm, 0]
            self.references_proton_labels = ['H2O', 'Manual']
            
            self.references_other = []
            self.references_other_labels = []

            # Determine frequency of other nuclei at 0ppm based on water ppm
            self.sfrq0 = self.nucleus_frequency_direct/(1+self.water_ppm*1e-6)
            self.dfrq_13C = self.sfrq0*0.251449530 
            self.dfrq_15N = self.sfrq0*0.101329118
            self.dfrq_P31 = self.sfrq0*0.4048064954
            self.dfrq_F19 = self.sfrq0*0.941

            # Use dfrq2 and dfrq3 to calculate chemical shifts of the indirect dimensions
            
            # Determine order of nucleus labels based on dfrq2/dfrq3
            self.ordered_indirect_labels = []
            self.sw_ordered = []
            found_match = False
            for i, value in enumerate(self.nucleus_frequencies_indirect):
                if(np.abs(value-self.dfrq_13C)/self.dfrq_13C * 100 < 1):
                    self.ordered_indirect_labels.append('C13')
                    self.references_other.append((value-self.dfrq_13C)/self.dfrq_13C * 1E6)
                    self.references_other_labels.append('C13 (Referenced to H2O)')
                    found_match = True
                elif(np.abs(value-self.dfrq_15N)/self.dfrq_15N * 100 < 1):
                    self.ordered_indirect_labels.append('N15')
                    self.references_other.append((value-self.dfrq_15N)/self.dfrq_15N * 1E6)
                    self.references_other_labels.append('N15 (Referenced to H2O)')
                    found_match = True
                elif(np.abs(value-self.dfrq_P31)/self.dfrq_P31 * 100 < 1):
                    self.ordered_indirect_labels.append('P31')
                    self.references_other.append((value-self.dfrq_P31)/self.dfrq_P31 * 1E6)
                    self.references_other_labels.append('P31 (Referenced to H2O)')
                    found_match = True
                elif(np.abs(value-self.dfrq_F19)/self.dfrq_F19 * 100 < 1):
                    self.ordered_indirect_labels.append('F19')
                    self.references_other.append((value-self.dfrq_F19)/self.dfrq_F19 * 1E6)
                    self.references_other_labels.append('F19 (Referenced to H2O)')
                    found_match = True
                elif(np.abs(value-self.sfrq0)/self.sfrq0 * 100 < 1):
                    self.ordered_indirect_labels.append('H1')
                    self.references_other.append(self.water_ppm)
                    self.references_other_labels.append('H2O')
                    found_match = True
                if(found_match==True):
                    # self.sw_ordered.append(self.sw_indirect[i])
                    found_match = False
        
        else:
            self.references_proton = [0.0]
            self.references_proton_labels = ['Manual']

            self.references_other = [0.0]
            self.references_other_labels = ['Manual']






    
    def find_axes_pseudo_varian(self):
        # Search through procpar to find array values 
        self.phase = False
        self.phase2 = False
        self.other_params = False

        self.found_array = False
        for i, line in enumerate(self.procpar_file_lines):
            if(self.found_array==True):
                array = line.split()[1]
                self.found_array = False
            if(line.split()[0] == 'array'):
                self.found_array = True

        array = array.split('"')[1].split(',')
        if(array==['']):
            array = []
        
        self.reorder_array = False
        count=0
        if('phase' in array):
            self.phase = True
            count+=1
        if('phase2' in array):
            self.phase2 = True
            count+=1
        if(count!=len(array)):
            self.other_params = True

        # Delete phase and phase2 from array
        self.found_arrayed_parameter = False
        self.number_of_arrayed_parameters = 0
        self.phases = []
        if(self.other_params==True):
            if(self.phase==True):
                array.remove('phase')
                self.phases.append('phase')
            if(self.phase2==True):
                array.remove('phase2')
                self.phases.append('phase2')
            self.arrayed_parameter = array[0]


        
            # Search through procpar file to find the number of values of the arrayed parameter
            
            for i, line in enumerate(self.procpar_file_lines):
                if(self.found_arrayed_parameter==True):
                    self.number_of_arrayed_parameters = int(line.split()[0])
                    self.found_arrayed_parameter = False
                if(line.split()[0] == self.arrayed_parameter):
                    self.found_arrayed_parameter = True

            
        for i, line in enumerate(self.procpar_file_lines):
            if(self.found_array==True):
                array = line.split()[1]
                self.found_array = False
            if(line.split()[0] == 'array'):
                self.found_array = True

        array = array.split('"')[1].split(',')
        if(array==['']):
            array = []
        
        self.reorder_array = False
        count=0
        if('phase' in array and 'phase2' in array):
            self.phase = True
            count+=1
            if(array[-1]=='phase' or array[-1]=='phase2' and len(array)>2):
                self.reorder_array = True
        elif('phase2' in array):
            self.phase2 = True
            count+=1
            if(array[-1]=='phase2' and len(array)>1):
                self.reorder_array = True
        elif('phase' in array):
            self.phase = True
            count+=1
            if(array[-1]=='phase' and len(array)>1):
                self.reorder_array = True
        self.reverse_acquisition_order = False
        if('phase' in array and 'phase2' in array):
            if(array==['phase','phase2']):
                self.reverse_acquisition_order = True
            else:
                self.reverse_acquisition_order = False

        if(self.reverse_acquisition_order == False):
            self.size_indirect.reverse()


        if(0 in self.size_indirect):
            self.size_indirect.remove(0)
            self.size_indirect.append(0)
           
      



        
            

        

        


    def find_sw_bruker(self):
        for i in range(len(self.size_indirect)+1):
            if(i==0):
                for j in range(len(self.acqus_file_lines)):
                    if('##$SW_h=' in self.acqus_file_lines[j]):
                        line = self.acqus_file_lines[j].split()
                        self.sw_direct = float(line[1])
                        break
            if(i==1):
                try:
                    self.sw_indirect = []
                    file = open('acqu2s','r')
                    file_lines = file.readlines()
                    file.close()
                    for j in range(len(file_lines)):
                        if('##$SW_h=' in file_lines[j]):
                            line = file_lines[j].split()
                            self.sw_indirect.append(float(line[1]))
                            break
                except:
                    self.sw_indirect.append(0)
            if(i==2):
                try:
                    file = open('acqu3s','r')
                    file_lines = file.readlines()
                    file.close()
                    for j in range(len(file_lines)):
                        if('##$SW_h=' in file_lines[j]):
                            line = file_lines[j].split()
                            self.sw_indirect.append(float(line[1]))
                            break
                except:
                    self.sw_indirect.append(0)
            if(i==3):
                try:
                    file = open('acqu4s','r')
                    file_lines = file.readlines()
                    file.close()
                    for j in range(len(file_lines)):
                        if('##$SW_h=' in file_lines[j]):
                            line = file_lines[j].split()
                            self.sw_indirect.append(float(line[1]))
                            break
                except:
                    self.sw_indirect.append(0)

            if(i>3):
                dlg = wx.MessageDialog(self.tempframe,'Error: Only able to convert data with up to 4 indirect dimensions. Unable to convert data to NMRPipe format. Please check the acqus file and try again.', 'Error', wx.OK | wx.ICON_ERROR)
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                dlg.ShowModal()
                dlg.Destroy()
                exit()


    def find_nucleus_frequencies_bruker(self):
        # Find the nucleus frequencies in the acqus file
        self.nucleus_frequencies = []
        for i in range(len(self.size_indirect)+1):
            if(i==0):
                for j in range(len(self.acqus_file_lines)):
                    if('##$SFO1=' in self.acqus_file_lines[j]):
                        line = self.acqus_file_lines[j].split()
                        self.nucleus_frequencies.append(float(line[1]))
                        break
            if(i==1):
                try:
                    file = open('acqu2s','r')
                    file_lines = file.readlines()
                    file.close()
                    for j in range(len(file_lines)):
                        if('##$SFO1=' in file_lines[j]):
                            line = file_lines[j].split()
                            self.nucleus_frequencies.append(float(line[1]))
                            break
                except:
                    self.nucleus_frequencies.append(0)
            if(i==2):
                try:
                    file = open('acqu3s','r')
                    file_lines = file.readlines()
                    file.close()
                    for j in range(len(file_lines)):
                        if('##$SFO1=' in file_lines[j]):
                            line = file_lines[j].split()
                            self.nucleus_frequencies.append(float(line[1]))
                            break
                except:
                    self.nucleus_frequencies.append(0)
            if(i==3):
                try:
                    file = open('acqu4s','r')
                    file_lines = file.readlines()
                    file.close()
                    for j in range(len(file_lines)):
                        if('##$SFO1=' in file_lines[j]):
                            line = file_lines[j].split()
                            self.nucleus_frequencies.append(float(line[1]))
                            break
                except:
                    self.nucleus_frequencies.append(0)



    def find_labels_bruker(self):
        # Find the labels in the acqus file
        self.labels = []
        for j in range(len(self.acqus_file_lines)):
            if('##$NUC' in self.acqus_file_lines[j] and '##$NUCLEUS' not in self.acqus_file_lines[j] and '##$NUCLEI' not in self.acqus_file_lines[j]):
                line = self.acqus_file_lines[j].split()
                self.labels.append(line[1].split('<')[1].split('>')[0])
  


        # If the labels are off, then change them to ID
        for i in range(len(self.labels)):
            if(self.labels[i] == 'off'):
                self.labels[i] = 'ID'


        # Work out what the correct order of labels is (i.e. if the nuclei match their frequencies etc)
        self.find_gamma_bruker()
        self.nuc1 = self.labels[0]
        # Find out the field strength using gamma and the carrier frequency
        field = self.nucleus_frequencies[0]/self.gamma[self.nuc1]
        # Find the correct order of labels
        self.labels_correct_order = []
        self.labels_correct_order.append(self.nuc1)
        for i in range(len(self.nucleus_frequencies)):
            if(i==0):
                continue
            else:
                if(self.nucleus_frequencies[i]==0):
                    self.labels_correct_order.append('ID')
                else:
                    for key in self.gamma:
                        if(abs(abs(field*self.gamma[key])-self.nucleus_frequencies[i])<1):
                            if(key in self.labels):
                                if(key in self.labels_correct_order):
                                    self.labels_correct_order.append(key+'_1')
                                else:
                                    self.labels_correct_order.append(key)
            


    

    def find_gamma_bruker(self):
        # Write out gyromagnetic ratio values known for nuclei (1H, 19F, 13C, 15N, 31P, 23Na, 2H) etc
        self.gamma = {}
        self.gamma['1H'] = 267.5153151E6   #267.5221877E6
        self.gamma['19F'] = 251.6628277E6
        self.gamma['13C'] = 67.262E6
        self.gamma['14N'] = 1.93297E7
        self.gamma['15N'] = -27.116E6
        self.gamma['31P'] = 108.282E6
        self.gamma['23Na'] = 70.882E6
        self.gamma['25Mg'] = -1.639E7
        self.gamma['39K'] = 1.2498E7
        self.gamma['41K'] = 0.686E7
        self.gamma['43Ca'] = -1.8025E7
        self.gamma['2H'] = 41.065E6
        self.gamma['7Li'] = 103.962E6
        self.gamma['17O'] = -36.264E6
        self.gamma['10B'] = 2.87471E7
        self.gamma['11B'] = 8.58406E7
        self.gamma['27Al'] = 6.97594E7
        self.gamma['29Si'] = -5.3146E7
        self.gamma['35Cl'] = 2.62401E7
        self.gamma['37Cl'] = 2.18428E7
        self.gamma['50V'] = 2.67164E7
        self.gamma['51V'] = 7.04578E7
        self.gamma['55Mn'] = 6.59777E7
        self.gamma['57Fe'] = 0.86399E7
        self.gamma['59Co'] = 6.3472E7
        self.gamma['63Cu'] = 7.0965E7
        self.gamma['65Cu'] = 7.6018E7
        self.gamma['67Zn'] = 16.767E6 
        self.gamma['69Ga'] = 6.43685E7
        self.gamma['71Ga'] = 8.180163E7
        self.gamma['77Se'] = 5.115E7
        self.gamma['79Br'] = 6.70186E7
        self.gamma['81Br'] = 7.22421E7
        self.gamma['103Rh'] = -0.84579E7
        self.gamma['107Ag'] = -1.08718E7
        self.gamma['109Ag'] = -1.25001E7
        self.gamma['111Cd'] = -5.69259E7
        self.gamma['113Cd'] = -5.95504E7
        self.gamma['117Sn'] = -9.57865E7
        self.gamma['119Sn'] = -10.01926
        self.gamma['123Te'] = -7.04893E7
        self.gamma['125Te'] = -8.49722E7
        self.gamma['127I'] = 5.37937E7
        self.gamma['129Xe'] = -7.44069E7
        self.gamma['131Xe'] = 2.20564E7
        self.gamma['183W'] = 1.12070E7
        self.gamma['195Pt'] = 5.80466E7
        self.gamma['197Au'] = 0.4692E7
        self.gamma['199Hg'] = 4.81519E6
        self.gamma['201Hg'] = -1.77748E7
        self.gamma['203Tl'] = 15.43599E7
        self.gamma['205Tl'] = 15.58829E7
        self.gamma['207Pb'] = 5.64661E7
        



    


                



        





    def find_pseudo_axis_bruker(self):   
        # Read pulseprogram.precomp (or pulseprogram) file to see if the data has a pseudo axis. Will have a QF in the acqusition mode of that dimension
        self.acqusition_modes = []
        if(self.data_dimensions>1 or self.size_indirect!=[]):
            self.pseudo_flag = 0
            try:
                self.pulseprogram_file = open('pulseprogram.precomp','r')
                self.pulseprogram_file_lines = self.pulseprogram_file.readlines()
                self.pulseprogram_file.close()
            except:
                try:
                    self.pulseprogram_file = open('pulseprogram','r')
                    self.pulseprogram_file_lines = self.pulseprogram_file.readlines()
                    self.pulseprogram_file.close()
                except:
                    pass
            
            count = 0
            try:
                
                for i in range(len(self.pulseprogram_file_lines)):
                    if('AQ_mode' in self.pulseprogram_file_lines[i]):
                        line = self.pulseprogram_file_lines[i].split()
                        # Find the index of AQ_mode and delete everything before it
                        for j in range(len(line)):
                            if(line[j] == 'AQ_mode'):
                                index = j
                        line = line[index+1:]
                        # Remove any terms containing brackets
                        new_line = []
                        for j in range(len(line)):
                            if('(' not in line[j]):
                                new_line.append(line[j])
                        break
            
                if(new_line.count('QF')>0):
                    # Then at least one of the dimensions is a pseudo axis
                    if(new_line.count('QF') == 1):
                        # Then there is only one pseudo axis
                        self.pseudo_flag = 1
                    else:
                        # Then there are multiple pseudo axes
                        for i in range(len(new_line)):
                            if(new_line[i] == 'QF'):
                                count += 1
                        self.pseudo_flag = count


                self.acqusition_modes = new_line
                

            except:
                try:
                    self.pseudo_flag = 0
                    # Read pulseprogram file to see if there is an FnMode option to input the acqusition mode
                    for i in range(len(self.pulseprogram_file_lines)):
                        if('FnMODE' in self.pulseprogram_file_lines[i]):
                            line = self.pulseprogram_file_lines[i].split('\n')[0].split('FnMODE:')[1].strip()
                            if(len(line.split())>1):
                                self.acqusition_modes = line.split()
                            elif(line.lower() == 'echo-antiecho'):
                                self.acqusition_modes = ['Echo-Antiecho']
                            elif(line.lower() == 'complex'):
                                self.acqusition_modes = ['Complex']
                            elif(line.lower() == 'tppi'):
                                self.acqusition_modes = ['TPPI']
                            elif(line.lower() == 'qf'):
                                self.acqusition_modes = ['QF']

                            
                            
            
                except:
                    try:
                        self.pseudo_flag = 0
                        self.acqusition_modes = []
                        for line in self.pulseprogram_file_lines:
                            if('define list<frequency>' in line or 'define list<delay>' in line):
                                self.pseudo_flag = 1
                                self.acqusition_modes = ['QF']

                                break

                    
                    except:
                        self.pseudo_flag = 0
                        self.acqusition_modes = ['Complex']

            
            while(len(self.acqusition_modes) < len(self.size_indirect)):
                for i in range(len(self.size_indirect)-len(self.acqusition_modes)):
                    self.acqusition_modes.append('Complex')


            if(len(self.size_indirect) == 1):
                if(self.pseudo_flag==0):
                    if(self.size_indirect[0]==1):
                        self.pseudo_flag = 1
                        self.acqusition_modes[-1] = 'QF'
            elif(len(self.size_indirect) == 2):
                if(self.pseudo_flag==0):
                    if(self.size_indirect[0]==1):
                        self.pseudo_flag += 1
                        self.acqusition_modes[-2] = 'QF'
                    if(self.size_indirect[1]==1):
                        self.pseudo_flag += 1
                        self.acqusition_modes[-1] = 'QF'

            

            # If there is a pseudo axis, then the number of real points in that dimension is the same as the number of complex points
            if(self.pseudo_flag == 1):
                sizes = []
                keys = []
                for key in self.udic:
                    if(key == 'ndim'):
                        continue
                    else:
                        keys.append(key)
                        sizes.append(self.udic[key]['size'])
                
                # Find the index of the minimum size (this is default set to the pseudo axis)
                index = sizes.index(min(sizes))
                self.udic[keys[index]]['complex'] = False
                self.udic[keys[index]]['encoding'] = 'Real'
                self.udic[keys[index]]['obs'] = 0
                self.udic[keys[index]]['sw'] = 0
                self.udic[keys[index]]['car'] = 0


            elif(self.pseudo_flag == 2):
                sizes = []
                keys = []
                for key in self.udic:
                    if(key == 'ndim'):
                        continue
                    else:
                        keys.append(key)
                        sizes.append(self.udic[key]['size'])
                
                # Find the index of the minimum size (this is default set to the pseudo axis)
                index = sizes.index(min(sizes))
                self.udic[keys[index]]['encoding'] = 'Real'
                self.udic[keys[index]]['complex'] = False
                self.udic[keys[index]]['obs'] = 1
                self.udic[keys[index]]['sw'] = 1
                self.udic[keys[index]]['car'] = 1
                # Then find the next smallest size and set this to the pseudo axis
                sizes[index] = max(sizes)
                index = sizes.index(min(sizes))
                self.udic[keys[index]]['encoding'] = 'Real'
                self.udic[keys[index]]['complex'] = False
                self.udic[keys[index]]['obs'] = 1
                self.udic[keys[index]]['sw'] = 1
                self.udic[keys[index]]['car'] = 1 

        



        
    def find_temperature_bruker(self):
        # Find the temperature the spectrum was recorded at
        self.temperature = 298.15 # Default temperature is 298.15K
        for i in range(len(self.acqus_file_lines)):
            if('##$TE=' in self.acqus_file_lines[i]):
                line = self.acqus_file_lines[i].split()
                self.temperature = float(line[1])
                break

    def calculate_carrier_frequency_bruker(self):
        # Calculate the carrier frequency for each dimension
        # For direct dimension, if proton calculate based on water and include water referencing as an option

        if(len(self.size_indirect) == 0):
            #return water chemical shift in range 0-100oC
            if(self.labels_correct_order[0] == '1H' or self.labels_correct_order[0] == 'H1' or self.labels_correct_order[0] == 'H'):
                self.water_ppm = 7.83 - self.temperature/96.9
                
                # Use O1/BF1 to calculate a second carrier frequency in case not centred on water
                for j in range(len(self.acqus_file_lines)):
                    if('##$O1=' in self.acqus_file_lines[j]):
                        line = self.acqus_file_lines[j].split()
                        self.O1 = float(line[1])
                        break
                for j in range(len(self.acqus_file_lines)):
                    if('##$BF1=' in self.acqus_file_lines[j]):
                        line = self.acqus_file_lines[j].split()
                        self.BF1 = float(line[1])
                        break
                self.carrier_frequency_1 = self.O1/self.BF1


                self.references_proton = [self.water_ppm, self.carrier_frequency_1]
                self.references_proton_labels = ['H2O', 'O1/BF1']
            else:
                # Use O1/BF1 to calculate a second carrier frequency in case not centred on water
                for j in range(len(self.acqus_file_lines)):
                    if('##$O1=' in self.acqus_file_lines[j]):
                        line = self.acqus_file_lines[j].split()
                        self.O1 = float(line[1])
                        break
                for j in range(len(self.acqus_file_lines)):
                    if('##$BF1=' in self.acqus_file_lines[j]):
                        line = self.acqus_file_lines[j].split()
                        self.BF1 = float(line[1])
                        break
                self.carrier_frequency_1 = self.O1/self.BF1


                self.references_proton = [self.carrier_frequency_1]
                self.references_proton_labels = ['O1/BF1']


                    
   
        else:
            if(self.labels_correct_order[0] == '1H' or self.labels_correct_order[0] == 'H1' or self.labels_correct_order[0] == 'H'):
                #return water chemical shift in range 0-100oC
                self.water_ppm = 7.83 - self.temperature/96.9
            
            # Use O1/BF1 to calculate a second carrier frequency in case not centred on water
            for j in range(len(self.acqus_file_lines)):
                if('##$O1=' in self.acqus_file_lines[j]):
                    line = self.acqus_file_lines[j].split()
                    self.O1 = float(line[1])
                    break
            for j in range(len(self.acqus_file_lines)):
                if('##$BF1=' in self.acqus_file_lines[j]):
                    line = self.acqus_file_lines[j].split()
                    self.BF1 = float(line[1])
                    break
            self.carrier_frequency_1 = self.O1/self.BF1

            # Calculate carrier frequencies based on O2/BF2, O3/BF3

            for j in range(len(self.acqus_file_lines)):
                if('##$O2=' in self.acqus_file_lines[j]):
                    line = self.acqus_file_lines[j].split()
                    self.O2 = float(line[1])
                    break
            for j in range(len(self.acqus_file_lines)):
                if('##$O3=' in self.acqus_file_lines[j]):
                    line = self.acqus_file_lines[j].split()
                    self.O3 = float(line[1])
                    break
            for j in range(len(self.acqus_file_lines)):
                if('##$BF2=' in self.acqus_file_lines[j]):
                    line = self.acqus_file_lines[j].split()
                    self.BF2 = float(line[1])
                    break
            for j in range(len(self.acqus_file_lines)):
                if('##$BF3=' in self.acqus_file_lines[j]):
                    line = self.acqus_file_lines[j].split()
                    self.BF3 = float(line[1])
                    break
            # Calculate the carrier frequency
            self.carrier_frequency_2 = self.O2/self.BF2
            self.carrier_frequency_3 = self.O3/self.BF3

            self.ppms_referenced = []
            self.ppms_referenced_labels = []

            if(self.labels_correct_order[0] == '1H' or self.labels_correct_order[0] == 'H1' or self.labels_correct_order[0] == 'H'):
                # Calculate referenced carrier frequencies based on water chemical shifts
                
                self.sfrq0 = self.nucleus_frequencies[0]/(1+self.water_ppm*1e-6)
                self.dfrq_13C = self.sfrq0*0.251449530 
                self.dfrq_15N = self.sfrq0*0.101329118
                self.dfrq_P31 = self.sfrq0*0.4048064954
                self.dfrq_F19 = self.sfrq0*0.9412866605363297

                

                for i, label in enumerate(self.labels_correct_order):
                    # if(i+1>len(self.nucleus_frequencies)):
                    #     break
                    if(label=='15N' or label=='N15'):
                        self.ppms_referenced.append(((self.nucleus_frequencies[i]-self.dfrq_15N)/self.dfrq_15N)*1e6)
                        self.ppms_referenced_labels.append('15N (Referenced to H2O)')
                    if(label=='13C' or label=='C13'):
                        self.ppms_referenced.append(((self.nucleus_frequencies[i]-self.dfrq_13C)/self.dfrq_13C)*1e6)
                        self.ppms_referenced_labels.append('13C (Referenced to H2O)')
                    if(label=='31P' or label=='P31'):
                        self.ppms_referenced.append(((self.nucleus_frequencies[i]-self.dfrq_P31)/self.dfrq_P31)*1e6)
                        self.ppms_referenced_labels.append('31P (Referenced to H2O)')
                    if(label=='19F' or label=='F19'):
                        self.ppms_referenced.append(((self.nucleus_frequencies[i]-self.dfrq_F19)/self.dfrq_F19)*1e6)
                        self.ppms_referenced_labels.append('19F (Referenced to H2O)')






            if(self.labels_correct_order[0]=='1H' or self.labels_correct_order[0]=='H1' or self.labels_correct_order[0]=='H'):
                self.references_proton = [self.water_ppm, self.carrier_frequency_1, self.carrier_frequency_2, self.carrier_frequency_3]
                self.references_proton_labels = ['H2O', 'O1/BF1', 'O2/BF2', 'O3/BF3']
            else:
                self.references_proton = [self.carrier_frequency_1, self.carrier_frequency_2, self.carrier_frequency_3]
                self.references_proton_labels = ['O1/BF1', 'O2/BF2', 'O3/BF3']

            self.references_other = []
            self.references_other_labels = []
            for i, ppm in enumerate(self.ppms_referenced):
                self.references_other.append(ppm)
                self.references_other_labels.append(self.ppms_referenced_labels[i])
            self.references_other.append(self.carrier_frequency_1)
            self.references_other.append(self.carrier_frequency_2)
            self.references_other.append(self.carrier_frequency_3)
            self.references_other_labels.append('O1/BF1')
            self.references_other_labels.append('O2/BF2')
            self.references_other_labels.append('O3/BF3')




            

            

            


# # This class creates the GUI

class MyApp(wx.Frame):
    def __init__(self):
        # Get the monitor size and set the window size to 85% of the monitor size
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 0.6*self.monitorWidth
        self.height = 0.6*self.monitorHeight
        self.app_frame = wx.Frame.__init__(self, None, wx.ID_ANY,'SpinConverter',wx.DefaultPosition, size=(int(self.width), int(self.height)))
        
        self.file_parser = False
        # Read the NMR data in the current directory
        
        self.set_variables()
        self.create_canvas()
        self.create_menu_bar()


        self.Bind(wx.EVT_CLOSE, self.OnClose)
        

        self.Show()
        self.Centre()





    def OnClose(self, event):
        self.Destroy()
        sys.exit()

    def set_variables(self):
        self.pseudo_flag = []
        self.threeDprocessing = False

    def create_canvas(self):

        # Create the main sizer
        self.main_sizer=wx.BoxSizer(wx.VERTICAL)


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


        # Get the initialise the NMR data class
        self.nmrdata = Converter(self)




        

    
    def create_menu_bar(self):
        # Have a row of buttons/TxtCtrl boxes at the top of the screen
        self.parameters_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.menu_bar = wx.BoxSizer(wx.VERTICAL)
        self.extra_sizers = wx.BoxSizer(wx.VERTICAL)
        self.parameters_sizer.Add(self.menu_bar)

        if(self.nmrdata.spectrometer == 'Bruker'):
            self.input_sizes_bruker()
            self.input_acquisition_modes_bruker()
            self.input_sweep_widths_bruker()
            self.get_nuclei_frequency_bruker()
            self.get_nuclei_labels_bruker()
            self.get_carrier_frequencies_bruker()
            if(len(self.N_complex_boxes) > 1):
                self.acquisition_2D_mode_combo_box()
            self.create_temperature_box()
            self.find_bruker_digital_filter_parameters()
            self.create_bruker_digital_filter_box()
            self.determine_byte_order()
            self.determine_byte_size()
            self.create_make_fid_conversion_box()
            self.create_other_options_box()
            self.create_intensity_scaling_box()
            if(self.nmrdata.size_indirect != []):
                self.find_nus_file()
                self.input_NUS_list_box()
            else:
                self.include_NUS = False
        elif(self.nmrdata.spectrometer == 'Varian'):
            self.input_sizes_varian()
            self.input_acquisition_modes_varian()
            self.input_sweep_widths_varian()
            self.get_nuclei_frequency_varian()
            self.get_nuclei_labels_varian()
            self.get_carrier_frequencies_varian()
            if(len(self.N_complex_boxes) > 1):
                self.acquisition_2D_mode_combo_box()
            self.create_temperature_box()
            self.create_make_fid_conversion_box()
            self.create_intensity_scaling_box()
            if(self.nmrdata.phase != False or self.nmrdata.phase2 != False):
                self.find_nus_file()
                self.input_NUS_list_box()
            else:
                self.include_NUS = False


        self.main_sizer.Add(self.parameters_sizer, 0, wx.CENTER)
        self.main_sizer.Add(self.extra_sizers, 0, wx.CENTER)

        self.SetSizerAndFit(self.main_sizer)

        # Get the width and height of the main_sizer
        self.width, self.height = self.main_sizer.GetSize()
        self.SetSize((int(self.width*1.25), int(self.height*1.25)))
        self.Centre()



    
    def input_sizes_varian(self):
        # Create the boxes for the complex point numbers
        self.title_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.title_sizer.AddSpacer(250)
        self.N_complex_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.N_complex_txt = wx.StaticText(self, label='Number complex points:')
        self.N_complex_boxes = []
        if(self.nmrdata.other_params == False):
            if(self.nmrdata.phase==False and self.nmrdata.phase2==False):
                for i in range(1):
                    self.title_sizer.Add(wx.StaticText(self, label='Dimension '+str(i+1)))
                    self.title_sizer.AddSpacer(145)
                    if(i==0):
                        self.N_complex_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.size_direct), size=(200,20)))
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==False):
                for i in range(2):
                    self.title_sizer.Add(wx.StaticText(self, label='Dimension '+str(i+1)))
                    self.title_sizer.AddSpacer(145)
                    if(i==0):
                        self.N_complex_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.size_direct), size=(200,20)))
                    if(i==1):
                        self.N_complex_boxes.append(wx.TextCtrl(self, value=str(int(self.nmrdata.size_indirect[i-1])*2), size=(200,20)))
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==True):
                for i in range(3):
                    self.title_sizer.Add(wx.StaticText(self, label='Dimension '+str(i+1)))
                    self.title_sizer.AddSpacer(145)
                    if(i==0):
                        self.N_complex_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.size_direct), size=(200,20)))
                    if(i==1):
                        self.N_complex_boxes.append(wx.TextCtrl(self, value=str(int(self.nmrdata.size_indirect[i-1])*2), size=(200,20)))
                    if(i==2):
                        self.N_complex_boxes.append(wx.TextCtrl(self, value=str(int(self.nmrdata.size_indirect[i-1])*2), size=(200,20)))
            elif(self.nmrdata.phase==False and self.nmrdata.phase2==True):
                for i in range(2):
                    self.title_sizer.Add(wx.StaticText(self, label='Dimension '+str(i+1)))
                    self.title_sizer.AddSpacer(145)
                    if(i==0):
                        self.N_complex_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.size_direct), size=(200,20)))
                    if(i==1):
                        self.N_complex_boxes.append(wx.TextCtrl(self, value=str(int(self.nmrdata.size_indirect[i-1])*2), size=(200,20)))
        else:
            if(self.nmrdata.phase==False and self.nmrdata.phase2==False):
                for i in range(2):
                    self.title_sizer.Add(wx.StaticText(self, label='Dimension '+str(i+1)))
                    self.title_sizer.AddSpacer(145)
                    if(i==0):
                        self.N_complex_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.size_direct), size=(200,20)))
                    if(i==1):
                        self.N_complex_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.number_of_arrayed_parameters), size=(200,20)))

            elif(self.nmrdata.phase==True and self.nmrdata.phase2==False):
                for i in range(3):
                    self.title_sizer.Add(wx.StaticText(self, label='Dimension '+str(i+1)))
                    self.title_sizer.AddSpacer(145)
                    if(i==0):
                        self.N_complex_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.size_direct), size=(200,20)))
                    if(i==1):
                        self.N_complex_boxes.append(wx.TextCtrl(self, value=str(int(self.nmrdata.size_indirect[i-1])*2), size=(200,20)))
                    if(i==2):
                        self.N_complex_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.number_of_arrayed_parameters), size=(200,20)))
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==True):
                for i in range(4):
                    self.title_sizer.Add(wx.StaticText(self, label='Dimension '+str(i+1)))
                    self.title_sizer.AddSpacer(145)
                    if(i==0):
                        self.N_complex_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.size_direct), size=(200,20)))
                    if(i==1):
                        self.N_complex_boxes.append(wx.TextCtrl(self, value=str(int(self.nmrdata.size_indirect[i-1])*2), size=(200,20)))
                    if(i==2):
                        self.N_complex_boxes.append(wx.TextCtrl(self, value=str(int(self.nmrdata.size_indirect[i-1])*2), size=(200,20)))
                    if(i==3):
                        self.N_complex_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.number_of_arrayed_parameters), size=(200,20)))
            elif(self.nmrdata.phase==False and self.nmrdata.phase2==True):
                for i in range(3):
                    self.title_sizer.Add(wx.StaticText(self, label='Dimension '+str(i+1)))
                    self.title_sizer.AddSpacer(145)
                    if(i==0):
                        self.N_complex_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.size_direct), size=(200,20)))
                    if(i==1):
                        self.N_complex_boxes.append(wx.TextCtrl(self, value=str(int(self.nmrdata.size_indirect[i-1])*2), size=(200,20)))
                    if(i==2):
                        self.N_complex_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.number_of_arrayed_parameters), size=(200,20)))
        
        self.N_complex_sizer.AddSpacer(20)
        self.N_complex_sizer.Add(self.N_complex_txt)
        self.N_complex_sizer.AddSpacer(20)
        for i in range(0,len(self.N_complex_boxes)):
            self.N_complex_sizer.Add(self.N_complex_boxes[i])
            self.N_complex_sizer.AddSpacer(20)
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.title_sizer)
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.N_complex_sizer)

        # Create the boxes for the real point numbers
        self.N_real_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.N_real_txt = wx.StaticText(self, label='Number real points:        ')
        self.N_real_boxes = []
        if(self.nmrdata.other_params == False):
            for i in range(len(self.N_complex_boxes)):
                if(i==0):
                    self.N_real_boxes.append(wx.TextCtrl(self, value=str(int(self.nmrdata.size_direct/2)), size=(200,20)))
                else:
                    if(i==1):
                        if(self.nmrdata.size_indirect[i-1] == 1):
                            continue
                        else:
                            self.N_real_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.size_indirect[i-1]), size=(200,20)))
                    if(i==2):
                        if(self.nmrdata.size_indirect[i-1] == 1):
                            continue
                        else:
                            self.N_real_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.size_indirect[i-1]), size=(200,20)))

        else:
            if(self.nmrdata.phase==False and self.nmrdata.phase2==False):
                for i in range(2):
                    if(i==0):
                        self.N_real_boxes.append(wx.TextCtrl(self, value=str(int(self.nmrdata.size_direct/2)), size=(200,20)))
                    if(i==1):
                        self.N_real_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.number_of_arrayed_parameters), size=(200,20)))

            elif(self.nmrdata.phase==True and self.nmrdata.phase2==False):
                for i in range(3):
                    if(i==0):
                        self.N_real_boxes.append(wx.TextCtrl(self, value=str(int(self.nmrdata.size_direct/2)), size=(200,20)))
                    if(i==1):
                        self.N_real_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.size_indirect[i-1]), size=(200,20)))
                    if(i==2):
                        self.N_real_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.number_of_arrayed_parameters), size=(200,20)))
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==True):
                for i in range(4):
                    if(i==0):
                        self.N_real_boxes.append(wx.TextCtrl(self, value=str(int(self.nmrdata.size_direct/2)), size=(200,20)))
                    if(i==1):
                        self.N_real_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.size_indirect[i-1]), size=(200,20)))
                    if(i==2):
                        self.N_real_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.size_indirect[i-1]), size=(200,20)))
                    if(i==3):
                        self.N_real_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.number_of_arrayed_parameters), size=(200,20)))
            elif(self.nmrdata.phase==False and self.nmrdata.phase2==True):
                for i in range(3):
                    if(i==0):
                        self.N_real_boxes.append(wx.TextCtrl(self, value=str(int(self.nmrdata.size_direct/2)), size=(200,20)))
                    if(i==1):
                        self.N_real_boxes.append(wx.TextCtrl(self, value=str(int(self.nmrdata.size_indirect[i-1]/2)), size=(200,20)))
                    if(i==2):
                        self.N_real_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.number_of_arrayed_parameters), size=(200,20)))


        
        self.N_real_sizer.AddSpacer(20)
        self.N_real_sizer.Add(self.N_real_txt)
        self.N_real_sizer.AddSpacer(20)
        for i in range(0,len(self.N_real_boxes)):
            self.N_real_sizer.Add(self.N_real_boxes[i])
            self.N_real_sizer.AddSpacer(20)
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.N_real_sizer)



    def input_sizes_bruker(self):
        # Create the boxes for the complex point numbers
        self.title_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.title_sizer.AddSpacer(250)
        self.N_complex_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.N_complex_txt = wx.StaticText(self, label='Number complex points:')
        self.N_complex_boxes = []
        for i in range(len(self.nmrdata.size_indirect)+1):
            self.title_sizer.Add(wx.StaticText(self, label='Dimension '+str(i+1)))
            self.title_sizer.AddSpacer(145)
            if(i==0):
                self.N_complex_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.size_direct_complex), size=(200,20)))
            else:
                try:
                    if(self.nmrdata.labels_correct_order[i] == 'ID'):
                        labelval = 'off'
                    else:
                        labelval = self.nmrdata.labels_correct_order[i]
                    size = str(self.nmrdata.indirect_sizes_dict[labelval])
                except:
                    size = str(self.nmrdata.size_indirect[i-1])
                self.N_complex_boxes.append(wx.TextCtrl(self, value=size, size=(200,20)))
        self.N_complex_sizer.AddSpacer(20)
        self.N_complex_sizer.Add(self.N_complex_txt)
        self.N_complex_sizer.AddSpacer(20)
        for i in range(0,len(self.N_complex_boxes)):
            self.N_complex_sizer.Add(self.N_complex_boxes[i])
            self.N_complex_sizer.AddSpacer(20)
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.title_sizer)
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.N_complex_sizer)

        # Create the boxes for the real point numbers
        self.N_real_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.N_real_txt = wx.StaticText(self, label='Number real points:        ')
        self.N_real_boxes = []
 
        if(self.pseudo_flag == []):
            for i in range(len(self.nmrdata.size_indirect)+1):
                if(i==0):
                    self.N_real_boxes.append(wx.TextCtrl(self, value=str(int(self.nmrdata.size_direct/2)), size=(200,20)))
                else:
                    
                    try:
                        if(self.nmrdata.labels_correct_order[i] == 'ID'):
                            labelval = 'off'
                            size = str(int(self.nmrdata.indirect_sizes_dict[labelval]))
                        else:
                            labelval = self.nmrdata.labels_correct_order[i]
                            size = str(int(self.nmrdata.indirect_sizes_dict[labelval]/2))

                    except:
                        size = str(int(self.nmrdata.size_indirect[i-1]/2))
                    self.N_real_boxes.append(wx.TextCtrl(self, value=size, size=(200,20)))

        else:
            for i in range(self.nmrdata.data_dimensions):
                if(i==0):
                    self.N_real_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.size_direct), size=(200,20)))
                else:
                    if(i in self.pseudo_flag):
                        try:
                            if(self.nmrdata.labels_correct_order[i] == 'ID'):
                                labelval = 'off'
                            else:
                                labelval = self.nmrdata.labels_correct_order[i]
                            size = str(self.nmrdata.indirect_sizes_dict[labelval])
                        except: 
                            size = str(self.nmrdata.size_indirect[i-1])
                        self.N_real_boxes.append(wx.TextCtrl(self, value=size), size=(200,20))
                    else:
                        try:
                            if(self.nmrdata.labels_correct_order[i] == 'ID'):
                                labelval = 'off'
                            else:
                                labelval = self.nmrdata.labels_correct_order[i]
                            size = str(self.nmrdata.indirect_sizes_dict[labelval])
                        except:
                            size = str(int(self.nmrdata.size_indirect[i-1]/2))
                        self.N_real_boxes.append(wx.TextCtrl(self, value=str(int(self.nmrdata.nmr_data.shape[self.nmrdata.data_dimensions-1-i]/2)), size=(200,20)))
        
        self.N_real_sizer.AddSpacer(20)
        self.N_real_sizer.Add(self.N_real_txt)
        self.N_real_sizer.AddSpacer(20)
        for i in range(0,len(self.N_real_boxes)):
            self.N_real_sizer.Add(self.N_real_boxes[i])
            self.N_real_sizer.AddSpacer(20)
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.N_real_sizer)


    
    def input_acquisition_modes_varian(self):
        # Create drop down options for the acquisition modes
        self.acquisition_mode_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.acquisition_mode_txt = wx.StaticText(self, label='Acquisition mode:           ')
        self.acquisition_mode_options_direct = ['Complex', 'Sequential', 'Real', 'DQD']
        self.acquisition_mode_options_indirect = ['Complex', 'States-TPPI', 'Rance-Kay','Echo-AntiEcho', 'TPPI', 'States', 'Real']

        self.acqusition_combo_boxes = []
        if(self.nmrdata.other_params == False):
            if(self.nmrdata.phase == False and self.nmrdata.phase2 == False):
                self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_direct[0], choices=self.acquisition_mode_options_direct, size=(200,20)))
                self.acqusition_combo_boxes[0].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)
            elif(self.nmrdata.phase == True and self.nmrdata.phase2 == False):
                for i in range(2):
                    if(i==0):
                        self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_direct[0], choices=self.acquisition_mode_options_direct, size=(200,20)))
                        self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)
                    if(i==1):
                        self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_indirect[0], choices=self.acquisition_mode_options_indirect, size=(200,20)))
                        self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)
            elif(self.nmrdata.phase == True and self.nmrdata.phase2 == True):
                for i in range(3):
                    if(i==0):
                        self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_direct[0], choices=self.acquisition_mode_options_direct, size=(200,20)))
                        self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)
                    if(i==1):
                        self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_indirect[0], choices=self.acquisition_mode_options_indirect, size=(200,20)))
                        self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)
                    if(i==2):
                        self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_indirect[0], choices=self.acquisition_mode_options_indirect, size=(200,20)))
                        self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)
        else:
            if(self.nmrdata.phase == False and self.nmrdata.phase2 == False):
                for i in range(2):
                    if(i==0):
                        self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_direct[0], choices=self.acquisition_mode_options_direct, size=(200,20)))
                        self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)
                    if(i==1):
                        self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_indirect[6], choices=self.acquisition_mode_options_indirect, size=(200,20)))
                        self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)
            elif(self.nmrdata.phase == True and self.nmrdata.phase2 == False):
                for i in range(3):
                    if(i==0):
                        self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_direct[0], choices=self.acquisition_mode_options_direct, size=(200,20)))
                        self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)
                    if(i==1):
                        self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_indirect[0], choices=self.acquisition_mode_options_indirect, size=(200,20)))
                        self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)
                    if(i==2):
                        self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_indirect[6], choices=self.acquisition_mode_options_indirect, size=(200,20)))
                        self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)
            elif(self.nmrdata.phase == True and self.nmrdata.phase2 == True):
                for i in range(4):
                    if(i==0):
                        self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_direct[0], choices=self.acquisition_mode_options_direct, size=(200,20)))
                        self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)
                    if(i==1):
                        self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_indirect[0], choices=self.acquisition_mode_options_indirect, size=(200,20)))
                        self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)
                    if(i==2):
                        self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_indirect[0], choices=self.acquisition_mode_options_indirect, size=(200,20)))
                        self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)
                    if(i==3):
                        self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_indirect[6], choices=self.acquisition_mode_options_indirect, size=(200,20)))
                        self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)
            elif(self.nmrdata.phase == False and self.nmrdata.phase2 == True):
                for i in range(3):
                    if(i==0):
                        self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_direct[0], choices=self.acquisition_mode_options_direct, size=(200,20)))
                        self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)
                    if(i==1):
                        self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_indirect[0], choices=self.acquisition_mode_options_indirect, size=(200,20)))
                        self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)
                    if(i==2):
                        self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_indirect[6], choices=self.acquisition_mode_options_indirect, size=(200,20)))
                        self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)


        self.acquisition_mode_sizer.AddSpacer(20)
        self.acquisition_mode_sizer.Add(self.acquisition_mode_txt)
        self.acquisition_mode_sizer.AddSpacer(20)
        for i in range(0,len(self.acqusition_combo_boxes)):
            self.acquisition_mode_sizer.Add(self.acqusition_combo_boxes[i])
            self.acquisition_mode_sizer.AddSpacer(20)
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.acquisition_mode_sizer)


    
    def input_acquisition_modes_bruker(self):
        # Create drop down options for the acquisition modes
        self.acquisition_mode_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.acquisition_mode_txt = wx.StaticText(self, label='Acquisition mode:           ')
        self.acquisition_mode_options_direct = ['DQD', 'Complex', 'Sequential','Real']
        self.acquisition_mode_options_indirect = ['Complex', 'States-TPPI', 'Echo-AntiEcho', 'TPPI', 'States', 'Real']

        self.acqusition_combo_boxes = []
        for i in range(len(self.nmrdata.size_indirect)+1):
            if(i==0):
                self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_direct[0], choices=self.acquisition_mode_options_direct, size=(200,20)))
                self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)
            else:
                if(self.nmrdata.acqusition_modes[i-1] == 'QF'):
                    self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_indirect[5], choices=self.acquisition_mode_options_indirect, size=(200,20)))
                    self.N_real_boxes[i].SetValue(str(self.N_complex_boxes[i].GetValue()))
                    self.nmrdata.sw_indirect[i-1] = 0
                    self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)
                else:
                    determined_value = self.nmrdata.acqusition_modes[i-1]
                    detected_value = ''
                    for j in range(len(self.acquisition_mode_options_indirect)):
                        if(self.acquisition_mode_options_indirect[j].upper() == determined_value.upper()):
                            detected_value = self.acquisition_mode_options_indirect[j]
                            break
                    if(detected_value == ''):
                        self.acqusition_combo_boxes.append(wx.ComboBox(self, value=self.acquisition_mode_options_indirect[0], choices=self.acquisition_mode_options_indirect, size=(200,20)))
                        self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)
                    else:
                        self.acqusition_combo_boxes.append(wx.ComboBox(self, value=detected_value, choices=self.acquisition_mode_options_indirect, size=(200,20)))
                        self.acqusition_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_acquisition_mode_change)

        self.acquisition_mode_sizer.AddSpacer(20)
        self.acquisition_mode_sizer.Add(self.acquisition_mode_txt)
        self.acquisition_mode_sizer.AddSpacer(20)
        for i in range(0,len(self.acqusition_combo_boxes)):
            self.acquisition_mode_sizer.Add(self.acqusition_combo_boxes[i])
            self.acquisition_mode_sizer.AddSpacer(20)
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.acquisition_mode_sizer)

    
    def on_acquisition_mode_change(self, event):
        # If combobox is changed to real, then change the number of real points to the same as the number of complex points, also set label to ID, sw to 0
        # Find out which combobox has been changed
        for i in range(len(self.acqusition_combo_boxes)):
            if(self.acqusition_combo_boxes[i] == event.GetEventObject()):
                index = i
                break
        if(self.acqusition_combo_boxes[index].GetValue() == 'Real'):
            if(index == 0):
                self.N_real_boxes[index].SetValue(str(self.nmrdata.size_direct))
                # Set label to ID
                self.nucleus_type_boxes[index].SetValue('ID')
                # Set sw to 0
                self.sweep_width_boxes[index].SetValue(str(1))
                # Set nucleus frequency to 0
                self.nuclei_frequency_boxes[index].SetValue(str(1))
                # Set carrier frequency to 0
                self.carrier_frequency_boxes[index].SetValue(str(1))
                # Set combobox to N/A
                self.carrier_combo_boxes[index].SetValue('N/A')

            else:
                self.N_real_boxes[index].SetValue(self.N_complex_boxes[index].GetValue())
                # Set label to ID
                self.nucleus_type_boxes[index].SetValue('ID')
                # Set sw to 0
                self.sweep_width_boxes[index].SetValue(str(1))
                # Set nucleus frequency to 0
                self.nuclei_frequency_boxes[index].SetValue(str(1))
                # Set carrier frequency to 0
                self.carrier_frequency_boxes[index].SetValue(str(1))
                # Set combobox to N/A
                self.carrier_combo_boxes[index].SetValue('N/A')
        else:
            if(index == 0):
                self.N_real_boxes[index].SetValue(str(int(self.nmrdata.size_direct/2)))
                # Set label back to original label
                self.nucleus_type_boxes[index].SetValue(self.nmrdata.labels_correct_order[index])
                # Set sw back to original sw
                self.sweep_width_boxes[index].SetValue(str(self.nmrdata.sw_direct))
                # Set nucleus frequency back to original nucleus frequency
                self.nuclei_frequency_boxes[index].SetValue(str(self.nmrdata.nucleus_frequencies[index]))
                # Set carrier frequency back to original carrier frequency
                self.carrier_frequency_boxes[index].SetValue(str(self.nmrdata.references_proton[index]))
                # Set combobox back to original combobox
                self.carrier_combo_boxes[index].SetValue(self.nmrdata.references_proton_labels[index])
            else:
                if(self.nmrdata.spectrometer=='Bruker'):
                    self.N_real_boxes[index].SetValue(str(int(int(self.N_complex_boxes[index].GetValue())/2)))
                else:
                    self.N_real_boxes[index].SetValue(str(int(self.nmrdata.size_indirect[index-1])))
                # Set label back to original label
                if(self.nmrdata.spectrometer=='Bruker'):
                    self.nucleus_type_boxes[index].SetValue(self.nmrdata.labels_correct_order[index])
                else:
                    self.nucleus_type_boxes[index].SetValue(self.nmrdata.labels_correct_order[index-1])
                # Set sw back to original sw
                if(self.nmrdata.spectrometer == 'Bruker'):
                    self.sweep_width_boxes[index].SetValue(str(self.nmrdata.sw_indirect[index-1]))
                else:
                    self.sweep_width_boxes[index].SetValue(str(self.nmrdata.sw_indirect['sw'+str(index)]))
                # Set nucleus frequency back to original nucleus frequency
                self.nuclei_frequency_boxes[index].SetValue(str(self.nmrdata.nucleus_frequencies[index]))

                # Set carrier frequency back to original carrier frequency
                self.carrier_frequency_boxes[index].SetValue(str(self.nmrdata.references_other[index-1]))
                # Set combobox back to original combobox
                self.carrier_combo_boxes[index].SetValue(self.nmrdata.references_other_labels[index-1])



    def input_sweep_widths_bruker(self):
        # Sweep width sizer
        self.sweep_width_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.sweep_width_txt = wx.StaticText(self, label='Sweep width (Hz):          ')
        self.sweep_width_boxes = []
        choices = [str(self.nmrdata.sw_direct)]
        for j in range(len(self.nmrdata.size_indirect)):
            choices.append(str(self.nmrdata.sw_indirect[j]))
        for i in range(len(self.nmrdata.size_indirect)+1):
            if(i==0):
                self.sweep_width_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.sw_direct),choices=choices, size=(200,20)))
            else:
                self.sweep_width_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.sw_indirect[i-1]), choices = choices, size=(200,20)))
        self.sweep_width_sizer.AddSpacer(20)
        self.sweep_width_sizer.Add(self.sweep_width_txt)
        self.sweep_width_sizer.AddSpacer(20)
        for i in range(0,len(self.sweep_width_boxes)):
            self.sweep_width_sizer.Add(self.sweep_width_boxes[i])
            self.sweep_width_sizer.AddSpacer(20)
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.sweep_width_sizer)


    def input_sweep_widths_varian(self):
        # Sweep width sizer
        self.sweep_width_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.sweep_width_txt = wx.StaticText(self, label='Sweep width (Hz):          ')
        self.sweep_width_boxes = []
        options_sw = [str(self.nmrdata.sw_direct)]
        try: 
            options_sw = options_sw + [str(self.nmrdata.sw_indirect['sw1'])]
            try:
                options_sw = options_sw + [str(self.nmrdata.sw_indirect['sw2'])]
                try:
                    options_sw = options_sw + [str(self.nmrdata.sw_indirect['sw3'])] + ['0.0']
                except:
                    options_sw = options_sw + ['0.0']
            except:
                options_sw = options_sw + ['0.0']
        except:
            options_sw = options_sw + ['0.0'] 
        if(self.nmrdata.other_params == False):
            if(self.nmrdata.phase==False and self.nmrdata.phase2==False):
                self.sweep_width_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.sw_direct),choices=options_sw, size=(200,20)))
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==False):
                for i in range(2):
                    if(i==0):
                        self.sweep_width_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.sw_direct),choices=options_sw, size=(200,20)))
                    if(i==1):
                        self.sweep_width_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.sw_indirect['sw1']),choices=options_sw, size=(200,20)))
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==True):
                for i in range(3):
                    if(i==0):
                        self.sweep_width_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.sw_direct),choices=options_sw, size=(200,20)))
                    if(i==1):
                        self.sweep_width_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.sw_indirect['sw1']),choices=options_sw, size=(200,20)))
                    if(i==2):
                        self.sweep_width_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.sw_indirect['sw2']),choices=options_sw, size=(200,20)))
            elif(self.nmrdata.phase==False and self.nmrdata.phase2==True):
                for i in range(2):
                    if(i==0):
                        self.sweep_width_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.sw_direct),choices=options_sw, size=(200,20)))
                    if(i==1):
                        self.sweep_width_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.sw_indirect['sw2']),choices=options_sw, size=(200,20)))
        else:
            if(self.nmrdata.phase==False and self.nmrdata.phase2==False):
                for i in range(2):
                    if(i==0):
                        self.sweep_width_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.sw_direct), choices=options_sw,size=(200,20)))
                    if(i==1):
                        self.sweep_width_boxes.append(wx.ComboBox(self, value=str(0),choices=options_sw, size=(200,20)))
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==False):
                for i in range(3):
                    if(i==0):
                        self.sweep_width_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.sw_direct),choices=options_sw, size=(200,20)))
                    if(i==1):
                        self.sweep_width_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.sw_indirect['sw1']),choices=options_sw, size=(200,20)))
                    if(i==2):
                        self.sweep_width_boxes.append(wx.ComboBox(self, value=str(0),choices=options_sw, size=(200,20)))
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==True):
                for i in range(4):
                    if(i==0):
                        self.sweep_width_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.sw_direct),choices=options_sw, size=(200,20)))
                    if(i==1):
                        self.sweep_width_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.sw_indirect['sw1']), choices=options_sw,size=(200,20)))
                    if(i==2):
                        self.sweep_width_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.sw_indirect['sw2']), choices=options_sw,size=(200,20)))
                    if(i==3):
                        self.sweep_width_boxes.append(wx.ComboBox(self, value=str(0),choices=options_sw, size=(200,20)))
            elif(self.nmrdata.phase==False and self.nmrdata.phase2==True):
                for i in range(3):
                    if(i==0):
                        self.sweep_width_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.sw_direct),choices=options_sw, size=(200,20)))
                    if(i==1):
                        self.sweep_width_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.sw_indirect['sw2']),choices=options_sw, size=(200,20)))
                    if(i==2):
                        self.sweep_width_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.number_of_arrayed_parameters),choices=options_sw, size=(200,20)))
        self.sweep_width_sizer.AddSpacer(20)
        self.sweep_width_sizer.Add(self.sweep_width_txt)
        self.sweep_width_sizer.AddSpacer(20)
        for i in range(0,len(self.sweep_width_boxes)):
            self.sweep_width_sizer.Add(self.sweep_width_boxes[i])
            self.sweep_width_sizer.AddSpacer(20)
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.sweep_width_sizer)



    def get_nuclei_frequency_bruker(self):
        # Nuclei frequency sizer (in MHz)
        # If Bruker, the frequencies are contained in the acqus file and called SFO1, SFO2, SFO3, etc. If Varian, the frequencies are contained in the procpar file and called sfrq, sfrq2, sfrq3, etc.
        self.nuclei_frequency_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.nuclei_frequency_txt = wx.StaticText(self, label='Nuclei frequency (MHz):')
        self.nuclei_frequency_boxes = []
        for i in range(len(self.nmrdata.size_indirect)+1):
            if(i==0):
                self.nuclei_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.nucleus_frequencies[0]), size=(200,20)))
            elif(self.acqusition_combo_boxes[i].GetValue() == 'Real'):
                self.nuclei_frequency_boxes.append(wx.TextCtrl(self, value=str(1), size=(200,20)))
                self.nmrdata.labels_correct_order[i] = 'ID'
            else:
                self.nuclei_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.nucleus_frequencies[i]), size=(200,20)))
            
        self.nuclei_frequency_sizer.AddSpacer(20)
        self.nuclei_frequency_sizer.Add(self.nuclei_frequency_txt)
        self.nuclei_frequency_sizer.AddSpacer(19)
        for i in range(0,len(self.nuclei_frequency_boxes)):
            self.nuclei_frequency_sizer.Add(self.nuclei_frequency_boxes[i])
            self.nuclei_frequency_sizer.AddSpacer(20)
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.nuclei_frequency_sizer)


    def get_nuclei_frequency_varian(self):
        self.nuclei_frequency_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.nuclei_frequency_txt = wx.StaticText(self, label='Nuclei frequency (MHz):')
        self.nuclei_frequency_boxes = []
        options = [str(self.nmrdata.nucleus_frequency_direct)] + [str(freq) for freq in self.nmrdata.nucleus_frequencies_indirect]
        if(self.nmrdata.other_params == False):
            if(self.nmrdata.phase==False and self.nmrdata.phase2==False):
                self.nuclei_frequency_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.nucleus_frequency_direct), choices = options, size=(200,20)))
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==False):
                for i in range(2):
                    if(i==0):
                        self.nuclei_frequency_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.nucleus_frequency_direct), choices = options, size=(200,20)))
                    if(i==1):
                        if(self.nmrdata.nucleus_frequencies_indirect[0] == 0):
                            try:
                                self.nuclei_frequency_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.nucleus_frequencies_indirect[1]), choices=options,size=(200,20)))
                            except:
                                self.nuclei_frequency_boxes.append(wx.ComboBox(self, value=str(1),choices=options, size=(200,20)))
                        else:
                            self.nuclei_frequency_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.nucleus_frequencies_indirect[0]),choices=options, size=(200,20)))
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==True):
                for i in range(3):
                    if(i==0):
                        self.nuclei_frequency_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.nucleus_frequency_direct),choices=options, size=(200,20)))
                    if(i==1):
                        self.nuclei_frequency_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.nucleus_frequencies_indirect[0]),choices=options, size=(200,20)))
                    if(i==2):
                        self.nuclei_frequency_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.nucleus_frequencies_indirect[1]), choices=options,size=(200,20)))
            elif(self.nmrdata.phase==False and self.nmrdata.phase2==True):
                for i in range(2):
                    if(i==0):
                        self.nuclei_frequency_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.nucleus_frequency_direct),choices=options, size=(200,20)))
                    if(i==1):
                        self.nuclei_frequency_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.nucleus_frequencies_indirect[0]),choices=options, size=(200,20)))
        else:
            if(self.nmrdata.phase==False and self.nmrdata.phase2==False):
                for i in range(2):
                    if(i==0):
                        self.nuclei_frequency_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.nucleus_frequency_direct),choices=options, size=(200,20)))
                    if(i==1):
                        self.nuclei_frequency_boxes.append(wx.ComboBox(self, value='0.0',choices=options, size=(200,20)))
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==False):
                for i in range(3):
                    if(i==0):
                        self.nuclei_frequency_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.nucleus_frequency_direct),choices=options, size=(200,20)))
                    if(i==1):
                        self.nuclei_frequency_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.nucleus_frequencies_indirect[0]),choices=options, size=(200,20)))
                    if(i==2):
                        self.nuclei_frequency_boxes.append(wx.ComboBox(self, value='0.0',choices=options, size=(200,20)))
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==True):
                for i in range(4):
                    if(i==0):
                        self.nuclei_frequency_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.nucleus_frequency_direct),choices=options, size=(200,20)))
                    if(i==1):
                        self.nuclei_frequency_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.nucleus_frequencies_indirect[0]),choices=options, size=(200,20)))
                    if(i==2):
                        self.nuclei_frequency_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.nucleus_frequencies_indirect[1]),choices=options, size=(200,20)))
                    if(i==3):
                        self.nuclei_frequency_boxes.append(wx.ComboBox(self, value='0.0',choices=options, size=(200,20)))
            elif(self.nmrdata.phase==False and self.nmrdata.phase2==True):
                for i in range(3):
                    if(i==0):
                        self.nuclei_frequency_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.nucleus_frequency_direct),choices=options, size=(200,20)))
                    if(i==1):
                        self.nuclei_frequency_boxes.append(wx.ComboBox(self, value=str(self.nmrdata.nucleus_frequencies_indirect[0]),choices=options, size=(200,20)))
                    if(i==2):
                        self.nuclei_frequency_boxes.append(wx.ComboBox(self, value='0.0',choices=options, size=(200,20)))

        self.nuclei_frequency_sizer.AddSpacer(20)
        self.nuclei_frequency_sizer.Add(self.nuclei_frequency_txt)
        self.nuclei_frequency_sizer.AddSpacer(19)
        for i in range(0,len(self.nuclei_frequency_boxes)):
            self.nuclei_frequency_sizer.Add(self.nuclei_frequency_boxes[i])
            self.nuclei_frequency_sizer.AddSpacer(20)
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.nuclei_frequency_sizer)
            
                      


    def get_nuclei_labels_bruker(self):
         # Add a section for axis labels
        self.nucleus_type_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.nucleus_type_txt = wx.StaticText(self, label='Label:                              ')
        self.nucleus_type_boxes = []
        j = 0
        for i in range(len(self.nmrdata.size_indirect)+1):
            if(i==0):
                self.nucleus_type_boxes.append(wx.TextCtrl(self, value=self.nmrdata.labels_correct_order[0], size=(200,20)))
            elif(self.acqusition_combo_boxes[i].GetValue() == 'Real'):
                self.nucleus_type_boxes.append(wx.TextCtrl(self, value='ID', size=(200,20)))
                j += 1
            else:
                self.nucleus_type_boxes.append(wx.TextCtrl(self, value=self.nmrdata.labels_correct_order[i], size=(200,20)))
            
        self.nucleus_type_sizer.AddSpacer(20)
        self.nucleus_type_sizer.Add(self.nucleus_type_txt)
        self.nucleus_type_sizer.AddSpacer(22)
        for i in range(0,len(self.nucleus_type_boxes)):
            self.nucleus_type_sizer.Add(self.nucleus_type_boxes[i])
            self.nucleus_type_sizer.AddSpacer(20)
        
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.nucleus_type_sizer)


    def get_nuclei_labels_varian(self):
        # Add a section for axis labels
        self.nucleus_type_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.nucleus_type_txt = wx.StaticText(self, label='Label:                              ')
        self.nucleus_type_boxes = []
        self.nmrdata.labels_indirect.remove('')
        options = [self.nmrdata.label_direct] + [self.nmrdata.labels_indirect[i] for i in range(len(self.nmrdata.labels_indirect))] + ['ID']
        if(self.nmrdata.other_params == True):
            options = options + [self.nmrdata.arrayed_parameter]
        j = 0

        if(self.nmrdata.other_params == False):
            if(self.nmrdata.phase==False and self.nmrdata.phase2==False):
                    self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.label_direct,choices=options, size=(200,20)))
                    
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==False):       
                for i in range(2):
                    if(i==0):
                        self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.label_direct,choices=options, size=(200,20)))
                    if(i==1):
                        self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.labels_indirect[0],choices=options, size=(200,20)))
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==True):
                for i in range(3):
                    if(i==0):
                        self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.label_direct,choices=options, size=(200,20)))
                    if(i==1):
                        self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.labels_indirect[0],choices=options, size=(200,20)))
                    if(i==2):
                        self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.labels_indirect[1],choices=options, size=(200,20)))
            elif(self.nmrdata.phase==False and self.nmrdata.phase2==True):
                for i in range(2):
                    if(i==0):
                        self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.label_direct,choices=options, size=(200,20)))
                    if(i==1):
                        self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.labels_indirect[1],choices=options, size=(200,20)))
            else:
                for i in range(2):
                    if(i==0):
                        self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.label_direct,choices=options, size=(200,20)))
                    if(i==1):
                        self.nucleus_type_boxes.append(wx.ComboBox(self, value='ID',choices=options, size=(200,20)))
        else:
            if(self.nmrdata.phase==False and self.nmrdata.phase2==False):
                for i in range(2):
                    if(i==0):
                        self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.label_direct,choices=options, size=(200,20)))
                    if(i==1):
                        try:
                            self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.arrayed_parameter,choices=options, size=(200,20)))
                        except:
                            self.nucleus_type_boxes.append(wx.ComboBox(self, value='ID',choices=options, size=(200,20)))
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==False):
                for i in range(3):
                    if(i==0):
                        self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.label_direct,choices=options, size=(200,20)))
                    if(i==1):
                        self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.labels_indirect[0],choices=options, size=(200,20)))
                    if(i==2):
                        try:
                            self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.arrayed_parameter,choices=options, size=(200,20)))
                        except:
                            self.nucleus_type_boxes.append(wx.ComboBox(self, value='ID',choices=options, size=(200,20)))
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==True):
                for i in range(4):
                    if(i==0):
                        self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.label_direct,choices=options, size=(200,20)))
                    if(i==1):
                        self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.labels_indirect[0],choices=options, size=(200,20)))
                    if(i==2):
                        self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.labels_indirect[1],choices=options, size=(200,20)))
                    if(i==3):
                        try:
                            self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.arrayed_parameter,choices=options, size=(200,20)))
                        except:
                            self.nucleus_type_boxes.append(wx.ComboBox(self, value='ID',choices=options, size=(200,20)))
            elif(self.nmrdata.phase==False and self.nmrdata.phase2==True):
                for i in range(3):
                    if(i==0):
                        self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.label_direct,choices=options, size=(200,20)))
                    if(i==1):
                        self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.labels_indirect[1],choices=options, size=(200,20)))
                    if(i==2):
                        try:
                            self.nucleus_type_boxes.append(wx.ComboBox(self, value=self.nmrdata.arrayed_parameter,choices=options, size=(200,20)))
                        except:
                            self.nucleus_type_boxes.append(wx.ComboBox(self, value='ID',choices=options, size=(200,20)))


            
        self.nucleus_type_sizer.AddSpacer(20)
        self.nucleus_type_sizer.Add(self.nucleus_type_txt)
        self.nucleus_type_sizer.AddSpacer(22)
        for i in range(0,len(self.nucleus_type_boxes)):
            self.nucleus_type_sizer.Add(self.nucleus_type_boxes[i])
            self.nucleus_type_sizer.AddSpacer(20)
        
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.nucleus_type_sizer)


    def get_carrier_frequencies_bruker(self):
        # Add a section for the carrier frequencies
        self.carrier_frequency_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.carrier_frequency_txt = wx.StaticText(self, label='Carrier frequency (ppm):')
        self.carrier_frequency_boxes = []
        self.carrier_combo_boxes = []
        if(self.nmrdata.size_indirect == []):
            self.options_other_dimensions = []
        else:
            self.options_other_dimensions = self.nmrdata.references_other_labels

        id_columns = 0
        for i in range(len(self.nmrdata.size_indirect)+1):
            if(i==0):
                if(self.nmrdata.labels_correct_order[i] == '1H' or self.nmrdata.labels_correct_order[i] == 'H1' or self.nmrdata.labels_correct_order[i] == 'H'):
                    self.carrier_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.water_ppm), size=(200,20)))
                    self.carrier_combo_boxes.append(wx.ComboBox(self, value=self.nmrdata.references_proton_labels[0], choices=self.nmrdata.references_proton_labels, size=(200,20)))
                else:
                    self.carrier_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.references_proton[i]), size=(200,20)))
                    self.carrier_combo_boxes.append(wx.ComboBox(self, value=self.nmrdata.references_proton_labels[0], choices=self.nmrdata.references_proton_labels, size=(200,20)))
            elif(self.nmrdata.labels_correct_order[i] == 'ID'):
                self.carrier_frequency_boxes.append(wx.TextCtrl(self, value='0.00', size=(200,20)))
                self.carrier_combo_boxes.append(wx.ComboBox(self, value = 'N/A', choices = self.options_other_dimensions, size = (200,20)))
                id_columns += 1
            else:
                self.carrier_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.references_other[i-1]), size=(200,20)))
                self.carrier_combo_boxes.append(wx.ComboBox(self, value=self.nmrdata.references_other_labels[i-1], choices=self.options_other_dimensions, size=(200,20)))

        self.carrier_frequency_sizer.AddSpacer(20)      
        self.carrier_frequency_sizer.Add(self.carrier_frequency_txt)
        self.carrier_frequency_sizer.AddSpacer(16)
        for i in range(0,len(self.carrier_frequency_boxes)):
            self.carrier_frequency_sizer.Add(self.carrier_frequency_boxes[i])
            self.carrier_frequency_sizer.AddSpacer(20)
        

        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.carrier_frequency_sizer)

        self.carrier_combo_box_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.carrier_combo_box_sizer.AddSpacer(20)
        self.carrier_combo_box_txt = wx.StaticText(self, label='Referencing Mode:          ')
        self.carrier_combo_box_sizer.Add(self.carrier_combo_box_txt)
        self.carrier_combo_box_sizer.AddSpacer(17)
        for i in range(0,len(self.carrier_combo_boxes)):
            if(self.carrier_combo_boxes[i] == 200):
                self.carrier_combo_box_sizer.AddSpacer(200)
                self.carrier_combo_box_sizer.AddSpacer(20)
            else:
                self.carrier_combo_box_sizer.Add(self.carrier_combo_boxes[i])
                self.carrier_combo_box_sizer.AddSpacer(20)
                self.carrier_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_carrier_combo_box_change)



        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.carrier_combo_box_sizer)




    def on_carrier_combo_box_change(self, event):
        # Find the index of the combo box that was changed
        index = self.carrier_combo_boxes.index(event.GetEventObject())
        # Find the value of the combo box
        value = event.GetEventObject().GetValue()
        index_selection = self.carrier_combo_boxes[index].GetSelection()

        if(self.nmrdata.labels_correct_order[index] == '1H' or self.nmrdata.labels_correct_order[index] == 'H1' or self.nmrdata.labels_correct_order[index] == 'H'):
            if(value == 'H2O'):
                self.carrier_frequency_boxes[index].SetValue(str(self.nmrdata.references_proton[0]))
            elif(value == 'O1/BF1'):
                self.carrier_frequency_boxes[index].SetValue(str(self.nmrdata.references_proton[1]))
            else:
                self.carrier_frequency_boxes[index].SetValue(str(self.nmrdata.references_other[index_selection]))
        else:
            if(value == 'O1/BF1'):
                self.carrier_frequency_boxes[index].SetValue(str(self.nmrdata.references_proton[0]))
            else:
                self.carrier_frequency_boxes[index].SetValue(str(self.nmrdata.references_other[index_selection]))


    def get_carrier_frequencies_varian(self):
        # Add a section for the carrier frequencies
        self.carrier_frequency_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.carrier_frequency_txt = wx.StaticText(self, label='Carrier frequency (ppm):')
        self.carrier_frequency_boxes = []
        self.carrier_combo_boxes = []
        if(self.nmrdata.label_direct == '1H' or self.nmrdata.label_direct == 'H1'):
            self.options_proton = ['H2O', 'Other']
        else:
            self.options_proton = ['Manual']
        if(self.nmrdata.size_indirect == []):
            self.options_other_dimensions = ['Other']
        else:
            self.options_other_dimensions = self.nmrdata.references_other_labels + ['Other']

        id_columns = 0

        if(self.nmrdata.other_params == False):
            if(self.nmrdata.phase==False and self.nmrdata.phase2==False):
                self.carrier_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.references_proton[0]), size=(200,20)))
                self.carrier_combo_boxes.append(wx.ComboBox(self, value=self.options_proton[0], choices=self.options_proton, size=(200,20)))

            elif(self.nmrdata.phase==True and self.nmrdata.phase2==False):
                for i in range(2):
                    if(i==0):
                        self.carrier_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.references_proton[0]), size=(200,20)))
                        self.carrier_combo_boxes.append(wx.ComboBox(self, value=self.options_proton[0], choices=self.options_proton, size=(200,20)))
                    if(i==1):
                            self.carrier_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.references_other[0]), size=(200,20)))
                            self.carrier_combo_boxes.append(wx.ComboBox(self, value=self.options_other_dimensions[0], choices=self.options_other_dimensions, size=(200,20)))
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==True):
                for i in range(3):
                    if(i==0):
                        self.carrier_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.references_proton[0]), size=(200,20)))
                        self.carrier_combo_boxes.append(wx.ComboBox(self, value=self.nmrdata.references_proton_labels[0], choices=self.options_proton, size=(200,20)))
                    if(i==1):
                            self.carrier_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.references_other[0]), size=(200,20)))
                            self.carrier_combo_boxes.append(wx.ComboBox(self, value=self.options_other_dimensions[0], choices=self.options_other_dimensions, size=(200,20)))
                    if(i==2):
                            self.carrier_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.references_other[1]), size=(200,20)))
                            self.carrier_combo_boxes.append(wx.ComboBox(self, value=self.options_other_dimensions[1], choices=self.options_other_dimensions, size=(200,20)))
            elif(self.nmrdata.phase==False and self.nmrdata.phase2==True):
                for i in range(2):
                    if(i==0):
                        self.carrier_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.references_proton[0]), size=(200,20)))
                        self.carrier_combo_boxes.append(wx.ComboBox(self, value=self.nmrdata.references_proton_labels[0], choices=self.options_proton, size=(200,20)))
                    if(i==1):
                            self.carrier_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.references_other[1]), size=(200,20)))
                            self.carrier_combo_boxes.append(wx.ComboBox(self, value=self.options_other_dimensions[1], choices=self.options_other_dimensions, size=(200,20)))

        else:
            if(self.nmrdata.phase==False and self.nmrdata.phase2==False):
                for i in range(2):
                    if(i==0):
                        self.carrier_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.references_proton[0]), size=(200,20)))
                        self.carrier_combo_boxes.append(wx.ComboBox(self, value=self.nmrdata.references_proton_labels[0], choices=self.options_proton, size=(200,20)))
                    if(i==1):
                        try:
                            self.carrier_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.references_other[0]), size=(200,20)))
                            self.carrier_combo_boxes.append(wx.ComboBox(self, value=self.options_other_dimensions[0], choices=self.options_other_dimensions, size=(200,20)))
                        except:
                            self.carrier_frequency_boxes.append(wx.TextCtrl(self, value='0.0', size=(200,20)))
                            self.carrier_combo_boxes.append(wx.ComboBox(self, value='ID', choices=self.options_other_dimensions, size=(200,20)))
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==False):
                for i in range(3):
                    if(i==0):
                        self.carrier_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.references_proton[0]), size=(200,20)))
                        self.carrier_combo_boxes.append(wx.ComboBox(self, value=self.nmrdata.references_proton_labels[0], choices=self.options_proton, size=(200,20)))
                    if(i==1):
                            self.carrier_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.references_other[0]), size=(200,20)))
                            self.carrier_combo_boxes.append(wx.ComboBox(self, value=self.options_other_dimensions[0], choices=self.options_other_dimensions, size=(200,20)))
                    if(i==2):
                            self.carrier_frequency_boxes.append(wx.TextCtrl(self, value='0.0', size=(200,20)))
                            self.carrier_combo_boxes.append(wx.ComboBox(self, value='Other', choices=self.options_other_dimensions, size=(200,20)))
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==True):
                for i in range(4):
                    if(i==0):
                        self.carrier_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.references_proton[0]), size=(200,20)))
                        self.carrier_combo_boxes.append(wx.ComboBox(self, value=self.nmrdata.references_proton_labels[0], choices=self.options_proton, size=(200,20)))
                    if(i==1):
                            self.carrier_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.references_other[0]), size=(200,20)))
                            self.carrier_combo_boxes.append(wx.ComboBox(self, value=self.options_other_dimensions[0], choices=self.options_other_dimensions, size=(200,20)))
                    if(i==2):
                            self.carrier_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.references_other[1]), size=(200,20)))
                            self.carrier_combo_boxes.append(wx.ComboBox(self, value=self.options_other_dimensions[1], choices=self.options_other_dimensions, size=(200,20)))
                    if(i==3):
                            self.carrier_frequency_boxes.append(wx.TextCtrl(self, value='0.0', size=(200,20)))
                            self.carrier_combo_boxes.append(wx.ComboBox(self, value='Other', choices=self.options_other_dimensions, size=(200,20)))
            elif(self.nmrdata.phase==False and self.nmrdata.phase2==True):
                for i in range(3):
                    if(i==0):
                        self.carrier_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.references_proton[0]), size=(200,20)))
                        self.carrier_combo_boxes.append(wx.ComboBox(self, value=self.nmrdata.references_proton_labels[0], choices=self.options_proton, size=(200,20)))
                    if(i==1):
                            self.carrier_frequency_boxes.append(wx.TextCtrl(self, value=str(self.nmrdata.references_other[1]), size=(200,20)))
                            self.carrier_combo_boxes.append(wx.ComboBox(self, value=self.options_other_dimensions[1], choices=self.options_other_dimensions, size=(200,20)))
                    if(i==2):
                            self.carrier_frequency_boxes.append(wx.TextCtrl(self, value='0.0', size=(200,20)))
                            self.carrier_combo_boxes.append(wx.ComboBox(self, value='Other', choices=self.options_other_dimensions, size=(200,20)))

        self.carrier_frequency_sizer.AddSpacer(20)      
        self.carrier_frequency_sizer.Add(self.carrier_frequency_txt)
        self.carrier_frequency_sizer.AddSpacer(16)
        for i in range(0,len(self.carrier_frequency_boxes)):
            self.carrier_frequency_sizer.Add(self.carrier_frequency_boxes[i])
            self.carrier_frequency_sizer.AddSpacer(20)
        

        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.carrier_frequency_sizer)

        self.carrier_combo_box_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.carrier_combo_box_sizer.AddSpacer(20)
        self.carrier_combo_box_txt = wx.StaticText(self, label='Referencing Mode:          ')
        self.carrier_combo_box_sizer.Add(self.carrier_combo_box_txt)
        self.carrier_combo_box_sizer.AddSpacer(17)
        for i in range(0,len(self.carrier_combo_boxes)):
            if(self.carrier_combo_boxes[i] == 200):
                self.carrier_combo_box_sizer.AddSpacer(200)
                self.carrier_combo_box_sizer.AddSpacer(20)
            else:
                self.carrier_combo_box_sizer.Add(self.carrier_combo_boxes[i])
                self.carrier_combo_box_sizer.AddSpacer(20)
                self.carrier_combo_boxes[i].Bind(wx.EVT_COMBOBOX, self.on_carrier_combo_box_change_varian)
        
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.carrier_combo_box_sizer)

    
    def on_carrier_combo_box_change_varian(self, event):
        # Find the index of the combo box that was changed
        index = self.carrier_combo_boxes.index(event.GetEventObject())
        # Find the value of the combo box
        value = event.GetEventObject().GetValue()
        index_selection = self.carrier_combo_boxes[index].GetSelection()

        if(index==0):
            if(self.nmrdata.label_direct == '1H' or self.nmrdata.label_direct == 'H1'):
                if(value == 'H2O'):
                    self.carrier_frequency_boxes[index].SetValue(str(self.nmrdata.references_proton[0]))
                elif(value == 'Other'):
                    self.carrier_frequency_boxes[index].SetValue(str(self.nmrdata.references_proton[1]))
        else:
            if(self.nmrdata.label_direct == '1H' or self.nmrdata.label_direct == 'H1'):
                if(value == 'Other'):
                    self.carrier_frequency_boxes[index].SetValue(str(0.0))
                else:
                    self.carrier_frequency_boxes[index].SetValue(str(self.nmrdata.references_other[index_selection]))


    def acquisition_2D_mode_combo_box(self):
        self.acquisition_2D_mode_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.acquisition_2D_mode_txt = wx.StaticText(self, label='Acquisition Mode 2D:       ')
        self.acquisition_mode_2D_choices = ['States', 'TPPI', 'Magnitude', 'Real', 'Complex', 'Image']
        if(self.nmrdata.spectrometer == 'Bruker'):
            val = self.nmrdata.pseudo_flag
        else:
            val = self.nmrdata.other_params

        if(val == False):
            self.acquisition_2D_mode_box = wx.ComboBox(self, value=self.acquisition_mode_2D_choices[0], choices=self.acquisition_mode_2D_choices, size=(200,20))
        else:
            if(len(self.nmrdata.size_indirect)==1):
                self.acquisition_2D_mode_box = wx.ComboBox(self, value=self.acquisition_mode_2D_choices[3], choices=self.acquisition_mode_2D_choices, size=(200,20))
            else:
                self.acquisition_2D_mode_box = wx.ComboBox(self, value=self.acquisition_mode_2D_choices[0], choices=self.acquisition_mode_2D_choices, size=(200,20))

        self.acquisition_2D_mode_sizer.AddSpacer(20)
        self.acquisition_2D_mode_sizer.Add(self.acquisition_2D_mode_txt)
        self.acquisition_2D_mode_sizer.AddSpacer(11)
        self.acquisition_2D_mode_sizer.AddSpacer(220)
        self.acquisition_2D_mode_sizer.Add(self.acquisition_2D_mode_box)
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.acquisition_2D_mode_sizer)



    def create_temperature_box(self):
        self.extra_boxes_total = wx.BoxSizer(wx.VERTICAL)
        self.extra_boxes_0 = wx.BoxSizer(wx.HORIZONTAL)
        self.extra_boxes_0_total = wx.BoxSizer(wx.HORIZONTAL)
        self.extra_boxes = wx.BoxSizer(wx.HORIZONTAL)
        self.extra_boxes_total_1 = wx.BoxSizer(wx.HORIZONTAL)
        self.extra_boxes.AddSpacer(10)
        self.extra_boxes_0.AddSpacer(20)
        self.extra_boxes_0_total.Add(self.extra_boxes_0, 0, wx.CENTER)
        self.extra_boxes_total.Add(self.extra_boxes_0_total, 0, wx.CENTER)
        self.extra_boxes_total.AddSpacer(20)
        self.extra_boxes_total_1.Add(self.extra_boxes, 0, wx.CENTER)
        self.extra_boxes_total.Add(self.extra_boxes_total_1, 0, wx.CENTER)
        self.temperature_box_label = wx.StaticBox(self, -1, label = 'Experiment Temperature:')
        self.temperature_box = wx.StaticBoxSizer(self.temperature_box_label, wx.VERTICAL)
        self.text = wx.StaticText(self, label='Temp (K):')
        self.temperature_value = str(self.nmrdata.temperature)
        self.row_1 = wx.BoxSizer(wx.HORIZONTAL)
        self.temperature_parameter = wx.StaticText(self, label=self.temperature_value)
        self.row_1.AddSpacer(5)
        self.row_1.Add(self.text)
        self.row_1.AddSpacer(5)
        self.row_1.Add(self.temperature_parameter)
        self.temperature_input_button = wx.Button(self, label='Change Temperature')
        self.temperature_input_button.Bind(wx.EVT_BUTTON,self.on_temperature_change_button)
        self.temperature_box.Add(self.row_1)
        self.temperature_box.AddSpacer(10)
        self.temperature_box.Add(self.temperature_input_button)

        self.extra_boxes_0.AddSpacer(10)
        self.extra_boxes_0.Add(self.temperature_box)
        self.extra_boxes_0.AddSpacer(20)

        self.extra_sizers.AddSpacer(20)
        self.extra_sizers.Add(self.extra_boxes_total)



    def on_temperature_change_button(self,event):
        # Create a popout window to change the temperature
        self.temperature_change_window = wx.Frame(self, wx.ID_ANY, 'Change Temperature', wx.DefaultPosition, size=(300,120))
        # Change background colour depending of in darkmode or not
        if(darkdetect.isDark() == True and platform!='windows'):
            self.temperature_change_window.SetBackgroundColour((53, 53, 53, 255))
        else:
            self.temperature_change_window.SetBackgroundColour('White')
        self.temperature_change_window_sizer = wx.BoxSizer(wx.VERTICAL)
        self.temperature_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.temperature_label = wx.StaticText(self.temperature_change_window, label='New Temperature (K):')
        self.temperature_input = wx.TextCtrl(self.temperature_change_window, value=str(self.nmrdata.temperature), size=(100,20))
        self.temperature_sizer.AddSpacer(20)
        self.temperature_sizer.Add(self.temperature_label)
        self.temperature_sizer.AddSpacer(5)
        self.temperature_sizer.Add(self.temperature_input)
        self.temperature_change_window_sizer.AddSpacer(10)
        self.temperature_change_window_sizer.Add(self.temperature_sizer)
        self.temperature_change_window_sizer.AddSpacer(20)
        self.save_and_exit_temperature_change = wx.Button(self.temperature_change_window, label='Save', size=(100,20))
        self.save_and_exit_temperature_change.Bind(wx.EVT_BUTTON, self.on_save_and_exit_temperature_change)
        self.save_button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.save_button_sizer.AddSpacer(100)
        self.save_button_sizer.Add(self.save_and_exit_temperature_change, 0, wx.CENTER)
        self.save_button_sizer.AddSpacer(100)
        self.temperature_change_window_sizer.Add(self.save_button_sizer)
        self.temperature_change_window_sizer.AddSpacer(20)
        self.temperature_change_window.SetSizer(self.temperature_change_window_sizer)
        self.temperature_change_window.Show()


    def on_save_and_exit_temperature_change(self,event):
        self.nmrdata.temperature = float(self.temperature_input.GetValue())
        if(self.nmrdata.spectrometer == 'Bruker'):
            self.nmrdata.calculate_carrier_frequency_bruker()
        else:
            self.nmrdata.calculate_carrier_frequency_varian()
        
        self.temperature_parameter.SetLabel(str(self.nmrdata.temperature))
        self.temperature_change_window.Close()

        if(self.nmrdata.spectrometer=='Bruker'):
            # Get all the indexes of the current direct/indirect dimension carrier frequency boxes
            j=0
            for i in range(len(self.carrier_combo_boxes)):
                if(self.carrier_combo_boxes[i] == 200):
                    j += 1
                else:
                    label = self.carrier_combo_boxes[i].GetValue()
                    index = 0
                    
                    if(i==0):
                        for k in range(len(self.nmrdata.references_proton_labels)):
                            if(label == self.nmrdata.references_proton_labels[k]):
                                index = k
                                break
                        self.carrier_frequency_boxes[i].SetValue(str(self.nmrdata.references_proton[index]))

                    else:
                        for k in range(len(self.nmrdata.references_other_labels)):
                            if(label == self.nmrdata.references_other_labels[k]):
                                index = k
                                break
                        self.carrier_frequency_boxes[i].SetValue(str(self.nmrdata.references_other[index]))
        else:
            # Get all the indexes of the current direct/indirect dimension carrier frequency boxes
            j=0
            for i in range(len(self.carrier_combo_boxes)):
                if(self.carrier_combo_boxes[i].GetValue() == 'ID' or self.carrier_combo_boxes[i].GetValue() == 'Other' or self.carrier_combo_boxes[i].GetValue() == 'Manual' or self.acqusition_combo_boxes[i].GetValue() == 'Real'):
                    j += 1
                else:
                    label = self.carrier_combo_boxes[i].GetValue()
                    index = 0
                    
                    if(i==0):
                        for k in range(len(self.nmrdata.references_proton_labels)):
                            if(label == self.nmrdata.references_proton_labels[k]):
                                index = k
                                break
                        self.carrier_frequency_boxes[i].SetValue(str(self.nmrdata.references_proton[index]))

                    else:
                        for k in range(len(self.nmrdata.references_other_labels)):
                            if(label == self.nmrdata.references_other_labels[k]):
                                index = k
                                break
                        self.carrier_frequency_boxes[i].SetValue(str(self.nmrdata.references_other[index]))
        
  
        




    def find_nus_file(self):
        self.nuslist_found = False
        for file in os.listdir():
            if(file == 'nuslist'):
                self.nuslist_found = True
                self.nusfile = file
                self.include_NUS = True
                break
        
        if(self.nuslist_found == False):
            self.include_NUS = False
            self.nusfile = ''

        

        



    def input_NUS_list_box(self):
        self.NUS_label = wx.StaticBox(self,-1,label='NUS Information')
        self.NUS_box = wx.StaticBoxSizer(self.NUS_label,wx.VERTICAL)
        self.NUS_box_2 = wx.BoxSizer(wx.VERTICAL)

        # Tick box to see if user wants to perform NUS reconstruction
        self.NUS_tickbox = wx.CheckBox(self,-1,label='Include NUS Reconstruction')
        self.NUS_tickbox.Bind(wx.EVT_CHECKBOX, self.On_NUS_CheckBox)
        self.NUS_box.Add(self.NUS_tickbox)
        self.NUS_box_2 = wx.BoxSizer(wx.VERTICAL)

        if(self.include_NUS == True):
            self.NUS_offset = 0
            self.NUS_sample_count = 0
            self.nusfile = ''
            for file in os.listdir():
                if(file == 'nuslist'):
                    self.nuslist_found = True
                    self.nusfile = file
                    break
            if(self.nuslist_found == True):
                lines = open('nuslist','r').readlines()
                self.NUS_sample_count = len(lines)
                if(lines[0].split()[0] == '0'):
                    self.NUS_offset = 0
                else:
                    self.NUS_offset = 1
            

            self.NUS_tickbox.SetValue(True)

            if(self.nuslist_found == True):
                self.digital_filter_radio_box.SetSelection(0)


            self.nus_details_box()


        

        self.NUS_box.Add(self.NUS_box_2)


        self.extra_boxes.AddSpacer(20)
        self.extra_boxes.Add(self.NUS_box)


    def nus_details_box(self):
        self.NUS_sample_count_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.NUS_sample_count_txt = wx.StaticText(self, label='NUS Sample Count:')
        self.NUS_sample_count_box = wx.TextCtrl(self, value=str(self.NUS_sample_count), size=(50,20))
        self.NUS_sample_count_sizer.AddSpacer(5)
        self.NUS_sample_count_sizer.Add(self.NUS_sample_count_txt)
        self.NUS_sample_count_sizer.AddSpacer(5)
        self.NUS_sample_count_sizer.Add(self.NUS_sample_count_box)

        self.NUS_box_2.AddSpacer(10)
        self.NUS_box_2.Add(self.NUS_sample_count_sizer)

        self.nusfile_label = wx.StaticText(self, label="NUS schedule:")
        self.nusfile_input = wx.TextCtrl(self,value = self.nusfile,size=(120,20))
        self.find_nus_file = wx.Button(self, label="...", size=(25,20))
        self.find_nus_file.Bind(wx.EVT_BUTTON, lambda evt: self.on_find_nus_file(evt, self.nusfile))
        self.nusfile_input.SetValue(self.nusfile)



        self.nus_extras_box = wx.BoxSizer(wx.HORIZONTAL)
        self.nus_offset_label = wx.StaticText(self, label='NUS Offset:')
        self.nus_offset_box = wx.TextCtrl(self, value=str(self.NUS_offset), size=(30,20))
        self.nus_extras_box.AddSpacer(5)
        self.nus_extras_box.Add(self.nus_offset_label)
        self.nus_extras_box.AddSpacer(5)
        self.nus_extras_box.Add(self.nus_offset_box)

        self.reverse_NUS_tickbox = wx.CheckBox(self,-1,label='Reverse NUS Schedule')
        self.reverse_NUS_tickbox.SetValue(False)
        self.nus_extras_box.AddSpacer(10)
        self.nus_extras_box.Add(self.reverse_NUS_tickbox)


        

        # Put all the NUS sizers into one box
        self.nusfile_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.nusfile_sizer.AddSpacer(5)
        self.nusfile_sizer.Add(self.nusfile_label)
        self.nusfile_sizer.AddSpacer(5)
        self.nusfile_sizer.Add(self.nusfile_input)
        self.nusfile_sizer.AddSpacer(10)
        self.nusfile_sizer.Add(self.find_nus_file)
        self.nusfile_sizer.AddSpacer(5)

        self.NUS_box_2.AddSpacer(10)
        self.NUS_box_2.Add(self.nusfile_sizer)
        self.NUS_box_2.AddSpacer(10)
        self.NUS_box_2.Add(self.nus_extras_box)

    def on_find_nus_file(self, e, textBox):
        #get dialog box here
        cwd = os.getcwd()
        dlg = wx.FileDialog(self, message="Choose a file", defaultDir=cwd, defaultFile=self.nusfile, style=wx.FD_OPEN | wx.FD_CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            if(cwd in path):
                splitPath = path.split(cwd)
                textBox =  '.'+splitPath[1]
            else:
                textBox = path

            self.nusfile_input.SetValue(textBox)
            self.Layout()
            self.Refresh()

        dlg.Destroy()
        
        try:
            self.read_nus_file()
        except:
            # Give popout saying NUS file not read correctly
            dlg = wx.MessageDialog(self, 'NUS file not read correctly', 'Error', wx.OK | wx.ICON_ERROR)
            self.Raise()
            self.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            self.nusfile_input.SetValue('')
            self.nusfile = ''
            self.NUS_sample_count = 0
            self.NUS_sample_count_box.SetValue(str(self.NUS_sample_count))
            self.nus_offset_box.SetValue('0')
            self.NUS_offset = 0
            self.Layout()
            self.Refresh()


        


    def read_nus_file(self):
        # Read the nus file
        self.nusfile = self.nusfile_input.GetValue()
        lines = open(self.nusfile,'r').readlines()
        self.NUS_sample_count = len(lines)
        self.NUS_sample_count_box.SetValue(str(self.NUS_sample_count))
        if(lines[0].split()[0] == '0'):
            self.NUS_offset = 0
        else:
            self.NUS_offset = 1
        self.nus_offset_box.SetValue(str(self.NUS_offset))
        self.Layout()
        self.Refresh()


    def On_NUS_CheckBox(self, event):
        
        if(self.NUS_tickbox.GetValue() == True):
            self.include_NUS = True
            self.NUS_box_2.Clear(True)
            self.NUS_box.Clear(True)
            self.extra_boxes.Remove(self.NUS_box)
            self.extra_boxes.Detach(len(self.extra_boxes.GetChildren())-1)
            self.input_NUS_list_box()
            self.Layout()
            self.Refresh()
        else:
            self.include_NUS = False
            self.NUS_box_2.Clear(True)
            self.NUS_box.Clear(True)
            self.extra_boxes.Remove(self.NUS_box)
            self.extra_boxes.Detach(len(self.extra_boxes.GetChildren())-1)
            self.input_NUS_list_box()
            self.Layout()
            self.Refresh()
            
        


    def find_bruker_digital_filter_parameters(self):
        # Search through the acqus file to find the parameters
        acqus_file_lines = open(self.nmrdata.parameter_file,'r').readlines()

        # Find the decim, dspfvs and grpdly parameters
        self.decim = 0
        self.dspfvs = 0
        self.grpdly = 0
        try:
            for i in range(len(acqus_file_lines)):
                if('##$DECIM=' in acqus_file_lines[i]):
                    line = acqus_file_lines[i].split()
                    self.decim = float(line[1])
                if('##$DSPFVS=' in acqus_file_lines[i]):
                    line = acqus_file_lines[i].split()
                    self.dspfvs = int(line[1])
                if('##$GRPDLY=' in acqus_file_lines[i]):
                    line = acqus_file_lines[i].split()
                    self.grpdly = float(line[1])
            self.include_digital_filter = True

        except:
            self.decim = 0
            self.dspfvs = 0
            self.grpdly = 0
            self.include_digital_filter = False



    def create_bruker_digital_filter_box(self):
        self.digital_filter_box = wx.StaticBox(self, -1, label = 'Digital Filter')
        self.digital_filter_box_sizer_total = wx.StaticBoxSizer(self.digital_filter_box, wx.VERTICAL)
        self.digital_filter_box_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.digital_filter_box_sizer_total.AddSpacer(11)
        self.digital_filter_box_sizer_total.Add(self.digital_filter_box_sizer)
        self.digital_filter_box_sizer_total.AddSpacer(11)
        self.digital_filter_checkbox = wx.CheckBox(self, -1, label='Remove Digital Filter')
        self.digital_filter_checkbox.Bind(wx.EVT_CHECKBOX, self.on_digital_filter_checkbox)
        self.digital_filter_box_sizer.Add(self.digital_filter_checkbox)
        self.digital_filter_box_sizer.AddSpacer(10)
        if(self.include_digital_filter==True):
            self.digital_filter_checkbox.SetValue(True)
            self.create_bruker_digital_filter_options()
        else:
            self.digital_filter_box_sizer_total.AddSpacer(5)



        self.extra_boxes_0.Add(self.digital_filter_box_sizer_total)
 
        
        
        
        
    
    def create_bruker_digital_filter_options(self):
        self.decim_box = wx.BoxSizer(wx.HORIZONTAL)
        self.decim_text = wx.StaticText(self, label='Decimation Rate:')
        self.decim_value = str(self.decim)
        self.decim_textbox = wx.TextCtrl(self, value=str(self.decim_value))
        self.decim_box.AddSpacer(5)
        self.decim_box.Add(self.decim_text)
        self.decim_box.AddSpacer(5)
        self.decim_box.Add(self.decim_textbox)
        self.digital_filter_box_sizer.Add(self.decim_box)

        self.dspfvs_box = wx.BoxSizer(wx.HORIZONTAL)
        self.dspfvs_text = wx.StaticText(self, label='DSP Firmware Version:')
        self.dspfvs_value = str(self.dspfvs)
        self.dspfvs_textbox = wx.TextCtrl(self, value=str(self.dspfvs_value))
        self.dspfvs_box.AddSpacer(5)
        self.dspfvs_box.Add(self.dspfvs_text)
        self.dspfvs_box.AddSpacer(5)
        self.dspfvs_box.Add(self.dspfvs_textbox)
        self.digital_filter_box_sizer.AddSpacer(10)
        self.digital_filter_box_sizer.Add(self.dspfvs_box)

        self.grpdly_box = wx.BoxSizer(wx.HORIZONTAL)
        self.grpdly_text = wx.StaticText(self, label='Group Delay:')
        self.grpdly_value = str(self.grpdly)
        self.grpdly_textbox = wx.TextCtrl(self, value=str(self.grpdly_value))
        self.grpdly_box.AddSpacer(5)
        self.grpdly_box.Add(self.grpdly_text)
        self.grpdly_box.AddSpacer(5)
        self.grpdly_box.Add(self.grpdly_textbox)
        self.digital_filter_box_sizer.AddSpacer(10)
        self.digital_filter_box_sizer.Add(self.grpdly_box)

        # Have a radiobox to either remove the digital filter before processing or during processing
        self.digital_filter_radio_box = wx.RadioBox(self, choices=['Remove Before Processing', 'Remove During Processing'], majorDimension=1, style=wx.RA_SPECIFY_ROWS)
        self.digital_filter_radio_box.SetSelection(1)
        self.digital_filter_box_sizer_total.AddSpacer(10)
        self.digital_filter_box_sizer_total.Add(self.digital_filter_radio_box)



    def on_digital_filter_checkbox(self,event):
        if(self.digital_filter_checkbox.GetValue() == True):
            self.include_digital_filter = True
        else:
            self.include_digital_filter = False


        self.digital_filter_box_sizer.Clear(True)
        self.digital_filter_box_sizer_total.Clear(True)
        self.extra_boxes_0.Remove(self.digital_filter_box_sizer_total)
        # self.extra_boxes_0.Detach(len(self.extra_boxes_0.GetChildren())-1)
        self.create_bruker_digital_filter_box()
        self.Layout()
        self.Refresh()

    def find_bruker_scaling_parameters(self):
        # Search through the acqus file to find the NS and NC parameters
        acqus_file_lines = open(self.nmrdata.parameter_file,'r').readlines()

        # Find the NS and NC parameters
        try:
            for i in range(len(acqus_file_lines)):
                if('##$NS=' in acqus_file_lines[i]):
                    line = acqus_file_lines[i].split()
                    self.NS = int(line[1])
                if('##$NC=' in acqus_file_lines[i]):
                    line = acqus_file_lines[i].split()
                    self.NC = int(line[1])
            self.include_scaling = True
        except:
            self.NS = 0
            self.NC = 0
            self.include_scaling = False

    def find_varian_scaling_parameters(self):
        found_nt = False
        self.include_scaling = False
        for i,line in enumerate(self.nmrdata.procpar_file_lines):
            if(found_nt == True):
                self.NS = int(line.split()[1])
                self.include_scaling = True
                break
            elif(line.split()[0] == 'nt'):
                self.NS = int(line.split()[1])
                found_nt = True
        
            

    
    def create_intensity_scaling_box(self):
        if(self.nmrdata.spectrometer=='Bruker'):
            self.find_bruker_scaling_parameters()
        else:
            self.find_varian_scaling_parameters()
        self.scaling_box = wx.StaticBox(self, -1, label = 'Intensity Scaling')
        self.scaling_box_sizer_total = wx.StaticBoxSizer(self.scaling_box, wx.VERTICAL)
        self.scaling_box_sizer = wx.BoxSizer(wx.HORIZONTAL)
        if(self.nmrdata.spectrometer=='Bruker'):
            self.scaling_box_sizer_total.AddSpacer(2)
            self.scaling_box_sizer_total.Add(self.scaling_box_sizer)
            self.scaling_box_sizer_total.AddSpacer(2)
        else:
            self.scaling_box_sizer_total.AddSpacer(13)
            self.scaling_box_sizer_total.Add(self.scaling_box_sizer)
            self.scaling_box_sizer_total.AddSpacer(13)
        # Create tick boxes for intensity scaling
        self.scaling_NS_checkbox = wx.CheckBox(self, -1, label='1/NS')  # Normalise by number of scans
        if(self.nmrdata.spectrometer=='Bruker'):
            self.scaling_NC = wx.CheckBox(self, -1, label='2^NC')   # Normalise by bruker normalisation constant
        self.scaling_by_number = wx.CheckBox(self, -1, label='x1000')
        self.scaling_NS_checkbox.Bind(wx.EVT_CHECKBOX, self.on_scaling_checkbox)
        if(self.nmrdata.spectrometer=='Bruker'):
            self.scaling_NC.Bind(wx.EVT_CHECKBOX, self.on_scaling_checkbox)
        self.scaling_by_number.Bind(wx.EVT_CHECKBOX, self.on_scaling_checkbox)

        if(self.include_scaling == True):
            self.scaling_NS_checkbox.SetValue(True)
            if(self.nmrdata.spectrometer=='Bruker'):
                self.scaling_NC.SetValue(True)
        self.scaling_by_number.SetValue(True)
        if(self.nmrdata.spectrometer=='Bruker'):
            self.scaling_factor = (1/self.NS)*(2**self.NC)*1000
        else:
            self.scaling_factor = (1/self.NS)*1000
        self.scaling_text = wx.StaticText(self, label='Scaling Factor:')
        if(self.nmrdata.spectrometer=='Bruker'):
            self.scaling_number = wx.TextCtrl(self, value='{:.4E}'.format(self.scaling_factor), size=(50,20))
        else:
            self.scaling_number = wx.TextCtrl(self, value='{:.4E}'.format(self.scaling_factor), size=(80,20))
        self.scaling_box_sizer.Add(self.scaling_NS_checkbox)
        self.scaling_box_sizer.AddSpacer(20)
        if(self.nmrdata.spectrometer=='Bruker'):
            self.scaling_box_sizer.Add(self.scaling_NC)
            self.scaling_box_sizer.AddSpacer(20)
        self.scaling_box_sizer.Add(self.scaling_by_number)
        self.scaling_box_sizer.AddSpacer(20)
        self.scaling_box_sizer.Add(self.scaling_text)
        self.scaling_box_sizer.AddSpacer(5)
        self.scaling_box_sizer.Add(self.scaling_number)

        if(self.nmrdata.spectrometer=='Bruker'):
            self.bottom_left_box.AddSpacer(10)
            self.bottom_left_box.Add(self.scaling_box_sizer_total)
        else:
            self.extra_boxes_0.AddSpacer(20)
            self.extra_boxes_0.Add(self.scaling_box_sizer_total)


    def on_scaling_checkbox(self,event):

        value = 1
        if(self.scaling_NS_checkbox.GetValue() == True):
            value = value*(1/self.NS)
        if(self.nmrdata.spectrometer=='Bruker'):
            if(self.scaling_NC.GetValue() == True):
                value = value*(2**self.NC)
        if(self.scaling_by_number.GetValue() == True):
            value = value*1000
        
        self.scaling_factor = value
        self.scaling_number.SetValue('{:.4E}'.format(self.scaling_factor))
        

    def create_make_fid_conversion_box(self):
        # Have a button for make fid.com
        if(self.nmrdata.spectrometer=='Bruker'):
            self.fid_conversion_box = wx.BoxSizer(wx.VERTICAL)
        else:
            self.fid_conversion_box = wx.BoxSizer(wx.HORIZONTAL)

        self.make_fid_button = wx.Button(self, label='Make Conversion Script', size=(175,20))
        if(self.nmrdata.spectrometer=='Bruker'):
            self.make_fid_button.Bind(wx.EVT_BUTTON, self.on_make_fid_button_bruker)
        else:
            self.make_fid_button.Bind(wx.EVT_BUTTON, self.on_make_fid_button_varian)
        self.fid_conversion_box.AddSpacer(20)
        self.fid_conversion_box.Add(self.make_fid_button)

        # Have a button for show fid.com
        self.show_fid_button = wx.Button(self, label='Show Conversion Script', size=(175,20))
        self.show_fid_button.Bind(wx.EVT_BUTTON, self.on_show_fid_button)
        self.fid_conversion_box.AddSpacer(20)
        self.fid_conversion_box.Add(self.show_fid_button)

        # Have a button for convert 
        self.convert_button = wx.Button(self, label='Convert Data (nmrPipe)', size=(175,20))
        self.convert_button.Bind(wx.EVT_BUTTON, self.on_convert_button)
        self.fid_conversion_box.AddSpacer(20)
        self.fid_conversion_box.Add(self.convert_button)

        # Have a button for convert 
        self.convert_button2 = wx.Button(self, label='Convert Data (nmrglue)', size=(175,20))
        self.convert_button2.Bind(wx.EVT_BUTTON, self.on_convert_button2)
        self.fid_conversion_box.AddSpacer(20)
        self.fid_conversion_box.Add(self.convert_button2)

        self.extra_boxes.AddSpacer(20)
        self.extra_boxes.Add(self.fid_conversion_box)

    
    def determine_byte_order(self):
        try:
            with open(self.nmrdata.parameter_file) as file:
                lines = file.readlines()
                for line in lines:
                    if('##$BYTORDA') in line:
                        self.byte_order = line.split('\n')[0].split()[-1]
                        break
        except:
            self.byte_order = 0
        

    def determine_byte_size(self):
        try:
            with open(self.nmrdata.parameter_file) as file:
                lines = file.readlines()
                for line in lines:
                    if('##$DTYPA') in line:
                        self.d_type = line.split('\n')[0].split()[-1]
                        break
        except:
            self.d_type = 0


    def on_make_fid_button_bruker(self,event):

        # Check to see if any of the labels have duplicate names, if so, give a warning
        total_labels = []
        for i,label in enumerate(self.nucleus_type_boxes):
            total_labels.append(label.GetValue())
        if(len(total_labels) != len(set(total_labels))):
            dlg = wx.MessageDialog(self, 'Duplicate labels found. Please rename labels so each dimension has a different label, then try again.', 'Warning', wx.OK | wx.ICON_WARNING)
            self.Raise()
            self.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            return

        # Check to see if the fid.com file already exists
        if(self.file_parser == True):
            os.chdir(self.path)
        if('fid.com' in os.listdir()):
            dlg = wx.MessageDialog(self, 'fid.com already exists. Do you want to overwrite it?', 'Warning', wx.YES_NO | wx.ICON_WARNING)
            self.Raise()
            self.SetFocus()
            if(dlg.ShowModal() == wx.ID_NO):
                dlg.Destroy()
                if(self.file_parser == True):
                    os.chdir(self.cwd)
                return
            dlg.Destroy()
        # Create the fid.com file
        fid_file = open('fid.com','w')
        fid_file.write('#!/bin/csh\n\n')

        if(self.include_NUS == True):
            fid_file.write('nusExpand.tcl -mode bruker -sampleCount ' + str(self.NUS_sample_count) + ' -off ' + str(self.NUS_offset) + ' \\\n' + ' -in ./' +self.nmrdata.files[0] + ' -out ./ser_full -sample ' + self.nusfile)
            if(self.reverse_NUS_tickbox.GetValue() == True):
                fid_file.write(' -rev')

            fid_file.write('\n\n')

            spectype = './ser_full'
        
        else:
            spectype = self.nmrdata.files[0]

        if(int(self.byte_order)==1):
            self.byte_order_text = '-noaswap'
        else:
            self.byte_order_text = '-aswap'

        if(self.d_type=='2'):
            self.byte_size_text = ' -ws 8 -noi2f '
        else:
            self.byte_size_text = ''


        
        if(self.nmrdata.remove_acquisition_padding==True):
            if(self.include_digital_filter == False):
                fid_file.write('bruk2pipe -verb -in ' + spectype + ' \\\n   -bad ' + str(self.nmrdata.bad_point_threshold) + ' -ext ' + self.byte_order_text +  self.byte_size_text + '     \\\n')
            else:
                if(self.digital_filter_radio_box.GetSelection() == 0):
                    self.AMX_vs_DMX = '-DMX'
                else:
                    self.AMX_vs_DMX = '-AMX'
                fid_file.write('bruk2pipe -verb -in ' + spectype +  ' \\\n   -bad ' + str(self.nmrdata.bad_point_threshold) + ' -ext ' + self.byte_order_text + ' ' + self.AMX_vs_DMX + ' -decim ' + str(self.decim) + ' -dspfvs ' + str(self.dspfvs) + ' -grpdly ' + str(self.grpdly) + self.byte_size_text +   '\\\n')

        else:
            if(self.include_digital_filter == False):
                fid_file.write('bruk2pipe -verb -in ' + spectype + ' \\\n' +'   '+self.byte_order_text + self.byte_size_text + '\\\n')
            else:
                if(self.digital_filter_radio_box.GetSelection() == 0):
                    self.AMX_vs_DMX = '-DMX'
                else:
                    self.AMX_vs_DMX = '-AMX'
                fid_file.write('bruk2pipe -verb -in ' + spectype + ' \\\n'+ '   ' + self.byte_order_text + ' ' + self.AMX_vs_DMX + ' -decim ' + str(self.decim) + ' -dspfvs ' + str(self.dspfvs) + ' -grpdly ' + str(self.grpdly) + self.byte_size_text + ' \\\n')

        if(self.nmrdata.size_indirect==[]):
            file_array = []
            file_array.append(['-xN',self.N_complex_boxes[0].GetValue().strip(),' \\\n'])
            file_array.append(['-xT',self.N_real_boxes[0].GetValue().strip(),'  \\\n'])
            file_array.append(['-xMODE',self.acqusition_combo_boxes[0].GetValue().strip(),' \\\n'])
            file_array.append(['-xSW',self.sweep_width_boxes[0].GetValue().strip(),' \\\n'])
            file_array.append(['-xOBS',self.nuclei_frequency_boxes[0].GetValue().strip(),' \\\n'])
            file_array.append(['-xCAR', self.carrier_frequency_boxes[0].GetValue().strip(), '\\\n'])
            file_array.append(['-xLAB',self.nucleus_type_boxes[0].GetValue().strip(),' \\\n'])
            file_array.append(['-ndim', '1', ' \\\n'])
            for row in file_array:
                fid_file.write("{:>10} {:>20} {:>5}".format(*row))
            if(self.scaling_factor!=1):
                fid_file.write('| nmrPipe -fn MULT -c ' + str(self.scaling_factor) + ' \\\n')
            fid_file.write(' -ov -out ./test.fid\n')
        elif(len(self.nmrdata.size_indirect)==1):
            file_array = []
            file_array.append(['-xN',self.N_complex_boxes[0].GetValue().strip(),'-yN',self.N_complex_boxes[1].GetValue().strip(),'\\\n'])
            file_array.append(['-xT',self.N_real_boxes[0].GetValue().strip(),'-yT',self.N_real_boxes[1].GetValue().strip(),'\\\n'])
            file_array.append(['-xMODE',self.acqusition_combo_boxes[0].GetValue().strip(),'-yMODE',self.acqusition_combo_boxes[1].GetValue().strip(),'\\\n'])
            file_array.append(['-xSW',self.sweep_width_boxes[0].GetValue().strip(),'-ySW',self.sweep_width_boxes[1].GetValue().strip(),'\\\n'])
            file_array.append(['-xOBS',self.nuclei_frequency_boxes[0].GetValue().strip(),'-yOBS',self.nuclei_frequency_boxes[1].GetValue().strip(),'\\\n'])
            file_array.append(['-xCAR', self.carrier_frequency_boxes[0].GetValue().strip(), '-yCAR', self.carrier_frequency_boxes[1].GetValue().strip(), '\\\n'])
            file_array.append(['-xLAB',self.nucleus_type_boxes[0].GetValue().strip(),'-yLAB',self.nucleus_type_boxes[1].GetValue().strip(),'\\\n'])
            file_array.append(['-ndim', '2', '-aq2D',self.acquisition_2D_mode_box.GetValue().strip(),'\\\n'])
            for row in file_array:
                fid_file.write("{:>10} {:>20} {:>10} {:>20} {:>5}".format(*row))
            if(self.scaling_factor!=1):
                fid_file.write('| nmrPipe -fn MULT -c ' + str(self.scaling_factor) + ' \\\n')
            fid_file.write(' -ov -out ./test.fid\n')
        elif(len(self.nmrdata.size_indirect)==2):
            file_array = []
            file_array.append(['-xN',self.N_complex_boxes[0].GetValue().strip(),'-yN',self.N_complex_boxes[1].GetValue().strip(),'-zN',self.N_complex_boxes[2].GetValue().strip(),'\\\n'])
            file_array.append(['-xT',self.N_real_boxes[0].GetValue().strip(),'-yT',self.N_real_boxes[1].GetValue().strip(),'-zT',self.N_real_boxes[2].GetValue().strip(),'\\\n'])
            file_array.append(['-xMODE',self.acqusition_combo_boxes[0].GetValue().strip(),'-yMODE',self.acqusition_combo_boxes[1].GetValue().strip(),'-zMODE',self.acqusition_combo_boxes[2].GetValue().strip(),'\\\n'])
            file_array.append(['-xSW',self.sweep_width_boxes[0].GetValue().strip(),'-ySW',self.sweep_width_boxes[1].GetValue().strip(),'-zSW',self.sweep_width_boxes[2].GetValue().strip(),'\\\n'])
            file_array.append(['-xOBS',self.nuclei_frequency_boxes[0].GetValue().strip(),'-yOBS',self.nuclei_frequency_boxes[1].GetValue().strip(),'-zOBS',self.nuclei_frequency_boxes[2].GetValue().strip(),'\\\n'])
            file_array.append(['-xCAR', self.carrier_frequency_boxes[0].GetValue().strip(), '-yCAR', self.carrier_frequency_boxes[1].GetValue().strip(), '-zCAR', self.carrier_frequency_boxes[2].GetValue().strip(), '\\\n'])
            file_array.append(['-xLAB',self.nucleus_type_boxes[0].GetValue().strip(),'-yLAB',self.nucleus_type_boxes[1].GetValue().strip(),'-zLAB',self.nucleus_type_boxes[2].GetValue().strip(),'\\\n'])
            for row in file_array:
                fid_file.write("{:>10} {:>20} {:>10} {:>20} {:>10} {:>20} {:>5}".format(*row))

            
            row = ['-ndim', '3','-aq2D',self.acquisition_2D_mode_box.GetValue().strip(),'', '', '\\\n']
            fid_file.write("{:>10} {:>20} {:>10} {:>20} {:>10} {:>20} {:>5}".format(*row))
            if(self.scaling_factor!=1):
                fid_file.write('| nmrPipe -fn MULT -c ' + str(self.scaling_factor) + ' \\\n')
            fid_file.write(' -ov -out ./test.fid\n')


        fid_file.close()

        if(self.file_parser == True):
            os.chdir(self.cwd)
    
    def on_make_fid_button_varian(self,event):
         # Check to see if any of the labels have duplicate names, if so, give a warning
        total_labels = []
        for i,label in enumerate(self.nucleus_type_boxes):
            total_labels.append(label.GetValue())
        if(len(total_labels) != len(set(total_labels))):
            dlg = wx.MessageDialog(self, 'Duplicate labels found. Please rename labels so each dimension has a different label, then try again.', 'Warning', wx.OK | wx.ICON_WARNING)
            self.Raise()
            self.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            return
        if(self.file_parser == True):
            os.chdir(self.path)
        # Check to see if the fid.com file already exists
        if('fid.com' in os.listdir()):
            dlg = wx.MessageDialog(self, 'fid.com already exists. Do you want to overwrite it?', 'Warning', wx.YES_NO | wx.ICON_WARNING)
            self.Raise()
            self.SetFocus()
            if(dlg.ShowModal() == wx.ID_NO):
                dlg.Destroy()
                if(self.file_parser == True):
                    os.chdir(self.cwd)
                return
            dlg.Destroy()
        
        # Create the fid.com file for the varian data
        fid_file = open('fid.com','w')
        fid_file.write('#!/bin/csh\n\n')

        if(self.nmrdata.reorder_array==True):
            # Need to run RelaxFix.out to reorder the data
            # Give a message to the user saying that data reshuffling is necessary and that RelaxFix.out will be run, make sure RelaxFix.out is in the correct path
            dlg = wx.MessageDialog(self, 'Data reshuffling is necessary. RelaxFix.out will be run to reorder the data', 'Warning', wx.YES_NO | wx.ICON_WARNING)
            self.Raise()
            self.SetFocus()
            if(dlg.ShowModal() == wx.ID_NO):
                dlg.Destroy()
                if(self.file_parser == True):
                    os.chdir(self.cwd)
                return
            # Try the RelaxFix.out command, if it fails, give a message to the user saying RelaxFix.out not found/failed to run and exit
            try:
                p = subprocess.Popen('RelaxFix.out',stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                p.wait()
                out, err = p.communicate()
                if(out==['']):
                    dlg = wx.MessageDialog(self, 'RelaxFix.out not found or failed to run', 'Error', wx.OK | wx.ICON_ERROR)
                    self.Raise()
                    self.SetFocus()
                    dlg.ShowModal()
                    dlg.Destroy()
                    if(self.file_parser == True):
                        os.chdir(self.cwd)
                    return
                else:
                    pass
            except:
                dlg = wx.MessageDialog(self, 'RelaxFix.out not found or failed to run', 'Error', wx.OK | wx.ICON_ERROR)
                self.Raise()
                self.SetFocus()
                dlg.ShowModal()
                dlg.Destroy()
                if(self.file_parser == True):
                    os.chdir(self.cwd)
                return

            fid_file.write('RelaxFix.out {} {} {} 0 fid.final {} \n\n'.format(int(self.nmrdata.size_direct),int(self.nmrdata.size_indirect[0]),int(self.nmrdata.number_of_arrayed_parameters),self.nmrdata.files[0]))

        if(self.include_NUS == True):
            fid_file.write('nusExpand.tcl -mode varian -sampleCount ' + str(self.NUS_sample_count) + ' -off ' + str(self.NUS_offset) + ' \\\n' + ' -in ./' + self.nmrdata.files[0]+ ' -out ./fid_full -sample ' + self.nusfile)
            if(self.reverse_NUS_tickbox.GetValue() == True):
                fid_file.write(' -rev')

            fid_file.write('\n\n')

            spectype = './fid_full'
        
        else:
            spectype = self.nmrdata.files[0]

        if(self.nmrdata.reorder_array==True):
            spectype = './fid.final'

        if(self.nmrdata.reverse_acquisition_order==False):
            fid_file.write('var2pipe -verb -in ' + spectype + ' \\\n' + '-noaswap' + '\\\n')
        else:
            fid_file.write('var2pipe -verb -in ' + spectype + ' \\\n' + '-noaswap -aqORD 1'+ '\\\n')

        if(self.nmrdata.other_params==False):
            try:
                self.nmrdata.size_indirect.remove(1)
            except:
                pass
            try:
                self.nmrdata.size_indirect.remove(1)
            except:
                pass
            if(self.nmrdata.size_indirect==[]):
                file_array = []
                file_array.append(['-xN',self.N_complex_boxes[0].GetValue().strip(),' \\\n'])
                file_array.append(['-xT',self.N_real_boxes[0].GetValue().strip(),'  \\\n'])
                file_array.append(['-xMODE',self.acqusition_combo_boxes[0].GetValue().strip(),' \\\n'])
                file_array.append(['-xSW',self.sweep_width_boxes[0].GetValue().strip(),' \\\n'])
                file_array.append(['-xOBS',self.nuclei_frequency_boxes[0].GetValue().strip(),' \\\n'])
                file_array.append(['-xCAR', self.carrier_frequency_boxes[0].GetValue().strip(), '\\\n'])
                file_array.append(['-xLAB',self.nucleus_type_boxes[0].GetValue().strip(),' \\\n'])
                file_array.append(['-ndim', '1', ' \\\n'])
                for row in file_array:
                    fid_file.write("{:>10} {:>20} {:>5}".format(*row))
                if(self.scaling_factor!=1):
                    fid_file.write('| nmrPipe -fn MULT -c ' + str(self.scaling_factor) + ' \\\n')
                fid_file.write(' -ov -out ./test.fid\n')
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==False):
                file_array = []
                file_array.append(['-xN',self.N_complex_boxes[0].GetValue().strip(),'-yN',self.N_complex_boxes[1].GetValue().strip(),'\\\n'])
                file_array.append(['-xT',self.N_real_boxes[0].GetValue().strip(),'-yT',self.N_real_boxes[1].GetValue().strip(),'\\\n'])
                file_array.append(['-xMODE',self.acqusition_combo_boxes[0].GetValue().strip(),'-yMODE',self.acqusition_combo_boxes[1].GetValue().strip(),'\\\n'])
                file_array.append(['-xSW',self.sweep_width_boxes[0].GetValue().strip(),'-ySW',self.sweep_width_boxes[1].GetValue().strip(),'\\\n'])
                file_array.append(['-xOBS',self.nuclei_frequency_boxes[0].GetValue().strip(),'-yOBS',self.nuclei_frequency_boxes[1].GetValue().strip(),'\\\n'])
                file_array.append(['-xCAR', self.carrier_frequency_boxes[0].GetValue().strip(), '-yCAR', self.carrier_frequency_boxes[1].GetValue().strip(), '\\\n'])
                file_array.append(['-xLAB',self.nucleus_type_boxes[0].GetValue().strip(),'-yLAB',self.nucleus_type_boxes[1].GetValue().strip(),'\\\n'])
                file_array.append(['-ndim', '2', '-aq2D',self.acquisition_2D_mode_box.GetValue().strip(),'\\\n'])
                for row in file_array:
                    fid_file.write("{:>10} {:>20} {:>10} {:>20} {:>5}".format(*row))
                if(self.scaling_factor!=1):
                    fid_file.write('| nmrPipe -fn MULT -c ' + str(self.scaling_factor) + ' \\\n')
                fid_file.write(' -ov -out ./test.fid\n')
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==True):
                file_array = []
                # print(self.nucleus_type_boxes[0].GetValue().strip())
                # print(self.nucleus_type_boxes[1].GetValue().strip())
                # print(self.nucleus_type_boxes[2].GetValue().strip())
                file_array.append(['-xN',self.N_complex_boxes[0].GetValue().strip(),'-yN',self.N_complex_boxes[1].GetValue().strip(),'-zN',self.N_complex_boxes[2].GetValue().strip(),'\\\n'])
                file_array.append(['-xT',self.N_real_boxes[0].GetValue().strip(),'-yT',self.N_real_boxes[1].GetValue().strip(),'-zT',self.N_real_boxes[2].GetValue().strip(),'\\\n'])
                file_array.append(['-xMODE',self.acqusition_combo_boxes[0].GetValue().strip(),'-yMODE',self.acqusition_combo_boxes[1].GetValue().strip(),'-zMODE',self.acqusition_combo_boxes[2].GetValue().strip(),'\\\n'])
                file_array.append(['-xSW',self.sweep_width_boxes[0].GetValue().strip(),'-ySW',self.sweep_width_boxes[1].GetValue().strip(),'-zSW',self.sweep_width_boxes[2].GetValue().strip(),'\\\n'])
                file_array.append(['-xOBS',self.nuclei_frequency_boxes[0].GetValue().strip(),'-yOBS',self.nuclei_frequency_boxes[1].GetValue().strip(),'-zOBS',self.nuclei_frequency_boxes[2].GetValue().strip(),'\\\n'])
                file_array.append(['-xCAR', self.carrier_frequency_boxes[0].GetValue().strip(), '-yCAR', self.carrier_frequency_boxes[1].GetValue().strip(), '-zCAR', self.carrier_frequency_boxes[2].GetValue().strip(), '\\\n'])
                file_array.append(['-xLAB',self.nucleus_type_boxes[0].GetValue().strip(),'-yLAB',self.nucleus_type_boxes[1].GetValue().strip(),'-zLAB',self.nucleus_type_boxes[2].GetValue().strip(),'\\\n'])
                for row in file_array:
                    fid_file.write("{:>10} {:>20} {:>10} {:>20} {:>10} {:>20} {:>5}".format(*row))

                
                row = ['-ndim', '3','-aq2D',self.acquisition_2D_mode_box.GetValue().strip(),'', '', '\\\n']
                fid_file.write("{:>10} {:>20} {:>10} {:>20} {:>10} {:>20} {:>5}".format(*row))
                if(self.scaling_factor!=1):
                    fid_file.write('| nmrPipe -fn MULT -c ' + str(self.scaling_factor) + ' \\\n')
                fid_file.write(' -ov -out ./test.fid\n')
            elif(self.nmrdata.phase==False and self.nmrdata.phase2==True):
                file_array = []
                file_array.append(['-xN',self.N_complex_boxes[0].GetValue().strip(),'-yN',self.N_complex_boxes[1].GetValue().strip(),'\\\n'])
                file_array.append(['-xT',self.N_real_boxes[0].GetValue().strip(),'-yT',self.N_real_boxes[1].GetValue().strip(),'\\\n'])
                file_array.append(['-xMODE',self.acqusition_combo_boxes[0].GetValue().strip(),'-yMODE',self.acqusition_combo_boxes[1].GetValue().strip(),'\\\n'])
                file_array.append(['-xSW',self.sweep_width_boxes[0].GetValue().strip(),'-ySW',self.sweep_width_boxes[1].GetValue().strip(),'\\\n'])
                file_array.append(['-xOBS',self.nuclei_frequency_boxes[0].GetValue().strip(),'-yOBS',self.nuclei_frequency_boxes[1].GetValue().strip(),'\\\n'])
                file_array.append(['-xCAR', self.carrier_frequency_boxes[0].GetValue().strip(), '-yCAR', self.carrier_frequency_boxes[1].GetValue().strip(), '\\\n'])
                file_array.append(['-xLAB',self.nucleus_type_boxes[0].GetValue().strip(),'-yLAB',self.nucleus_type_boxes[1].GetValue().strip(),'\\\n'])
                file_array.append(['-ndim', '2', '-aq2D',self.acquisition_2D_mode_box.GetValue().strip(),'\\\n'])
                for row in file_array:
                    fid_file.write("{:>10} {:>20} {:>10} {:>20} {:>5}".format(*row))
                if(self.scaling_factor!=1):
                    fid_file.write('| nmrPipe -fn MULT -c ' + str(self.scaling_factor) + ' \\\n')
                fid_file.write(' -ov -out ./test.fid\n')
        else:
            if(self.nmrdata.phase==False and self.nmrdata.phase2==False):
                file_array = []
                file_array.append(['-xN',self.N_complex_boxes[0].GetValue().strip(),'-yN',self.N_complex_boxes[1].GetValue().strip(),'\\\n'])
                file_array.append(['-xT',self.N_real_boxes[0].GetValue().strip(),'-yT',self.N_real_boxes[1].GetValue().strip(),'\\\n'])
                file_array.append(['-xMODE',self.acqusition_combo_boxes[0].GetValue().strip(),'-yMODE',self.acqusition_combo_boxes[1].GetValue().strip(),'\\\n'])
                file_array.append(['-xSW',self.sweep_width_boxes[0].GetValue().strip(),'-ySW',self.sweep_width_boxes[1].GetValue().strip(),'\\\n'])
                file_array.append(['-xOBS',self.nuclei_frequency_boxes[0].GetValue().strip(),'-yOBS',self.nuclei_frequency_boxes[1].GetValue().strip(),'\\\n'])
                file_array.append(['-xCAR', self.carrier_frequency_boxes[0].GetValue().strip(), '-yCAR', self.carrier_frequency_boxes[1].GetValue().strip(), '\\\n'])
                file_array.append(['-xLAB',self.nucleus_type_boxes[0].GetValue().strip(),'-yLAB',self.nucleus_type_boxes[1].GetValue().strip(),'\\\n'])
                file_array.append(['-ndim', '2', '-aq2D','Real','\\\n'])
                for row in file_array:
                    fid_file.write("{:>10} {:>20} {:>10} {:>20} {:>5}".format(*row))
                if(self.scaling_factor!=1):
                    fid_file.write('| nmrPipe -fn MULT -c ' + str(self.scaling_factor) + ' \\\n')
                fid_file.write(' -ov -out ./test.fid\n')
            elif(self.nmrdata.phase==True and self.nmrdata.phase2==False):
                file_array = []
                file_array.append(['-xN',self.N_complex_boxes[0].GetValue().strip(),'-yN',self.N_complex_boxes[1].GetValue().strip(),'-zN',self.N_complex_boxes[2].GetValue().strip(),'\\\n'])
                file_array.append(['-xT',self.N_real_boxes[0].GetValue().strip(),'-yT',self.N_real_boxes[1].GetValue().strip(),'-zT',self.N_real_boxes[2].GetValue().strip(),'\\\n'])
                file_array.append(['-xMODE',self.acqusition_combo_boxes[0].GetValue().strip(),'-yMODE',self.acqusition_combo_boxes[1].GetValue().strip(),'-zMODE',self.acqusition_combo_boxes[2].GetValue().strip(),'\\\n'])
                file_array.append(['-xSW',self.sweep_width_boxes[0].GetValue().strip(),'-ySW',self.sweep_width_boxes[1].GetValue().strip(),'-zSW',self.sweep_width_boxes[2].GetValue().strip(),'\\\n'])
                file_array.append(['-xOBS',self.nuclei_frequency_boxes[0].GetValue().strip(),'-yOBS',self.nuclei_frequency_boxes[1].GetValue().strip(),'-zOBS',self.nuclei_frequency_boxes[2].GetValue().strip(),'\\\n'])
                file_array.append(['-xCAR', self.carrier_frequency_boxes[0].GetValue().strip(), '-yCAR', self.carrier_frequency_boxes[1].GetValue().strip(), '-zCAR', self.carrier_frequency_boxes[2].GetValue().strip(), '\\\n'])
                file_array.append(['-xLAB',self.nucleus_type_boxes[0].GetValue().strip(),'-yLAB',self.nucleus_type_boxes[1].GetValue().strip(),'-zLAB',self.nucleus_type_boxes[2].GetValue().strip(),'\\\n'])
                file_array.append(['-ndim', '3', '-aq2D',self.acquisition_2D_mode_box.GetValue().strip(),'','','\\\n'])
                for row in file_array:
                    fid_file.write("{:>10} {:>20} {:>10} {:>20} {:>10} {:>20} {:>5}".format(*row))
                if(self.scaling_factor!=1):
                    fid_file.write('| nmrPipe -fn MULT -c ' + str(self.scaling_factor) + ' \\\n')
                fid_file.write(' -ov -out ./test.fid\n')
                self.threeDprocessing = True

            elif(self.nmrdata.phase==True and self.nmrdata.phase2==True):
                file_array = []
                file_array.append(['-xN',self.N_complex_boxes[0].GetValue().strip(),'-yN',self.N_complex_boxes[1].GetValue().strip(),'-zN',self.N_complex_boxes[2].GetValue().strip(),'-aN',self.N_complex_boxes[3].GetValue().strip(),'\\\n'])
                file_array.append(['-xT',self.N_real_boxes[0].GetValue().strip(),'-yT',self.N_real_boxes[1].GetValue().strip(),'-zT',self.N_real_boxes[2].GetValue().strip(),'-aT',self.N_real_boxes[3].GetValue.strip(),'\\\n'])
                file_array.append(['-xMODE',self.acqusition_combo_boxes[0].GetValue().strip(),'-yMODE',self.acqusition_combo_boxes[1].GetValue().strip(),'-zMODE',self.acqusition_combo_boxes[2].GetValue().strip(),'-aMODE',self.acqusition_combo_boxes[3].GetValue().strip(),'\\\n'])
                file_array.append(['-xSW',self.sweep_width_boxes[0].GetValue().strip(),'-ySW',self.sweep_width_boxes[1].GetValue().strip(),'-zSW',self.sweep_width_boxes[2].GetValue().strip(),'-aSW',self.sweep_width_boxes[3].GetValue().strip(),'\\\n'])
                file_array.append(['-xOBS',self.nuclei_frequency_boxes[0].GetValue().strip(),'-yOBS',self.nuclei_frequency_boxes[1].GetValue().strip(),'-zOBS',self.nuclei_frequency_boxes[2].GetValue().strip(),'-aOBS',self.nuclei_frequency_boxes[3].GetValue().strip(),'\\\n'])
                file_array.append(['-xCAR', self.carrier_frequency_boxes[0].GetValue().strip(), '-yCAR', self.carrier_frequency_boxes[1].GetValue().strip(), '-zCAR', self.carrier_frequency_boxes[2].GetValue().strip(), '-aCAR', self.carrier_frequency_boxes[3].GetValue().strip(), '\\\n'])
                file_array.append(['-xLAB',self.nucleus_type_boxes[0].GetValue().strip(),'-yLAB',self.nucleus_type_boxes[1].GetValue().strip(),'-zLAB',self.nucleus_type_boxes[2].GetValue().strip(),'-aLAB',self.nucleus_type_boxes[3].GetValue().strip(),'\\\n'])
                file_array.append(['-ndim', '4', '-aq2D',self.acquisition_2D_mode_box.GetValue().strip(),'','','','','\\\n'])
                for row in file_array:
                    fid_file.write("{:>10} {:>20} {:>10} {:>20} {:>10} {:>20} {:>10} {:>20} {:>5}".format(*row))
                if(self.scaling_factor!=1):
                    fid_file.write('| nmrPipe -fn MULT -c ' + str(self.scaling_factor) + ' \\\n')
                fid_file.write(' -ov -out ./test.fid\n')
                self.threeDprocessing = True
            
            elif(self.nmrdata.phase==False and self.nmrdata.phase2==True):
                file_array = []
                file_array.append(['-xN',self.N_complex_boxes[0].GetValue().strip(),'-yN',self.N_complex_boxes[1].GetValue().strip(),'-zN',self.N_complex_boxes[2].GetValue().strip(),'\\\n'])
                file_array.append(['-xT',self.N_real_boxes[0].GetValue().strip(),'-yT',self.N_real_boxes[1].GetValue().strip(),'-zT',self.N_real_boxes[2].GetValue().strip(),'\\\n'])
                file_array.append(['-xMODE',self.acqusition_combo_boxes[0].GetValue().strip(),'-yMODE',self.acqusition_combo_boxes[1].GetValue().strip(),'-zMODE',self.acqusition_combo_boxes[2].GetValue().strip(),'\\\n'])
                file_array.append(['-xSW',self.sweep_width_boxes[0].GetValue().strip(),'-ySW',self.sweep_width_boxes[1].GetValue().strip(),'-zSW',self.sweep_width_boxes[2].GetValue().strip(),'\\\n'])
                file_array.append(['-xOBS',self.nuclei_frequency_boxes[0].GetValue().strip(),'-yOBS',self.nuclei_frequency_boxes[1].GetValue().strip(),'-zOBS',self.nuclei_frequency_boxes[2].GetValue().strip(),'\\\n'])
                file_array.append(['-xCAR', self.carrier_frequency_boxes[0].GetValue().strip(), '-yCAR', self.carrier_frequency_boxes[1].GetValue().strip(), '-zCAR', self.carrier_frequency_boxes[2].GetValue().strip(), '\\\n'])
                file_array.append(['-xLAB',self.nucleus_type_boxes[0].GetValue().strip(),'-yLAB',self.nucleus_type_boxes[1].GetValue().strip(),'-zLAB',self.nucleus_type_boxes[2].GetValue().strip(),'\\\n'])
                file_array.append(['-ndim', '3', '-aq2D',self.acquisition_2D_mode_box.GetValue().strip(),'','','\\\n'])
                for row in file_array:
                    fid_file.write("{:>10} {:>20} {:>10} {:>20} {:>10} {:>20} {:>5}".format(*row))
                if(self.scaling_factor!=1):
                    fid_file.write('| nmrPipe -fn MULT -c ' + str(self.scaling_factor) + ' \\\n')
                fid_file.write(' -ov -out ./test.fid\n')
                self.threeDprocessing = True
        fid_file.close()

        if(self.file_parser == True):
            os.chdir(self.cwd)


        



    def on_show_fid_button(self,event):
        # Open up a new window containing a textcontrol box that a user can manually edit 
        self.show_fid_frame = wx.Frame(self,wx.ID_ANY,'fid.com',wx.DefaultPosition, size = (700,350))
        # Create the main sizer
        self.main_sizer_fid_com=wx.BoxSizer(wx.VERTICAL)
        


        if(darkdetect.isDark() == True and platform!='windows'):
            self.show_fid_frame.SetBackgroundColour((53, 53, 53, 255))
            matplotlib.rc('axes',edgecolor='white')
            matplotlib.rc('xtick',color='white')
            matplotlib.rc('ytick',color='white')
            matplotlib.rc('axes',labelcolor='white')
            matplotlib.rc('axes',facecolor='#282A36')
            matplotlib.rc('legend', labelcolor='white')


        else:
            self.show_fid_frame.SetBackgroundColour('White')
            matplotlib.rc('axes',edgecolor='black')
            matplotlib.rc('xtick',color='black')
            matplotlib.rc('ytick',color='black')
            matplotlib.rc('axes',labelcolor='black')
            matplotlib.rc('axes',facecolor='white')
            matplotlib.rc('legend', labelcolor='black')

        
        self.fid_com_sizer = wx.BoxSizer(wx.HORIZONTAL)

        if(self.file_parser == True):
            # print(self.path)
            # print(os.getcwd())
            os.chdir(self.path)

        try:
            fid_file = open('fid.com','r').readlines()
        except:
            dlg = wx.MessageDialog(self, 'The fid.com file does not exist. Please create the fid.com file first by clicking on the make conversion script button.', 'Warning', wx.OK | wx.ICON_WARNING)
            self.Raise()
            self.SetFocus()
            dlg.ShowModal()
            if(self.file_parser == True):
                os.chdir(self.cwd)
            dlg.Destroy()
            return
        


        self.fid_text = wx.StaticText(self.show_fid_frame, id=-1,label = ''.join(fid_file), size = (680,300))
        
        self.fid_com_sizer.Add(self.fid_text, wx.CENTER | wx.GROW)
        self.main_sizer_fid_com.AddSpacer(5)
        self.main_sizer_fid_com.Add(self.fid_com_sizer, wx.CENTER | wx.GROW)
        self.main_sizer_fid_com.AddSpacer(5)

        self.show_fid_frame.SetSizer(self.main_sizer_fid_com)
        self.show_fid_frame.Show()

        if(self.file_parser == True):
            os.chdir(self.cwd)


        
    def on_convert_button2(self,event):
        self.on_convert_button_total(conversion_platform='windows')


    def on_convert_button(self,event):
        if(platform=='windows'):
            # Outputting a message saying that nmrPipe conversion is not possible on windows
            dlg = wx.MessageDialog(self, 'nmrPipe not installed. Please use nmrglue convert button instead.', 'Warning', wx.OK | wx.ICON_WARNING)
            self.Raise()
            self.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            return
        else:
            self.on_convert_button_total(conversion_platform='mac')
        

    def on_convert_button_total(self,conversion_platform):
        # Check to see if the fid.com file exists
        if(self.file_parser == True):
            os.chdir(self.path)
        if('fid.com' not in os.listdir()):
            dlg = wx.MessageDialog(self, 'The fid.com file does not exist. Please create the fid.com file first by clicking on the make conversion script button.', 'Warning', wx.OK | wx.ICON_WARNING)
            self.Raise()
            self.SetFocus()
            dlg.ShowModal()
            if(self.file_parser == True):
                os.chdir(self.cwd)
            dlg.Destroy()
            return
        
        # Check to see if pipe2xyz in fid.com
        f = open('fid.com','r')
        lines = f.readlines()
        f.close()
        found_pipe2xyz = False
        for line in lines:
            if('pipe2xyz' in line):
                found_pipe2xyz = True

    
        # Check to see if test.fid already exists in the current directory and ask the user if they want to overwrite it
        if(found_pipe2xyz==False):
            if('test.fid' in os.listdir()):
                dlg = wx.MessageDialog(self, 'The test.fid file already exists. Do you want to overwrite it?', 'Warning', wx.YES_NO | wx.ICON_WARNING)
                self.Raise()
                self.SetFocus()
                if(dlg.ShowModal() == wx.ID_NO):
                    dlg.Destroy()
                    if(self.file_parser == True):
                        os.chdir(self.cwd)
                    return
                dlg.Destroy()
        else:
            if('fids' in os.listdir()):
                dlg = wx.MessageDialog(self, 'The \'fids\' directory already exists. Do you want to overwrite it?', 'Warning', wx.YES_NO | wx.ICON_WARNING)
                self.Raise()
                self.SetFocus()
                if(dlg.ShowModal() == wx.ID_NO):
                    dlg.Destroy()
                    if(self.file_parser == True):
                        os.chdir(self.cwd)
                    return
                dlg.Destroy()


        


        if(conversion_platform == 'windows'):
            print('converting nmrglue')
            # Check to see if the NUS flag is ticked
            if(self.include_NUS == True):
                # Give an error saying that NUS conversion is not currently supported on windows, please use a linux/mac containing nmrPipe to complete conversion
                dlg = wx.MessageDialog(self, 'NUS conversion is not currently supported on windows or using nmrglue conversion. Please use a linux/mac containing nmrPipe to complete conversion', 'Warning', wx.OK | wx.ICON_WARNING)
                self.Raise()
                self.SetFocus()
                dlg.ShowModal()
                if(self.file_parser == True):
                    os.chdir(self.cwd)
                dlg.Destroy()
                return
            
            # Unable to perform nmrPipe commands on windows 
            self.perform_conversion_windows()
            if(self.file_parser == True):
                os.chdir(self.cwd)
        else:
            print('converting nmrPipe')
            # Add the necessary permissions to the fid.com file
            os.system('chmod +x fid.com')

            # Run the fid.com file
            command = 'csh fid.com'
            p = subprocess.Popen(command, shell=True)
            p.wait()

            # Check to see if the output file exists
            if(os.path.exists('test.fid') == False):
                dlg = wx.MessageDialog(self, 'The converted FID file (test.fid) file cannot be found in the current directory. Conversion unsuccessful. Ensure that nmrPipe has been downloaded and added to the path.', 'Warning', wx.OK | wx.ICON_WARNING)
                self.Raise()
                self.SetFocus()
                dlg.ShowModal()
                if(self.file_parser == True):
                    os.chdir(self.cwd)
                return
            else:
                dlg = wx.MessageDialog(self, 'The FID file has been successfully converted to nmrPipe format (test.fid)', 'Success', wx.OK | wx.ICON_INFORMATION)
                self.Raise()
                self.SetFocus()
                dlg.ShowModal()
                if(self.file_parser == True):
                    os.chdir(self.cwd)
                return


    def perform_conversion_windows(self):
        # Convert the udic object read from the raw fid into nmrpipe format

        # Set Rance-Kay/echo-antiecho to false initially
        rance_kay = False

        if(self.nmrdata.spectrometer=='Bruker'):
            C = ng.convert.converter()



            if(len(self.N_complex_boxes)==1):
                u = {'ndim': 1, 0: {'sw': 0, 'complex': True, 'obs': 0, 'car': 0, 'size': 0, 'label': '', 'encoding': 'direct', 'time': True, 'freq': False}}
                u['ndim'] = 1
                u[0]['size'] = int(int(self.N_complex_boxes[0].GetValue().strip())/2)
                if(self.acqusition_combo_boxes[0].GetValue().strip() == 'Real'):
                    u[0]['complex'] = False
                else:
                    u[0]['complex'] = True
                u[0]['encoding'] = 'direct'
                u[0]['sw'] = float(self.sweep_width_boxes[0].GetValue().strip())
                u[0]['obs'] = float(self.nuclei_frequency_boxes[0].GetValue().strip())
                u[0]['car'] = float(self.carrier_frequency_boxes[0].GetValue().strip())*u[0]['obs']
                u[0]['label'] = self.nucleus_type_boxes[0].GetValue().strip()
    
            elif(len(self.N_complex_boxes)==2):
                u = {'ndim': 2, 0: {'sw': 0, 'complex': True, 'obs': 0, 'car': 0, 'size': 0, 'label': '', 'encoding': 'direct', 'time': True, 'freq': False}, 1: {'sw': 0, 'complex': True, 'obs': 0, 'car': 0, 'size': 0, 'label': '', 'encoding': 'direct', 'time': True, 'freq': False}}
                u[1]['size'] = int(int(self.N_complex_boxes[0].GetValue().strip())/2)
                if(self.acqusition_combo_boxes[0].GetValue().strip() == 'Real'):
                    u[1]['complex'] = False
                else:
                    u[1]['complex'] = True
                u[1]['encoding'] = 'direct'
                u[1]['sw'] = float(self.sweep_width_boxes[0].GetValue().strip())
                u[1]['obs'] = float(self.nuclei_frequency_boxes[0].GetValue().strip())
                u[1]['car'] = float(self.carrier_frequency_boxes[0].GetValue().strip())*u[1]['obs']
                u[1]['label'] = self.nucleus_type_boxes[0].GetValue().strip()
    
                u[0]['size'] = int(self.N_complex_boxes[1].GetValue().strip())
                if(self.acqusition_combo_boxes[1].GetValue().strip() == 'Real'):
                    u[0]['complex'] = False
                else:
                    u[0]['complex'] = True
                if(self.acqusition_combo_boxes[1].GetValue().strip() == 'Real'):
                    u[0]['encoding'] = 'real'
                elif(self.acqusition_combo_boxes[1].GetValue().strip() == 'Complex'):
                    u[0]['encoding'] = 'complex'
                elif(self.acqusition_combo_boxes[1].GetValue().strip() == 'States'):
                    u[0]['encoding'] = 'states'
                elif(self.acqusition_combo_boxes[1].GetValue().strip() == 'TPPI'):
                    u[0]['encoding'] = 'tppi'
                elif(self.acqusition_combo_boxes[1].GetValue().strip() == 'States-TPPI'):
                    u[0]['encoding'] = 'states-tppi'
                elif(self.acqusition_combo_boxes[1].GetValue().strip() == 'Echo-Antiecho' or self.acqusition_combo_boxes[1].GetValue().strip() == 'Echo-AntiEcho' or self.acqusition_combo_boxes[1].GetValue().strip() == 'Rance-Kay'):
                    u[0]['encoding'] = 'complex'
                    rance_kay = True
                u[0]['sw'] = float(self.sweep_width_boxes[1].GetValue().strip())
                u[0]['obs'] = float(self.nuclei_frequency_boxes[1].GetValue().strip())
                u[0]['car'] = float(self.carrier_frequency_boxes[1].GetValue().strip())*u[0]['obs']
                u[0]['label'] = self.nucleus_type_boxes[1].GetValue().strip()

            else:
                u = {'ndim': 3, 0: {'sw': 0, 'complex': True, 'obs': 0, 'car': 0, 'size': 0, 'label': '', 'encoding': 'direct', 'time': True, 'freq': False}, 1: {'sw': 0, 'complex': True, 'obs': 0, 'car': 0, 'size': 0, 'label': '', 'encoding': 'direct', 'time': True, 'freq': False},2: {'sw': 0, 'complex': True, 'obs': 0, 'car': 0, 'size': 0, 'label': '', 'encoding': 'direct', 'time': True, 'freq': False}}
                u[2]['size'] = int(int(self.N_complex_boxes[0].GetValue().strip())/2)
                if(self.acqusition_combo_boxes[0].GetValue().strip() == 'Real'):
                    u[2]['complex'] = False
                else:
                    u[2]['complex'] = True
                u[2]['encoding'] = 'direct'
                u[2]['sw'] = float(self.sweep_width_boxes[0].GetValue().strip())
                u[2]['obs'] = float(self.nuclei_frequency_boxes[0].GetValue().strip())
                u[2]['car'] = float(self.carrier_frequency_boxes[0].GetValue().strip())*u[2]['obs']
                u[2]['label'] = self.nucleus_type_boxes[0].GetValue().strip()
    
                u[1]['size'] = int(self.N_complex_boxes[1].GetValue().strip())
                if(self.acqusition_combo_boxes[1].GetValue().strip() == 'Real'):
                    u[1]['complex'] = False
                else:
                    u[1]['complex'] = True
                if(self.acqusition_combo_boxes[1].GetValue().strip() == 'Real'):
                    u[1]['encoding'] = 'real'
                elif(self.acqusition_combo_boxes[1].GetValue().strip() == 'Complex'):
                    u[1]['encoding'] = 'complex'
                elif(self.acqusition_combo_boxes[1].GetValue().strip() == 'States'):
                    u[1]['encoding'] = 'states'
                elif(self.acqusition_combo_boxes[1].GetValue().strip() == 'TPPI'):
                    u[1]['encoding'] = 'tppi'
                elif(self.acqusition_combo_boxes[1].GetValue().strip() == 'States-TPPI'):
                    u[1]['encoding'] = 'states-tppi'
                elif(self.acqusition_combo_boxes[1].GetValue().strip() == 'Echo-Antiecho' or self.acqusition_combo_boxes[1].GetValue().strip() == 'Echo-AntiEcho' or self.acqusition_combo_boxes[1].GetValue().strip() == 'Rance-Kay'):
                    u[1]['encoding'] = 'complex'
                    rance_kay = True
                u[1]['sw'] = float(self.sweep_width_boxes[1].GetValue().strip())
                u[1]['obs'] = float(self.nuclei_frequency_boxes[1].GetValue().strip())
                u[1]['car'] = float(self.carrier_frequency_boxes[1].GetValue().strip())*u[1]['obs']
                u[1]['label'] = self.nucleus_type_boxes[1].GetValue().strip()

                u[0]['size'] = int(self.N_complex_boxes[2].GetValue().strip())
                if(self.acqusition_combo_boxes[2].GetValue().strip() == 'Real'):
                    u[0]['complex'] = False
                else:
                    u[0]['complex'] = True
                if(self.acqusition_combo_boxes[2].GetValue().strip() == 'Real'):
                    u[0]['encoding'] = 'real'
                elif(self.acqusition_combo_boxes[2].GetValue().strip() == 'Complex'):
                    u[0]['encoding'] = 'complex'
                elif(self.acqusition_combo_boxes[2].GetValue().strip() == 'States'):
                    u[0]['encoding'] = 'states'
                elif(self.acqusition_combo_boxes[2].GetValue().strip() == 'TPPI'):
                    u[0]['encoding'] = 'tppi'
                elif(self.acqusition_combo_boxes[2].GetValue().strip() == 'States-TPPI'):
                    u[0]['encoding'] = 'states-tppi'
                elif(self.acqusition_combo_boxes[2].GetValue().strip() == 'Echo-Antiecho' or self.acqusition_combo_boxes[2].GetValue().strip() == 'Echo-AntiEcho' or self.acqusition_combo_boxes[2].GetValue().strip() == 'Rance-Kay'):
                    u[0]['encoding'] = 'complex'
                    rance_kay = True
                u[0]['sw'] = float(self.sweep_width_boxes[2].GetValue().strip())
                u[0]['obs'] = float(self.nuclei_frequency_boxes[2].GetValue().strip())
                u[0]['car'] = float(self.carrier_frequency_boxes[2].GetValue().strip())*u[0]['obs']
                u[0]['label'] = self.nucleus_type_boxes[2].GetValue().strip()


                # # 3D data conversion is currently not supported. Please use a linux/mac containing nmrpipe to complete conversion
                # dlg = wx.MessageDialog(self, '3D data conversion is currently not supported. Please use a linux/mac machine containing nmrPipe to complete conversion', 'Warning', wx.OK | wx.ICON_WARNING)
                # self.Raise()
                # self.SetFocus()
                # dlg.ShowModal()
                # if(self.file_parser == True):
                #     os.chdir(self.cwd)
                # dlg.Destroy()
                # return

            if(self.d_type == 0):
                dtypa = False
            else:
                dtypa = True
            
            if(self.byte_order==0):
                big = False
            else:
                big = True

            cplex = True


            dic,data = ng.bruker.read('./')  

            if(len(self.N_real_boxes)==2):
                # If have 2D data but nmrglue has read it in as a 1D, need to split it up
                if(len(data.shape)==1):
                    data = np.array(np.split(data,int(self.N_real_boxes[-1].GetValue())))

            data = ng.proc_base.fft(data)
            data_proc = ng.bruker.remove_digital_filter(dic,data,post_proc=True)
            data = ng.proc_base.ifft(data_proc)


                # Rance-Kay/Echo-Antiecho reshuffling taken from https://github.com/jjhelmus/nmrglue/issues/149 
            if(rance_kay == True):
                summ, diff = data[0::2], data[1::2]
                A = summ - diff
                B = -1j * (summ + diff)
                shuffled_data = np.empty(data.shape, data.dtype)
                shuffled_data[0::2] = A
                shuffled_data[1::2] = B


            C.from_bruker(dic, data,u)
            pdic, pdata = C.to_pipe()
            pdic["FDPIPEFLAG"] = 1.0
            

            if(u[0]['encoding']=='real' and len(self.N_complex_boxes)==2):
                pdic['FDF1TDSIZE'] = u[0]['size'] 
                pdic['FDF1FTSIZE'] = u[0]['size'] 
                pdic['FDF1APOD'] = u[0]['size'] 
                pdic['FDF1QUADFLAG'] = 1.0
                pdic['FDF1OBS'] = 1.0
                pdic['FDF1SW'] = 1.0
                pdic['FDF1ORIG'] = 1.0
                pdic['FD2DPHASE'] = 0


        else:
            C = ng.convert.converter()
            dic,data = ng.fileio.varian.read('./')

            if(len(self.N_complex_boxes)==1):
                u = {}
                u['ndim'] = 1
                u[0] = {}
                u[0]['size'] = int(int(self.N_complex_boxes[0].GetValue().strip())/2)
                if(self.acqusition_combo_boxes[0].GetValue().strip() == 'Real'):
                    u[0]['complex'] = False
                else:
                    u[0]['complex'] = True
                u[0]['encoding'] = 'direct'
                u[0]['sw'] = float(self.sweep_width_boxes[0].GetValue().strip())
                u[0]['obs'] = float(self.nuclei_frequency_boxes[0].GetValue().strip())
                u[0]['car'] = float(self.carrier_frequency_boxes[0].GetValue().strip())*u[0]['obs']
                u[0]['label'] = self.nucleus_type_boxes[0].GetValue().strip()
                u[0]['time'] = True
                u[0]['freq'] = False
    
            elif(len(self.N_complex_boxes)==2):
                u = {}
                u['ndim'] = 2
                u[0] = {}
                u[1] = {}
                u[1]['size'] = int(int(self.N_complex_boxes[0].GetValue().strip())/2)
                if(self.acqusition_combo_boxes[0].GetValue().strip() == 'Real'):
                    u[1]['complex'] = False
                else:
                    u[1]['complex'] = True
                u[1]['encoding'] = 'direct'
                u[1]['sw'] = float(self.sweep_width_boxes[0].GetValue().strip())
                u[1]['obs'] = float(self.nuclei_frequency_boxes[0].GetValue().strip())
                u[1]['car'] = float(self.carrier_frequency_boxes[0].GetValue().strip())*u[1]['obs']
                u[1]['label'] = self.nucleus_type_boxes[0].GetValue().strip()
                u[1]['time'] = True
                u[1]['freq'] = False
    
                u[0]['size'] = int(self.N_complex_boxes[1].GetValue().strip())
                if(self.acqusition_combo_boxes[1].GetValue().strip() == 'Real'):
                    u[0]['complex'] = False
                else:
                    u[0]['complex'] = True
                if(self.acqusition_combo_boxes[1].GetValue().strip() == 'Real'):
                    u[0]['encoding'] = 'real'
                elif(self.acqusition_combo_boxes[1].GetValue().strip() == 'Complex'):
                    u[0]['encoding'] = 'complex'
                elif(self.acqusition_combo_boxes[1].GetValue().strip() == 'States'):
                    u[0]['encoding'] = 'states'
                elif(self.acqusition_combo_boxes[1].GetValue().strip() == 'TPPI'):
                    u[0]['encoding'] = 'tppi'
                elif(self.acqusition_combo_boxes[1].GetValue().strip() == 'States-TPPI'):
                    u[0]['encoding'] = 'states-tppi'
                elif(self.acqusition_combo_boxes[1].GetValue().strip() == 'Echo-Antiecho' or self.acqusition_combo_boxes[1].GetValue().strip() == 'Echo-AntiEcho' or self.acqusition_combo_boxes[1].GetValue().strip() == 'Rance-Kay'):
                    u[0]['encoding'] = 'complex'
                    rance_kay = True
                u[0]['sw'] = float(self.sweep_width_boxes[1].GetValue().strip())
                u[0]['obs'] = float(self.nuclei_frequency_boxes[1].GetValue().strip())
                u[0]['car'] = float(self.carrier_frequency_boxes[1].GetValue().strip())*u[0]['obs']
                u[0]['label'] = self.nucleus_type_boxes[1].GetValue().strip()
                u[0]['time'] = False
                u[0]['freq'] = False
                if(u[0]['obs']==0.0):
                    u[0]['obs'] = 1.0
            
            else:
                # 3D data conversion is currently not supported. Please use a linux/mac containing nmrpipe to complete conversion
                dlg = wx.MessageDialog(self, '3D data conversion is currently not supported. Please use a linux/mac machine containing nmrPipe to complete conversion', 'Warning', wx.OK | wx.ICON_WARNING)
                self.Raise()
                self.SetFocus()
                dlg.ShowModal()
                if(self.file_parser == True):
                    os.chdir(self.cwd)
                dlg.Destroy()
                return
            

            # Rance-Kay/Echo-Antiecho reshuffling taken from https://github.com/jjhelmus/nmrglue/issues/149 
            if(rance_kay == True):
                summ, diff = data[0::2], data[1::2]
                A = summ - diff
                B = -1j * (summ + diff)
                shuffled_data = np.empty(data.shape, data.dtype)
                shuffled_data[0::2] = A
                shuffled_data[1::2] = B
                data = shuffled_data
            

            C.from_varian(dic,data,u)
            pdic, pdata = C.to_pipe()

            if(u[0]['encoding']=='real' and len(self.N_complex_boxes)==2):
                pdic['FDF1TDSIZE'] = u[0]['size'] 
                pdic['FDF1FTSIZE'] = u[0]['size'] 
                pdic['FDF1APOD'] = u[0]['size'] 
                pdic['FDF1QUADFLAG'] = 1.0
                pdic['FDF1OBS'] = 0.0
                pdic['FDF1SW'] = 0.0
                pdic['FDF1ORIG'] = 0.0
                pdic['FD2DPHASE'] = 0

 
        ng.pipe.write('test.fid', pdic,pdata, overwrite=True)

        # Check to see if the output file exists
        if(os.path.exists('test.fid') == False):
            dlg = wx.MessageDialog(self, 'The converted FID file (test.fid) file cannot be found in the current directory. Conversion unsuccessful.', 'Warning', wx.OK | wx.ICON_WARNING)
            self.Raise()
            self.SetFocus()
            dlg.ShowModal()
            if(self.file_parser == True):
                os.chdir(self.cwd)
            return
        else:
            dlg = wx.MessageDialog(self, 'The FID file has been successfully converted to nmrPipe format (test.fid)', 'Success', wx.OK | wx.ICON_INFORMATION)
            self.Raise()
            self.SetFocus()
            dlg.ShowModal()
            if(self.file_parser == True):
                os.chdir(self.cwd)
            return

        

    def create_other_options_box(self):
        self.bottom_left_box = wx.BoxSizer(wx.VERTICAL)
        self.other_options_box = wx.StaticBox(self, -1, label = 'Other Options')
        self.other_options_box_sizer_total = wx.StaticBoxSizer(self.other_options_box, wx.VERTICAL)
        self.other_options_box_sizer = wx.BoxSizer(wx.VERTICAL)
        self.other_options_box_sizer_total.AddSpacer(2)
        self.other_options_box_sizer_total.Add(self.other_options_box_sizer)
        self.other_options_box_sizer_total.AddSpacer(2)

        self.row_1 = wx.BoxSizer(wx.HORIZONTAL)

        # Bad point removal threshold
        self.bad_point_threshold_box = wx.BoxSizer(wx.HORIZONTAL)
        self.nmrdata.bad_point_threshold = 0.0
        self.bad_point_threshold_text = wx.StaticText(self, label='Bad Point Threshold:')
        self.bad_point_threshold_value = str(self.nmrdata.bad_point_threshold)
        self.bad_point_threshold_textbox = wx.TextCtrl(self, value=str(self.bad_point_threshold_value))
        self.bad_point_threshold_box.AddSpacer(5)
        self.bad_point_threshold_box.Add(self.bad_point_threshold_text)
        self.bad_point_threshold_box.AddSpacer(5)
        self.bad_point_threshold_box.Add(self.bad_point_threshold_textbox)
        self.row_1.Add(self.bad_point_threshold_box)

        # Remove acquisition padding
        self.remove_acquisition_padding_box = wx.BoxSizer(wx.HORIZONTAL)
        self.nmrdata.remove_acquisition_padding = True
        self.remove_acquisition_padding_checkbox = wx.CheckBox(self, label='Remove Acquisition Padding')
        self.remove_acquisition_padding_checkbox.Bind(wx.EVT_CHECKBOX, self.on_remove_acquisition_padding_checkbox)
        self.remove_acquisition_padding_checkbox.SetValue(True)
        self.remove_acquisition_padding_box.AddSpacer(5)
        self.remove_acquisition_padding_box.Add(self.remove_acquisition_padding_checkbox)
        self.row_1.AddSpacer(10)
        self.row_1.Add(self.remove_acquisition_padding_box)

        self.other_options_box_sizer.Add(self.row_1)

        self.bottom_left_box.Add(self.other_options_box_sizer_total)
        self.extra_boxes.AddSpacer(10)
        self.extra_boxes.Add(self.bottom_left_box)


    def on_remove_acquisition_padding_checkbox(self,event):
        if(self.remove_acquisition_padding_checkbox.GetValue() == True):
            self.nmrdata.remove_acquisition_padding = True
        else:
            self.nmrdata.remove_acquisition_padding = False
    
    

def main():
    app = wx.App()
    frame = MyApp()
    app.MainLoop()


if __name__ == '__main__':
    main()