import wx # type: ignore
import os 
import nmrglue as ng # type: ignore
import numpy as np
from natsort import natsorted

from SpinExplorer.SpinExplorer_CL_tools.processSpec import FindingParameters
from SpinExplorer.SpinExplorer_CL_tools.pulse_sequence_parsing import PulseSequenceParser
from SpinExplorer.SpinExplorer_CL_tools.make_parameter_file_cl import parameter_write_cl
from SpinExplorer.SpinExplorer_CL_tools.convert_nmrglue_cl import Convert_nmrglue
from SpinExplorer.SpinExplorer_CL_tools.config_register import registry

# This class reads in the NMRPipe data
class GetData:
    def __init__(self, app, file=""):

        # Create a hidden frame to be used as a parent for popout messages
        self.tempframe = wx.Frame(None, title="Temporary Parent", size=(1, 1))
        self.tempframe.Hide()  # Hide the frame since we don't need it to be visible

        self.app = app

        self.file = file
        self.path = os.getcwd()
        if self.file == "":
            self.get_filename()
        self.read_data()
        self.dim = self.get_dimensions()
        if self.file != ".":
            # NMRPipe data
            self.get_axislabels_nmrglue()
        else:
            # Bruker Topspin processed data
            self.generic_labels_bruker()

    # Get the filename of the NMRPipe data file
    def get_filename(self):
        self.found_file = False
        current_directory = os.getcwd()
        files = os.listdir(current_directory)
        spectrum_file = []
        self.brukerdata = False
        for file in files:
            if file.endswith(".ft"):
                spectrum_file.append(file)
            if file.endswith(".ft2"):
                spectrum_file.append(file)
            if file.endswith(".ft3"):
                spectrum_file.append(file)
            if file.endswith(".pipe"):
                spectrum_file.append(file)
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
                # Topspin processed Bruker data is present
                spectrum_file.append(".")
                break

        if len(spectrum_file) == 0:
            try:
                dlg = wx.MessageDialog(
                    self.tempframe,
                    "No NMRPipe or Bruker data files in current directory. We will attempt to auto-analyse.",
                    "Information",
                    wx.OK | wx.ICON_INFORMATION,
                )
                input_dat = FindingParameters()
                pp_parser = PulseSequenceParser()
                sequence = pp_parser.parse()

                config = registry.get_default_config(sequence)

                nmr_glue_conv = Convert_nmrglue(input_dat.params, input_dat)

                params = parameter_write_cl(nmr_glue_conv, config)
                params.write_out_dict(params.dictionary)
        
                config.process_data(pseudo_flag=nmr_glue_conv.params.pseudo_flag)
                self.file = config.ft_name
            except:
                dlg = wx.MessageDialog(
                    self.tempframe,
                    "No NMRPipe or Bruker data files in current directory and automatic processing failed.",
                    "Error",
                    wx.OK | wx.ICON_INFORMATION,
                )
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                dlg.ShowModal()
                dlg.Destroy()
                self.app.Destroy()
        if len(spectrum_file) == 1:
            self.file = spectrum_file[0]
        if len(spectrum_file) > 1:
            res = ChooseFile(spectrum_file, self)
            res.Raise()
            res.SetFocus()
            res.ShowModal()
            res.Destroy()

    # Read in the NMRPipe data file
    def read_data(self):
        self.found_file = False
        try:
            if self.file != ".":
                self.dic, self.data = ng.pipe.read(self.file)
                print(self.dic)
                if('nmrglue' in self.dic['FDCOMMENT']):
                    self.nmrglue_flag = True
                else:
                    self.nmrglue_flag = False
                if('pseudo' in self.dic['FDCOMMENT']):
                    self.pseudo_flag = True
                else:
                    self.pseudo_flag = False
            else:
                self.dic, self.data = ng.bruker.read_pdata(self.file)
            if len(self.data) == 0:
                # Give a popout saying the NMRPipe file has not been read properly. Retry processing
                dlg = wx.MessageDialog(
                    self.tempframe,
                    "Data file was read but data array is empty. Ensure raw data is downloaded to the local device.",
                    "Error",
                    wx.OK | wx.ICON_INFORMATION,
                )
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                dlg.ShowModal()
                dlg.Destroy()
                self.found_file = True
                self.app.Destroy()


        except:
            if self.found_file == False:
                # Give a popout saying the NMRPipe file has not been read properly. Retry processing
                dlg = wx.MessageDialog(
                    self.tempframe,
                    "NMRPipe file not read properly. Please retry processing the file then try again.",
                    "Error",
                    wx.OK | wx.ICON_INFORMATION,
                )
                self.tempframe.Raise()
                self.tempframe.SetFocus()
                dlg.ShowModal()
                dlg.Destroy()
                self.app.Destroy()

    # Work out NMR spectrum dimensions in order to get the plotting correct (need contour plot for 2D/3D but not for 1D)
    def get_dimensions(self):
        if (
            type(self.data[0]) == np.float32
            or type(self.data[0]) == np.float64
            or type(self.data[0]) == np.complex64
        ):
            return 1
        if len(self.data.shape) == 2:
            pseudo = False
            for val in self.data.shape:
                if val == 1:
                    pseudo = True
            if pseudo == True:
                self.data = self.data[0]
                return 1
            else:
                return 2
        if len(self.data.shape) == 3:
            pseudo = False
            for val in self.data.shape:
                if val == 1:
                    pseudo = True
            if pseudo == True:
                self.data_new = []
                for i, val2 in enumerate(self.data):
                    if self.data.shape[i] != 1:
                        self.data_new.append(val2)
                self.data = self.data_new
                return 3
            else:
                return 3

    def read_labels_file(self):
        file = open("labels.txt", "r")
        label = file.readlines()
        for i, line in enumerate(label):
            if i == 0:
                line = line.split("\n")[0].split(",")
                self.axislabels = line
        file.close()


    def get_axislabels_nmrglue(self):
        """
        Reading the nmrglue dictionary (self.dic) to obtain the correct axis
        labels associated with the data.
        """

        try:
            # If the user has already opened and customised the labels they will be in the labels.txt file
            self.read_labels_file()
        except:

            self.axislabels = []


            if self.dim == 1:
                # If 1D take FDF1LABEL
                self.axislabels.append(self.dic["FDF1LABEL"])
            elif self.dim == 2:
                # If 2D take FDF2LABEL as direct and FDF1LABEL as indirect
                if(self.pseudo_flag == False):
                    self.axislabels.append(self.dic["FDF1LABEL"])
                    self.axislabels.append(self.dic["FDF2LABEL"])
                else:
                    self.axislabels.append(self.dic["FDF2LABEL"])
                    self.axislabels.append(self.dic["FDF1LABEL"])
            else:
                # If 3D take FDF3LABEL as direct, FDF1LABEL as indirect1 and FDF2LABEL as indirect3
                if(self.pseudo_flag==True and self.nmrglue_flag==True):
                    self.axislabels.append(self.dic["FDF3LABEL"])
                    self.axislabels.append(self.dic["FDF2LABEL"])
                    self.axislabels.append(self.dic["FDF1LABEL"])
                else:
                    self.axislabels.append(self.dic["FDF1LABEL"])
                    self.axislabels.append(self.dic["FDF2LABEL"])
                    self.axislabels.append(self.dic["FDF3LABEL"])

    def generic_labels_bruker(self):
        """
        Input generic dim1, dim2, dim3 axis labels for Topspin
        processed data. This is temporary and should be updated
        to include correct labels from the Bruker dictionary.
        """
        self.labels = []

        if self.dim == 1:
            self.labels = ["dim1"]
        if self.dim == 2:
            self.labels = ["dim1", "dim2"]
        if self.dim == 3:
            self.labels = ["dim1", "dim2", "dim3"]
        self.axislabels = self.labels


class ChooseFile(wx.Dialog):
    def __init__(self, spectrum_file, parent, session_choice=False):
        if session_choice == False:
            name = "Select NMRPipe Data File"
        else:
            name = "Select Session File"
        dialog = wx.Dialog.__init__(
            self,
            None,
            wx.ID_ANY,
            name,
            wx.DefaultPosition,
            size=(300, 200),
            style=wx.DEFAULT_DIALOG_STYLE,
        )
        self.spectrum_file = spectrum_file
        self.parent = parent
        self.session_choice = session_choice
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.main_sizer.AddSpacer(10)
        if self.session_choice == False:
            self.message = wx.StaticText(
                self,
                label="Multiple NMRPipe data files in current directory. Please select an NMRPipe file to show.\n",
            )
        else:
            self.message = wx.StaticText(
                self,
                label="Multiple session files in current directory. Please select a session file to load.\n",
            )
        self.main_sizer.Add(self.message, 0, wx.ALL, 5)
        self.file_combobox = wx.ComboBox(
            self, choices=natsorted(spectrum_file), style=wx.CB_READONLY
        )
        self.main_sizer.Add(
            self.file_combobox, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 5
        )
        self.ok_button = wx.Button(self, label="OK")
        self.ok_button.Bind(wx.EVT_BUTTON, self.OnOK)
        self.main_sizer.Add(self.ok_button, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.ALL, 5)
        self.SetSizer(self.main_sizer)
        self.Centre()
        self.Show()

    def OnOK(self, event):
        file_selection = self.file_combobox.GetSelection()
        self.parent.file = self.spectrum_file[file_selection]
        self.parent.session_file = self.spectrum_file[file_selection]
        self.Close()
        if self.session_choice == False:
            self.parent.read_data()
            self.parent.dim = self.parent.get_dimensions()
            self.parent.get_axislabels_nmrglue()
