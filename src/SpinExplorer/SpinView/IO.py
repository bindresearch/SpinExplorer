import wx # type: ignore
import os
import traceback 
import nmrglue as ng # type: ignore
import numpy as np
from natsort import natsorted

from SpinExplorer.SpinExplorer_CL_tools.processSpec import FindingParameters
from SpinExplorer.SpinExplorer_CL_tools.pulse_sequence_parsing import (
    PulseSequenceParser,
    PulseSequenceError,
    PulseProgramNotFoundError,
    NoConfigurationError,
)
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

        if getattr(self, "auto_failed", False) == True:
            # Why there is nothing to show has been said already
            return

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
            if file.endswith(".ft1"):
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
                self.file = self.auto_process()
            except Exception as problem:
                self.report_auto_failure(problem)

        if len(spectrum_file) == 1:
            self.file = spectrum_file[0]
        if len(spectrum_file) > 1:
            res = ChooseFile(spectrum_file, self)
            res.Raise()
            res.SetFocus()
            res.ShowModal()
            res.Destroy()

    def auto_process(self) -> str:
        """
        Convert and process the raw data in this directory, which is done when
        there is no processed spectrum to show. The recipe is the one registered
        for the pulse programme which was run.
        """
        print("No processed data found, so the raw data will be processed.")

        parameters = FindingParameters()

        self.auto_sequence = PulseSequenceParser().parse()
        print("Pulse programme: {}".format(self.auto_sequence))

        config = registry.get_default_config(self.auto_sequence)

        converted = Convert_nmrglue(parameters.params, parameters)

        # The sizes are what the conversion depends on, so they are worth saying
        # out loud: a dimension which was not read from the parameter files
        # leaves the data with fewer dimensions than the experiment has
        print("Direct dimension: {} points, indirect: {}".format(
            parameters.params.size_direct, parameters.params.size_indirect))
        print("Dimensions being converted: {}".format(converted.complex_sizes))

        written = parameter_write_cl(converted, config)
        written.write_out_dict(written.dictionary)

        if converted.params.remove_filter_before_processing == True:
            remove_filter = False
        else:
            remove_filter = True

        config.process_data(
            pseudo_flag=converted.params.pseudo_flag, filter_removal=remove_filter
        )

        if os.path.exists(config.ft_name) == False:
            raise RuntimeError(
                "the processing finished without writing " + str(config.ft_name)
            )

        return config.ft_name

    def report_auto_failure(self, problem):
        """
        Say why the raw data could not be processed. The whole of what went
        wrong is printed to the terminal, and the window says which step it was,
        as each one fails for its own reasons and needs something different
        doing about it.
        """
        print(traceback.format_exc())

        # The reason is given here, so reading the data is not tried again and
        # does not add a second message of its own
        self.auto_failed = True

        sequence = getattr(self, "auto_sequence", "")

        if isinstance(problem, NoConfigurationError) == True:
            message = (
                "The pulse programme which was run ({}) is not one of the "
                "experiments which can be processed automatically. Process the "
                "data with SpinProcess instead.".format(sequence)
            )
        elif isinstance(problem, PulseProgramNotFoundError) == True:
            message = (
                "No processed data was found in this directory, and the "
                "pulseprogram file which says which experiment was run is not "
                "here either, so the raw data cannot be processed "
                "automatically. Process the data with SpinProcess instead."
            )
        elif isinstance(problem, PulseSequenceError) == True:
            message = (
                "No processed data was found in this directory, and the name of "
                "the experiment could not be read from the pulseprogram file. "
                "Process the data with SpinProcess instead."
            )
        else:
            message = (
                "No processed data was found in this directory, and processing "
                "the raw data did not work. The reason was:\n\n{}\n\nThe whole "
                "of what went wrong has been printed to the terminal. Process "
                "the data with SpinProcess to look at it step by "
                "step.".format(problem)
            )

        dlg = wx.MessageDialog(
            self.tempframe, message, "Error", wx.OK | wx.ICON_INFORMATION
        )
        self.tempframe.Raise()
        self.tempframe.SetFocus()
        dlg.ShowModal()
        dlg.Destroy()
        self.app.Destroy()

    # Read in the NMRPipe data file
    def read_data(self):
        if getattr(self, "auto_failed", False) == True:
            # There is nothing to read, and why has already been said
            return

        self.found_file = False
        try:
            if self.file != ".":
                self.dic, self.data = ng.pipe.read(self.file)
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
                    "NMRPipe file not read properly. Ensure raw data is downloaded to the local device or please retry processing the data then try again.",
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

            # FDDIMORDER says which dimension of the spectrum each axis of the
            # data holds: the last axis is FDDIMORDER[0], the one before it is
            # FDDIMORDER[1] and so on. Taking the labels from it means they
            # always describe the data, however the spectrum was processed.
            def label_of_dimension(dimension):
                return self.dic["FDF{}LABEL".format(int(self.dic["FDDIMORDER"][dimension]))]

            if self.dim == 1 or self.dim == 2:
                # 1D and 2D data are labelled with the last axis of the data
                # first (the direct dimension of an ordinary spectrum)
                for dimension in range(self.dim):
                    self.axislabels.append(label_of_dimension(dimension))
            else:
                # 3D data is labelled in the order of the axes of the data
                for axis in range(3):
                    self.axislabels.append(label_of_dimension(2 - axis))

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
        file_selection = self.file_combobox.GetValue()
        self.parent.file = file_selection
        self.parent.session_file = file_selection
        self.Close()
        if self.session_choice == False:
            self.parent.read_data()
            self.parent.dim = self.parent.get_dimensions()
            self.parent.get_axislabels_nmrglue()
