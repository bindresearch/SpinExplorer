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

print("-------------------------------------------------------------")
print("                        SpinConverter                        ")
print("-------------------------------------------------------------")
print("               (version 1.3) 17th February 2026              ")
print(" (c) 2025 James Eaton, Andrew Baldwin (University of Oxford) ")
print("                  2025-2026, Bind Research                   ")
print("                        MIT License                          ")
print("-------------------------------------------------------------")
print("            Converting NMR data to nmrPipe format            ")
print("-------------------------------------------------------------")
print(" Video tutorials at:")
print(" https://www.youtube.com/@BindResearch")
print("-------------------------------------------------------------")
print("")


# Import relevant external modules
import sys
import wx
import subprocess
import darkdetect
import warnings
import pathlib
import wx.adv
from appdirs import user_data_dir
import os

# Importing internal classes
from SpinExplorer.SpinConverter.FindingParameters.parameters import FindingParameters
from SpinExplorer.SpinConverter.FormattingGUI.bruker_formatting import (
    FormatParametersBruker,
)
from SpinExplorer.SpinConverter.FormattingGUI.varian_formatting import (
    FormatParametersVarian,
)
from SpinExplorer.SpinConverter.FormattingGUI.shared_formatting import (
    SharedFormatting,
)
from SpinExplorer.SpinConverter.StoringParameters.save_parameters import Save_json
from SpinExplorer.SpinConverter.StoringParameters.read_parameters import Read_json
from SpinExplorer.SpinConverter.Conversion.convert_pipe import Convert_pipe
from SpinExplorer.SpinConverter.Conversion.convert_nmrglue import Convert_nmrglue
from SpinExplorer.SpinConverter.Conversion.add_fid import Add_fid

from SpinExplorer.SpinExpLogo import SpinExpLogo


warnings.simplefilter("ignore", UserWarning)

    # except:
    #     platform = "windows"


# James Eaton, 10/06/2025, University of Oxford
# James Eaton, 25/09/2025, Bind Research
# This program is designed to allow the user to convert NMR data from Bruker/Varian into NMRPipe format so it can be viewed using
# SpinView.py. It is designed to be used with the SpinProcess.py program used to process the converted nmrPipe FID to produce
# an NMR spectrum. These spectra can then be viewed using SpinView.py, a GUI for viewing NMR data.


# task bar dock icon adapted from https://wiki.wxpython.org/Custom%20Mac%20OsX%20Dock%20Bar%20Icon
class TaskBarIcon(wx.adv.TaskBarIcon):
    TBMENU_CLOSE   = wx.NewId()
    TBMENU_CHANGE  = wx.NewId()
    TBMENU_REMOVE  = wx.NewId()
   
    def __init__(self, frame):
        wx.adv.TaskBarIcon.__init__(self, iconType=wx.adv.TBI_DOCK)
        self.frame = frame

        # Set the image
        icon = self.MakeIcon(SpinExpLogo.GetImage())
        self.SetIcon(icon, "SpinConverter")
        self.imgidx = 1      
        # bind some events
        self.Bind(wx.EVT_MENU, self.OnTaskBarClose, id=self.TBMENU_CLOSE)

    def CreatePopupMenu(self):
        """
        This method is called by the base class when it needs to popup
        the menu for the default EVT_RIGHT_DOWN event.
        """
        menu = wx.Menu()
        menu.Append(self.TBMENU_CLOSE,   "Close SpinConverter")
        return menu

    def MakeIcon(self, img):
        icon = wx.Icon()
        icon.CopyFromBitmap(img.ConvertToBitmap())
        return icon

    def OnTaskBarClose(self, evt):
        self.frame.Destroy()
        sys.exit()




class SpinConverter(wx.Frame):
    def __init__(self, explorer=False):
        """
        This class creates the GUI showing the found parameters with scope
        for changing the parameters.
        """
        # Check the platform and if nmrpipe is installed
        self.check_platform()

        # Get the title for the panels
        self.title = self.GetTitle()
        
        # Setup the dock/task bar with the logo
        if(explorer==False):
            self.tbicon = TaskBarIcon(self)

        # Get the monitor size and set the window size to 85% of the monitor size
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 0.6 * self.monitorWidth
        self.height = 0.6 * self.monitorHeight
        self.app_frame = wx.Frame.__init__(
            self,
            None,
            wx.ID_ANY,
            self.title,
            wx.DefaultPosition,
            size=(int(self.width), int(self.height)),
        )

        self.file_parser = False

        # Initialise the NMR data and parameter class
        self.nmrdata = FindingParameters(self)

        # Creating a canvas and formatting the app on it
        self.create_canvas()
        self.format_app()

        # Reading previously saved parameters if present
        read = Read_json(self.nmrdata.params, self.nmrdata, self)

        self.Bind(wx.EVT_CLOSE, self.OnClose)

        
        self.SetIcon(SpinExpLogo.GetIcon())
        self.Show()
        self.Centre()


    def check_platform(self):

        appname = 'SpinExplorer'
        appauthor = "James Eaton"
        data_dir = user_data_dir(appname, appauthor)
        os.makedirs(data_dir, exist_ok=True)
        logfile_location = os.path.join(data_dir, 'logfile.txt')


        # Check to see if using mac, linux or windows
        if sys.platform == "darwin":
            platform = "mac"
        elif sys.platform == "linux":
            platform = "linux"
        else:
            platform = "windows"


        # See if the nmrPipe command works, if not set the platform to windows
        if platform == "mac" or platform == "linux":

                p = subprocess.Popen(["csh", "-c", 'nmrPipe'], stdout=subprocess.DEVNULL,
                     stderr=subprocess.PIPE)
                out, err = p.communicate()
                
                if "NMRPipe System Version" in str(out) or "NMRPipe System Version" in str(err):
                    platform = platform
                else:
                    platform = "windows"


        self.platform = platform

    def GetTitle(self):
        """
        Finding an appropriate title for the panel. 
        The title for the panel is:
        SpinConverter + the current working directory (last 3 elements)
        + the title + pulseprogram (for Bruker data)
        """
        title = 'SpinConverter: '
        p = pathlib.Path.cwd()
        dirs = p.parts[-3:]
        last_directories_path = pathlib.Path(*dirs)
        title = title + "/" + str(last_directories_path)
        
        # If pdata/1/title exists, add this title too
        try:
            with open('pdata/1/title') as file:
                line = file.readlines()[0]
                title_extra = ''
                line = line.split('\n')[0]
                title_extra+= line 

                
                title = title + '(' + title_extra + ')'
        except:
            pass

        return title


    def OnClose(self, event):
        """
        Ensuring the application is closed after pressing close
        """
        self.Destroy()
        # sys.exit()

    def create_canvas(self) -> None:
        """
        Creating a canvas for the application
        """
        # Create the main sizer
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)


    def create_sizers(self) -> None:
        self.parameters_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.menu_bar = wx.BoxSizer(wx.VERTICAL)
        self.extra_sizers = wx.BoxSizer(wx.VERTICAL)
        self.parameters_sizer.Add(self.menu_bar)
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

    def format_app(self) -> None:
        """
        Adding the TextControl sizers and buttons to the app and
        populating them with parameters.
        """

        self.create_sizers()

        if self.nmrdata.spectrometer == "Bruker":
            self.format = FormatParametersBruker(
                self, self.nmrdata.params, self.nmrdata
            )
            self.shared_format = SharedFormatting(
                self, self.nmrdata.params, self.nmrdata
            )
            self.format.input_sizes_bruker()
            self.format.input_acquisition_modes_bruker()
            self.format.input_sweep_widths_bruker()
            self.format.get_nuclei_frequency_bruker()
            self.format.get_nuclei_labels_bruker()
            self.format.get_carrier_frequencies_bruker()
            if len(self.format.N_complex_boxes) > 1:
                self.shared_format.acquisition_2D_mode_combo_box()
            self.shared_format.create_temperature_box()
            self.format.create_bruker_digital_filter_box()
            self.shared_format.create_conversion_box()
            self.format.create_other_options_box()
            self.shared_format.create_intensity_scaling_box()
            if self.nmrdata.params.size_indirect != []:
                self.shared_format.find_nus_file()
                self.shared_format.input_NUS_list_box()
            else:
                self.shared_format.include_NUS = False
        elif self.nmrdata.spectrometer == "Varian":
            self.format = FormatParametersVarian(self, self.nmrdata.params, self.nmrdata)
            self.shared_format = SharedFormatting(
                self, self.nmrdata.params, self.nmrdata
            )
            self.format.input_sizes_varian()
            self.format.input_acquisition_modes_varian()
            self.format.input_sweep_widths_varian()
            self.format.get_nuclei_frequency_varian()
            self.format.get_nuclei_labels_varian()
            self.format.get_carrier_frequencies_varian()
            if len(self.format.N_complex_boxes) > 1:
                self.shared_format.acquisition_2D_mode_combo_box()
            self.shared_format.create_temperature_box()
            self.shared_format.create_conversion_box()
            self.shared_format.create_intensity_scaling_box()
            if self.nmrdata.params.phase != False or self.nmrdata.params.phase2 != False:
                self.shared_format.find_nus_file()
                self.shared_format.input_NUS_list_box()
            else:
                self.shared_format.include_NUS = False

        self.main_sizer.Add(self.parameters_sizer, 0, wx.CENTER)
        self.main_sizer.Add(self.extra_sizers, 0, wx.CENTER)

        self.SetSizerAndFit(self.main_sizer)

        # Get the width and height of the main_sizer
        self.width, self.height = self.main_sizer.GetSize()
        self.SetSize((int(self.width * 1.25), int(self.height * 1.25)))
        self.Centre()

    def on_save_parameters(self, event) -> None:
        """
        Saving the current SpinConverter parameters to parameters.json
        """
        save = Save_json(self.nmrdata.params, self.nmrdata, self)

    def on_convert_pipe(self, event) -> None:
        """
        Checking to see that nmrPipe is installed and then performing
        nmrPipe conversion.
        """
        
        if self.platform == "windows":
            # Outputting a message saying that nmrPipe conversion is not possible on windows
            dlg = wx.MessageDialog(
                self,
                "It seems like nmrPipe is not installed. Please use the nmrglue convert button instead.",
                "Warning",
                wx.OK | wx.ICON_WARNING,
            )
            self.Raise()
            self.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            return
        else:
            # Saving conversion parameters
            self.on_save_parameters(wx.EVT_BUTTON)
            # Performing nmrpipe conversion
            pipe_conversion = Convert_pipe(self, self.nmrdata.params, self.nmrdata)

    def on_convert_glue(self, event) -> None:
        # Saving conversion parameters
        self.on_save_parameters(wx.EVT_BUTTON)
        # Performing nmrglue conversion
        glue_conversion = Convert_nmrglue(self, self.nmrdata.params, self.nmrdata)

    
    def on_add_fids(self, event) -> None:
        """
        Creating a popout where users can select different fids to add up (.fid files)
        and can also remove selected FIDs too. When the Add button is pressed, the user
        chooses a new directory to call the folder (with added fid saved as test.fid).
        If the addition fails, an error message pops out and no addition takes place
        """
        add_fid = Add_fid(self, self.nmrdata)



def main():
    app = wx.App()
    frame = SpinConverter()
    app.MainLoop()


if __name__ == "__main__":
    main()
