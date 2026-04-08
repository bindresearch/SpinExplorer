#!/usr/bin/env python3

"""
MIT License

Copyright (c) 2025 James Eaton, Andrew Baldwin (University of Oxford)
              2025-2026, Bind Research


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
print("                          SpinView                           ")
print("-------------------------------------------------------------")
print("               (version 1.4) 13th March 2026                 ")
print(" (c) 2025 James Eaton, Andrew Baldwin (University of Oxford) ")
print("                  2025-2026, Bind Research                   ")
print("                        MIT License                          ")
print("-------------------------------------------------------------")
print("              Viewing and analysing NMR spectra              ")
print("-------------------------------------------------------------")
print(" Video tutorials at:")
print(" https://www.youtube.com/@BindResearch")
print("-------------------------------------------------------------")
print("")


import sys
import os
import wx # type: ignore
import wx.adv # type: ignore

# Import relevant modules
import numpy as np
from appdirs import user_data_dir # type: ignore

appname = 'SpinExplorer'
appauthor = "James Eaton"
data_dir = user_data_dir(appname, appauthor)
os.makedirs(data_dir, exist_ok=True)

import matplotlib

mpl_cache = os.path.join(data_dir, "mpl-cache")
os.makedirs(mpl_cache, exist_ok=True)
matplotlib.get_cachedir = lambda: mpl_cache
import pathlib

matplotlib.use("wxAgg")

from matplotlib.backend_bases import MouseEvent as MPLMouseEvent
import sys
from SpinExplorer.SpinExpLogo import SpinExpLogo

from SpinExplorer.SpinView.Viewers.loading import GetData, ReadSession, ChooseFile
from SpinExplorer.SpinView.Viewers.oned_view import OneDViewer
from SpinExplorer.SpinView.Viewers.twod_view import TwoDViewer
from SpinExplorer.SpinView.Viewers.threed_view import ThreeDViewer

from SpinExplorer.SpinView.config import colours, twoD_colours
from SpinExplorer.SpinView.config import reference_range_values, multiply_range_values, vertical_range_values


matplotlib.rcParams["font.sans-serif"] = "Arial"
matplotlib.rcParams["font.family"] = "sans-serif"


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
        self.SetIcon(icon, "SpinView")
        self.imgidx = 1      
        # bind some events
        self.Bind(wx.EVT_MENU, self.OnTaskBarClose, id=self.TBMENU_CLOSE)

    def CreatePopupMenu(self):
        """
        This method is called by the base class when it needs to popup
        the menu for the default EVT_RIGHT_DOWN event.
        """
        menu = wx.Menu()
        menu.Append(self.TBMENU_CLOSE,   "Close SpinView")
        return menu

    def MakeIcon(self, img):
        icon = wx.Icon()
        icon.CopyFromBitmap(img.ConvertToBitmap())
        return icon

    def OnTaskBarClose(self, evt):
        self.frame.Destroy()
        sys.exit()




# This class creates the GUI main frame
class SpinView(wx.Frame):
    def __init__(self, explorer=False, session_file='', fid_data=[]):
        # Get the monitor size and set the window size to 85% of the monitor size
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]

        self.width = int(1.0 * sizes[0][0])
        self.height = int(0.875 * sizes[0][1])
        self.reprocess = False

        # Get the title for the panels
        self.title = self.GetTitle()

        # Setup the dock/task bar with the logo
        if(explorer==False):
            try:
                self.tbicon = TaskBarIcon(self)
            except:
                pass

        self.app_frame = wx.Frame.__init__(
            self,
            None,
            wx.ID_ANY,
            self.title,
            wx.DefaultPosition,
            size=(int(self.width), int(self.height)),
        )

        self.main_sizer = wx.BoxSizer(wx.VERTICAL)

        if(len(fid_data)==0):
            self.build_app(session_file)
        else:
            self.build_fid_app(fid_data)


    def build_fid_app(self, fid_data):

        self.main_sizer.AddSpacer(10)
        self.display_index_current = wx.Display.GetFromWindow(self)

        # Variables needed to set the correct path so code can be used with unidecFile parser
        self.path = ""
        self.cwd = ""
        self.file_parser = False

        nmrdata = fid_data[0]
        if nmrdata.dim == 1:
            self.viewer = OneDViewer(parent=self, nmrdata=nmrdata, fid_viewer=True)
            self.main_sizer.Add(self.viewer, 1, wx.EXPAND)
        elif nmrdata.dim == 2:
            self.viewer = TwoDViewer(parent=self, nmrdata=nmrdata, fid_viewer=True)
            self.main_sizer.Add(self.viewer, 1, wx.EXPAND)
        elif nmrdata.dim == 3:
            self.viewer = ThreeDViewer(parent=self, nmrdata=nmrdata, fid_viewer=True)
            self.main_sizer.Add(self.viewer, 1, wx.EXPAND)


        self.SetSizer(self.main_sizer)

        # Make negative contour lines solid
        matplotlib.rc("contour", negative_linestyle="solid")

        self.Show()
        self.Centre()

        # Bind method to check/resize the window when the frame is moved
        self.Bind(wx.EVT_MOVE, self.OnMoveFrame)

        # Bind method to resize the window when the frame is resized
        self.Bind(wx.EVT_SIZE, self.OnSizeFrame)

        self.Bind(wx.EVT_CLOSE, self.OnClose)

    def build_app(self, session_file=''):

        self.main_sizer.AddSpacer(10)
        self.display_index_current = wx.Display.GetFromWindow(self)

        # Variables needed to set the correct path so code can be used with unidecFile parser
        self.path = ""
        self.cwd = ""
        self.file_parser = False

        # Find if there are any sessions saved in the current directory
        if(session_file==''):
            self.find_sessions()
        else:
            self.session_file=session_file
        if self.session_file != "":
            ReadSession(self, self.session_file)
        else:
            self.nmrdata = GetData(self)
            if self.nmrdata.dim == 1:
                self.viewer = OneDViewer(parent=self, nmrdata=self.nmrdata)
                self.main_sizer.Add(self.viewer, 1, wx.EXPAND)
            elif self.nmrdata.dim == 2:
                self.viewer = TwoDViewer(parent=self, nmrdata=self.nmrdata)
                self.main_sizer.Add(self.viewer, 1, wx.EXPAND)
            elif self.nmrdata.dim == 3:
                self.viewer = ThreeDViewer(parent=self, nmrdata=self.nmrdata)
                self.main_sizer.Add(self.viewer, 1, wx.EXPAND)

        if(self.session_file!=''):
            self.title += ' Session: ' + pathlib.Path(self.session_file).name
            self.SetTitle(self.title)

        self.SetSizer(self.main_sizer)

        # Make negative contour lines solid
        matplotlib.rc("contour", negative_linestyle="solid")

        self.Show()
        self.Centre()

        # Bind method to check/resize the window when the frame is moved
        self.Bind(wx.EVT_MOVE, self.OnMoveFrame)

        # Bind method to resize the window when the frame is resized
        self.Bind(wx.EVT_SIZE, self.OnSizeFrame)

        self.Bind(wx.EVT_CLOSE, self.OnClose)

    def GetTitle(self):
        """
        Finding an appropriate title for the panel.
        The title for the panel is:
        SpinView + the current working directory (last 3 elements)
        + the title + pulseprogram (for Bruker data)
        """
        title = "SpinView: "
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
        # Save the session file if the user wants to save it
        if self.reprocess == True:
            return
        else:
            try:
                if self.nmrdata.dim == 1 or self.nmrdata.dim == 2:
                    dlg = wx.MessageDialog(
                        None,
                        "Do you want to save the session file?",
                        "Save Session File",
                        wx.YES_NO | wx.ICON_INFORMATION,
                    )
                    self.Raise()
                    self.SetFocus()
                    result = dlg.ShowModal()
                    if result == wx.ID_YES:
                        dlg.Destroy()
                        try:
                            if self.nmrdata.dim == 1:
                                self.viewer.OnSaveSessionButton(wx.EVT_BUTTON)
                            else:
                                self.viewer.OnSaveSessionButton2D(wx.EVT_BUTTON)
                        except:
                            dlg2 = wx.MessageDialog(
                                None,
                                "Session file not saved properly.",
                                "Save Session File",
                                wx.OK | wx.ICON_INFORMATION,
                            )
                            self.Raise()
                            self.SetFocus()
                            result2 = dlg2.ShowModal()
                            dlg2.Destroy()

                    dlg.Destroy()
                self.Destroy()
                # sys.exit()
            except:
                self.Destroy()
                # sys.exit()

    def OnMoveFrame(self, event):
        # Get the new default display if the frame is moved
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        display_index = wx.Display.GetFromWindow(self)
        if display_index != self.display_index_current:
            self.width = int(1.0 * sizes[display_index][0])
            self.height = int(0.875 * sizes[display_index][1])
            self.SetSize((self.width, self.height))
            self.viewer.canvas.SetSize(
                (
                    self.width * 0.0104,
                    (self.height - self.viewer.bottom_sizer.GetMinSize()[1] - 100)
                    * 0.0104,
                )
            )
            self.viewer.fig.set_size_inches(
                self.width * 0.0104,
                (self.height - self.viewer.bottom_sizer.GetMinSize()[1] - 100) * 0.0104,
            )
            self.viewer.UpdateFrame()
            self.display_index_current = display_index

        event.Skip()

    def OnSizeFrame(self, event):
        # Get the new frame size
        self.width, self.height = self.GetSize()
        self.SetSize((self.width, self.height))
        self.viewer.canvas.SetSize(
            (
                self.width * 0.0104,
                (self.height - self.viewer.bottom_sizer.GetMinSize()[1] - 100) * 0.0104,
            )
        )
        self.viewer.fig.set_size_inches(
            self.width * 0.0104,
            (self.height - self.viewer.bottom_sizer.GetMinSize()[1] - 100) * 0.0104,
        )
        self.viewer.UpdateFrame()
        event.Skip()

    def find_sessions(self, ask_user=True):
        self.sessions = []
        if self.path != "":
            os.chdir(self.path)
            files = os.listdir(self.path)
        else:
            files = os.listdir()
        for file in files:
            if file.endswith(".session"):
                self.sessions.append(file)
        if self.path != "":
            os.chdir(self.cwd)

        # If there are no found sessions then session flag needs to be set to False
        if len(self.sessions) == 0:
            self.session_file = ""
        elif len(self.sessions) == 1:
            # Ask the user if they want to load the session file
            if(ask_user==True):
                dlg = wx.MessageDialog(
                    None,
                    "Session file found ({}). Do you want to load the session file?".format(
                        self.sessions[0]
                    ),
                    "Session File Found",
                    wx.YES_NO | wx.ICON_INFORMATION,
                )
                self.Raise()
                self.SetFocus()
                result = dlg.ShowModal()
                if result == wx.ID_YES:
                    self.session_file = self.sessions[0]
                else:
                    self.session_file = ""
                dlg.Destroy()
            else:
                self.session_file = self.sessions[0]
        else:
            # Asking the user if they wish to load a session file
            if(ask_user==True):
                dlg = wx.MessageDialog(
                    None,
                    "Multiple session files found. Do you want to load a session file?",
                    "Session Files Found",
                    wx.YES_NO | wx.ICON_INFORMATION,
                )
                self.Raise()
                self.SetFocus()
                result = dlg.ShowModal()
                if result == wx.ID_YES:
                    dlg.Destroy()
                    # Asking the user to select the session file they want to load
                    res = ChooseFile(self.sessions, self, session_choice=True)
                    res.ShowModal()
                    res.Destroy()
                else:
                    self.session_file = ""
                    dlg.Destroy()
            else:
                # Asking the user to select the session file they want to load
                res = ChooseFile(self.sessions, self, session_choice=True)
                res.ShowModal()
                res.Destroy()

    # Initialising global app variables variables variables
    def set_variables(self):
        # Colours for 1D lines
        self.colours = colours
        self.colour_value = self.colours[0]

        # Initial 1D slice colour for 2D/3D spectra is set to navy
        self.colour_slice = "navy"

        # List of cmap colours for when overlaying multiple spectra
        self.cmap = "#e41a1c"
        self.cmap_neg = "#377eb8"
        self.twoD_colours = twoD_colours
        self.twoD_label_colours = self.twoD_colours

        self.twoD_slices_horizontal = []
        self.twoD_slices_vertical = []

        # Range of the sliders to for moving spectra left/right/up/down
        self.reference_range_values = reference_range_values
        self.reference_range = float(self.reference_range_values[0])
        self.reference_rangeX = float(self.reference_range_values[0])
        self.reference_rangeY = float(self.reference_range_values[0])

        # Range of the sliders to for moving spectra up/down in 1D spectra
        self.vertical_range_values = vertical_range_values

        # Range of the sliders to for multiplying 1D spectra
        self.multiply_range_values = multiply_range_values

        # Initial x,y movements for referencing are set to zero
        self.x_movement = 0
        self.y_movement = 0

        # Multiplot mode is initially set to off
        self.multiplot_mode = False

        # Dictionary to store the values of the sliders for each spectrum in multiplot mode
        self.values_dictionary = {}

        # Initial multiply factor is 1
        self.multiply_factor = 1

        # 1D slice color of 2D spectra is initially set to navy
        self.slice_colour = "navy"

        # Initial colour/reference/vertical index from list of colours is set to 0
        self.index = 0
        self.ref_index = 0
        self.vertical_index = 0

        # List to hold the multiple 2D spectra in multiplot mode
        self.twoD_spectra = []

        self.linewidth = 1.0
        self.linewidth1D = 1.5

        self.x_difference = 0
        self.y_difference = 0

        # Initially set the transpose flag to False
        self.transpose = False
        self.transposed2D = False

        # Default options for pivot point for P1 phasing
        self.pivot_x_default = 0
        self.pivot_x = self.pivot_x_default

        self.pivot_y_default = 0
        self.pivot_y = self.pivot_y_default

        self.slice_mode = None

        # Suppressing complex warning from numpy - prevents the complex warning from being printed to terminal when phasing
        import warnings

        # warnings.simplefilter("ignore", np.ComplexWarning)  # For old numpy versions
        warnings.simplefilter(
            "ignore", np.exceptions.ComplexWarning
        )  # For new numpy versions

    def UpdateFrame(self):
        self.canvas.draw()
        self.canvas.Refresh()
        self.canvas.Update()
        self.panel.Refresh()
        self.panel.Update()

    def OnSysColourChanged(self, event):

        pass


def main():
    app = wx.App()
    frame = SpinView()
    app.MainLoop()


if __name__ == "__main__":
    main()
