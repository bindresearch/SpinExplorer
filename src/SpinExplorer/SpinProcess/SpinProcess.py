#!/usr/bin/env python3

"""MIT License

Copyright (c) 2025 James Eaton, Andrew Baldwin
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
print("                         SpinProcess                         ")
print("-------------------------------------------------------------")
print("                (version 2.0) 29th July 2025                 ")
print(" (c) 2025 James Eaton, Andrew Baldwin (University of Oxford) ")
print("                     2025, Bind Research                     ")
print("                        MIT License                          ")
print("-------------------------------------------------------------")
print("                     Processing NMR Data                     ")
print("-------------------------------------------------------------")
print(" Documentation at:")
print(" https://github.com/james-eaton-1/SpinExplorer")
print("-------------------------------------------------------------")
print("")


import sys

import wx
import wx.adv

# Import relevant modules
import numpy as np
import matplotlib

matplotlib.use("WXAgg")

import pathlib

matplotlib.rcParams["font.sans-serif"] = "Arial"
matplotlib.rcParams["font.family"] = "sans-serif"

# Suppress complex warning from numpy
import warnings

# warnings.simplefilter("ignore", np.ComplexWarning)  # For old numpy versions
warnings.simplefilter("ignore", np.exceptions.ComplexWarning)  # For new numpy versions


# Importing SpinProcess modules
from SpinExplorer.SpinProcess.ReadingData.read_fid import ReadFID
from SpinExplorer.SpinProcess.FormattingGUI.notebook import NotebookProcess
from SpinExplorer.SpinExpLogo import SpinExpLogo

# Find out the version of operating system being used (Mac, Linux, Windows)
if sys.platform == "linux":
    platform = "linux"
    height = 30
elif sys.platform == "darwin":
    platform = "mac"
    height = 16
else:
    platform = "windows"
    height = 30


# James Eaton, 10/06/2025, University of Oxford
# James Eaton, 25/09/2025, Bind Research
# This program is designed to allow the user to process NMR FID data that has been converted to nmrPipe format.


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
        self.SetIcon(icon, "SpinProcess")
        self.imgidx = 1      
        # bind some events
        self.Bind(wx.EVT_MENU, self.OnTaskBarClose, id=self.TBMENU_CLOSE)

    def CreatePopupMenu(self):
        """
        This method is called by the base class when it needs to popup
        the menu for the default EVT_RIGHT_DOWN event.
        """
        menu = wx.Menu()
        menu.Append(self.TBMENU_CLOSE,   "Close SpinProcess")
        return menu

    def MakeIcon(self, img):
        icon = wx.Icon()
        icon.CopyFromBitmap(img.ConvertToBitmap())
        return icon

    def OnTaskBarClose(self, evt):
        self.frame.Destroy()
        sys.exit()


class SpinProcess(wx.Frame):
    def __init__(
        self, original_frame=None, file_parser=False, path="", cwd="", reprocess=False
    ):
        # Get the monitor size and set the window size to 85% of the monitor size
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 0.8 * self.monitorWidth
        self.height = 0.75 * self.monitorHeight

        # Initially set the reprocessing flag to False
        self.reprocess = reprocess
        self.original_frame = original_frame
        self.file_parser = file_parser
        self.path = path
        self.cwd = cwd

        # Get the title for the panels
        self.title = self.GetTitle()

        # Setup the dock/task bar with the logo
        self.tbicon = TaskBarIcon(self)

        # Create the main window
        self.main_window = wx.Frame.__init__(
            self, None, title=self.title, size=(self.width, self.height)
        )

        # Read the NMR data in the current directory
        self.nmr_data = ReadFID(self)

        self.notebook = NotebookProcess(self, self.nmr_data)
        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.main_sizer.AddSpacer(10)
        self.main_sizer.Add(self.notebook, 1, wx.EXPAND)
        self.notebook.create_buttons(parent=self)

        self.SetSizerAndFit(self.main_sizer)
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
    
    def GetTitle(self):
        """
        Finding an appropriate title for the panel. 
        The title for the panel is:
        SpinProcess + the current working directory (last 3 elements)
        + the title + pulseprogram (for Bruker data)
        """
        title = 'SpinProcess: '
        p = pathlib.Path.cwd()
        dirs = p.parts[-3:]
        last_directories_path = pathlib.Path(*dirs)
        title = title + "/" + str(last_directories_path)
        
        # If pdata/1/title exists, add this title too
        try:
            with open('pdata/1/title') as file:
                lines = file.readlines()
                title = title + ' ('
                for line in lines:
                    line = line.split('\n')[0]
                    line = line + ' '
                    title+= line 
                
                title +=')'
        except:
            pass

        return title


def main():
    app = wx.App()
    frame = SpinProcess()
    app.MainLoop()


if __name__ == "__main__":
    main()
