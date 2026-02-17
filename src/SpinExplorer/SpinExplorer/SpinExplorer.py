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
print("                        SpinExplorer                         ")
print("-------------------------------------------------------------")
print("               (version 1.2) 17th February 2026              ")
print(" (c) 2025 James Eaton, Andrew Baldwin (University of Oxford) ")
print("                  2025-2026, Bind Research                   ")
print("                        MIT License                          ")
print("-------------------------------------------------------------")
print(" Video tutorials at:")
print(" https://www.youtube.com/@BindResearch")
print("-------------------------------------------------------------")
print("")


import sys

import wx
import wx.adv
import os
import json
import wx.lib.agw.hyperlink as hl
from appdirs import user_data_dir

# cache_dir = os.path.expanduser("~/.SpinExplorer_mpl_cache")
# os.makedirs(cache_dir, exist_ok=True)

appname = 'SpinExplorer'
appauthor = "James Eaton"
data_dir = user_data_dir(appname, appauthor)
os.makedirs(data_dir, exist_ok=True)


import matplotlib
mpl_cache = os.path.join(data_dir, "mpl-cache")
os.makedirs(mpl_cache, exist_ok=True)
matplotlib.get_cachedir = lambda: mpl_cache

import pathlib


# Importing buttons and images
from SpinExplorer.SpinExpLogo import SpinExpLogo
from SpinExplorer.SpinExplorer.SpinExplorerHeader import SpinExplorerHeader
from SpinExplorer.SpinExplorer.SpinConverterButton import SpinConverterButton
from SpinExplorer.SpinExplorer.SpinProcessButton import SpinProcessButton
from SpinExplorer.SpinExplorer.SpinViewButton import SpinViewButton
from SpinExplorer.SpinExplorer.Logo import Logo

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


class SpinExplorer(wx.Frame):
    def __init__(
        self
    ):
        # Get the monitor size and set the window size to 85% of the monitor size
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        self.width = 840
        self.height = 800


        # Setup the dock/task bar with the logo
        try:
            self.tbicon = TaskBarIcon(self)
        except:
            pass

        self.main_window = wx.Frame.__init__(
            self, None, wx.ID_ANY,'SpinExplorer', wx.DefaultPosition, size=(self.width, self.height)
        )

        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.main_sizer.AddSpacer(10)
        self.create_main_sizer()

        self.package_directory = user_data_dir(appname, appauthor)
        os.makedirs(self.package_directory, exist_ok=True)
        self.recently_opened_file = os.path.join(self.package_directory, "recently_opened.json")
        self.recent_directories = self.load_recent_directories()

        self.SetSizer(self.main_sizer)
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self.SetIcon(SpinExpLogo.GetIcon())
        self.Show()
        self.Centre()

        


    def create_main_sizer(self):
        
        bmp = SpinExplorerHeader.GetBitmap()
        top = wx.StaticBitmap(self, -1, bitmap=bmp)
        self.main_sizer.Add(top, 0, wx.ALIGN_CENTER_HORIZONTAL, 10)


        button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        bmp1 = SpinConverterButton.GetBitmap()
        left = wx.BitmapButton(self, -1, bitmap=bmp1)
        left.Bind(wx.EVT_BUTTON, self.OnClickSpinConverter)
        button_sizer.Add(left, 0, wx.ALIGN_CENTER_VERTICAL, 10)

        button_sizer.AddSpacer(20)

        bmp2 = SpinProcessButton.GetBitmap()
        middle = wx.BitmapButton(self, -1, bitmap=bmp2)
        middle.Bind(wx.EVT_BUTTON, self.OnClickSpinProcess)
        button_sizer.Add(middle, 0, wx.ALIGN_CENTER_VERTICAL, 10)

        button_sizer.AddSpacer(20)

        bmp3 = SpinViewButton.GetBitmap()
        right = wx.BitmapButton(self, -1, bitmap=bmp3)
        right.Bind(wx.EVT_BUTTON, self.OnClickSpinView)
        button_sizer.Add(right, 0, wx.ALIGN_CENTER_VERTICAL, 10)

        self.main_sizer.Add(button_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL, 10)


        # Add a listbox which shows recently opened directories

        self.load_spectra_box_label = wx.StaticBox(self, -1, "Select data directory:")
        self.load_spectra_box = wx.StaticBoxSizer(self.load_spectra_box_label, wx.HORIZONTAL)

        self.listbox = wx.ListBox(self, style=wx.LB_SINGLE, size=(300,100))

        self.listbox.Bind(wx.EVT_LISTBOX, self.OnListBoxPress)

        self.left_box = wx.BoxSizer(wx.VERTICAL)
        text = wx.StaticText(self, label="Recently opened:")

        self.left_box.Add(text)
        self.left_box.AddSpacer(5)
        self.left_box.Add(self.listbox)

        self.load_spectra_box.Add(self.left_box, 1, wx.ALIGN_CENTER_VERTICAL)
        self.load_spectra_box.AddSpacer(20)

        self.rightbox = wx.BoxSizer(wx.VERTICAL)

        self.file_button = wx.Button(self, -1, "Open File Browser", size=(300,20))
        self.file_button.Bind(wx.EVT_BUTTON, self.OnOpenFileBrowser)
        self.rightbox.AddSpacer(10)
        self.rightbox.Add(self.file_button)

        self.rightbox.AddSpacer(10)
        text2 = wx.StaticText(self, label="Selected Directory:")
        self.rightbox.Add(text2)
        self.directory_box = wx.TextCtrl(self, -1, '', size=(300,20), style=wx.CB_READONLY)
        self.rightbox.AddSpacer(10)
        self.rightbox.Add(self.directory_box)

        self.rightbox.AddSpacer(10)

        self.text3 = wx.StaticText(self, label="")
        
        self.rightbox.Add(self.text3)

        self.title_text = wx.StaticText(self, label="")

        self.rightbox.AddSpacer(5)
        self.rightbox.Add(self.title_text)


        self.load_spectra_box.Add(self.rightbox, 1, wx.ALIGN_CENTER_VERTICAL)



        self.main_sizer.AddSpacer(20)
        
        self.main_sizer.Add(self.load_spectra_box, 1, wx.ALIGN_CENTER_HORIZONTAL)



        # Add a bottom box with citations etc
        self.bottom_box = wx.BoxSizer(wx.VERTICAL)

        self.extra_info = wx.StaticText(self, label="This package was developed by Bind Research and the Baldwin Group at the University of Oxford.")

        self.video_link = hl.HyperLinkCtrl(
            self,
            -1,
            "Video Tutorials",
            URL="https://www.youtube.com/@BindResearch",
        )

        self.source_code_link = hl.HyperLinkCtrl(
            self,
            -1,
            "Source code",
            URL="https://github.com/bindresearch/SpinExplorer",
        )

        
        self.citation_info = wx.StaticText(self, label="If you use SpinExplorer, please cite the following:")

        self.citation1 = hl.HyperLinkCtrl(
            self,
            -1,
            "SpinExplorer",
            URL="xxx",
        )

        self.citation2 = hl.HyperLinkCtrl(
            self,
            -1,
            "nmrglue",
            URL="https://doi.org/10.1007/s10858-013-9718-x",
        )

        self.citation3 = hl.HyperLinkCtrl(
            self,
            -1,
            "nmrPipe",
            URL="https://doi.org/10.1007/BF00197809",
        )

        self.citation_box = wx.BoxSizer(wx.HORIZONTAL)

        self.copyright_statement = wx.StaticText(self, label="(c) 2025, James Eaton, Andrew Baldwin - University of Oxford")
        self.copyright_statement1 = wx.StaticText(self, label="(c) 2025-2026, Bind Research (Version 1.2)")

        self.citation_box.Add(self.citation1)
        self.citation_box.AddSpacer(10)
        self.citation_box.Add(self.citation2)
        self.citation_box.AddSpacer(10)
        self.citation_box.Add(self.citation3)

        self.link_box = wx.BoxSizer(wx.HORIZONTAL)
        self.link_box.Add(self.video_link)
        self.link_box.AddSpacer(10)
        self.link_box.Add(self.source_code_link)


        self.bottom_box.Add(self.extra_info, 1, wx.ALIGN_CENTER_HORIZONTAL)
        self.bottom_box.Add(self.link_box, 1, wx.ALIGN_CENTER_HORIZONTAL)
        self.bottom_box.Add(self.citation_info, 1, wx.ALIGN_CENTER_HORIZONTAL)
        self.bottom_box.Add(self.citation_box, 1, wx.ALIGN_CENTER_HORIZONTAL)

        self.bottom_box.Add(self.copyright_statement, 1, wx.ALIGN_CENTER_HORIZONTAL)
        self.bottom_box.Add(self.copyright_statement1, 1, wx.ALIGN_CENTER_HORIZONTAL)

        bmp = Logo.GetBitmap()
        logo = wx.StaticBitmap(self, -1, bitmap=bmp)
        

        self.main_sizer.AddSpacer(5)
        self.main_sizer.Add(self.bottom_box, 1, wx.ALIGN_CENTER_HORIZONTAL)
        self.main_sizer.Add(logo, 0, wx.ALIGN_CENTER_HORIZONTAL)


    def OnOpenFileBrowser(self, event):
        with wx.DirDialog(self, "Choose a data directory (containing the FID file)", style=wx.DD_DEFAULT_STYLE) as dirDialog:
            if dirDialog.ShowModal() == wx.ID_CANCEL:
                return 
            
            try:
                dir = dirDialog.GetPath()
                self.directory = pathlib.Path(dir)
                dirs = self.directory.parts[-3:]
                try:
                    last_directories_path = str(pathlib.Path(*dirs))
                except:
                    last_directories_path = self.directory

                
                self.directory_box.SetValue(last_directories_path)
                self.find_title()
                # Change the current directory to the one which was found
                os.chdir(self.directory)
                self.add_recent_dir(self.directory)

            except:
                # Give the user an error to say that the directory selected was not loaded correctly
                dlg = wx.MessageDialog(
                self,
                "The selected directory was not loaded correctly. Please select a different directory and try again.",
                "Warning",
                wx.OK | wx.ICON_WARNING,
                )
                dlg.ShowModal()
                dlg.Destroy()
                return



    def find_title(self):

        try:

            # Read the pdata/1/title to get the name
            title_file = self.directory / 'pdata/1/title'

            title = ''

            with open(title_file) as file:
                line = file.readlines()[0]
                title_extra = ''
                line = line.split('\n')[0]
                title_extra+= line 
                
            title = title + title_extra

            self.title_text.SetLabel(title)

            # Update the text in text3
            self.text3.SetLabel('Title for selected experiment:')

        except:
            # Unable to read the title
            self.title_text.SetLabel('No title found')
            self.text3.SetLabel('Title for selected experiment:')
            pass


    def OnListBoxPress(self, event):
        """
        Update the current selected directionary with the selected value
        """
        index = self.listbox.GetSelection()
        self.directory = pathlib.Path(self.recent_directories[index])



        self.directory_box.SetValue(self.recent_directories_short[index])
        self.find_title()
        # Change the current directory to the one which was found
        os.chdir(self.directory)
        self.add_recent_dir(self.directory)



    def load_recent_directories(self):

        if os.path.exists(self.recently_opened_file):
            with open(self.recently_opened_file, "r") as f:
                self.recent_directories = json.load(f)
                self.recent_directories_short = []
                for val in self.recent_directories:
                    dir = pathlib.Path(val)
                    dirs = dir.parts[-3:]
                    self.recent_directories_short.append(str(pathlib.Path(*dirs)))

                
                self.listbox.SetItems(self.recent_directories_short)
                return self.recent_directories
        else:
            return []
        

    def save_recent_dirs(self, recent_dirs):
        with open(self.recently_opened_file, "w") as f:
            json.dump(recent_dirs, f, indent=4)

    def add_recent_dir(self, path):
        path = os.path.abspath(path)
        recent_directories = self.recent_directories

        if path in recent_directories:
            recent_directories.remove(path)

        recent_directories.insert(0, path)
        recent_directories = recent_directories[:20]

        self.save_recent_dirs(recent_directories)

        self.recent_directories_short = []
        for val in recent_directories:
            dir = pathlib.Path(val)
            dirs = dir.parts[-3:]
            self.recent_directories_short.append(str(pathlib.Path(*dirs)))

        
        self.listbox.SetItems(self.recent_directories_short)


    def OnClickSpinConverter(self, event):
        from SpinExplorer.SpinConverter.SpinConverter import SpinConverter
        converting_frame = SpinConverter(explorer=True)


    def OnClickSpinProcess(self, event):
        from SpinExplorer.SpinProcess.SpinProcess import SpinProcess
        processing_frame = SpinProcess(explorer=True)
        

    def OnClickSpinView(self, event):
        from SpinExplorer.SpinView.SpinView import SpinView
        viewing_frame = SpinView(explorer=True)
    

    def OnClose(self, event):
        self.Destroy()
        sys.exit()

    


def main():
    app = wx.App()
    frame = SpinExplorer()
    app.MainLoop()


if __name__ == "__main__":
    main()
