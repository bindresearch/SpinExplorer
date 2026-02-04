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
import pathlib
import numpy as np
import nmrglue as ng


class Add_fid(wx.Frame):
    def __init__(self, app, nmrdata) -> None:
        """
        This class will open up a frame to select fids to add together
        """
        self.app = app
        self.nmrdata = nmrdata

        # Setting some initial variables
        self.fid_paths_short = []
        self.fid_paths_long = []
        self.selected_index = ''

        self.add_fid_frame = wx.Frame.__init__(
            self,
            None,
            wx.ID_ANY,
            "Add FIDs",
            wx.DefaultPosition,
            size=(400,400))
        
        self.create_sizer()

        self.SetSizer(self.main_sizer)
        self.Show()



    def create_sizer(self) -> None:
        """
        Creating the sizer for the add fid window
        """

        self.main_sizer = wx.BoxSizer(wx.VERTICAL)

        self.select_fid_text = wx.StaticText(self, -1, "Select multiple converted FIDs (.fid) and then press 'Add FIDs' to continue")
        self.select_fid_button = wx.Button(self, -1, 'Select FID (.fid)')
        self.select_fid_button.Bind(wx.EVT_BUTTON, self.on_select_fid)

        self.remove_fid_button = wx.Button(self, -1, 'Remove FID')
        self.remove_fid_button.Bind(wx.EVT_BUTTON, self.on_remove_fid)

        self.main_sizer.AddSpacer(10)
        self.main_sizer.Add(self.select_fid_text, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.main_sizer.AddSpacer(10)

        self.top_sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.top_sizer.Add(self.select_fid_button)
        self.top_sizer.AddSpacer(10)
        self.top_sizer.Add(self.remove_fid_button)

        self.main_sizer.Add(self.top_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)


        self.listbox = wx.ListBox(self, style=wx.LB_SINGLE, size=(300,100))
        self.listbox.Bind(wx.EVT_LISTBOX, self.on_listbox_press)

        self.main_sizer.AddSpacer(10)
        self.main_sizer.Add(self.listbox, 0, wx.ALIGN_CENTER_HORIZONTAL)


        self.add_fid_button = wx.Button(self, -1, 'Add FIDs')
        self.add_fid_button.Bind(wx.EVT_BUTTON, self.on_add_fid)

        self.main_sizer.AddSpacer(10)
        self.main_sizer.Add(self.add_fid_button, 0, wx.ALIGN_CENTER_HORIZONTAL)

        self.add_fid_text = wx.StaticText(self, -1, "Added FID will be saved in the current directory as test_add.fid")

        self.main_sizer.AddSpacer(10)
        self.main_sizer.Add(self.add_fid_text, 0, wx.ALIGN_CENTER_HORIZONTAL)


    def on_select_fid(self, event):
        with wx.FileDialog(self, "Choose an fid file (.fid)", style=wx.DD_DEFAULT_STYLE, wildcard='.fid') as FileDialog:
            if FileDialog.ShowModal() == wx.ID_CANCEL:
                return 
            
            try:
                file = FileDialog.GetPath()
                self.directory = pathlib.Path(file)
                dirs = self.directory.parts[-3:]
                try:
                    last_directories_path = str(pathlib.Path(*dirs))
                except:
                    last_directories_path = self.directory

                self.fid_paths_long.append(self.directory)
                self.fid_paths_short.append(last_directories_path)

                self.listbox.SetItems(self.fid_paths_short)
                


            except:
                # Give the user an error to say that the directory selected was not loaded correctly
                dlg = wx.MessageDialog(
                self,
                "The selected fid was not loaded correctly. Please select a different fid and try again.",
                "Warning",
                wx.OK | wx.ICON_WARNING,
                )
                dlg.ShowModal()
                dlg.Destroy()
                return

    def on_remove_fid(self, event):
        if(self.selected_index!=''):
            self.fid_paths_long.pop(self.selected_index)
            self.fid_paths_short.pop(self.selected_index)
            self.listbox.SetItems(self.fid_paths_short)


    def on_listbox_press(self, event):
        self.selected_index = self.listbox.GetSelection()

    def on_add_fid(self, event):
        """
        The FID addition works as follows:
        - First the code checks if the FIDs are of the same shape (they need to be the same to be added together)
        - Then the FIDs will be added together
        - Then the user is prompted to choose a directory and name to save the added FID as
        """

        if(len(self.fid_paths_long)<2):
            dlg = wx.MessageDialog(
                self,
                "Less than 2 FIDs have been selected. Please select multiple FIDs and then try again.",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return

        dic_list = []
        data_list = []
        data_list_shapes = []
        for file in self.fid_paths_long:
            dic, data = ng.pipe.read(file)
            dic_list.append(dic)
            data_list.append(data)
            data_list_shapes.append(data.shape)

        sets = set(data_list_shapes)

        if(len(sets)==1):
            # All FIDs have the same shape so can continue
            for i, data in enumerate(data_list):
                if(i==0):
                    data_added = data
                else:
                    data_added = data_added+data

        else:
            message = "Not all datasets have the same shape. Please ensure converted FIDs have identical shapes before adding them together.\n"
            for i, fid_path in enumerate(self.fid_paths_short):
                shape = np.array_str(data_list_shapes[i])
                message += fid_path + ': ' + shape + '\n'
            dlg = wx.MessageDialog(
                self,
                message,
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return
        
        if(os.path.exists('./test_add.fid')):
            dlg = wx.MessageDialog(
                self,
                "The file test_add.fid already exists. Continuing will overwrite this, would you like to continue?",
                "Warning",
                wx.YES_NO | wx.ICON_WARNING,
            )
            result = dlg.ShowModal()
            if(result==wx.ID_NO):
                dlg.Destroy()
                return
            dlg.Destroy()

        

        ng.pipe.write('./test_add.fid', dic, data, overwrite=True)
        

        

 