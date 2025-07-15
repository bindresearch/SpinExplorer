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


import wx
import warnings

warnings.simplefilter("ignore", UserWarning)


class FormatParametersVarian:
    def __init__(self, app, params, nmrdata) -> None:
        """
        This class contains the functions to populate the GUI with the
        determined parameters (stored in params) for Varian NMR data
        """
        self.app = app
        self.params = params
        self.nmrdata = nmrdata

    def input_sizes_varian(self):
        """
        Creating TextCtrl boxes to input the complex and read data
        dimension sizes
        """
        self.title_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.title_sizer.AddSpacer(250)
        self.N_complex_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.N_complex_txt = wx.StaticText(self, label="Number complex points:")
        self.N_complex_boxes = []
        if self.nmrdata.other_params == False:
            if self.nmrdata.phase == False and self.nmrdata.phase2 == False:
                for i in range(1):
                    self.title_sizer.Add(
                        wx.StaticText(self, label="Dimension " + str(i + 1))
                    )
                    self.title_sizer.AddSpacer(145)
                    if i == 0:
                        self.N_complex_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(self.nmrdata.size_direct),
                                size=(200, 20),
                            )
                        )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == False:
                for i in range(2):
                    self.title_sizer.Add(
                        wx.StaticText(self, label="Dimension " + str(i + 1))
                    )
                    self.title_sizer.AddSpacer(145)
                    if i == 0:
                        self.N_complex_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(self.nmrdata.size_direct),
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.N_complex_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(int(self.nmrdata.size_indirect[i - 1]) * 2),
                                size=(200, 20),
                            )
                        )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == True:
                for i in range(3):
                    self.title_sizer.Add(
                        wx.StaticText(self, label="Dimension " + str(i + 1))
                    )
                    self.title_sizer.AddSpacer(145)
                    if i == 0:
                        self.N_complex_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(self.nmrdata.size_direct),
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.N_complex_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(int(self.nmrdata.size_indirect[i - 1]) * 2),
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.N_complex_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(int(self.nmrdata.size_indirect[i - 1]) * 2),
                                size=(200, 20),
                            )
                        )
            elif self.nmrdata.phase == False and self.nmrdata.phase2 == True:
                for i in range(2):
                    self.title_sizer.Add(
                        wx.StaticText(self, label="Dimension " + str(i + 1))
                    )
                    self.title_sizer.AddSpacer(145)
                    if i == 0:
                        self.N_complex_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(self.nmrdata.size_direct),
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.N_complex_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(int(self.nmrdata.size_indirect[i - 1]) * 2),
                                size=(200, 20),
                            )
                        )
        else:
            if self.nmrdata.phase == False and self.nmrdata.phase2 == False:
                for i in range(2):
                    self.title_sizer.Add(
                        wx.StaticText(self, label="Dimension " + str(i + 1))
                    )
                    self.title_sizer.AddSpacer(145)
                    if i == 0:
                        self.N_complex_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(self.nmrdata.size_direct),
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.N_complex_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(self.nmrdata.number_of_arrayed_parameters),
                                size=(200, 20),
                            )
                        )

            elif self.nmrdata.phase == True and self.nmrdata.phase2 == False:
                for i in range(3):
                    self.title_sizer.Add(
                        wx.StaticText(self, label="Dimension " + str(i + 1))
                    )
                    self.title_sizer.AddSpacer(145)
                    if i == 0:
                        self.N_complex_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(self.nmrdata.size_direct),
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.N_complex_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(int(self.nmrdata.size_indirect[i - 1]) * 2),
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.N_complex_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(self.nmrdata.number_of_arrayed_parameters),
                                size=(200, 20),
                            )
                        )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == True:
                for i in range(4):
                    self.title_sizer.Add(
                        wx.StaticText(self, label="Dimension " + str(i + 1))
                    )
                    self.title_sizer.AddSpacer(145)
                    if i == 0:
                        self.N_complex_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(self.nmrdata.size_direct),
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.N_complex_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(int(self.nmrdata.size_indirect[i - 1]) * 2),
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.N_complex_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(int(self.nmrdata.size_indirect[i - 1]) * 2),
                                size=(200, 20),
                            )
                        )
                    if i == 3:
                        self.N_complex_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(self.nmrdata.number_of_arrayed_parameters),
                                size=(200, 20),
                            )
                        )
            elif self.nmrdata.phase == False and self.nmrdata.phase2 == True:
                for i in range(3):
                    self.title_sizer.Add(
                        wx.StaticText(self, label="Dimension " + str(i + 1))
                    )
                    self.title_sizer.AddSpacer(145)
                    if i == 0:
                        self.N_complex_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(self.nmrdata.size_direct),
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.N_complex_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(int(self.nmrdata.size_indirect[i - 1]) * 2),
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.N_complex_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(self.nmrdata.number_of_arrayed_parameters),
                                size=(200, 20),
                            )
                        )

        self.N_complex_sizer.AddSpacer(20)
        self.N_complex_sizer.Add(self.N_complex_txt)
        self.N_complex_sizer.AddSpacer(20)
        for i in range(0, len(self.N_complex_boxes)):
            self.N_complex_sizer.Add(self.N_complex_boxes[i])
            self.N_complex_sizer.AddSpacer(20)
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.title_sizer)
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.N_complex_sizer)

        # Create the boxes for the real point numbers
        self.N_real_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.N_real_txt = wx.StaticText(self, label="Number real points:        ")
        self.N_real_boxes = []
        if self.nmrdata.other_params == False:
            for i in range(len(self.N_complex_boxes)):
                if i == 0:
                    self.N_real_boxes.append(
                        wx.TextCtrl(
                            self,
                            value=str(int(self.nmrdata.size_direct / 2)),
                            size=(200, 20),
                        )
                    )
                else:
                    if i == 1:
                        if self.nmrdata.size_indirect[i - 1] == 1:
                            continue
                        else:
                            self.N_real_boxes.append(
                                wx.TextCtrl(
                                    self,
                                    value=str(self.nmrdata.size_indirect[i - 1]),
                                    size=(200, 20),
                                )
                            )
                    if i == 2:
                        if self.nmrdata.size_indirect[i - 1] == 1:
                            continue
                        else:
                            self.N_real_boxes.append(
                                wx.TextCtrl(
                                    self,
                                    value=str(self.nmrdata.size_indirect[i - 1]),
                                    size=(200, 20),
                                )
                            )

        else:
            if self.nmrdata.phase == False and self.nmrdata.phase2 == False:
                for i in range(2):
                    if i == 0:
                        self.N_real_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(int(self.nmrdata.size_direct / 2)),
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.N_real_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(self.nmrdata.number_of_arrayed_parameters),
                                size=(200, 20),
                            )
                        )

            elif self.nmrdata.phase == True and self.nmrdata.phase2 == False:
                for i in range(3):
                    if i == 0:
                        self.N_real_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(int(self.nmrdata.size_direct / 2)),
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.N_real_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(self.nmrdata.size_indirect[i - 1]),
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.N_real_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(self.nmrdata.number_of_arrayed_parameters),
                                size=(200, 20),
                            )
                        )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == True:
                for i in range(4):
                    if i == 0:
                        self.N_real_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(int(self.nmrdata.size_direct / 2)),
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.N_real_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(self.nmrdata.size_indirect[i - 1]),
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.N_real_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(self.nmrdata.size_indirect[i - 1]),
                                size=(200, 20),
                            )
                        )
                    if i == 3:
                        self.N_real_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(self.nmrdata.number_of_arrayed_parameters),
                                size=(200, 20),
                            )
                        )
            elif self.nmrdata.phase == False and self.nmrdata.phase2 == True:
                for i in range(3):
                    if i == 0:
                        self.N_real_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(int(self.nmrdata.size_direct / 2)),
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.N_real_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(int(self.nmrdata.size_indirect[i - 1] / 2)),
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.N_real_boxes.append(
                            wx.TextCtrl(
                                self,
                                value=str(self.nmrdata.number_of_arrayed_parameters),
                                size=(200, 20),
                            )
                        )

        self.N_real_sizer.AddSpacer(20)
        self.N_real_sizer.Add(self.N_real_txt)
        self.N_real_sizer.AddSpacer(20)
        for i in range(0, len(self.N_real_boxes)):
            self.N_real_sizer.Add(self.N_real_boxes[i])
            self.N_real_sizer.AddSpacer(20)
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.N_real_sizer)

    def input_acquisition_modes_varian(self):
        """
        Create drop down options for the acquisition modes
        """
        self.acquisition_mode_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.acquisition_mode_txt = wx.StaticText(
            self, label="Acquisition mode:           "
        )
        self.acquisition_mode_options_direct = ["Complex", "Sequential", "Real", "DQD"]
        self.acquisition_mode_options_indirect = [
            "Complex",
            "States-TPPI",
            "Rance-Kay",
            "Echo-AntiEcho",
            "TPPI",
            "States",
            "Real",
        ]

        self.acqusition_combo_boxes = []
        if self.nmrdata.other_params == False:
            if self.nmrdata.phase == False and self.nmrdata.phase2 == False:
                self.acqusition_combo_boxes.append(
                    wx.ComboBox(
                        self,
                        value=self.acquisition_mode_options_direct[0],
                        choices=self.acquisition_mode_options_direct,
                        size=(200, 20),
                    )
                )
                self.acqusition_combo_boxes[0].Bind(
                    wx.EVT_COMBOBOX, self.on_acquisition_mode_change
                )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == False:
                for i in range(2):
                    if i == 0:
                        self.acqusition_combo_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.acquisition_mode_options_direct[0],
                                choices=self.acquisition_mode_options_direct,
                                size=(200, 20),
                            )
                        )
                        self.acqusition_combo_boxes[i].Bind(
                            wx.EVT_COMBOBOX, self.on_acquisition_mode_change
                        )
                    if i == 1:
                        self.acqusition_combo_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.acquisition_mode_options_indirect[0],
                                choices=self.acquisition_mode_options_indirect,
                                size=(200, 20),
                            )
                        )
                        self.acqusition_combo_boxes[i].Bind(
                            wx.EVT_COMBOBOX, self.on_acquisition_mode_change
                        )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == True:
                for i in range(3):
                    if i == 0:
                        self.acqusition_combo_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.acquisition_mode_options_direct[0],
                                choices=self.acquisition_mode_options_direct,
                                size=(200, 20),
                            )
                        )
                        self.acqusition_combo_boxes[i].Bind(
                            wx.EVT_COMBOBOX, self.on_acquisition_mode_change
                        )
                    if i == 1:
                        self.acqusition_combo_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.acquisition_mode_options_indirect[0],
                                choices=self.acquisition_mode_options_indirect,
                                size=(200, 20),
                            )
                        )
                        self.acqusition_combo_boxes[i].Bind(
                            wx.EVT_COMBOBOX, self.on_acquisition_mode_change
                        )
                    if i == 2:
                        self.acqusition_combo_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.acquisition_mode_options_indirect[0],
                                choices=self.acquisition_mode_options_indirect,
                                size=(200, 20),
                            )
                        )
                        self.acqusition_combo_boxes[i].Bind(
                            wx.EVT_COMBOBOX, self.on_acquisition_mode_change
                        )
        else:
            if self.nmrdata.phase == False and self.nmrdata.phase2 == False:
                for i in range(2):
                    if i == 0:
                        self.acqusition_combo_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.acquisition_mode_options_direct[0],
                                choices=self.acquisition_mode_options_direct,
                                size=(200, 20),
                            )
                        )
                        self.acqusition_combo_boxes[i].Bind(
                            wx.EVT_COMBOBOX, self.on_acquisition_mode_change
                        )
                    if i == 1:
                        self.acqusition_combo_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.acquisition_mode_options_indirect[6],
                                choices=self.acquisition_mode_options_indirect,
                                size=(200, 20),
                            )
                        )
                        self.acqusition_combo_boxes[i].Bind(
                            wx.EVT_COMBOBOX, self.on_acquisition_mode_change
                        )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == False:
                for i in range(3):
                    if i == 0:
                        self.acqusition_combo_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.acquisition_mode_options_direct[0],
                                choices=self.acquisition_mode_options_direct,
                                size=(200, 20),
                            )
                        )
                        self.acqusition_combo_boxes[i].Bind(
                            wx.EVT_COMBOBOX, self.on_acquisition_mode_change
                        )
                    if i == 1:
                        self.acqusition_combo_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.acquisition_mode_options_indirect[0],
                                choices=self.acquisition_mode_options_indirect,
                                size=(200, 20),
                            )
                        )
                        self.acqusition_combo_boxes[i].Bind(
                            wx.EVT_COMBOBOX, self.on_acquisition_mode_change
                        )
                    if i == 2:
                        self.acqusition_combo_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.acquisition_mode_options_indirect[6],
                                choices=self.acquisition_mode_options_indirect,
                                size=(200, 20),
                            )
                        )
                        self.acqusition_combo_boxes[i].Bind(
                            wx.EVT_COMBOBOX, self.on_acquisition_mode_change
                        )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == True:
                for i in range(4):
                    if i == 0:
                        self.acqusition_combo_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.acquisition_mode_options_direct[0],
                                choices=self.acquisition_mode_options_direct,
                                size=(200, 20),
                            )
                        )
                        self.acqusition_combo_boxes[i].Bind(
                            wx.EVT_COMBOBOX, self.on_acquisition_mode_change
                        )
                    if i == 1:
                        self.acqusition_combo_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.acquisition_mode_options_indirect[0],
                                choices=self.acquisition_mode_options_indirect,
                                size=(200, 20),
                            )
                        )
                        self.acqusition_combo_boxes[i].Bind(
                            wx.EVT_COMBOBOX, self.on_acquisition_mode_change
                        )
                    if i == 2:
                        self.acqusition_combo_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.acquisition_mode_options_indirect[0],
                                choices=self.acquisition_mode_options_indirect,
                                size=(200, 20),
                            )
                        )
                        self.acqusition_combo_boxes[i].Bind(
                            wx.EVT_COMBOBOX, self.on_acquisition_mode_change
                        )
                    if i == 3:
                        self.acqusition_combo_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.acquisition_mode_options_indirect[6],
                                choices=self.acquisition_mode_options_indirect,
                                size=(200, 20),
                            )
                        )
                        self.acqusition_combo_boxes[i].Bind(
                            wx.EVT_COMBOBOX, self.on_acquisition_mode_change
                        )
            elif self.nmrdata.phase == False and self.nmrdata.phase2 == True:
                for i in range(3):
                    if i == 0:
                        self.acqusition_combo_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.acquisition_mode_options_direct[0],
                                choices=self.acquisition_mode_options_direct,
                                size=(200, 20),
                            )
                        )
                        self.acqusition_combo_boxes[i].Bind(
                            wx.EVT_COMBOBOX, self.on_acquisition_mode_change
                        )
                    if i == 1:
                        self.acqusition_combo_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.acquisition_mode_options_indirect[0],
                                choices=self.acquisition_mode_options_indirect,
                                size=(200, 20),
                            )
                        )
                        self.acqusition_combo_boxes[i].Bind(
                            wx.EVT_COMBOBOX, self.on_acquisition_mode_change
                        )
                    if i == 2:
                        self.acqusition_combo_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.acquisition_mode_options_indirect[6],
                                choices=self.acquisition_mode_options_indirect,
                                size=(200, 20),
                            )
                        )
                        self.acqusition_combo_boxes[i].Bind(
                            wx.EVT_COMBOBOX, self.on_acquisition_mode_change
                        )

        self.acquisition_mode_sizer.AddSpacer(20)
        self.acquisition_mode_sizer.Add(self.acquisition_mode_txt)
        self.acquisition_mode_sizer.AddSpacer(20)
        for i in range(0, len(self.acqusition_combo_boxes)):
            self.acquisition_mode_sizer.Add(self.acqusition_combo_boxes[i])
            self.acquisition_mode_sizer.AddSpacer(20)
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.acquisition_mode_sizer)

    def input_sweep_widths_varian(self):
        """
        Creating ComboBoxes for the spectral sweep widths for each dimension.
        """
        self.sweep_width_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.sweep_width_txt = wx.StaticText(self, label="Sweep width (Hz):          ")
        self.sweep_width_boxes = []
        options_sw = [str(self.nmrdata.sw_direct)]
        try:
            options_sw = options_sw + [str(self.nmrdata.sw_indirect["sw1"])]
            try:
                options_sw = options_sw + [str(self.nmrdata.sw_indirect["sw2"])]
                try:
                    options_sw = (
                        options_sw + [str(self.nmrdata.sw_indirect["sw3"])] + ["0.0"]
                    )
                except:
                    options_sw = options_sw + ["0.0"]
            except:
                options_sw = options_sw + ["0.0"]
        except:
            options_sw = options_sw + ["0.0"]
        if self.nmrdata.other_params == False:
            if self.nmrdata.phase == False and self.nmrdata.phase2 == False:
                self.sweep_width_boxes.append(
                    wx.ComboBox(
                        self,
                        value=str(self.nmrdata.sw_direct),
                        choices=options_sw,
                        size=(200, 20),
                    )
                )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == False:
                for i in range(2):
                    if i == 0:
                        self.sweep_width_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.sw_direct),
                                choices=options_sw,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.sweep_width_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.sw_indirect["sw1"]),
                                choices=options_sw,
                                size=(200, 20),
                            )
                        )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == True:
                for i in range(3):
                    if i == 0:
                        self.sweep_width_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.sw_direct),
                                choices=options_sw,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.sweep_width_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.sw_indirect["sw1"]),
                                choices=options_sw,
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.sweep_width_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.sw_indirect["sw2"]),
                                choices=options_sw,
                                size=(200, 20),
                            )
                        )
            elif self.nmrdata.phase == False and self.nmrdata.phase2 == True:
                for i in range(2):
                    if i == 0:
                        self.sweep_width_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.sw_direct),
                                choices=options_sw,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.sweep_width_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.sw_indirect["sw2"]),
                                choices=options_sw,
                                size=(200, 20),
                            )
                        )
        else:
            if self.nmrdata.phase == False and self.nmrdata.phase2 == False:
                for i in range(2):
                    if i == 0:
                        self.sweep_width_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.sw_direct),
                                choices=options_sw,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.sweep_width_boxes.append(
                            wx.ComboBox(
                                self, value=str(0), choices=options_sw, size=(200, 20)
                            )
                        )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == False:
                for i in range(3):
                    if i == 0:
                        self.sweep_width_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.sw_direct),
                                choices=options_sw,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.sweep_width_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.sw_indirect["sw1"]),
                                choices=options_sw,
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.sweep_width_boxes.append(
                            wx.ComboBox(
                                self, value=str(0), choices=options_sw, size=(200, 20)
                            )
                        )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == True:
                for i in range(4):
                    if i == 0:
                        self.sweep_width_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.sw_direct),
                                choices=options_sw,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.sweep_width_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.sw_indirect["sw1"]),
                                choices=options_sw,
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.sweep_width_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.sw_indirect["sw2"]),
                                choices=options_sw,
                                size=(200, 20),
                            )
                        )
                    if i == 3:
                        self.sweep_width_boxes.append(
                            wx.ComboBox(
                                self, value=str(0), choices=options_sw, size=(200, 20)
                            )
                        )
            elif self.nmrdata.phase == False and self.nmrdata.phase2 == True:
                for i in range(3):
                    if i == 0:
                        self.sweep_width_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.sw_direct),
                                choices=options_sw,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.sweep_width_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.sw_indirect["sw2"]),
                                choices=options_sw,
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.sweep_width_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.number_of_arrayed_parameters),
                                choices=options_sw,
                                size=(200, 20),
                            )
                        )
        self.sweep_width_sizer.AddSpacer(20)
        self.sweep_width_sizer.Add(self.sweep_width_txt)
        self.sweep_width_sizer.AddSpacer(20)
        for i in range(0, len(self.sweep_width_boxes)):
            self.sweep_width_sizer.Add(self.sweep_width_boxes[i])
            self.sweep_width_sizer.AddSpacer(20)
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.sweep_width_sizer)

    def get_nuclei_frequency_varian(self):
        """
        Inputting the nucleus frequencies of each dimension
        """
        self.nuclei_frequency_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.nuclei_frequency_txt = wx.StaticText(self, label="Nuclei frequency (MHz):")
        self.nuclei_frequency_boxes = []
        options = [str(self.nmrdata.nucleus_frequency_direct)] + [
            str(freq) for freq in self.nmrdata.nucleus_frequencies_indirect
        ]
        if self.nmrdata.other_params == False:
            if self.nmrdata.phase == False and self.nmrdata.phase2 == False:
                self.nuclei_frequency_boxes.append(
                    wx.ComboBox(
                        self,
                        value=str(self.nmrdata.nucleus_frequency_direct),
                        choices=options,
                        size=(200, 20),
                    )
                )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == False:
                for i in range(2):
                    if i == 0:
                        self.nuclei_frequency_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.nucleus_frequency_direct),
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        if self.nmrdata.nucleus_frequencies_indirect[0] == 0:
                            try:
                                self.nuclei_frequency_boxes.append(
                                    wx.ComboBox(
                                        self,
                                        value=str(
                                            self.nmrdata.nucleus_frequencies_indirect[1]
                                        ),
                                        choices=options,
                                        size=(200, 20),
                                    )
                                )
                            except:
                                self.nuclei_frequency_boxes.append(
                                    wx.ComboBox(
                                        self,
                                        value=str(1),
                                        choices=options,
                                        size=(200, 20),
                                    )
                                )
                        else:
                            self.nuclei_frequency_boxes.append(
                                wx.ComboBox(
                                    self,
                                    value=str(
                                        self.nmrdata.nucleus_frequencies_indirect[0]
                                    ),
                                    choices=options,
                                    size=(200, 20),
                                )
                            )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == True:
                for i in range(3):
                    if i == 0:
                        self.nuclei_frequency_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.nucleus_frequency_direct),
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.nuclei_frequency_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.nucleus_frequencies_indirect[0]),
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.nuclei_frequency_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.nucleus_frequencies_indirect[1]),
                                choices=options,
                                size=(200, 20),
                            )
                        )
            elif self.nmrdata.phase == False and self.nmrdata.phase2 == True:
                for i in range(2):
                    if i == 0:
                        self.nuclei_frequency_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.nucleus_frequency_direct),
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.nuclei_frequency_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.nucleus_frequencies_indirect[0]),
                                choices=options,
                                size=(200, 20),
                            )
                        )
        else:
            if self.nmrdata.phase == False and self.nmrdata.phase2 == False:
                for i in range(2):
                    if i == 0:
                        self.nuclei_frequency_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.nucleus_frequency_direct),
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.nuclei_frequency_boxes.append(
                            wx.ComboBox(
                                self, value="0.0", choices=options, size=(200, 20)
                            )
                        )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == False:
                for i in range(3):
                    if i == 0:
                        self.nuclei_frequency_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.nucleus_frequency_direct),
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.nuclei_frequency_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.nucleus_frequencies_indirect[0]),
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.nuclei_frequency_boxes.append(
                            wx.ComboBox(
                                self, value="0.0", choices=options, size=(200, 20)
                            )
                        )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == True:
                for i in range(4):
                    if i == 0:
                        self.nuclei_frequency_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.nucleus_frequency_direct),
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.nuclei_frequency_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.nucleus_frequencies_indirect[0]),
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.nuclei_frequency_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.nucleus_frequencies_indirect[1]),
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 3:
                        self.nuclei_frequency_boxes.append(
                            wx.ComboBox(
                                self, value="0.0", choices=options, size=(200, 20)
                            )
                        )
            elif self.nmrdata.phase == False and self.nmrdata.phase2 == True:
                for i in range(3):
                    if i == 0:
                        self.nuclei_frequency_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.nucleus_frequency_direct),
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.nuclei_frequency_boxes.append(
                            wx.ComboBox(
                                self,
                                value=str(self.nmrdata.nucleus_frequencies_indirect[0]),
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.nuclei_frequency_boxes.append(
                            wx.ComboBox(
                                self, value="0.0", choices=options, size=(200, 20)
                            )
                        )

        self.nuclei_frequency_sizer.AddSpacer(20)
        self.nuclei_frequency_sizer.Add(self.nuclei_frequency_txt)
        self.nuclei_frequency_sizer.AddSpacer(19)
        for i in range(0, len(self.nuclei_frequency_boxes)):
            self.nuclei_frequency_sizer.Add(self.nuclei_frequency_boxes[i])
            self.nuclei_frequency_sizer.AddSpacer(20)
        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.nuclei_frequency_sizer)

    def get_nuclei_labels_varian(self):
        """
        Creating boxes for the labels in each dimension
        """
        self.nucleus_type_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.nucleus_type_txt = wx.StaticText(
            self, label="Label:                              "
        )
        self.nucleus_type_boxes = []
        self.nmrdata.labels_indirect.remove("")
        options = (
            [self.nmrdata.label_direct]
            + [
                self.nmrdata.labels_indirect[i]
                for i in range(len(self.nmrdata.labels_indirect))
            ]
            + ["ID"]
        )
        if self.nmrdata.other_params == True:
            options = options + [self.nmrdata.arrayed_parameter]
        j = 0

        if self.nmrdata.other_params == False:
            if self.nmrdata.phase == False and self.nmrdata.phase2 == False:
                self.nucleus_type_boxes.append(
                    wx.ComboBox(
                        self,
                        value=self.nmrdata.label_direct,
                        choices=options,
                        size=(200, 20),
                    )
                )

            elif self.nmrdata.phase == True and self.nmrdata.phase2 == False:
                for i in range(2):
                    if i == 0:
                        self.nucleus_type_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.nmrdata.label_direct,
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.nucleus_type_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.nmrdata.labels_indirect[0],
                                choices=options,
                                size=(200, 20),
                            )
                        )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == True:
                for i in range(3):
                    if i == 0:
                        self.nucleus_type_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.nmrdata.label_direct,
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.nucleus_type_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.nmrdata.labels_indirect[0],
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.nucleus_type_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.nmrdata.labels_indirect[1],
                                choices=options,
                                size=(200, 20),
                            )
                        )
            elif self.nmrdata.phase == False and self.nmrdata.phase2 == True:
                for i in range(2):
                    if i == 0:
                        self.nucleus_type_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.nmrdata.label_direct,
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.nucleus_type_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.nmrdata.labels_indirect[1],
                                choices=options,
                                size=(200, 20),
                            )
                        )
            else:
                for i in range(2):
                    if i == 0:
                        self.nucleus_type_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.nmrdata.label_direct,
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.nucleus_type_boxes.append(
                            wx.ComboBox(
                                self, value="ID", choices=options, size=(200, 20)
                            )
                        )
        else:
            if self.nmrdata.phase == False and self.nmrdata.phase2 == False:
                for i in range(2):
                    if i == 0:
                        self.nucleus_type_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.nmrdata.label_direct,
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        try:
                            self.nucleus_type_boxes.append(
                                wx.ComboBox(
                                    self,
                                    value=self.nmrdata.arrayed_parameter,
                                    choices=options,
                                    size=(200, 20),
                                )
                            )
                        except:
                            self.nucleus_type_boxes.append(
                                wx.ComboBox(
                                    self, value="ID", choices=options, size=(200, 20)
                                )
                            )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == False:
                for i in range(3):
                    if i == 0:
                        self.nucleus_type_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.nmrdata.label_direct,
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.nucleus_type_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.nmrdata.labels_indirect[0],
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        try:
                            self.nucleus_type_boxes.append(
                                wx.ComboBox(
                                    self,
                                    value=self.nmrdata.arrayed_parameter,
                                    choices=options,
                                    size=(200, 20),
                                )
                            )
                        except:
                            self.nucleus_type_boxes.append(
                                wx.ComboBox(
                                    self, value="ID", choices=options, size=(200, 20)
                                )
                            )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == True:
                for i in range(4):
                    if i == 0:
                        self.nucleus_type_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.nmrdata.label_direct,
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.nucleus_type_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.nmrdata.labels_indirect[0],
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.nucleus_type_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.nmrdata.labels_indirect[1],
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 3:
                        try:
                            self.nucleus_type_boxes.append(
                                wx.ComboBox(
                                    self,
                                    value=self.nmrdata.arrayed_parameter,
                                    choices=options,
                                    size=(200, 20),
                                )
                            )
                        except:
                            self.nucleus_type_boxes.append(
                                wx.ComboBox(
                                    self, value="ID", choices=options, size=(200, 20)
                                )
                            )
            elif self.nmrdata.phase == False and self.nmrdata.phase2 == True:
                for i in range(3):
                    if i == 0:
                        self.nucleus_type_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.nmrdata.label_direct,
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.nucleus_type_boxes.append(
                            wx.ComboBox(
                                self,
                                value=self.nmrdata.labels_indirect[1],
                                choices=options,
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        try:
                            self.nucleus_type_boxes.append(
                                wx.ComboBox(
                                    self,
                                    value=self.nmrdata.arrayed_parameter,
                                    choices=options,
                                    size=(200, 20),
                                )
                            )
                        except:
                            self.nucleus_type_boxes.append(
                                wx.ComboBox(
                                    self, value="ID", choices=options, size=(200, 20)
                                )
                            )

        self.nucleus_type_sizer.AddSpacer(20)
        self.nucleus_type_sizer.Add(self.nucleus_type_txt)
        self.nucleus_type_sizer.AddSpacer(22)
        for i in range(0, len(self.nucleus_type_boxes)):
            self.nucleus_type_sizer.Add(self.nucleus_type_boxes[i])
            self.nucleus_type_sizer.AddSpacer(20)

        self.menu_bar.AddSpacer(10)
        self.menu_bar.Add(self.nucleus_type_sizer)

    def get_carrier_frequencies_varian(self):
        """
        Creating boxes for the carrier frequencies (ppm) in each
        dimension.
        """
        self.carrier_frequency_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.carrier_frequency_txt = wx.StaticText(
            self.app, label="Carrier frequency (ppm):"
        )
        self.carrier_frequency_boxes = []
        self.carrier_combo_boxes = []
        if self.nmrdata.label_direct == "1H" or self.nmrdata.label_direct == "H1":
            self.options_proton = ["H2O", "Other"]
        else:
            self.options_proton = ["Manual"]
        if self.nmrdata.size_indirect == []:
            self.options_other_dimensions = ["Other"]
        else:
            self.options_other_dimensions = self.nmrdata.references_other_labels + [
                "Other"
            ]

        id_columns = 0

        if self.nmrdata.other_params == False:
            if self.nmrdata.phase == False and self.nmrdata.phase2 == False:
                self.carrier_frequency_boxes.append(
                    wx.TextCtrl(
                        self.app,
                        value=str(self.nmrdata.references_proton[0]),
                        size=(200, 20),
                    )
                )
                self.carrier_combo_boxes.append(
                    wx.ComboBox(
                        self.app,
                        value=self.options_proton[0],
                        choices=self.options_proton,
                        size=(200, 20),
                    )
                )

            elif self.nmrdata.phase == True and self.nmrdata.phase2 == False:
                for i in range(2):
                    if i == 0:
                        self.carrier_frequency_boxes.append(
                            wx.TextCtrl(
                                self.app,
                                value=str(self.nmrdata.references_proton[0]),
                                size=(200, 20),
                            )
                        )
                        self.carrier_combo_boxes.append(
                            wx.ComboBox(
                                self.app,
                                value=self.options_proton[0],
                                choices=self.options_proton,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.carrier_frequency_boxes.append(
                            wx.TextCtrl(
                                self.app,
                                value=str(self.nmrdata.references_other[0]),
                                size=(200, 20),
                            )
                        )
                        self.carrier_combo_boxes.append(
                            wx.ComboBox(
                                self.app,
                                value=self.options_other_dimensions[0],
                                choices=self.options_other_dimensions,
                                size=(200, 20),
                            )
                        )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == True:
                for i in range(3):
                    if i == 0:
                        self.carrier_frequency_boxes.append(
                            wx.TextCtrl(
                                self.app,
                                value=str(self.nmrdata.references_proton[0]),
                                size=(200, 20),
                            )
                        )
                        self.carrier_combo_boxes.append(
                            wx.ComboBox(
                                self.app,
                                value=self.nmrdata.references_proton_labels[0],
                                choices=self.options_proton,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.carrier_frequency_boxes.append(
                            wx.TextCtrl(
                                self.app,
                                value=str(self.nmrdata.references_other[0]),
                                size=(200, 20),
                            )
                        )
                        self.carrier_combo_boxes.append(
                            wx.ComboBox(
                                self.app,
                                value=self.options_other_dimensions[0],
                                choices=self.options_other_dimensions,
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.carrier_frequency_boxes.append(
                            wx.TextCtrl(
                                self.app,
                                value=str(self.nmrdata.references_other[1]),
                                size=(200, 20),
                            )
                        )
                        self.carrier_combo_boxes.append(
                            wx.ComboBox(
                                self.app,
                                value=self.options_other_dimensions[1],
                                choices=self.options_other_dimensions,
                                size=(200, 20),
                            )
                        )
            elif self.nmrdata.phase == False and self.nmrdata.phase2 == True:
                for i in range(2):
                    if i == 0:
                        self.carrier_frequency_boxes.append(
                            wx.TextCtrl(
                                self.app,
                                value=str(self.nmrdata.references_proton[0]),
                                size=(200, 20),
                            )
                        )
                        self.carrier_combo_boxes.append(
                            wx.ComboBox(
                                self.app,
                                value=self.nmrdata.references_proton_labels[0],
                                choices=self.options_proton,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.carrier_frequency_boxes.append(
                            wx.TextCtrl(
                                self.app,
                                value=str(self.nmrdata.references_other[1]),
                                size=(200, 20),
                            )
                        )
                        self.carrier_combo_boxes.append(
                            wx.ComboBox(
                                self.app,
                                value=self.options_other_dimensions[1],
                                choices=self.options_other_dimensions,
                                size=(200, 20),
                            )
                        )

        else:
            if self.nmrdata.phase == False and self.nmrdata.phase2 == False:
                for i in range(2):
                    if i == 0:
                        self.carrier_frequency_boxes.append(
                            wx.TextCtrl(
                                self.app,
                                value=str(self.nmrdata.references_proton[0]),
                                size=(200, 20),
                            )
                        )
                        self.carrier_combo_boxes.append(
                            wx.ComboBox(
                                self.app,
                                value=self.nmrdata.references_proton_labels[0],
                                choices=self.options_proton,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        try:
                            self.carrier_frequency_boxes.append(
                                wx.TextCtrl(
                                    self.app,
                                    value=str(self.nmrdata.references_other[0]),
                                    size=(200, 20),
                                )
                            )
                            self.carrier_combo_boxes.append(
                                wx.ComboBox(
                                    self.app,
                                    value=self.options_other_dimensions[0],
                                    choices=self.options_other_dimensions,
                                    size=(200, 20),
                                )
                            )
                        except:
                            self.carrier_frequency_boxes.append(
                                wx.TextCtrl(self.app, value="0.0", size=(200, 20))
                            )
                            self.carrier_combo_boxes.append(
                                wx.ComboBox(
                                    self.app,
                                    value="ID",
                                    choices=self.options_other_dimensions,
                                    size=(200, 20),
                                )
                            )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == False:
                for i in range(3):
                    if i == 0:
                        self.carrier_frequency_boxes.append(
                            wx.TextCtrl(
                                self.app,
                                value=str(self.nmrdata.references_proton[0]),
                                size=(200, 20),
                            )
                        )
                        self.carrier_combo_boxes.append(
                            wx.ComboBox(
                                self.app,
                                value=self.nmrdata.references_proton_labels[0],
                                choices=self.options_proton,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.carrier_frequency_boxes.append(
                            wx.TextCtrl(
                                self.app,
                                value=str(self.nmrdata.references_other[0]),
                                size=(200, 20),
                            )
                        )
                        self.carrier_combo_boxes.append(
                            wx.ComboBox(
                                self.app,
                                value=self.options_other_dimensions[0],
                                choices=self.options_other_dimensions,
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.carrier_frequency_boxes.append(
                            wx.TextCtrl(self.app, value="0.0", size=(200, 20))
                        )
                        self.carrier_combo_boxes.append(
                            wx.ComboBox(
                                self.app,
                                value="Other",
                                choices=self.options_other_dimensions,
                                size=(200, 20),
                            )
                        )
            elif self.nmrdata.phase == True and self.nmrdata.phase2 == True:
                for i in range(4):
                    if i == 0:
                        self.carrier_frequency_boxes.append(
                            wx.TextCtrl(
                                self.app,
                                value=str(self.nmrdata.references_proton[0]),
                                size=(200, 20),
                            )
                        )
                        self.carrier_combo_boxes.append(
                            wx.ComboBox(
                                self.app,
                                value=self.nmrdata.references_proton_labels[0],
                                choices=self.options_proton,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.carrier_frequency_boxes.append(
                            wx.TextCtrl(
                                self.app,
                                value=str(self.nmrdata.references_other[0]),
                                size=(200, 20),
                            )
                        )
                        self.carrier_combo_boxes.append(
                            wx.ComboBox(
                                self.app,
                                value=self.options_other_dimensions[0],
                                choices=self.options_other_dimensions,
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.carrier_frequency_boxes.append(
                            wx.TextCtrl(
                                self.app,
                                value=str(self.nmrdata.references_other[1]),
                                size=(200, 20),
                            )
                        )
                        self.carrier_combo_boxes.append(
                            wx.ComboBox(
                                self.app,
                                value=self.options_other_dimensions[1],
                                choices=self.options_other_dimensions,
                                size=(200, 20),
                            )
                        )
                    if i == 3:
                        self.carrier_frequency_boxes.append(
                            wx.TextCtrl(self.app, value="0.0", size=(200, 20))
                        )
                        self.carrier_combo_boxes.append(
                            wx.ComboBox(
                                self.app,
                                value="Other",
                                choices=self.options_other_dimensions,
                                size=(200, 20),
                            )
                        )
            elif self.nmrdata.phase == False and self.nmrdata.phase2 == True:
                for i in range(3):
                    if i == 0:
                        self.carrier_frequency_boxes.append(
                            wx.TextCtrl(
                                self.app,
                                value=str(self.nmrdata.references_proton[0]),
                                size=(200, 20),
                            )
                        )
                        self.carrier_combo_boxes.append(
                            wx.ComboBox(
                                self.app,
                                value=self.nmrdata.references_proton_labels[0],
                                choices=self.options_proton,
                                size=(200, 20),
                            )
                        )
                    if i == 1:
                        self.carrier_frequency_boxes.append(
                            wx.TextCtrl(
                                self.app,
                                value=str(self.nmrdata.references_other[1]),
                                size=(200, 20),
                            )
                        )
                        self.carrier_combo_boxes.append(
                            wx.ComboBox(
                                self.app,
                                value=self.options_other_dimensions[1],
                                choices=self.options_other_dimensions,
                                size=(200, 20),
                            )
                        )
                    if i == 2:
                        self.carrier_frequency_boxes.append(
                            wx.TextCtrl(self.app, value="0.0", size=(200, 20))
                        )
                        self.carrier_combo_boxes.append(
                            wx.ComboBox(
                                self.app,
                                value="Other",
                                choices=self.options_other_dimensions,
                                size=(200, 20),
                            )
                        )

        self.carrier_frequency_sizer.AddSpacer(20)
        self.carrier_frequency_sizer.Add(self.carrier_frequency_txt)
        self.carrier_frequency_sizer.AddSpacer(16)
        for i in range(0, len(self.carrier_frequency_boxes)):
            self.carrier_frequency_sizer.Add(self.carrier_frequency_boxes[i])
            self.carrier_frequency_sizer.AddSpacer(20)

        self.app.menu_bar.AddSpacer(10)
        self.app.menu_bar.Add(self.carrier_frequency_sizer)

        self.carrier_combo_box_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.carrier_combo_box_sizer.AddSpacer(20)
        self.carrier_combo_box_txt = wx.StaticText(
            self.app, label="Referencing Mode:          "
        )
        self.carrier_combo_box_sizer.Add(self.carrier_combo_box_txt)
        self.carrier_combo_box_sizer.AddSpacer(17)
        for i in range(0, len(self.carrier_combo_boxes)):
            if self.carrier_combo_boxes[i] == 200:
                self.carrier_combo_box_sizer.AddSpacer(200)
                self.carrier_combo_box_sizer.AddSpacer(20)
            else:
                self.carrier_combo_box_sizer.Add(self.carrier_combo_boxes[i])
                self.carrier_combo_box_sizer.AddSpacer(20)
                self.carrier_combo_boxes[i].Bind(
                    wx.EVT_COMBOBOX, self.on_carrier_combo_box_change_varian
                )

        self.app.menu_bar.AddSpacer(10)
        self.app.menu_bar.Add(self.carrier_combo_box_sizer)

    def on_carrier_combo_box_change_varian(self, event):
        """
        If the carrier combo box is changed, update the carrier
        accordingly
        """
        # Find the index of the combo box that was changed
        index = self.carrier_combo_boxes.index(event.GetEventObject())
        # Find the value of the combo box
        value = event.GetEventObject().GetValue()
        index_selection = self.carrier_combo_boxes[index].GetSelection()

        if index == 0:
            if self.nmrdata.label_direct == "1H" or self.nmrdata.label_direct == "H1":
                if value == "H2O":
                    self.carrier_frequency_boxes[index].SetValue(
                        str(self.nmrdata.references_proton[0])
                    )
                elif value == "Other":
                    self.carrier_frequency_boxes[index].SetValue(
                        str(self.nmrdata.references_proton[1])
                    )
        else:
            if self.nmrdata.label_direct == "1H" or self.nmrdata.label_direct == "H1":
                if value == "Other":
                    self.carrier_frequency_boxes[index].SetValue(str(0.0))
                else:
                    self.carrier_frequency_boxes[index].SetValue(
                        str(self.nmrdata.references_other[index_selection])
                    )
