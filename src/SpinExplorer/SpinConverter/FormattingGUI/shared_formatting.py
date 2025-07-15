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


import os
import wx
import warnings

warnings.simplefilter("ignore", UserWarning)


class SharedFormatting:
    def __init__(self, app, params, nmrdata) -> None:
        """
        This class contains the functions to populate the GUI.
        These functions are general for both Bruker and Varian
        data formats.
        """
        self.app = app
        self.params = params
        self.nmrdata = nmrdata

    def create_temperature_box(self):
        self.temperature_box_label = wx.StaticBox(
            self.app, -1, label="Experiment Temperature:"
        )
        self.temperature_box = wx.StaticBoxSizer(
            self.temperature_box_label, wx.VERTICAL
        )
        self.text = wx.StaticText(self.app, label="Temp (K):")
        self.temperature_value = str(self.params.temperature)
        self.row_1 = wx.BoxSizer(wx.HORIZONTAL)
        self.temperature_parameter = wx.StaticText(
            self.app, label=self.temperature_value
        )
        self.row_1.AddSpacer(5)
        self.row_1.Add(self.text)
        self.row_1.AddSpacer(5)
        self.row_1.Add(self.temperature_parameter)
        self.temperature_input_button = wx.Button(self.app, label="Change Temperature")
        self.temperature_input_button.Bind(
            wx.EVT_BUTTON, self.on_temperature_change_button
        )
        self.temperature_box.Add(self.row_1)
        self.temperature_box.AddSpacer(10)
        self.temperature_box.Add(self.temperature_input_button)

        self.app.extra_boxes_0.AddSpacer(10)
        self.app.extra_boxes_0.Add(self.temperature_box)
        self.app.extra_boxes_0.AddSpacer(20)

        self.app.extra_sizers.AddSpacer(20)
        self.app.extra_sizers.Add(self.app.extra_boxes_total)

    def on_temperature_change_button(self, event):
        # Create a popout window to change the temperature
        self.temperature_change_window = wx.Frame(
            self.app,
            wx.ID_ANY,
            "Change Temperature",
            wx.DefaultPosition,
            size=(300, 120),
        )

        self.temperature_change_window_sizer = wx.BoxSizer(wx.VERTICAL)
        self.temperature_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.temperature_label = wx.StaticText(
            self.temperature_change_window, label="New Temperature (K):"
        )
        self.temperature_input = wx.TextCtrl(
            self.temperature_change_window,
            value=str(self.params.temperature),
            size=(100, 20),
        )
        self.temperature_sizer.AddSpacer(20)
        self.temperature_sizer.Add(self.temperature_label)
        self.temperature_sizer.AddSpacer(5)
        self.temperature_sizer.Add(self.temperature_input)
        self.temperature_change_window_sizer.AddSpacer(10)
        self.temperature_change_window_sizer.Add(self.temperature_sizer)
        self.temperature_change_window_sizer.AddSpacer(20)
        self.save_and_exit_temperature_change = wx.Button(
            self.temperature_change_window, label="Save", size=(100, 20)
        )
        self.save_and_exit_temperature_change.Bind(
            wx.EVT_BUTTON, self.on_save_and_exit_temperature_change
        )
        self.save_button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.save_button_sizer.AddSpacer(100)
        self.save_button_sizer.Add(self.save_and_exit_temperature_change, 0, wx.CENTER)
        self.save_button_sizer.AddSpacer(100)
        self.temperature_change_window_sizer.Add(self.save_button_sizer)
        self.temperature_change_window_sizer.AddSpacer(20)
        self.temperature_change_window.SetSizer(self.temperature_change_window_sizer)
        self.temperature_change_window.Show()

    def on_save_and_exit_temperature_change(self, event):
        self.params.temperature = float(self.temperature_input.GetValue())
        if self.nmrdata.spectrometer == "Bruker":
            self.params.calculate_carrier_frequency_bruker()
        else:
            self.params.calculate_carrier_frequency_varian()

        self.temperature_parameter.SetLabel(str(self.params.temperature))
        self.temperature_change_window.Close()

        if self.nmrdata.spectrometer == "Bruker":
            # Get all the indexes of the current direct/indirect dimension carrier frequency boxes
            j = 0
            for i in range(len(self.app.format.carrier_combo_boxes)):
                if self.app.format.carrier_combo_boxes[i] == 200:
                    j += 1
                else:
                    label = self.app.format.carrier_combo_boxes[i].GetValue()
                    index = 0

                    if i == 0:
                        for k in range(len(self.params.references_proton_labels)):
                            if label == self.params.references_proton_labels[k]:
                                index = k
                                break
                        self.app.format.carrier_frequency_boxes[i].SetValue(
                            str(self.params.references_proton[index])
                        )

                    else:
                        for k in range(len(self.params.references_other_labels)):
                            if label == self.params.references_other_labels[k]:
                                index = k
                                break
                        self.app.format.carrier_frequency_boxes[i].SetValue(
                            str(self.params.references_other[index])
                        )
        else:
            # Get all the indexes of the current direct/indirect dimension carrier frequency boxes
            j = 0
            for i in range(len(self.app.format.carrier_combo_boxes)):
                if (
                    self.app.format.carrier_combo_boxes[i].GetValue() == "ID"
                    or self.app.format.carrier_combo_boxes[i].GetValue() == "Other"
                    or self.app.format.carrier_combo_boxes[i].GetValue() == "Manual"
                    or self.app.format.acqusition_combo_boxes[i].GetValue() == "Real"
                ):
                    j += 1
                else:
                    label = self.app.format.carrier_combo_boxes[i].GetValue()
                    index = 0

                    if i == 0:
                        for k in range(len(self.params.references_proton_labels)):
                            if label == self.params.references_proton_labels[k]:
                                index = k
                                break
                        self.app.format.carrier_frequency_boxes[i].SetValue(
                            str(self.params.references_proton[index])
                        )

                    else:
                        for k in range(len(self.params.references_other_labels)):
                            if label == self.params.references_other_labels[k]:
                                index = k
                                break
                        self.app.format.carrier_frequency_boxes[i].SetValue(
                            str(self.params.references_other[index])
                        )

    def create_intensity_scaling_box(self):
        self.scaling_box = wx.StaticBox(self.app, -1, label="Intensity Scaling")
        self.scaling_box_sizer_total = wx.StaticBoxSizer(self.scaling_box, wx.VERTICAL)
        self.scaling_box_sizer = wx.BoxSizer(wx.HORIZONTAL)
        if self.nmrdata.spectrometer == "Bruker":
            self.scaling_box_sizer_total.AddSpacer(2)
            self.scaling_box_sizer_total.Add(self.scaling_box_sizer)
            self.scaling_box_sizer_total.AddSpacer(2)
        else:
            self.scaling_box_sizer_total.AddSpacer(13)
            self.scaling_box_sizer_total.Add(self.scaling_box_sizer)
            self.scaling_box_sizer_total.AddSpacer(13)
        # Create tick boxes for intensity scaling
        self.scaling_NS_checkbox = wx.CheckBox(
            self.app, -1, label="1/NS"
        )  # Normalise by number of scans
        if self.nmrdata.spectrometer == "Bruker":
            self.scaling_NC = wx.CheckBox(
                self.app, -1, label="2^NC"
            )  # Normalise by bruker normalisation constant
        self.scaling_by_number = wx.CheckBox(self.app, -1, label="x1000")
        self.scaling_NS_checkbox.Bind(wx.EVT_CHECKBOX, self.on_scaling_checkbox)
        if self.nmrdata.spectrometer == "Bruker":
            self.scaling_NC.Bind(wx.EVT_CHECKBOX, self.on_scaling_checkbox)
        self.scaling_by_number.Bind(wx.EVT_CHECKBOX, self.on_scaling_checkbox)

        if self.params.include_scaling == True:
            self.scaling_NS_checkbox.SetValue(True)
            if self.nmrdata.spectrometer == "Bruker":
                self.scaling_NC.SetValue(True)
        self.scaling_by_number.SetValue(True)
        if self.nmrdata.spectrometer == "Bruker":
            self.params.scaling_factor = (
                (1 / self.params.NS) * (2**self.params.NC) * 1000
            )
        else:
            self.params.scaling_factor = (1 / self.params.NS) * 1000
        self.scaling_text = wx.StaticText(self.app, label="Scaling Factor:")
        if self.nmrdata.spectrometer == "Bruker":
            self.scaling_number = wx.TextCtrl(
                self.app,
                value="{:.2E}".format(self.params.scaling_factor),
                size=(50, 20),
            )
        else:
            self.scaling_number = wx.TextCtrl(
                self.app,
                value="{:.2E}".format(self.params.scaling_factor),
                size=(80, 20),
            )
        self.scaling_box_sizer.Add(self.scaling_NS_checkbox)
        self.scaling_box_sizer.AddSpacer(20)
        if self.nmrdata.spectrometer == "Bruker":
            self.scaling_box_sizer.Add(self.scaling_NC)
            self.scaling_box_sizer.AddSpacer(20)
        self.scaling_box_sizer.Add(self.scaling_by_number)
        self.scaling_box_sizer.AddSpacer(20)
        self.scaling_box_sizer.Add(self.scaling_text)
        self.scaling_box_sizer.AddSpacer(5)
        self.scaling_box_sizer.Add(self.scaling_number)

        if self.nmrdata.spectrometer == "Bruker":
            self.app.bottom_left_box.AddSpacer(10)
            self.app.bottom_left_box.Add(self.scaling_box_sizer_total)
        else:
            self.app.extra_boxes_0.AddSpacer(20)
            self.app.extra_boxes_0.Add(self.scaling_box_sizer_total)

    def on_scaling_checkbox(self, event):

        value = 1
        if self.scaling_NS_checkbox.GetValue() == True:
            value = value * (1 / self.params.NS)
        if self.nmrdata.spectrometer == "Bruker":
            if self.scaling_NC.GetValue() == True:
                value = value * (2**self.params.NC)
        if self.scaling_by_number.GetValue() == True:
            value = value * 1000

        self.params.scaling_factor = value
        self.scaling_number.SetValue("{:.2E}".format(self.params.scaling_factor))

    def on_acquisition_mode_change(self, event) -> None:
        """
        Perform an action if the acquisition mode of a dimension is changed
        by the user to real. # If combobox is changed to real, then change
        the number of real points to the same as the number of complex points,
        also set label to ID, sw to 0.
        """

        # Find out which combobox has been changed
        for i in range(len(self.app.format.acqusition_combo_boxes)):
            if self.app.format.acqusition_combo_boxes[i] == event.GetEventObject():
                index = i
                break
        if self.app.format.acqusition_combo_boxes[index].GetValue() == "Real":
            if index == 0:
                self.app.format.N_real_boxes[index].SetValue(
                    str(self.params.size_direct)
                )
                # Set label to ID
                self.app.format.nucleus_type_boxes[index].SetValue("ID")
                # Set sw to 0
                self.app.format.sweep_width_boxes[index].SetValue(str(1))
                # Set nucleus frequency to 0
                self.app.format.nuclei_frequency_boxes[index].SetValue(str(1))
                # Set carrier frequency to 0
                self.app.format.carrier_frequency_boxes[index].SetValue(str(1))
                # Set combobox to N/A
                self.app.format.carrier_combo_boxes[index].SetValue("N/A")

            else:
                self.app.format.N_real_boxes[index].SetValue(
                    self.app.format.N_complex_boxes[index].GetValue()
                )
                # Set label to ID
                self.app.format.nucleus_type_boxes[index].SetValue("ID")
                # Set sw to 0
                self.app.format.sweep_width_boxes[index].SetValue(str(1))
                # Set nucleus frequency to 0
                self.app.format.nuclei_frequency_boxes[index].SetValue(str(1))
                # Set carrier frequency to 0
                self.app.format.carrier_frequency_boxes[index].SetValue(str(1))
                # Set combobox to N/A
                self.app.format.carrier_combo_boxes[index].SetValue("N/A")
        else:
            if index == 0:
                self.app.format.N_real_boxes[index].SetValue(
                    str(int(self.params.size_direct / 2))
                )
                # Set label back to original label
                self.app.format.nucleus_type_boxes[index].SetValue(
                    self.params.labels_correct_order[index]
                )
                # Set sw back to original sw
                self.app.format.sweep_width_boxes[index].SetValue(
                    str(self.params.sw_direct)
                )
                # Set nucleus frequency back to original nucleus frequency
                self.app.format.nuclei_frequency_boxes[index].SetValue(
                    str(self.params.nucleus_frequencies[index])
                )
                # Set carrier frequency back to original carrier frequency
                self.app.format.carrier_frequency_boxes[index].SetValue(
                    str(self.params.references_proton[index])
                )
                # Set combobox back to original combobox
                self.app.format.carrier_combo_boxes[index].SetValue(
                    self.params.references_proton_labels[index]
                )
            else:
                if self.nmrdata.spectrometer == "Bruker":
                    self.app.format.N_real_boxes[index].SetValue(
                        str(
                            int(
                                int(self.app.format.N_complex_boxes[index].GetValue())
                                / 2
                            )
                        )
                    )
                else:
                    self.app.format.N_real_boxes[index].SetValue(
                        str(int(self.params.size_indirect[index - 1]))
                    )
                # Set label back to original label
                if self.nmrdata.spectrometer == "Bruker":
                    self.app.format.nucleus_type_boxes[index].SetValue(
                        self.params.labels_correct_order[index]
                    )
                else:
                    self.app.format.nucleus_type_boxes[index].SetValue(
                        self.params.labels_correct_order[index - 1]
                    )
                # Set sw back to original sw
                if self.nmrdata.spectrometer == "Bruker":
                    self.app.format.sweep_width_boxes[index].SetValue(
                        str(self.params.sw_indirect[index - 1])
                    )
                else:
                    self.app.format.sweep_width_boxes[index].SetValue(
                        str(self.params.sw_indirect["sw" + str(index)])
                    )
                # Set nucleus frequency back to original nucleus frequency
                self.app.format.nuclei_frequency_boxes[index].SetValue(
                    str(self.params.nucleus_frequencies[index])
                )

                # Set carrier frequency back to original carrier frequency
                self.app.format.carrier_frequency_boxes[index].SetValue(
                    str(self.params.references_other[index - 1])
                )
                # Set combobox back to original combobox
                self.app.format.carrier_combo_boxes[index].SetValue(
                    self.params.references_other_labels[index - 1]
                )

    def acquisition_2D_mode_combo_box(self) -> None:
        """
        A combobox for the 2D acquisition mode
        """
        self.acquisition_2D_mode_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.acquisition_2D_mode_txt = wx.StaticText(
            self.app, label="Acquisition Mode 2D:       "
        )
        self.acquisition_mode_2D_choices = [
            "States",
            "TPPI",
            "Magnitude",
            "Real",
            "Complex",
            "Image",
        ]

        if self.nmrdata.params.pseudo_flag == 0:
            self.acquisition_2D_mode_box = wx.ComboBox(
                self.app,
                value=self.acquisition_mode_2D_choices[0],
                choices=self.acquisition_mode_2D_choices,
                size=(200, 20),
            )
        else:
            if len(self.params.size_indirect) == 1:
                self.acquisition_2D_mode_box = wx.ComboBox(
                    self.app,
                    value=self.acquisition_mode_2D_choices[3],
                    choices=self.acquisition_mode_2D_choices,
                    size=(200, 20),
                )
            else:
                self.acquisition_2D_mode_box = wx.ComboBox(
                    self.app,
                    value=self.acquisition_mode_2D_choices[0],
                    choices=self.acquisition_mode_2D_choices,
                    size=(200, 20),
                )

        self.acquisition_2D_mode_sizer.AddSpacer(20)
        self.acquisition_2D_mode_sizer.Add(self.acquisition_2D_mode_txt)
        self.acquisition_2D_mode_sizer.AddSpacer(11)
        self.acquisition_2D_mode_sizer.AddSpacer(220)
        self.acquisition_2D_mode_sizer.Add(self.acquisition_2D_mode_box)
        self.app.menu_bar.AddSpacer(10)
        self.app.menu_bar.Add(self.acquisition_2D_mode_sizer)

    def create_conversion_box(self):
        # Have a button for make fid.com
        if self.nmrdata.spectrometer == "Bruker":
            self.fid_conversion_box = wx.BoxSizer(wx.VERTICAL)
        else:
            self.fid_conversion_box = wx.BoxSizer(wx.HORIZONTAL)

        self.save_parameters_button = wx.Button(
            self.app, label="Save Parameters", size=(175, 20)
        )
        self.save_parameters_button.Bind(wx.EVT_BUTTON, self.app.on_save_parameters)
        self.fid_conversion_box.AddSpacer(20)
        self.fid_conversion_box.Add(self.save_parameters_button)

        # Have a button for convert
        self.convert_button = wx.Button(
            self.app, label="Convert Data (nmrPipe)", size=(175, 20)
        )
        self.convert_button.Bind(wx.EVT_BUTTON, self.app.on_convert_pipe)
        self.fid_conversion_box.AddSpacer(20)
        self.fid_conversion_box.Add(self.convert_button)

        # Have a button for convert
        self.convert_button2 = wx.Button(
            self.app, label="Convert Data (nmrglue)", size=(175, 20)
        )
        self.convert_button2.Bind(wx.EVT_BUTTON, self.app.on_convert_glue)
        self.fid_conversion_box.AddSpacer(20)
        self.fid_conversion_box.Add(self.convert_button2)

        self.app.extra_boxes.AddSpacer(20)
        self.app.extra_boxes.Add(self.fid_conversion_box)

    def find_nus_file(self):
        self.nuslist_found = False
        for file in os.listdir():
            if file == "nuslist":
                self.nuslist_found = True
                self.nusfile = file
                self.include_NUS = True
                break

        if self.nuslist_found == False:
            self.include_NUS = False
            self.nusfile = ""

    def input_NUS_list_box(self):
        self.NUS_label = wx.StaticBox(self.app, -1, label="NUS Information")
        self.NUS_box = wx.StaticBoxSizer(self.NUS_label, wx.VERTICAL)
        self.NUS_box_2 = wx.BoxSizer(wx.VERTICAL)

        # Tick box to see if user wants to perform NUS reconstruction
        self.NUS_tickbox = wx.CheckBox(self.app, -1, label="Include NUS Reconstruction")
        self.NUS_tickbox.Bind(wx.EVT_CHECKBOX, self.On_NUS_CheckBox)
        self.NUS_box.Add(self.NUS_tickbox)
        self.NUS_box_2 = wx.BoxSizer(wx.VERTICAL)

        if self.include_NUS == True:
            self.NUS_offset = 0
            self.NUS_sample_count = 0
            self.nusfile = ""
            for file in os.listdir():
                if file == "nuslist":
                    self.nuslist_found = True
                    self.nusfile = file
                    break
            if self.nuslist_found == True:
                lines = open("nuslist", "r").readlines()
                self.NUS_sample_count = len(lines)
                if lines[0].split()[0] == "0":
                    self.NUS_offset = 0
                else:
                    self.NUS_offset = 1

            self.NUS_tickbox.SetValue(True)

            self.nus_details_box()

        self.NUS_box.Add(self.NUS_box_2)

        self.app.extra_boxes.AddSpacer(20)
        self.app.extra_boxes.Add(self.NUS_box)

    def nus_details_box(self):
        self.NUS_sample_count_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.NUS_sample_count_txt = wx.StaticText(self.app, label="NUS Sample Count:")
        self.NUS_sample_count_box = wx.TextCtrl(
            self.app, value=str(self.NUS_sample_count), size=(50, 20)
        )
        self.NUS_sample_count_sizer.AddSpacer(5)
        self.NUS_sample_count_sizer.Add(self.NUS_sample_count_txt)
        self.NUS_sample_count_sizer.AddSpacer(5)
        self.NUS_sample_count_sizer.Add(self.NUS_sample_count_box)

        self.NUS_box_2.AddSpacer(10)
        self.NUS_box_2.Add(self.NUS_sample_count_sizer)

        self.nusfile_label = wx.StaticText(self.app, label="NUS schedule:")
        self.nusfile_input = wx.TextCtrl(self.app, value=self.nusfile, size=(120, 20))
        self.find_nus_file = wx.Button(self.app, label="...", size=(25, 20))
        self.find_nus_file.Bind(
            wx.EVT_BUTTON, lambda evt: self.on_find_nus_file(evt, self.nusfile)
        )
        self.nusfile_input.SetValue(self.nusfile)

        self.nus_extras_box = wx.BoxSizer(wx.HORIZONTAL)
        self.nus_offset_label = wx.StaticText(self.app, label="NUS Offset:")
        self.nus_offset_box = wx.TextCtrl(
            self.app, value=str(self.NUS_offset), size=(30, 20)
        )
        self.nus_extras_box.AddSpacer(5)
        self.nus_extras_box.Add(self.nus_offset_label)
        self.nus_extras_box.AddSpacer(5)
        self.nus_extras_box.Add(self.nus_offset_box)

        self.reverse_NUS_tickbox = wx.CheckBox(
            self.app, -1, label="Reverse NUS Schedule"
        )
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
        # get dialog box here
        cwd = os.getcwd()
        dlg = wx.FileDialog(
            self.app,
            message="Choose a file",
            defaultDir=cwd,
            defaultFile=self.nusfile,
            style=wx.FD_OPEN | wx.FD_CHANGE_DIR,
        )
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            if cwd in path:
                splitPath = path.split(cwd)
                textBox = "." + splitPath[1]
            else:
                textBox = path

            self.nusfile_input.SetValue(textBox)
            self.app.Layout()
            self.app.Refresh()

        dlg.Destroy()

        try:
            self.read_nus_file()
        except:
            # Give popout saying NUS file not read correctly
            dlg = wx.MessageDialog(
                self.app, "NUS file not read correctly", "Error", wx.OK | wx.ICON_ERROR
            )
            self.Raise()
            self.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            self.nusfile_input.SetValue("")
            self.nusfile = ""
            self.NUS_sample_count = 0
            self.NUS_sample_count_box.SetValue(str(self.NUS_sample_count))
            self.nus_offset_box.SetValue("0")
            self.NUS_offset = 0
            self.app.Layout()
            self.app.Refresh()

    def read_nus_file(self):
        # Read the nus file
        self.nusfile = self.nusfile_input.GetValue()
        lines = open(self.nusfile, "r").readlines()
        self.NUS_sample_count = len(lines)
        self.NUS_sample_count_box.SetValue(str(self.NUS_sample_count))
        if lines[0].split()[0] == "0":
            self.NUS_offset = 0
        else:
            self.NUS_offset = 1
        self.nus_offset_box.SetValue(str(self.NUS_offset))
        self.app.Layout()
        self.app.Refresh()

    def On_NUS_CheckBox(self, event):

        if self.NUS_tickbox.GetValue() == True:
            self.include_NUS = True
            self.NUS_box_2.Clear(True)
            self.NUS_box.Clear(True)
            self.app.extra_boxes.Remove(self.NUS_box)
            self.app.extra_boxes.Detach(len(self.app.extra_boxes.GetChildren()) - 1)
            self.input_NUS_list_box()
            self.app.Layout()
            self.app.Refresh()
        else:
            self.include_NUS = False
            self.NUS_box_2.Clear(True)
            self.NUS_box.Clear(True)
            self.app.extra_boxes.Remove(self.NUS_box)
            self.app.extra_boxes.Detach(len(self.app.extra_boxes.GetChildren()) - 1)
            self.input_NUS_list_box()
            self.app.Layout()
            self.app.Refresh()
