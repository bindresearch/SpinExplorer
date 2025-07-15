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


class FormatParametersBruker:
    def __init__(self, app, params, nmrdata) -> None:
        """
        This class contains the functions to populate the GUI with the
        determined parameters (stored in params)
        """
        self.app = app
        self.params = params
        self.nmrdata = nmrdata

    def input_sizes_bruker(self) -> None:
        """
        Creating and inputting boxes for the complex and real dimension sizes
        """

        # Create the boxes for the complex point numbers
        self.title_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.title_sizer.AddSpacer(250)
        self.N_complex_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.N_complex_txt = wx.StaticText(self.app, label="Number complex points:")
        self.N_complex_boxes = []
        for i in range(len(self.params.size_indirect) + 1):
            self.title_sizer.Add(
                wx.StaticText(self.app, label="Dimension " + str(i + 1))
            )
            self.title_sizer.AddSpacer(145)
            if i == 0:
                self.N_complex_boxes.append(
                    wx.TextCtrl(
                        self.app,
                        value=str(self.params.size_direct_complex),
                        size=(200, 20),
                    )
                )
            else:
                try:
                    if self.params.labels_correct_order[i] == "ID":
                        labelval = "off"
                    else:
                        labelval = self.params.labels_correct_order[i]
                    size = str(self.params.indirect_sizes_dict[labelval])
                except:
                    size = str(self.params.size_indirect[i - 1])
                self.N_complex_boxes.append(
                    wx.TextCtrl(self.app, value=size, size=(200, 20))
                )
        self.N_complex_sizer.AddSpacer(20)
        self.N_complex_sizer.Add(self.N_complex_txt)
        self.N_complex_sizer.AddSpacer(20)
        for i in range(0, len(self.N_complex_boxes)):
            self.N_complex_sizer.Add(self.N_complex_boxes[i])
            self.N_complex_sizer.AddSpacer(20)
        self.app.menu_bar.AddSpacer(10)
        self.app.menu_bar.Add(self.title_sizer)
        self.app.menu_bar.AddSpacer(10)
        self.app.menu_bar.Add(self.N_complex_sizer)

        # Create the boxes for the real point numbers
        self.N_real_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.N_real_txt = wx.StaticText(self.app, label="Number real points:        ")
        self.N_real_boxes = []

        if self.params.pseudo_flag == 0:
            for i in range(len(self.params.size_indirect) + 1):
                if i == 0:
                    self.N_real_boxes.append(
                        wx.TextCtrl(
                            self.app,
                            value=str(int(self.params.size_direct / 2)),
                            size=(200, 20),
                        )
                    )
                else:

                    try:
                        if self.params.labels_correct_order[i] == "ID":
                            labelval = "off"
                            size = str(int(self.params.indirect_sizes_dict[labelval]))
                        else:
                            labelval = self.params.labels_correct_order[i]
                            size = str(
                                int(self.params.indirect_sizes_dict[labelval] / 2)
                            )

                    except:
                        size = str(int(self.params.size_indirect[i - 1] / 2))
                    self.N_real_boxes.append(
                        wx.TextCtrl(self.app, value=size, size=(200, 20))
                    )

        else:
            for i in range(self.nmrdata.data_dimensions):
                if i == 0:
                    self.N_real_boxes.append(
                        wx.TextCtrl(
                            self.app, value=str(self.params.size_direct), size=(200, 20)
                        )
                    )
                else:
                    if i in self.params.pseudo_flag:
                        try:
                            if self.params.labels_correct_order[i] == "ID":
                                labelval = "off"
                            else:
                                labelval = self.params.labels_correct_order[i]
                            size = str(self.params.indirect_sizes_dict[labelval])
                        except:
                            size = str(self.params.size_indirect[i - 1])
                        self.N_real_boxes.append(
                            wx.TextCtrl(self.app, value=size), size=(200, 20)
                        )
                    else:
                        try:
                            if self.params.labels_correct_order[i] == "ID":
                                labelval = "off"
                            else:
                                labelval = self.params.labels_correct_order[i]
                            size = str(self.params.indirect_sizes_dict[labelval])
                        except:
                            size = str(int(self.params.size_indirect[i - 1] / 2))
                        self.N_real_boxes.append(
                            wx.TextCtrl(
                                self.app,
                                value=str(
                                    int(
                                        self.nmrdata.nmr_data.shape[
                                            self.nmrdata.data_dimensions - 1 - i
                                        ]
                                        / 2
                                    )
                                ),
                                size=(200, 20),
                            )
                        )

        self.N_real_sizer.AddSpacer(20)
        self.N_real_sizer.Add(self.N_real_txt)
        self.N_real_sizer.AddSpacer(20)
        for i in range(0, len(self.N_real_boxes)):
            self.N_real_sizer.Add(self.N_real_boxes[i])
            self.N_real_sizer.AddSpacer(20)
        self.app.menu_bar.AddSpacer(10)
        self.app.menu_bar.Add(self.N_real_sizer)

    def input_acquisition_modes_bruker(self) -> None:
        """
        Creating a drop down combobox with all the standard acquisition modes
        for Bruker direct and indirect dimensions. Initial selection is based
        on the acqusition modes found by ParameterExtractorBruker class
        (params)
        """
        self.acquisition_mode_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.acquisition_mode_txt = wx.StaticText(
            self.app, label="Acquisition mode:           "
        )
        self.acquisition_mode_options_direct = ["DQD", "Complex", "Sequential", "Real"]
        self.acquisition_mode_options_indirect = [
            "Complex",
            "States-TPPI",
            "Echo-AntiEcho",
            "TPPI",
            "States",
            "Real",
        ]

        self.acqusition_combo_boxes = []
        for i in range(len(self.params.size_indirect) + 1):
            if i == 0:
                self.acqusition_combo_boxes.append(
                    wx.ComboBox(
                        self.app,
                        value=self.acquisition_mode_options_direct[0],
                        choices=self.acquisition_mode_options_direct,
                        size=(200, 20),
                        style=wx.CB_READONLY,
                    )
                )
                self.acqusition_combo_boxes[i].Bind(
                    wx.EVT_COMBOBOX, self.app.shared_format.on_acquisition_mode_change
                )
            else:
                if self.params.acqusition_modes[i - 1] == "QF":
                    self.acqusition_combo_boxes.append(
                        wx.ComboBox(
                            self.app,
                            value=self.acquisition_mode_options_indirect[5],
                            choices=self.acquisition_mode_options_indirect,
                            size=(200, 20),
                            style=wx.CB_READONLY,
                        )
                    )
                    self.N_real_boxes[i].SetValue(
                        str(self.N_complex_boxes[i].GetValue())
                    )
                    self.params.sw_indirect[i - 1] = 0
                    self.acqusition_combo_boxes[i].Bind(
                        wx.EVT_COMBOBOX, self.shared_format.on_acquisition_mode_change
                    )
                else:
                    determined_value = self.params.acqusition_modes[i - 1]
                    detected_value = ""
                    for j in range(len(self.acquisition_mode_options_indirect)):
                        if (
                            self.acquisition_mode_options_indirect[j].upper()
                            == determined_value.upper()
                        ):
                            detected_value = self.acquisition_mode_options_indirect[j]
                            break
                    if detected_value == "":
                        self.acqusition_combo_boxes.append(
                            wx.ComboBox(
                                self.app,
                                value=self.acquisition_mode_options_indirect[0],
                                choices=self.acquisition_mode_options_indirect,
                                size=(200, 20),
                                style=wx.CB_READONLY,
                            )
                        )
                        self.acqusition_combo_boxes[i].Bind(
                            wx.EVT_COMBOBOX,
                            self.app.shared_format.on_acquisition_mode_change,
                        )
                    else:
                        self.acqusition_combo_boxes.append(
                            wx.ComboBox(
                                self.app,
                                value=detected_value,
                                choices=self.acquisition_mode_options_indirect,
                                size=(200, 20),
                                style=wx.CB_READONLY,
                            )
                        )
                        self.acqusition_combo_boxes[i].Bind(
                            wx.EVT_COMBOBOX,
                            self.app.shared_format.on_acquisition_mode_change,
                        )

        self.acquisition_mode_sizer.AddSpacer(20)
        self.acquisition_mode_sizer.Add(self.acquisition_mode_txt)
        self.acquisition_mode_sizer.AddSpacer(20)
        for i in range(0, len(self.acqusition_combo_boxes)):
            self.acquisition_mode_sizer.Add(self.acqusition_combo_boxes[i])
            self.acquisition_mode_sizer.AddSpacer(20)
        self.app.menu_bar.AddSpacer(10)
        self.app.menu_bar.Add(self.acquisition_mode_sizer)

    def input_sweep_widths_bruker(self) -> None:
        """
        Creating textcontrol boxes for the sweep widths in each dimension.
        Will set the values to those found by ParameterExtractorBruker
        class (params)
        """
        self.sweep_width_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.sweep_width_txt = wx.StaticText(
            self.app, label="Sweep width (Hz):          "
        )
        self.sweep_width_boxes = []
        choices = [str(self.params.sw_direct)]
        for j in range(len(self.params.size_indirect)):
            choices.append(str(self.params.sw_indirect[j]))
        for i in range(len(self.params.size_indirect) + 1):
            if i == 0:
                self.sweep_width_boxes.append(
                    wx.ComboBox(
                        self.app,
                        value=str(self.params.sw_direct),
                        choices=choices,
                        size=(200, 20),
                    )
                )
            else:
                self.sweep_width_boxes.append(
                    wx.ComboBox(
                        self.app,
                        value=str(self.params.sw_indirect[i - 1]),
                        choices=choices,
                        size=(200, 20),
                    )
                )
        self.sweep_width_sizer.AddSpacer(20)
        self.sweep_width_sizer.Add(self.sweep_width_txt)
        self.sweep_width_sizer.AddSpacer(20)
        for i in range(0, len(self.sweep_width_boxes)):
            self.sweep_width_sizer.Add(self.sweep_width_boxes[i])
            self.sweep_width_sizer.AddSpacer(20)
        self.app.menu_bar.AddSpacer(10)
        self.app.menu_bar.Add(self.sweep_width_sizer)

    def get_nuclei_frequency_bruker(self) -> None:
        """
        Creating a text box to house the spectrometer frequencies in MHz.
        Initial values are set to those found in the ParameterExtractorBruker
        class (params).
        """
        self.nuclei_frequency_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.nuclei_frequency_txt = wx.StaticText(
            self.app, label="Nuclei frequency (MHz):"
        )
        self.nuclei_frequency_boxes = []
        for i in range(len(self.params.size_indirect) + 1):
            if i == 0:
                self.nuclei_frequency_boxes.append(
                    wx.TextCtrl(
                        self.app,
                        value=str(self.params.nucleus_frequencies[0]),
                        size=(200, 20),
                    )
                )
            elif self.acqusition_combo_boxes[i].GetValue() == "Real":
                self.nuclei_frequency_boxes.append(
                    wx.TextCtrl(self.app, value=str(1), size=(200, 20))
                )
                self.params.labels_correct_order[i] = "ID"
            else:
                self.nuclei_frequency_boxes.append(
                    wx.TextCtrl(
                        self.app,
                        value=str(self.params.nucleus_frequencies[i]),
                        size=(200, 20),
                    )
                )

        self.nuclei_frequency_sizer.AddSpacer(20)
        self.nuclei_frequency_sizer.Add(self.nuclei_frequency_txt)
        self.nuclei_frequency_sizer.AddSpacer(19)
        for i in range(0, len(self.nuclei_frequency_boxes)):
            self.nuclei_frequency_sizer.Add(self.nuclei_frequency_boxes[i])
            self.nuclei_frequency_sizer.AddSpacer(20)
        self.app.menu_bar.AddSpacer(10)
        self.app.menu_bar.Add(self.nuclei_frequency_sizer)

    def get_nuclei_labels_bruker(self) -> None:
        """
        Text controls for spectrum dimension labels. By default these will be
        set to those found for each dimension in the ParamaterExtractorBruker
        class (params)
        """
        # Add a section for axis labels
        self.nucleus_type_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.nucleus_type_txt = wx.StaticText(
            self.app, label="Label:                              "
        )
        self.nucleus_type_boxes = []
        j = 0
        for i in range(len(self.params.size_indirect) + 1):
            if i == 0:
                self.nucleus_type_boxes.append(
                    wx.TextCtrl(
                        self.app,
                        value=self.params.labels_correct_order[0],
                        size=(200, 20),
                    )
                )
            elif self.acqusition_combo_boxes[i].GetValue() == "Real":
                self.nucleus_type_boxes.append(
                    wx.TextCtrl(self.app, value="ID", size=(200, 20))
                )
                j += 1
            else:
                self.nucleus_type_boxes.append(
                    wx.TextCtrl(
                        self.app,
                        value=self.params.labels_correct_order[i],
                        size=(200, 20),
                    )
                )

        self.nucleus_type_sizer.AddSpacer(20)
        self.nucleus_type_sizer.Add(self.nucleus_type_txt)
        self.nucleus_type_sizer.AddSpacer(22)
        for i in range(0, len(self.nucleus_type_boxes)):
            self.nucleus_type_sizer.Add(self.nucleus_type_boxes[i])
            self.nucleus_type_sizer.AddSpacer(20)

        self.app.menu_bar.AddSpacer(10)
        self.app.menu_bar.Add(self.nucleus_type_sizer)

    def get_carrier_frequencies_bruker(self) -> None:
        """
        TextControls for the carrier frequencies for each dimension.
        By default these will be set to those found using the
        ParameterExtractorBruker class (params).
        """
        self.carrier_frequency_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.carrier_frequency_txt = wx.StaticText(
            self.app, label="Carrier frequency (ppm):"
        )
        self.carrier_frequency_boxes = []
        self.carrier_combo_boxes = []
        if self.params.size_indirect == []:
            self.options_other_dimensions = []
        else:
            self.options_other_dimensions = self.params.references_other_labels

        id_columns = 0
        for i in range(len(self.params.size_indirect) + 1):
            if i == 0:
                if (
                    self.params.labels_correct_order[i] == "1H"
                    or self.params.labels_correct_order[i] == "H1"
                    or self.params.labels_correct_order[i] == "H"
                ):
                    self.carrier_frequency_boxes.append(
                        wx.TextCtrl(
                            self.app, value=str(self.params.water_ppm), size=(200, 20)
                        )
                    )
                    self.carrier_combo_boxes.append(
                        wx.ComboBox(
                            self.app,
                            value=self.params.references_proton_labels[0],
                            choices=self.params.references_proton_labels,
                            size=(200, 20),
                            style=wx.CB_READONLY,
                        )
                    )
                else:
                    self.carrier_frequency_boxes.append(
                        wx.TextCtrl(
                            self.app,
                            value=str(self.params.references_proton[i]),
                            size=(200, 20),
                        )
                    )
                    self.carrier_combo_boxes.append(
                        wx.ComboBox(
                            self.app,
                            value=self.params.references_proton_labels[0],
                            choices=self.params.references_proton_labels,
                            size=(200, 20),
                            style=wx.CB_READONLY,
                        )
                    )
            elif self.params.labels_correct_order[i] == "ID":
                self.carrier_frequency_boxes.append(
                    wx.TextCtrl(self.app, value="0.00", size=(200, 20))
                )
                self.carrier_combo_boxes.append(
                    wx.ComboBox(
                        self.app,
                        value="N/A",
                        choices=self.options_other_dimensions,
                        size=(200, 20),
                        style=wx.CB_READONLY,
                    )
                )
                id_columns += 1
            else:
                self.carrier_frequency_boxes.append(
                    wx.TextCtrl(
                        self.app,
                        value=str(self.params.references_other[i - 1]),
                        size=(200, 20),
                    )
                )
                self.carrier_combo_boxes.append(
                    wx.ComboBox(
                        self.app,
                        value=self.params.references_other_labels[i - 1],
                        choices=self.options_other_dimensions,
                        size=(200, 20),
                        style=wx.CB_READONLY,
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
                    wx.EVT_COMBOBOX, self.on_carrier_combo_box_change
                )

        self.app.menu_bar.AddSpacer(10)
        self.app.menu_bar.Add(self.carrier_combo_box_sizer)

    def on_carrier_combo_box_change(self, event) -> None:
        """
        If a combo box for the carrier is changed, such as from water
        to O1/BF1, then change the carrier in the textcontrol to this
        value.
        """
        # Find the index of the combo box that was changed
        index = self.carrier_combo_boxes.index(event.GetEventObject())
        # Find the value of the combo box
        value = event.GetEventObject().GetValue()
        index_selection = self.carrier_combo_boxes[index].GetSelection()

        if (
            self.params.labels_correct_order[index] == "1H"
            or self.params.labels_correct_order[index] == "H1"
            or self.params.labels_correct_order[index] == "H"
        ):
            if value == "H2O":
                self.carrier_frequency_boxes[index].SetValue(
                    str(self.params.references_proton[0])
                )
            elif value == "O1/BF1":
                self.carrier_frequency_boxes[index].SetValue(
                    str(self.params.references_proton[1])
                )
            else:
                self.carrier_frequency_boxes[index].SetValue(
                    str(self.params.references_other[index_selection])
                )
        else:
            if value == "O1/BF1":
                self.carrier_frequency_boxes[index].SetValue(
                    str(self.params.references_proton[0])
                )
            else:
                self.carrier_frequency_boxes[index].SetValue(
                    str(self.params.references_other[index_selection])
                )

    def create_bruker_digital_filter_box(self):
        """
        Creating a box to insert all of the Bruker digital filter information
        """
        self.digital_filter_box = wx.StaticBox(self.app, -1, label="Digital Filter")
        self.digital_filter_box_sizer_total = wx.StaticBoxSizer(
            self.digital_filter_box, wx.VERTICAL
        )
        self.digital_filter_box_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.digital_filter_box_sizer_total.AddSpacer(11)
        self.digital_filter_box_sizer_total.Add(self.digital_filter_box_sizer)
        self.digital_filter_box_sizer_total.AddSpacer(11)
        self.digital_filter_checkbox = wx.CheckBox(
            self.app, -1, label="Remove Digital Filter"
        )
        self.digital_filter_checkbox.Bind(
            wx.EVT_CHECKBOX, self.on_digital_filter_checkbox
        )
        self.digital_filter_box_sizer.Add(self.digital_filter_checkbox)
        self.digital_filter_box_sizer.AddSpacer(10)
        if self.params.include_digital_filter == True:
            self.digital_filter_checkbox.SetValue(True)
            self.create_bruker_digital_filter_options()
        else:
            self.digital_filter_box_sizer_total.AddSpacer(5)

        self.app.extra_boxes_0.Add(self.digital_filter_box_sizer_total)

    def create_bruker_digital_filter_options(self):
        """
        Creating boxes for all the Bruker digital filter options
        """
        self.decim_box = wx.BoxSizer(wx.HORIZONTAL)
        self.decim_text = wx.StaticText(self.app, label="Decimation Rate:")
        self.decim_value = str(self.params.decim)
        self.decim_textbox = wx.TextCtrl(self.app, value=str(self.decim_value))
        self.decim_box.AddSpacer(5)
        self.decim_box.Add(self.decim_text)
        self.decim_box.AddSpacer(5)
        self.decim_box.Add(self.decim_textbox)
        self.digital_filter_box_sizer.Add(self.decim_box)

        self.dspfvs_box = wx.BoxSizer(wx.HORIZONTAL)
        self.dspfvs_text = wx.StaticText(self.app, label="DSP Firmware Version:")
        self.dspfvs_value = str(self.params.dspfvs)
        self.dspfvs_textbox = wx.TextCtrl(self.app, value=str(self.dspfvs_value))
        self.dspfvs_box.AddSpacer(5)
        self.dspfvs_box.Add(self.dspfvs_text)
        self.dspfvs_box.AddSpacer(5)
        self.dspfvs_box.Add(self.dspfvs_textbox)
        self.digital_filter_box_sizer.AddSpacer(10)
        self.digital_filter_box_sizer.Add(self.dspfvs_box)

        self.grpdly_box = wx.BoxSizer(wx.HORIZONTAL)
        self.grpdly_text = wx.StaticText(self.app, label="Group Delay:")
        self.grpdly_value = str(self.params.grpdly)
        self.grpdly_textbox = wx.TextCtrl(self.app, value=str(self.grpdly_value))
        self.grpdly_box.AddSpacer(5)
        self.grpdly_box.Add(self.grpdly_text)
        self.grpdly_box.AddSpacer(5)
        self.grpdly_box.Add(self.grpdly_textbox)
        self.digital_filter_box_sizer.AddSpacer(10)
        self.digital_filter_box_sizer.Add(self.grpdly_box)

        # Have a radiobox to either remove the digital filter before processing or during processing
        self.digital_filter_radio_box = wx.RadioBox(
            self.app,
            choices=[
                "Remove before fourier transform",
                "Remove after fourier transform",
            ],
            majorDimension=1,
            style=wx.RA_SPECIFY_ROWS,
        )
        self.digital_filter_radio_box.SetSelection(1)
        self.digital_filter_box_sizer_total.AddSpacer(10)
        self.digital_filter_box_sizer_total.Add(self.digital_filter_radio_box)

    def on_digital_filter_checkbox(self, event):
        """
        Change the formatting of the application canvas when digital filter
        checkbox is toggled on or off.
        """
        if self.digital_filter_checkbox.GetValue() == True:
            self.params.include_digital_filter = True
        else:
            self.params.include_digital_filter = False

        self.digital_filter_box_sizer.Clear(True)
        self.digital_filter_box_sizer_total.Clear(True)
        self.app.extra_boxes_0.Remove(self.digital_filter_box_sizer_total)
        # self.extra_boxes_0.Detach(len(self.extra_boxes_0.GetChildren())-1)
        self.create_bruker_digital_filter_box()
        self.app.Layout()
        self.app.Refresh()

    def create_other_options_box(self):
        """
        Creating a box with other options, these are implemented only for
        nmrPipe based conversions.
        """
        self.app.bottom_left_box = wx.BoxSizer(wx.VERTICAL)
        self.app.other_options_box = wx.StaticBox(self.app, -1, label="Other Options")
        self.app.other_options_box_sizer_total = wx.StaticBoxSizer(
            self.app.other_options_box, wx.VERTICAL
        )
        self.app.other_options_box_sizer = wx.BoxSizer(wx.VERTICAL)
        self.app.other_options_box_sizer_total.AddSpacer(2)
        self.app.other_options_box_sizer_total.Add(self.app.other_options_box_sizer)
        self.app.other_options_box_sizer_total.AddSpacer(2)

        self.row_1 = wx.BoxSizer(wx.HORIZONTAL)

        # Bad point removal threshold
        self.bad_point_threshold_box = wx.BoxSizer(wx.HORIZONTAL)
        self.nmrdata.params.bad_point_threshold = 0.0
        self.bad_point_threshold_text = wx.StaticText(
            self.app, label="Bad Point Threshold:"
        )
        self.bad_point_threshold_value = str(self.nmrdata.params.bad_point_threshold)
        self.bad_point_threshold_textbox = wx.TextCtrl(
            self.app, value=str(self.bad_point_threshold_value)
        )
        self.bad_point_threshold_box.AddSpacer(5)
        self.bad_point_threshold_box.Add(self.bad_point_threshold_text)
        self.bad_point_threshold_box.AddSpacer(5)
        self.bad_point_threshold_box.Add(self.bad_point_threshold_textbox)
        self.row_1.Add(self.bad_point_threshold_box)

        # Remove acquisition padding
        self.remove_acquisition_padding_box = wx.BoxSizer(wx.HORIZONTAL)
        self.nmrdata.params.remove_acquisition_padding = True
        self.remove_acquisition_padding_checkbox = wx.CheckBox(
            self.app, label="Remove Acquisition Padding"
        )
        self.remove_acquisition_padding_checkbox.Bind(
            wx.EVT_CHECKBOX, self.on_remove_acquisition_padding_checkbox
        )
        self.remove_acquisition_padding_checkbox.SetValue(True)
        self.remove_acquisition_padding_box.AddSpacer(5)
        self.remove_acquisition_padding_box.Add(
            self.remove_acquisition_padding_checkbox
        )
        self.row_1.AddSpacer(10)
        self.row_1.Add(self.remove_acquisition_padding_box)

        self.app.other_options_box_sizer.Add(self.row_1)

        self.app.bottom_left_box.Add(self.app.other_options_box_sizer_total)
        self.app.extra_boxes.AddSpacer(10)
        self.app.extra_boxes.Add(self.app.bottom_left_box)

    def on_remove_acquisition_padding_checkbox(self, event):
        """
        When the checkbox is clicked, the remove_acquisition_padding flag
        is changed.
        """
        if self.remove_acquisition_padding_checkbox.GetValue() == True:
            self.nmrdata.params.remove_acquisition_padding = True
        else:
            self.nmrdata.params.remove_acquisition_padding = False
