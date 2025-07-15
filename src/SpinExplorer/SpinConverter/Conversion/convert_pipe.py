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
import os
import subprocess


class Convert_pipe:
    def __init__(self, app, params, nmrdata) -> None:
        """
        This class will perform the conversion of the NMR data
        to nmrPipe format using nmrPipe.
        """
        self.app = app
        self.params = params
        self.nmrdata = nmrdata

        duplicate_names = self.check_duplicate_names()
        if duplicate_names == True:
            return

        stopping = self.check_fid_com_file()
        if stopping == True:
            return

        # Checking to see if the nmrPipe command works (i.e. if nmrPipe is installed)
        conversion_script = WritePipe(self.app, self.params, self.nmrdata)

        # Running the conversion script
        self.run_conversion_script()

    def check_duplicate_names(self) -> bool:
        """
        Check to see if any of the labels have duplicate names, if so, give a
        warning
        """
        total_labels = []
        for i, label in enumerate(self.app.format.nucleus_type_boxes):
            total_labels.append(label.GetValue())
        if len(total_labels) != len(set(total_labels)):
            dlg = wx.MessageDialog(
                self.app,
                "Duplicate labels found. Please rename labels so each dimension has a different label, then try again.",
                "Warning",
                wx.OK | wx.ICON_WARNING,
            )
            self.app.Raise()
            self.app.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()
            return True

        return False

    def check_fid_com_file(self) -> bool:
        """
        Check to see if the fid.com file already exists
        """
        if self.app.file_parser == True:
            os.chdir(self.app.path)
        if "fid.com" in os.listdir():
            dlg = wx.MessageDialog(
                self.app,
                "fid.com already exists. Do you want to overwrite it?",
                "Warning",
                wx.YES_NO | wx.ICON_WARNING,
            )
            self.app.Raise()
            self.app.SetFocus()
            if dlg.ShowModal() == wx.ID_NO:
                dlg.Destroy()
                if self.app.file_parser == True:
                    os.chdir(self.app.cwd)
                return True
            dlg.Destroy()

        return False

    def testing_nmrpipe(self):
        """
        Testing to see if nmrPipe is actually installed and
        working correctly.
        """

    def run_conversion_script(self):
        """
        Trying to run NMRPipe conversion script
        """
        if "test.fid" in os.listdir():
            dlg = wx.MessageDialog(
                self.app,
                "The test.fid file already exists. Do you want to overwrite it?",
                "Warning",
                wx.YES_NO | wx.ICON_WARNING,
            )
            self.app.Raise()
            self.app.SetFocus()
            if dlg.ShowModal() == wx.ID_NO:
                dlg.Destroy()
                if self.app.file_parser == True:
                    os.chdir(self.cwd)
                return
            dlg.Destroy()

        print("converting nmrPipe")
        # Add the necessary permissions to the fid.com file
        os.system("chmod +x fid.com")
        # Run the fid.com file
        command = "csh fid.com"
        p = subprocess.Popen(command, shell=True)
        p.wait()


class WritePipe:
    def __init__(self, app, params, nmrdata) -> None:
        """
        This class reads the current SpinConverter parameters to
        write the relevant nmrpipe conversion script (fid.com)
        """
        self.app = app
        self.params = params
        self.nmrdata = nmrdata

        # Create the fid.com file
        fid_file = open("fid.com", "w")
        fid_file.write("#!/bin/csh\n\n")

        self.spectype = self.nmrdata.files[0]

        if len(self.app.format.N_complex_boxes) > 1:
            if self.app.shared_format.NUS_tickbox.GetValue() == True:
                fid_file = self.write_nus(fid_file)

        if self.nmrdata.spectrometer == "Bruker":
            fid_file = self.write_bruk2pipe_line(fid_file)

        else:
            fid_file = self.write_var2pipe_line(fid_file)

        if self.params.size_indirect == []:
            fid_file = self.write_spectral_parameters_1D(fid_file)
        elif len(self.params.size_indirect) == 1:
            fid_file = self.write_spectral_parameters_2D(fid_file)
        elif len(self.params.size_indirect) == 2:
            fid_file = self.write_spectral_parameters_3D(fid_file)

        fid_file.close()
        if self.app.file_parser == True:
            os.chdir(self.app.cwd)

    def write_bruk2pipe_line(self, fid_file):
        """
        Writing the bruk2pipe line detailing information about data byte
        size, big/little endian data and bruker digital filter information
        """
        if int(self.params.byte_order) == 1:
            self.byte_order_text = "-noaswap"
        else:
            self.byte_order_text = "-aswap"

        if self.params.d_type == "2":
            self.byte_size_text = " -ws 8 -noi2f "
        else:
            self.byte_size_text = ""

        if self.params.remove_acquisition_padding == True:
            if self.app.format.digital_filter_checkbox.GetValue() == False:
                fid_file.write(
                    "bruk2pipe -verb -in "
                    + self.spectype
                    + " \\\n   -bad "
                    + str(self.params.bad_point_threshold)
                    + " -ext "
                    + self.byte_order_text
                    + self.byte_size_text
                    + "     \\\n"
                )
            else:
                if self.app.format.digital_filter_radio_box.GetSelection() == 0:
                    self.AMX_vs_DMX = "-DMX"
                else:
                    self.AMX_vs_DMX = "-AMX"
                fid_file.write(
                    "bruk2pipe -verb -in "
                    + self.spectype
                    + " \\\n   -bad "
                    + str(self.params.bad_point_threshold)
                    + " -ext "
                    + self.byte_order_text
                    + " "
                    + self.AMX_vs_DMX
                    + " -decim "
                    + str(self.params.decim)
                    + " -dspfvs "
                    + str(self.params.dspfvs)
                    + " -grpdly "
                    + str(self.params.grpdly)
                    + self.byte_size_text
                    + "\\\n"
                )

        else:
            if self.app.format.digital_filter_checkbox.GetValue() == False:
                fid_file.write(
                    "bruk2pipe -verb -in "
                    + self.spectype
                    + " \\\n"
                    + "   "
                    + self.byte_order_text
                    + self.byte_size_text
                    + "\\\n"
                )
            else:
                if self.app.format.digital_filter_radio_box.GetSelection() == 0:
                    self.AMX_vs_DMX = "-DMX"
                else:
                    self.AMX_vs_DMX = "-AMX"
                fid_file.write(
                    "bruk2pipe -verb -in "
                    + self.spectype
                    + " \\\n"
                    + "   "
                    + self.byte_order_text
                    + " "
                    + self.AMX_vs_DMX
                    + " -decim "
                    + str(self.decim)
                    + " -dspfvs "
                    + str(self.dspfvs)
                    + " -grpdly "
                    + str(self.grpdly)
                    + self.byte_size_text
                    + " \\\n"
                )

        return fid_file

    def write_nus(self, fid_file):
        """
        Writing the script to reshuffle the data so it is in
        the correct order with missing values replaced with
        zero.
        """
        sample_count = str(self.app.shared_format.NUS_sample_count_box.GetValue())
        offset = str(self.app.shared_format.nus_offset_box.GetValue())
        nusfile = str(self.app.shared_format.nusfile_input.GetValue())
        if self.nmrdata.spectrometer == "Bruker":
            fid_file.write(
                "nusExpand.tcl -mode bruker -sampleCount "
                + str(sample_count)
                + " -off "
                + str(offset)
                + " \\\n"
                + " -in ./"
                + self.nmrdata.files[0]
                + " -out ./ser_full -sample "
                + nusfile
            )
            if self.app.shared_format.reverse_NUS_tickbox.GetValue() == True:
                fid_file.write(" -rev")
            fid_file.write("\n\n")
            self.spectype = "./ser_full"

        else:
            fid_file.write(
                "nusExpand.tcl -mode varian -sampleCount "
                + str(sample_count)
                + " -off "
                + str(offset)
                + " \\\n"
                + " -in ./"
                + self.nmrdata.files[0]
                + " -out ./fid_full -sample "
                + nusfile
            )
            if self.app.shared_format.reverse_NUS_tickbox.GetValue() == True:
                fid_file.write(" -rev")

            fid_file.write("\n\n")

            self.spectype = "./fid_full"

        return fid_file

    def write_var2pipe_line(self, fid_file):

        if self.params.reverse_acquisition_order == False:
            fid_file.write(
                "var2pipe -verb -in " + self.spectype + " \\\n" + "-noaswap" + "\\\n"
            )
        else:
            fid_file.write(
                "var2pipe -verb -in "
                + self.spectype
                + " \\\n"
                + "-noaswap -aqORD 1"
                + "\\\n"
            )

        return fid_file

    def write_spectral_parameters_1D(self, fid_file):
        """
        Writing all the number of complex/real points, frequencies
        etc for a 1D spectrum
        """

        file_array = []
        file_array.append(
            ["-xN", self.app.format.N_complex_boxes[0].GetValue().strip(), " \\\n"]
        )
        file_array.append(
            ["-xT", self.app.format.N_real_boxes[0].GetValue().strip(), "  \\\n"]
        )
        file_array.append(
            [
                "-xMODE",
                self.app.format.acqusition_combo_boxes[0].GetValue().strip(),
                " \\\n",
            ]
        )
        file_array.append(
            [
                "-xSW",
                self.app.format.sweep_width_boxes[0].GetValue().strip(),
                " \\\n",
            ]
        )
        file_array.append(
            [
                "-xOBS",
                self.app.format.nuclei_frequency_boxes[0].GetValue().strip(),
                " \\\n",
            ]
        )
        file_array.append(
            [
                "-xCAR",
                self.app.format.carrier_frequency_boxes[0].GetValue().strip(),
                "\\\n",
            ]
        )
        file_array.append(
            [
                "-xLAB",
                self.app.format.nucleus_type_boxes[0].GetValue().strip(),
                " \\\n",
            ]
        )
        file_array.append(["-ndim", "1", " \\\n"])
        for row in file_array:
            fid_file.write("{:>10} {:>20} {:>5}".format(*row))
        if float(self.app.shared_format.scaling_number.GetValue()) != 1:
            fid_file.write(
                "| nmrPipe -fn MULT -c "
                + str(self.app.shared_format.scaling_number.GetValue())
                + " \\\n"
            )
        fid_file.write(" -ov -out ./test.fid\n")

        return fid_file

    def write_spectral_parameters_2D(self, fid_file):
        """
        Writing the necessary spectral parameters for a 2D
        spectrum
        """
        file_array = []
        file_array.append(
            [
                "-xN",
                self.app.format.N_complex_boxes[0].GetValue().strip(),
                "-yN",
                self.app.format.N_complex_boxes[1].GetValue().strip(),
                "\\\n",
            ]
        )
        file_array.append(
            [
                "-xT",
                self.app.format.N_real_boxes[0].GetValue().strip(),
                "-yT",
                self.app.format.N_real_boxes[1].GetValue().strip(),
                "\\\n",
            ]
        )
        file_array.append(
            [
                "-xMODE",
                self.app.format.acqusition_combo_boxes[0].GetValue().strip(),
                "-yMODE",
                self.app.format.acqusition_combo_boxes[1].GetValue().strip(),
                "\\\n",
            ]
        )
        file_array.append(
            [
                "-xSW",
                self.app.format.sweep_width_boxes[0].GetValue().strip(),
                "-ySW",
                self.app.format.sweep_width_boxes[1].GetValue().strip(),
                "\\\n",
            ]
        )
        file_array.append(
            [
                "-xOBS",
                self.app.format.nuclei_frequency_boxes[0].GetValue().strip(),
                "-yOBS",
                self.app.format.nuclei_frequency_boxes[1].GetValue().strip(),
                "\\\n",
            ]
        )
        file_array.append(
            [
                "-xCAR",
                self.app.format.carrier_frequency_boxes[0].GetValue().strip(),
                "-yCAR",
                self.app.format.carrier_frequency_boxes[1].GetValue().strip(),
                "\\\n",
            ]
        )
        file_array.append(
            [
                "-xLAB",
                self.app.format.nucleus_type_boxes[0].GetValue().strip(),
                "-yLAB",
                self.app.format.nucleus_type_boxes[1].GetValue().strip(),
                "\\\n",
            ]
        )
        file_array.append(
            [
                "-ndim",
                "2",
                "-aq2D",
                self.app.shared_format.acquisition_2D_mode_box.GetValue().strip(),
                "\\\n",
            ]
        )
        for row in file_array:
            fid_file.write("{:>10} {:>20} {:>10} {:>20} {:>5}".format(*row))
        if float(self.app.shared_format.scaling_number.GetValue()) != 1:
            fid_file.write(
                "| nmrPipe -fn MULT -c "
                + str(self.app.shared_format.scaling_number.GetValue())
                + " \\\n"
            )
        fid_file.write(" -ov -out ./test.fid\n")

        return fid_file

    def write_spectral_parameters_3D(self, fid_file):
        """
        Writing the necessary spectral parameters for a 3D
        spectrum
        """
        file_array = []
        file_array.append(
            [
                "-xN",
                self.app.format.N_complex_boxes[0].GetValue().strip(),
                "-yN",
                self.app.format.N_complex_boxes[1].GetValue().strip(),
                "-zN",
                self.app.format.N_complex_boxes[2].GetValue().strip(),
                "\\\n",
            ]
        )
        file_array.append(
            [
                "-xT",
                self.app.format.N_real_boxes[0].GetValue().strip(),
                "-yT",
                self.app.format.N_real_boxes[1].GetValue().strip(),
                "-zT",
                self.app.format.N_real_boxes[2].GetValue().strip(),
                "\\\n",
            ]
        )
        file_array.append(
            [
                "-xMODE",
                self.app.format.acqusition_combo_boxes[0].GetValue().strip(),
                "-yMODE",
                self.app.format.acqusition_combo_boxes[1].GetValue().strip(),
                "-zMODE",
                self.app.format.acqusition_combo_boxes[2].GetValue().strip(),
                "\\\n",
            ]
        )
        file_array.append(
            [
                "-xSW",
                self.app.format.sweep_width_boxes[0].GetValue().strip(),
                "-ySW",
                self.app.format.sweep_width_boxes[1].GetValue().strip(),
                "-zSW",
                self.app.format.sweep_width_boxes[2].GetValue().strip(),
                "\\\n",
            ]
        )
        file_array.append(
            [
                "-xOBS",
                self.app.format.nuclei_frequency_boxes[0].GetValue().strip(),
                "-yOBS",
                self.app.format.nuclei_frequency_boxes[1].GetValue().strip(),
                "-zOBS",
                self.app.format.nuclei_frequency_boxes[2].GetValue().strip(),
                "\\\n",
            ]
        )
        file_array.append(
            [
                "-xCAR",
                self.app.format.carrier_frequency_boxes[0].GetValue().strip(),
                "-yCAR",
                self.app.format.carrier_frequency_boxes[1].GetValue().strip(),
                "-zCAR",
                self.app.format.carrier_frequency_boxes[2].GetValue().strip(),
                "\\\n",
            ]
        )
        file_array.append(
            [
                "-xLAB",
                self.app.format.nucleus_type_boxes[0].GetValue().strip(),
                "-yLAB",
                self.app.format.nucleus_type_boxes[1].GetValue().strip(),
                "-zLAB",
                self.app.format.nucleus_type_boxes[2].GetValue().strip(),
                "\\\n",
            ]
        )
        for row in file_array:
            fid_file.write(
                "{:>10} {:>20} {:>10} {:>20} {:>10} {:>20} {:>5}".format(*row)
            )
        row = [
            "-ndim",
            "3",
            "-aq2D",
            self.app.shared_format.acquisition_2D_mode_box.GetValue().strip(),
            "",
            "",
            "\\\n",
        ]
        fid_file.write("{:>10} {:>20} {:>10} {:>20} {:>10} {:>20} {:>5}".format(*row))
        if float(self.app.shared_format.scaling_number.GetValue()) != 1:
            fid_file.write(
                "| nmrPipe -fn MULT -c "
                + str(self.app.shared_format.scaling_number.GetValue())
                + " \\\n"
            )
        fid_file.write(" -ov -out ./test.fid\n")

        return fid_file
