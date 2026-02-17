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


import nmrglue as ng
import os
import warnings

warnings.simplefilter("ignore", UserWarning)
from typing import List


class ParameterExtractorBruker:
    def __init__(self, nmrdata) -> None:
        """
        This class contains functions to search through the bruker acqus files
        to find relevant parameters needed for conversion to nmrPipe format.
        """
        self.nmrdata = nmrdata

        self.acqus_file = open(self.nmrdata.parameter_file, "r")
        self.acqus_file_lines = self.acqus_file.readlines()
        self.acqus_file.close()
        self.udic = ng.bruker.guess_udic(self.nmrdata.nmr_dic, self.nmrdata.nmr_data)

    def find_size_bruker(self) -> None:
        """
        Search through Bruker acqus file to search for NMR spectrum sizes.
        The TD parameter denotes the direct dimension size and TD_INDIRECT
        for the indirect dimension sizes (at least as a first guess).
        Sometimes this varies depending on pulse sequence.
        """
        # Look for TD in the acqus file
        self.size_direct, self.size_direct_complex = self.find_direct_bruker()
        # Create an array of initial guess for the indirect dimension sizes
        self.find_indirect_bruker()
        # Checking that the direct dimension TD value agrees with the data size
        self.checking_direct_bruker()
        # # Checking the indirect dimension
        # self.checking_indirect_bruker()

    def find_direct_bruker(self) -> List[int]:
        """
        Reading through Bruker acqus file lines to find the TD entry for
        direct dimension size.
        """
        for i in range(len(self.acqus_file_lines)):
            if "##$TD=" in self.acqus_file_lines[i]:
                line = self.acqus_file_lines[i].split()
                size_direct = int(line[1])
                size_direct_complex = int(size_direct)
                break

        return size_direct, size_direct_complex

    def checking_direct_bruker(self) -> None:
        try:
            self.size_direct
            # Find the max value in nmr_data.shape
            if len(self.nmrdata.nmr_data.shape) == 1:
                self.size_1 = max(self.nmrdata.nmr_data.shape)
                self.size_2 = self.size_1
                for i in range(len(self.size_indirect_non_reduced)):
                    self.size_2 = self.size_2 / self.size_indirect_non_reduced[i]
                if self.size_direct_complex < int(self.size_2 * 2):
                    self.size_direct_complex = int(self.size_2 * 2)
                if(self.size_direct_complex > 1.5*self.size_direct):
                    self.size_direct_complex = self.size_direct
            else:
                # Sometimes have the issue where stored complex data size is larger than TD, this ensures that the direct dimension
                # size is altered to the larger size to correctly be read
                self.size_1 = self.nmrdata.nmr_data.shape[-1]
                self.size_2 = self.size_1
                if self.size_1 * 2 > self.size_direct:
                    self.size_direct_complex = int(self.size_1 * 2)
                self.size_direct = self.size_1 * 2

        except:
            print("error in finding dimensions sizes")


    def find_indirect_bruker(self) -> None:
        """
        If have not found TD_INDIRECT values, check through acqu2s and acqu3s
        so see if there are any values present. If no other values can be
        found, the spectrum will be assumed to be 1D.
        """
        # Look to see if there is an acqu2s/acqu3s file
        self.size_indirect = []
        self.size_indirect_non_reduced = (
            []
        )  # values without substracting -1 in the case of odd numbers
        # Checking other spectrometer files if TD_INDIRECT is empty
        if "acqu2s" in os.listdir("./"):
            self.sizes_dim2 = []
            self.sizes_dim2_nus = []
            # Look for TD in the acqu2s file
            with open("acqu2s", "r") as file:
                file_lines = file.readlines()
                for i in range(len(file_lines)):
                    if "##$NUSTD=" in file_lines[i] or "##$NusTD=" in file_lines[i]:
                        line = file_lines[i].split("\n")[0].split()
                        self.sizes_dim2_nus = int(line[1])
                    if "##$TD=" in file_lines[i]:
                        line = file_lines[i].split()
                        self.sizes_dim2 = int(line[1])
            if (
                self.sizes_dim2 != []
                and self.sizes_dim2 < self.sizes_dim2_nus
                and self.sizes_dim2 != 1
            ):
                if self.sizes_dim2 % 2 != 0:
                    self.sizes_dim2 -= 1
                    self.size_indirect_non_reduced.append(self.sizes_dim2 + 1)
                self.size_indirect.append(self.sizes_dim2)
            else:
                self.size_indirect.append(self.sizes_dim2_nus)
                self.size_indirect_non_reduced.append(self.sizes_dim2_nus)

        if "acqu3s" in os.listdir("./"):
            self.sizes_dim3 = []
            self.sizes_dim3_nus = []
            # Look for TD in the acqu3s file
            with open("acqu3s", "r") as file:
                file_lines = file.readlines()
                for i in range(len(file_lines)):
                    if "##$NUSTD=" in file_lines[i] or "##$NusTD=" in file_lines[i]:
                        line = file_lines[i].split("\n")[0].split()
                        self.sizes_dim3_nus = int(line[1])
                    if "##$TD=" in file_lines[i]:
                        line = file_lines[i].split()
                        self.sizes_dim3 = int(line[1])
            if (
                self.sizes_dim3 != []
                and self.sizes_dim3 < self.sizes_dim3_nus
                and self.sizes_dim3 != 1
            ):
                if self.sizes_dim3 % 2 != 0:
                    self.sizes_dim3 -= 1
                    self.size_indirect_non_reduced.append(self.sizes_dim3 + 1)
                self.size_indirect.append(self.sizes_dim3)
            else:
                self.size_indirect.append(self.sizes_dim3_nus)
                self.size_indirect_non_reduced.append(self.sizes_dim3)

        try:
            self.size_indirect
        except:
            print("Error: TD_INDIRECT not found in acqus file. Unable to determine size of data for indirect dimension. Unable to convert data to NMRPipe format. Please check the acqus file and try again.")
    

        if self.size_indirect != [] and "acqu2s" not in os.listdir("./"):
            self.size_indirect = []

        # Remove values from size indirect if they are equal to 1
        sizes_new = []
        for size in self.size_indirect:
            if size != 1:
                sizes_new.append(size)

        self.size_indirect = sizes_new

        # Try to go through acqu2s and acqu3s and find the nucleus labels and corresponding TD values
        self.indirect_sizes_dict = {}
        if "acqu2s" in os.listdir("./"):
            # Look for TD in the acqu2s file
            nus = False
            td = False
            with open("acqu2s", "r") as file:
                file_lines = file.readlines()

                for i in range(len(file_lines)):
                    if "##$NUSTD=" in file_lines[i] or "##$NusTD=" in file_lines[i]:
                        nus = True
                        line = file_lines[i].split("\n")[0].split()
                        self.sizes_dim2_nus = int(line[1])
                    if "##$TD=" in file_lines[i]:
                        td = True
                        line = file_lines[i].split()
                        self.sizes_dim2 = int(line[1])
                    if "##$NUC1=" in file_lines[i]:
                        line = file_lines[i].split("\n")[0]
                        nuc = line.split("<")[1].split(">")[0]

            try:
                if (
                    td == True
                    and self.sizes_dim2 < self.sizes_dim2_nus
                    and self.sizes_dim2 != 1
                ):
                    if self.sizes_dim2 % 2 != 0:
                        self.sizes_dim2 -= 1
                    self.indirect_sizes_dict[nuc] = self.sizes_dim2
                else:
                    self.indirect_sizes_dict[nuc] = self.sizes_dim2_nus
            except:
                pass

        if "acqu3s" in os.listdir("./"):
            # Look for TD in the acqu2s file
            nus = False
            td = False
            with open("acqu3s", "r") as file:
                file_lines = file.readlines()
                for i in range(len(file_lines)):
                    if "##$NUSTD=" in file_lines[i] or "##$NusTD=" in file_lines[i]:
                        nus = True
                        line = file_lines[i].split("\n")[0].split()
                        self.sizes_dim3_nus = int(line[1])
                    if "##$TD=" in file_lines[i]:
                        td = True
                        line = file_lines[i].split()
                        self.sizes_dim3 = int(line[1])
                    if "##$NUC1=" in file_lines[i]:
                        line = file_lines[i].split("\n")[0]
                        nuc = line.split("<")[1].split(">")[0]

            try:
                if (
                    td == True
                    and self.sizes_dim3 < self.sizes_dim3_nus
                    and self.sizes_dim3 != 1
                ):
                    if nuc not in self.indirect_sizes_dict.keys():
                        if self.sizes_dim3 % 2 != 0:
                            self.sizes_dim3 -= 1
                        self.indirect_sizes_dict[nuc] = self.sizes_dim3
                    else:
                        if self.sizes_dim3 % 2 != 0:
                            self.sizes_dim3 -= 1
                        self.indirect_sizes_dict[nuc + "_1"] = self.sizes_dim3

                else:
                    if nuc not in self.indirect_sizes_dict.keys():
                        self.indirect_sizes_dict[nuc] = self.sizes_dim3_nus
                    else:
                        self.indirect_sizes_dict[nuc + "_1"] = self.sizes_dim3_nus
            except:
                pass

    def find_sw_bruker(self) -> None:
        """
        Reading through the Bruker acqus file to find the sweep widths
        for each spectrum dimension.
        """
        self.sw_indirect = []
        for i in range(len(self.size_indirect) + 1):
            if i == 0:
                for j in range(len(self.acqus_file_lines)):
                    if "##$SW_h=" in self.acqus_file_lines[j]:
                        line = self.acqus_file_lines[j].split()
                        self.sw_direct = float(line[1])
                        break
            if i == 1:
                try:
                    file = open("acqu2s", "r")
                    file_lines = file.readlines()
                    file.close()
                    for j in range(len(file_lines)):
                        if "##$SW_h=" in file_lines[j]:
                            line = file_lines[j].split()
                            self.sw_indirect.append(float(line[1]))
                            break
                except:
                    self.sw_indirect.append(0)
            if i == 2:
                try:
                    file = open("acqu3s", "r")
                    file_lines = file.readlines()
                    file.close()
                    for j in range(len(file_lines)):
                        if "##$SW_h=" in file_lines[j]:
                            line = file_lines[j].split()
                            self.sw_indirect.append(float(line[1]))
                            break
                except:
                    self.sw_indirect.append(0)
            if i == 3:
                try:
                    file = open("acqu4s", "r")
                    file_lines = file.readlines()
                    file.close()
                    for j in range(len(file_lines)):
                        if "##$SW_h=" in file_lines[j]:
                            line = file_lines[j].split()
                            self.sw_indirect.append(float(line[1]))
                            break
                except:
                    self.sw_indirect.append(0)

            if i > 3:
                print("Error: Only able to convert data with up to 4 indirect dimensions. Unable to convert data to NMRPipe format. Please check the acqus file and try again.")


    def find_nucleus_frequencies_bruker(self) -> None:
        """
        Searching through the Bruker acqus file for the Larmor frequency
        of each nucleus recorded.
        """
        self.nucleus_frequencies = []
        nuclei = []
        for i in range(len(self.size_indirect) + 1):
            if i == 0:
                for j in range(len(self.acqus_file_lines)):
                    if "##$SFO1=" in self.acqus_file_lines[j]:
                        line = self.acqus_file_lines[j].split()
                        self.nucleus_frequencies.append(float(line[1]))
                        break
            if i >= 1:
                try:
                    file = open("acqu"+str(i+1) + "s", "r")
                    file_lines = file.readlines()
                    file.close()
                    nucleus=''
                    param = ''
                    for j in range(len(file_lines)):
                        # There is a bug where sometimes O1 is set to 0 in acqu2s even though the
                        # this was not the case. Check if O1=0, if it is, will need to get reference
                        # from the acqus file
                        if("##$NUC1=" in file_lines[j]):
                            line = file_lines[j].split()
                            nucleus = line[1]
                            nuclei.append(nucleus)

                        # if("##$O1=" in file_lines[j]):
                        #     line = file_lines[j].split()
                        #     o1 = float(line[1])
                        #     if(o1==0.0):
                            # search acqus file for nucleus
                        if(nucleus!=''):
                            for j in range(len(self.acqus_file_lines)):
                                if('##$NUC' in self.acqus_file_lines[j]):
                                    if(nucleus in self.acqus_file_lines[j]):
                                        # count of this nucleus already in 
                                        channel = self.acqus_file_lines[j].split('##$NUC')[1].split('=')[0]
                                        param = '##$SFO' + channel + '='
                                
                                if(param!=''):
                                    if param in self.acqus_file_lines[j]:
                                        line = self.acqus_file_lines[j].split()
                                        self.nucleus_frequencies.append(float(line[1]))
                        

                        if "##$SFO1=" in file_lines[j] and param=='':
                            line = file_lines[j].split()
                            # Checking that sfo1_acqus is not equal to bf1
                            self.nucleus_frequencies.append(float(line[1]))
                            break
                except:
                    self.nucleus_frequencies.append(0)
            # if i == 2:
            #     try:
            #         file = open("acqu3s", "r")
            #         file_lines = file.readlines()
            #         file.close()
            #         for j in range(len(file_lines)):
            #             if "##$SFO1=" in file_lines[j]:
            #                 line = file_lines[j].split()
            #                 self.nucleus_frequencies.append(float(line[1]))
            #                 break
            #     except:
            #         self.nucleus_frequencies.append(0)
            # if i == 3:
            #     try:
            #         file = open("acqu4s", "r")
            #         file_lines = file.readlines()
            #         file.close()
            #         for j in range(len(file_lines)):
            #             if "##$SFO1=" in file_lines[j]:
            #                 line = file_lines[j].split()
            #                 self.nucleus_frequencies.append(float(line[1]))
            #                 break
            #     except:
            #         self.nucleus_frequencies.append(0)

    def find_labels_bruker(self) -> None:
        """
        Search through the Bruker acqus file to find the labels for each
        dimension.
        """
        self.labels = []
        for j in range(len(self.acqus_file_lines)):
            if (
                "##$NUC" in self.acqus_file_lines[j]
                and "##$NUCLEUS" not in self.acqus_file_lines[j]
                and "##$NUCLEI" not in self.acqus_file_lines[j]
            ):
                line = self.acqus_file_lines[j].split()
                self.labels.append(line[1].split("<")[1].split(">")[0])

        # If the labels are off, then change them to ID
        for i in range(len(self.labels)):
            if self.labels[i] == "off":
                self.labels[i] = "ID"

        # Work out what the correct order of labels is (i.e. if the nuclei match their frequencies etc)
        self.find_gamma_bruker()
        self.nuc1 = self.labels[0]
        # Find out the field strength using gamma and the carrier frequency
        field = self.nucleus_frequencies[0] / self.gamma[self.nuc1]
        # Find the correct order of labels
        self.labels_correct_order = []
        self.labels_correct_order.append(self.nuc1)

        for i in range(len(self.nucleus_frequencies)):
            if i == 0:
                continue
            else:
                added_nuc = 0
                if self.nucleus_frequencies[i] == 0:
                    self.labels_correct_order.append("ID")
                    added_nuc = 1
                else:
                    for key in self.gamma:
                        if (
                            abs(
                                abs(field * self.gamma[key])
                                - self.nucleus_frequencies[i]
                            )
                            < 0.1
                        ):
                            if key in self.labels_correct_order:
                                self.labels_correct_order.append(key + "_1")
                                added_nuc = 1
                            else:
                                self.labels_correct_order.append(key)
                                added_nuc = 1
                if(added_nuc==0):
                    self.labels_correct_order.append('ID')


    # def find_nucleus_frequencies_bruker(self) -> None:
    #     """
    #     Searching through the Bruker acqus file for the Larmor frequency
    #     of each nucleus recorded.
    #     """
    #     self.nucleus_frequencies = []
    #     for i in range(len(self.size_indirect) + 1):
    #         if i == 0:
    #             for j in range(len(self.acqus_file_lines)):
    #                 if "##$SFO1=" in self.acqus_file_lines[j]:
    #                     line = self.acqus_file_lines[j].split()
    #                     self.nucleus_frequencies.append(float(line[1]))
    #                     break
    #         if i == 1:
    #             try:
    #                 file = open("acqu2s", "r")
    #                 file_lines = file.readlines()
    #                 file.close()
    #                 for j in range(len(file_lines)):
    #                     if "##$SFO1=" in file_lines[j]:
    #                         line = file_lines[j].split()
    #                         self.nucleus_frequencies.append(float(line[1]))
    #                         break
    #             except:
    #                 self.nucleus_frequencies.append(0)
    #         if i == 2:
    #             try:
    #                 file = open("acqu3s", "r")
    #                 file_lines = file.readlines()
    #                 file.close()
    #                 for j in range(len(file_lines)):
    #                     if "##$SFO1=" in file_lines[j]:
    #                         line = file_lines[j].split()
    #                         self.nucleus_frequencies.append(float(line[1]))
    #                         break
    #             except:
    #                 self.nucleus_frequencies.append(0)
    #         if i == 3:
    #             try:
    #                 file = open("acqu4s", "r")
    #                 file_lines = file.readlines()
    #                 file.close()
    #                 for j in range(len(file_lines)):
    #                     if "##$SFO1=" in file_lines[j]:
    #                         line = file_lines[j].split()
    #                         self.nucleus_frequencies.append(float(line[1]))
    #                         break
    #             except:
    #                 self.nucleus_frequencies.append(0)

    # def find_labels_bruker(self) -> None:
    #     """
    #     Search through the Bruker acqus file to find the labels for each
    #     dimension.
    #     """
    #     self.labels = []
    #     for j in range(len(self.acqus_file_lines)):
    #         if (
    #             "##$NUC" in self.acqus_file_lines[j]
    #             and "##$NUCLEUS" not in self.acqus_file_lines[j]
    #             and "##$NUCLEI" not in self.acqus_file_lines[j]
    #         ):
    #             line = self.acqus_file_lines[j].split()
    #             self.labels.append(line[1].split("<")[1].split(">")[0])

    #     # If the labels are off, then change them to ID
    #     for i in range(len(self.labels)):
    #         if self.labels[i] == "off":
    #             self.labels[i] = "ID"

    #     # Work out what the correct order of labels is (i.e. if the nuclei match their frequencies etc)
    #     self.find_gamma_bruker()
    #     self.nuc1 = self.labels[0]
    #     # Find out the field strength using gamma and the carrier frequency
    #     field = self.nucleus_frequencies[0] / self.gamma[self.nuc1]
    #     # Find the correct order of labels
    #     self.labels_correct_order = []
    #     self.labels_correct_order.append(self.nuc1)
    #     for i in range(len(self.nucleus_frequencies)):
    #         if i == 0:
    #             continue
    #         else:
    #             if self.nucleus_frequencies[i] == 0:
    #                 self.labels_correct_order.append("ID")
    #             else:
    #                 for key in self.gamma:
    #                     if (
    #                         abs(
    #                             abs(field * self.gamma[key])
    #                             - self.nucleus_frequencies[i]
    #                         )
    #                         < 1
    #                     ):
    #                         if key in self.labels:
    #                             if key in self.labels_correct_order:
    #                                 self.labels_correct_order.append(key + "_1")
    #                             else:
    #                                 self.labels_correct_order.append(key)

    def find_aqseq(self) -> None:
        """
        The aqseq parameter determines the acquisition order. For 3D or pseudo3D
        data this is essential because the order of acquisition matters for the
        conversion of the data. It can be either 321 or 312 for 3D data and can
        be found in the pulseprogram or pulseprogram.precomp file.
        """

        aqseq_value = 0
        try:
            with open("pulseprogram") as file:
                lines = file.readlines()
                for line in lines:
                    if "aqseq" in line:
                        if line.split()[0] == "aqseq":
                            aqseq_value = line.split()[1]
        except:
            try:
                with open("pulseprogram.precomp") as file:
                    lines = file.readlines()
                    for line in lines:
                        if "aqseq" in line:
                            aqseq_value = line.split()[1]

            except:
                print("Unable to find acquisition order. Using standard ordering.")

        if aqseq_value == "312":
            # The indirect dimensions need reversing
            self.reverse_indirect_dimension_parameters()

    def reverse_indirect_dimension_parameters(self) -> None:
        """
        If aqseq = 312 then the indirect dimension parameters need to
        be reversed.
        """

        self.size_indirect.reverse()
        self.sw_indirect.reverse()
        freq1 = self.nucleus_frequencies[-1]
        freq2 = self.nucleus_frequencies[-2]
        self.nucleus_frequencies[-1] = freq2
        self.nucleus_frequencies[-2] = freq1
        self.acqusition_modes_indirect.reverse()
        lab1 = self.labels_correct_order[-1]
        lab2 = self.labels_correct_order[-2]
        self.labels_correct_order[-1] = lab2
        self.labels_correct_order[-2] = lab1

    def find_gamma_bruker(self) -> None:
        """
        Producing a dictionary containing the gyromagnetic ratio of most NMR
        active nuclei.
        """
        self.gamma = {}
        self.gamma["1H"] = 267.5153151e6  # 267.5221877E6
        self.gamma["19F"] = 251.6628277e6
        self.gamma["13C"] = 67.262e6
        self.gamma["14N"] = 1.93297e7
        self.gamma["15N"] = -27.116e6
        self.gamma["31P"] = 108.282e6
        self.gamma["23Na"] = 70.882e6
        self.gamma["25Mg"] = -1.639e7
        self.gamma["39K"] = 1.2498e7
        self.gamma["41K"] = 0.686e7
        self.gamma["43Ca"] = -1.8025e7
        self.gamma["2H"] = 41.065e6
        self.gamma["7Li"] = 103.962e6
        self.gamma["17O"] = -36.264e6
        self.gamma["10B"] = 2.87471e7
        self.gamma["11B"] = 8.58406e7
        self.gamma["27Al"] = 6.97594e7
        self.gamma["29Si"] = -5.3146e7
        self.gamma["35Cl"] = 2.62401e7
        self.gamma["37Cl"] = 2.18428e7
        self.gamma["50V"] = 2.67164e7
        self.gamma["51V"] = 7.04578e7
        self.gamma["55Mn"] = 6.59777e7
        self.gamma["57Fe"] = 0.86399e7
        self.gamma["59Co"] = 6.3472e7
        self.gamma["63Cu"] = 7.0965e7
        self.gamma["65Cu"] = 7.6018e7
        self.gamma["67Zn"] = 16.767e6
        self.gamma["69Ga"] = 6.43685e7
        self.gamma["71Ga"] = 8.180163e7
        self.gamma["77Se"] = 5.115e7
        self.gamma["79Br"] = 6.70186e7
        self.gamma["81Br"] = 7.22421e7
        self.gamma["103Rh"] = -0.84579e7
        self.gamma["107Ag"] = -1.08718e7
        self.gamma["109Ag"] = -1.25001e7
        self.gamma["111Cd"] = -5.69259e7
        self.gamma["113Cd"] = -5.95504e7
        self.gamma["117Sn"] = -9.57865e7
        self.gamma["119Sn"] = -10.01926
        self.gamma["123Te"] = -7.04893e7
        self.gamma["125Te"] = -8.49722e7
        self.gamma["127I"] = 5.37937e7
        self.gamma["129Xe"] = -7.44069e7
        self.gamma["131Xe"] = 2.20564e7
        self.gamma["183W"] = 1.12070e7
        self.gamma["195Pt"] = 5.80466e7
        self.gamma["197Au"] = 0.4692e7
        self.gamma["199Hg"] = 4.81519e6
        self.gamma["201Hg"] = -1.77748e7
        self.gamma["203Tl"] = 15.43599e7
        self.gamma["205Tl"] = 15.58829e7
        self.gamma["207Pb"] = 5.64661e7

    def find_acquisition_modes_bruker(self) -> None:
        """
        Function to try to determine the acquisition modes of the bruker
        data using the acqus files from the AQ_Mod for the direct dimension
        in acqus and FnMode for the indirect dimensions in acqu2s and acqu3s
        etc.

        For the direct dimension acquisition mode:
        AQ_mod = 0, 1, 2, 3 for QF, QSEQ, QSIM, DQD

        For the indirect dimension acqusition modes:
        FnMode:
        0 = undefined (real),
        1 = QF (real),
        2 = QSED,
        3 = TPPI,
        4 = States,
        5 = States-TPPI,
        6 = Echo-AntiEcho
        """

        self.acqusition_mode_direct = 3  # setting the default value
        self.acqusition_modes_indirect = []
        self.pseudo_flag = 0
        fn_mode = 0

        for line in self.acqus_file_lines:
            if "##$AQ_mod=" in line:
                line = line.split()[1]
                self.acqusition_mode_direct = int(line)
            if "##$FnMode=" in line:
                line = line.split()[1]
                fn_mode = int(line)

        if self.nmrdata.data_dimensions > 1 or self.size_indirect != []:
            self.pseudo_flag = 0

            # Try to read through acqu2s and acqu3s to find the FnMODE parameter

            with open("acqu2s", "r") as file:
                file_lines = file.readlines()
                for line in file_lines:
                    if "##$FnMODE=" in line:
                        line = line.split()[1]
                        if(int(line) < 7):
                            val = line
                        else:
                            val = fn_mode
                        
                        self.acqusition_modes_indirect.append(int(val))
                        if int(val) == 0 or int(val) == 1:
                            self.pseudo_flag += 1
                        break

            if len(self.size_indirect) > 1:
                with open("acqu3s", "r") as file:
                    file_lines = file.readlines()
                    for line in file_lines:
                        if "##$FnMODE=" in line:
                            line = line.split()[1]
                            if(int(line) < 7):
                                val = line
                            else:
                                val = fn_mode
                            self.acqusition_modes_indirect.append(int(val))
                            if int(val) == 0 or int(val) == 1:
                                self.pseudo_flag += 1
                            break

    def find_temperature_bruker(self) -> None:
        """
        Find the temperature the spectrum was recorded at
        """
        self.temperature = 298.15  # Default temperature is 298.15K
        for i in range(len(self.acqus_file_lines)):
            if "##$TE=" in self.acqus_file_lines[i]:
                line = self.acqus_file_lines[i].split()
                self.temperature = float(line[1])
                break

    def calculate_carrier_frequency_bruker(self) -> None:
        """
        Calculate the carrier frequency for each dimension. If the direct
        dimension is proton, calculate the carrier based on water and include
        water referencing as an option.
        """

        if len(self.size_indirect) == 0:
            # return water chemical shift in range 0-100oC
            if (
                self.labels_correct_order[0] == "1H"
                or self.labels_correct_order[0] == "H1"
                or self.labels_correct_order[0] == "H"
            ):
                self.water_ppm = 7.83 - self.temperature / 96.9

                # Use O1/BF1 to calculate a second carrier frequency in case not centred on water
                for j in range(len(self.acqus_file_lines)):
                    if "##$O1=" in self.acqus_file_lines[j]:
                        line = self.acqus_file_lines[j].split()
                        self.O1 = float(line[1])
                        break
                for j in range(len(self.acqus_file_lines)):
                    if "##$BF1=" in self.acqus_file_lines[j]:
                        line = self.acqus_file_lines[j].split()
                        self.BF1 = float(line[1])
                        break
                self.carrier_frequency_1 = self.O1 / self.BF1

                self.references_proton = [self.water_ppm, self.carrier_frequency_1]
                self.references_proton_labels = ["H2O", "O1/BF1"]
            else:
                # Use O1/BF1 to calculate a second carrier frequency in case not centred on water
                for j in range(len(self.acqus_file_lines)):
                    if "##$O1=" in self.acqus_file_lines[j]:
                        line = self.acqus_file_lines[j].split()
                        self.O1 = float(line[1])
                        break
                for j in range(len(self.acqus_file_lines)):
                    if "##$BF1=" in self.acqus_file_lines[j]:
                        line = self.acqus_file_lines[j].split()
                        self.BF1 = float(line[1])
                        break
                self.carrier_frequency_1 = self.O1 / self.BF1

                self.references_proton = [self.carrier_frequency_1]
                self.references_proton_labels = ["O1/BF1"]

        else:
            if (
                self.labels_correct_order[0] == "1H"
                or self.labels_correct_order[0] == "H1"
                or self.labels_correct_order[0] == "H"
            ):
                # return water chemical shift in range 0-100oC
                self.water_ppm = 7.83 - self.temperature / 96.9

            # Use O1/BF1 to calculate a second carrier frequency in case not centred on water
            for j in range(len(self.acqus_file_lines)):
                if "##$O1=" in self.acqus_file_lines[j]:
                    line = self.acqus_file_lines[j].split()
                    self.O1 = float(line[1])
                    break
            for j in range(len(self.acqus_file_lines)):
                if "##$BF1=" in self.acqus_file_lines[j]:
                    line = self.acqus_file_lines[j].split()
                    self.BF1 = float(line[1])
                    break
            self.carrier_frequency_1 = self.O1 / self.BF1

            # Calculate carrier frequencies based on O2/BF2, O3/BF3

            # Search through acqu2s and acqu3s files to find O1 and BF1 values to get chemical shifts
            try:
                with open("acqu2s", "r") as file:
                    file_lines = file.readlines()
                    for i in range(len(file_lines)):
                        if "##$BF1=" in file_lines[i]:
                            line = file_lines[i].split("\n")[0].split()
                            self.BF2 = float(line[1])
                        if "##$O1=" in file_lines[i]:
                            line = file_lines[i].split("\n")[0].split()
                            self.O2 = float(line[1])

            except:
                for j in range(len(self.acqus_file_lines)):
                    if "##$O2=" in self.acqus_file_lines[j]:
                        line = self.acqus_file_lines[j].split()
                        self.O2 = float(line[1])
                    if "##$BF2=" in self.acqus_file_lines[j]:
                        line = self.acqus_file_lines[j].split()
                        self.BF2 = float(line[1])

            try:
                with open("acqu3s", "r") as file:
                    file_lines = file.readlines()
                    for i in range(len(file_lines)):
                        if "##$BF1=" in file_lines[i]:
                            line = file_lines[i].split("\n")[0].split()
                            self.BF3 = float(line[1])
                        if "##$O1=" in file_lines[i]:
                            line = file_lines[i].split("\n")[0].split()
                            self.O3 = float(line[1])

            except:
                for j in range(len(self.acqus_file_lines)):
                    if "##$O3=" in self.acqus_file_lines[j]:
                        line = self.acqus_file_lines[j].split()
                        self.O3 = float(line[1])
                    if "##$BF3=" in self.acqus_file_lines[j]:
                        line = self.acqus_file_lines[j].split()
                        self.BF3 = float(line[1])

            # for j in range(len(self.acqus_file_lines)):
            #     if "##$O3=" in self.acqus_file_lines[j]:
            #         line = self.acqus_file_lines[j].split()
            #         self.O3 = float(line[1])
            #         break
            # for j in range(len(self.acqus_file_lines)):
            #     if "##$BF2=" in self.acqus_file_lines[j]:
            #         line = self.acqus_file_lines[j].split()
            #         self.BF2 = float(line[1])
            #         break
            # for j in range(len(self.acqus_file_lines)):
            #     if "##$BF3=" in self.acqus_file_lines[j]:
            #         line = self.acqus_file_lines[j].split()
            #         self.BF3 = float(line[1])
            #         break
            # Calculate the carrier frequency
            self.carrier_frequency_2 = self.O2 / self.BF2
            self.carrier_frequency_3 = self.O3 / self.BF3

            self.ppms_referenced = []
            self.ppms_referenced_labels = []

            if (
                self.labels_correct_order[0] == "1H"
                or self.labels_correct_order[0] == "H1"
                or self.labels_correct_order[0] == "H"
            ):
                # Calculate referenced carrier frequencies based on water chemical shifts
                self.sfrq0 = self.nucleus_frequencies[0] / (1 + self.water_ppm * 1e-6)
                # Frequency for nucleus A = sfrq0 * gammaH/gammaA
                self.dfrq_13C = self.sfrq0 * 0.251449530
                self.dfrq_15N = self.sfrq0 * 0.101329118
                self.dfrq_P31 = self.sfrq0 * 0.4048064954
                self.dfrq_F19 = self.sfrq0 * 0.9412866605363297

                for i, label in enumerate(self.labels_correct_order):
                    # if(i+1>len(self.nucleus_frequencies)):
                    #     break
                    if label == "15N" or label == "N15":
                        self.ppms_referenced.append(
                            (
                                (self.nucleus_frequencies[i] - self.dfrq_15N)
                                / self.dfrq_15N
                            )
                            * 1e6
                        )
                        self.ppms_referenced_labels.append("15N (Referenced to H2O)")
                    if label == "13C" or label == "C13":
                        self.ppms_referenced.append(
                            (
                                (self.nucleus_frequencies[i] - self.dfrq_13C)
                                / self.dfrq_13C
                            )
                            * 1e6
                        )
                        self.ppms_referenced_labels.append("13C (Referenced to H2O)")
                    if label == "31P" or label == "P31":
                        self.ppms_referenced.append(
                            (
                                (self.nucleus_frequencies[i] - self.dfrq_P31)
                                / self.dfrq_P31
                            )
                            * 1e6
                        )
                        self.ppms_referenced_labels.append("31P (Referenced to H2O)")
                    if label == "19F" or label == "F19":
                        self.ppms_referenced.append(
                            (
                                (self.nucleus_frequencies[i] - self.dfrq_F19)
                                / self.dfrq_F19
                            )
                            * 1e6
                        )
                        self.ppms_referenced_labels.append("19F (Referenced to H2O)")

            if (
                self.labels_correct_order[0] == "1H"
                or self.labels_correct_order[0] == "H1"
                or self.labels_correct_order[0] == "H"
            ):
                self.references_proton = [
                    self.water_ppm,
                    self.carrier_frequency_1,
                    self.carrier_frequency_2,
                    self.carrier_frequency_3,
                ]
                self.references_proton_labels = ["H2O", "O1/BF1", "O2/BF2", "O3/BF3"]
            else:
                self.references_proton = [
                    self.carrier_frequency_1,
                    self.carrier_frequency_2,
                    self.carrier_frequency_3,
                ]
                self.references_proton_labels = ["O1/BF1", "O2/BF2", "O3/BF3"]

            self.references_other = []
            self.references_other_labels = []
            for i, ppm in enumerate(self.ppms_referenced):
                self.references_other.append(ppm)
                self.references_other_labels.append(self.ppms_referenced_labels[i])
            self.references_other.append(self.carrier_frequency_1)
            self.references_other.append(self.carrier_frequency_2)
            self.references_other.append(self.carrier_frequency_3)
            self.references_other_labels.append("O1/BF1")
            self.references_other_labels.append("O2/BF2")
            self.references_other_labels.append("O3/BF3")

    def find_bruker_digital_filter_parameters(self):
        # Search through the acqus file to find the parameters
        acqus_file_lines = open(self.nmrdata.parameter_file, "r").readlines()

        # Find the decim, dspfvs and grpdly parameters
        self.decim = 0
        self.dspfvs = 0
        self.grpdly = 0
        self.remove_filter_before_processing=False
        try:
            for i in range(len(acqus_file_lines)):
                if "##$DECIM=" in acqus_file_lines[i]:
                    line = acqus_file_lines[i].split()
                    self.decim = float(line[1])
                if "##$DSPFVS=" in acqus_file_lines[i]:
                    line = acqus_file_lines[i].split()
                    self.dspfvs = int(line[1])
                if "##$GRPDLY=" in acqus_file_lines[i]:
                    line = acqus_file_lines[i].split()
                    self.grpdly = float(line[1])

                # Finding out what the digitiser mode is (need baseopt - 4) for nmrglue filter removal before processing
                if "##$DIGMOD=" in acqus_file_lines[i]:
                    line = acqus_file_lines[i].split()
                    digimod = int(line[1])
                    if(digimod==4):
                        self.remove_filter_before_processing=True
                        
            
            self.include_digital_filter = True

            

        except:
            self.decim = 0
            self.dspfvs = 0
            self.grpdly = 0
            self.remove_filter_before_processing=False
            self.include_digital_filter = False

    def find_bruker_scaling_parameters(self):
        # Search through the acqus file to find the NS and NC parameters
        acqus_file_lines = open(self.nmrdata.parameter_file, "r").readlines()

        # Find the NS and NC parameters
        try:
            for i in range(len(acqus_file_lines)):
                if "##$NS=" in acqus_file_lines[i]:
                    line = acqus_file_lines[i].split()
                    self.NS = int(line[1])
                if "##$NC=" in acqus_file_lines[i]:
                    line = acqus_file_lines[i].split()
                    self.NC = int(line[1])
            self.include_scaling = True
        except:
            self.NS = 0
            self.NC = 0
            self.include_scaling = False

    def determine_byte_order(self):
        try:
            with open(self.nmrdata.parameter_file) as file:
                lines = file.readlines()
                for line in lines:
                    if ("##$BYTORDA") in line:
                        self.byte_order = line.split("\n")[0].split()[-1]
                        break
        except:
            self.byte_order = 0

    def determine_byte_size(self):
        try:
            with open(self.nmrdata.parameter_file) as file:
                lines = file.readlines()
                for line in lines:
                    if ("##$DTYPA") in line:
                        self.d_type = line.split("\n")[0].split()[-1]
                        break
        except:
            self.d_type = 0
