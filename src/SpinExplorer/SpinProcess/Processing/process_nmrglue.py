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


"""
    The functions listed below were originally obtained from nmrglue,
    followed by editing and customisation:
    remove_digital_filter, rm_digital_filter, sol_general_nd, sol_general,
    suppress_solvent_3d, ext, base, base2, lp, lp2, zf
    
    Copyright Notice and Statement for the nmrglue Project
    Copyright (c) 2010-2015 Jonathan J. Helmus
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:


    a. Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.


    b. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the
    distribution.


    c. Neither the name of the author nor the names of contributors may
    be used to endorse or promote products derived from this software
    without specific prior written permission.


    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    """


from SpinExplorer.SpinConverter.Conversion.convert_nmrglue import Convert_nmrglue


import wx
import numpy as np
import nmrglue as ng
import os
import json
import copy
import traceback
import shutil


class ProcessNMRGlue:

    def __init__(
        self,
        notebook,
        dimension_tabs,
        nmr_data,
        interactive_phasing=False,
    ):
        """
        This class will perform nmrglue processing of the nmr data
        based on current SpinProcess parameters.
        """
        self.notebook = notebook
        self.dimension_tabs = dimension_tabs
        self.nmr_data = nmr_data

        if interactive_phasing == True:
            # Just initialising the class for bruker digital filter removal
            return

        try:
            # # If Bruker, need to re-convert the data using the correct phasing in order to remove frowns in the direct dimension
            # converter = Convert_nmrglue()
            
            # Now need to process the data. Phasing in the direct dimension can be ignored if Bruker because this has 
            # been applied during conversion
            self.on_run_processing_nmrglue()
            self.success_output_message()

        except:
            message = traceback.format_exc()
            self.fail_output_message(message)

        if self.notebook.parent.original_frame != None:
            self.notebook.parent.Destroy()

    def success_output_message(self):
        """
        Provides an output message to the user to say that the
        conversion is complete.
        """
        dlg = wx.MessageDialog(
            self.notebook,
            "Data processing using nmrglue is complete.",
            "Complete",
            wx.OK,
        )
        self.notebook.Raise()
        self.notebook.SetFocus()
        dlg.ShowModal()
        dlg.Destroy()



    def fail_output_message(self, message):
        """
        Provides an output message to the user to say that the
        conversion did not work correctly
        """
        dlg = wx.MessageDialog(
            self.notebook,
            "Data processing using nmrglue did not complete correctly. Traceback:\n\n" + message,
            "Error",
            wx.OK,
        )
        self.notebook.Raise()
        self.notebook.SetFocus()
        dlg.ShowModal()
        dlg.Destroy()

    def on_run_processing_nmrglue(self):
        """
        Apply the processing parameters to the data
        """
        self.apply_processing_parameters()

    def apply_processing_parameters(self):
        # Process the data according to the user inputted processing parameters

        # Initial FID data
        dic, data = copy.deepcopy(self.nmr_data.dic), copy.deepcopy(self.nmr_data.data)

        # Checking whether processing is required for second/third dimensions
        include_dim2, include_dim3 = self.checking_dimensions()

        # Set the comment to nmrglue so that the fact nmrglue processing was used is noted in the processed spectrum header
        dic['FDCOMMENT'] = 'nmrglue'
        if self.nmr_data.pseudo_axis == True:
            # Adding the fact that there is a pseudo axis to the FDCOMMENT
            dic['FDCOMMENT'] += '_pseudo'

        # For the direct dimension, apply the processing functions
        dic, data = self.apply_dimension_processing(
            dic, data, 0, self.dimension_tabs[0]
        )

        # Adding processing lines for the first complex indirect dimension
        if include_dim2:
            check_nus = self.check_nus(1)
            if check_nus != [0]:
                # Give an error saying that NUS reconstruction is not currently supported using nmrglue processing.
                self.nus_nmrglue_error()
                return
            else:
                if self.nmr_data.pseudo_axis == False and include_dim3 == False:
                    dic, data = ng.pipe_proc.tp(dic, data)
                elif self.nmr_data.pseudo_axis == True:
                    if self.nmr_data.index == 2:
                        dic, data = ng.pipe_proc.tp(dic, data)
                    elif self.nmr_data.index == 1:
                        # If the pseudo axis is the central axis then need to move the third axis
                        dic, data = self.transpose_3d(dic, data, auto=True)
                        dic, data = self.zero_transpose_3d(dic, data)
                else:
                    dic, data = self.transpose_3d(dic, data, auto=True)

            dic, data = self.apply_dimension_processing(
                dic, data, 1, self.dimension_tabs[1]
            )

        if include_dim3:
            dic, data = self.zero_transpose_3d(dic, data)
            dic, data = self.apply_dimension_processing(
                dic, data, 2, self.dimension_tabs[2]
            )
            dic, data = self.zero_transpose_3d(dic, data)
            dic1 = copy.deepcopy(dic)
            data1 = copy.deepcopy(data)
            self.create_3D_projections(dic1, data1)

        if (
            include_dim2 == True
            and self.nmr_data.pseudo_axis == True
            and self.nmr_data.index == 1
        ):
            dic, data = self.zero_transpose_3d(dic, data)
            dic1 = copy.deepcopy(dic)
            data1 = copy.deepcopy(data)
            self.create_3D_projections(dic1, data1)

        self.write_output(dic, data)

        original_frame = []
        if self.notebook.parent.original_frame != None:
            self.update_spinview_frame()

    def apply_dimension_processing(self, dic, data, dimension, dimension_tab):
        """
        Applying NMR processing to a given dimension using nmrglue
        functions. May need to be altered for 3D datasets
        """

        if dimension == 0:
            dic, data = self.add_solvent_suppression(
                dic, data, dimension, dimension_tab
            )
        dic, data = self.add_linear_prediction(dic, data, dimension, dimension_tab)
        dic, data = self.add_apodization(dic, data, dimension, dimension_tab)
        dic, data = self.add_zero_filling(dic, data, dimension, dimension_tab)
        dic, data = self.add_fourier_transform(dic, data, dimension, dimension_tab)
        dic, data = self.add_phasing(dic, data, dimension, dimension_tab)
        dic, data = self.add_extraction(dic, data, dimension, dimension_tab)
        dic, data = self.add_baseline_correction(dic, data, dimension, dimension_tab)

        return dic, data

    def write_output(self, dic, data):
        """
        Writing the processed data to the nmrfile.
        Also setting the dictionary quadrature flags
        to real as fourier transforms should have been
        applied to all dimensions.
        """

        dim = self.nmr_data.dim
        pseudo_axis = self.nmr_data.pseudo_axis

        if dim == 1:
            nmrfile = "test.ft"
        elif dim == 2 and pseudo_axis == False:
            nmrfile = "test.ft2"
        elif dim == 2 and pseudo_axis == True:
            nmrfile = "test.ft"
        elif dim == 3 and pseudo_axis == True:
            nmrfile = "test.ft2"
        elif dim == 3 and pseudo_axis == False:
            nmrfile = "test.ft3"
        else:
            nmrfile = "test.ft"

        # Data is now processed to can set all QUAD flags to 1 (Real)
        dic["FDF1QUADFLAG"] = 1.0
        dic["FDF2QUADFLAG"] = 1.0
        dic["FDF3QUADFLAG"] = 1.0
        dic["FDQUADFLAG"] = 1.0
        dic["FDFILECOUNT"] = 1

        data = data.real
        data = data.astype(np.float32)

        ng.pipe.write(nmrfile, dic, data, overwrite=True)

    def check_nus(self, dimension) -> list[int]:
        """
        Checking if NUS reconstruction has been selected.
        Output:
        [0] - no NUS reconstruction for any indirect dimension
        [1] - NUS reconstruction for the first indirect dimension
        [2] - NUS reconstruction for the second indirect dimension
              (if present)
        [1,2] - NUS reconstruction for the first and second indirect
                dimensions.

        The code will also check the nusfile to ensure that the length
        of the rows in nusfile are constistent with the length of the
        output array. This prevents errors where  if there are 2 NUS
        dimensions, but the user has only selected to reconstruct one
        dimension. The user will get a warning to either turn off
        NUS reconstruction or to select both dimensions for NUS
        reconstruction in the graphical interface.
        """
        if (
            self.dimension_tabs[
                dimension
            ].linear_prediction.linear_prediction_radio_box_indirect.GetSelection()
            == 2
        ):
            try:
                if (
                    self.dimension_tabs[
                        dimension + 1
                    ].linear_prediction.linear_prediction_radio_box_indirect.GetSelection()
                    == 2
                ):
                    nus = [1, 2]
                else:
                    nus = [1]
            except:
                nus = [1]
        else:
            nus = [0]

        return nus

    def checking_dimensions(self) -> list[bool]:
        """
        Based on the number of SpinProcess tabs, determine whether
        second/third dimension processing is required.
        """
        include_dim2 = False
        include_dim3 = False
        if self.nmr_data.dim == 2 and self.nmr_data.pseudo_axis == False:
            include_dim2 = True
        elif self.nmr_data.dim == 3:
            include_dim2 = True
        else:
            include_dim2 = False
        if self.nmr_data.dim == 3 and self.nmr_data.pseudo_axis == False:
            include_dim3 = True
        else:
            include_dim3 = False

        return include_dim2, include_dim3

    def nus_nmrglue_error(self):
        """
        Outputting an error informing the user that NUS reconstruction is not
        currently supported for nmrglue processing.
        """

        dlg = wx.MessageDialog(
            self.notebook,
            "NUS reconstruction is not currently supported for nmrglue processing. Process the data using the nmrPipe processing function.",
            "Warning",
            wx.OK | wx.ICON_WARNING,
        )
        self.notebook.Raise()
        self.notebook.SetFocus()
        result = dlg.ShowModal()

    def update_spinview_frame(self):
        """
        Updating the SpinView frame after processing using the
        re-process button.
        """
        original_frame = True
        self.notebook.parent.original_frame.Enable()
        path = self.notebook.parent.original_frame.parent.path
        cwd = self.notebook.parent.original_frame.parent.cwd
        self.notebook.parent.original_frame.parent.reprocess = True
        self.notebook.parent.original_frame.parent.Destroy()
        if self.notebook.parent.original_frame.parent.path != "":
            os.chdir(self.notebook.parent.original_frame.parent.path)
        from SpinExplorer.SpinView.SpinView import SpinView

        app = SpinView()
        if self.notebook.parent.original_frame.parent.cwd != "":
            app.path = path
            app.cwd = cwd
        

    def create_3D_projections(self, dic, data):
        """
        This function will form skyline projections over the data along a given
        axis. e.g. a HNCO will have H-N, H-CO, N-CO planes.
        """
        # Move all existing .dat files to a folder called OldProjections

        current_dir = os.getcwd()
        old_dir = os.path.join(current_dir, "OldProjections")
        
        # Create 'Old' directory if it doesn't exist
        os.makedirs(old_dir, exist_ok=True)
        
        # Loop through files in the current directory
        for filename in os.listdir(current_dir):
            if filename.endswith(".dat") and os.path.isfile(filename):
                source = os.path.join(current_dir, filename)
                destination = os.path.join(old_dir, filename)
                shutil.move(source, destination)

        data0 = np.max(data, axis=0)
        dic0 = copy.deepcopy(dic)

        dim_0 = dic["FDDIMORDER"][2]
        fn = "FDF" + str(int(dim_0))
        dic0["FDDIMCOUNT"] = 2
        dic0[fn + "SIZE"] = 0
        dic0[fn + "TDSIZE"] = 0
        dic0[fn + "FTSIZE"] = 0
        dic0[fn + "APOD"] = 0
        dic0[fn + "APODSIZE"] = 0
        dic0[fn + "SW"] = 0
        dic0[fn + "CENTER"] = 0
        # dic0[fn + "LABEL"] = ""

        name = dic0["FDF1LABEL"] + "." + dic0["FDF2LABEL"] + ".dat"

        ng.pipe.write(name, dic0, data0, overwrite=True)

        dic1, data1 = self.zero_transpose_3d(dic, data)
        data1_1 = np.max(data1, axis=0)
        dic1 = copy.deepcopy(dic1)

        dim_1 = dic1["FDDIMORDER"][2]
        fn = "FDF" + str(int(dim_1))
        dic1["FDDIMCOUNT"] = 2
        dic1[fn + "SIZE"] = 0
        dic1[fn + "TDSIZE"] = 0
        dic1[fn + "FTSIZE"] = 0
        dic1[fn + "APOD"] = 0
        dic1[fn + "APODSIZE"] = 0
        dic1[fn + "SW"] = 0
        dic1[fn + "CENTER"] = 0
        # dic1[fn + "LABEL"] = ""

        name = dic1["FDF3LABEL"] + "." + dic1["FDF2LABEL"] + ".dat"

        ng.pipe.write(name, dic1, data1_1, overwrite=True)

        dic2, data2 = self.transpose_3d(dic, data)
        dic2, data2 = self.zero_transpose_3d(dic2, data2)

        data2_1 = np.max(data2, axis=0)
        dic2 = copy.deepcopy(dic2)

        dim_2 = dic2["FDDIMORDER"][2]
        fn = "FDF" + str(int(dim_2))
        dic2["FDDIMCOUNT"] = 2
        dic2[fn + "SIZE"] = 0
        dic2[fn + "TDSIZE"] = 0
        dic2[fn + "FTSIZE"] = 0
        dic2[fn + "APOD"] = 0
        dic2[fn + "APODSIZE"] = 0
        dic2[fn + "SW"] = 0
        dic2[fn + "CENTER"] = 0
        # dic2[fn + "LABEL"] = ""

        name = dic2["FDF1LABEL"] + "." + dic2["FDF3LABEL"] + ".dat"

        ng.pipe.write(name, dic2, data2_1, overwrite=True)

        # front = np.max(data, axis=1)
        # side = np.max(data, axis=2)

    def add_solvent_suppression(self, dic, data, dimension, dimension_tab):
        """
        1 - checking if the solvent suppression checkbox is ticked
        2 - if it is ticked, then perform solvent suppression on the
            direct dimension
        """
        tab = dimension_tab.solvent_suppression

        if tab.solvent_suppression_checkbox.GetValue() == True:
            # Apply solvent suppression
            data_orgiginal = data
            if tab.solvent_suppression_filter_selection == 0:
                filter_size = int(
                    tab.solvent_suppression_filter_length
                )  # Larger filter in time domain is larger filter in the frequency domain
                if int(tab.solvent_suppression_lowpass_shape_selection) + 1 == 1:
                    from scipy.signal.windows import boxcar

                    filter = boxcar(filter_size)
                elif int(tab.solvent_suppression_lowpass_shape_selection) + 1 == 2:
                    filter = np.cos(np.pi * np.linspace(-0.5, 0.5, filter_size))
                else:
                    filter = np.cos(np.pi * np.linspace(-0.5, 0.5, filter_size)) ** 2

                data = self.sol_general(data, filter, w=filter_size, mode="same")

            else:
                # Give an error saying that the selected solvent suppression filter is not supported for windows processing, please use a machine containing nmrPipe
                dlg = wx.MessageDialog(
                    self.notebook,
                    "The selected solvent suppression filter is not supported for windows processing. Continuing without solvent suppression. Please change to low bandpass filter or use a machine containing nmrPipe.",
                    "Warning",
                    wx.OK | wx.ICON_WARNING,
                )
                self.Raise()
                self.SetFocus()
                result = dlg.ShowModal()

        return dic, data

    def add_linear_prediction(self, dic, data, dimension, dimension_tab):
        """
        1 - checking if the linear prediction checkbox is ticked
            for the direct dimension or if the radio box is on
            linear prediction for the indirect dimension.
        2 - if it is ticked, then perform linear prediction
        """
        tab = dimension_tab.linear_prediction

        if dimension == 0:
            check = tab.linear_prediction_checkbox.GetValue()
            if check == False:
                return dic, data
            option = tab.linear_prediction_options_selection
            coefficient_option = tab.linear_prediction_coefficients_selection
        else:
            selection = tab.linear_prediction_radio_box_indirect.GetSelection()
            if selection != 1:
                return dic, data
            option = tab.linear_prediction_indirect_options_selection
            coefficient_option = tab.linear_prediction_indirect_coefficients_selection

        # Apply linear prediction

        if option == 0:
            append = "after"
        else:
            append = "before"
        if coefficient_option == 0:
            mode = "f"
        elif coefficient_option == 1:
            mode = "b"
        else:
            mode = "fb"
        dic, data = self.lp(dic, data, pred="default", mode=mode, append=append)

        return dic, data

    def add_apodization(self, dic, data, dimension, dimension_tab):
        """
        1 - checking if the apodization check box is ticked
        2 - if it is ticked, then perform the selected apodization
        """

        tab = dimension_tab.apodization

        if tab.apodization_checkbox.GetValue() == True:
            # Apply apodization
            if tab.apodization_combobox_selection == 0:
                # No apodization, just 1st point scaling
                dic, data = ng.pipe_proc.em(
                    dic,
                    data,
                    lb=0.0,
                    c=float(tab.apodization_first_point_scaling),
                )
            elif tab.apodization_combobox_selection == 1:
                # Exponential line broadening
                dic, data = ng.pipe_proc.em(
                    dic,
                    data,
                    lb=float(tab.exponential_line_broadening),
                    c=float(tab.apodization_first_point_scaling),
                )
            elif tab.apodization_combobox_selection == 2:
                # Lorentz to gauss apodization
                dic, data = ng.pipe_proc.gm(
                    dic,
                    data,
                    g1=float(tab.g1),
                    g2=float(tab.g2),
                    g3=float(tab.g3),
                    c=float(tab.apodization_first_point_scaling),
                )
            elif tab.apodization_combobox_selection == 3:
                # Sinebell apodization
                dic, data = ng.pipe_proc.sp(
                    dic,
                    data,
                    off=float(tab.offset),
                    end=float(tab.end),
                    pow=int(tab.power),
                    c=float(tab.apodization_first_point_scaling),
                )
            elif tab.apodization_combobox_selection == 4:
                # Gaussian broadening apodization
                dic, data = ng.pipe_proc.gmb(
                    dic,
                    data,
                    lb=float(tab.a),
                    gb=float(tab.b),
                    c=float(tab.apodization_first_point_scaling),
                )
            elif tab.apodization_combobox_selection == 5:
                # Trapezoid apodization
                dic, data = ng.pipe_proc.tp(
                    dic,
                    data,
                    t1=float(tab.t1),
                    t2=float(tab.t2),
                    c=float(tab.apodization_first_point_scaling),
                )
            elif tab.apodization_combobox_selection == 6:
                # Triangle apodization
                dic, data = ng.pipe_proc.tri(
                    dic,
                    data,
                    loc=float(tab.loc),
                    c=float(tab.apodization_first_point_scaling),
                )

        return dic, data

    def add_zero_filling(self, dic, data, dimension, dimension_tab):
        """
        1 - checking if the zero filling check box is ticked
        2 - if it is ticked, then perform the selected zero filling
        """

        tab = dimension_tab.zero_filling

        if tab.zero_filling_checkbox.GetValue() == True:
            if tab.zero_filling_round_checkbox.GetValue() == True:
                round = True
            else:
                round = False
            if tab.zero_filling_combobox_selection == 0:
                dic, data = self.zf(
                    dic,
                    data,
                    zf=int(tab.zero_filling_value_doubling_times),
                    auto=round,
                )
            elif tab.zero_filling_combobox_selection == 1:
                dic, data = self.zf(
                    dic,
                    data,
                    pad=int(tab.zero_filling_value_zeros_to_add),
                    auto=round,
                )
            elif tab.zero_filling_combobox_selection == 2:
                dic, data = self.zf(
                    dic,
                    data,
                    size=int(tab.zero_filling_value_final_data_size),
                    auto=round,
                )

        return dic, data

    def add_fourier_transform(self, dic, data, dimension, dimension_tab):
        """
        1 - checking if the fourier transform check box is ticked
        2 - if it is ticked, then perform the selected fourier transform
        3 - if the digital filter removal post-FT option was chosen,
            remove the digital filter
        """

        tab = dimension_tab.fourier_transform

        if tab.fourier_transform_checkbox.GetValue() == True:
            if tab.ft_method_selection == 0:
                dic, data = ng.pipe_proc.ft(dic, data)
            if tab.ft_method_selection == 1:
                dic, data = ng.pipe_proc.ft(dic, data, auto=True)
            elif tab.ft_method_selection == 2:
                dic, data = ng.pipe_proc.ft(dic, data, real=True)
            elif tab.ft_method_selection == 3:
                dic, data = ng.pipe_proc.ft(dic, data, inv=True)
            elif tab.ft_method_selection == 4:
                dic, data = ng.pipe_proc.ft(dic, data, alt=True)
            elif tab.ft_method_selection == 5:
                dic, data = ng.pipe_proc.ft(dic, data, neg=True)
            elif tab.ft_method_selection == 6:
                dic, data = ng.pipe_proc.ft(dic, data, alt=True, neg=True)

        if dimension == 0:
            digital_filter_removal = self.check_digital_filter_removal()
            if digital_filter_removal == True:
                dic_bruker, dat_bruker = ng.bruker.read("./")
                data = self.remove_digital_filter(dic_bruker, data, truncate=False)

        return dic, data

    def check_digital_filter_removal(self) -> bool:
        """
        Try to search through the conversion parameter dictionary to see if the
        digital filter removal was selected before or after fourier transform
        in the SpinConverter page. Default output will be True.
        """

        try:
            with open("parameters.json", "r") as file:
                parameter_dictionary = json.load(file)["conversion"]
            if (
                parameter_dictionary["digital filter parameters"][
                    "Remove Digital Filter"
                ]
                == True
            ):
                if (
                    parameter_dictionary["digital filter parameters"][
                        "Remove Before/After Fourier Transform"
                    ]
                    == "After"
                ):
                    return True
                else:
                    return False
            return False

        except:
            return True

    def add_phasing(self, dic, data, dimension, dimension_tab):
        """
        1 - checking if the phasing check box is ticked
        2 - if it is ticked, then perform the selected phasing
        """
        tab = dimension_tab.phasing

        if dimension == 0:
            check = tab.phase_correction_checkbox.GetValue()
            p0 = float(tab.phase_correction_p0_textcontrol.GetValue())
            p1 = float(tab.phase_correction_p1_textcontrol.GetValue())
        else:
            check = tab.phase_correction_checkbox_indirect.GetValue()
            p0 = float(tab.phase_correction_p0_textcontrol_indirect.GetValue())
            p1 = float(tab.phase_correction_p1_textcontrol_indirect.GetValue())


        if check == True:
            dic, data = ng.pipe_proc.ps(
                dic,
                data,
                p0=p0,
                p1=p1,
            )

        if dimension == 0:
            # Magnitude mode for the direct dimension
            if tab.magnitude_mode_checkbox.GetValue() == True:
                dic, data = ng.pipe_proc.mc(dic, data)

        dic, data = ng.pipe_proc.di(dic, data)

        return dic, data

    def add_extraction(self, dic, data, dimension, dimension_tab):
        """
        1 - checking if the extraction check box is ticked
        2 - if it is ticked, then perform the selected extraction
        """

        tab = dimension_tab.extraction

        if tab.extraction_checkbox.GetValue() == True:
            # Find the indexes of the ppm values selected
            # Get the ppm values from the data
            ppm_values = ng.pipe.make_uc(dic, data, dim=len(data.shape) - 1)
            ppm_values = ppm_values.ppm_scale()
            x_initial = np.abs(
                ppm_values - float(tab.extraction_ppm_start_textcontrol.GetValue())
            ).argmin()
            x_final = np.abs(
                ppm_values - float(tab.extraction_ppm_end_textcontrol.GetValue())
            ).argmin()
            if x_initial > x_final:
                x_initial, x_final = x_final, x_initial
            # Change x_initial and x_final so that the difference is an even number
            if (x_final - x_initial + 1) % 2 != 0:
                x_final += 1
            dic, data = self.ext(dic, data, x1=x_initial, xn=x_final, sw=True)
        return dic, data

    def add_baseline_correction(self, dic, data, dimension, dimension_tab):
        """
        1 - checking if the baseline correction check box is ticked
        2 - if it is ticked, then perform the selected extraction
        """

        tab = dimension_tab.baseline_correction

        if tab.baseline_correction_checkbox.GetValue() == True:
            if tab.baseline_correction_radio_box_selection == 1:
                # If POLY baseline correction is selected, this is not currently supported on nmrglue
                message = "The selected baseline correction method is not supported for nmrglue processing. Continuing without baselining. Please use a machine containing nmrPipe or use a linear baselining method."
                dlg = wx.MessageDialog(
                    self.notebook, message, "Warning", wx.OK | wx.ICON_WARNING
                )
                self.Raise()
                self.SetFocus()
                result = dlg.ShowModal()
                return dic, data

            # Split the node list
            node_list = tab.baseline_correction_node_list_textcontrol.GetValue()
            node_list = node_list.split(",")
            node_list_final = []
            for node in node_list:
                node_list_final.append(float(node))

            # Convert nodes into points
            node_list_final = np.array(node_list_final)
            node_list_final = (node_list_final / 100) * data.shape[-1]
            node_list_final = node_list_final.astype(int)
            # Replace any zeros with a number greater than 1 to allow the nmrglue baselining routines to work correctly
            node_list_final[node_list_final == 0] = (
                int(tab.baseline_correction_nodes_textcontrol.GetValue()) + 1
            )

            dic, data = self.base(
                dic,
                data,
                nl=node_list_final,
                nw=int(tab.baseline_correction_nodes_textcontrol.GetValue()),
            )

        return dic, data

    

    """
    Obtained from nmrglue followed by customisation
    
    Copyright Notice and Statement for the nmrglue Project
    Copyright (c) 2010-2015 Jonathan J. Helmus
    All rights reserved.
    """

    def remove_digital_filter(self, dic, data, truncate=True):
        """
        Remove the digital filter from Bruker data.

        Parameters
        ----------
        dic : dict
            Dictionary of Bruker parameters.
        data : ndarray
            Array of NMR data to remove digital filter from.
        truncate : bool, optional
            True to truncate the phase shift prior to removing the digital filter.
            This typically produces a better looking spectrum but may remove
            useful data.  False uses a non-truncated phase.
        post_proc : bool, optional
            True if the digital filter is to be removed post processing, i.e after
            fourier transformation. The corrected FID will not be returned, only a
            corrected spectrum in the frequency dimension will be returned

        Returns
        -------
        ndata : ndarray
            Array of NMR data with digital filter removed

        See Also
        ---------
        rm_dig_filter : Remove digital filter by specifying parameters.

        """
        if "acqus" not in dic:
            raise ValueError("dictionary does not contain acqus parameters")

        if "DECIM" not in dic["acqus"]:
            raise ValueError("dictionary does not contain DECIM parameter")
        decim = dic["acqus"]["DECIM"]

        if "DSPFVS" not in dic["acqus"]:
            raise ValueError("dictionary does not contain DSPFVS parameter")
        dspfvs = dic["acqus"]["DSPFVS"]

        if "GRPDLY" not in dic["acqus"]:
            grpdly = 0
        else:
            grpdly = dic["acqus"]["GRPDLY"]

        return self.rm_dig_filter(data, decim, dspfvs, grpdly, truncate)
    

    """
    Obtained from nmrglue followed by customisation
    
    Copyright Notice and Statement for the nmrglue Project
    Copyright (c) 2010-2015 Jonathan J. Helmus
    All rights reserved.
    """

    def rm_dig_filter(self, data, decim, dspfvs, grpdly=0, truncate_grpdly=True):
        """
        Remove the digital filter from Bruker data.

        Parameters
        ----------
        data : ndarray
            Array of NMR data to remove digital filter from.
        decim : int
            Decimation rate (Bruker DECIM parameter).
        dspfvs : int
            Firmware version (Bruker DSPFVS parameter).
        grpdly : float, optional
            Group delay. (Bruker GRPDLY parameter). When non-zero decim and
            dspfvs are ignored.
        truncate_grpdly : bool, optional
            True to truncate the value of grpdly provided or determined from
            the decim and dspfvs parameters before removing the digital filter.
            This typically produces a better looking spectrum but may remove useful
            data.  False uses a non-truncated grpdly value.
        post_proc : bool, optional
            True if the digital filter is to be removed post processing, i.e after
            fourier transformation. The corrected time domain data will not be
            returned, only the corrected spectrum in the frequency dimension will
            be returned

        Returns
        -------
        ndata : ndarray
            Array of NMR data with digital filter removed.

        See Also
        --------
        remove_digital_filter : Remove digital filter using Bruker dictionary.

        """
        #    A first order phase correction equal to 2*PI*GRPDLY is applied to the
        #    data and the time-corrected FT data is returned

        # The frequency dimension will have the same number of points as the
        # original time domain data, but the time domain data will remain
        # uncorrected
        # -----------------------------------------------------------------------

        if grpdly > 0:  # use group delay value if provided (not 0 or -1)
            phase = grpdly

        # determine the phase correction
        else:
            if dspfvs >= 14:  # DSPFVS greater than 14 give no phase correction.
                phase = 0.0
            else:  # loop up the phase in the table
                bruker_dsp_table = {
                    10: {
                        2: 44.75,
                        3: 33.5,
                        4: 66.625,
                        6: 59.083333333333333,
                        8: 68.5625,
                        12: 60.375,
                        16: 69.53125,
                        24: 61.020833333333333,
                        32: 70.015625,
                        48: 61.34375,
                        64: 70.2578125,
                        96: 61.505208333333333,
                        128: 70.37890625,
                        192: 61.5859375,
                        256: 70.439453125,
                        384: 61.626302083333333,
                        512: 70.4697265625,
                        768: 61.646484375,
                        1024: 70.48486328125,
                        1536: 61.656575520833333,
                        2048: 70.492431640625,
                    },
                    11: {
                        2: 46.0,
                        3: 36.5,
                        4: 48.0,
                        6: 50.166666666666667,
                        8: 53.25,
                        12: 69.5,
                        16: 72.25,
                        24: 70.166666666666667,
                        32: 72.75,
                        48: 70.5,
                        64: 73.0,
                        96: 70.666666666666667,
                        128: 72.5,
                        192: 71.333333333333333,
                        256: 72.25,
                        384: 71.666666666666667,
                        512: 72.125,
                        768: 71.833333333333333,
                        1024: 72.0625,
                        1536: 71.916666666666667,
                        2048: 72.03125,
                    },
                    12: {
                        2: 46.0,
                        3: 36.5,
                        4: 48.0,
                        6: 50.166666666666667,
                        8: 53.25,
                        12: 69.5,
                        16: 71.625,
                        24: 70.166666666666667,
                        32: 72.125,
                        48: 70.5,
                        64: 72.375,
                        96: 70.666666666666667,
                        128: 72.5,
                        192: 71.333333333333333,
                        256: 72.25,
                        384: 71.666666666666667,
                        512: 72.125,
                        768: 71.833333333333333,
                        1024: 72.0625,
                        1536: 71.916666666666667,
                        2048: 72.03125,
                    },
                    13: {
                        2: 2.75,
                        3: 2.8333333333333333,
                        4: 2.875,
                        6: 2.9166666666666667,
                        8: 2.9375,
                        12: 2.9583333333333333,
                        16: 2.96875,
                        24: 2.9791666666666667,
                        32: 2.984375,
                        48: 2.9895833333333333,
                        64: 2.9921875,
                        96: 2.9947916666666667,
                    },
                }
                if dspfvs not in bruker_dsp_table:
                    raise ValueError("dspfvs not in lookup table")
                if decim not in bruker_dsp_table[dspfvs]:
                    raise ValueError("decim not in lookup table")
                phase = bruker_dsp_table[dspfvs][decim]

        if truncate_grpdly:  # truncate the phase
            phase = np.floor(phase)

        s = data.shape[-1]
        pdata = data * np.exp(-2.0j * np.pi * phase * np.arange(s) / s)
        pdata = pdata.astype(data.dtype)
        return pdata

    """
    Obtained from nmrglue followed by customisation
    
    Copyright Notice and Statement for the nmrglue Project
    Copyright (c) 2010-2015 Jonathan J. Helmus
    All rights reserved.
    """

    def sol_general_nd(self, data, filter, axis=-1, mode="same"):
        """
        Generalized solvent suppression filter for N-D data.

        Applies solvent suppression along a specific axis using convolution.

        Parameters
        ----------
        data : ndarray
            N-D array of NMR data.
        filter : ndarray
            1D filter array to convolve with.
        axis : int, optional
            Axis along which to apply the filter (default: last axis).
        mode : {'valid', 'same', 'full'}, optional
            Convolution mode (usually 'same').

        Returns
        -------
        ndata : ndarray
            Filtered NMR data.
        """
        A = filter.sum()
        if A == 0:
            raise ValueError("Filter sum cannot be zero.")

        filtered_data = np.zeros(shape=data.shape)

        # Apply filter to each trace
        import scipy

        if len(data.shape) == 1:
            filtered_data = data - scipy.signal.convolve(data, filter, mode=mode) / A
        elif len(data.shape) == 2:
            for i, dat in enumerate(data):
                filtered_data[i] = (
                    dat - scipy.signal.convolve(dat, filter, mode=mode) / A
                )
        elif len(data.shape) == 3:
            for j, dat in enumerate(data):
                for k, dat2 in enumerate(dat):
                    filtered_data[j][k] = (
                        dat2 - scipy.signal.convolve(dat2, filter, mode=mode) / A
                    )

        return filtered_data

    """
    Obtained from nmrglue followed by customisation
    
    Copyright Notice and Statement for the nmrglue Project
    Copyright (c) 2010-2015 Jonathan J. Helmus
    All rights reserved.
    """

    def sol_general(self, data, filter, w=16, mode="same"):
        """
        Solvent filter with generic filter.

        Algorithm described in: Marion et al. JMR 1989 84 425-430

        Parameters
        ----------
        data : 1D or 2D ndarray
            Array of 1D or 2D NMR data.
        filter : ndarray
            Filter to convolve with data.  Not used in solvent filter functions
            which specific the filter, e.g. sol_boxcar.
        w : int, optional
            Filter length.  Not used here but is used in solent filter functions
            which specify the filter, e.g. sol_boxcar.
        mode : {'valid', 'same', 'full'}, optional
            Convolution mode, 'same' should be used.

        Returns
        -------
        ndata : 1D or 2D ndarray
            NMR data with solvent filter applied

        """
        import scipy

        A = filter.sum()
        if data.ndim == 2:
            filter = filter.reshape((1, -1))  # apply along axis=1
        elif data.ndim == 3:
            filter = filter.reshape((1, 1, -1))  # apply along axis=2
        return data - scipy.signal.convolve(data, filter, mode=mode) / A

    """
    Obtained from nmrglue followed by customisation
    
    Copyright Notice and Statement for the nmrglue Project
    Copyright (c) 2010-2015 Jonathan J. Helmus
    All rights reserved.
    """

    def suppress_solvent_3d(self, data, filt, axis=-1, mode="same"):
        """
        Applies 1D solvent suppression filter along one axis of 3D data.

        Parameters:
            data: np.ndarray (3D or higher)
            filt: 1D array-like (your boxcar or other filter)
            axis: axis to apply filter on (e.g., -1 for last)
            mode: convolution mode, usually 'same'
        """
        A = np.sum(filt)
        # Move target axis to last for convenience
        data = np.moveaxis(data, axis, -1)
        from scipy.signal import fftconvolve

        # Apply 1D convolution along last axis
        filtered = fftconvolve(
            data, filt[None, None, :], mode=mode
        )  # shape must match broadcasting
        result = data - filtered / A

        # Move axis back to original position
        result = np.moveaxis(result, -1, axis)
        return result

    def zero_transpose_3d(self, dic, data):
        """
        Transpose NMRPipe-style data from (X, Y, Z) to (Z, Y, X),
        including correct updates to the NMRPipe dictionary.

        Parameters:
            dic (dict): NMRPipe dictionary
            data (ndarray): NMRPipe data, assumed shape (X, Y, Z)

        Returns:
            new_dic (dict): Transposed dictionary
            new_data (ndarray): Transposed data, shape (Z, Y, X)
        """
        # Transpose data from (X, Y, Z) to (Z, Y, X)
        new_data = np.transpose(data, axes=(2, 1, 0))

        fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
        fn3 = "FDF" + str(int(dic["FDDIMORDER"][2]))  # F1, F2, etc

        # Create new dictionary
        new_dic = dic.copy()

        # for i, new_i in enumerate(
        #     new_axis_order
        # ):  # i = new dim index, new_i = old dim index
        #     for key in dic:
        #         if key.startswith(axis_keys[new_i]):
        #             # e.g., FDF1SW -> FDF1SW, becomes FDF3SW when i == 0 (Z)
        #             suffix = key[len(axis_keys[new_i]) :]  # e.g. SW, ORIG
        #             new_key = axis_keys[i] + suffix  # FDF1SW, FDF2SW, etc.
        #             new_dic[new_key] = dic[key]

        # swapping the FDDIMORDER1 and FDDIMORDER3 values
        order1 = dic["FDDIMORDER1"]
        order3 = dic["FDDIMORDER3"]
        new_dic["FDDIMORDER1"] = order3
        new_dic["FDDIMORDER3"] = order1
        new_dic["FDDIMORDER"][0] = order3
        new_dic["FDDIMORDER"][2] = order1

        new_dic["FDDIMORDER"] = [
            new_dic["FDDIMORDER1"],
            new_dic["FDDIMORDER2"],
            new_dic["FDDIMORDER3"],
            new_dic["FDDIMORDER4"],
        ]

        new_dic["FDSLICECOUNT"] = new_data.shape[-2]
        new_dic["FDSPECNUM"] = new_dic["FDSLICECOUNT"]
        new_dic["FDSIZE"] = new_data.shape[-1]

        if dic[fn3 + "QUADFLAG"] != 1:
            # unpack complex as needed
            new_data = np.array(ng.proc_base.c2ri(new_data), dtype="complex64")
            new_dic[fn3 + "SIZE"] = int(new_dic[fn3 + "SIZE"] / 2)

        return new_dic, new_data

    """
    Obtained from nmrglue followed by customisation
    
    Copyright Notice and Statement for the nmrglue Project
    Copyright (c) 2010-2015 Jonathan J. Helmus
    All rights reserved.
    """

    def transpose_3d(
        self, dic, data, hyper=False, nohyper=False, auto=False, nohdr=False
    ):
        """
        Transpose data (2D).

        Parameters
        ----------
        dic : dict
            Dictionary of NMRPipe parameters.
        data : ndarray
            Array of NMR data.
        hyper : bool
            True to perform hypercomplex transpose.
        nohyper : bool
            True to suppress hypercomplex transpose.
        auto : bool
            True to choose transpose mode automatically.
        nohdr : bool
            True to not update the transpose parameters in ndic.

        Returns
        -------
        ndic : dict
            Dictionary of updated NMRPipe parameters.
        ndata : ndarray
            Array of NMR data which has been transposed.

        """
        # XXX test if works with TPPI
        if nohyper:
            hyper = False

        fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
        fn2 = "FDF" + str(int(dic["FDDIMORDER"][1]))  # F1, F2, etc

        if auto:
            if (dic[fn + "QUADFLAG"] != 1) and (dic[fn2 + "QUADFLAG"] != 1):
                hyper = True
            else:
                hyper = False

        if hyper:  # Hypercomplex transpose need type recast
            data = np.array(ng.proc_base.tp_hyper(data), dtype="complex64")
        else:
            data = np.transpose(data, axes=(0, 2, 1))
            if dic[fn2 + "QUADFLAG"] != 1 and nohyper is False:
                # unpack complex as needed
                data = np.array(ng.proc_base.c2ri(data), dtype="complex64")

        # update the dimensionality and order
        dic["FDSLICECOUNT"] = data.shape[-2]
        if (data.dtype == "float32") and (nohyper is True):
            # when nohyper is True and the new last dimension was complex
            # prior to transposing then FDSIZE is set as if the dimension was
            # converted to complex data, that is half the actual size.
            dic["FDSIZE"] = data.shape[-1] / 2
        else:
            dic["FDSIZE"] = data.shape[-1]

        dic["FDSPECNUM"] = dic["FDSLICECOUNT"]
        dic["FDDIMORDER1"], dic["FDDIMORDER2"] = (
            dic["FDDIMORDER2"],
            dic["FDDIMORDER1"],
        )
        dic["FDDIMORDER"] = [
            dic["FDDIMORDER1"],
            dic["FDDIMORDER2"],
            dic["FDDIMORDER3"],
            dic["FDDIMORDER4"],
        ]

        if dic["FD2DPHASE"] == 0:
            dic["FDF1QUADFLAG"], dic["FDF2QUADFLAG"] = (
                dic["FDF2QUADFLAG"],
                dic["FDF1QUADFLAG"],
            )

        if nohdr is not True:
            dic["FDTRANSPOSED"] = (dic["FDTRANSPOSED"] + 1) % 2

        dic = ng.pipe_proc.clean_minmax(dic)
        return dic, data

    """
    Obtained from nmrglue followed by customisation
    
    Copyright Notice and Statement for the nmrglue Project
    Copyright (c) 2010-2015 Jonathan J. Helmus
    All rights reserved.
    """

    def ext(self, dic, data, x1, xn, sw):
        """
        Extract a region. Adapted from nmrglue

        Parameters
        ----------
        dic : dict
            Dictionary of NMRPipe parameters.
        data : ndarray
            Array of NMR data.
        x1 : int or 'default'
            Starting point of the X-axis extraction. 'default' will start the
            extraction at the first point.
        xn : int or 'default'
            Ending point of the X-axis extraction. 'default' will stop the
            extraction at the last point.
        sw : bool
            True to update the sweep width and ppm calibration parameters,
            recommended.

        Returns
        -------
        ndic : dict
            Dictionary of updated NMRPipe parameters.
        ndata : ndarray
            Extracted region of NMR data.

        """

        # store old sizes
        old_x = float(data.shape[-1])

        # slice find limits
        if x1 == "default":
            x_min = 0
        else:
            x_min = np.round(x1) - 1

        if xn == "default":
            x_max = data.shape[-1]
        else:
            x_max = np.round(xn)

        r_x = 1
        fn = "FDF" + str(int(dic["FDDIMORDER"][0]))

        # round size to be multiple of r_x when axis is cut
        if x1 != "default" or xn != "default":
            remain_x = (x_min - x_max) % r_x  # -len_x%r_x
            x_min = x_min - np.floor(remain_x / 2)
            x_max = x_max + remain_x - np.floor(remain_x / 2)

        if x_min < 0:
            x_max = x_max - x_min
            x_min = 0.0

        if x_max > data.shape[-1]:
            x_min = x_min - (x_max - data.shape[-1])
            x_max = data.shape[-1]

        no_of_dimensions = len(data.shape)
        if no_of_dimensions == 1:  # 1D Array
            data = data[int(x_min) : int(x_max)]
            dic["FDSIZE"] = x_max - x_min
            dic[fn + "SIZE"] = x_max - x_min
            dic[fn + "FTSIZE"] = x_max - x_min
            dic[fn + "TDSIZE"] = x_max - x_min
            dic[fn + "APODSIZE"] = x_max - x_min
        else:
            data = data[..., int(x_min) : int(x_max)]
            dic["FDSIZE"] = x_max - x_min
            dic[fn + "SIZE"] = x_max - x_min
            dic[fn + "FTSIZE"] = x_max - x_min
            dic[fn + "TDSIZE"] = x_max - x_min
            dic[fn + "APODSIZE"] = x_max - x_min

        # adjust sweep width and ppm calibration
        if sw:
            fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
            s = data.shape[-1]

            if dic[fn + "FTFLAG"] == 0:  # time domain
                dic[fn + "CENTER"] = float(int(s / 2.0 + 1))
                dic[fn + "APOD"] = s
                dic[fn + "TDSIZE"] = s
                dic = ng.pipe_proc.recalc_orig(dic, data, fn)
            else:  # freq domain
                dic[fn + "X1"] = x_min + 1
                dic[fn + "XN"] = x_max
                dic[fn + "APOD"] = np.floor(dic[fn + "APOD"] * s / old_x)
                dic[fn + "CENTER"] = dic[fn + "CENTER"] - x_min
                dic[fn + "SW"] = dic[fn + "SW"] * s / old_x
                dic = ng.pipe_proc.recalc_orig(dic, data, fn)

        dic = ng.pipe_proc.update_minmax(dic, data)
        return dic, data

    """
    Obtained from nmrglue followed by customisation
    
    Copyright Notice and Statement for the nmrglue Project
    Copyright (c) 2010-2015 Jonathan J. Helmus
    All rights reserved.
    """

    def base(self, dic, data, nl=None, nw=0):
        """
        Linear baseline correction.

        Parameters
        ----------
        dic : dict
            Dictionary of NMRPipe parameters.
        data : ndarray
            Array of NMR data.
        nl : list
            List of baseline node points.
        nw : int
            Node width in points.

        Returns
        -------
        ndic : dict
            Dictionary of updated NMRPipe parameters.
        ndata : ndarray
            Array of NMR data with a linear baseline correction applied.

        """

        # change values in node list to start at 0
        nl = [i - 1 for i in nl]

        data = self.base2(data, nl, nw)
        dic = ng.pipe_proc.update_minmax(dic, data)
        return dic, data

    """
    Obtained from nmrglue followed by customisation
    
    Copyright Notice and Statement for the nmrglue Project
    Copyright (c) 2010-2015 Jonathan J. Helmus
    All rights reserved.
    """

    def base2(self, data, nl, nw=0):
        """
        Linear (first-order) baseline correction based on node list.

        Parameters
        ----------
        data : 1D or 2D ndarray
            Array of 1D or 2D NMR data.
        nl : list
            List of baseline nodes.
        nw : float, optional
            Node half-width in points.

        Returns
        -------
        ndata : ndarray
            NMR data with first order baseline correction applied.  For 2D data
            baseline correction is applied for each trace along the last
            dimension.

        """
        if data.ndim == 1:
            data = data - ng.proc_bl.calc_bl_linear(data, nl, nw)
        elif data.ndim == 2:  # for 2D array loop over traces
            for i, vec in enumerate(data):
                data[i] = data[i] - ng.proc_bl.calc_bl_linear(vec, nl, nw)

        else:
            for i, vec in enumerate(data):
                for j, vec2 in enumerate(vec):
                    data[i][j] = data[i][j] - ng.proc_bl.calc_bl_linear(vec2, nl, nw)

        return data

    """
    Obtained from nmrglue followed by customisation
    
    Copyright Notice and Statement for the nmrglue Project
    Copyright (c) 2010-2015 Jonathan J. Helmus
    All rights reserved.
    """

    def lp(
        self,
        dic,
        data,
        pred="default",
        x1="default",
        xn="default",
        ord=8,
        mode="f",
        append="after",
        bad_roots="auto",
        mirror=None,
        fix_mode="on",
        method="tls",
    ):
        """
        Linear Prediction

        Parameters
        ----------
        dic : dict
            Dictionary of NMRPipe parameters.
        data : ndarray
            Array of NMR data.
        pred : int
            Number of points to predict, "default" chooses the vector size for
            forward prediction, 1 for backward prediction
        x1 : int or 'default'
            First point in 1D vector to use to extract LP filter. 'default' will
            use the first or last point depending on the mode.
        xn : int or 'default'
            Last point in 1D vector to use to extract LP filter. 'default' will use
            the first or last point depending on the mode.
        ord : int
            Prediction order, number of LP coefficients used in prediction.
        mode : {'f', 'b', 'fb'}
            Mode to generate LP filter, 'f' for forward,'b' for backward, 'fb' for
            forward-backward.
        append : {'before' or 'after'}
            Location to append predicted data, 'before' or 'after' the existing
            data.
        bad_roots {'incr', 'decr', None, 'auto'} :
            Type of roots which are will be marked as bad and stabilized. Choices
            are 'incr' for increasing roots, 'decr' for decreasing roots, or None
            for not root stabilization. The default 'auto' will set this parameter
            based upon the LP `mode` parameter: 'f' and 'fb' will results in an
            'incr' parameter. 'b' in 'decr'.
        mirror : {'90-180', '0-0', None}
            Mirror mode, option are '90-180' for a one point shifted mirror image,
            '0-0' for an exact mirror image, and None for no mirror imaging of the
            data.
        fix_mode : {'on', 'reflect'}
            Method used to stabilize bad roots, 'on' moves bad roots onto the unit
            circle, 'reflect' reflect bad roots across the unit circle.
        method : {'svd', 'qr', 'choleskey', 'tls'}
            Method to use to calculate the LP filter.

        Notes
        -----
        The results from this function do not match NMRPipe's LP function.  Also
        some additional parameter and different parameter in this function.

        Returns
        -------
        ndic : dict
            Dictionary of updated NMRPipe parameters.
        ndata : ndarray
            Array of NMR data with linear prediction applied.

        """
        # check parameter
        if mirror not in [None, "90-180", "0-0"]:
            raise ValueError("mirror must be None, '90-180' or '0-0'")

        # pred default values
        if pred == "default":
            if mode == "f" or mode == "fb":
                pred = data.shape[-1]  # double the number of points
            else:
                pred = 1  # predict 1 point before the data

        # remove first pred points if appending before data
        if append == "before":
            data = data[..., pred:]

        # create slice object
        if x1 == "default":
            x_min = 0
        elif mode == "before":
            x_min = x1 - pred - 1
        else:
            x_min = x1 - 1

        if xn == "default":
            x_max = data.shape[-1]
        else:
            x_max = xn - 1
        sl = slice(x_min, x_max)

        # mirror mode (remap to proc_lp names
        mirror = {None: None, "90-180": "180", "0-0": "0"}[mirror]

        # mode, append, bad_roots, fix_mode, and method are passed unchanged
        # use LP-TLS for best results
        data = self.lp2(
            data, pred, sl, ord, mode, append, bad_roots, fix_mode, mirror, method
        )

        # calculation for dictionary updates
        fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc
        s = data.shape[-1]
        s2 = s / 2.0 + 1

        # update the dictionary
        dic[fn + "CENTER"] = s2
        if dic["FD2DPHASE"] == 1 and fn != "FDF2":  # TPPI data
            dic[fn + "CENTER"] = np.round(s2 / 2.0 + 0.001)
        dic = ng.pipe_proc.recalc_orig(dic, data, fn)
        dic["FDSIZE"] = s
        dic[fn + "SIZE"] = s
        dic[fn + "APOD"] = s
        dic[fn + "TDSIZE"] = s

        dic = ng.pipe_proc.update_minmax(dic, data)
        return dic, data

    """
    Obtained from nmrglue followed by customisation
    
    Copyright Notice and Statement for the nmrglue Project
    Copyright (c) 2010-2015 Jonathan J. Helmus
    All rights reserved.
    """

    def lp2(
        self,
        data,
        pred=1,
        slice=slice(None),
        order=8,
        mode="f",
        append="after",
        bad_roots="auto",
        fix_mode="on",
        mirror=None,
        method="svd",
    ):
        """
        Linear prediction extrapolation of 1D or 2D data.

        Parameters
        ----------
        data : ndarray
            1D or 2D NMR data with the last (-1) axis in the time domain.
        pred : int
            Number of points to predict along the last axis.
        slice : slice object, optional
            Slice object which selects the region along the last axis to use in LP
            equation.  The default (slice(None)) will use all points.
        order : int
            Prediction order, number of LP coefficients calculated.
        mode : {'f', 'b', 'fb' or 'bf'}
            Mode to generate LP filter. 'f' for forward,'b' for backward, fb for
            'forward-backward and 'bf' for backward-forward.
        append : {'before', 'after'}
            Location to append the data, either 'before' the current data, or
            'after' the existing data. This is independent of the `mode` parameter.
        bad_roots : {'incr', 'decr', None, 'auto'}
            Type of roots which to consider bad and to stabilize.  Option are those
            with increasing signals 'incr' or decreasing signals 'decr'.  None will
            perform no root stabilizing.  The default ('auto') will set the
            parameter based on the `mode` parameter.  'f' or 'fb' `mode` will
            results in a 'incr' `bad_roots` parameter, 'b' or 'bf` in 'decr'
        fix_mode : {'on', 'reflect'}
            Method used to stabilize bad roots, 'on' to move the roots onto the
            unit circle, 'reflect' to reflect bad roots across the unit circle.
            This parameter is ignored when `bad_roots` is None.
        mirror : {None, '0', '180'}
            Mode to form mirror image of data before processing.  None will
            process the data trace as provided (no mirror image). '0' or '180'
            forms a mirror image of the sliced trace to calculate the LP filter.
            '0' should be used with data with no delay, '180' with data
            with an initial half-point delay.
        method : {'svd', 'qr', 'choleskey', 'tls'}
            Method to use to calculate the LP filter. Choices are a SVD ('svd'), QR
            ('qr'), or Choleskey ('choleskey') decomposition, or Total Least
            Squares ('tls').

        Returns
        -------
        ndata : ndarray
            NMR data with `pred` number of points linear predicted and appended to
            the original data.

        Notes
        -----
        When given 2D data a series of 1D linear predictions are made to
        each row in the array, extending each by pred points. To perform a 2D
        linear prediction using a 2D prediction matrix use :py:func:`lp2d`.

        In forward-backward or backward-forward mode root stabilizing is done
        on both sets of signal roots as calculated in the first mode direction.
        After averaging the coefficient the roots are again stabilized.

        When the append parameter does not match the LP mode, for example
        if a backward linear prediction (mode='b') is used to predict points
        after the trace (append='after'), any root fixing is done before reversing
        the filter.

        """
        if data.ndim == 1:
            return ng.proc_lp.lp_1d(
                data,
                pred,
                slice,
                order,
                mode,
                append,
                bad_roots,
                fix_mode,
                mirror,
                method,
            )
        elif data.ndim == 2:
            # create empty array to hold output
            s = list(data.shape)
            s[-1] = s[-1] + pred
            new = np.empty(s, dtype=data.dtype)
            # vector-wise 1D LP
            for i, trace in enumerate(data):
                new[i] = ng.proc_lp.lp_1d(
                    trace,
                    pred,
                    slice,
                    order,
                    mode,
                    append,
                    bad_roots,
                    fix_mode,
                    mirror,
                    method,
                )
            return new
        else:
            # create empty array to hold output
            s = list(data.shape)
            s[-1] = s[-1] + pred
            new = np.empty(s, dtype=data.dtype)
            # vector-wise 1D LP
            for i, trace in enumerate(data):
                for j, trace2 in enumerate(trace):
                    new[i][j] = ng.proc_lp.lp_1d(
                        trace2,
                        pred,
                        slice,
                        order,
                        mode,
                        append,
                        bad_roots,
                        fix_mode,
                        mirror,
                        method,
                    )
            return new

    """
    Obtained from nmrglue followed by customisation
    
    Copyright Notice and Statement for the nmrglue Project
    Copyright (c) 2010-2015 Jonathan J. Helmus
    All rights reserved.
    """

    def zf(
        self,
        dic,
        data,
        zf=1,
        pad="auto",
        size="auto",
        mid=False,
        inter=False,
        auto=False,
        inv=False,
    ):
        """
        Zero fill

        Parameters
        ----------
        dic : dict
            Dictionary of NMRPipe parameters.
        data : ndarray
            Array of NMR data.
        zf : int, optional.
            Number of times to double the current dimensions size.
        pad : int
            Number of zeros to pad the data with.
        size : int
            Desired final size of the current dimension.
        mid : bool
            True to zero fill in the middle of the current dimension
        inter : bool
            True to zero fill between points.
        auto : bool
            True to round final size to nearest power of two.
        inv : bool
            True to extract the time domain data (remove zero filling).

        Returns
        -------
        ndic : dict
            Dictionary of updated NMRPipe parameters.
        ndata : ndarray
            Array of NMR data which has zero filled.

        Notes
        -----
        Only one of the `zf`, `pad` and `size` parameter should be used, the other
        should be left as the default value.  If any of the `mid`, `inter`, `auto`
        and `inv` parameters are True other parameter may be ignored.

        """
        fn = "FDF" + str(int(dic["FDDIMORDER"][0]))  # F1, F2, etc

        if inv:  # recover original time domain points
            # calculation for dictionary updates
            s = dic[fn + "TDSIZE"]
            s2 = s / 2.0 + 1

            # update the dictionary
            dic[fn + "ZF"] = -1.0 * s
            dic[fn + "CENTER"] = s2
            dic = ng.pipe_proc.recalc_orig(dic, data, fn)
            dic["FDSIZE"] = s
            return dic, data[..., : int(s)]

        if inter:  # zero filling between points done first
            data = ng.proc_base.zf_inter(data, zf)
            dic[fn + "SW"] = dic[fn + "SW"] * (zf + 1)
            zf = 0
            pad = 0  # NMRPipe ignores pad after a inter zf

        # set zpad, the number of zeros to be padded
        zpad = data.shape[-1] * 2**zf - data.shape[-1]

        if pad != "auto":
            zpad = pad
        if size != "auto":
            zpad = size - data.shape[-1]

        # auto is applied on top of other parameters:
        if auto:
            fsize = data.shape[-1] + zpad
            fsize = 2 ** (np.ceil(np.log(fsize) / np.log(2)))
            zpad = fsize - data.shape[-1]

        if zpad < 0:
            zpad = 0

        data = ng.proc_base.zf_pad(data, pad=zpad, mid=mid)

        # calculation for dictionary updates
        s = data.shape[-1]
        s2 = s / 2.0 + 1

        # update the dictionary
        dic[fn + "ZF"] = -1.0 * s
        dic[fn + "SIZE"] = s
        dic[fn + "TDSIZE"] = s
        dic[fn + "APODSIZE"] = s
        dic[fn + "CENTER"] = s2
        if dic["FD2DPHASE"] == 1 and fn != "FDF2":  # TPPI data
            dic[fn + "CENTER"] = np.round(s2 / 2.0 + 0.001)
        dic = ng.pipe_proc.recalc_orig(dic, data, fn)
        dic["FDSIZE"] = s
        dic = ng.pipe_proc.update_minmax(dic, data)
        return dic, data
