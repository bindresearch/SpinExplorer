#!/usr/bin/env python3

"""MIT License

Copyright (c) 2025 James Eaton, Andrew Baldwin (University of Oxford)

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
import numpy as np
import nmrglue as ng
import os
import json


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

        self.on_run_processing_nmrglue()

    def on_run_processing_nmrglue(self):
        """
        Apply the processing parameters to the data
        """
        self.apply_processing_parameters()

    def apply_processing_parameters(self):
        # Process the data according to the user inputted processing parameters

        # Initial FID data
        dic, data = self.nmr_data.dic, self.nmr_data.data

        # Checking whether processing is required for second/third dimensions
        include_dim2, include_dim3 = self.checking_dimensions()

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
                if self.nmr_data.pseudo_axis == False:
                    # data = data.T
                    # dic = self.update_dictionary_from_transpose(dic)
                    dic, data = ng.pipe_proc.tp(dic, data)
                elif self.nmr_data.pseudo_axis == True:
                    if self.nmr_data.index == 2:
                        dic, data = ng.pipe_proc.tp(dic, data)
                    elif self.nmr_data.index == 1:
                        # If the pseudo axis is the central axis then need to move the third axis
                        dic, data = ng.pipe_proc.tp(dic, data)
                        dic, data = self.zero_transpose_3d(dic, data)
            dic, data = self.apply_dimension_processing(
                dic, data, 1, self.dimension_tabs[1]
            )
            # dic, data = ng.pipe_proc.tp(dic, data)

        if include_dim3:
            dic, data = self.zero_transpose_3d(dic, data)
            dic, data = self.apply_dimension_processing(
                dic, data, 2, self.dimension_tabs[2]
            )

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
        dic["FDF3QUADFLAG"] = 1.0
        dic["FDQUADFLAG"] = 1.0

        data = data.real
        data = data.astype(np.float32)

        ng.pipe.write(nmrfile, dic, data, overwrite=True)

    def update_dictionary_from_transpose(self, dic):

        import copy

        new_dic = copy.deepcopy(dic)

        # Step 1: Swap FDF1* and FDF2* keys
        for key in list(dic.keys()):
            if key.startswith("FDF1"):
                suffix = key[4:]
                f1_key = f"FDF1{suffix}"
                f2_key = f"FDF2{suffix}"
                if f2_key in dic:
                    new_dic[f1_key], new_dic[f2_key] = dic[f2_key], dic[f1_key]

        # Step 2: Swap entries in FDDIMORDER
        if "FDDIMORDER" in dic:
            new_dic["FDDIMORDER"] = [
                2.0 if x == 1.0 else 1.0 if x == 2.0 else x for x in dic["FDDIMORDER"]
            ]

        # Step 3: Swap individual FDDIMORDERn values
        for i in range(1, 5):
            key = f"FDDIMORDER{i}"
            if key in dic:
                val = dic[key]
                if val == 1.0:
                    new_dic[key] = 2.0
                elif val == 2.0:
                    new_dic[key] = 1.0

        # Step 4: Set transpose flag
        new_dic["FDTRANSPOSED"] = 1.0

        return new_dic

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
        self.notebook.parent.original_frame.parent.Close()
        if self.notebook.parent.original_frame.parent.path != "":
            os.chdir(self.notebook.parent.original_frame.parent.path)
        from SpinExplorer.SpinView.SpinView import MyApp

        app = MyApp()
        if self.notebook.parent.original_frame.parent.cwd != "":
            app.path = path
            app.cwd = cwd

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
        else:
            selection = tab.linear_prediction_radio_box_indirect.GetSelection()
            if selection != 1:
                return dic, data

        # Apply linear prediction

        if tab.linear_prediction_options_selection == 0:
            append = "after"
        else:
            append = "before"
        if tab.linear_prediction_coefficients_selection == 0:
            mode = "f"
        elif tab.linear_prediction_coefficients_selection == 1:
            mode = "b"
        else:
            mode = "fb"
        dic, data = ng.pipe_proc.lp(dic, data, pred="default", mode=mode, append=append)

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
                dic, data = ng.pipe_proc.zf(
                    dic,
                    data,
                    zf=int(tab.zero_filling_value_doubling_times),
                    auto=round,
                )
            elif tab.zero_filling_combobox_selection == 1:
                dic, data = ng.pipe_proc.zf(
                    dic,
                    data,
                    pad=int(tab.zero_filling_value_zeros_to_add),
                    auto=round,
                )
            elif tab.zero_filling_combobox_selection == 2:
                dic, data = ng.pipe_proc.zf(
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
                dic, data = ng.pipe_proc.ft(dic, data, auto=True)
            elif tab.ft_method_selection == 1:
                dic, data = ng.pipe_proc.ft(dic, data, real=True)
            elif tab.ft_method_selection == 2:
                dic, data = ng.pipe_proc.ft(dic, data, inv=True)
            elif tab.ft_method_selection == 3:
                dic, data = ng.pipe_proc.ft(dic, data, alt=True)

        if dimension == 0:
            digital_filter_removal = self.check_digital_filter_removal()
            if digital_filter_removal == True:
                dic_bruker, dat_bruker = ng.bruker.read("./")
                data = self.remove_digital_filter(dic_bruker, data)

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
            ppm_values = ng.pipe.make_uc(dic, data, dim=dimension)
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
            dic, data = ng.pipe_proc.ext(dic, data, x1=x_initial, xn=x_final, sw=True)

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
                message = "The selected baseline correction method is not supported for nmrglue processing. Continuing without baselining. Please use a machine containing nmrPipe or use a linear baselining method for polynomial baselining."
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
            node_list_final = (node_list_final / 100) * len(data)
            node_list_final = node_list_final.astype(int)
            # Replace any zeros with a number greater than 1 to allow the nmrglue baselining routines to work correctly
            node_list_final[node_list_final == 0] = (
                int(tab.baseline_correction_nodes_textcontrol.GetValue()) + 1
            )

            dic, data = ng.pipe_proc.base(
                dic,
                data,
                nl=node_list_final,
                nw=int(tab.baseline_correction_nodes_textcontrol.GetValue()),
            )

        return dic, data

    """
    Obtained from the nmrglue code nmrglue/nmrglue/fileio/bruker.py for customisation
    
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
    Obtained from the nmrglue code nmrglue/nmrglue/fileio/bruker.py for customisation
    
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

    # def process_dimension_1(self, dic, data, dim):
    #     # Process the first dimension
    #     if self.tabDim1.solvent_suppression_checkbox.GetValue() == True:
    #         # Apply solvent suppression
    #         data_orgiginal = data
    #         if self.tabDim1.solvent_suppression_filter_selection == 0:
    #             filter_size = int(
    #                 self.tabDim1.solvent_suppression_filter_length
    #             )  # Larger filter in time domain is larger filter in the frequency domain
    #             if (
    #                 int(self.tabDim1.solvent_suppression_lowpass_shape_selection) + 1
    #                 == 1
    #             ):
    #                 from scipy.signal.windows import boxcar

    #                 filter = boxcar(filter_size)
    #             elif (
    #                 int(self.tabDim1.solvent_suppression_lowpass_shape_selection) + 1
    #                 == 2
    #             ):
    #                 filter = np.cos(np.pi * np.linspace(-0.5, 0.5, filter_size))
    #             else:
    #                 filter = np.cos(np.pi * np.linspace(-0.5, 0.5, filter_size)) ** 2

    #             data = self.sol_general(data, filter, w=filter_size, mode="same")
    #             # data = self.suppress_solvent_3d(data,filter)
    #             # dic,data = ng.pipe_proc.sol(dic,data_orgiginal, mode = 'low', fs=int(self.tabDim1.solvent_suppression_lowpass_shape_selection)+1)
    #             # fl = 16
    #             # fs = int(self.tabDim1.solvent_suppression_lowpass_shape_selection)+1
    #             # import scipy.signal.windows
    #             # if(fs == 1):
    #             #     filter = scipy.signal.windows.boxcar(33)
    #             #     data =  self.sol_general_nd(data, filter)

    #             # try:
    #             #     dic,data = ng.pipe_proc.sol(dic,data_orgiginal, mode = 'low', fs=int(self.tabDim1.solvent_suppression_lowpass_shape_selection)+1)
    #             # except:
    #             #     try:
    #             #         fl = 16
    #             #         fs = int(self.tabDim1.solvent_suppression_lowpass_shape_selection)+1
    #             #         import scipy.signal.windows
    #             #         if(fs == 1):
    #             #             filter = scipy.signal.windows.boxcar(33)
    #             #             data =  ng.proc_bl.sol_general(data, filter, w=33, mode='same')
    #             #         else:
    #             #             # Give a message saying the solvent suppression did not work correctly, the spectrum was processed without solvent suppression
    #             #             dlg = wx.MessageDialog(self, 'The solvent suppression filter did not work correctly. Continuing without digital solvent suppression.', 'Warning', wx.OK | wx.ICON_WARNING)
    #             #             self.Raise()
    #             #             self.SetFocus()
    #             #             result = dlg.ShowModal()
    #             #             data = data_orgiginal
    #             #     except:
    #             #         # Give a message saying the solvent suppression did not work correctly, the spectrum was processed without solvent suppression
    #             #         dlg = wx.MessageDialog(self, 'The solvent suppression filter did not work correctly. Continuing without digital solvent suppression.', 'Warning', wx.OK | wx.ICON_WARNING)
    #             #         self.Raise()
    #             #         self.SetFocus()
    #             #         result = dlg.ShowModal()
    #             #         data = data_orgiginal

    #         else:
    #             # Give an error saying that the selected solvent suppression filter is not supported for windows processing, please use a machine containing nmrPipe
    #             dlg = wx.MessageDialog(
    #                 self,
    #                 "The selected solvent suppression filter is not supported for windows processing. Please change to low bandpass filter or use a machine containing nmrPipe.",
    #                 "Warning",
    #                 wx.OK | wx.ICON_WARNING,
    #             )
    #             self.Raise()
    #             self.SetFocus()
    #             result = dlg.ShowModal()
    #             return
    #     if self.tabDim1.linear_prediction_checkbox.GetValue() == True:
    #         # Apply linear prediction
    #         if self.tabDim1.linear_prediction_options_selection == 0:
    #             if self.tabDim1.linear_prediction_options_selection == 0:
    #                 append = "after"
    #             else:
    #                 append = "before"
    #             if self.tabDim1.linear_prediction_coefficients_selection == 0:
    #                 mode = "f"
    #             elif self.tabDim1.linear_prediction_coefficients_selection == 1:
    #                 mode = "b"
    #             else:
    #                 mode = "fb"
    #             dic, data = ng.pipe_proc.lp(
    #                 dic, data, pred="default", mode=mode, append=append
    #             )

    #     if self.tabDim1.apodization_checkbox.GetValue() == True:
    #         # Apply apodization
    #         if self.tabDim1.apodization_combobox_selection == 0:
    #             dic, data = ng.pipe_proc.em(
    #                 dic,
    #                 data,
    #                 lb=0.0,
    #                 c=float(self.tabDim1.apodization_first_point_scaling),
    #             )
    #         elif self.tabDim1.apodization_combobox_selection == 1:
    #             dic, data = ng.pipe_proc.em(
    #                 dic,
    #                 data,
    #                 lb=float(self.tabDim1.exponential_line_broadening),
    #                 c=float(self.tabDim1.apodization_first_point_scaling),
    #             )
    #         elif self.tabDim1.apodization_combobox_selection == 2:
    #             dic, data = ng.pipe_proc.gm(
    #                 dic,
    #                 data,
    #                 g1=float(self.tabDim1.g1),
    #                 g2=float(self.tabDim1.g2),
    #                 g3=float(self.tabDim1.g3),
    #                 c=float(self.tabDim1.apodization_first_point_scaling),
    #             )
    #         elif self.tabDim1.apodization_combobox_selection == 3:
    #             dic, data = ng.pipe_proc.sp(
    #                 dic,
    #                 data,
    #                 off=float(self.tabDim1.offset),
    #                 end=float(self.tabDim1.end),
    #                 pow=int(self.tabDim1.power),
    #                 c=float(self.tabDim1.apodization_first_point_scaling),
    #             )
    #         elif self.tabDim1.apodization_combobox_selection == 4:
    #             dic, data = ng.pipe_proc.gmb(
    #                 dic,
    #                 data,
    #                 lb=float(self.tabDim1.a),
    #                 gb=float(self.tabDim1.b),
    #                 c=float(self.tabDim1.apodization_first_point_scaling),
    #             )
    #         elif self.tabDim1.apodization_combobox_selection == 5:
    #             dic, data = ng.pipe_proc.tp(
    #                 dic,
    #                 data,
    #                 t1=float(self.tabDim1.t1),
    #                 t2=float(self.tabDim1.t2),
    #                 c=float(self.tabDim1.apodization_first_point_scaling),
    #             )
    #         elif self.tabDim1.apodization_combobox_selection == 6:
    #             dic, data = ng.pipe_proc.tri(
    #                 dic,
    #                 data,
    #                 loc=float(self.tabDim1.loc),
    #                 c=float(self.tabDim1.apodization_first_point_scaling),
    #             )

    #     if self.tabDim1.zero_filling_checkbox.GetValue() == True:
    #         if self.tabDim1.zero_filling_round_checkbox.GetValue() == True:
    #             round = True
    #         else:
    #             round = False
    #         if self.tabDim1.zero_filling_combobox_selection == 0:
    #             dic, data = ng.pipe_proc.zf(
    #                 dic,
    #                 data,
    #                 zf=int(self.tabDim1.zero_filling_value_doubling_times),
    #                 auto=round,
    #             )
    #         elif self.tabDim1.zero_filling_combobox_selection == 1:
    #             dic, data = ng.pipe_proc.zf(
    #                 dic,
    #                 data,
    #                 pad=int(self.tabDim1.zero_filling_value_zeros_to_add),
    #                 auto=round,
    #             )
    #         elif self.tabDim1.zero_filling_combobox_selection == 2:
    #             dic, data = ng.pipe_proc.zf(
    #                 dic,
    #                 data,
    #                 size=int(self.tabDim1.zero_filling_value_final_data_size),
    #                 auto=round,
    #             )

    #     if self.tabDim1.fourier_transform_checkbox.GetValue() == True:
    #         if self.tabDim1.ft_method_selection == 0:
    #             dic, data = ng.pipe_proc.ft(dic, data, auto=True)
    #         elif self.tabDim1.ft_method_selection == 1:
    #             dic, data = ng.pipe_proc.ft(dic, data, real=True)
    #         elif self.tabDim1.ft_method_selection == 2:
    #             dic, data = ng.pipe_proc.ft(dic, data, inv=True)
    #         elif self.tabDim1.ft_method_selection == 3:
    #             dic, data = ng.pipe_proc.ft(dic, data, alt=True)

    #     dic_bruker, dat_bruker = ng.bruker.read("./")
    #     data = self.remove_digital_filter(dic_bruker, data)

    #     if self.tabDim1.phase_correction_checkbox.GetValue() == True:
    #         dic, data = ng.pipe_proc.ps(
    #             dic,
    #             data,
    #             p0=float(self.tabDim1.phase_correction_p0_textcontrol.GetValue()),
    #             p1=float(self.tabDim1.phase_correction_p1_textcontrol.GetValue()),
    #         )

    #     if self.tabDim1.magnitude_mode_checkbox.GetValue() == True:
    #         dic, data = ng.pipe_proc.mc(dic, data)

    #     if self.tabDim1.extraction_checkbox.GetValue() == True:
    #         # Find the indexes of the ppm values selected
    #         # Get the ppm values from the data
    #         ppm_values = ng.pipe.make_uc(dic, data, dim=dim)
    #         ppm_values = ppm_values.ppm_scale()
    #         x_initial = np.abs(
    #             ppm_values
    #             - float(self.tabDim1.extraction_ppm_start_textcontrol.GetValue())
    #         ).argmin()
    #         x_final = np.abs(
    #             ppm_values
    #             - float(self.tabDim1.extraction_ppm_end_textcontrol.GetValue())
    #         ).argmin()
    #         if x_initial > x_final:
    #             x_initial, x_final = x_final, x_initial
    #         # Change x_initial and x_final so that the difference is an even number
    #         if (x_final - x_initial + 1) % 2 != 0:
    #             x_final += 1
    #         dic, data = ng.pipe_proc.ext(dic, data, x1=x_initial, xn=x_final, sw=True)

    #     if self.tabDim1.baseline_correction_checkbox.GetValue() == True:
    #         if self.tabDim1.baseline_correction_radio_box_selection == 1:
    #             # If POLY baseline correction is selected, this is not currently supported on windows without nmrPipe
    #             message = "The selected baseline correction method is not supported for windows processing. Please use a machine containing nmrPipe or use a linear baselining method."
    #             dlg = wx.MessageDialog(
    #                 self, message, "Warning", wx.OK | wx.ICON_WARNING
    #             )
    #             self.Raise()
    #             self.SetFocus()
    #             result = dlg.ShowModal()
    #             return

    #         # Split the node list
    #         node_list = (
    #             self.tabDim1.baseline_correction_node_list_textcontrol.GetValue()
    #         )
    #         node_list = node_list.split(",")
    #         node_list_final = []
    #         for node in node_list:
    #             node_list_final.append(float(node))

    #         # Convert nodes into points
    #         node_list_final = np.array(node_list_final)
    #         node_list_final = (node_list_final / 100) * len(data)
    #         node_list_final = node_list_final.astype(int)
    #         # Replace any zeros with a number greater than 1 to allow the nmrglue baselining routines to work correctly
    #         node_list_final[node_list_final == 0] = (
    #             int(self.tabDim1.baseline_correction_nodes_textcontrol.GetValue()) + 1
    #         )

    #         dic, data = ng.pipe_proc.base(
    #             dic,
    #             data,
    #             nl=node_list_final,
    #             nw=int(self.tabDim1.baseline_correction_nodes_textcontrol.GetValue()),
    #         )

    #     return dic, data

    # def process_dimension_2(self, dic, data):
    #     # Transpose to the second dimension
    #     dic, data = ng.pipe_proc.tp(dic, data)

    #     # Process the second dimension
    #     if self.tabDim2.linear_prediction_radio_box_dim2_selection == 1:
    #         # Apply linear prediction
    #         if self.tabDim2.linear_prediction_dim2_options_selection == 0:
    #             append = "after"
    #         else:
    #             append = "before"
    #         if self.tabDim2.linear_prediction_dim2_coefficients_selection == 0:
    #             mode = "f"
    #         elif self.tabDim2.linear_prediction_dim2_coefficients_selection == 1:
    #             mode = "b"
    #         else:
    #             mode = "fb"
    #         dic, data = ng.pipe_proc.lp(
    #             dic, data, pred="default", mode=mode, append=append
    #         )

    #     if self.tabDim2.apodization_checkbox_dim2.GetValue() == True:
    #         # Apply apodization
    #         if self.tabDim2.apodization_dim2_combobox_selection == 0:
    #             dic, data = ng.pipe_proc.em(
    #                 dic,
    #                 data,
    #                 lb=0.0,
    #                 c=float(self.tabDim2.apodization_first_point_scaling_dim2),
    #             )
    #         elif self.tabDim2.apodization_dim2_combobox_selection == 1:
    #             dic, data = ng.pipe_proc.em(
    #                 dic,
    #                 data,
    #                 lb=float(self.tabDim2.exponential_line_broadening_dim2),
    #                 c=float(self.tabDim2.apodization_first_point_scaling_dim2),
    #             )
    #         elif self.tabDim2.apodization_dim2_combobox_selection == 2:
    #             dic, data = ng.pipe_proc.gm(
    #                 dic,
    #                 data,
    #                 g1=float(self.tabDim2.g1_dim2),
    #                 g2=float(self.tabDim2.g2_dim2),
    #                 g3=float(self.tabDim2.g3_dim2),
    #                 c=float(self.tabDim2.apodization_first_point_scaling_dim2),
    #             )
    #         elif self.tabDim2.apodization_dim2_combobox_selection == 3:
    #             dic, data = ng.pipe_proc.sp(
    #                 dic,
    #                 data,
    #                 off=float(self.tabDim2.offset_dim2),
    #                 end=float(self.tabDim2.end_dim2),
    #                 pow=int(self.tabDim2.power_dim2),
    #                 c=float(self.tabDim2.apodization_first_point_scaling_dim2),
    #             )
    #         elif self.tabDim2.apodization_dim2_combobox_selection == 4:
    #             dic, data = ng.pipe_proc.gmb(
    #                 dic,
    #                 data,
    #                 lb=float(self.tabDim2.a_dim2),
    #                 gb=float(self.tabDim2.b_dim2),
    #                 c=float(self.tabDim2.apodization_first_point_scaling_dim2),
    #             )
    #         elif self.tabDim2.apodization_dim2_combobox_selection == 5:
    #             dic, data = ng.pipe_proc.tp(
    #                 dic,
    #                 data,
    #                 t1=float(self.tabDim2.t1_dim2),
    #                 t2=float(self.tabDim2.t2_dim2),
    #                 c=float(self.tabDim2.apodization_first_point_scaling_dim2),
    #             )
    #         elif self.tabDim2.apodization_dim2_combobox_selection == 6:
    #             dic, data = ng.pipe_proc.tri(
    #                 dic,
    #                 data,
    #                 loc=float(self.tabDim2.loc_dim2),
    #                 c=float(self.tabDim2.apodization_first_point_scaling_dim2),
    #             )

    #     if self.tabDim2.zero_filling_checkbox_dim2.GetValue() == True:
    #         if self.tabDim2.zero_filling_round_checkbox_dim2.GetValue() == True:
    #             round = True
    #         else:
    #             round = False

    #         if self.tabDim2.zero_filling_dim2_combobox_selection == 0:
    #             dic, data = ng.pipe_proc.zf(
    #                 dic,
    #                 data,
    #                 zf=int(self.tabDim2.zero_filling_dim2_value_doubling_times),
    #                 auto=round,
    #             )
    #         elif self.tabDim2.zero_filling_dim2_combobox_selection == 1:
    #             dic, data = ng.pipe_proc.zf(
    #                 dic,
    #                 data,
    #                 pad=int(self.tabDim2.zero_filling_dim2_value_zeros_to_add),
    #                 auto=round,
    #             )
    #         elif self.tabDim2.zero_filling_dim2_combobox_selection == 2:
    #             dic, data = ng.pipe_proc.zf(
    #                 dic,
    #                 data,
    #                 size=int(self.tabDim2.zero_filling_dim2_value_final_data_size),
    #                 auto=round,
    #             )

    #     if self.tabDim2.fourier_transform_checkbox_dim2.GetValue() == True:
    #         if self.tabDim2.ft_method_selection_dim2 == 0:
    #             dic, data = ng.pipe_proc.ft(dic, data, auto=True)
    #         elif self.tabDim2.ft_method_selection_dim2 == 1:
    #             dic, data = ng.pipe_proc.ft(dic, data, real=True)
    #         elif self.tabDim2.ft_method_selection_dim2 == 2:
    #             dic, data = ng.pipe_proc.ft(dic, data, inv=True)
    #         elif self.tabDim2.ft_method_selection_dim2 == 3:
    #             dic, data = ng.pipe_proc.ft(dic, data, alt=True)

    #     if self.tabDim2.phase_correction_checkbox_dim2.GetValue() == True:
    #         dic, data = ng.pipe_proc.ps(
    #             dic,
    #             data,
    #             p0=float(self.tabDim2.phase_correction_p0_textcontrol_dim2.GetValue()),
    #             p1=float(self.tabDim2.phase_correction_p1_textcontrol_dim2.GetValue()),
    #         )

    #     if self.tabDim2.extraction_checkbox_dim2.GetValue() == True:
    #         # Find the indexes of the ppm values selected
    #         # Get the ppm values from the data
    #         ppm_values = ng.pipe.make_uc(dic, data, dim=1)
    #         ppm_values = ppm_values.ppm_scale()
    #         x_initial = np.abs(
    #             ppm_values
    #             - float(self.tabDim2.extraction_ppm_start_textcontrol_dim2.GetValue())
    #         ).argmin()
    #         x_final = np.abs(
    #             ppm_values
    #             - float(self.tabDim2.extraction_ppm_end_textcontrol_dim2.GetValue())
    #         ).argmin()
    #         if x_initial > x_final:
    #             x_initial, x_final = x_final, x_initial

    #         if (x_final - x_initial + 1) % 2 != 0:
    #             x_final += 1
    #         dic, data = ng.pipe_proc.ext(dic, data, x1=x_initial, xn=x_final, sw=True)

    #     if self.tabDim2.baseline_correction_checkbox_dim2.GetValue() == True:
    #         if self.tabDim2.baseline_correction_radio_box_selection_dim2 == 1:
    #             # If POLY baseline correction is selected, this is not currently supported on windows without nmrPipe
    #             message = "The selected baseline correction method is not supported for windows processing. Please use a machine containing nmrPipe or use a linear baselining method."
    #             dlg = wx.MessageDialog(
    #                 self, message, "Warning", wx.OK | wx.ICON_WARNING
    #             )
    #             self.Raise()
    #             self.SetFocus()
    #             result = dlg.ShowModal()
    #             return

    #         # Split the node list
    #         node_list = (
    #             self.tabDim2.baseline_correction_node_list_textcontrol_dim2.GetValue()
    #         )
    #         node_list = node_list.split(",")
    #         node_list_final = []
    #         for node in node_list:
    #             node_list_final.append(float(node))

    #         # Convert nodes into points
    #         node_list_final = np.array(node_list_final)
    #         node_list_final = (node_list_final / 100) * len(data)
    #         node_list_final = node_list_final.astype(int)
    #         # Replace any zeros with a number greater than 1 to allow the nmrglue baselining routines to work correctly
    #         node_list_final[node_list_final == 0] = (
    #             int(self.tabDim2.baseline_correction_nodes_textcontrol_dim2.GetValue())
    #             + 1
    #         )

    #         dic, data = ng.pipe_proc.base(
    #             dic,
    #             data,
    #             nl=node_list_final,
    #             nw=int(
    #                 self.tabDim2.baseline_correction_nodes_textcontrol_dim2.GetValue()
    #             ),
    #         )

    #     return dic, data

    # def process_dimension_3(self, dic, data):

    #     # Process the second dimension
    #     if self.tabDim3.linear_prediction_radio_box_dim3_selection == 1:
    #         # Apply linear prediction
    #         if self.tabDim3.linear_prediction_dim3_options_selection == 0:
    #             append = "after"
    #         else:
    #             append = "before"
    #         if self.tabDim3.linear_prediction_dim3_coefficients_selection == 0:
    #             mode = "f"
    #         elif self.tabDim3.linear_prediction_dim3_coefficients_selection == 1:
    #             mode = "b"
    #         else:
    #             mode = "fb"
    #         dic, data = ng.pipe_proc.lp(
    #             dic, data, pred="default", mode=mode, append=append
    #         )

    #     if self.tabDim3.apodization_checkbox_dim3.GetValue() == True:
    #         # Apply apodization
    #         if self.tabDim3.apodization_dim3_combobox_selection == 0:
    #             dic, data = ng.pipe_proc.em(
    #                 dic,
    #                 data,
    #                 lb=0.0,
    #                 c=float(self.tabDim3.apodization_first_point_scaling_dim3),
    #             )
    #         elif self.tabDim3.apodization_dim3_combobox_selection == 1:
    #             dic, data = ng.pipe_proc.em(
    #                 dic,
    #                 data,
    #                 lb=float(self.tabDim3.exponential_line_broadening_dim3),
    #                 c=float(self.tabDim3.apodization_first_point_scaling_dim3),
    #             )
    #         elif self.tabDim3.apodization_dim3_combobox_selection == 2:
    #             dic, data = ng.pipe_proc.gm(
    #                 dic,
    #                 data,
    #                 g1=float(self.tabDim3.g1_dim3),
    #                 g2=float(self.tabDim3.g2_dim3),
    #                 g3=float(self.tabDim3.g3_dim3),
    #                 c=float(self.tabDim3.apodization_first_point_scaling_dim3),
    #             )
    #         elif self.tabDim3.apodization_dim3_combobox_selection == 3:
    #             dic, data = ng.pipe_proc.sp(
    #                 dic,
    #                 data,
    #                 off=float(self.tabDim3.offset_dim3),
    #                 end=float(self.tabDim3.end_dim3),
    #                 pow=int(self.tabDim3.power_dim3),
    #                 c=float(self.tabDim3.apodization_first_point_scaling_dim3),
    #             )
    #         elif self.tabDim3.apodization_dim3_combobox_selection == 4:
    #             dic, data = ng.pipe_proc.gmb(
    #                 dic,
    #                 data,
    #                 lb=float(self.tabDim3.a_dim3),
    #                 gb=float(self.tabDim3.b_dim3),
    #                 c=float(self.tabDim3.apodization_first_point_scaling_dim3),
    #             )
    #         elif self.tabDim3.apodization_dim3_combobox_selection == 5:
    #             dic, data = ng.pipe_proc.tp(
    #                 dic,
    #                 data,
    #                 t1=float(self.tabDim3.t1_dim3),
    #                 t2=float(self.tabDim3.t2_dim3),
    #                 c=float(self.tabDim3.apodization_first_point_scaling_dim3),
    #             )
    #         elif self.tabDim3.apodization_dim3_combobox_selection == 6:
    #             dic, data = ng.pipe_proc.tri(
    #                 dic,
    #                 data,
    #                 loc=float(self.tabDim3.loc_dim3),
    #                 c=float(self.tabDim3.apodization_first_point_scaling_dim3),
    #             )

    #     if self.tabDim3.zero_filling_checkbox_dim3.GetValue() == True:
    #         if self.tabDim3.zero_filling_round_checkbox_dim3.GetValue() == True:
    #             round = True
    #         else:
    #             round = False

    #         if self.tabDim3.zero_filling_dim3_combobox_selection == 0:
    #             dic, data = ng.pipe_proc.zf(
    #                 dic,
    #                 data,
    #                 zf=int(self.tabDim3.zero_filling_dim3_value_doubling_times),
    #                 auto=round,
    #             )
    #         elif self.tabDim3.zero_filling_dim3_combobox_selection == 1:
    #             dic, data = ng.pipe_proc.zf(
    #                 dic,
    #                 data,
    #                 pad=int(self.tabDim3.zero_filling_dim3_value_zeros_to_add),
    #                 auto=round,
    #             )
    #         elif self.tabDim3.zero_filling_dim3_combobox_selection == 2:
    #             dic, data = ng.pipe_proc.zf(
    #                 dic,
    #                 data,
    #                 size=int(self.tabDim3.zero_filling_dim3_value_final_data_size),
    #                 auto=round,
    #             )

    #     if self.tabDim3.fourier_transform_checkbox_dim3.GetValue() == True:
    #         if self.tabDim3.ft_method_selection_dim3 == 0:
    #             dic, data = ng.pipe_proc.ft(dic, data, auto=True)
    #         elif self.tabDim3.ft_method_selection_dim3 == 1:
    #             dic, data = ng.pipe_proc.ft(dic, data, real=True)
    #         elif self.tabDim3.ft_method_selection_dim3 == 2:
    #             dic, data = ng.pipe_proc.ft(dic, data, inv=True)
    #         elif self.tabDim3.ft_method_selection_dim3 == 3:
    #             dic, data = ng.pipe_proc.ft(dic, data, alt=True)

    #     if self.tabDim3.phase_correction_checkbox_dim3.GetValue() == True:
    #         dic, data = ng.pipe_proc.ps(
    #             dic,
    #             data,
    #             p0=float(self.tabDim3.phase_correction_p0_textcontrol_dim3.GetValue()),
    #             p1=float(self.tabDim3.phase_correction_p1_textcontrol_dim3.GetValue()),
    #         )

    #     if self.tabDim3.extraction_checkbox_dim3.GetValue() == True:
    #         # Find the indexes of the ppm values selected
    #         # Get the ppm values from the data
    #         ppm_values = ng.pipe.make_uc(dic, data, dim=0)
    #         ppm_values = ppm_values.ppm_scale()
    #         x_initial = np.abs(
    #             ppm_values
    #             - float(self.tabDim3.extraction_ppm_start_textcontrol_dim3.GetValue())
    #         ).argmin()
    #         x_final = np.abs(
    #             ppm_values
    #             - float(self.tabDim3.extraction_ppm_end_textcontrol_dim3.GetValue())
    #         ).argmin()
    #         if x_initial > x_final:
    #             x_initial, x_final = x_final, x_initial

    #         if (x_final - x_initial + 1) % 2 != 0:
    #             x_final += 1
    #         dic, data = ng.pipe_proc.ext(dic, data, x1=x_initial, xn=x_final, sw=True)

    #     if self.tabDim3.baseline_correction_checkbox_dim3.GetValue() == True:
    #         if self.tabDim3.baseline_correction_radio_box_selection_dim3 == 1:
    #             # If POLY baseline correction is selected, this is not currently supported on windows without nmrPipe
    #             message = "The selected baseline correction method is not supported for windows processing. Please use a machine containing nmrPipe or use a linear baselining method."
    #             dlg = wx.MessageDialog(
    #                 self, message, "Warning", wx.OK | wx.ICON_WARNING
    #             )
    #             self.Raise()
    #             self.SetFocus()
    #             result = dlg.ShowModal()
    #             return

    #         # Split the node list
    #         node_list = (
    #             self.tabDim3.baseline_correction_node_list_textcontrol_dim3.GetValue()
    #         )
    #         node_list = node_list.split(",")
    #         node_list_final = []
    #         for node in node_list:
    #             node_list_final.append(float(node))

    #         # Convert nodes into points
    #         node_list_final = np.array(node_list_final)
    #         node_list_final = (node_list_final / 100) * len(data)
    #         node_list_final = node_list_final.astype(int)
    #         # Replace any zeros with a number greater than 1 to allow the nmrglue baselining routines to work correctly
    #         node_list_final[node_list_final == 0] = (
    #             int(self.tabDim3.baseline_correction_nodes_textcontrol_dim3.GetValue())
    #             + 1
    #         )

    #         dic, data = ng.pipe_proc.base(
    #             dic,
    #             data,
    #             nl=node_list_final,
    #             nw=int(
    #                 self.tabDim3.baseline_correction_nodes_textcontrol_dim3.GetValue()
    #             ),
    #         )

    #     # dic, data = self.ztp(dic,data)

    #     # dic['FDF1TDSIZE'] = data.shape[0]
    #     # dic['FDF1FTSIZE'] = data.shape[0]
    #     # dic['FDF1QUADFLAG'] = 1.0

    #     # dic, data = self.ztp(dic,data)
    #     # dic, data = ng.pipe_proc.tp(dic,data)

    #     return dic, data

    # def zero_transpose_3d(self, dic, data):
    #     # Transpose axes 0 and 1 in the 3D array
    #     new_data = data.swapaxes(0, 1)

    #     # Deep copy of the dictionary
    #     import copy

    #     new_dic = copy.deepcopy(dic)

    #     # Swap all FDF1 and FDF2 values
    #     for key in list(dic.keys()):
    #         if key.startswith("FDF1"):
    #             f1_key = key
    #             f2_key = "FDF2" + key[4:]
    #             if f2_key in dic:
    #                 new_dic[f1_key], new_dic[f2_key] = dic[f2_key], dic[f1_key]

    #     # Swap dimension order values
    #     if "FDDIMORDER" in dic:
    #         new_dic["FDDIMORDER"] = [
    #             2.0 if x == 1.0 else 1.0 if x == 2.0 else x for x in dic["FDDIMORDER"]
    #         ]

    #     for i in range(1, 5):
    #         key = f"FDDIMORDER{i}"
    #         if key in dic:
    #             val = dic[key]
    #             if val == 1.0:
    #                 new_dic[key] = 2.0
    #             elif val == 2.0:
    #                 new_dic[key] = 1.0

    #     # Update the FDTRANSPOSED flag
    #     new_dic["FDTRANSPOSED"] = 1.0

    #     return new_dic, new_data

    # def ztp(self, dic, data, nohdr=False):
    #     """
    #     Z-axis transpose (ZTP) for 3D+ NMRPipe data.
    #     Moves the last axis to the front, updates headers.

    #     Parameters
    #     ----------
    #     dic : dict
    #         Dictionary of NMRPipe parameters.
    #     data : ndarray
    #         NMR data array.
    #     nohdr : bool, optional
    #         If True, do not update header metadata.

    #     Returns
    #     -------
    #     dic : dict
    #         Updated NMRPipe dictionary.
    #     data : ndarray
    #         Transposed data.
    #     """
    #     print(dic)
    #     ndim = data.ndim
    #     if ndim < 3:
    #         raise ValueError("ZTP requires at least 3D data.")

    #     # Rotate last axis to front
    #     data = np.transpose(data, axes=(ndim - 1,) + tuple(range(ndim - 1)))

    #     # --- Update FDDIMORDER (dimension order metadata)
    #     # Shift dimensions left, move last to first
    #     dim_keys = [dic.get(f"FDDIMORDER{i+1}", 0) for i in range(3)]
    #     dim_keys = [dim_keys[-1]] + dim_keys[:-1]
    #     dim_keys.append(4.0)

    #     for i, key in enumerate(dim_keys):
    #         dic[f"FDDIMORDER{i+1}"] = key
    #     dic["FDDIMORDER"] = dim_keys

    #     # --- Update QUADFLAGs
    #     quad_keys = [dic.get(f"FDF{i+1}QUADFLAG", 0) for i in range(4)]
    #     quad_keys = [quad_keys[-1]] + quad_keys[:-1]
    #     for i, q in enumerate(quad_keys):
    #         dic[f"FDF{i+1}QUADFLAG"] = q

    #     # --- Update size metadata
    #     dic["FDSIZE"] = data.shape[1]  # second axis after transpose
    #     dic["FDSLICECOUNT"] = data.shape[0]
    #     dic["FDSPECNUM"] = data.shape[0]

    #     if not nohdr:
    #         dic["FDTRANSPOSED"] = (dic.get("FDTRANSPOSED", 0) + 1) % 2

    #     dic = ng.pipe_proc.clean_minmax(dic)
    #     print("\n\n\n")
    #     print(dic)
    #     return dic, data
