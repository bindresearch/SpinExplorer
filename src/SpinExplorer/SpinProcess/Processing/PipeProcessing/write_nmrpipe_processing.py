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
import os
import subprocess


class WriteNMRPipeProcessing:
    def __init__(self, notebook, dimension_tabs, nmr_data):
        """
        This class contains functions to return strings for
        each nmrpipe processing line
        """
        self.notebook = notebook
        self.dimension_tabs = dimension_tabs
        self.nmr_data = nmr_data

    def solvent_suppression(self, dimension, nmrproc_com):
        """
        This function will write the processing line for
        nmrPipe digital solvent suppression.
        """

        dimension_tab = self.dimension_tabs[dimension]

        solvent_suppression_line = "| nmrPipe -fn SOL"
        if dimension_tab.solvent_suppression.solvent_suppression_filter_selection == 0:
            solvent_suppression_line += " -mode 1"
            if (
                dimension_tab.solvent_suppression.solvent_suppression_lowpass_shape_selection
                == 0
            ):
                solvent_suppression_line += " -fs 1"
            elif (
                dimension_tab.solvent_suppression.solvent_suppression_lowpass_shape_selection
                == 1
            ):
                solvent_suppression_line += " -fs 2"
            elif (
                dimension_tab.solvent_suppression.solvent_suppression_lowpass_shape_selection
                == 2
            ):
                solvent_suppression_line += " -fs 3"
            solvent_suppression_line += " -fl " + str(
                int(dimension_tab.solvent_suppression.solvent_suppression_filter_length)
            )
        elif (
            dimension_tab.solvent_suppression.solvent_suppression_filter_selection == 1
        ):
            solvent_suppression_line += " -mode 2"
        elif (
            dimension_tab.solvent_suppression.solvent_suppression_filter_selection == 2
        ):
            solvent_suppression_line += " -mode 3"

        nmrproc_com.write(solvent_suppression_line + " \\\n")
        return nmrproc_com

    def linear_prediction(self, dimension, nmrproc_com):
        """
        This function will write the processing line for
        nmrPipe linear prediction.
        """

        dimension_tab = self.dimension_tabs[dimension]
        linear_prediction_line = "| nmrPipe -fn LP"
        if dimension_tab.linear_prediction.linear_prediction_options_selection == 0:
            linear_prediction_line += " -after"
        elif dimension_tab.linear_prediction.linear_prediction_options_selection == 1:
            linear_prediction_line += " -before"
        if (
            dimension_tab.linear_prediction.linear_prediction_coefficients_selection
            == 0
        ):
            linear_prediction_line += " -f"
        elif (
            dimension_tab.linear_prediction.linear_prediction_coefficients_selection
            == 1
        ):
            linear_prediction_line += " -b"
        elif (
            dimension_tab.linear_prediction.linear_prediction_coefficients_selection
            == 2
        ):
            linear_prediction_line += " -fb"
        nmrproc_com.write(linear_prediction_line + " \\\n")

        return nmrproc_com

    def apodization(self, dimension, nmrproc_com):
        """
        This function will write the processing line for
        nmrPipe apodization.
        """
        dimension_tab = self.dimension_tabs[dimension]
        apodization_line = "| nmrPipe -fn"
        if dimension_tab.apodization.apodization_combobox_selection == 0:
            apodization_line += " EM"
            # Input a line broadening of 0 and first point scaling
            apodization_line += " -lb 0 -c 0.5"
        elif dimension_tab.apodization.apodization_combobox_selection == 1:
            apodization_line += " EM"
            apodization_line += (
                " -lb "
                + str(dimension_tab.apodization.exponential_line_broadening)
                + " -c "
                + str(dimension_tab.apodization.apodization_first_point_scaling)
            )
        elif dimension_tab.apodization.apodization_combobox_selection == 2:
            apodization_line += " GM"
            apodization_line += (
                " -g1 "
                + str(dimension_tab.apodization.g1)
                + " -g2 "
                + str(dimension_tab.apodization.g2)
                + " -g3 "
                + str(dimension_tab.apodization.g3)
                + " -c "
                + str(dimension_tab.apodization.apodization_first_point_scaling)
            )
        elif dimension_tab.apodization.apodization_combobox_selection == 3:
            apodization_line += " SP"
            apodization_line += (
                " -off "
                + str(dimension_tab.apodization.offset)
                + " -end "
                + str(dimension_tab.apodization.end)
                + " -pow "
                + str(dimension_tab.apodization.power)
                + " -c "
                + str(dimension_tab.apodization.apodization_first_point_scaling)
            )
        elif dimension_tab.apodization.apodization_combobox_selection == 4:
            apodization_line += " GMB"
            apodization_line += (
                " -lb "
                + str(dimension_tab.apodization.a)
                + " -gb "
                + str(dimension_tab.apodization.b)
                + " -c "
                + str(dimension_tab.apodization.apodization_first_point_scaling)
            )
        elif dimension_tab.apodization.apodization_combobox_selection == 5:
            apodization_line += " TM"
            apodization_line += (
                " -t1 "
                + str(dimension_tab.apodization.t1)
                + " -t2 "
                + str(dimension_tab.apodization.t2)
                + " -c "
                + str(dimension_tab.apodization.apodization_first_point_scaling)
            )
        elif dimension_tab.apodization.apodization_combobox_selection == 6:
            apodization_line += " TRI"
            apodization_line += (
                " -loc "
                + str(dimension_tab.apodization.loc)
                + " -c "
                + str(dimension_tab.apodization.apodization_first_point_scaling)
            )
        nmrproc_com.write(apodization_line + " \\\n")

        return nmrproc_com

    def zero_filling(self, dimension, nmrproc_com):
        """
        This function will write the processing line for
        nmrPipe zero filling.
        """

        dimension_tab = self.dimension_tabs[dimension]
        zero_filling_line = "| nmrPipe -fn ZF"
        if dimension_tab.zero_filling.zero_filling_combobox_selection == 0:
            zero_filling_line += " -zf " + str(
                dimension_tab.zero_filling.zero_filling_value_doubling_times
            )
        if dimension_tab.zero_filling.zero_filling_combobox_selection == 1:
            zero_filling_line += " -pad " + str(
                dimension_tab.zero_filling.zero_filling_value_zeros_to_add
            )
        elif dimension_tab.zero_filling.zero_filling_combobox_selection == 2:
            zero_filling_line += " -size " + str(
                dimension_tab.zero_filling.zero_filling_value_final_data_size
            )
        if dimension_tab.zero_filling.zero_filling_round_checkbox.GetValue() == True:
            zero_filling_line += " -auto"
        nmrproc_com.write(zero_filling_line + " \\\n")

        return nmrproc_com

    def fourier_transform(self, dimension, nmrproc_com):
        """
        This function will write the processing line for
        nmrPipe fourier transform.
        """
        dimension_tab = self.dimension_tabs[dimension]
        fourier_transform_line = "| nmrPipe -fn FT"
        if dimension_tab.fourier_transform.ft_method_selection == 0:
            fourier_transform_line += " -auto"
        elif dimension_tab.fourier_transform.ft_method_selection == 1:
            fourier_transform_line += " -real"
        elif dimension_tab.fourier_transform.ft_method_selection == 2:
            fourier_transform_line += " -inv"
        elif dimension_tab.fourier_transform.ft_method_selection == 3:
            fourier_transform_line += " -alt"
        elif dimension_tab.fourier_transform.ft_method_selection == 4:
            fourier_transform_line += " -neg"
        nmrproc_com.write(fourier_transform_line + " \\\n")

        return nmrproc_com

    def phasing(self, dimension, nmrproc_com):
        """
        This function will write the processing line for
        nmrPipe phasing.
        """

        dimension_tab = self.dimension_tabs[dimension]

        phase_correction_line = "| nmrPipe -fn PS"
        if dimension == 0:
            phase_correction_line += (
                " -p0 "
                + str(dimension_tab.phasing.phase_correction_p0_textcontrol.GetValue())
                + " -p1 "
                + str(dimension_tab.phasing.phase_correction_p1_textcontrol.GetValue())
                + " -di "
            )
        else:
            phase_correction_line += (
                " -p0 "
                + str(
                    dimension_tab.phasing.phase_correction_p0_textcontrol_indirect.GetValue()
                )
                + " -p1 "
                + str(
                    dimension_tab.phasing.phase_correction_p1_textcontrol_indirect.GetValue()
                )
                + " -di "
            )
        nmrproc_com.write(phase_correction_line + "\\\n")

        return nmrproc_com

    def magnitude_mode(self, nmrproc_com):
        """
        This function will write the processing line for
        nmrPipe magnitude mode.
        """
        mc_line = "| nmrPipe -fn MC"
        nmrproc_com.write(mc_line + "\\\n")

        return nmrproc_com

    def extraction(self, dimension, nmrproc_com):
        """
        This function will write the processing line for
        nmrPipe chemical shift extraction.
        """

        dimension_tab = self.dimension_tabs[dimension]
        extraction_line = "| nmrPipe -fn EXT"
        extraction_line += (
            " -x1 "
            + str(dimension_tab.extraction.extraction_ppm_start_textcontrol.GetValue())
            + "ppm -xn "
            + str(dimension_tab.extraction.extraction_ppm_end_textcontrol.GetValue())
            + "ppm -sw "
        )
        nmrproc_com.write(extraction_line + " \\\n")
        return nmrproc_com

    def baseline_correction(self, dimension, nmrproc_com):
        """
        This function will write the processing line for
        nmrPipe chemical shift extraction.
        """

        dimension_tab = self.dimension_tabs[dimension]

        if (
            dimension_tab.baseline_correction.baseline_correction_radio_box_selection
            == 0
        ):
            baseline_correction_line = "| nmrPipe -fn BASE"
        else:
            baseline_correction_line = "| nmrPipe -fn POLY"
        # add the node width and node list to the baseline correction line
        baseline_correction_line += (
            " -nw "
            + str(
                dimension_tab.baseline_correction.baseline_correction_nodes_textcontrol.GetValue()
            )
            + " -nl "
        )
        node_list = dimension_tab.baseline_correction.baseline_correction_node_list_textcontrol.GetValue().split(
            ","
        )
        node_list_final = []
        for node in node_list:
            node_list_final.append(float(node))
        for node in node_list_final:
            baseline_correction_line += str(node) + "% "

        if (
            dimension_tab.baseline_correction.baseline_correction_radio_box_selection
            == 1
        ):
            baseline_correction_line += " -ord " + str(
                dimension_tab.baseline_correction.baseline_correction_polynomial_order_textcontrol.GetValue()
            )

        nmrproc_com.write(baseline_correction_line + "\\\n")
        return nmrproc_com

    def transpose_line(self, nmrproc_com):
        """
        This function will write the processing line for
        nmrPipe transposing.
        """
        nmrproc_com.write("| nmrPipe -fn TP \\\n")
        return nmrproc_com

    def zero_transpose_line(self, nmrproc_com):
        """
        This function will write the processing line for
        nmrPipe zero transposing.
        """
        nmrproc_com.write("| nmrPipe -fn ZTP \\\n")
        return nmrproc_com

    def apply_NUS_SMILE(self, nmrproc_com, nus_dimensions):
        """
        Apply the relevent lines for NUS reconstruction using
        SMILE.
        """
        if nus_dimensions == [1, 2]:
            dimension_tab1 = self.dimension_tabs[1]
            dimension_tab2 = self.dimension_tabs[2]
            nmrproc_com.write("| nmrPipe -fn ZTP \\\n")
            nmrproc_com.write("| nmrPipe -fn TP \\\n")

        elif nus_dimensions == [1]:
            dimension_tab1 = self.dimension_tabs[1]
        elif nus_dimensions == [2]:
            dimension_tab2 = self.dimension_tabs[2]

        smile_line = (
            "| nmrPipe -fn SMILE -nDim "
            + str(int(self.nmr_data.dim))
            + " -nThread "
            + str(
                dimension_tab1.linear_prediction.smile_nus_cpu_textcontrol_indirect.GetValue()
            )
            + " -report 1 -sample "
            + str(dimension_tab1.linear_prediction.nuslist_name_indirect)
            + "\\\n"
            + "              -maxIter "
            + str(
                int(
                    dimension_tab1.linear_prediction.smile_nus_iterations_textcontrol_indirect.GetValue()
                )
            )
            + "        \\\n"
        )

        if 1 in nus_dimensions:
            smile_line += "       -xP0 {} -xP1 {} -xT {} ".format(
                dimension_tab1.phasing.phase_correction_p0_textcontrol_indirect.GetValue(),
                dimension_tab1.phasing.phase_correction_p1_textcontrol_indirect.GetValue(),
                dimension_tab1.linear_prediction.smile_nus_extension_textcontrol_indirect.GetValue(),
            )
            nmrproc_com.write(smile_line + " \\\n")
            if dimension_tab1.apodization.apodization_checkbox_value == True:
                if dimension_tab1.apodization.apodization_combobox_selection == 0:
                    smile_line_apod = " -xApod EM -xQ1 0.0 -xQ2 0.0 -xQ3 0.0 \\\n"
                elif dimension_tab1.apodization.apodization_combobox_selection == 1:
                    smile_line_apod = (
                        " -xApod EM -xQ1 {} -xQ2 0.0 -xQ3 0.0 \\\n".format(
                            str(dimension_tab1.apodization.exponential_line_broadening)
                        )
                    )
                elif dimension_tab1.apodization.apodization_combobox_selection == 2:
                    smile_line_apod = " -xApod GM -xQ1 {} -xQ2 {} -xQ3 {} \\\n".format(
                        str(dimension_tab1.apodization.g1),
                        str(dimension_tab1.apodization.g2),
                        str(dimension_tab1.apodization.g3),
                    )
                elif dimension_tab1.apodization.apodization_combobox_selection == 3:
                    smile_line_apod = " -xApod SP -xQ1 {} -xQ2 {} -xQ3 {} \\\n".format(
                        str(dimension_tab1.apodization.offset),
                        str(dimension_tab1.apodization.end),
                        str(dimension_tab1.apodization.power),
                    )
                elif dimension_tab1.apodization.apodization_combobox_selection == 4:
                    smile_line_apod = (
                        " -xApod GMB -xQ1 {} -xQ2 {} -xQ3 0.0 \\\n".format(
                            str(dimension_tab1.apodization.a),
                            str(dimension_tab1.apodization.b),
                        )
                    )
                elif dimension_tab1.apodization.apodization_combobox_selection == 5:
                    smile_line_apod = " -xApod TM -xQ1 {} -xQ2 {} -xQ3 0.0 \\\n".format(
                        str(dimension_tab1.apodization.t1),
                        str(dimension_tab1.apodization.t2),
                    )
                elif dimension_tab1.apodization.apodization_combobox_selection == 5:
                    smile_line_apod = (
                        " -xApod TRI -xQ1 {} -xQ2 0.0 -xQ3 0.0 \\\n".format(
                            str(dimension_tab1.apodization.loc)
                        )
                    )
            nmrproc_com.write(smile_line_apod)

        if nus_dimensions == [2]:
            smile_line2 = (
                "| nmrPipe -fn SMILE -nDim "
                + str(int(self.nmr_data.dim))
                + " -nThread "
                + str(
                    dimension_tab1.linear_prediction.smile_nus_cpu_textcontrol_indirect.GetValue()
                )
                + " -report 1 -sample "
                + str(dimension_tab1.linear_prediction.nuslist_name_indirect)
                + "\\\n"
                + "              -maxIter "
                + str(
                    int(
                        dimension_tab1.linear_prediction.smile_nus_iterations_textcontrol_indirect.GetValue()
                    )
                )
                + "        \\\n"
            )
        else:
            smile_line2 = ""

        if 2 in nus_dimensions:
            smile_line2 += "       -yP0 {} -yP1 {} -yT {} ".format(
                dimension_tab2.phasing.phase_correction_p0_textcontrol_indirect.GetValue(),
                dimension_tab2.phasing.phase_correction_p1_textcontrol_indirect.GetValue(),
                dimension_tab2.linear_prediction.smile_nus_extension_textcontrol_indirect.GetValue(),
            )
            nmrproc_com.write(smile_line2 + " \\\n")
            if dimension_tab2.apodization.apodization_checkbox_value == True:
                if dimension_tab2.apodization.apodization_combobox_selection == 0:
                    smile_line_apod = " -yApod EM -yQ1 0.0 -yQ2 0.0 -yQ3 0.0 \\\n"
                elif dimension_tab2.apodization.apodization_combobox_selection == 1:
                    smile_line_apod = (
                        " -yApod EM -yQ1 {} -yQ2 0.0 -yQ3 0.0 \\\n".format(
                            str(dimension_tab2.apodization.exponential_line_broadening)
                        )
                    )
                elif dimension_tab2.apodization.apodization_combobox_selection == 2:
                    smile_line_apod = " -yApod GM -yQ1 {} -yQ2 {} -yQ3 {} \\\n".format(
                        str(dimension_tab2.apodization.g1),
                        str(dimension_tab2.apodization.g2),
                        str(dimension_tab2.apodization.g3),
                    )
                elif dimension_tab2.apodization.apodization_combobox_selection == 3:
                    smile_line_apod = " -yApod SP -yQ1 {} -yQ2 {} -yQ3 {} \\\n".format(
                        str(dimension_tab2.apodization.offset),
                        str(dimension_tab2.apodization.end),
                        str(dimension_tab2.apodization.power),
                    )
                elif dimension_tab2.apodization.apodization_combobox_selection == 4:
                    smile_line_apod = (
                        " -yApod GMB -yQ1 {} -yQ2 {} -yQ3 0.0 \\\n".format(
                            str(dimension_tab2.apodization.a),
                            str(dimension_tab2.apodization.b),
                        )
                    )
                elif dimension_tab2.apodization.apodization_combobox_selection == 5:
                    smile_line_apod = " -yApod TM -yQ1 {} -yQ2 {} -yQ3 0.0 \\\n".format(
                        str(dimension_tab2.apodization.t1),
                        str(dimension_tab2.apodization.t2),
                    )
                elif dimension_tab2.apodization.apodization_combobox_selection == 5:
                    smile_line_apod = (
                        " -yApod TRI -yQ1 {} -yQ2 0.0 -yQ3 0.0 \\\n".format(
                            str(dimension_tab2.apodization.loc)
                        )
                    )
            nmrproc_com.write(smile_line_apod)

        if nus_dimensions == [1, 2]:
            nmrproc_com.write("| nmrPipe -fn TP \\\n")
            nmrproc_com.write("| nmrPipe -fn ZTP \\\n")

        return nmrproc_com
