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


class CheckingParameters:
    def __init__(self, notebook, dimension_tabs):
        """
        This class will check the SpinProcess TextControl input
        values to make sure that they are valid
        """
        self.notebook = notebook
        self.dimension_tabs = dimension_tabs

        self.check = True  # Start off having passed the check
        self.check_parameter_validity()

    def check_parameter_validity(self):
        """
        This function will run through all dimensions and perform
        checks on parameters in each dimension.
        If after any dimension the check failed (is False),
        the function will return False and not continue.
        """
        for i, dimension_tab in enumerate(self.dimension_tabs):
            self.check = self.check_parameters(i)
            if self.check == False:
                return self.check

        return self.check

    def check_parameters(self, dimension: int) -> bool:
        """
        Checking that each of the processing parameters is suitable
        for processing. e.g. if a user has entered a word rather than
        a number into a textcontrol.
        Dimension can be 0 (direct) or 1 (indirect). If dimension is indirect, then use different
        checks as the functions are slightly different.

        This function will return True if all the checks are passed,
        otherwise it will return False
        """

        dimension_tab = self.dimension_tabs[dimension]

        check = True  # Initially set to pass (True)

        # Check parameters that could have been edited by a user
        check = self.check_extraction(dimension_tab, dimension)
        if check == False:
            return check
        check = self.check_phasing(dimension_tab, dimension)
        if check == False:
            return check
        check = self.check_baseline(dimension_tab, dimension)
        if check == False:
            return check
        if dimension == 1:
            check = self.check_nus(dimension_tab, dimension)
        return check

    def check_extraction(self, dimension_tab, dimension: int) -> bool:
        """
        This function will check the extraction parameters for
        a given dimension tab.
        """
        if dimension_tab.extraction.extraction_checkbox.GetValue() == True:
            try:
                float(
                    dimension_tab.extraction.extraction_ppm_start_textcontrol.GetValue()
                )
                float(
                    dimension_tab.extraction.extraction_ppm_end_textcontrol.GetValue()
                )
            except:
                dlg = wx.MessageDialog(
                    self.notebook,
                    "Extraction error (Dimension {}): The extraction values must be numbers".format(
                        dimension + 1
                    ),
                    "Warning",
                    wx.OK | wx.ICON_WARNING,
                )
                self.notebook.Raise()
                self.notebook.SetFocus()
                result = dlg.ShowModal()
                # self.change_to_cwd()
                return False
            if float(
                dimension_tab.extraction.extraction_ppm_start_textcontrol.GetValue()
            ) >= float(
                dimension_tab.extraction.extraction_ppm_end_textcontrol.GetValue()
            ):
                dlg = wx.MessageDialog(
                    self.notebook,
                    "Extraction error (Dimension {}): The extraction start value (ppm) must be less than the extraction end value (ppm)".format(
                        dimension + 1
                    ),
                    "Warning",
                    wx.OK | wx.ICON_WARNING,
                )
                self.notebook.Raise()
                self.notebook.SetFocus()
                result = dlg.ShowModal()
                # self.change_to_cwd()
                return False
        return True

    def check_phasing(self, dimension_tab, dimension: int) -> bool:
        """
        Check to see that the phasing values are valid
        """
        try:
            if dimension == 0:
                if dimension_tab.phasing.phase_correction_checkbox.GetValue() == True:
                    float(
                        dimension_tab.phasing.phase_correction_p0_textcontrol.GetValue()
                    )
                    float(
                        dimension_tab.phasing.phase_correction_p1_textcontrol.GetValue()
                    )
            else:
                if (
                    dimension_tab.phasing.phase_correction_checkbox_indirect.GetValue()
                    == True
                ):
                    float(
                        dimension_tab.phasing.phase_correction_p0_textcontrol_indirect.GetValue()
                    )
                    float(
                        dimension_tab.phasing.phase_correction_p1_textcontrol_indirect.GetValue()
                    )
        except:
            dlg = wx.MessageDialog(
                self.notebook,
                "Phasing error (Dimension {}): The phase correction values must be numbers".format(
                    dimension + 1
                ),
                "Warning",
                wx.OK | wx.ICON_WARNING,
            )
            self.notebook.Raise()
            self.notebook.SetFocus()
            result = dlg.ShowModal()
            # self.change_to_cwd()
            return False

        return True

    def check_baseline(self, dimension_tab, dimension: int) -> bool:
        """
        Checking the validity of user inputted parameters for baseline
        correction sections.
        """
        if (
            dimension_tab.baseline_correction.baseline_correction_checkbox.GetValue()
            == True
        ):
            try:
                int(
                    dimension_tab.baseline_correction.baseline_correction_nodes_textcontrol.GetValue()
                )
            except:
                dlg = wx.MessageDialog(
                    self.notebook,
                    "Baseline error (Dimension {}): The node width must be an integer number of points".format(
                        dimension + 1
                    ),
                    "Warning",
                    wx.OK | wx.ICON_WARNING,
                )
                self.notebook.Raise()
                self.notebook.SetFocus()
                result = dlg.ShowModal()
                # self.change_to_cwd()
                return False
            try:
                node_list = dimension_tab.baseline_correction.baseline_correction_node_list_textcontrol.GetValue().split(
                    ","
                )
                node_list_final = []
                for node in node_list:
                    node_list_final.append(float(node))

                if len(node_list_final) == 0:
                    dlg = wx.MessageDialog(
                        self,
                        "Baseline error (Dimension {}): The node list must contain at least one value".format(
                            dimension + 1
                        ),
                        "Warning",
                        wx.OK | wx.ICON_WARNING,
                    )
                    self.notebook.Raise()
                    self.notebook.SetFocus()
                    result = dlg.ShowModal()
                    # self.change_to_cwd()
                    return False

            except:
                dlg = wx.MessageDialog(
                    self.notebook,
                    "Baseline error (Dimension {}): The node list must be a list of comma separated numbers".format(
                        dimension + 1
                    ),
                    "Warning",
                    wx.OK | wx.ICON_WARNING,
                )
                self.notebook.Raise()
                self.notebook.SetFocus()
                result = dlg.ShowModal()
                # self.change_to_cwd()
                return False
            # If polynomial is selected, check to see that the polynomial order is valid
            if (
                dimension_tab.baseline_correction.baseline_correction_radio_box_selection
                == 1
            ):
                try:
                    int(
                        dimension_tab.baseline_correction.baseline_correction_polynomial_order_textcontrol.GetValue()
                    )
                except:
                    dlg = wx.MessageDialog(
                        self.notebook,
                        "Baseline error (Dimension {}): The polynomial order must be an integer for polynomial baselining".format(
                            dimension + 1
                        ),
                        "Warning",
                        wx.OK | wx.ICON_WARNING,
                    )
                    self.notebook.Raise()
                    self.notebook.SetFocus()
                    result = dlg.ShowModal()
                    # self.change_to_cwd()
                    return False

        return True

    def check_nus(self, dimension_tab, dimension: int) -> bool:
        # If SMILE processing is selected, check to see that the SMILE file exists
        if (
            dimension_tab.linear_prediction.linear_prediction_radio_box_indirect.GetSelection()
            == 2
        ):
            # List the files in the current directory
            files = os.listdir()
            # Check to see if the SMILE file exists
            if (
                dimension_tab.linear_prediction.smile_nus_file_textcontrol_indirect.GetValue()
                not in files
            ):
                message = (
                    "SMILE NUS reconstruction error (dimension {}): The NUS file ".format(
                        dimension + 1
                    )
                    + dimension_tab.linear_prediction.smile_nus_file_textcontrol_indirect.GetValue()
                    + " cannot be found in the current directory"
                )
                dlg = wx.MessageDialog(
                    self.notebook, message, "Warning", wx.OK | wx.ICON_WARNING
                )
                self.notebook.Raise()
                self.notebook.SetFocus()
                result = dlg.ShowModal()
                # self.change_to_cwd()
                return False

            # Check that the number of CPUs is an integer
            try:
                int(
                    dimension_tab.linear_prediction.smile_nus_cpu_textcontrol_indirect.GetValue()
                )
            except:
                dlg = wx.MessageDialog(
                    self,
                    "SMILE NUS reconstruction error (dimension {}): The number of CPUs must be an integer".format(
                        dimension + 1
                    ),
                    "Warning",
                    wx.OK | wx.ICON_WARNING,
                )
                self.notebook.Raise()
                self.notebook.SetFocus()
                result = dlg.ShowModal()
                # self.change_to_cwd()
                return False

            # Check that the number of iterations is an integer
            try:
                int(
                    dimension_tab.linear_prediction.smile_nus_iterations_textcontrol_indirect.GetValue()
                )
            except:
                dlg = wx.MessageDialog(
                    self.notebook,
                    "SMILE NUS reconstruction error (dimension {}): The number of iterations must be an integer".format(
                        dimension + 1
                    ),
                    "Warning",
                    wx.OK | wx.ICON_WARNING,
                )
                self.notebook.Raise()
                self.notebook.SetFocus()
                result = dlg.ShowModal()
                # self.change_to_cwd()
                return False

        return True
