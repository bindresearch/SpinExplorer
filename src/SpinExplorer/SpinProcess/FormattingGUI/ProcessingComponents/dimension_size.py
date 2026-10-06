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
import numpy as np


def update_dimension_size(parent):
    """
    Update the dimension size text of a processing tab. This is called by the
    processing options which change the number of points in the dimension.
    """
    dimension_size = getattr(parent, "dimension_size", None)
    if dimension_size != None:
        dimension_size.update_dimension_size()


class DimensionSize:

    def __init__(self, app, nmr_data, parent, dimension):
        """
        This class shows the number of points in the current dimension
        along with the number of points the dimension will have once
        linear prediction, NUS data extension and zero filling have
        been applied.
        """

        self.app = app
        self.nmr_data = nmr_data
        self.parent = parent
        self.dimension_index = dimension

        # Which dimension of the nmr data this tab refers to. This is set by
        # the notebook once all of the tabs have been created as it is not
        # always the same as the tab number (see set_data_dimension)
        self.data_dimension_index = dimension

        self.create_dimension_size_sizer(parent)

    def create_dimension_size_sizer(self, parent):
        """
        Create a box showing the size of the current dimension
        """
        self.dimension_size_box = wx.StaticBox(parent, -1, "Dimension Size")
        self.dimension_size_sizer = wx.StaticBoxSizer(
            self.dimension_size_box, wx.HORIZONTAL
        )
        self.dimension_size_label = wx.StaticText(
            self.dimension_size_box, -1, self.find_dimension_size_text()
        )
        self.dimension_size_sizer.AddSpacer(10)
        self.dimension_size_sizer.Add(
            self.dimension_size_label, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.dimension_size_sizer.AddSpacer(10)

        parent.sizer_1.Add(self.dimension_size_sizer)
        parent.sizer_1.AddSpacer(10)

    def update_dimension_size(self):
        """
        Update the text showing the size of the dimension. This is called
        whenever a processing option which changes the number of points
        is altered.
        """
        try:
            self.dimension_size_label.SetLabel(self.find_dimension_size_text())
            self.dimension_size_sizer.Layout()
            self.parent.Layout()
        except (RuntimeError, AttributeError):
            # The label is no longer shown (the interface is being rebuilt)
            pass

    def set_data_dimension(self, dimension: int):
        """
        Set which dimension of the nmr data this tab refers to. This is the
        same as the tab number apart from pseudo 3D datasets, where the pseudo
        axis is not given a processing tab, and the final indirect tab of a 3D,
        which is created using the same class as the first indirect tab.
        """
        self.data_dimension_index = dimension
        self.update_dimension_size()

    def find_current_size(self) -> int:
        """
        The number of complex points currently in this dimension. The direct
        dimension is stored as complex points whereas the indirect dimensions
        are stored as interleaved real/imaginary points.
        """
        dimension = self.data_dimension_index
        points = int(self.nmr_data.number_of_points[dimension])

        if dimension == 0:
            return points

        return int(points / 2)

    def find_truncation_size(self, size: int):
        """
        The number of complex points after truncation, where fewer points are
        processed than were recorded.
        """
        truncation = getattr(self.parent, "truncation", None)
        if truncation == None:
            return size, False

        if truncation.truncation_checkbox_value == False:
            return size, False

        try:
            points = int(truncation.find_truncation_points())
        except (ValueError, TypeError):
            return size, False

        if points < 1 or points >= size:
            return size, False

        return points, True

    def find_linear_prediction_size(self, size: int):
        """
        The number of complex points after linear prediction. Predicting
        points after the FID doubles the number of points (or adds a single
        point when only backward coefficients are used), predicting points
        before the FID replaces the first points and so does not change
        the size.
        """
        linear_prediction = getattr(self.parent, "linear_prediction", None)
        if linear_prediction == None:
            return size, False

        if self.dimension_index == 0:
            if linear_prediction.linear_prediction_checkbox_value == False:
                return size, False
            append = linear_prediction.linear_prediction_options_selection
            coefficients = linear_prediction.linear_prediction_coefficients_selection
        else:
            if linear_prediction.linear_prediction_radio_box_indirect_selection != 1:
                return size, False
            append = linear_prediction.linear_prediction_indirect_options_selection
            coefficients = (
                linear_prediction.linear_prediction_indirect_coefficients_selection
            )

        if append != 0:
            # Points are predicted before the FID, the size is unchanged
            return size, True

        if coefficients == 1:
            # Backward coefficients only predict a single point
            return size + 1, True

        return size * 2, True

    def find_extension_size(self, size: int):
        """
        The number of complex points after the NUS data extension. Only the
        indirect dimensions can be extended.
        """
        linear_prediction = getattr(self.parent, "linear_prediction", None)
        if linear_prediction == None or self.dimension_index == 0:
            return size, False

        selection = linear_prediction.linear_prediction_radio_box_indirect_selection
        if selection == 2:
            extension = linear_prediction.smile_data_extension_number_indirect
        elif selection == 3:
            extension = linear_prediction.ist_data_extension_number_indirect
        else:
            return size, False

        try:
            extension = int(extension)
        except (ValueError, TypeError):
            return size, False

        if extension <= 0:
            return size, False

        return size + extension, True

    def find_zero_filling_size(self, size: int):
        """
        The number of complex points after zero filling. This follows the
        same calculation as the zero filling applied during processing.
        """
        zero_filling = getattr(self.parent, "zero_filling", None)
        if zero_filling == None:
            return size, False

        if zero_filling.zero_filling_checkbox_value == False:
            return size, False

        selection = zero_filling.zero_filling_combobox_selection
        try:
            if selection == 0:
                final_size = size * 2 ** int(
                    zero_filling.zero_filling_value_doubling_times
                )
            elif selection == 1:
                final_size = size + int(zero_filling.zero_filling_value_zeros_to_add)
            else:
                final_size = int(zero_filling.zero_filling_value_final_data_size)
        except (ValueError, TypeError):
            return size, False

        if final_size < size:
            # Zero filling never removes points
            final_size = size

        if zero_filling.zero_filling_round_checkbox_value == True:
            final_size = int(2 ** (np.ceil(np.log(final_size) / np.log(2))))

        if final_size == size:
            return size, False

        return final_size, True

    def find_dimension_size_text(self) -> str:
        """
        Create the text showing the current size of the dimension and, if
        any of the processing options change the number of points, the size
        of the dimension after they have been applied.
        """
        current_size = self.find_current_size()

        try:
            size, truncation = self.find_truncation_size(current_size)
            size, linear_prediction = self.find_linear_prediction_size(size)
            size, extension = self.find_extension_size(size)
            size, zero_filling = self.find_zero_filling_size(size)
        except (RuntimeError, AttributeError):
            # A processing section is in the middle of being rebuilt, the text
            # is updated again once the interface has been recreated
            size, truncation, linear_prediction, extension, zero_filling = (
                current_size,
                False,
                False,
                False,
                False,
            )

        text = "Current size: {} complex points ({} real and imaginary points)".format(
            current_size, 2 * current_size
        )

        applied = []
        if truncation == True:
            applied.append("truncation")
        if linear_prediction == True:
            applied.append("linear prediction")
        if extension == True:
            applied.append("data extension")
        if zero_filling == True:
            applied.append("zero filling")

        if len(applied) == 0:
            return text

        if len(applied) == 1:
            applied_text = applied[0]
        else:
            applied_text = ", ".join(applied[:-1]) + " and " + applied[-1]

        text += "\nSize after {}: {} complex points ({} real and imaginary points)".format(
            applied_text, size, 2 * size
        )

        return text
