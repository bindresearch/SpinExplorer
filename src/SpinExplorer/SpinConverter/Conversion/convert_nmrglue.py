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

import numpy as np
import nmrglue as ng
from typing import Dict
from numpy.typing import NDArray


class Convert_nmrglue:
    def __init__(self, app, params, nmrdata) -> None:
        """
        This class will perform the conversion of the NMR data
        to nmrPipe format using nmrglue.
        """
        self.app = app
        self.params = params
        self.nmrdata = nmrdata

        sizes = []
        for i, box in enumerate(self.app.format.N_complex_boxes):
            size = int(box.GetValue())
            if i == 0:
                size = int(size / 2)
            sizes.append(size)

        sizes.reverse()

        C = ng.convert.converter()
        # Obtain first guesses of dictionary values
        if self.nmrdata.spectrometer == "Bruker":
            dic, data = ng.fileio.bruker.read("./", shape=tuple(sizes))
        else:
            dic, data = ng.fileio.varian.read("./", shape=tuple(sizes))

        u = self.create_conversion_dictionary()

        self.perform_conversion(C, u, dic, data)

        # Give an output to say that the conversion was successful
        self.success_output_message()

    def success_output_message(self):
        """
        Provides an output message to the user to say that the
        conversion is complete.
        """
        dlg = wx.MessageDialog(
            self.app,
            "Data conversion to nmrPipe format using nmrglue is complete.",
            "Complete",
            wx.OK,
        )
        self.app.Raise()
        self.app.SetFocus()
        dlg.ShowModal()
        dlg.Destroy()

    def perform_conversion(self, C, u, dic, data):
        """
        Performing any necessary data reshuffling and then
        performing the data conversion to nmrPipe format before
        saving as test.fid
        """

        if len(self.app.format.N_real_boxes) == 2:
            # If have 2D data but nmrglue has read it in as a 1D, need to split it up
            if len(data.shape) == 1:
                data = np.array(
                    np.split(data, int(self.app.format.N_real_boxes[-1].GetValue()))
                )

        if len(self.app.format.N_complex_boxes) > 1:
            # Check to see if the NUS flag is ticked
            if self.app.shared_format.NUS_tickbox.GetValue() == True:
                data = self.reshape_nus_data(data)

        # Rance-Kay/Echo-Antiecho reshuffling
        if self.rance_kay == True:
            dic, data = self.rancekay_shuffling(dic, data, u)

        if self.nmrdata.spectrometer == "Bruker":
            if self.app.format.digital_filter_checkbox.GetValue() == True:
                if self.app.format.digital_filter_radio_box.GetSelection() == 0:
                    # Removing Bruker digital filter pre-processing
                    # i.e. before Fourier transform
                    data = self.remove_digital_filter_fid(data)

            C.from_bruker(dic, data, u)
        else:
            C.from_varian(dic, data, u)
        pdic, pdata = C.to_pipe()

        pdata = self.add_intensity_scaling(pdata)

        pdic["FDPIPEFLAG"] = 1.0  # Setting the pipe flag to true

        # For pseudo2D spectra it is necessary to update the dictionary accordingly
        if u[0]["encoding"] == "real" and len(self.app.format.N_complex_boxes) == 2:
            pdic["FDF1TDSIZE"] = u[0]["size"]
            pdic["FDF1FTSIZE"] = u[0]["size"]
            pdic["FDF1APOD"] = u[0]["size"]
            pdic["FDF1QUADFLAG"] = 1.0
            pdic["FDF1OBS"] = 1.0
            pdic["FDF1SW"] = 1.0
            pdic["FDF1ORIG"] = 1.0
            pdic["FD2DPHASE"] = 0

        ng.pipe.write("test.fid", pdic, pdata, overwrite=True)

    def create_conversion_dictionary(self) -> Dict:
        """
        Create a conversion dictionary based on the relevant current parameters
        in the SpinConverter GUI
        """

        # Initially, the Rance-Kay (Echo-AntiEcho) flag is set to falso
        self.rance_kay = False

        if len(self.app.format.N_complex_boxes) == 1:
            # Create a generic 1D conversion dictionary
            u = {}
            u["ndim"] = 1
            u[0] = {}
            u[0]["time"] = True
            u[0]["freq"] = False

            # Input values into the dictionary
            u = self.populate_conversion_dictionary(u["ndim"], u, 0)

        elif len(self.app.format.N_complex_boxes) == 2:
            # Create a generic 2D conversion dictionary
            u = {}
            u["ndim"] = 2
            for i in range(2):
                u[i] = {}
                u[i]["time"] = True
                u[i]["freq"] = False

            # Input conversion parameters into the dictionary
            u = self.populate_conversion_dictionary(u["ndim"], u, 0)
            u = self.populate_conversion_dictionary(u["ndim"], u, 1)

        else:
            # Create a generic 3D conversion dictionary
            u = {}
            u["ndim"] = 3
            for i in range(3):
                u[i] = {}
                u[i]["time"] = True
                u[i]["freq"] = False

            # Input conversion parameters into the dictionary
            u = self.populate_conversion_dictionary(u["ndim"], u, 0)
            u = self.populate_conversion_dictionary(u["ndim"], u, 1)
            u = self.populate_conversion_dictionary(u["ndim"], u, 2)

        return u

    def populate_conversion_dictionary(
        self, dimensions: int, u: Dict, index: int
    ) -> dict:
        """
        This function populates the conversion dictionary (u) for
        the current dimension index
        """

        dict_index = dimensions - 1 - index

        if index == 0:
            u[dict_index]["size"] = int(
                int(self.app.format.N_complex_boxes[index].GetValue().strip()) / 2
            )
            if (
                self.app.format.acqusition_combo_boxes[index].GetValue().strip()
                == "Real"
            ):
                u[dict_index]["complex"] = False
            else:
                u[dict_index]["complex"] = True

            u[dict_index]["encoding"] = "direct"
        else:
            u[dict_index]["size"] = int(
                self.app.format.N_complex_boxes[index].GetValue().strip()
            )
            if (
                self.app.format.acqusition_combo_boxes[index].GetValue().strip()
                == "Real"
            ):
                u[dict_index]["encoding"] = "real"
                u[dict_index]["complex"] = False
            elif (
                self.app.format.acqusition_combo_boxes[index].GetValue().strip()
                == "Complex"
            ):
                u[dict_index]["encoding"] = "complex"
                u[dict_index]["complex"] = True
            elif (
                self.app.format.acqusition_combo_boxes[index].GetValue().strip()
                == "States"
            ):
                u[dict_index]["encoding"] = "states"
                u[dict_index]["complex"] = True
            elif (
                self.app.format.acqusition_combo_boxes[index].GetValue().strip()
                == "TPPI"
            ):
                u[dict_index]["encoding"] = "tppi"
                u[dict_index]["complex"] = True
            elif (
                self.app.format.acqusition_combo_boxes[index].GetValue().strip()
                == "States-TPPI"
            ):
                u[dict_index]["encoding"] = "states-tppi"
                u[dict_index]["complex"] = True
            elif (
                self.app.format.acqusition_combo_boxes[index].GetValue().strip()
                == "Echo-Antiecho"
                or self.app.format.acqusition_combo_boxes[index].GetValue().strip()
                == "Echo-AntiEcho"
                or self.app.format.acqusition_combo_boxes[index].GetValue().strip()
                == "Rance-Kay"
            ):
                u[dict_index]["encoding"] = "complex"
                u[dict_index]["complex"] = True
                self.rance_kay = True

        u[dict_index]["sw"] = float(
            self.app.format.sweep_width_boxes[index].GetValue().strip()
        )
        u[dict_index]["obs"] = float(
            self.app.format.nuclei_frequency_boxes[index].GetValue().strip()
        )
        u[dict_index]["car"] = (
            float(self.app.format.carrier_frequency_boxes[index].GetValue().strip())
            * u[dict_index]["obs"]
        )
        u[dict_index]["label"] = (
            self.app.format.nucleus_type_boxes[index].GetValue().strip()
        )

        return u

    def add_intensity_scaling(self, pdata: NDArray) -> NDArray:
        """
        If the intensity scaling box is not equal to 1 then the FID data
        needs to be scaled by the scaling box number
        """
        try:
            scaling_number = float(self.app.shared_format.scaling_number.GetValue())
            pdata = scaling_number * pdata
            return pdata
        except:
            # Multiplication by scaling number did not work
            return pdata

    def reshape_nus_data(self, data: NDArray) -> NDArray:
        """
        Reshaping the NUS FID to the correct order and inserting
        zeros into the missing gaps.
        """
        # Need to reshape the data
        shape = []
        for k, value in enumerate(self.app.format.N_complex_boxes):
            if k == 0:
                # Taking the real size for the direct dimension
                shape.append(int(self.app.format.N_real_boxes[0].GetValue()))
            else:
                shape.append(int(value.GetValue()))
        shape.reverse()
        shape = tuple(shape)
        nuslist_tuple = ng.bruker.read_nuslist(
            fname=self.app.shared_format.nusfile_input.GetValue()
        )
        data = ng.proc_base.expand_nus(data, shape, nuslist_tuple)
        return data

    def remove_digital_filter_fid(self, data: NDArray) -> NDArray:
        """
        Removing the Bruker digital filter before Fourier transform
        (post_proc=False). This amounts to a circular shift of the
        data to account for the group delay.
        """
        decim = float(self.app.format.decim_textbox.GetValue())
        dspfvs = int(self.app.format.dspfvs_textbox.GetValue())
        grpdly = float(self.app.format.grpdly_textbox.GetValue())
        data = ng.bruker.rm_dig_filter(data, decim, dspfvs, grpdly, post_proc=False)
        return data

    def rancekay_shuffling(self, dic, data, udic, rotate_phase=True, **kwargs):
        """
        Reshuffle the data according to the Rance-Kay quadrature scheme.

        Adapted from NMRglue to make general for >2D spectra or spectra
        with multiple Rance-Kay dimensions.

        Parameters
        ----------
        dic : dict
            Dictionary of NMRPipe parameters
        data : ndarray
            Array of NMR data.
        rotate_phase : bool, optional
            Remove the requirement for a 90 degree zero-order phase correction

        Returns
        -------
        ndic : dict
            Dictionary of updated NMRPipe parameters.
        ndata : ndarray
            Array of NMR data which has been reshuffled according to the
            Rance-Kay scheme.

        """

        # Finding which dimensions in udic are Rance-Kay
        rance_kay_dimensions = []
        for i, box in enumerate(self.app.format.acqusition_combo_boxes):
            box = box.GetValue().strip()
            if box == "Echo-AntiEcho" or box == "Rance-Kay":
                rance_kay_dimensions.append((len(data.shape) - 1) - i)

        # Creating an empty array to store the reshuffled data
        shuffled_data = np.empty(data.shape, data.dtype)
        # If final dimension is Rance-Kay/Echo-AntiEcho
        if rance_kay_dimensions == [0]:
            for i in range(0, data.shape[0], 2):
                shuffled_data[i] = (
                    1.0 * (data[i].real - data[i + 1].real)
                    + 1.0 * (data[i].imag - data[i + 1].imag) * 1j
                )
                if rotate_phase is True:
                    shuffled_data[i + 1] = (
                        -1.0 * (data[i].imag + data[i + 1].imag)
                        + 1.0 * (data[i].real + data[i + 1].real) * 1j
                    )
                else:
                    shuffled_data[i + 1] = (
                        1.0 * (data[i].real + data[i + 1].real)
                        + 1.0 * (data[i].imag + data[i + 1].imag) * 1j
                    )

        # If second to last dimension is Rance-Kay/Echo-AntiEcho
        elif rance_kay_dimensions == [1]:
            if len(data.shape) == 3:
                for i in range(0, data.shape[1], 2):
                    shuffled_data[:, i, :] = (
                        1.0 * (data[:, i, :].real - data[:, i + 1, :].real)
                        + 1.0 * (data[:, i, :].imag - data[:, i + 1, :].imag) * 1j
                    )
                    if rotate_phase is True:
                        shuffled_data[:, i + 1, :] = (
                            -1.0 * (data[:, i, :].imag + data[:, i + 1, :].imag)
                            + 1.0 * (data[:, i, :].real + data[:, i + 1, :].real) * 1j
                        )
                    else:
                        shuffled_data[:, i + 1, :] = (
                            1.0 * (data[:, i, :].real + data[:, i + 1, :].real)
                            + 1.0 * (data[:, i, :].imag + data[:, i + 1, :].imag) * 1j
                        )

        return dic, shuffled_data
