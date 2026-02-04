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
import nmrglue as ng # type: ignore
from typing import Dict
from numpy.typing import NDArray
import traceback


class Convert_nmrglue:
    def __init__(self, params, nmrdata) -> None:
        """
        This class will perform the conversion of the NMR data
        to nmrPipe format using nmrglue.
        """

        # Note we are currently defaulting to remove 
        # digital filter after fourier transform

        self.params = params
        self.nmrdata = nmrdata
        self.params.spectrometer = self.nmrdata.spectrometer
        self.params.bad_point_threshold = 0.0
        self.ndim = len(self.params.size_indirect)+1
        self.get_scaling()
        self.get_complex_and_real_sizes()
        self.sizes = [int(self.complex_sizes[0]/2), *self.complex_sizes[1:]]
        C = ng.convert.converter()
        # Obtain first guesses of dictionary values
        if self.nmrdata.spectrometer == "Bruker":
            dic, data = ng.fileio.bruker.read("./", shape=tuple(self.sizes[::-1]))
        else:
            dic, data = ng.fileio.varian.read("./", shape=tuple(self.sizes[::-1]))

       
        
        self.get_acquisition_modes()
        self.sweep_widths = [self.params.sw_direct, *self.params.sw_indirect]
        self.get_nuclei_frequencies()
        self.get_carrier_frequencies()
        self.get_nucelus_type()
        self.NUS_tick = False
        
        u = self.create_conversion_dictionary()

        try:
            self.perform_conversion(C, u, dic, data)
            # Give an output to say that the conversion was successful
            self.success_output_message()

        except:
            message = traceback.format_exc()
            self.fail_output_message(message)
    
    def get_scaling(self)->None:
        if self.nmrdata.spectrometer == "Bruker":
            self.params.scaling_factor = (1 / self.params.NS) * (2**self.params.NC) * 1000
        else:
            self.params.scaling_factor = (1 / self.params.NS) * 1000

    

    def get_complex_and_real_sizes(self)->None:
        self.complex_sizes = []
        self.real_sizes = []

        for i in range(self.ndim):
            if i == 0:
                self.complex_sizes.append(self.params.size_direct_complex)
            else:
                try:
                    if self.params.labels_correct_order[i]=="ID":
                        labelval = "off"
                    else:
                        labelval = self.params.labels_correct_order[i]
                    self.complex_sizes.append(self.params.indirect_sizes_dict[labelval])
                except:
                    self.complex_sizes.append(self.params.size_indirect[i-1])

        if self.params.pseudo_flag == 0:
            for i in range(self.ndim):
                if i == 0:
                    self.real_sizes.append(self.params.size_direct/2)
                else:
                    try:
                        if self.params.labels_correct_order[i]=="ID":
                            labelval = "off"
                            size = int(self.params.indirect_sizes_dict[labelval])
                        else:
                            labelval = self.params.labels_correct_order[i]
                            size = int(self.params.indirect_sizes_dict[labelval] / 2)
                    except:
                        size = int(self.params.size_indirect[i - 1] / 2)

                    self.real_sizes.append(size)
        
        else:
            for i in range(self.ndim):
                if i == 0:
                    self.real_sizes.append(int(self.params.size_direct/2))
                else:
                    if (
                        self.params.acqusition_modes_indirect[i - 1] == 0
                        or self.params.acqusition_modes_indirect[i - 1] == 1
                    ):
                        # pseudo (real) axis
                        real_size = str(self.params.size_indirect[i - 1])
                    else:
                        real_size = str(int(self.params.size_indirect[i - 1] / 2))
                    self.real_sizes.append(real_size)


    def get_nuclei_frequencies(self)->None:
        self.nuclei_frequencies = []
        
        for i in range(self.ndim):
            if i == 0:
                self.nuclei_frequencies.append(self.params.nucleus_frequencies[0])
            elif self.acq_modes[i]=="Real":
                self.nuclei_frequencies.append(1)
            else:
                self.nuclei_frequencies.append(self.params.nucleus_frequencies[i])

    def get_nucelus_type(self)->None:
        self.nucleus_type = []
        for i in range(self.ndim):
            if i == 0: 
                self.nucleus_type.append(self.params.labels_correct_order[0])
            elif self.acq_modes[i]=='Real':
                self.nucleus_type.append('ID')
            else:
                self.nucleus_type.append(self.params.labels_correct_order[i])


    def get_carrier_frequencies(self)->None:
        self.carrier_frequencies = []
        for i in range(self.ndim):
            if i == 0:
                if (
                    self.params.labels_correct_order[i] == "1H"
                    or self.params.labels_correct_order[i] == "H1"
                    or self.params.labels_correct_order[i] == "H"
                ):
                    self.carrier_frequencies.append(self.params.water_ppm)
                else:
                    self.carrier_frequencies.append(self.params.references_proton[i])
            elif self.params.labels_correct_order[i] == 'ID':
                self.carrier_frequencies.append(0.00)
            else:
                nuc_type = self.params.labels_correct_order[i]
                main_index = i-1 
                for k, reference_label in enumerate(self.params.references_other_labels):
                    if nuc_type in reference_label:
                        main_index = k 
                self.carrier_frequencies.append(self.params.references_other[main_index])


    def get_acquisition_modes(self)->None:
        """
        pull together acuqisition mode from parameter file
        
        :param self: Description
        """
        acquisition_mode_options_direct = ["DQD", "Complex", "Sequential", "Real"]
        direct_mapping = {3:0, 1:2, 2:1, 4:1}
        acquisition_mode_options_indirect = [
                "Complex",
                "States-TPPI",
                "Echo-AntiEcho",
                "TPPI",
                "States",
                "Real",
            ]
        indirect_mapping = {0: 5, 1: 5, 3: 3, 4: 4, 5: 1, 6: 2}

        direct_acq_mode = acquisition_mode_options_direct[direct_mapping[self.params.acqusition_mode_direct]]
        indirect_acq_mode = [acquisition_mode_options_indirect[indirect_mapping[val]]
                                                              for val in self.params.acqusition_modes_indirect]
        self.acq_modes = [direct_acq_mode,*indirect_acq_mode]

    def success_output_message(self):
        """
        Provides an output message to the user to say that the
        conversion is complete.
        """
        
        print("Data conversion to nmrPipe format using nmrglue is complete.")


    def fail_output_message(self, message):
        """
        Provides an output message to the user to say that the
        conversion did not complete correctly
        """

        print("Data conversion to nmrPipe format using nmrglue did not complete correctly. The following error was reported:\n\n" + message)


    def perform_conversion(self, C, u, dic, data):
        """
        Performing any necessary data reshuffling and then
        performing the data conversion to nmrPipe format before
        saving as test.fid
        """

        if len(self.real_sizes) == 2:
            # If have 2D data but nmrglue has read it in as a 1D, need to split it up
            if len(data.shape) == 1:
                data = np.array(
                    np.split(data, int(self.real_sizes[-1]))
                )

        if len(self.complex_sizes) > 1:
            # Check to see if the NUS flag is ticked
            try:
                self.NUS_tick
                nusbox = True
            except:
                nusbox = False

            if(nusbox == True):
                if self.NUS_tick == True:
                    data = self.reshape_nus_data(data)

        # Rance-Kay/Echo-Antiecho reshuffling
        if self.rance_kay == True:
            dic, data = self.rancekay_shuffling(dic, data, u)

        if self.nmrdata.spectrometer == "Bruker":
            if self.params.include_digital_filter:
                # default to remove digital filter before fourier transform
                dic, data = self.remove_digital_filter_fid(dic, data)

            C.from_bruker(dic, data, u)
        else:
            C.from_varian(dic, data, u)
        pdic, pdata = C.to_pipe()
        
        pdata = self.add_intensity_scaling(pdata)

        pdic["FDPIPEFLAG"] = 1.0  # Setting the pipe flag to true

        # For pseudo2D spectra it is necessary to update the dictionary accordingly
        if u[0]["encoding"] == "real" and len(self.sizes) == 2:
            pdic["FDF1TDSIZE"] = u[0]["size"]
            pdic["FDF1FTSIZE"] = u[0]["size"]
            pdic["FDF1APOD"] = u[0]["size"]
            pdic["FDF1QUADFLAG"] = 1.0
            pdic["FDF1OBS"] = 1.0
            pdic["FDF1SW"] = 1.0
            pdic["FDF1ORIG"] = 1.0
            pdic["FD2DPHASE"] = 0

        pdic['FDCOMMENT'] = 'nmrglue'
        ng.pipe.write("test.fid", pdic, pdata, overwrite=True)

    def create_conversion_dictionary(self) -> Dict:
        """
        Create a conversion dictionary based on the relevant current parameters
        in the SpinConverter GUI
        """

        # Initially, the Rance-Kay (Echo-AntiEcho) flag is set to falso
        self.rance_kay = False

        if len(self.complex_sizes) == 1:
            # Create a generic 1D conversion dictionary
            u = {}
            u["ndim"] = 1
            u[0] = {}
            u[0]["time"] = True
            u[0]["freq"] = False

            # Input values into the dictionary
            u = self.populate_conversion_dictionary(u["ndim"], u, 0)

        elif len(self.complex_sizes) == 2:
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
            u[dict_index]["size"] = int(self.complex_sizes[index] / 2)
            
            if self.acq_modes[0] == "Real":
                u[dict_index]["complex"] = False
            else:
                u[dict_index]["complex"] = True

            u[dict_index]["encoding"] = "direct"
        else:
            u[dict_index]["size"] = int(
                self.complex_sizes[index])
            
            if (
                self.acq_modes[index]== "Real"
            ):
                u[dict_index]["encoding"] = "real"
                u[dict_index]["complex"] = False
            elif (
                self.acq_modes[index]== "Complex"
            ):
                u[dict_index]["encoding"] = "complex"
                u[dict_index]["complex"] = True
            elif (
                self.acq_modes[index]== "States"
            ):
                u[dict_index]["encoding"] = "states"
                u[dict_index]["complex"] = True
            elif (
                self.acq_modes[index]== "TPPI"
            ):
                u[dict_index]["encoding"] = "tppi"
                u[dict_index]["complex"] = True
            elif (
                self.acq_modes[index]== "States-TPPI"
            ):
                u[dict_index]["encoding"] = "states-tppi"
                u[dict_index]["complex"] = True
            elif (
                self.acq_modes[index]== "Echo-Antiecho"
                or self.acq_modes[index]== "Echo-AntiEcho"
                or self.acq_modes[index]== "Rance-Kay"
            ):
                u[dict_index]["encoding"] = "complex"
                u[dict_index]["complex"] = True
                self.rance_kay = True

        u[dict_index]["sw"] = float(
            self.sweep_widths[index]
        )
        
        obs = self.nuclei_frequencies[index]

        if(obs!=0.0):
            u[dict_index]["obs"] = obs
        else:
            # obs cannot be equal to zero in nmrglue as it performs car/obs which is undefined if obs is 0
            u[dict_index]['obs'] = 1.0

        u[dict_index]["car"] = (self.carrier_frequencies[index]*u[dict_index]["obs"])
        u[dict_index]["label"] = (self.nucleus_type[index])

        return u
    


    def add_intensity_scaling(self, pdata: NDArray) -> NDArray:
        """
        If the intensity scaling box is not equal to 1 then the FID data
        needs to be scaled by the scaling box number
        """
        try:
            scaling_number = self.params.scaling_factor
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
        for k, value in enumerate(self.complex_sizes):
            if k == 0:
                # Taking the real size for the direct dimension
                shape.append(self.real_sizes[0])
            else:
                shape.append(value)
        shape.reverse()
        shape = tuple(shape)
        nuslist_tuple = ng.bruker.read_nuslist(
            fname=self.nusfile
        )
        
        data = ng.proc_base.expand_nus(data, shape, nuslist_tuple)
        return data
    
    
    
    def find_bruker_initial_point(self, fid, start):
        """
        Estimate initial complex amplitude and phase of the fid
        """
        c = 1j
        amp = np.abs(fid[start])
        ph = 0.0
        n = 0

        for i in range(start - 2, start // 2, -2):
            val = fid[i]
            if np.abs(val) > 0.0:
                ph += np.angle(val)
                n += 1

        if n > 0:
            ph /= n
            c = amp * np.exp(1j * ph)

        return c

    def remove_digital_filter_fid(self, dic, data: NDArray) -> NDArray:
        """
        Removing the Bruker digital filter before Fourier transform
        (post_proc=False). This amounts to a circular shift of the
        data to account for the group delay.
        """
        decim = self.params.decim
        dspfvs = self.params.dspfvs
        grpdly = self.params.grpdly

        # fid = the first slice
        fid = data if data.ndim == 1 else data[0]

        if grpdly == 0.0:
            return data

        # Estimating the phase correction
        start = int(grpdly + 0.5)
        fid = np.asarray(fid, dtype=np.complex128)
        ph0 = 0.0
        init_pt = self.find_bruker_initial_point(fid, start)
        p1 = np.rad2deg(np.angle(init_pt))
        ph0 += p1

        data = ng.proc_base.zf_double(data, 1)

        # ax.plot(np.linspace(0,1,len(data[0])), data[0])

        data = ng.proc_base.fft(data)
        data = ng.proc_base.ps(data, p0=-ph0)
        data = ng.bruker.remove_digital_filter(dic, data, post_proc=True)
        from scipy.signal import hilbert # type: ignore
        # ax.plot(np.linspace(0,1,len(data[0])), data[0].imag)
        data = ng.proc_base.ht(data, data.shape[-1])
        
        data_real = data.real
        data_imag = ng.proc_base.ps(data, p0=180).imag
        data = data_real + 1j*data_imag
        # ax.plot(np.linspace(0,1,len(data[0])), data[0].imag)
        # data = hilbert(data_real, data.shape[0])
        data = ng.proc_base.ifft(data)

        # data = data[..., start:start + data.shape[-1]]
        midpoint = int(data.shape[-1]/2)

        data = data[...,:midpoint:]

        # ax.plot(np.linspace(0,1,len(data[0])), data[0])

        return dic,data



    """
    The functions shown below was originally obtained from nmrglue, 
    followed by customisation.
    
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
        for i, val in enumerate(self.acq_modes):
            if val == "Echo-AntiEcho" or val == "Rance-Kay":
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
