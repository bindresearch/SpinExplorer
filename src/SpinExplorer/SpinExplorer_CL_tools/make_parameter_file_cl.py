import json
from SpinExplorer.SpinExplorer_CL_tools.convert_nmrglue_cl import Convert_nmrglue
from SpinExplorer.SpinExplorer_CL_tools.experiment_config import ExperimentConfigStore

class parameter_write_cl:
    def __init__(self, conv_object: Convert_nmrglue, processing_config: ExperimentConfigStore):
        self.conv = conv_object
        self.dictionary = self.create_dictionary()
        self.dictionary = self.add_general(self.dictionary)
        self.dictionary = self.add_complex_and_real_sizes(self.dictionary)
        self.dictionary = self.add_acqusition_mode(self.dictionary)
        self.dictionary = self.add_sweep_widths(self.dictionary)
        self.dictionary = self.add_nuclei_frequency(self.dictionary)
        self.dictionary = self.add_labels(self.dictionary)
        self.dictionary = self.add_carrier_frequency(self.dictionary)
        self.dictionary = self.add_scaling_information(self.dictionary)
        self.dictionary = self.add_digital_filter_information(self.dictionary)
        self.dictionary = self.add_other_options_information(self.dictionary)
        self.dictionary = self.add_nus_information(self.dictionary)
        self.dictionary.update(processing_config.make_processing_dictionary())

        self.write_out_dict(self.dictionary)


    def create_dictionary(self) -> dict:
            """
            Create a dictionary to hold the relevant parameters for the
            conversion (note: for bruker data extra information such as
            digital filter parameters also need to be saved).
            """
            if self.conv.params.spectrometer == "Bruker":
                dictionary = {
                    "conversion": {
                        "general": {"spectrometer": "Bruker"},
                        "spectral parameters": {},
                        "intensity scaling": {},
                        "digital filter parameters": {},
                        "other parameters": {},
                    }
                }
            else:
                dictionary = {
                    "conversion": {
                        "general": {"spectrometer": "Varian"},
                        "spectral parameters": {},
                        "intensity scaling": {},
                    }
                }
            return dictionary

    def add_general(self, dictionary: dict) -> dict:
        """
        Adding general information such as spectrometer type
        and temperature and number of scans
        """
        dictionary["conversion"]["general"]["temperature"] = str(self.conv.params.temperature)
        dictionary["conversion"]["general"]["number of scans (NS)"] = str(self.conv.params.NS)
        return dictionary

    def add_complex_and_real_sizes(self, dictionary: dict) -> dict:
        """
        Adding the complex dimension sizes to the dictionary
        """
        dictionary["conversion"]["spectral parameters"]["sizes"] = {}
        dictionary["conversion"]["spectral parameters"]["sizes"]["complex"] = [str(int(size)) for size in self.conv.complex_sizes]
        dictionary["conversion"]["spectral parameters"]["sizes"]["real"] = [str(int(size)) for size in self.conv.real_sizes]
        return dictionary

    def add_acqusition_mode(self, dictionary: dict) -> dict:
        """
        Adding the acqusition mode for each dimension
        (both name and selection index)
        """
        indirect_modes = []
        indirect_mode_indexes = []
        dictionary["conversion"]["spectral parameters"]["acqusition modes"] = {}
        dictionary["conversion"]["spectral parameters"]["acqusition modes"]["direct"] = {}
        dictionary["conversion"]["spectral parameters"]["acqusition modes"]["indirect"] = {}
        for i, acq_mode in enumerate(self.conv.acq_modes):
            if i == 0:
                dictionary["conversion"]["spectral parameters"]["acqusition modes"]["direct"][
                    "mode"
                ] = acq_mode
                dictionary["conversion"]["spectral parameters"]["acqusition modes"]["direct"][
                    "index"
                ] = str(i)

            else:
                indirect_modes.append(acq_mode)
                indirect_mode_indexes.append(str(i))

        if indirect_modes != []:
            dictionary["conversion"]["spectral parameters"]["acqusition modes"]["indirect"][
                "mode"
            ] = indirect_modes
            dictionary["conversion"]["spectral parameters"]["acqusition modes"]["indirect"][
                "index"
            ] = indirect_mode_indexes

        return dictionary
    
    def add_sweep_widths(self, dictionary: dict) -> dict:
        """
        Adding the sweep width from each dimension to the dictionary
        """
        sweep_widths = []
        sweep_width_indexes = []
        for i,val in enumerate(self.conv.sweep_widths):
            sweep_widths.append(str(val))
            sweep_width_indexes.append(str(i))

        dictionary["conversion"]["spectral parameters"]["sweep widths"] = {}
        dictionary["conversion"]["spectral parameters"]["sweep widths"]["values"] = sweep_widths
        dictionary["conversion"]["spectral parameters"]["sweep widths"]["indexes"] = sweep_width_indexes

        return dictionary
    
    def add_nuclei_frequency(self, dictionary: dict) -> dict:
        """
        Adding the nucleus frequency from each dimension to the dictionary
        """
        dictionary["conversion"]["spectral parameters"]["nuclei frequencies"] = [str(freq) for freq in self.conv.nuclei_frequencies]

        return dictionary
    
    def add_labels(self, dictionary: dict) -> dict:
        """
        Adding the labels from each dimension to the dictionary
        """
        dictionary["conversion"]["spectral parameters"]["labels"] = [nuc for nuc in self.conv.nucleus_type]
        return dictionary
    
    def add_carrier_frequency(self, dictionary: dict) -> dict:
        """
        Adding the carrier frequencies from each dimension to the dictionary
        """

        dictionary["conversion"]["spectral parameters"]["carrier frequencies"] = {}
        dictionary["conversion"]["spectral parameters"]["carrier frequencies"][
            "frequency"
        ] = self.conv.carrier_frequencies
        dictionary["conversion"]["spectral parameters"]["carrier frequencies"][
            "combobox"
        ] = self.conv.params.labels_correct_order
        dictionary["conversion"]["spectral parameters"]["carrier frequencies"][
            "combobox index"
        ] = [i for i in range(self.conv.ndim)]

        return dictionary
    
    def add_scaling_information(self, dictionary: dict) -> dict:
        """
        Adding the spectrum scaling information to the dictionary
        """

        dictionary["conversion"]["intensity scaling"][
            "Scale by number of scans (NS)"
        ] = True
        if(dictionary["conversion"]["general"]["spectrometer"] == "Bruker"):
            dictionary["conversion"]["intensity scaling"][
                "Scale by Bruker normalisation constant (NC)"
            ] = True
        dictionary["conversion"]["intensity scaling"]["Scale by 1000"] = True
        dictionary["conversion"]["intensity scaling"]["Scaling number"] = str(self.conv.params.scaling_factor)

        return dictionary
    
    def add_digital_filter_information(self, dictionary: dict) -> dict:
            
        """
        Adding the digital filter information to the dictionary
        from the Bruker data
        """
        decim = self.conv.params.decim
        dspfvs = self.conv.params.dspfvs
        grpdly = self.conv.params.grpdly

        # Remove digital filter checkbox
        remove_digital_filter = True

        if(self.conv.params.remove_filter_before_processing==True):
            remove='Before'
        else:
            remove='After'


        dictionary["conversion"]["digital filter parameters"]["Remove Digital Filter"] = remove_digital_filter
        dictionary["conversion"]["digital filter parameters"]["Remove Before/After Fourier Transform"] = remove
        dictionary["conversion"]["digital filter parameters"]["Decimation Rate (decim)"] = str(decim)
        dictionary["conversion"]["digital filter parameters"]["DSP Firmware Version (dspfvs)"] = str(dspfvs)
        dictionary["conversion"]["digital filter parameters"]["Group Delay (grpdly)"] = str(grpdly)

        return dictionary
    
    def add_other_options_information(self, dictionary: dict) -> dict:
        """
        Adding extra information to the dictionary for Bruker data
        """
        dictionary["conversion"]["other parameters"]["remove acqusition padding"] = True
        dictionary["conversion"]["other parameters"]["bad point threshold"] = self.conv.params.bad_point_threshold
        return dictionary

    def add_nus_information(self, dictionary: dict) -> dict:
        """
        Adding the current NUS information to the dictionary
        """
        dictionary["conversion"]["NUS information"] = "N/A"

        return dictionary
    

    def write_out_dict(self, dictionary)->None:
        import json

        with open('parameters.json','w') as outy:
            json.dump(dictionary, outy, indent = 4)