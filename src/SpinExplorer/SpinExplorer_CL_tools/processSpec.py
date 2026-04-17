import nmrglue as ng # type: ignore
import os
from SpinExplorer.SpinExplorer_CL_tools.bruker_params_cl import ParameterExtractorBruker
from SpinExplorer.SpinExplorer_CL_tools.convert_nmrglue_cl import Convert_nmrglue
from SpinExplorer.SpinExplorer_CL_tools.make_parameter_file_cl import parameter_write_cl
from SpinExplorer.SpinExplorer_CL_tools.pulse_sequence_parsing import PulseSequenceParser
from SpinExplorer.SpinExplorer_CL_tools.config_register import registry
from pathlib import Path

#from SpinExplorer.SpinConverter.FindingParameters.bruker_parameters import ParameterExtractorBruker

class FindingParameters:

    def __init__(self) -> None:
        """
        Get parameters from folder (focus on Bruker data here)
        """

        found_fid = self.find_nmr_files()

        if not found_fid:
            print("could not find NMR files")
            return 
        self.read_nmr_data()
        self.find_parameters()
    
    def find_nmr_files(self) -> bool:
        """
        Check the current folder for the relevant NMR files 
        required to do conversions.
        """
        if "acqus" in os.listdir("."):
            self.spectrometer = "Bruker"
            self.parameter_file = "acqus"

        elif "procpar" in os.listdir("."):
            self.spectrometer = "Varian"
            self.parameter_file = "propcar"

        elif "acqu" in os.listdir("."):
            self.spectrometer = "Bruker"
            self.parameter_file = "acqu"
            print("we are defaulting to the acqu file")
        
        else:
            print('no parameter file found in this folder - cannot analyse data')
            return False 
        
        if self.spectrometer == "Bruker":
            self.files = [f for f in os.listdir(".") if f in ("ser","fid")]
        
        elif self.spectrometer == "Varian":
            self.files = [f for f in os.listdir(".") if f in ("fid", "origfig")]

        if len(self.files)==0:
            print("No raw FID found to analyse")
            return False 
        
        if len(self.files)>1:
            print("more than one potential raw data file found")
            return False 

        self.file = self.files[0]

        return True
    
    def read_nmr_data(self) -> None:
        """
        Read in NMR data and get info on data dimensions
        """
        if self.spectrometer == "Bruker":
            if os.path.isdir("pdata") and not os.listdir("pdata"):
                os.rename("pdata", "pdata_original")

            try:
                self.nmr_dic, self.nmr_data = ng.bruker.read(dir="./", bin_file = self.file)
            except:
                print("Error: Unable to read the Bruker data in")
        
        elif self.spectrometer == "Varian":
            try:
                self.nmr_dic, self.nmr_data = ng.varian.read(dir="./", bin_file = self.file)
            except:
                print("Error: Unable to read the Varian data in")
        
        self.data_dimensions = len(self.nmr_data.shape)

    def find_parameters(self)->None:
        """
        Finding the relevant parameters from the spectrometer parameter files
        """
        if self.spectrometer == "Bruker":
            # Initialising the class to extract relevant Bruker parameters
            self.params = ParameterExtractorBruker(self)
            # Extracting spectrum dimensions
            self.params.find_size_bruker()
            # Determine if any axis is a pseudo (non-complex) axis
            self.params.find_acquisition_modes_bruker()
            # Finding the spectrum sweep widths
            self.params.find_sw_bruker()
            # Finding the nucleus spectrometer frequencies
            self.params.find_nucleus_frequencies_bruker()
            # Finding the dimension labels
            self.params.find_labels_bruker()
            # Find the acquisition order
            self.params.find_aqseq()
            # Finding the temperature the experiment was performed at
            self.params.find_temperature_bruker()
            # Finding the carrier frequency of each dimension
            self.params.calculate_carrier_frequency_bruker()
            # Finding the digital filter parameters for bruker spectra
            self.params.find_bruker_digital_filter_parameters()
            # Find bruker scaling parameters
            self.params.find_bruker_scaling_parameters()
            # Finding bruker byte order and byte size
            self.params.determine_byte_order()
            self.params.determine_byte_size()

        else:
            # Spectrometer is Varian
            # Initialising the class to extract relevant varian parameters
            self.params = ParameterExtractorVarian(self)
            # Extracting spectrum dimensions
            self.params.find_size_varian()
            # Determine if any axis is a pseudo (non-complex) axis
            self.params.find_axes_pseudo_varian()
            # Finding the spectrum sweep widths
            self.params.find_sw_varian()
            # Finding the nucleus spectrometer frequencies
            self.params.find_nucleus_frequencies_varian()
            # Finding the dimension labels
            self.params.find_labels_varian()
            # Finding the temperature the experiment was performed at
            self.params.find_temperature_varian()
            # Finding the carrier frequency of each dimension
            self.params.calculate_carrier_frequency_varian()
            # Find varian scaling parameters
            self.params.find_varian_scaling_parameters()
    
    def create_dictionary(self) -> dict:
        """
        Create a dictionary to hold the relevant parameters for the
        conversion (note: for bruker data extra information such as
        digital filter parameters also need to be saved).
        """
        if self.spectrometer == "Bruker":
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



def main():
    input_dat = FindingParameters()
    pp_parser = PulseSequenceParser()
    sequence = pp_parser.parse()

    # TODO: introduce possibility of an overwite so we can choose a specific config file
    config = registry.get_default_config(sequence)
    

    nmr_glue_conv = Convert_nmrglue(input_dat.params, input_dat)
    
    params = parameter_write_cl(nmr_glue_conv, config)
    params.write_out_dict(params.dictionary)

    
    config.process_data(pseudo_flag=nmr_glue_conv.params.pseudo_flag)
    
    import subprocess

    # Basic usage - waits for completion
    subprocess.Popen(
    ['SpinView'],
    stdin=subprocess.DEVNULL,
    stdout=subprocess.DEVNULL, 
    stderr=subprocess.DEVNULL,
    start_new_session=True)
    print('opened Spinview...')

if __name__ == "__main__":
    main()