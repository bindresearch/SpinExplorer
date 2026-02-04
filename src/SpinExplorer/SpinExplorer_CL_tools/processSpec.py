import nmrglue as ng # type: ignore
import os
from bruker_params_cl import ParameterExtractorBruker
from convert_nmrglue_cl import Convert_nmrglue
from make_parameter_file_cl import parameter_write_cl
from pulse_sequence_parsing import PulseSequenceParser
from config_register import registry
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


def interactive_plot_make(dic, data):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider, CheckButtons
    uc0 = ng.pipe.make_uc(dic, data, dim=1)
    uc1 = ng.pipe.make_uc(dic, data, dim=0)
    
    ppms_0 = uc0.ppm_scale()
    ppms_1 = uc1.ppm_scale()
    xx,yy = np.meshgrid(ppms_0,ppms_1)

    fig,ax = plt.subplots()
    one_d_proj = data[128]

    line, = ax.plot(ppms_0, np.real(one_d_proj))
    ax.invert_xaxis()

    ax_phase = plt.axes([0.15, 0.02, 0.65, 0.02])

    slider_phase = Slider(ax_phase, 'Phase', -180,180, valinit=0.0)
    
    def update(val):
        phase_val = slider_phase.val 
        adjusted_spec=one_d_proj*np.exp(1.0j*phase_val*np.pi/180.)
        line.set_ydata(np.real(adjusted_spec))
        #ax.relim()
        #ax.autoscale_view()
        fig.canvas.draw_idle() 
    
    
    slider_phase.on_changed(update)
    plt.show()

def write_parameter_file():
    pass

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np
    input_dat = FindingParameters()
    pp_parser = PulseSequenceParser()
    sequence = pp_parser.parse()

    # TODO: introduce possibility of an overwite so we can choose a specific config file
    config = registry.get_default_config(sequence)
    

    nmr_glue_conv = Convert_nmrglue(input_dat.params, input_dat)
    
    params = parameter_write_cl(nmr_glue_conv, config)
    params.write_out_dict(params.dictionary)
    
    # TODO: write out appropriate filename
    config.process_data('test.fid','test.ft2')
    
    import subprocess

    # Basic usage - waits for completion
    subprocess.Popen(
    ['python', '/Users/gogs/Bind/NMR/source/SpinExplorer/src/SpinExplorer/SpinView/SpinView.py'],
    stdin=subprocess.DEVNULL,
    stdout=subprocess.DEVNULL, 
    stderr=subprocess.DEVNULL,
    start_new_session=True)
    print('opened Spinview...')

    # interactive_plot_make(dic,data)
    # one_d_proj = np.sum(data, axis=0)
    # print(data.shape)
    
    # uc0 = ng.pipe.make_uc(dic, data, dim=1)
    # uc1 = ng.pipe.make_uc(dic, data, dim=0)
    
    # ppms_0 = uc0.ppm_scale()
    # ppms_1 = uc1.ppm_scale()
    # plt.plot(ppms_0,one_d_proj)
    # plt.gca().invert_xaxis()
    # plt.show()
    
    # print('hello')
    # print(data.shape)
    # cl = [np.max(data*0.1)*1.2**x for x in range(14)]
    # xx,yy = np.meshgrid(ppms_0,ppms_1)
    # fig,ax = plt.subplots()
    # print(xx.shape)
    # print(yy.shape)
    # ax.contour(xx,yy,data,cl)
    # ax.set_xlim(max(ppms_0), min(ppms_0))
    # ax.set_ylim(max(ppms_1), min(ppms_1))
    # plt.show()