A python-based graphical interface package to perform NMR data processing and analysis.

Full processing functionality is enabled through nmrPipe installation: https://www.ibbr.umd.edu/nmrpipe/install.html  
In the absence of nmrPipe (e.g. windows systems), NMR processing is performed using nmrglue.

Full documentation provided at "https://github.com/james-eaton-1/SpinExplorer" in Documentation.pdf

Once installed, the commands SpinConverter, SpinProcess and SpinView can be ran from a terminal in a directory containing raw NMR data to perform NMR data conversion, processing, and viewing/analysis, respectively. 

Installation (macOS/Linux):
- Ensure python3 version is greater than 3.10 (package not tested on python<3.10)
- It is recommended to create a python3 virtual environment in the terminal using the command:  
"python3 -m venv ~/SpinEnv"    
More details at https://docs.python.org/3/tutorial/venv.html
- Alternatively conda environments can be used if preferred
- Activate the virtual environment using the command:  
"source ~/SpinEnv/bin/activate"  
This will need to be ran every time a new terminal is created, or the following line can be
added to the ~/.zshrc or ~/.bashrc file: "source ~/SpinEnv/bin/activate"
- Install the SpinExplorer package and all dependencies using:
"python3 -m pip install SpinExplorer"

Note: the wxPython dependency might not install correctly on Linux systems.  
More details and fixes at https://wxpython.org/pages/downloads/  
Currently, wxPython needs to be build from source for Linux systems without x86_64 architecture.  
More details at https://wxpython.org/blog/2017-08-17-builds-for-linux-with-pip/index.html

Installation (Windows):
- Using py (py-launcher) is recommended for windows systems  
If py is not present, download latest python from https://www.python.org, and ensure install py-launcher is ticked
- If desired, a virtual environment for the package can be created according to https://docs.python.org/3/tutorial/venv.html
- Install the SpinExplorer package and all dependencies using the following command in command prompt/power shell:  
"py -m pip install SpinExplorer"




