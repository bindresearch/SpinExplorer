A python-based graphical interface package to perform NMR data processing and analysis.

Processing can be performed using nmrglue or nmrPipe. All processing except SMILE NUS reconstruction can be performed using nmrglue.

Video tutorials will be created shortly and made available at "https://www.youtube.com/@BindResearch".

Written documentation is available (SpinExplorerDocumentation.pdf)

Releases
--------
Each release will have application bundles for the app to make the installation as simple as possible. In the current pre-release, only MacOS (apple silicon and x86_64) architectures are supported. Windows on x86_64 architectures is also supported. For Linux users, follow the installation for GitHub instructions below (Linux users may also need to compile wxPython from source).


Installation from GitHub
------------------------
If users prefer to view/edit the source code directly or add extra functionality, cloning the repository from GitHub is possible.

- Clone the GitHub repository
- Create a virtual environment containing python3>=3.11
- Activate the virtual environment
- From the main package directory (containing README.md) run the command "pip install ."

Once installed from GitHub, the command SpinExplorer can be run in a terminal to open the application in the same manner as the releases. In addition, the commands SpinConverter, SpinProcess and SpinView can be ran from a terminal in a directory containing raw NMR data to perform NMR data conversion, processing, and viewing/analysis, respectively. This is useful if a user prefers to use the command line for NMR processing/analysis.





