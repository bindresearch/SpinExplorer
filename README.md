A python-based graphical interface package to perform NMR data processing and analysis.

Processing can be performed using nmrglue or nmrPipe. All processing except SMILE NUS reconstruction can be performed using nmrglue.

Video tutorials will be created shortly and made available at "https://www.youtube.com/@BindResearch".

Written documentation is available (SpinExplorerDocumentation.pdf).

Test NMR data is available in the downloadable TestData.zip folder.

Releases
--------
Each release will have application bundles for the app to make the installation as simple as possible. In the current pre-release, only MacOS (apple silicon and x86_64) architectures are supported. Windows on x86_64 architectures is also supported. For Linux users, follow the installation for GitHub instructions below (Linux users may also need to compile wxPython from source).


Opening the application
------------------------

- Download the compressed application folder (e.g. SpinExplorer_MacOS_AppleSilicon.zip) to your machine and extract the folder.
- Move the extracted downloaded folder to your applications folder
- For **windows** computers you may need to also download the visual C++ package (https://aka.ms/vc14/vc_redist.x64.exe). For further information visit https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist?view=msvc-170#latest-supported-redistributable-version
- Double click on the SpinExplorer application executable (within the downloaded SpinExplorer folder) to open the application. If you are opening the application for the first time, it may take a while to open.

Security
--------
On MacOS and Windows systems, it is possible that your computer may prevent the application from opening for security reasons. We are trying to become recognised as official Apple/Microsoft developers to prevent this issue for future releases.

To override these warning and open the application:
- **MacOS:** Go to "system settings" and click on "Privacy & Security". Scroll down to the "Security" section and click "allow anyway" on the SpinExplorer application.
- **Windows:** When you first open the application the following pop-out message may occur "Windows protected your PC. Microsoft Defender SmartScreen prevented an unrecognised app from starting. Running this app might put your PC at risk."
Click on more info and then click on “Run anyway” to open the application.


Build information
-----------------
The application was built on the following systems:

- SpinExplorer_MacOS_AppleSilicon - MacBook Pro (2024), M4 Pro, macOS Tahoe 26.3
- SpinExplorer_MacOS_x86_64 - MacBook Air (early 2015), 1.6GHz Dual-Core Intel i5, macOS Monterey 12.7.6
- SpinExplorer_Windows_AMD64 - HP Z2 Tower Workstation, Intel(R) Core(TM) i7-14700K 3.40GHz, Windows 11 Pro 23H2



Installation from GitHub
------------------------
If users prefer to view/edit the source code directly or add extra functionality, cloning the repository from GitHub is possible.

- Clone the GitHub repository
- Create a virtual environment containing python3>=3.11
- Activate the virtual environment
- From the main package directory (containing README.md) run the command "pip install ."

Once installed from GitHub, the command SpinExplorer can be run in a terminal to open the application in the same manner as the releases. In addition, the commands SpinConverter, SpinProcess and SpinView can be run from a terminal in a directory containing raw NMR data to perform NMR data conversion, processing, and viewing/analysis, respectively. This is useful if a user prefers to use the command line for NMR processing/analysis.





