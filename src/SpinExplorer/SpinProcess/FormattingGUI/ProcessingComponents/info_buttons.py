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
import wx.lib.agw.hyperlink as hl


class InfoButtons:

    def __init__(self, app):
        """
        This class contains pop out messages containing information
        for all the different NMR processing options.
        """

        self.app = app

        # Popout text colour
        self.colour = "BLUE"

    def on_solvent_suppression_info(self, event):
        """
        Include information on how the low-pass, splines
        and polynomial solvent suppression works. In nmrglue
        only low-pass filter solvent suppression is implemented
        """

        general_information = """Aqueous NMR spectra often suffer from
        substantial solvent signals reducing spectral quality. \n If 
        NMR spectra have been run with the solvent set at the carrier
          frequency there will low frequency oscillations in the FID\n 
          (in the rotating frame) corresponding to the solvent signal. 
          These oscillations can be removed by applying a low-pass filter
            to the FID."""
        general_information = " ".join(general_information.split()) + "\n\n"

        lowpass_filter = """Low-pass filter: \n A low-pass filter is a filter
          that passes signals with a frequency lower than a certain cutoff 
          frequency and attenuates signals with frequencies higher than the
            cutoff frequency."""
        lowpass_filter = " ".join(lowpass_filter.split()) + "\n\n"

        lowpass_filter_shape = """Low-pass filter shape: \n The shape of 
        the low-pass filter can be set to a boxcar, sine, or sine squared
          shape."""
        lowpass_filter_shape = " ".join(lowpass_filter_shape.split()) + "\n\n"

        additional_advanced_options = """Additional advanced options using
          splines or polynomial digital solvent suppression can also be 
          selected. Further advanced options such as filter length (-fl) 
          can be added manually to the nmrproc.com file to further refine
            the solvent suppression."""
        additional_advanced_options = (
            " ".join(additional_advanced_options.split()) + "\n\n"
        )

        # Create a popup window with the information
        self.solvent_suppression_info_window = wx.Frame(
            self.app, -1, "Solvent Suppression Information", size=(450, 600)
        )

        self.solvent_suppression_info_window_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.solvent_suppression_info_window.SetSizer(
            self.solvent_suppression_info_window_sizer
        )

        self.solvent_suppression_info_window_sizer.AddSpacer(10)

        # Create a sizer for the solvent suppression information
        self.solvent_suppression_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.solvent_suppression_info_sizer.AddSpacer(10)

        self.solvent_suppression_info_window_sizer.AddSpacer(10)

        # Have a text control for the information
        self.solvent_suppression_info_textcontrol = wx.StaticText(
            self.solvent_suppression_info_window,
            label=general_information
            + lowpass_filter
            + lowpass_filter_shape
            + additional_advanced_options,
            size=(400, 500),
        )

        self.solvent_suppression_info_sizer.Add(
            self.solvent_suppression_info_textcontrol, 0, wx.ALIGN_CENTER_HORIZONTAL
        )

        # Add a hyperlink to the sizer for the NMRPipe SOL help page
        self.solvent_suppression_info_sizer.AddSpacer(10)

        self.sol_hyperlink = hl.HyperLinkCtrl(
            self.solvent_suppression_info_window,
            -1,
            "NMRPipe Solvent Suppression Help Page",
            URL="http://www.nmrscience.com/ref/nmrpipe/sol.html",
        )

        self.sol_hyperlink.SetColours(self.colour, self.colour, self.colour)
        self.sol_hyperlink.SetUnderlines(False, False, False)
        self.sol_hyperlink.UpdateLink()

        self.solvent_suppression_info_sizer.Add(
            self.sol_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL
        )

        self.solvent_suppression_info_sizer.AddSpacer(10)

        self.solvent_suppression_info_window_sizer.Add(
            self.solvent_suppression_info_sizer, 0, wx.ALIGN_CENTER_VERTICAL
        )
        self.solvent_suppression_info_sizer.AddSpacer(10)
        self.solvent_suppression_info_window.SetSizer(
            self.solvent_suppression_info_window_sizer
        )
        self.solvent_suppression_info_window.Show()

    def on_linear_prediction_info(self, event):
        """
        Create a popout window with information about linear prediction
        """

        # Create a new frame
        self.linear_prediction_info_frame = wx.Frame(
            self.app, -1, "Linear Prediction Information", size=(500, 300)
        )

        # Create a sizer to hold the box
        self.linear_prediction_info_sizer_window = wx.BoxSizer(wx.VERTICAL)
        self.linear_prediction_info_sizer_window.AddSpacer(10)

        # Create a sizer to hold the text
        self.linear_prediction_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.linear_prediction_info_sizer.AddSpacer(10)

        # Create a text box with the information
        # Linear prediction information
        linear_prediction_information = """Linear prediction is a method 
        used to increase the resolution of NMR spectra. It is used to predict
        the points of truncated FIDs (especially in indirect dimensions) and
        increase signal resolution.\n\n The linear prediction coefficients 
        can be predicted using the forward FID data, backward data or an 
        average of both directions. Then these can be used to add predicted 
        points either before or after the current FID.\n\n Note that 
        advanced options such as  -pred (number of predicted points) and
        -ord (number of predicted coefficients) can be implemented by 
        manually added them to the nmrproc.com file."""

        linear_prediction_information = (
            " ".join(linear_prediction_information.split()) + "\n"
        )

        self.linear_prediction_info_text = wx.StaticText(
            self.linear_prediction_info_frame,
            -1,
            linear_prediction_information,
            size=(450, 200),
            style=wx.ALIGN_CENTER_HORIZONTAL,
        )

        # Add the text to the sizer
        self.linear_prediction_info_sizer.Add(
            self.linear_prediction_info_text, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.linear_prediction_info_sizer.AddSpacer(10)

        # Add a url to the nmrPipe help page
        url = "http://www.nmrscience.com/ref/nmrpipe/lp.html"
        self.linear_prediction_info_url = hl.HyperLinkCtrl(
            self.linear_prediction_info_frame,
            -1,
            "NMRPipe Help Page for Linear Prediction",
            URL=url,
        )
        self.linear_prediction_info_url.SetColours(
            self.colour, self.colour, self.colour
        )
        self.linear_prediction_info_url.SetUnderlines(False, False, False)
        self.linear_prediction_info_url.UpdateLink()

        # Add url to the sizer
        self.linear_prediction_info_sizer.Add(
            self.linear_prediction_info_url, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.linear_prediction_info_sizer.AddSpacer(10)

        # Add the sizer to the window sizer
        self.linear_prediction_info_sizer_window.Add(
            self.linear_prediction_info_sizer, 0, wx.ALIGN_CENTER
        )
        self.linear_prediction_info_sizer_window.AddSpacer(10)

        # Add the window sizer to the frame
        self.linear_prediction_info_frame.SetSizer(
            self.linear_prediction_info_sizer_window
        )

        # Show the frame
        self.linear_prediction_info_frame.Show()

    def on_apodization_info(self, event):
        """
        Include information on how the apodization functions work
        and when might be best to use each one.
        """

        # Make a pop out window with all of the apodization information
        # Create a new frame
        self.apodization_info_frame = wx.Frame(
            self.app, -1, "Apodization Information", size=(800, 700)
        )

        # Create a sizer to hold the box
        self.apodization_info_sizer_window = wx.BoxSizer(wx.VERTICAL)
        self.apodization_info_sizer_window.AddSpacer(10)

        general_information = """Apodization (window functions) can be 
        multiplied to the FID to increase signal to noise/increase resolution. 
        Furthermore, window functions can be applied to data whose FID has 
        been truncated (often in the indirect dimension) in order to reduce 
        the presence of Sinc wiggles in the NMR spectrum. Below is an 
        explanation of common apodization functions\n and when to use them."""

        general_information = " ".join(general_information.split())

        first_point_information = """First point scaling: In data where there 
        is no first order (p1) phase correction required a value of 0.5 should 
        be used, otherwise a value of 1.0 should be used."""

        first_point_information = " ".join(first_point_information.split())

        exponential_information = """Exponential: The exponential window 
        function is used to apply an exponential decay to the FID. This 
        suppresses the noise at the end of an FID enhancing signal to noise, 
        but leads to a reduction in resolution. The line broadening term 
        dictates the strength of the decay of the exponential function applied 
        to the FID."""

        exponential_information = " ".join(exponential_information.split())
        exponential_info_hyperlink = hl.HyperLinkCtrl(
            self.apodization_info_frame,
            -1,
            "Further information on the nmrPipe exponential window function",
            URL="http://www.nmrscience.com/ref/nmrpipe/em.html",
        )
        exponential_info_hyperlink.SetColours(self.colour, self.colour, self.colour)
        exponential_info_hyperlink.SetUnderlines(False, False, False)
        exponential_info_hyperlink.UpdateLink()

        lorentz_to_gauss_information = """Lorentz to Gauss: The Lorentz to 
        Gauss window function is used to apply a Lorentzian to Gaussian 
        transformation to the FID. The inverse Lorentzian parameter can be 
        tuned to remove all Lorentzian character from peaks in the NMR spectrum 
        leaving them with a pure Gaussian peak shape. This can enhance signal 
        resolution and is advantageous for peak picking routines which are often 
        more accurate with Gaussian peak shapes. The Gaussian broadening value can 
        be increased further to enhance the signal to noise at the cost of 
        reduced resolution. As a rule of thumb, the Gaussian broadening parameter 
        should be 3x larger than the inverse Lorentzian value."""

        lorentz_to_gauss_information = (
            " ".join(lorentz_to_gauss_information.split()) + "\n"
        )

        lorentz_to_gauss_info_hyperlink = hl.HyperLinkCtrl(
            self.apodization_info_frame,
            -1,
            "Further information on the nmrPipe Lorentz to Gauss window function",
            URL="http://www.nmrscience.com/ref/nmrpipe/gm.html",
        )
        lorentz_to_gauss_info_hyperlink.SetColours(
            self.colour, self.colour, self.colour
        )
        lorentz_to_gauss_info_hyperlink.SetUnderlines(False, False, False)
        lorentz_to_gauss_info_hyperlink.UpdateLink()

        sinebell_information = """Sinebell: \n The Sinebell window function 
        is used to apply a sinebell function to the FID. The offset value 
        adjusts the phase of the sinebell function, the end value adjusts the 
        end of the sinebell function, and the power value adjusts the power of 
        the sinebell function."""

        sinebell_information = " ".join(sinebell_information.split()) + "\n"

        sinebell_info_hyperlink = hl.HyperLinkCtrl(
            self.apodization_info_frame,
            -1,
            "Further information on the nmrPipe Sinebell window function",
            URL="http://www.nmrscience.com/ref/nmrpipe/sp.html",
        )
        sinebell_info_hyperlink.SetColours(self.colour, self.colour, self.colour)
        sinebell_info_hyperlink.SetUnderlines(False, False, False)
        sinebell_info_hyperlink.UpdateLink()

        gauss_broadening_information = """Gauss broadening: \n The Gauss 
        broadening window function is used to apply a Decaying Exponential 
        / Gaussian function to the FID. The line broadening value adjusts 
        the decay of the exponential function, whilst the gaussian broadening 
        term adjusts the decay of the gaussian function."""

        gauss_broadening_information = (
            " ".join(gauss_broadening_information.split()) + "\n"
        )

        gauss_broadening_hyperlink = hl.HyperLinkCtrl(
            self.apodization_info_frame,
            -1,
            "Further information on the nmrPipe Gaussian broadening window function",
            URL="http://www.nmrscience.com/ref/nmrpipe/gmb.html",
        )
        gauss_broadening_hyperlink.SetColours(self.colour, self.colour, self.colour)
        gauss_broadening_hyperlink.SetUnderlines(False, False, False)
        gauss_broadening_hyperlink.UpdateLink()

        trapeziod_broadening_information = """Trapezoid: A trapezoid function 
        is applied to the FID. The number of ramp up and ramp down points can 
        be adjusted to change the shape of the trapezoid"""

        trapeziod_broadening_information = " ".join(
            trapeziod_broadening_information.split()
        )

        trapezoid_hyperlink = hl.HyperLinkCtrl(
            self.apodization_info_frame,
            -1,
            "Further information on the nmrPipe trapezoid window function",
            URL="http://www.nmrscience.com/ref/nmrpipe/tm.html",
        )
        trapezoid_hyperlink.SetColours(self.colour, self.colour, self.colour)
        trapezoid_hyperlink.SetUnderlines(False, False, False)
        trapezoid_hyperlink.UpdateLink()

        triangle_information = """Triangle: A triangle function is applied 
        to the FID. The location of the maximum of the triangle can be adjusted
          between 0 (first point) and 1 (last point)."""
        triangle_information = " ".join(triangle_information.split())

        triangle_hyperlink = hl.HyperLinkCtrl(
            self.apodization_info_frame,
            -1,
            "Further information on the nmrPipe triangle window function",
            URL="http://www.nmrscience.com/ref/nmrpipe/tri.html",
        )
        triangle_hyperlink.SetColours(self.colour, self.colour, self.colour)
        triangle_hyperlink.SetUnderlines(False, False, False)
        triangle_hyperlink.UpdateLink()

        # Create a sizer to hold the text
        self.apodization_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.apodization_info_sizer.AddSpacer(10)

        self.apodization_info_text1 = wx.StaticText(
            self.apodization_info_frame,
            -1,
            general_information,
        )
        self.apodization_info_text1.Wrap(700)
        self.apodization_info_text2 = wx.StaticText(
            self.apodization_info_frame,
            -1,
            first_point_information,
        )
        self.apodization_info_text2.Wrap(700)
        self.apodization_info_text3 = wx.StaticText(
            self.apodization_info_frame,
            -1,
            exponential_information,
        )
        self.apodization_info_sizer.Add(
            self.apodization_info_text1, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.apodization_info_sizer.AddSpacer(10)
        self.apodization_info_sizer.Add(
            self.apodization_info_text2, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.apodization_info_sizer.AddSpacer(10)
        self.apodization_info_sizer.Add(
            self.apodization_info_text3, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.apodization_info_text3.Wrap(700)
        self.apodization_info_sizer.Add(
            exponential_info_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.apodization_info_sizer.AddSpacer(10)
        text_lorentz_to_gauss = wx.StaticText(
            self.apodization_info_frame,
            -1,
            lorentz_to_gauss_information,
        )
        text_lorentz_to_gauss.Wrap(700)

        self.apodization_info_sizer.Add(
            text_lorentz_to_gauss, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.apodization_info_sizer.Add(
            lorentz_to_gauss_info_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.apodization_info_sizer.AddSpacer(10)
        text_sinebell = wx.StaticText(
            self.apodization_info_frame,
            -1,
            sinebell_information,
        )
        text_sinebell.Wrap(700)
        self.apodization_info_sizer.Add(text_sinebell, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.apodization_info_sizer.Add(
            sinebell_info_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.apodization_info_sizer.AddSpacer(10)
        text_gauss_broadening = wx.StaticText(
            self.apodization_info_frame,
            -1,
            gauss_broadening_information,
        )
        text_gauss_broadening.Wrap(700)
        self.apodization_info_sizer.Add(
            text_gauss_broadening,
            0,
        )
        self.apodization_info_sizer.Add(
            gauss_broadening_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.apodization_info_sizer.AddSpacer(10)
        text_trapezoid = wx.StaticText(
            self.apodization_info_frame,
            -1,
            trapeziod_broadening_information,
        )
        text_trapezoid.Wrap(700)
        self.apodization_info_sizer.Add(text_trapezoid, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.apodization_info_sizer.Add(
            trapezoid_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.apodization_info_sizer.AddSpacer(10)
        text_triangle = wx.StaticText(
            self.apodization_info_frame,
            -1,
            triangle_information,
        )
        text_triangle.Wrap(700)
        self.apodization_info_sizer.Add(text_triangle, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.apodization_info_sizer.Add(
            triangle_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.apodization_info_sizer.AddSpacer(10)

        # Add the sizer to the window sizer
        self.apodization_info_sizer_window.Add(
            self.apodization_info_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.apodization_info_sizer_window.AddSpacer(10)

        # Add the window sizer to the frame
        self.apodization_info_frame.SetSizer(self.apodization_info_sizer_window)

        # Show the frame
        self.apodization_info_frame.Center()
        self.apodization_info_frame.Show()

    def on_zero_fill_info(self, event):
        """
        Create a popout window with information about zero filling
        """

        # Create a new frame
        self.zero_fill_info_frame = wx.Frame(
            self.app, -1, "Zero Filling Information", size=(500, 300)
        )

        # Create a sizer to hold the box
        self.zero_fill_info_sizer_window = wx.BoxSizer(wx.VERTICAL)
        self.zero_fill_info_sizer_window.AddSpacer(10)

        # Create a sizer to hold the text
        self.zero_fill_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.zero_fill_info_sizer.AddSpacer(10)

        # Create a text box with the information
        # Zero filling information
        zero_fill_information = """Zero filling is a method used to 
        increase the resolution of NMR spectra. It is used to add 
        zeros to the end of the FID to increase the number of points.\n\n 
        It is important that the size of the data is at least doubled to 
        prevent loss of resolution when the imaginary component of the complex
        data is deleted. In addition it is advised that the data is rounded 
        to the nearest power of 2 in order to speed up the Fast Fourier 
        Transform process.\n\nFurther advanced zero filling options can be 
        added manually to the nmrproc.com file if performing nmrPipe
        processing."""

        zero_fill_information = " ".join(zero_fill_information.split())

        self.zero_fill_info_text = wx.StaticText(
            self.zero_fill_info_frame,
            -1,
            zero_fill_information,
            size=(450, 150),
            style=wx.ALIGN_CENTER,
        )

        # Add the text to the sizer
        self.zero_fill_info_sizer.Add(
            self.zero_fill_info_text, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.zero_fill_info_sizer.AddSpacer(10)

        # Have a url to the nmrPipe help page for zero filling
        url = "http://www.nmrscience.com/ref/nmrpipe/zf.html"
        self.zero_fill_info_url = hl.HyperLinkCtrl(
            self.zero_fill_info_frame, -1, "NMRPipe Help Page for Zero Filling", URL=url
        )
        self.zero_fill_info_url.SetColours(self.colour, self.colour, self.colour)
        self.zero_fill_info_url.SetUnderlines(False, False, False)
        self.zero_fill_info_url.UpdateLink()

        # Add the url to the sizer
        self.zero_fill_info_sizer.Add(
            self.zero_fill_info_url, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.zero_fill_info_sizer.AddSpacer(10)

        # Add the sizer to the window sizer
        self.zero_fill_info_sizer_window.Add(
            self.zero_fill_info_sizer, 0, wx.ALIGN_CENTER
        )
        self.zero_fill_info_sizer_window.AddSpacer(10)

        # Add the window sizer to the frame
        self.zero_fill_info_frame.SetSizer(self.zero_fill_info_sizer_window)

        # Show the frame
        self.zero_fill_info_frame.Show()

    def on_fourier_transform_info(self, event):
        """
        Creating a popout showing fourier transform information
        """
        ft_text = """The fourier transform applies a complex fourier 
        transform to the FID to convert it to a frequency domain spectrum.
        \n Further information can be found using the link below."""

        ft_text = " ".join(ft_text.split())

        # Create a popup window with the information
        self.fourier_transform_info_window = wx.Frame(
            self.app, -1, "Fourier Transform Information", size=(450, 150)
        )

        self.fourier_transform_info_window_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.fourier_transform_info_window.SetSizer(
            self.fourier_transform_info_window_sizer
        )

        self.fourier_transform_info_window_sizer.AddSpacer(10)

        # Create a sizer for the fourier transform information
        self.fourier_transform_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.fourier_transform_info_sizer.AddSpacer(10)
        self.fourier_transform_info_sizer.Add(
            wx.StaticText(
                self.fourier_transform_info_window, -1, ft_text, size=(400, 50)
            ),
            0,
            wx.ALIGN_CENTER,
        )
        self.fourier_transform_info_sizer.AddSpacer(10)

        # Have a hyperlink to the fourier transform information
        self.fourier_transform_info_hyperlink = hl.HyperLinkCtrl(
            self.fourier_transform_info_window,
            -1,
            "NMRPipe Help Page for Fourier Transform",
            URL="http://www.nmrscience.com/ref/nmrpipe/ft.html",
        )
        self.fourier_transform_info_hyperlink.SetColours(
            self.colour, self.colour, self.colour
        )
        self.fourier_transform_info_hyperlink.SetUnderlines(False, False, False)
        self.fourier_transform_info_hyperlink.SetBold(False)
        self.fourier_transform_info_hyperlink.UpdateLink()
        self.fourier_transform_info_sizer.Add(
            self.fourier_transform_info_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.fourier_transform_info_sizer.AddSpacer(10)

        self.fourier_transform_info_window_sizer.Add(
            self.fourier_transform_info_sizer, 0, wx.ALIGN_CENTER
        )

        self.fourier_transform_info_window.Show()

    def on_phase_correction_info(self, event):
        """
        Creating a popout showing phasing information for the direct
        dimension.
        """
        phase_correction_text = """Phase correction is a method to correct
          for phase errors in the FID. Zero order phase correction (p0) 
          is used to correct a phase offset that is applied equally across 
          the spectrum. However, a first order phase correction (p1) is used
            to correct the phasing in a spectrum where peaks in different 
            locations of the spectrum require a different phasing value. \n 
            Further information can be found using the link below."""

        phase_correction_text = " ".join(phase_correction_text.split())

        # Create a popup window with the information
        self.phase_correction_info_window = wx.Frame(
            self.app, -1, "Phase Correction Information", size=(450, 200)
        )
        self.phase_correction_info_window_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.phase_correction_info_window.SetSizer(
            self.phase_correction_info_window_sizer
        )
        self.phase_correction_info_window_sizer.AddSpacer(10)
        # Create a sizer for the phase correction information
        self.phase_correction_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.phase_correction_info_sizer.AddSpacer(10)
        self.phase_correction_info_sizer.Add(
            wx.StaticText(
                self.phase_correction_info_window,
                -1,
                phase_correction_text,
                size=(400, 100),
            ),
            0,
            wx.ALIGN_CENTER,
        )
        self.phase_correction_info_sizer.AddSpacer(10)
        # Have a hyperlink to the phase correction information
        self.phase_correction_info_hyperlink = hl.HyperLinkCtrl(
            self.phase_correction_info_window,
            -1,
            "NMRPipe Help Page for Phase Correction",
            URL="http://www.nmrscience.com/ref/nmrpipe/ps.html",
        )
        self.phase_correction_info_hyperlink.SetColours(
            self.colour, self.colour, self.colour
        )
        self.phase_correction_info_hyperlink.SetUnderlines(False, False, False)
        self.phase_correction_info_hyperlink.SetBold(False)
        self.phase_correction_info_hyperlink.UpdateLink()
        self.phase_correction_info_sizer.Add(
            self.phase_correction_info_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.phase_correction_info_sizer.AddSpacer(10)
        self.phase_correction_info_window_sizer.Add(
            self.phase_correction_info_sizer, 0, wx.ALIGN_CENTER
        )
        self.phase_correction_info_window.Show()

    def on_extraction_info(self, event):
        """
        Creating a popout with information on data extraction
        """
        extraction_text = """Extraction of data between two chemical 
        shift values can be used to extract a region of interest from 
        the spectrum. This can be useful for removing solvent signals 
        or other unwanted peaks. \n Further information can be found 
        using the link below."""

        extraction_text = " ".join(extraction_text.split())

        # Create a popup window with the information
        self.extraction_info_window = wx.Frame(
            self.app, -1, "Extraction Information", size=(450, 200)
        )

        self.extraction_info_window_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.extraction_info_window.SetSizer(self.extraction_info_window_sizer)
        self.extraction_info_window_sizer.AddSpacer(10)
        # Create a sizer for the extraction information
        self.extraction_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.extraction_info_sizer.AddSpacer(10)
        self.extraction_info_sizer.Add(
            wx.StaticText(
                self.extraction_info_window, -1, extraction_text, size=(400, 100)
            ),
            0,
            wx.ALIGN_CENTER,
        )
        self.extraction_info_sizer.AddSpacer(10)
        # Have a hyperlink to the extraction information
        self.extraction_info_hyperlink = hl.HyperLinkCtrl(
            self.extraction_info_window,
            -1,
            "NMRPipe Help Page for Extraction",
            URL="http://www.nmrscience.com/ref/nmrpipe/ext.html",
        )
        self.extraction_info_hyperlink.SetColours(self.colour, self.colour, self.colour)
        self.extraction_info_hyperlink.SetUnderlines(False, False, False)
        self.extraction_info_hyperlink.SetBold(False)
        self.extraction_info_hyperlink.UpdateLink()
        self.extraction_info_sizer.Add(
            self.extraction_info_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.extraction_info_sizer.AddSpacer(10)
        self.extraction_info_window_sizer.Add(
            self.extraction_info_sizer, 0, wx.ALIGN_CENTER
        )
        self.extraction_info_window.Show()

    def on_baseline_correction_info(self, event):
        """
        Creating a popout on baseline correction information
        """
        # Include information on linear and polynomial baseline correction
        linear_information = """Linear baseline correction: \n Linear baseline 
        correction is a method to correct for a linear baseline issue in the 
        spectrum. The linear baseline is removed by fitting a straight line to 
        the spectrum and subtracting it from the spectrum."""

        linear_information = " ".join(linear_information.split()) + "\n\n"

        polynomial_information = "Polynomial baseline correction: \n Polynomial baseline correction is a method to correct for a polynomial baseline issue in the spectrum. The polynomial baseline is removed by fitting a polynomial to the spectrum and subtracting it from the spectrum. \n\n"
        extra_information = "The node list is a list of percentages that are used to define the nodes (points which are expected to have 0 intensity) for the baseline correction. The node width is the number of points used to define the nodes. Further advanced options can be added to the processing file nmrproc.com file manually \n\n"

        # Create a popup window with the information
        self.baseline_correction_info_window = wx.Frame(
            self.app, -1, "Baseline Correction Information", size=(450, 450)
        )
        self.baseline_correction_info_window_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.baseline_correction_info_window.SetSizer(
            self.baseline_correction_info_window_sizer
        )
        self.baseline_correction_info_window_sizer.AddSpacer(10)
        # Create a sizer for the baseline correction information
        self.baseline_correction_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.baseline_correction_info_sizer.AddSpacer(10)
        self.baseline_correction_info_sizer.Add(
            wx.StaticText(
                self.baseline_correction_info_window,
                -1,
                linear_information + polynomial_information + extra_information,
                size=(400, 300),
            ),
            0,
            wx.ALIGN_CENTER,
        )
        self.baseline_correction_info_sizer.AddSpacer(10)
        # Have a hyperlink to the linear baseline correction information
        self.baseline_correction_info_hyperlink = hl.HyperLinkCtrl(
            self.baseline_correction_info_window,
            -1,
            "NMRPipe Help Page for Linear Baseline Correction",
            URL="http://www.nmrscience.com/ref/nmrpipe/base.html",
        )
        self.baseline_correction_info_hyperlink.SetColours(
            self.colour, self.colour, self.colour
        )
        self.baseline_correction_info_hyperlink.SetUnderlines(False, False, False)
        self.baseline_correction_info_hyperlink.SetBold(False)
        self.baseline_correction_info_hyperlink.UpdateLink()
        # Have a hyperlink to the polynomial baseline correction information
        self.baseline_correction_info_hyperlink_2 = hl.HyperLinkCtrl(
            self.baseline_correction_info_window,
            -1,
            "NMRPipe Help Page for Polynomial Baseline Correction",
            URL="http://www.nmrscience.com/ref/nmrpipe/poly.html",
        )
        self.baseline_correction_info_hyperlink_2.SetColours(
            self.colour, self.colour, self.colour
        )
        self.baseline_correction_info_hyperlink_2.SetUnderlines(False, False, False)
        self.baseline_correction_info_hyperlink_2.SetBold(False)
        self.baseline_correction_info_hyperlink_2.UpdateLink()
        self.baseline_correction_info_sizer.Add(
            self.baseline_correction_info_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.baseline_correction_info_sizer.AddSpacer(10)
        self.baseline_correction_info_sizer.Add(
            self.baseline_correction_info_hyperlink_2, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.baseline_correction_info_sizer.AddSpacer(10)
        self.baseline_correction_info_window_sizer.Add(
            self.baseline_correction_info_sizer, 0, wx.ALIGN_CENTER
        )
        self.baseline_correction_info_window.Show()

    def on_linear_prediction_info_indirect(self, event):
        """
        Create a popout window with information about linear prediction
        and non-uniform sampling (NUS)
        """
        # Create a new frame
        self.linear_prediction_info_frame = wx.Frame(
            self.app, -1, "Linear Prediction / SMILE Information", size=(500, 500)
        )

        # Create a sizer to hold the box
        self.linear_prediction_info_sizer_window = wx.BoxSizer(wx.VERTICAL)
        self.linear_prediction_info_sizer_window.AddSpacer(10)

        # Create a sizer to hold the text
        self.linear_prediction_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.linear_prediction_info_sizer.AddSpacer(10)

        # Create a text box with the information
        # Linear prediction information
        linear_prediction_information = """Linear prediction is a method 
        used to increase the resolution of NMR spectra. It is used to predict 
        the points of truncated FIDs (especially in indirect dimensions) and 
        increase signal resolution.\n\n The linear prediction coefficients can 
        be predicted using the forward FID data, backward data or an average of 
        both directions. Then these can be used to add predicted points either 
        before or after the current FID.\n\n Note that advanced options such as 
          -pred (number of predicted points) and -ord (number of predicted 
          coefficients) can be implemented by manually added them to the 
          nmrproc.com file for nmrPipe processing."""

        linear_prediction_information = " ".join(linear_prediction_information.split())

        self.linear_prediction_info_text = wx.StaticText(
            self.linear_prediction_info_frame,
            -1,
            linear_prediction_information,
            size=(450, 200),
            style=wx.ALIGN_CENTER_HORIZONTAL,
        )

        # Add the text to the sizer
        self.linear_prediction_info_sizer.Add(
            self.linear_prediction_info_text, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.linear_prediction_info_sizer.AddSpacer(10)

        # Add a url to the nmrPipe help page
        url = "http://www.nmrscience.com/ref/nmrpipe/lp.html"
        self.linear_prediction_info_url = hl.HyperLinkCtrl(
            self.linear_prediction_info_frame,
            -1,
            "NMRPipe Help Page for Linear Prediction",
            URL=url,
        )
        self.linear_prediction_info_url.SetColours(
            self.colour, self.colour, self.colour
        )
        self.linear_prediction_info_url.SetUnderlines(False, False, False)
        self.linear_prediction_info_url.UpdateLink()

        # Add url to the sizer
        self.linear_prediction_info_sizer.Add(
            self.linear_prediction_info_url, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.linear_prediction_info_sizer.AddSpacer(10)

        # Add the sizer to the window sizer
        self.linear_prediction_info_sizer_window.Add(
            self.linear_prediction_info_sizer, 0, wx.ALIGN_CENTER
        )
        self.linear_prediction_info_sizer_window.AddSpacer(10)

        # Have text to explain SMILE NUS reconstruction
        smile_nus_text = """SMILE NUS reconstruction is a method 
        used to reconstruct non-uniformly sampled data. The NUS 
        file is a list of points that have been sampled in the FID.
        \nThe number of CPU's (default=1) is the number of cores 
        that will be used to perform the reconstruction and the 
        number of iterations can be changed to improve the accuracy 
        (default=800).\n Furthermore, in order for accurate SMILE 
        reconstruction, the correct zero (p0) and first (p1) order 
        phase correction values need to be inputted."""

        smile_nus_text = " ".join(smile_nus_text.split())

        self.smile_nus_text = wx.StaticText(
            self.linear_prediction_info_frame,
            -1,
            smile_nus_text,
            size=(450, 200),
            style=wx.ALIGN_CENTER_HORIZONTAL,
        )
        self.linear_prediction_info_sizer_window.Add(
            self.smile_nus_text, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.linear_prediction_info_sizer_window.AddSpacer(10)

        # Add the window sizer to the frame
        self.linear_prediction_info_frame.SetSizer(
            self.linear_prediction_info_sizer_window
        )

        # Show the frame
        self.linear_prediction_info_frame.Show()

    def on_phase_correction_info_indirect(self, event):
        phase_correction_text = """Phase correction is a method to correct
          for phase errors in the FID. Zero order phase correction (p0) 
          is used to correct a phase offset that is applied equally across
            the spectrum. However, a first order phase correction (p1) is 
            used to correct the phasing in a spectrum where peaks in different
              locations of the spectrum require a different phasing value. 
              For the indirect dimension, it is often the case that the 
              acquisition is delayed by an exact time so that the resulting 
              spectrum can be phased using the phase values of: p0=-90, p1=180. 
              This is often termed F1180. \n Further information can be found 
              using the link below."""

        phase_correction_text = " ".join(phase_correction_text.split())

        # Create a popup window with the information
        self.phase_correction_info_window = wx.Frame(
            self.app, -1, "Phase Correction Information", size=(450, 300)
        )

        self.phase_correction_info_window_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.phase_correction_info_window.SetSizer(
            self.phase_correction_info_window_sizer
        )

        self.phase_correction_info_window_sizer.AddSpacer(10)

        # Create a sizer for the phase correction information
        self.phase_correction_info_sizer = wx.BoxSizer(wx.VERTICAL)
        self.phase_correction_info_sizer.AddSpacer(10)
        self.phase_correction_info_sizer.Add(
            wx.StaticText(
                self.phase_correction_info_window,
                -1,
                phase_correction_text,
                size=(400, 200),
            ),
            0,
            wx.ALIGN_CENTER,
        )
        self.phase_correction_info_sizer.AddSpacer(10)

        # Have a hyperlink to the phase correction information
        self.phase_correction_info_hyperlink = hl.HyperLinkCtrl(
            self.phase_correction_info_window,
            -1,
            "NMRPipe Help Page for Phase Correction",
            URL="http://www.nmrscience.com/ref/nmrpipe/ps.html",
        )
        self.phase_correction_info_hyperlink.SetColours(
            self.colour, self.colour, self.colour
        )
        self.phase_correction_info_hyperlink.SetUnderlines(False, False, False)
        self.phase_correction_info_hyperlink.SetBold(False)
        self.phase_correction_info_hyperlink.UpdateLink()
        self.phase_correction_info_sizer.Add(
            self.phase_correction_info_hyperlink, 0, wx.ALIGN_CENTER_HORIZONTAL
        )
        self.phase_correction_info_sizer.AddSpacer(10)

        self.phase_correction_info_window_sizer.Add(
            self.phase_correction_info_sizer, 0, wx.ALIGN_CENTER
        )

        self.phase_correction_info_window.Show()
