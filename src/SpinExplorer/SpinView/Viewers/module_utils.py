import wx # type: ignore
import numpy as np

class InputROI(wx.Frame):
    def __init__(self, title, parent):
        self.main_frame = parent
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        width = 400
        height = 200
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.panel_ROI_input = wx.Panel(self, -1)
        self.main_ROI_input = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_ROI_input)


        self.make_manual_ROI_input_sizer()
        self.Show()

    def make_manual_ROI_input_sizer(self):
        self.input_ROI_label = wx.StaticBox(
            self, -1, "Input region of interest chemical shift range (ppm):"
        )
        self.input_ROI_sizer = wx.StaticBoxSizer(
            self.input_ROI_label, wx.VERTICAL
        )

        self.row = wx.BoxSizer(wx.HORIZONTAL)

        self.min_text = wx.StaticText(self.input_ROI_label,-1,'Min:')

        self.min_box = wx.TextCtrl(
            self.input_ROI_label, -1, value='', size=(100,20), 
        )

        self.max_text = wx.StaticText(self.input_ROI_label,-1,'Max:')

        self.max_box = wx.TextCtrl(
            self.input_ROI_label, -1, value='', size=(100,20), 
        )

        self.row.AddSpacer(10)
        self.row.Add(self.min_text)
        self.row.AddSpacer(5)
        self.row.Add(self.min_box)
        self.row.AddSpacer(10)
        self.row.Add(self.max_text)
        self.row.AddSpacer(5)
        self.row.Add(self.max_box)
        self.row.AddSpacer(10)

        self.input_ROI_sizer.AddSpacer(10)
        self.input_ROI_sizer.Add(self.row, 0, wx.ALIGN_CENTER_HORIZONTAL)

        self.save_region_button = wx.Button(self.input_ROI_label, -1, "Add ROI")

        self.save_region_button.Bind(wx.EVT_BUTTON, self.OnSaveRegion)

        self.input_ROI_sizer.AddSpacer(10)
        self.input_ROI_sizer.Add(self.save_region_button, 0, wx.ALIGN_CENTER_HORIZONTAL)

        self.main_ROI_input.AddSpacer(10)
        self.main_ROI_input.Add(self.input_ROI_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)

    def OnSaveRegion(self, event):

        # check validity of values inputed
        try:
            float(self.min_box.GetValue())
            float(self.max_box.GetValue())
        except:
            # The inputed values are not numbers - changing this now
            message = "The inputed values are not numbers, please change these to numbers and try again."
            dlg = wx.MessageDialog(self, message, "Warning", wx.OK)
            dlg.ShowModal()
            dlg.Destroy()
            return
        
        if(float(self.min_box.GetValue()) > float(self.max_box.GetValue())):
            # The inputed values are not numbers - changing this now
            message = "The minimum value is not less than the maximum value. Please ensure the minimum value is less than the maximum value and try again."
            dlg = wx.MessageDialog(self, message, "Warning", wx.OK)
            dlg.ShowModal()
            dlg.Destroy()
            return

        input_success = self.main_frame.add_user_input_region(float(self.min_box.GetValue()), float(self.max_box.GetValue()))

        # close the window
        if(input_success):
            self.Destroy()

class DelaysManualInput(wx.Frame):
    def __init__(self, title, parent):
        self.main_frame = parent
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        width = int(0.3 * self.monitorWidth)
        height = int(0.5 * self.monitorHeight)
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.panel_delays_input = wx.Panel(self, -1)
        self.main_delays_input = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_delays_input)


        self.make_manual_delays_input_sizer()
        self.Show()

    def make_manual_delays_input_sizer(self):
        self.input_delays_label = wx.StaticBox(
            self, -1, "Input delays in seconds (one delay per line)"
        )
        self.input_delays_sizer = wx.StaticBoxSizer(
            self.input_delays_label, wx.VERTICAL
        )

        try:
            file = open("delays.txt")
            delays = []
            for line in file.readlines():
                line = line.split("\n")[0]
                delays.append(line)

            label = ""
            for delay in delays:
                label = label + delay + "\n"
            file.close()
        except:
            label = ""

        self.delay_box = wx.TextCtrl(
            self.input_delays_label, -1, value=label, size=(250, 400), style=wx.TE_MULTILINE
        )
        self.input_delays_sizer.AddSpacer(3)
        self.input_delays_sizer.Add(self.delay_box)
        self.input_delays_sizer.AddSpacer(10)

        self.save_delays_button = wx.Button(self.input_delays_label, -1, "Save delays")

        self.save_delays_button.Bind(wx.EVT_BUTTON, self.OnSaveDelays)

        self.input_delays_sizer.Add(self.save_delays_button)

        self.main_delays_input.Add(self.input_delays_sizer)

    def OnSaveDelays(self, event):
        # Ensure that there are not any blank lines in the delays
        new_lines = ""
        for i, line in enumerate(self.delay_box.GetValue().split("\n")):
            if line != "":
                new_lines = new_lines + line.rstrip() + "\n"

        if new_lines[-1] == "\n":
            new_lines = new_lines[:-1]
        self.delay_box.SetValue(new_lines)

        # Check that all the delays are numbers

        for delay in self.delay_box.GetValue().split("\n"):
            try:
                float(delay)
            except:
                error = (
                    "Please ensure that all delays entrered are numbers: "
                    + delay
                    + " is not a number"
                )
                msg = wx.MessageDialog(self, error, "Error", wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                return

        # Ensure that all delays are positive
        for delay in self.delay_box.GetValue().split("\n"):
            if float(delay) < 0:
                error = (
                    "Please ensure that all delays entrered are positive: "
                    + delay
                    + " is negative"
                )
                msg = wx.MessageDialog(self, error, "Error", wx.OK | wx.ICON_ERROR)
                msg.ShowModal()
                msg.Destroy()
                return

        # Ensure that the number of delays is the same as the number of slices
        if len(self.delay_box.GetValue().split("\n")) != len(self.main_frame.y_data):
            error = (
                "Please ensure that the number of delays is the same as the number of slices. There are "
                + str(len(self.main_frame.y_data))
                + " slices, but "
                + str(len(self.delay_box.GetValue().split("\n")))
                + " delays were given."
            )
            msg = wx.MessageDialog(self, error, "Error", wx.OK | wx.ICON_ERROR)
            msg.ShowModal()
            msg.Destroy()
            return

        file = open("delays.txt", "w")
        file.write(self.delay_box.GetValue())
        file.close()
        self.main_frame.delays = []
        for line in self.delay_box.GetValue().split("\n"):
            if line != "":
                self.main_frame.delays.append(float(line))

        self.main_frame.delays = np.array(self.main_frame.delays)


class GradientsManualTextInput(wx.Frame):
    def __init__(self, title, parent):
        self.main_frame = parent
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        width = int(0.3 * self.monitorWidth)
        height = int(0.5 * self.monitorHeight)
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.panel_gradients_text_input = wx.Panel(self, -1)
        self.main_gradients_text_input = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_gradients_text_input)


        self.make_manual_gradients_text_input_sizer()
        self.Show()

    def make_manual_gradients_text_input_sizer(self):
        if self.main_frame.spectrometer == "Bruker":
            self.input_gradient_text_label = wx.StaticBox(
                self, -1, "Input gradient percentages (one per line)"
            )
            self.input_gradient_text_sizer = wx.StaticBoxSizer(
                self.input_gradient_text_label, wx.VERTICAL
            )
        else:
            self.input_gradient_text_label = wx.StaticBox(
                self, -1, "Input gradient DAC values (one per line)"
            )
            self.input_gradient_text_sizer = wx.StaticBoxSizer(
                self.input_gradient_text_label, wx.VERTICAL
            )

        try:
            file = open("gradients.txt")
            gradients_percent = []
            for line in file.readlines():
                line = line.split("\n")[0]
                gradients_percent.append(line)

            label = ""
            for gradient in gradients_percent:
                label = label + gradient + "\n"
            file.close()
        except:
            label = ""

        self.gradient_box = wx.TextCtrl(
            self.input_gradient_text_label, -1, value=label, size=(250, 400), style=wx.TE_MULTILINE
        )
        self.input_gradient_text_sizer.AddSpacer(3)
        self.input_gradient_text_sizer.Add(self.gradient_box)
        self.input_gradient_text_sizer.AddSpacer(10)

        self.save_gradients_button = wx.Button(self.input_gradient_text_label, -1, "Save gradients")

        self.save_gradients_button.Bind(wx.EVT_BUTTON, self.OnSaveGradients)

        self.input_gradient_text_sizer.Add(self.save_gradients_button)

        self.main_gradients_text_input.Add(self.input_gradient_text_sizer)

    def OnSaveGradients(self, event):
        # Check the max gradient and integral factors are valid numbers
        try:
            self.gradient_integral_factor = float(self.main_frame.integral_factor_box.GetValue())
        except:
            self.main_frame.error_message('integral factor', self.main_frame.integral_factor_box.GetValue())
            return
        try:
            self.max_gradient = float(self.main_frame.max_gradient_box.GetValue())
        except:
            self.main_frame.error_message('max gradient', self.main_frame.max_gradient_box.GetValue())
            return
        # Remove all extra empty lines at the end of the text box
        if self.gradient_box.GetValue().split("\n")[-1] == "":
            old_value = self.gradient_box.GetValue()
            new_value = ""
            for i, line in enumerate(old_value.split("\n")):
                if line != "":
                    if i == len(old_value.split("\n")) - 1:
                        new_value = new_value + line.rstrip()
                    else:
                        new_value = new_value + line.rstrip() + "\n"
            # If the last element in the string contains a newline character, remove it
            if new_value[-1] == "\n":
                new_value = new_value[:-1]
            self.gradient_box.SetValue(new_value)
            self.panel_gradients_text_input.Layout()

        if self.gradient_box.GetValue() == "":
            # Give an error message saying that the gradients must be inputted first
            msg = wx.MessageDialog(
                self,
                "Please input the gradients first before saving",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        # Ensure all text are floats/integers
        for line in self.gradient_box.GetValue().split("\n"):
            if line != "":
                try:
                    float(line)
                except:
                    # Give an error message saying that all lines must be floats
                    error = (
                        "All gradient values must be numbers: "
                        + line
                        + " is not a number"
                    )
                    msg = wx.MessageDialog(self, error, "Error", wx.OK | wx.ICON_ERROR)
                    msg.ShowModal()
                    msg.Destroy()
                    return

        # Ensure there are no negative values
        for line in self.gradient_box.GetValue().split("\n"):
            if line != "":
                if float(line) < 0:
                    # Give an error message saying that there must be no negative values
                    msg = wx.MessageDialog(
                        self,
                        "There must be no negative values",
                        "Error",
                        wx.OK | wx.ICON_ERROR,
                    )
                    msg.ShowModal()
                    msg.Destroy()
                    return
        # If Bruker, ensure that all gradient percentages are between 0 and 100
        if self.main_frame.spectrometer == "Bruker":
            for line in self.gradient_box.GetValue().split("\n"):
                if line != "":
                    if float(line) > 100:
                        # Give an error message saying that there must be no values greater than 100
                        error = (
                            "There must be no gradient percentages greater than 100: "
                            + line
                            + " is greater than 100"
                        )
                        msg = wx.MessageDialog(
                            self, error, "Error", wx.OK | wx.ICON_ERROR
                        )
                        msg.ShowModal()
                        msg.Destroy()
                        return
        # If Varian, ensure that all gradient DAC values are between 0 and 30000
        else:
            for line in self.gradient_box.GetValue().split("\n"):
                if line != "":
                    if float(line) > 30000:
                        # Give an error message saying that there must be no values greater than 30000
                        error = (
                            "There must be no DAC gradient values greater than 30000: "
                            + line
                            + " is greater than 30000"
                        )
                        msg = wx.MessageDialog(
                            self, error, "Error", wx.OK | wx.ICON_ERROR
                        )
                        msg.ShowModal()
                        msg.Destroy()
                        return
            
            try:
                self.DAC_conversion = float(self.main_frame.dac_conversion_box.GetValue())
            except:
                self.main_frame.error_message('DAC conversion', self.main_frame.dac_conversion_box.GetValue())
                return
        # Ensure that the number of gradients entered is the same as the number of slices
        if len(self.gradient_box.GetValue().split("\n")) != len(self.main_frame.y_data):
            # Give an error message saying that the number of gradients must be the same as the number of slices
            error = (
                "The number of gradients must be the same as the number of slices in the data: "
                + str(len(self.gradient_box.GetValue().split("\n")))
                + " gradients entered, "
                + str(len(self.main_frame.y_data))
                + " slices in the data"
            )
            msg = wx.MessageDialog(self, error, "Error", wx.OK | wx.ICON_ERROR)
            msg.ShowModal()
            msg.Destroy()
            return
        file = open("gradients.txt", "w")
        file.write(self.gradient_box.GetValue())
        file.close()
        
        if self.main_frame.spectrometer == "Bruker":
            self.main_frame.gradients_percent = []
            self.main_frame.gradients = []
            for line in self.gradient_box.GetValue().split("\n"):
                if line != "":
                    self.main_frame.gradients_percent.append(float(line))

            self.main_frame.gradients_percent = np.array(
                self.main_frame.gradients_percent
            )
            self.main_frame.gradients = (
                (self.main_frame.gradients_percent / 100)
                * self.max_gradient
                * self.gradient_integral_factor
            )

        else:
            self.gradients_DAC = []
            for line in self.gradient_box.GetValue().split("\n"):
                if line != "":
                    self.gradients_DAC.append(float(line))
            self.gradients_DAC = np.array(self.gradients_DAC)

            self.gradients = (
                np.array(self.gradients_DAC) * self.DAC_conversion
            )
            self.gradients_percent = (
                np.array(self.gradients_DAC)
                * self.DAC_conversion
                / self.max_gradient
                * 100
            )
            self.main_frame.gradients_percent = self.gradients_percent
            self.main_frame.gradients = self.gradients * self.gradient_integral_factor


class DiffusionGradientManualInput(wx.Frame):
    def __init__(self, title, parent, spectrometer="Bruker"):
        self.main_frame = parent
        self.spectrometer = spectrometer
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        width = int(0.3 * self.monitorWidth)
        height = int(0.35 * self.monitorHeight)
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.panel_gradient_input = wx.Panel(self, -1)


        # Define initial default values
        self.number_of_gradients = 5
        self.min_gradient_percent = 10.0
        self.max_gradient_percent = 90.0
        self.min_gradient_DAC = 10000
        self.max_gradient_DAC = 30000
        self.gradient_distribution = "squared"

        self.main_sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.make_manual_gradient_input_sizer()
        self.make_manual_gradients_text_input_sizer()

        self.main_sizer.AddSpacer(5)
        self.main_sizer.Add(self.gradient_input_sizer)
        self.main_sizer.AddSpacer(20)
        self.main_sizer.Add(self.input_gradient_text_sizer)
        self.main_sizer.AddSpacer(5)

        self.SetSizer(self.main_sizer)
        self.Show()

    def make_manual_gradient_input_sizer(self):
        self.gradient_input_sizer_label = wx.StaticBox(self, -1, "Calculate Gradients")
        self.gradient_input_sizer = wx.StaticBoxSizer(
            self.gradient_input_sizer_label, wx.VERTICAL
        )
        self.gradient_input_sizer.AddSpacer(5)
        # TextCtrl for number of gradient values used
        self.number_of_gradients_label = wx.StaticText(self.gradient_input_sizer_label, -1, "Number of Gradients:")
        self.gradient_input_sizer.Add(self.number_of_gradients_label)
        self.gradient_input_sizer.AddSpacer(5)
        self.number_of_gradients = int(len(self.main_frame.y_data))
        self.number_of_gradients_box = wx.TextCtrl(
            self.gradient_input_sizer_label, -1, str(self.number_of_gradients), size=(50, -1)
        )
        self.gradient_input_sizer.Add(self.number_of_gradients_box)
        self.gradient_input_sizer.AddSpacer(5)

        if self.spectrometer == "Bruker":
            # TextCtrl for minimum gradient percentage
            self.min_gradient_label = wx.StaticText(self.gradient_input_sizer_label, -1, "Min Gradient (%):")
            self.gradient_input_sizer.Add(self.min_gradient_label)
            self.gradient_input_sizer.AddSpacer(5)
            self.min_gradient_box = wx.TextCtrl(
                self.gradient_input_sizer_label, -1, str(self.min_gradient_percent), size=(50, -1)
            )
            self.gradient_input_sizer.Add(self.min_gradient_box)
            self.gradient_input_sizer.AddSpacer(5)

            # TextCtrl for maximum gradient percentage
            self.max_gradient_label = wx.StaticText(self.gradient_input_sizer_label, -1, "Max Gradient (%):")
            self.gradient_input_sizer.Add(self.max_gradient_label)
            self.gradient_input_sizer.AddSpacer(5)
            self.max_gradient_box = wx.TextCtrl(
                self.gradient_input_sizer_label, -1, str(self.max_gradient_percent), size=(50, -1)
            )
            self.gradient_input_sizer.Add(self.max_gradient_box)
            self.gradient_input_sizer.AddSpacer(5)
        else:
            # TextCtrl for minimum gradient DAC value
            self.min_gradient_label = wx.StaticText(self.gradient_input_sizer_label, -1, "Min Gradient (DAC):")
            self.gradient_input_sizer.Add(self.min_gradient_label)
            self.gradient_input_sizer.AddSpacer(5)
            self.min_gradient_box = wx.TextCtrl(
                self.gradient_input_sizer_label, -1, str(self.min_gradient_DAC), size=(50, -1)
            )
            self.gradient_input_sizer.Add(self.min_gradient_box)
            self.gradient_input_sizer.AddSpacer(5)

            # TextCtrl for maximum gradient DAC value
            self.max_gradient_label = wx.StaticText(self.gradient_input_sizer_label, -1, "Max Gradient (DAC):")
            self.gradient_input_sizer.Add(self.max_gradient_label)
            self.gradient_input_sizer.AddSpacer(5)
            self.max_gradient_box = wx.TextCtrl(
                self.gradient_input_sizer_label, -1, str(self.max_gradient_DAC), size=(50, -1)
            )
            self.gradient_input_sizer.Add(self.max_gradient_box)
            self.gradient_input_sizer.AddSpacer(5)

        # RadioBox for gradient distribution
        self.gradient_distribution_label = wx.StaticText(self.gradient_input_sizer_label, -1, "Gradient Spacing:")
        self.gradient_input_sizer.Add(self.gradient_distribution_label)
        self.gradient_input_sizer.AddSpacer(5)
        self.gradient_distribution_choices = ["linear", "squared", "exponential"]
        self.gradient_distribution_radiobox = wx.RadioBox(
            self.gradient_input_sizer_label, -1, choices=self.gradient_distribution_choices, style=wx.RA_VERTICAL
        )
        self.gradient_distribution_radiobox.SetSelection(1)
        self.gradient_input_sizer.Add(self.gradient_distribution_radiobox)
        self.gradient_input_sizer.AddSpacer(5)

        # Button to confirm input
        self.confirm_button = wx.Button(self.gradient_input_sizer_label, -1, "Calculate")
        if self.spectrometer == "Bruker":
            self.confirm_button.Bind(wx.EVT_BUTTON, self.OnCalculateGradientsBruker)
        else:
            self.confirm_button.Bind(wx.EVT_BUTTON, self.OnCalculateGradientsVarian)
        self.gradient_input_sizer.Add(self.confirm_button)
        self.gradient_input_sizer.AddSpacer(5)

    def OnCalculateGradientsBruker(self, event):
        # Check number of gradients is an integer
        if self.number_of_gradients_box.GetValue().isdigit() == False:
            msg = wx.MessageDialog(
                self,
                "Number of gradients must be an integer",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return

        if (
            int(self.number_of_gradients_box.GetValue())
            - len(self.main_frame.deleted_slices)
        ) != len(self.main_frame.y_data):
            error = (
                "Number of gradients must be equal to the number of slices: "
                + str(len(self.main_frame.y_data))
                + " slices in the data, "
                + str(int(self.number_of_gradients_box.GetValue()))
                + " gradients entered"
            )
            msg = wx.MessageDialog(self, error, "Error", wx.OK | wx.ICON_ERROR)
            msg.ShowModal()
            msg.Destroy()
            return

        # Check min and max gradient percentages are between 0 and 100
        if (
            float(self.min_gradient_box.GetValue()) < 0
            or float(self.min_gradient_box.GetValue()) > 100
        ):
            msg = wx.MessageDialog(
                self,
                "Min gradient percentage must be between 0 and 100",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        if (
            float(self.max_gradient_box.GetValue()) < 0
            or float(self.max_gradient_box.GetValue()) > 100
        ):
            msg = wx.MessageDialog(
                self,
                "Max gradient percentage must be between 0 and 100",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        # Check min gradient percentage is less than max gradient percentage
        if float(self.min_gradient_box.GetValue()) > float(
            self.max_gradient_box.GetValue()
        ):
            msg = wx.MessageDialog(
                self,
                "Min gradient percentage must be less than max gradient percentage",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        self.number_of_gradients = int(self.number_of_gradients_box.GetValue())
        self.min_gradient_percent = float(self.min_gradient_box.GetValue())
        self.max_gradient_percent = float(self.max_gradient_box.GetValue())
        self.gradient_distribution = self.gradient_distribution_choices[
            self.gradient_distribution_radiobox.GetSelection()
        ]
        self.CalculateGradientsBruker()
        self.gradients = (
            self.gradients_percent
            * self.main_frame.max_gradient
            * self.main_frame.integral_factor
            / 100
        )
        self.main_frame.gradients_percent = self.gradients_percent
        self.main_frame.gradients = self.gradients
        string_of_gradients = ""
        for i, gradient in enumerate(self.gradients_percent):
            if i == len(self.gradients_percent) - 1:
                string_of_gradients = string_of_gradients + str(gradient)
            else:
                string_of_gradients = string_of_gradients + str(gradient) + "\n"

        # Save gradients to text file
        file = open("gradients.txt", "w")
        file.write(string_of_gradients)
        file.close()

        # Write the gradients to the text box
        self.gradient_box.SetValue(string_of_gradients)

        if len(self.main_frame.y_data) != len(self.main_frame.gradients):
            for i, deleted_slice in enumerate(self.main_frame.deleted_slices):
                self.main_frame.gradients = np.delete(
                    self.main_frame.gradients, deleted_slice, axis=0
                )
                self.main_frame.gradients_percent = np.delete(
                    self.main_frame.gradients_percent, deleted_slice, axis=0
                )

    def CalculateGradientsBruker(self):
        if self.gradient_distribution == "linear":
            self.gradients_percent = np.linspace(
                self.min_gradient_percent,
                self.max_gradient_percent,
                self.number_of_gradients,
            )
        elif self.gradient_distribution == "squared":
            self.gradients_percent = np.sqrt(
                np.linspace(
                    self.min_gradient_percent**2,
                    self.max_gradient_percent**2,
                    self.number_of_gradients,
                )
            )
        else:
            self.gradients_percent = np.log(
                np.linspace(
                    np.exp(self.min_gradient_percent),
                    np.exp(self.max_gradient_percent),
                    self.number_of_gradients,
                )
            )

    def OnCalculateGradientsVarian(self, event):
        # Check number of gradients is an integer
        if self.number_of_gradients_box.GetValue().isdigit() == False:
            msg = wx.MessageDialog(
                self,
                "Number of gradients must be an integer",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return

        if (
            int(self.number_of_gradients_box.GetValue())
            - len(self.main_frame.deleted_slices)
        ) != len(self.main_frame.y_data):
            error = (
                "Number of gradients must be equal to the number of slices: "
                + str(len(self.main_frame.y_data))
                + " slices in the data, "
                + str(int(self.number_of_gradients_box.GetValue()))
                + " gradients entered"
            )
            msg = wx.MessageDialog(self, error, "Error", wx.OK | wx.ICON_ERROR)
            msg.ShowModal()
            msg.Destroy()
            return

        # Check min gradient DAC is a number
        if self.min_gradient_box.GetValue().isdigit() == False:
            msg = wx.MessageDialog(
                self,
                "Min gradient DAC value must be a number",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        # Check min gradient DAC is a number above 0
        if float(self.min_gradient_box.GetValue()) <= 0:
            msg = wx.MessageDialog(
                self,
                "Min gradient DAC value must be above 0",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        # Check max gradient DAC is a number
        if self.max_gradient_box.GetValue().isdigit() == False:
            msg = wx.MessageDialog(
                self,
                "Max gradient DAC value must be a number",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        # Check max gradient DAC is a number above 0
        if float(self.max_gradient_box.GetValue()) <= 0:
            msg = wx.MessageDialog(
                self,
                "Max gradient DAC value must be above 0",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return

        # Check min gradient percentage is less than max gradient percentage
        if float(self.min_gradient_box.GetValue()) > float(
            self.max_gradient_box.GetValue()
        ):
            msg = wx.MessageDialog(
                self,
                "Min gradient DAC value must be less than max gradient DAC value",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        self.number_of_gradients = int(self.number_of_gradients_box.GetValue())
        self.min_gradient_DAC = float(self.min_gradient_box.GetValue())
        self.max_gradient_DAC = float(self.max_gradient_box.GetValue())
        self.gradient_distribution = self.gradient_distribution_choices[
            self.gradient_distribution_radiobox.GetSelection()
        ]
        self.CalculateGradientsVarian()
        self.gradients = np.array(self.gradients_DAC) * self.main_frame.DAC_conversion
        self.gradients_percent = (
            np.array(self.gradients_DAC)
            * self.main_frame.DAC_conversion
            / self.main_frame.max_gradient
            * 100
        )
        self.main_frame.gradients_percent = self.gradients_percent
        self.main_frame.gradients = self.gradients
        if len(self.main_frame.y_data) != len(self.main_frame.gradients):
            for i, deleted_slice in enumerate(self.main_frame.deleted_slices):
                self.main_frame.gradients = np.delete(
                    self.main_frame.gradients, deleted_slice, axis=0
                )
                self.main_frame.gradients_percent = np.delete(
                    self.main_frame.gradients_percent, deleted_slice, axis=0
                )

        string_of_gradients = ""
        for i, gradient in enumerate(self.gradients_DAC):
            if i == len(self.gradients_DAC) - 1:
                string_of_gradients = string_of_gradients + str(gradient)
            else:
                string_of_gradients = string_of_gradients + str(gradient) + "\n"

        # Save gradients to text file
        file = open("gradients.txt", "w")
        file.write(string_of_gradients)
        file.close()

        # Write the gradients to the text box
        self.gradient_box.SetValue(string_of_gradients)

    def CalculateGradientsVarian(self):
        if self.gradient_distribution == "linear":
            self.gradients_DAC = np.linspace(
                self.min_gradient_DAC, self.max_gradient_DAC, self.number_of_gradients
            )
        elif self.gradient_distribution == "squared":
            self.gradients_DAC = np.sqrt(
                np.linspace(
                    self.min_gradient_DAC**2,
                    self.max_gradient_DAC**2,
                    self.number_of_gradients,
                )
            )
        else:
            self.gradients_DAC = np.log(
                np.linspace(
                    np.exp(self.min_gradient_DAC),
                    np.exp(self.max_gradient_DAC),
                    self.number_of_gradients,
                )
            )

    def make_manual_gradients_text_input_sizer(self):
        if self.main_frame.spectrometer == "Bruker":
            self.input_gradient_text_label = wx.StaticBox(
                self, -1, "Input gradient percentages (one per line)"
            )
            self.input_gradient_text_sizer = wx.StaticBoxSizer(
                self.input_gradient_text_label, wx.VERTICAL
            )
        else:
            self.input_gradient_text_label = wx.StaticBox(
                self, -1, "Input gradient DAC values (one per line)"
            )
            self.input_gradient_text_sizer = wx.StaticBoxSizer(
                self.input_gradient_text_label, wx.VERTICAL
            )

        try:
            file = open("gradients.txt")
            gradients_percent = []
            for line in file.readlines():
                line = line.split("\n")[0]
                gradients_percent.append(line)

            label = ""
            for gradient in gradients_percent:
                label = label + gradient + "\n"
            file.close()
        except:
            label = ""

        self.gradient_box = wx.TextCtrl(
            self.input_gradient_text_label, -1, value=label, size=(250, 235), style=wx.TE_MULTILINE
        )
        self.input_gradient_text_sizer.AddSpacer(3)
        self.input_gradient_text_sizer.Add(self.gradient_box)
        self.input_gradient_text_sizer.AddSpacer(10)

        self.save_gradients_button = wx.Button(self.input_gradient_text_label, -1, "Save gradients")

        self.save_gradients_button.Bind(wx.EVT_BUTTON, self.OnSaveGradients)

        self.input_gradient_text_sizer.Add(self.save_gradients_button)

    def OnSaveGradients(self, event):
        # Remove all extra empty lines at the end of the text box
        if self.gradient_box.GetValue().split("\n")[-1] == "":
            old_value = self.gradient_box.GetValue()
            new_value = ""
            for i, line in enumerate(old_value.split("\n")):
                if line != "":
                    if i == len(old_value.split("\n")) - 1:
                        new_value = new_value + line.rstrip()
                    else:
                        new_value = new_value + line.rstrip() + "\n"
            # If the last element in the string contains a newline character, remove it
            if new_value[-1] == "\n":
                new_value = new_value[:-1]
            self.gradient_box.SetValue(new_value)
            self.panel_gradient_input.Layout()

        if self.gradient_box.GetValue() == "":
            # Give an error message saying that the gradients must be inputted first
            msg = wx.MessageDialog(
                self,
                "Please input the gradients first before saving",
                "Error",
                wx.OK | wx.ICON_ERROR,
            )
            msg.ShowModal()
            msg.Destroy()
            return
        # Ensure all text are floats/integers
        for line in self.gradient_box.GetValue().split("\n"):
            if line != "":
                try:
                    float(line)
                except:
                    # Give an error message saying that all lines must be floats
                    error = (
                        "All gradient values must be numbers: "
                        + line
                        + " is not a number"
                    )
                    msg = wx.MessageDialog(self, error, "Error", wx.OK | wx.ICON_ERROR)
                    msg.ShowModal()
                    msg.Destroy()
                    return

        # Ensure there are no negative values
        for line in self.gradient_box.GetValue().split("\n"):
            if line != "":
                if float(line) < 0:
                    # Give an error message saying that there must be no negative values
                    msg = wx.MessageDialog(
                        self,
                        "There must be no negative values",
                        "Error",
                        wx.OK | wx.ICON_ERROR,
                    )
                    msg.ShowModal()
                    msg.Destroy()
                    return
        # If Bruker, ensure that all gradient percentages are between 0 and 100
        if self.main_frame.spectrometer == "Bruker":
            for line in self.gradient_box.GetValue().split("\n"):
                if line != "":
                    if float(line) > 100:
                        # Give an error message saying that there must be no values greater than 100
                        error = (
                            "There must be no gradient percentages greater than 100: "
                            + line
                            + " is greater than 100"
                        )
                        msg = wx.MessageDialog(
                            self, error, "Error", wx.OK | wx.ICON_ERROR
                        )
                        msg.ShowModal()
                        msg.Destroy()
                        return
        # If Varian, ensure that all gradient DAC values are between 0 and 30000
        else:
            for line in self.gradient_box.GetValue().split("\n"):
                if line != "":
                    if float(line) > 30000:
                        # Give an error message saying that there must be no values greater than 30000
                        error = (
                            "There must be no DAC gradient values greater than 30000: "
                            + line
                            + " is greater than 30000"
                        )
                        msg = wx.MessageDialog(
                            self, error, "Error", wx.OK | wx.ICON_ERROR
                        )
                        msg.ShowModal()
                        msg.Destroy()
                        return
        # Ensure that the number of gradients entered is the same as the number of slices
        if len(self.gradient_box.GetValue().split("\n")) != len(self.main_frame.y_data):
            # Give an error message saying that the number of gradients must be the same as the number of slices
            error = (
                "The number of gradients must be the same as the number of slices in the data: "
                + str(len(self.gradient_box.GetValue().split("\n")))
                + " gradients entered, "
                + str(len(self.main_frame.y_data))
                + " slices in the data"
            )
            msg = wx.MessageDialog(self, error, "Error", wx.OK | wx.ICON_ERROR)
            msg.ShowModal()
            msg.Destroy()
            return
        file = open("gradients.txt", "w")
        file.write(self.gradient_box.GetValue())
        file.close()
        if self.main_frame.spectrometer == "Bruker":
            self.main_frame.gradients_percent = []
            self.main_frame.gradients = []
            for line in self.gradient_box.GetValue().split("\n"):
                if line != "":
                    self.main_frame.gradients_percent.append(float(line))

            self.main_frame.gradients_percent = np.array(
                self.main_frame.gradients_percent
            )
            self.main_frame.gradients = (
                (self.main_frame.gradients_percent / 100)
                * self.main_frame.max_gradient
                * self.main_frame.integral_factor
            )

        else:
            self.gradients_DAC = []
            for line in self.gradient_box.GetValue().split("\n"):
                if line != "":
                    self.gradients_DAC.append(float(line))
            self.gradients_DAC = np.array(self.gradients_DAC)

            self.gradients = (
                np.array(self.gradients_DAC) * self.main_frame.DAC_conversion
            )
            self.gradients_percent = (
                np.array(self.gradients_DAC)
                * self.main_frame.DAC_conversion
                / self.main_frame.max_gradient
                * 100
            )
            self.main_frame.gradients_percent = self.gradients_percent
            self.main_frame.gradients = self.gradients

