import wx # type: ignore
import os

class uSTA_Dialog(wx.Dialog):
    def __init__(self, title, parent):
        self.main_frame = parent
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        width = 450
        height = 150
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.panel_uSTAparams = wx.Panel(self, -1)
        self.main_uSTAparams = wx.BoxSizer(wx.HORIZONTAL)
        self.SetSizer(self.main_uSTAparams)

        self.make_uSTAparams_sizer()
        self.Show()

    def make_uSTAparams_sizer(self):
        # Make a sizer to hold the text box and button
        self.uSTAparams_sizer = wx.BoxSizer(wx.VERTICAL)
        self.uSTAparams_sizer.AddSpacer(30)

        self.uSTAparams_sizer2 = wx.BoxSizer(wx.HORIZONTAL)
        self.uSTAparams_sizer2.AddSpacer(5)

        self.mixing_time_label = wx.StaticText(self, -1, "Mixing time (s):")
        self.uSTAparams_sizer2.Add(self.mixing_time_label)

        self.uSTA_mixing_time = wx.TextCtrl(
            self, -1, value=self.main_frame.mixing_time, size=(100, 20)
        )
        self.uSTAparams_sizer2.AddSpacer(10)
        self.uSTAparams_sizer2.Add(self.uSTA_mixing_time)
        self.uSTAparams_sizer2.AddSpacer(10)

        self.power_level_label = wx.StaticText(self, -1, "Power level:")

        self.uSTAparams_sizer2.Add(self.power_level_label)
        self.uSTA_power_level = wx.TextCtrl(
            self, -1, value=self.main_frame.power_level, size=(100, 20)
        )
        self.uSTAparams_sizer2.AddSpacer(10)
        self.uSTAparams_sizer2.Add(self.uSTA_power_level)
        self.uSTAparams_sizer2.AddSpacer(10)

        self.uSTAparams_sizer.Add(self.uSTAparams_sizer2)

        self.uSTAparams_sizer.AddSpacer(30)

        # Have a button to confirm the selection
        self.confirm_button = wx.Button(self, -1, "Confirm")
        self.confirm_button.Bind(wx.EVT_BUTTON, self.OnConfirm)
        self.uSTAparams_sizer.Add(self.confirm_button, 0, wx.ALIGN_CENTER | wx.ALL, 10)
        self.uSTAparams_sizer.AddSpacer(10)

        self.main_uSTAparams.AddSpacer(10)
        self.main_uSTAparams.Add(self.uSTAparams_sizer)
        self.main_uSTAparams.AddSpacer(10)

    def OnConfirm(self, event):

        # Getting the ppm values
        ppm_values = self.main_frame.ppms_0
        if len(ppm_values) == 2:
            ppm_values = self.main_frame.ppms_1

        # Getting the current directory name
        current_dir_name = os.path.basename(os.getcwd())
        # Getting the intensities of on and on resonance spectra
        usta_data = self.main_frame.nmrdata.data
        if len(usta_data) > 2:
            usta_data = self.main_frame.nmrdata.data.T

        on_data = usta_data[0]
        off_data = usta_data[1]

        data_file_name = current_dir_name + ".data"
        data_file = open(data_file_name, "w")
        for i in range(len(ppm_values)):
            data_file.write(
                current_dir_name
                + "\t"
                + self.uSTA_mixing_time.GetValue()
                + "\t"
                + self.uSTA_power_level.GetValue()
                + "\t"
                + str(ppm_values[i])
                + "\t"
                + str(on_data[i])
                + "\t"
                + str(off_data[i])
                + "\n"
            )
        data_file.close()

        self.Destroy()