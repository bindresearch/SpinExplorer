import wx
import numpy as np
import nmrglue as ng
import matplotlib
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import (
    NavigationToolbar2WxAgg as NavigationToolbar,
)


class CESTOrder_Dialog(wx.Dialog):
    def __init__(self, title, parent):
        self.main_frame = parent
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        width = 300
        height = 150
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.panel_CESTOrder = wx.Panel(self, -1)
        self.main_CESTOrder = wx.BoxSizer(wx.HORIZONTAL)
        self.SetSizer(self.main_CESTOrder)

        self.make_CESTOrder_sizer()
        self.Show()

    def make_CESTOrder_sizer(self):
        # Make a sizer to hold the text box and button
        self.CESTOrder_sizer = wx.BoxSizer(wx.VERTICAL)
        self.CESTOrder_sizer.AddSpacer(5)

        self.CESTOrder_label = wx.StaticText(
            self, -1, "Interleaved order:", style=wx.ALIGN_LEFT
        )

        choices = ["OnResonance, OffResonance", "OffResonance, OnResonance"]
        self.CESTOrder_radiobox = wx.RadioBox(
            self, -1, choices=choices, style=wx.CB_READONLY | wx.RA_SPECIFY_ROWS
        )
        self.CESTOrder_sizer.AddSpacer(10)
        self.CESTOrder_sizer.Add(self.CESTOrder_radiobox, wx.ALIGN_LEFT)
        self.CESTOrder_sizer.AddSpacer(10)
        # Have a button to confirm the selection
        self.confirm_button = wx.Button(self, -1, "Confirm")
        self.confirm_button.Bind(wx.EVT_BUTTON, self.OnConfirm)
        self.CESTOrder_sizer.Add(self.confirm_button, wx.ALIGN_CENTER)
        self.CESTOrder_sizer.AddSpacer(10)

        self.main_CESTOrder.AddSpacer(10)
        self.main_CESTOrder.Add(self.CESTOrder_sizer, wx.ALIGN_CENTER)
        self.main_CESTOrder.AddSpacer(10)

    def OnConfirm(self, event):
        self.main_frame.CESTArrayOrder = self.CESTOrder_radiobox.GetSelection()
        self.main_frame.continue_deletion()
        self.Destroy()


class CESTFrame(wx.Frame):
    def __init__(self, title, parent=None, CESTArrayOrder=0):
        self.main_frame = parent
        # Get the monitor size and set the window size to 85% of the monitor size
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        self.display_index = wx.Display.GetFromWindow(parent)
        self.display_index_current = self.display_index
        self.width = int(1.0 * sizes[self.display_index][0])
        self.height = int(0.875 * sizes[self.display_index][1])
        wx.Frame.__init__(
            self, parent=parent, title=title, size=(self.width, self.height)
        )
        self.panel_CEST = wx.Panel(self, -1)
        self.main_CEST_sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.main_CEST_sizer)

        self.CESTArrayOrder = CESTArrayOrder

        self.fig_CEST = Figure()
        self.fig_CEST.tight_layout()
        self.canvas_CEST = FigCanvas(self, -1, self.fig_CEST)
        self.main_CEST_sizer.Add(self.canvas_CEST, 1, flag=wx.EXPAND | wx.ALL)
        self.toolbar_CEST = NavigationToolbar(self.canvas_CEST)
        self.main_CEST_sizer.Add(self.toolbar_CEST, 0, wx.EXPAND)


        self.make_CEST_sizer()

        self.titlecolor = "black"

        self.find_CEST_frequencies()
        self.organise_CEST_data()

        self.plot_CEST_data()
        self.Show()

        # Bind method to check/resize the window when the frame is moved
        self.Bind(wx.EVT_MOVE, self.OnMoveFrame)

        # Bind method to resize the window when the frame is resized
        self.Bind(wx.EVT_SIZE, self.OnSizeFrame)

    def find_CEST_frequencies(self):
        # Try to find procpar file in current directory, if can find it, work out frequencies based on tof_sel
        # If cannot find it will just have to the frequency indexes of 0 to n
        try:
            dic_v, data_v = ng.varian.read("./")

            # find the CEST offset values used
            offsets = dic_v["procpar"]["tof_sel"]["values"]
            self.offsets_Hz = []
            for i in range(len(offsets)):
                self.offsets_Hz.append(float(offsets[i]))

            # get ppm values for the offsets used
            self.tof = float(dic_v["procpar"]["tof"]["values"][0])

            # get the sfrq
            self.sfrq = float(dic_v["procpar"]["sfrq"]["values"][0])

            # open fid.com to find the carrier frequency
            with open("fid.com", "r") as file:
                lines = file.readlines()
                for line in lines:
                    if "-xCAR" in line:
                        self.carrier = float(line.split()[1])

            self.offsets_ppm = []
            for tof_sel in self.offsets_Hz:
                offset_ppm = (
                    tof_sel - self.tof
                ) / self.sfrq + self.carrier  # get actual ppm values
                self.offsets_ppm.append(offset_ppm)

        except:
            self.offsets_ppm = np.arange(0, len(self.main_frame.ppms_0), 1)

    def organise_CEST_data(self):
        # Get the CEST data from the main frame
        self.CEST_data = self.main_frame.nmrdata.data.T
        self.cest_on_data = []
        self.cest_off_data = []
        for i in range(len(self.CEST_data)):
            if i % 2 == 0:
                if self.CESTArrayOrder == 0:
                    self.cest_on_data.append(self.CEST_data[i])
                else:
                    self.cest_off_data.append(self.CEST_data[i])
            else:
                if self.CESTArrayOrder == 0:
                    self.cest_off_data.append(self.CEST_data[i])
                else:
                    self.cest_on_data.append(self.CEST_data[i])

        # Find the selected 1H chemical shift range in the main frame
        self.selected_shift = self.main_frame.line4.get_xdata()[0]
        self.selected_shift_index = np.argmin(
            np.abs(self.main_frame.ppms_1 - self.selected_shift)
        )

        self.non_normalized_cest_data = []
        self.normalized_cest_data = []
        for i in range(len(self.cest_on_data)):
            self.non_normalized_cest_data.append(
                self.cest_on_data[i][self.selected_shift_index]
            )
            self.normalized_cest_data.append(
                self.cest_on_data[i][self.selected_shift_index]
                / self.cest_off_data[i][self.selected_shift_index]
            )

    def OnMoveFrame(self, event):
        # Get the new default display if the frame is moved
        displays = (wx.Display(i) for i in range(wx.Display.GetCount()))
        sizes = [display.GetGeometry().GetSize() for display in displays]
        display_index = wx.Display.GetFromWindow(self)
        if display_index != self.display_index_current:
            self.display_index_current = display_index
            self.width = int(1.0 * sizes[display_index][0])
            self.height = int(0.875 * sizes[display_index][1])
            self.SetSize((self.width, self.height))
            self.canvas_CEST.SetSize(
                (
                    self.width * 0.0104,
                    (self.height - self.CEST_sizer.GetMinSize()[1] - 100) * 0.0104,
                )
            )
            self.fig_CEST.set_size_inches(
                self.width * 0.0104,
                (self.height - self.CEST_sizer.GetMinSize()[1] - 100) * 0.0104,
            )
            self.UpdateCESTFrame()
        event.Skip()

    def OnSizeFrame(self, event):
        # Get the new frame size
        self.width, self.height = self.GetSize()
        self.SetSize((self.width, self.height))
        self.canvas_CEST.SetSize(
            (
                self.width * 0.0104,
                (self.height - self.CEST_sizer.GetMinSize()[1] - 100) * 0.0104,
            )
        )
        self.fig_CEST.set_size_inches(
            self.width * 0.0104,
            (self.height - self.CEST_sizer.GetMinSize()[1] - 100) * 0.0104,
        )
        self.UpdateCESTFrame()
        event.Skip()

    def UpdateCESTFrame(self):
        self.canvas_CEST.draw()
        self.canvas_CEST.Refresh()
        self.canvas_CEST.Update()
        self.panel_CEST.Refresh()
        self.panel_CEST.Update()

    def make_CEST_sizer(self):
        self.CEST_sizer = wx.BoxSizer(wx.HORIZONTAL)

        # Create a button that opens a file for a user to input the delay times
        self.CEST_ppm_range_label = wx.StaticBox(
            self, -1, "3D plot chemical shift range"
        )
        self.CEST_ppm_sizer_total = wx.StaticBoxSizer(
            self.CEST_ppm_range_label, wx.VERTICAL
        )
        self.CEST_ppm_sizer = wx.BoxSizer(wx.HORIZONTAL)

        # ppm min value
        self.CEST_ppm_min_label = wx.StaticText(self, -1, "Min ppm:")
        self.min_val = min(self.main_frame.ppms_1)
        self.CEST_ppm_min_text = wx.TextCtrl(
            self, -1, str(self.min_val), style=wx.TE_PROCESS_ENTER
        )
        self.CEST_ppm_min_text.Bind(wx.EVT_TEXT_ENTER, self.OnCEST_ppm_change)
        self.CEST_ppm_sizer.Add(self.CEST_ppm_min_label, wx.ALIGN_CENTER)
        self.CEST_ppm_sizer.AddSpacer(5)
        self.CEST_ppm_sizer.Add(self.CEST_ppm_min_text)

        # ppm max value
        self.CEST_ppm_max_label = wx.StaticText(self, -1, "Max ppm:")
        self.max_val = max(self.main_frame.ppms_1)
        self.CEST_ppm_max_text = wx.TextCtrl(
            self, -1, str(self.max_val), style=wx.TE_PROCESS_ENTER
        )
        self.CEST_ppm_max_text.Bind(wx.EVT_TEXT_ENTER, self.OnCEST_ppm_change)
        self.CEST_ppm_sizer.AddSpacer(10)
        self.CEST_ppm_sizer.Add(self.CEST_ppm_max_label)
        self.CEST_ppm_sizer.AddSpacer(5)
        self.CEST_ppm_sizer.Add(self.CEST_ppm_max_text)

        self.CEST_ppm_sizer_total.Add(self.CEST_ppm_sizer)
        self.CEST_sizer.Add(self.CEST_ppm_sizer_total, wx.ALIGN_CENTER_HORIZONTAL)

        # Make a button which will save the normalised CEST data as a spectrum
        self.save_CEST_button = wx.Button(self, -1, "Save normalised CEST data")
        self.save_CEST_button.Bind(wx.EVT_BUTTON, self.OnSaveCESTData)
        self.CEST_sizer.AddSpacer(10)
        self.CEST_sizer.Add(self.save_CEST_button, wx.ALIGN_CENTER_HORIZONTAL)

        self.main_CEST_sizer.AddSpacer(10)
        self.main_CEST_sizer.Add(self.CEST_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.main_CEST_sizer.AddSpacer(10)

    def OnSaveCESTData(self, event):
        # The ppms are self.offsets_ppm
        # The data is self.normalized_cest_data
        # Save the data as an nmrPipe spectrum using nmrglue

        try:
            # Ask the user what the file name should be
            file_name = ""
            message = "Input the file name to save the CEST data as"
            dlg = wx.TextEntryDialog(None, message, "Save CEST data", "CEST_data.ft")
            if dlg.ShowModal() == wx.ID_OK:
                file_name = dlg.GetValue()
                dlg.Destroy()

            if file_name == "":
                # Give an error message to say that the user must input a file name
                message = "Error: No file name given"
                dlg = wx.MessageDialog(None, message, "Error", wx.OK | wx.ICON_ERROR)
                dlg.ShowModal()
                dlg.Destroy()
                return

            # Create the nmrPipe file
            data = np.flip(np.array(self.normalized_cest_data) * 100)
            data = data.astype(np.float32)

            obs = self.sfrq
            sw = max(self.offsets_Hz) - min(self.offsets_Hz)
            car = self.carrier
            size = len(data)
            label = "ppm"
            orig = car * obs - sw * (size / 2 - 1) / size
            center = self.main_frame.nmrdata.dic["FDF2CENTER"]
            udic = {
                "ndim": 1,
                0: {
                    "size": size,
                    "complex": False,
                    "encoding": "int",
                    "sw": sw,
                    "obs": obs,
                    "car": car,
                    "label": label,
                    "time": False,
                    "freq": True,
                },
            }

            dic = ng.pipe.create_dic(udic)

            dic["FDF2LABEL"] = label
            dic["FDF2OBS"] = obs
            dic["FDF2SW"] = sw
            dic["FDF2CAR"] = car
            dic["FDF2SIZE"] = size
            dic["FDF2ORIG"] = orig
            dic["FDF2CENTER"] = center
            ng.pipe.write(file_name, dic, data, overwrite=True)

            message = "CEST data saved as 'CEST_data.ft'"
            dlg = wx.MessageDialog(
                None, message, "Save successful", wx.OK | wx.ICON_INFORMATION
            )
            dlg.ShowModal()
            dlg.Destroy()

        except:
            message = "Error saving CEST data"
            dlg = wx.MessageDialog(None, message, "Error", wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()

    def plot_CEST_data(self):

        self.ax = self.fig_CEST.add_subplot(121, projection="3d")

        for i, y in enumerate(self.offsets_ppm):
            self.ax.plot(
                self.main_frame.ppms_1,
                np.full_like(self.main_frame.ppms_1, y),
                self.cest_on_data[i],
                color="tab:grey",
                linewidth=1.5,
                alpha=0.5,
            )
            self.ax.scatter(
                self.main_frame.ppms_1[self.selected_shift_index],
                y,
                self.non_normalized_cest_data[i],
                color="red",
                s=10,
            )

        self.ax.set_title("CEST data")
        self.ax.set_xlabel(self.main_frame.nmrdata.axislabels[1])
        self.ax.set_ylabel(self.main_frame.nmrdata.axislabels[0])

        self.ax.set_xlim([max(self.main_frame.ppms_1), min(self.main_frame.ppms_1)])
        self.ax.set_ylim([max(self.offsets_ppm), min(self.offsets_ppm)])
        self.ax.view_init(elev=10.0, azim=-45)

        self.ax2 = self.fig_CEST.add_subplot(122)
        self.ax2.plot(
            self.offsets_ppm, self.normalized_cest_data, color="tab:red", linewidth=1.5
        )
        self.ax2.set_title("Normalized CEST data")
        self.ax2.set_xlabel("CEST offset (ppm)")
        self.ax2.set_ylabel("Normalized CEST data")
        self.ax2.set_xlim([max(self.offsets_ppm), min(self.offsets_ppm)])
        self.ax2.set_ylim([-0.1, 1.1])

        self.UpdateCESTFrame()

    def OnCEST_ppm_change(self, event):
        min_val = float(self.CEST_ppm_min_text.GetValue())
        max_val = float(self.CEST_ppm_max_text.GetValue())

        # Find the indexes of the min and max values in ppms_1
        min_index = np.argmin(np.abs(self.main_frame.ppms_1 - min_val))
        max_index = np.argmin(np.abs(self.main_frame.ppms_1 - max_val))

        if min_index > max_index:
            min_index, max_index = max_index, min_index

        # Get the data between the min and max values
        self.cest_on_data_new = []
        self.cest_off_data_new = []
        for i in range(len(self.cest_on_data)):
            self.cest_on_data_new.append(
                self.cest_on_data[i].tolist()[min_index:max_index]
            )
            self.cest_off_data_new.append(self.cest_off_data[i][min_index:max_index])

        self.ax.clear()
        for i, y in enumerate(self.offsets_ppm):
            self.ax.plot(
                self.main_frame.ppms_1[min_index:max_index],
                np.full_like(self.main_frame.ppms_1[min_index:max_index], y),
                self.cest_on_data_new[i],
                color="tab:grey",
                linewidth=1.5,
                alpha=0.5,
            )
            if min_index <= self.selected_shift_index <= max_index:
                self.ax.scatter(
                    self.main_frame.ppms_1[self.selected_shift_index],
                    y,
                    self.non_normalized_cest_data[i],
                    color="red",
                    s=10,
                )

        self.ax.set_title("CEST data")
        self.ax.set_xlabel(self.main_frame.nmrdata.axislabels[1])
        self.ax.set_ylabel(self.main_frame.nmrdata.axislabels[0])

        self.ax.set_xlim(
            [
                max(self.main_frame.ppms_1[min_index:max_index]),
                min(self.main_frame.ppms_1[min_index:max_index]),
            ]
        )
        self.ax.set_ylim([max(self.offsets_ppm), min(self.offsets_ppm)])

        self.UpdateCESTFrame()
