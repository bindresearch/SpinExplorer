import wx # type: ignore


# This class creates sliders which can contain floating point values
# Source: https://stackoverflow.com/questions/4709087/wxslider-with-floating-point-values
class FloatSlider(wx.Slider):

    def __init__(
        self,
        parent,
        id,
        value,
        minval,
        maxval,
        res,
        size=wx.DefaultSize,
        style=wx.SL_HORIZONTAL,
        name="floatslider",
    ):
        self._value = value
        self._min = minval
        self._max = maxval
        self._res = res
        ival, imin, imax = [round(v / res) for v in (value, minval, maxval)]
        self._islider = super(FloatSlider, self)
        self._islider.__init__(
            parent, id, ival, imin, imax, size=size, style=style, name=name
        )
        self.Bind(wx.EVT_SCROLL, self._OnScroll)

    def _OnScroll(self, event):
        ival = self._islider.GetValue()
        imin = self._islider.GetMin()
        imax = self._islider.GetMax()
        if ival == imin:
            self._value = self._min
        elif ival == imax:
            self._value = self._max
        else:
            self._value = ival * self._res
        event.Skip()

    def GetValue(self):
        return self._value

    def GetMin(self):
        return self._min

    def GetMax(self):
        return self._max

    def GetRes(self):
        return self._res

    def SetValue(self, value):
        self._islider.SetValue(round(value / self._res))
        self._value = value

    def SetMin(self, minval):
        self._islider.SetMin(round(minval / self._res))
        self._min = minval

    def SetMax(self, maxval):
        self._islider.SetMax(round(maxval / self._res))
        self._max = maxval

    def SetRes(self, res):
        self._islider.SetRange(round(self._min / res), round(self._max / res))
        self._islider.SetValue(round(self._value / res))
        self._res = res

    def SetRange(self, minval, maxval):
        self._islider.SetRange(round(minval / self._res), round(maxval / self._res))
        self._min = minval
        self._max = maxval


class PhasingSliderRange(wx.Frame):
    def __init__(self, title, parent):
        self.main_frame = parent
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        width = 400
        height = 200
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.panel_slider_range = wx.Panel(self, -1)



        self.main_sizer = wx.BoxSizer(wx.VERTICAL)

        self.make_slider_range_sizer()


        self.SetSizer(self.main_sizer)
        self.Show()

    def make_slider_range_sizer(self):
        self.P0_label = wx.StaticBox(self, -1, "P0 range:")
        self.P0_sizer = wx.StaticBoxSizer(self.P0_label, wx.VERTICAL)
        
        self.coarse_label_p0 = wx.StaticText(self.P0_label, -1, "Coarse (+/-):")
        coarse_range_p0 = self.main_frame.P0_slider.GetMax()
        self.coarse_box = wx.TextCtrl(self.P0_label, -1, str(coarse_range_p0))

        self.fine_label_p0 = wx.StaticText(self.P0_label, -1, "Fine (+/-):")
        fine_range_p0 = self.main_frame.P0_slider_fine.GetMax()
        self.fine_box = wx.TextCtrl(self.P0_label, -1, str(fine_range_p0))

        self.p0_row1 = wx.BoxSizer(wx.HORIZONTAL)
        self.p0_row1.Add(self.coarse_label_p0)
        self.p0_row1.AddSpacer(5)
        self.p0_row1.Add(self.coarse_box)

        self.p0_row2 = wx.BoxSizer(wx.HORIZONTAL)
        self.p0_row2.Add(self.fine_label_p0)
        self.p0_row2.AddSpacer(5)
        self.p0_row2.Add(self.fine_box)

        self.P0_sizer.Add(self.p0_row1, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.P0_sizer.AddSpacer(10)
        self.P0_sizer.Add(self.p0_row2, 0, wx.ALIGN_CENTER_HORIZONTAL)

        self.P1_label = wx.StaticBox(self, -1, "P1 range:")
        self.P1_sizer = wx.StaticBoxSizer(self.P1_label, wx.VERTICAL)
        
        self.coarse_label_p1 = wx.StaticText(self.P1_label, -1, "Coarse (+/-):")
        coarse_range_p1 = self.main_frame.P1_slider.GetMax()
        self.coarse_box_p1 = wx.TextCtrl(self.P1_label, -1, str(coarse_range_p1))

        self.fine_label_p1 = wx.StaticText(self.P1_label, -1, "Fine (+/-):")
        fine_range_p1 = self.main_frame.P1_slider_fine.GetMax()
        self.fine_box_p1 = wx.TextCtrl(self.P1_label, -1, str(fine_range_p1))

        self.p1_row1 = wx.BoxSizer(wx.HORIZONTAL)
        self.p1_row1.Add(self.coarse_label_p1)
        self.p1_row1.AddSpacer(5)
        self.p1_row1.Add(self.coarse_box_p1)

        self.p1_row2 = wx.BoxSizer(wx.HORIZONTAL)
        self.p1_row2.Add(self.fine_label_p1)
        self.p1_row2.AddSpacer(5)
        self.p1_row2.Add(self.fine_box_p1)

        self.P1_sizer.Add(self.p1_row1, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.P1_sizer.AddSpacer(10)
        self.P1_sizer.Add(self.p1_row2, 0, wx.ALIGN_CENTER_HORIZONTAL)

        self.slider_range_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.slider_range_sizer.AddSpacer(10)
        self.slider_range_sizer.Add(self.P0_sizer, 0, wx.ALIGN_CENTER_VERTICAL)
        self.slider_range_sizer.AddSpacer(10)
        self.slider_range_sizer.Add(self.P1_sizer, 0, wx.ALIGN_CENTER_VERTICAL)

        self.main_sizer.AddSpacer(5)
        self.main_sizer.Add(self.slider_range_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)

        self.update_button = wx.Button(self, -1, "Update Sliders")
        self.update_button.Bind(wx.EVT_BUTTON, self.on_update_button)

        self.main_sizer.AddSpacer(10)
        self.main_sizer.Add(self.update_button, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.main_sizer.AddSpacer(10)

    def on_update_button(self, event):
        if(self.check_values()==True):
            # Update the slider ranges with the new values
            p0_coarse = float(self.coarse_box.GetValue())
            self.main_frame.P0_slider.SetRange(-1*p0_coarse, p0_coarse)
            p1_coarse = float(self.coarse_box_p1.GetValue())
            self.main_frame.P1_slider.SetRange(-1*p1_coarse, p1_coarse)
            p0_fine = float(self.fine_box.GetValue())
            self.main_frame.P0_slider_fine.SetRange(-1*p0_fine, p0_fine)
            p1_fine = float(self.fine_box_p1.GetValue())
            self.main_frame.P1_slider_fine.SetRange(-1*p1_fine, p1_fine)

            dlg = wx.MessageDialog(
                self,
                "Slider ranges have been updated successfully",
                "Warning",
                wx.OK
            )
            self.Raise()
            self.SetFocus()
            dlg.ShowModal()
            dlg.Destroy()


    def check_values(self):
        # check that all the values can be converted to floats
        try:
            p0_coarse = float(self.coarse_box.GetValue())
        except:
            self.error_message("p0 coarse", self.coarse_box.GetValue())
            return False

        try:
            p1_coarse = float(self.coarse_box_p1.GetValue())
        except:
            self.error_message("p1 coarse", self.coarse_box_p1.GetValue())
            return False

        try:
            p0_fine = float(self.fine_box.GetValue())
        except:
            self.error_message("p0 fine", self.fine_box.GetValue())
            return False

        try:
            p1_fine = float(self.fine_box_p1.GetValue())
        except:
            self.error_message("p1 fine", self.fine_box_p1.GetValue())
            return False
        
        return True


    def error_message(self, box, value):
        dlg = wx.MessageDialog(
                self,
                "The value {} cannot be converted to a float for the {} box. Please edit this to a valid number and try again.".format(value, box),
                "Warning",
                wx.OK | wx.ICON_WARNING,
            )
        self.Raise()
        self.SetFocus()
        dlg.ShowModal()
        dlg.Destroy()

        