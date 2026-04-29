    
import wx
import os
import numpy as np
import copy 
import matplotlib
matplotlib.use("wxAgg")
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import (
    NavigationToolbar2WxAgg as NavigationToolbar,
)




class analysis_frame(wx.Frame):

    def __init__(self, parent, df):
        
        self.df = copy.deepcopy(df)
        self.parent = parent

        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        width = int(self.monitorWidth/1.25)
        height = int(self.monitorHeight/1.25)
        wx.Frame.__init__(self, parent=parent, title='Analysis', size=(width, height))
        self.panel_fit = wx.Panel(self, -1)
        self.main_analysis_sizer = wx.BoxSizer(wx.VERTICAL)

        self.panel = wx.Panel(self)
        self.fig = Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        self.toolbar = NavigationToolbar(self.canvas)
        self.main_analysis_sizer.Add(self.canvas, 10, wx.EXPAND)
        self.main_analysis_sizer.Add(self.toolbar, 0, wx.EXPAND)
        self.create_window()
        self.get_data()
        self.produce_plot()

        self.SetSizer(self.main_analysis_sizer)
        self.Show()


    def create_window(self):

        self.main_sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.multiplication_factors_label = wx.StaticBox(self, -1, "CSP multiplication factors:")
        self.multiplication_factors_sizer = wx.StaticBoxSizer(self.multiplication_factors_label, wx.HORIZONTAL)

        self.factor1 = wx.StaticText(self, -1, 'Factor ({}):'.format(self.parent.main_frame.ax.get_xlabel()))
        self.factor1_input = wx.TextCtrl(self, value = str(self.parent.factor1), size = (50, 30), style = wx.TE_PROCESS_ENTER)
        self.factor1_input.Bind(wx.EVT_TEXT_ENTER, self.OnUpdateHeteroatomMultiplicationFactor)

        self.factor2 = wx.StaticText(self, -1, 'Factor ({}):'.format(self.parent.main_frame.ax.get_ylabel()))
        self.factor2_input = wx.TextCtrl(self, value = str(self.parent.factor2), size = (50, 30), style = wx.TE_PROCESS_ENTER)
        self.factor2_input.Bind(wx.EVT_TEXT_ENTER, self.OnUpdateHeteroatomMultiplicationFactor)

        self.multiplication_factors_sizer.AddSpacer(5)
        self.multiplication_factors_sizer.Add(self.factor1)
        self.multiplication_factors_sizer.AddSpacer(5)
        self.multiplication_factors_sizer.Add(self.factor1_input)
        self.multiplication_factors_sizer.AddSpacer(5)
        self.multiplication_factors_sizer.Add(self.factor2)
        self.multiplication_factors_sizer.AddSpacer(5)
        self.multiplication_factors_sizer.Add(self.factor2_input)


        self.threshold_label = wx.StaticBox(self, -1, "Threshold lines:")
        self.threshold_sizer = wx.StaticBoxSizer(self.threshold_label, wx.HORIZONTAL)

        # Initial CSP value is 0.005
        self.threshold1 = wx.StaticText(self, -1, 'CSP threshold (ppm):')
        self.threshold1_input = wx.TextCtrl(self, value = str(0.005), size = (50, 30), style = wx.TE_PROCESS_ENTER)
        self.threshold1_input.Bind(wx.EVT_TEXT_ENTER, self.OnUpdateCSPThreshold)

        # Initial intensity value is 0.1
        self.threshold2 = wx.StaticText(self, -1, 'Intensity Ratio Threshold:')
        self.threshold2_input = wx.TextCtrl(self, value = str(0.1), size = (50, 30), style = wx.TE_PROCESS_ENTER)
        self.threshold2_input.Bind(wx.EVT_TEXT_ENTER, self.OnUpdateIntensityRatioThreshold)

        self.threshold_sizer.AddSpacer(5)
        self.threshold_sizer.Add(self.threshold1)
        self.threshold_sizer.AddSpacer(5)
        self.threshold_sizer.Add(self.threshold1_input)
        self.threshold_sizer.AddSpacer(5)
        self.threshold_sizer.Add(self.threshold2)
        self.threshold_sizer.AddSpacer(5)
        self.threshold_sizer.Add(self.threshold2_input)


        # Download CSV button

        self.download_button = wx.Button(self, label='Download CSV')
        self.download_button.Bind(wx.EVT_BUTTON, self.download_csv)


        # Add all sizers together

        self.main_sizer.Add(self.multiplication_factors_sizer)
        self.main_sizer.AddSpacer(10)
        self.main_sizer.Add(self.threshold_sizer)
        self.main_sizer.AddSpacer(10)
        self.main_sizer.Add(self.download_button, 0, wx.ALIGN_CENTER_VERTICAL)

        self.main_analysis_sizer.Add(self.main_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.main_analysis_sizer.AddSpacer(10)


    def get_data(self):
        """
        Calculate CSPs based on the multiplication factors provided

        CSP = sqrt(factor1*(delta_shift1**2) + factor2*(delta_shift2**2))
        """

        factor1 = float(self.factor1_input.GetValue())
        factor2 = float(self.factor2_input.GetValue())

        delta_shift1 = self.df['shift1 1 (ppm)']-self.df['shift1 2 (ppm)']
        delta_shift2 = self.df['shift2 1 (ppm)']-self.df['shift2 2 (ppm)']

        self.df['CSP (ppm)'] = np.sqrt(factor1*(delta_shift1**2)+factor2*(delta_shift2)**2)


    def produce_plot(self):
        self.ax = self.fig.add_subplot(211)
        self.ax2 = self.fig.add_subplot(212)

        self.add_plot()

    def add_plot(self):

        self.residue_numbers = self.df['residue number']
        self.csps = self.df['CSP (ppm)']
        self.intensity_ratios = self.df['Intensity Ratio']

        self.names = self.df['peak name']

        
        self.csp_points = self.ax.scatter(self.residue_numbers, self.csps, picker=5)
        self.ax.plot(self.residue_numbers, self.csps)
        self.ax.set_ylabel('CSP (ppm)')
        self.ax.set_title(r'CSP = $\sqrt{f_1(\Delta\delta_1^2) + f_2(\Delta\delta_2^2)}$')
        self.ax.axhline(float(self.threshold1_input.GetValue()), color='k', linestyle = '--')

        
        self.intensity_ratio_points = self.ax2.scatter(self.residue_numbers, self.intensity_ratios, picker=5)
        self.ax2.plot(self.residue_numbers, self.intensity_ratios)
        self.ax2.set_xlabel('residue number')
        self.ax2.set_ylabel('Intensity Ratio')
        self.ax2.set_title('I(peaklist1)/I(peaklist2)')
        self.ax2.axhline(1+float(self.threshold2_input.GetValue()), color='k', linestyle = '--')
        self.ax2.axhline(1-float(self.threshold2_input.GetValue()), color='k', linestyle = '--')


        # Have the ability to hover and show which residues are under the cursor

        # Annotation for hover
        self.annotations_ax = self.ax.annotate(
                "",
                xy=(0, 0),
                xytext=(15, 15),
                textcoords="offset points",
                bbox=dict(boxstyle="round", fc="w"),
                arrowprops=dict(arrowstyle="->"),
            )
        self.annotations_ax.set_visible(False)

        self.annotations_ax2 = self.ax2.annotate(
                "",
                xy=(0, 0),
                xytext=(15, 15),
                textcoords="offset points",
                bbox=dict(boxstyle="round", fc="w"),
                arrowprops=dict(arrowstyle="->"),
            )
        self.annotations_ax2.set_visible(False)

        # Connect event
        self.hover_connect = self.canvas.mpl_connect(
            "motion_notify_event", self.on_hover
        )

        self.UpdateFrame()


    def on_hover(self, event):

        if event.inaxes != None:
            # Calculate distance from mouse to each point

            if(event.inaxes == self.ax):
                cont, ind = self.csp_points.contains(event)
                if cont:
                    # Show annotation
                    index = ind["ind"][0]  # first index found
                    name = self.names[index]
                    x = self.residue_numbers[index]
                    y = self.csps[index]
                    self.annotations_ax.xy = (x, y)
                    text = name
                    self.annotations_ax.set_text(text)
                    self.annotations_ax.set_color('k')
                    self.annotations_ax.set_position((36, 0))
                    self.annotations_ax.set_visible(True)
                    self.canvas.draw_idle()
                else:
                    if self.annotations_ax.get_visible():
                        self.annotations_ax.set_visible(False)
            else:
                if self.annotations_ax.get_visible():
                        self.annotations_ax.set_visible(False)

            if(event.inaxes == self.ax2):
                cont, ind = self.intensity_ratio_points.contains(event)
                if cont:
                    # Show annotation
                    index = ind["ind"][0]  # first index found
                    name = self.names[index]
                    x = self.residue_numbers[index]
                    y = self.intensity_ratios[index]
                    self.annotations_ax2.xy = (x, y)
                    text = name
                    self.annotations_ax2.set_text(text)
                    self.annotations_ax2.set_color('k')
                    self.annotations_ax2.set_position((36, 0))
                    self.annotations_ax2.set_visible(True)
                    self.canvas.draw_idle()
                else:
                    if self.annotations_ax2.get_visible():
                        self.annotations_ax2.set_visible(False)
            else:
                if self.annotations_ax2.get_visible():
                        self.annotations_ax2.set_visible(False)


            self.canvas.draw_idle()


    def UpdateFrame(self):

        self.canvas.draw()
        self.canvas.Refresh()
        self.canvas.Update()
        self.panel.Refresh()
        self.panel.Update()


    def OnUpdateHeteroatomMultiplicationFactor(self, event):
        """
        Recalculate the CSPs and update the plots
        """
        
        # Performing checks that the multiplication factors are numbers
        try:
            float(self.factor1_input.GetValue())
            float(self.factor2_input.GetValue())
        except:
            dlg = wx.MessageDialog(
                    self,
                    "One of the multiplication factors could not be converted to a number. Please ensure only numbers are entered and try again."
                    ,
                    "Warning",
                    wx.OK,
            )
            dlg.ShowModal()  
            dlg.Destroy()
            return


        # Clear the plots to be reloaded
        self.ax.clear()
        self.ax2.clear()


        self.get_data()
        self.add_plot()


    def OnUpdateCSPThreshold(self, event):
        """
        Update the horizontal line on the plot
        Work out which residues are outside this threshold (e.g. 0.1 would mean above
        0.1 is outside the threshold) and then update the plot highlighting
        these residue numbers as text next to their points
        """

        try:
            float(self.threshold1_input.GetValue())
        except:
            dlg = wx.MessageDialog(
                    self,
                    "The CSP threshold input could not be converted to a number. Please ensure only numbers are entered and try again."
                    ,
                    "Warning",
                    wx.OK,
            )
            dlg.ShowModal()  
            dlg.Destroy()
            return
        
        # Clear the plots to be reloaded
        self.ax.clear()
        self.ax2.clear()


        self.get_data()
        self.add_plot()

    def OnUpdateIntensityRatioThreshold(self, event):
        """
        Update the horizontal line on the plot
        Work out which residues are outside this threshold (e.g. 0.1 would mean below
        0.9 and above 1.1 are outside the threshold) and then update the plot highlighting
        these residue numbers as text next to their points
        """

        try:
            float(self.threshold2_input.GetValue())
        except:
            dlg = wx.MessageDialog(
                    self,
                    "The intensity threshold input could not be converted to a number. Please ensure only numbers are entered and try again."
                    ,
                    "Warning",
                    wx.OK,
            )
            dlg.ShowModal()  
            dlg.Destroy()
            return
        
        # Clear the plots to be reloaded
        self.ax.clear()
        self.ax2.clear()


        self.get_data()
        self.add_plot()


    def download_csv(self, event):
        """
        Download a csv showing the CSPs and intensity ratios calculated
        """

        dlg = wx.FileDialog(self, "Select the folder and name to save the dataframe as.", wildcard="", style=wx.FD_SAVE)
        dlg.SetDirectory(os.getcwd())
        dlg.SetFilename('CSP_Analysis.csv')
        if dlg.ShowModal() == wx.ID_OK:
            file = dlg.GetPath()
        else:
            dlg.Destroy()
            return
        
        self.df.to_csv(file)