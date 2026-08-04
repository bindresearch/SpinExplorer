import wx
import numpy as np
import copy 
import matplotlib
matplotlib.use("wxAgg")
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import (
    NavigationToolbar2WxAgg as NavigationToolbar,
)
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit
from sklearn.cluster import DBSCAN

from SpinExplorer.SpinView.UI_objects.UI_tools import FloatSlider


class fit_peaks():

    def __init__(self, peaklist_frame):

        self.peaklist_frame = peaklist_frame

        if(self.perform_fitting_checks()==None):
            return
        
        self.fit_success = self.run_fit()

    def perform_fitting_checks(self):
        """
        Performing checks to ensure that peak fitting can take place
        """
        if(self.peaklist_frame.active_select_peaks==False):
            # The select group feature must have been used to select the current peaks
            # in order to ensure that the area for peak fitting is defined
            dlg = wx.MessageDialog(
                    self.peaklist_frame,
                    "The select group feature must have been used to select the current peaks in order to ensure that the area for peak fitting is defined. Please try again."
                    ,
                    "Warning",
                    wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return None
        
        if(self.peaklist_frame.selected_area==[]):
            dlg = wx.MessageDialog(
                    self.peaklist_frame,
                    "Please select an area containing peaks using he select group feature and try again."
                    ,
                    "Warning",
                    wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return None
        
        if(self.peaklist_frame.selected_peak_indexes==[]):
            dlg = wx.MessageDialog(
                    self.peaklist_frame,
                    "Please ensure that at least one peak is selected using the select group feature and try again."
                    ,
                    "Warning",
                    wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return None
        
        # Passed the checks
        return True
    

    def get_data_for_fitting(self, indexes = []):

        if(self.peaklist_frame.main_frame.multiplot_mode==False):
            data = self.peaklist_frame.main_frame.nmrdata.data * self.peaklist_frame.main_frame.multiply_factor
            x_values = self.peaklist_frame.main_frame.new_x_ppms
            y_values = self.peaklist_frame.main_frame.new_y_ppms       
        else:
            data = self.peaklist_frame.main_frame.values_dictionary[self.peaklist_frame.main_frame.active_plot_index]["z_data"] * self.peaklist_frame.main_frame.values_dictionary[self.peaklist_frame.main_frame.active_plot_index]["multiply factor"]
            x_values = self.peaklist_frame.main_frame.values_dictionary[self.peaklist_frame.main_frame.active_plot_index]["new_x_ppms"]
            y_values = self.peaklist_frame.main_frame.values_dictionary[self.peaklist_frame.main_frame.active_plot_index]["new_y_ppms"]


        # Find the indexes of x_values and y_values that are within the xlimit and ylimits of self.selected_area

        if(indexes==[]):
            xmin, xmax, ymin, ymax = self.peaklist_frame.selected_area

            xmin_index = np.argmin(np.abs(x_values-xmin))
            xmax_index = np.argmin(np.abs(x_values-xmax))
            ymin_index = np.argmin(np.abs(y_values-ymin))
            ymax_index = np.argmin(np.abs(y_values-ymax))
        
        else:
            xmin_index, xmax_index, ymin_index, ymax_index = indexes

        x1,x2 = sorted([xmin_index, xmax_index])
        y1,y2 = sorted([ymin_index, ymax_index])

        self.selected_xvalues = x_values[x1:x2]
        self.selected_yvalues = y_values[y1:y2]
        self.selected_data_values = data[x1:x2,y1:y2]


        # Plot this data as a 2D contour briefly to see if it matches
        self.X,self.Y = np.meshgrid(self.selected_yvalues, self.selected_xvalues)

        xdata = (self.Y.ravel(), self.X.ravel())
        zdata = self.selected_data_values.ravel()

        return xdata, zdata
    

    def get_data_for_fitting_cluster(self, cluster_region):

        if(self.peaklist_frame.main_frame.multiplot_mode==False):
            data = self.peaklist_frame.main_frame.nmrdata.data * self.peaklist_frame.main_frame.multiply_factor
            x_values = self.peaklist_frame.main_frame.new_x_ppms
            y_values = self.peaklist_frame.main_frame.new_y_ppms       
        else:
            data = self.peaklist_frame.main_frame.values_dictionary[self.peaklist_frame.main_frame.active_plot_index]["z_data"] * self.peaklist_frame.main_frame.values_dictionary[self.peaklist_frame.main_frame.active_plot_index]["multiply factor"]
            x_values = self.peaklist_frame.main_frame.values_dictionary[self.peaklist_frame.main_frame.active_plot_index]["new_x_ppms"]
            y_values = self.peaklist_frame.main_frame.values_dictionary[self.peaklist_frame.main_frame.active_plot_index]["new_y_ppms"]


        # Find the indexes of x_values and y_values that are within the xlimit and ylimits of self.selected_area     

        xmin, xmax, ymin, ymax = cluster_region
        xmin_index = np.argmin(np.abs(x_values-xmin))
        xmax_index = np.argmin(np.abs(x_values-xmax))
        ymin_index = np.argmin(np.abs(y_values-ymin))
        ymax_index = np.argmin(np.abs(y_values-ymax))

        x1,x2 = sorted([xmin_index, xmax_index])
        y1,y2 = sorted([ymin_index, ymax_index])

        selected_xvalues = x_values[x1:x2]
        selected_yvalues = y_values[y1:y2]
        selected_data_values = data[x1:x2,y1:y2]

        # Plot this data as a 2D contour briefly to see if it matches
        X,Y = np.meshgrid(selected_yvalues, selected_xvalues)

        xdata = (Y.ravel(), X.ravel())
        zdata = selected_data_values.ravel()


        # Assume most of the data is noise
        data_noise_value = np.percentile(data,20)

        # mask the noise so that the fit can avoid this
        mask = (selected_data_values > (data_noise_value)).ravel()
        zdata_to_fit = selected_data_values.ravel()[mask]

        xdata_to_fit = (Y.ravel()[mask], X.ravel()[mask])

        return xdata, zdata, zdata_to_fit, xdata_to_fit
    

    def get_peaks_for_fitting(self):
        """
        For each peaks, get initial estimates for the fit parameters
        
        - Initial peak position is the current stored peak position

        - Lower and upper bounds for each parameter are determined
          Intensity: 0->10*initial_intensity
          Peak positions: initial_position+/-0.1*sigma
          Sigma: 0->infinity

        - sigma_x and sigma_y estimates are determined from a 1D lineshape analysis

        - The threshold for clustering in the x/y dimensions will be sigma_x_cluster
          and sigma_y_cluster
        """

        peaks = []
        lower_bounds = []
        upper_bounds = []

        sigma_x_values = []
        sigma_y_values = []

        

        for peak_index in self.peaklist_frame.selected_peak_indexes:
            peaklist_dict = self.peaklist_frame.peak_list_dictionary[self.peaklist_frame.current_peaklist_box.GetValue()]

            x_init = peaklist_dict['shift1'][peak_index]
            y_init = peaklist_dict['shift2'][peak_index]

            x_init_arg = np.argmin(np.abs(self.selected_xvalues-x_init))
            y_init_arg = np.argmin(np.abs(self.selected_yvalues-y_init))

            row = self.selected_data_values[x_init_arg, :]
            col = self.selected_data_values[:, y_init_arg]

            intensity = peaklist_dict['intensity'][peak_index]
            if(intensity == 0):
                try:
                    intensity = self.selected_data_values[x_init_arg, y_init_arg]
                except:
                    intensity = np.max(self.selected_data_values)

            

            sigma_y_init = self.estimate_sigma_1d_asymmetric(row, y_init_arg, self.selected_yvalues)
            sigma_x_init = self.estimate_sigma_1d_asymmetric(col, x_init_arg, self.selected_xvalues)

            sigma_x_values.append(sigma_x_init)
            sigma_y_values.append(sigma_y_init)

            peaks.append([intensity, x_init, y_init, sigma_x_init, sigma_y_init])

            lower_bounds.append([0,x_init-0.25*sigma_x_init, y_init-0.25*sigma_y_init, 0, 0])
            upper_bounds.append([intensity*10,x_init+0.25*sigma_x_init, y_init+0.25*sigma_y_init, sigma_x_init*2, sigma_y_init*2])


        sigma_x_cluster = np.median(np.array(sigma_x_values))
        sigma_y_cluster = np.median(np.array(sigma_y_values))

        return peaks, lower_bounds, upper_bounds, sigma_x_cluster, sigma_y_cluster, sigma_x_values, sigma_y_values
    
    
    def cluster_centres(self, peaks, lower_bounds, upper_bounds, sigma_x_cluster, sigma_y_cluster, sigma_x_values, sigma_y_values):
        """
        eps (set to 1 as the peak positions are normalised to sigma in each dimension)
        """
        centres = []
        for i, peak in enumerate(peaks):
            sigma_x = sigma_x_cluster
            sigma_y = sigma_y_cluster
                

            centres.append([peak[1]/(sigma_x), peak[2]/(sigma_y)])

        centres = np.asarray(centres)

        clustering = DBSCAN(eps=5, min_samples=1).fit(centres)
        labels = clustering.labels_


        clusters = []
        cluster_lower_bounds = []
        cluster_upper_bounds = []
        for k in np.unique(labels):
            clusters.append(np.array(peaks)[labels == k])
            cluster_lower_bounds.append(np.array(lower_bounds)[labels==k])
            cluster_upper_bounds.append(np.array(upper_bounds)[labels==k])

        core_points = centres[clustering.core_sample_indices_]
        core_labels = clustering.labels_[clustering.core_sample_indices_]

        return clusters, cluster_lower_bounds, cluster_upper_bounds, core_points, core_labels
    

    def perform_fit(self, clusters, cluster_lower_bounds, cluster_upper_bounds, sigma_x_cluster, sigma_y_cluster):

        self.popt = []


        for k, cluster in enumerate(clusters):
            p0 = []
            upper = []
            lower = []
            for c in cluster:
                p0+=list(c)
                if(sigma_x_cluster!=999):
                    if(len(cluster)>1):
                        p0[-2] = sigma_x_cluster
                        p0[-1] = sigma_y_cluster
            for l in cluster_lower_bounds[k]:
                lower+=list(l)
            for u in cluster_upper_bounds[k]:
                upper+=list(u)

            cluster_region = self.get_cluster_region(cluster)
            xdata, zdata, zdata_to_fit, xdata_to_fit= self.get_data_for_fitting_cluster(cluster_region)

            def model(xy, *params):
                return self.n_gaussians_2d(xy, *params, n=len(cluster))
            
            popt, pcov = curve_fit(model, xdata_to_fit, zdata_to_fit, p0=p0, maxfev=10000, bounds=(lower, upper))
            self.popt+=list(popt)
            zfit = self.n_gaussians_2d((self.Y.ravel(), self.X.ravel()), *popt, n=len(cluster)).reshape(self.X.shape)# your Gaussian fitting routine
            self.z_fit += zfit
    

    def get_cluster_region(self, cluster):
        """
        Determine a suitable region to perform the fitting on based on the cluster
        initial peak positions and linewidths\
        
        If the peak positions are at (x,y)

        The x region will be defined as min(x)-xthresh*sigma_x, max(x)+xthresh*sigma
        """

        x_thresh = 1.5
        y_thresh = 1.5


        x_values = []
        sigma_x_vals = []
        y_values = []
        sigma_y_vals = []

        for i in range(len(cluster)):
            x_values.append(cluster[i][1])
            y_values.append(cluster[i][2])
            sigma_x_vals.append(cluster[i][3])
            sigma_y_vals.append(cluster[i][4])


        sigma_x = np.median(np.array(sigma_x_vals))
        sigma_y = np.median(np.array(sigma_y_vals))

        xmin = np.min(np.array(x_values))-x_thresh*sigma_x
        xmax = np.max(np.array(x_values))+x_thresh*sigma_x

        ymin = np.min(np.array(y_values))-y_thresh*sigma_y
        ymax = np.max(np.array(y_values))+y_thresh*sigma_y
            
        return [xmin, xmax, ymin, ymax]

    def estimate_sigma_1d_asymmetric(self, profile, center_idx, x, peak=None, frac=0.606):
        """
        Estimate Gaussian sigma from a 1D profile using a centre-anchored
        threshold crossing method with linear interpolation.

        Parameters
        ----------
        x : array-like
            Coordinate values (must match profile axis)
        profile : array-like
            Intensity values
        center_idx : int
            Index of the peak / Gaussian centre
        peak : float or None
            Peak value. If None, uses profile[center_idx]
        frac : float
            Fraction of peak to use (0.606 ≈ 1 sigma for Gaussian)

        Returns
        -------
        sigma : float or None
            Estimated sigma in coordinate units
        """

        x = np.asarray(x)
        y = np.asarray(profile)

        if peak is None:
            peak = y[center_idx]

        target = frac * peak
        x0 = x[center_idx]

        def find_crossing(i_start, i_end, step):
            for i in range(i_start, i_end, step):
                if i + step < 0 or i + step >= len(y):
                    break

                y1, y2 = y[i], y[i + step]

                # check for crossing
                if (y1 - target) * (y2 - target) <= 0:
                    x1, x2 = x[i], x[i + step]

                    if y2 == y1:
                        return x1

                    t = (target - y1) / (y2 - y1)
                    return x1 + t * (x2 - x1)

            return None

        # search outward from centre
        left_cross = find_crossing(center_idx, -1, -1)
        right_cross = find_crossing(center_idx, len(y) - 1, +1)

        sigmas = []

        if left_cross is not None:
            sigmas.append(x0 - left_cross)

        if right_cross is not None:
            sigmas.append(right_cross - x0)

        if len(sigmas) == 0:
            return None
        
        return np.median(np.abs(sigmas))


    def n_gaussians_2d(self, xy, *params, n):
        """
        Sum of N 2D Gaussians.
        params layout (per Gaussian): [amplitude, x0, y0, sigma_x, sigma_y]
        Total params = N*5
        """
        x, y = xy
        result = np.zeros_like(x, dtype=float)
        
        for i in range(n):
            amplitude, x0, y0, sigma_x, sigma_y = params[i*5 : i*5 + 5]
            result += amplitude * np.exp(
                -((x - x0)**2 / (2*sigma_x**2) + (y - y0)**2 / (2*sigma_y**2))
            )
        
        return result
    
    def run_fit(self, indexes=[]):
        """
        Run 2D fitting
        """
        try:
            # Get the data into a format ready for fitting
            xdata, zdata = self.get_data_for_fitting(indexes)
                
            # Find initial parameters for the fit
            peaks, lower_bounds, upper_bounds, sigma_x_cluster, sigma_y_cluster, sigma_x_values, sigma_y_values = self.get_peaks_for_fitting()


            self.sigma_x_cluster = sigma_x_cluster
            self.sigma_y_cluster = sigma_y_cluster

            # Cluster selected peaks
            clusters, cluster_lower_bounds, cluster_upper_bounds, core_points, core_labels = self.cluster_centres(peaks, lower_bounds, upper_bounds, sigma_x_cluster, sigma_y_cluster, sigma_x_values, sigma_y_values)
            

            self.fit_clusters_to_data(clusters, cluster_lower_bounds, cluster_upper_bounds)
        
        except:
            # report that the fit failed
            dlg = wx.MessageDialog(
                    self.peaklist_frame,
                    "The fit of the selected peak or peaks failed. Please ensure that the peaks are moved close to their local maxima and try again."
                    ,
                    "Warning",
                    wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return None
        
        return self.z_fit, [clusters, cluster_lower_bounds, cluster_upper_bounds, core_points, core_labels]
    


    def fit_clusters_to_data(self, clusters, cluster_lower_bounds, cluster_upper_bounds):

        # Create an empty array to fill with the fitted data
        self.z_fit = np.zeros_like(self.selected_data_values)
        
        self.perform_fit(clusters, cluster_lower_bounds, cluster_upper_bounds, self.sigma_x_cluster, self.sigma_y_cluster)


    def add_2d_fit(self):
        
        if(self.main_frame.multiplot_mode==False):
            data = self.main_frame.nmrdata.data * self.main_frame.multiply_factor
            x_values = self.main_frame.new_x_ppms
            y_values = self.main_frame.new_y_ppms       
        else:
            data = self.main_frame.values_dictionary[self.main_frame.active_plot_index]["z_data"] * self.main_frame.values_dictionary[self.main_frame.active_plot_index]["multiply factor"]
            x_values = self.main_frame.values_dictionary[self.main_frame.active_plot_index]["new_x_ppms"]
            y_values = self.main_frame.values_dictionary[self.main_frame.active_plot_index]["new_y_ppms"]


        self.selected_area = np.min(x_values)-1, np.max(x_values)+1, np.min(y_values)-1, np.max(y_values)+1

        selected_area_indexes = 0, len(x_values), 0, len(y_values)

        self.selected_peak_indexes = range(len(self.peak_list_dictionary[self.current_peaklist_box.GetValue()]['shift1']))

        print('running fit 1')
        fit = self.run_fit(selected_area_indexes)
        print('fit 1 complete')

        clusters, cluster_lower_bounds, cluster_upper_bounds, core_points, core_labels = fit[1]

        # If any peak has a sigma greater than 3x standard deviation of sigmas, then add another point here and refit
        shifts1 = []
        shifts2 = []
        intensities = []
        sigma1 = []
        sigma2 = []
        for j, peak_index in enumerate(self.selected_peak_indexes):
            shifts1.append(self.popt[j*5+1])
            shifts2.append(self.popt[j*5+2])
            intensities.append(self.popt[j*5])
            sigma1.append(self.popt[j*5+3])
            sigma2.append(self.popt[j*5+4])

        median_sigma1 = np.median(sigma1)
        median_sigma2 = np.median(sigma2)

        stdev_sigma1 = np.std(sigma1)
        stdev_sigma2 = np.std(sigma2)
        

        number = len(self.peak_list_dictionary[self.current_peaklist_box.GetValue()]["peak_name"])

        count = 0
        points = 0
        for p, point in enumerate(self.selected_peak_indexes):
            points+=1
            if(sigma1[p]> median_sigma1+2*stdev_sigma1):
                add_point = True
            elif(sigma2[p] > median_sigma2+2*stdev_sigma2):
                add_point = True
            else:
                add_point = False

            if(add_point == True):
                name = str(number+1+count)
                points+=1
                self.peak_list_dictionary[self.current_peaklist_box.GetValue()]["peak_name"].append(name)
                self.peak_list_dictionary[self.current_peaklist_box.GetValue()]["shift1"].append(self.popt[p*5+1])
                self.peak_list_dictionary[self.current_peaklist_box.GetValue()]["shift2"].append(self.popt[p*5+2])
                self.peak_list_dictionary[self.current_peaklist_box.GetValue()]["intensity"].append(self.popt[p*5])
                count+=1




        new_points = []

        if(fit != None):
            data_minus_fit = data-fit[0]
            
            # Find the data regions where data_minus_fit is greater than the threshold
            indexes = np.argwhere(data_minus_fit> np.max(data)*float(self.peak_picking_threshold_box.GetValue())/100)

            labels = []
            print(len(indexes))
            for i, index in enumerate(indexes):
                xval = x_values[index[0]]/(5*self.sigma_x_cluster)
                yval = y_values[index[1]]/(5*self.sigma_y_cluster)

                # Find the cluster which already has a peak closest to the argwhere index
                new_point_values = [x_values[index[0]], y_values[index[1]]]
                assigned_cluster = self.assign_new_point([xval, yval], core_points, core_labels)
                if(assigned_cluster not in labels):
                    new_points.append(new_point_values)
                
                labels.append(assigned_cluster)

        print(len(new_points))
        for point in new_points:
            name = str(number+1+count)
            points+=1
            self.peak_list_dictionary[self.current_peaklist_box.GetValue()]["peak_name"].append(name)
            new_points.append(name)
            self.peak_list_dictionary[self.current_peaklist_box.GetValue()]["shift1"].append(point[0])
            self.peak_list_dictionary[self.current_peaklist_box.GetValue()]["shift2"].append(point[1])
            self.peak_list_dictionary[self.current_peaklist_box.GetValue()]["intensity"].append(0)
            count+=1


        self.selected_peak_indexes = range(len(self.peak_list_dictionary[self.current_peaklist_box.GetValue()]['shift1']))
        print('running fit 2')
        fit_new = self.run_fit(selected_area_indexes)
        print('fit 2 complete')

        for k in range(points):
            self.peak_list_dictionary[self.current_peaklist_box.GetValue()]["shift1"][k] = self.popt[k*5+1]
            self.peak_list_dictionary[self.current_peaklist_box.GetValue()]["shift2"][k] = self.popt[k*5+2]
            self.peak_list_dictionary[self.current_peaklist_box.GetValue()]["intensity"][k] = self.popt[k*5]


        self.selected_peak_indexes = []


    def assign_new_point(self, new_point, core_points, core_labels, eps=5.0):
        """
        Assign a new point to the nearest DBSCAN cluster,
        or label it as noise (-1) if no core point is within eps.
        """
        distances = np.linalg.norm(core_points - new_point, axis=1)
        nearest_idx = np.argmin(distances)

        if distances[nearest_idx] <= eps:
            return core_labels[nearest_idx]
        else:
            return -1  # Noise

class fit_peaks_2D_window(wx.Frame):
    def __init__(self, title, parent, fit_selected):
        """
        This class will pop out a window showing the result of the current peak
        fitting and it will ask the user if they wish to use the updated peak 
        positions and intensities
        """
        self.peaklist_frame = parent
        self.fit_selected = fit_selected
        self.monitorWidth, self.monitorHeight = wx.GetDisplaySize()
        width = int(self.monitorWidth/2)
        height = int(self.monitorHeight/1.25)
        wx.Frame.__init__(self, parent=parent, title=title, size=(width, height))
        self.panel_fit = wx.Panel(self, -1)
        self.main_fit_sizer = wx.BoxSizer(wx.VERTICAL)

        self.panel = wx.Panel(self)
        self.fig = Figure()
        self.canvas = FigCanvas(self, -1, self.fig)
        self.toolbar = NavigationToolbar(self.canvas)
        self.main_fit_sizer.Add(self.canvas, 10, wx.EXPAND)
        self.main_fit_sizer.Add(self.toolbar, 0, wx.EXPAND)
        self.make_fit_window()
        self.plot_fit()

        self.SetSizer(self.main_fit_sizer)

        
        self.Show()

        
    def make_fit_window(self):

        self.options_sizer = wx.BoxSizer(wx.HORIZONTAL)
        # Create a sizer for changing the contour levels of the spectrum
        self.contour_label = wx.StaticBox(self, -1, "Contour Start = max(data)/x")
        self.contour_sizer = wx.StaticBoxSizer(self.contour_label, wx.VERTICAL)
        self.csizer = wx.BoxSizer(wx.HORIZONTAL)
        self.x_val = 10.00
        self.contour2_label = wx.StaticText(self.contour_label, label="x:")
        self.contour_slider = FloatSlider(
            self.contour_label, id=-1, value=1, minval=0, maxval=3, res=0.1, size=(200, 25)
        )
        self.contour_slider.Bind(wx.EVT_SLIDER, self.OnMinContourFit)
        self.csizer.Add(self.contour2_label)
        self.csizer.AddSpacer(5)
        self.csizer.Add(self.contour_slider)

        self.contour_sizer.Add(self.csizer)
        self.options_sizer.Add(self.contour_sizer, 0, wx.ALIGN_CENTER_VERTICAL)


        # Producing buttons asking whether to accept changes to current peaklist or accept changes to a new peaklist
        self.accept_button = wx.Button(self, label='Accept updated positions/intensities')
        self.accept_button.Bind(wx.EVT_BUTTON, self.OnAcceptButton)

        self.options_sizer.AddSpacer(10)
        self.options_sizer.Add(self.accept_button, 0, wx.ALIGN_CENTER_VERTICAL)


        self.main_fit_sizer.Add(self.options_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)

    
    def OnAcceptButton(self, event):
        """
        Accept changes from the fit (peak locations and intensities) into the
        current peaklist
        """
        # Save the previous state of the peaklists (in case the user wants to undo the changes)
        if(len(self.peaklist_frame.previous_peaklists)>10):
            self.peaklist_frame.previous_peaklists.pop(0)
        self.peaklist_frame.previous_peaklists.append(copy.deepcopy(self.peaklist_frame.peak_list_dictionary))

        # Update the chemical shifts and intensities of these peaks in the current peaklist


        try:
            for j in range(len(self.selected_peak_indexes)):
                selected_peak_index = self.selected_peak_indexes[j]
                shift1 = self.shifts1[j]
                shift2 = self.shifts2[j]
                intensity = self.intensities[j]
                

                self.peaklist_frame.peak_list_dictionary[self.current_peaklist]['shift1'][selected_peak_index] = shift1
                self.peaklist_frame.peak_list_dictionary[self.current_peaklist]['shift2'][selected_peak_index] = shift2
                self.peaklist_frame.peak_list_dictionary[self.current_peaklist]['intensity'][selected_peak_index] = intensity

            self.peaklist_frame.AddToTable()
            self.peaklist_frame.main_frame.OnMinContour2D(wx.EVT_BUTTON, textcontrol=True)
            

        except:
            # Peak positions and intensities were not updated correctly
            dlg = wx.MessageDialog(
                    self,
                    "The peak positions and intensities in the peaklist could not be updated"
                    ,
                    "Warning",
                    wx.OK,
            )
            dlg.ShowModal()
            dlg.Destroy()
            return None

    def OnMinContourFit(self, event):
        self.x_val = 10 ** float(self.contour_slider.GetValue())

        xlim, ylim = self.ax.get_xlim(), self.ax.get_ylim()

        self.ax.clear()
        contour_num = 20  # number of contour levels
        contour_factor = 1.20  # scaling factor between contour levels
        contour_start = np.max(np.abs(self.data_values)) / self.x_val
        cl = contour_start * contour_factor ** np.arange(
            contour_num
        )

        contour1 = self.ax.contour(
            self.Y,
            self.X,
            self.data_values,
            cl,
            colors=self.cmap_data,
            linewidths=self.peaklist_frame.main_frame.linewidth,
            zorder=1,
            label='Fit'
        )
        contour1_fit = self.ax.contour(
            self.Y,
            self.X,
            self.z_fit,
            cl,
            colors=self.cmap_fit,
            linewidths=self.peaklist_frame.main_frame.linewidth,
            zorder=1,
            label='Fit'
        )


        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)

        self.ax.legend(self.custom_lines, self.custom_labels)

        self.ax.set_title('Peaklist = ' + self.current_peaklist)



        # Plot Peaklists
        
        cs = 'k'
        self.shifts1 = []
        self.shifts2 = []
        self.intensities = []
        for j, peak_index in enumerate(self.selected_peak_indexes):
            self.shifts1.append(self.popt[j*5+1])
            self.shifts2.append(self.popt[j*5+2])
            self.intensities.append(self.popt[j*5])

        self.points = self.ax.scatter(
                self.shifts1,
                self.shifts2,
                s=5,
                marker="o",
                c=cs,
                picker=5,
                zorder=2,
        )

        # Annotation for hover
        self.annotations= self.ax.annotate(
                    "",
                    xy=(0, 0),
                    xytext=(15, 15),
                    textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->"),
                )
        self.annotations.set_visible(False)
        self.hover_connect = self.canvas.mpl_connect(
            "motion_notify_event", self.on_hover)
            

        self.UpdateFrame()


    def on_hover(self, event):

        if event.inaxes != None:
            # Calculate distance from mouse to each point

                cont, ind = self.points.contains(event)

                if cont:
                    # Show annotation
                    index = ind["ind"][0]  # first index found
                    dictionary = self.peaklist_frame.peak_list_dictionary[self.current_peaklist]


                    peakname = (
                        dictionary["peak_name"][self.selected_peak_indexes[index]]
                    )
                    x = self.shifts1[index]
                    y = self.shifts2[index]
                    self.annotations.xy = (x, y)

                    text = peakname
                    self.annotations.set_text(text)
                    self.annotations.set_color('k')
                    self.annotations.set_position((36,0))
                    self.annotations.set_visible(True)
                    self.canvas.draw_idle()
                else:
                    if self.annotations.get_visible():
                        self.annotations.set_visible(False)

                self.canvas.draw_idle()

    def plot_fit(self):
        self.ax = self.fig.add_subplot(111)


        # Making copies of the main variables to keep this new frame independent of the main frame
        self.Y = copy.deepcopy(self.fit_selected.Y)
        self.X = copy.deepcopy(self.fit_selected.X)
        self.selected_peak_indexes = copy.deepcopy(self.peaklist_frame.selected_peak_indexes)
        self.data_values = copy.deepcopy(self.fit_selected.selected_data_values)
        self.z_fit = copy.deepcopy(self.fit_selected.z_fit)
        self.popt = copy.deepcopy(self.fit_selected.popt)
        self.current_peaklist = self.peaklist_frame.current_peaklist_box.GetValue()


        contour_num = 20  # number of contour levels
        contour_factor = 1.20  # scaling factor between contour levels
        contour_start = np.max(np.abs(self.data_values)) / self.x_val
        cl = contour_start * contour_factor ** np.arange(
            contour_num
        )

        if(self.peaklist_frame.main_frame.multiplot_mode==True):
            self.cmap_data = self.peaklist_frame.main_frame.twoD_colours[self.peaklist_frame.main_frame.active_plot_index]
        else:
            self.cmap_data = self.peaklist_frame.main_frame.twoD_colours[0]

        self.cmap_fit = 'grey'

        if(self.peaklist_frame.main_frame.multiplot_mode==True):
            self.data_label = self.peaklist_frame.main_frame.files.custom_labels[self.peaklist_frame.main_frame.active_plot_index]
        else:
            self.data_label = 'Data'


        self.custom_lines = [Line2D([0],[0], color=self.cmap_data, lw=1.5), Line2D([0],[0], color=self.cmap_fit, lw=1.5)]
        self.custom_labels = [self.data_label, 'Fit']

        contour1 = self.ax.contour(
            self.Y,
            self.X,
            self.data_values,
            cl,
            colors=self.cmap_data,
            linewidths=self.peaklist_frame.main_frame.linewidth,
            zorder=1,
        )
        contour1_fit = self.ax.contour(
            self.Y,
            self.X,
            self.z_fit,
            cl,
            colors=self.cmap_fit,
            linewidths=self.peaklist_frame.main_frame.linewidth,
            zorder=1,
        )

        xlim, ylim = self.ax.get_xlim(), self.ax.get_ylim()

        self.ax.set_xlim([max(xlim), min(xlim)])
        self.ax.set_ylim([max(ylim), min(ylim)])

        self.ax.legend(self.custom_lines, self.custom_labels)

        self.ax.set_title('Peaklist = '+ self.current_peaklist)


        self.OnMinContourFit(wx.EVT_SLIDER)


    def UpdateFrame(self):

        self.canvas.draw()
        self.canvas.Refresh()
        self.canvas.Update()
        self.panel.Refresh()
        self.panel.Update()