import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '/Users/alexhill/Documents/GitHub/NMR_Analysis')
from NMR_Tools import *
import matplotlib as mpl
import os 
import csv
import pandas as pd 
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

class LoadSpectra:

    def __init__(self):
        self.initial_ppm = []
        self.initial_amplitude = []
        self.water_sub_ppm = []        
        self.water_sub_amplitude = []
        self.fitting_function = []
        self.function_best_params = []
        self.fit = []
        self.filepath = []
        self.filename = []
        self.time_taken = []

    def ReadTextFile(self,
                    path = "/Users/alexhill/Desktop/Metabolomics/Data_Analysis/BrukerRerun3Data",
                    sample = 'D24',
                    pulse = 'zg30',
                    nscan = '256'):

        self.filepath = path
        self.filename = "%s_%s_ns%s.txt" % (sample, pulse, nscan)

        x, y = grab_data(path = self.filepath, file_name = self.filename)
        self.initial_ppm = x
        self.initial_amplitude = y
        self.water_sub_ppm = x
        self.water_sub_amplitude = y

    def ReadRawData(self,
                    path =  "/Users/alexhill/Desktop/Metabolomics/Rerun_Data/20240119_RERUN_3",
                    sample = 'D24',
                    pulse = 'zg30',
                    nscan = 256   
                    ):

        def process_acqus_file(folder, filename, linestart='##$NS'):
            file_path = os.path.join(folder, filename)
            try:
                with open(file_path, 'r') as file:
                    for line in file:
                        if line.startswith(linestart) and line[len(linestart)] in ('', ' ', '='):
                            _, value = line.split('=')
                            return value.strip()
                return None
            except Exception as e:
                print(f"Error processing file '{filename}': {e}")
                return None

        def extract_root_and_experiment(folder):
            # Split folder path into Root and Experiment
            root, experiment = os.path.split(folder)
            return root, experiment
            
        def search_and_save_to_csv(root_folder='.', output_csv='output.csv'):
            with open(output_csv, 'w', newline='') as csvfile:
                fieldnames = ['Root', 'Sample', 'Experiment', 'Filename', 'PULPROG', 'NS (Num. Scans)', 'TD0 (Loop count)', 'TD (Size of FID)', 'Time (seconds)']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()

                for foldername, subfolders, filenames in os.walk(root_folder):
                    for filename in filenames:
                        if filename.lower() == 'acqus':
                            ns_value = process_acqus_file(foldername, filename, linestart='##$NS')
                            td0_value = process_acqus_file(foldername, filename, linestart='##$TD0')
                            td_value = process_acqus_file(foldername, filename, linestart='##$TD')
                            pulprog_value = process_acqus_file(foldername, filename, linestart='##$PULPROG')
                            date_value = process_acqus_file(foldername, filename, linestart='##$DATE')
                            date_st_value = process_acqus_file(foldername, filename, linestart='##$DATE_START')

                            if ns_value is not None:
                                root_pr, experiment = extract_root_and_experiment(foldername)
                                root, sample = extract_root_and_experiment(root_pr)
                                writer.writerow({
                                    'Root': root,
                                    'Sample': sample,
                                    'Experiment': experiment,
                                    'Filename': filename,
                                    'PULPROG': pulprog_value,
                                    'NS (Num. Scans)': ns_value,
                                    'TD0 (Loop count)': td0_value,
                                    'TD (Size of FID)': td_value,
                                    'Time (seconds)': int(date_value) - int(date_st_value)
                                })

        def read_csv(file_path):
            """
            Read a CSV file into a pandas DataFrame.
            
            Parameters:
                file_path (str): Path to the CSV file.
                
            Returns:
                pandas.DataFrame: DataFrame containing the data from the CSV file.
            """
            df = pd.read_csv(file_path)
            return df

        # Example usage

        # Replace '.' with the root folder path if you want to start the search from a specific directory
        # Replace 'output.csv' with the desired CSV file name
        search_and_save_to_csv(path, 'temp.csv')
        df = read_csv('temp.csv')
        self.df = df
        time = df["Time (seconds)"].values
        pulses = df["PULPROG"].values
        sample_experiments = df["Sample"].values
        num_scan_per_loop = df["NS (Num. Scans)"].values
        loops = df["TD0 (Loop count)"].values
        tot_scans = loops * num_scan_per_loop

        sample_mask = np.core.defchararray.find(sample_experiments.astype(str), sample) != -1
        pulse_mask = np.core.defchararray.find(pulses.astype(str), pulse) != -1
        nscan_mask = (tot_scans == nscan)
        time_taken = time[sample_mask * pulse_mask * nscan_mask]

        if len(time_taken) == 1:
            self.time_taken = time_taken[0]
        elif len(time_taken) == 0:
            print('No match found in database')
        elif len(time_taken) > 1:
            print('Multiple matches found in database')

    def QuickPlot(self):
        "A function to quickly plot ppm and amplitude to see what's going on"

        fig, axs = plt.subplots(1,1, figsize = [6,6])
        axs.set_title(self.filename)
        axs.plot(self.initial_ppm, self.initial_amplitude, label = 'Initial data')
        axs.plot(self.water_sub_ppm, self.water_sub_amplitude, label = 'water_sub data')
        if len(self.fit) != 0:
            axs.plot(self.water_sub_ppm, self.fit, label = 'Fit', ls = '--')
        axs.set_xlabel('ppm')
        axs.set_ylabel('amplitude')
        axs.legend()
        axs.invert_xaxis()
        plt.show()
        plt.close()

    def SubtractWater(self,
                    functional_form = 'asymm_lorentzian'):


        def _lorentzian(x, A, x0, gamma):
            """
            Compute the Lorentzian function.

            Parameters:
            x : array-like
                The independent variable.
            A : float
                The amplitude of the peak.
            x0 : float
                The peak position.
            gamma : float
                The full width at half maximum (FWHM).

            Returns:
            array-like
                The value of the Lorentzian function at each point in x.
            """
            return A / (1 + ((x - x0) / (gamma / 2)) ** 2)

        def _asymm_lorentzian(x, A, x0, gamma_left, gamma_right):
            left_mask = x <= x0
            right_mask = ~left_mask  # Inverse of linear_mask

            #left_output = A / (1 + ((x[left_mask] - x0) / (gamma_left / 2)) ** 2)
            #right_output = A / (1 + ((x[right_mask] - x0) / (gamma_right / 2)) ** 2)

            left_output = _lorentzian(x[left_mask], A, x0, gamma_left)    
            right_output = _lorentzian(x[right_mask], A, x0, gamma_right)    

            output = np.empty_like(x)
            output[left_mask] = left_output
            output[right_mask] = right_output
            return output
        
        def _gaussian(x, amplitude, mu, sigma):
            return amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

        x, y = self.initial_ppm, self.initial_amplitude
        max_y = np.nanmax(y)
        x_where_max_y = x[y == max_y][0]
        amp = max_y
        mean = x_where_max_y
        fwhm = 1

        if functional_form == 'gaussian':
            initial_guess = [amp, mean, fwhm]
            popt, pcov = curve_fit(_gaussian, x, y, p0=initial_guess)
            output = _gaussian(x, popt[0], popt[1], popt[2])
            self.fitting_function = _gaussian
            self.water_sub_ppm = self.initial_ppm
            self.fit = self.fitting_function(self.water_sub_ppm, popt[0], popt[1], popt[2])

        elif functional_form == 'lorentzian':
            initial_guess = [amp, mean, fwhm]
            popt, pcov = curve_fit(_lorentzian, x, y, p0=initial_guess)
            output = _lorentzian(x, popt[0], popt[1], popt[2])
            self.fitting_function = _lorentzian
            self.water_sub_ppm = self.initial_ppm
            self.fit = self.fitting_function(self.water_sub_ppm, popt[0], popt[1], popt[2])

        elif functional_form == 'asymm_lorentzian':
            initial_guess = [amp, mean, fwhm, fwhm]
            popt, pcov = curve_fit(_asymm_lorentzian, x, y, p0=initial_guess)
            output = _asymm_lorentzian(x, popt[0], popt[1], popt[2], popt[3])
            self.fitting_function = _asymm_lorentzian
            self.water_sub_ppm = self.initial_ppm
            self.fit = self.fitting_function(self.water_sub_ppm, popt[0], popt[1], popt[2], popt[3])

        self.function_best_params = popt
        self.water_sub_amplitude = self.initial_amplitude - self.fit

class AnalyseSpectra:

    def __init__(self):
        self.initial_ppm = []
        self.initial_amplitude = []
        self.processed_ppm = []        
        self.processed_amplitude = []
        self.tsp_centre = []
        self.tsp_amplitude = []
        self.tsp_confidence = []
        self.tsp_centre = []
        self.tsp_amplitude = []
        self.tsp_integral = []
        self.tsp_confidence = []
        self.tsp_fit = []
        self.snr     = []
        self.spectra_signal = []
        self.spectra_noise  = []
        self.signal_bounds = []
        self.noise_bounds = []
    
    def InputData(self,
                x = [],
                y = []):
        self.initial_ppm = x
        self.initial_amplitude = y     
        self.processed_ppm = x
        self.processed_amplitude = y

    def FitTSP(self,
               fit_function = 'lorentzian',
               tsp_ppm_location = 0.0, 
               tsp_window = 0.2,
               return_uncertainty = True,
               plot_fit = False, 
               save_fig = False,
               plot_fnm = 'dummy_fit.png',
               title = None):

        def _lorentzian(x, A, x0, gamma):
            """
            Compute the Lorentzian function.

            Parameters:
            x : array-like
                The independent variable.
            A : float
                The amplitude of the peak.
            x0 : float
                The peak position.
            gamma : float
                The full width at half maximum (FWHM).

            Returns:
            array-like
                The value of the Lorentzian function at each point in x.
            """
            return A / (1 + ((x - x0) / (gamma / 2)) ** 2)

        def _asymm_lorentzian(x, A, x0, gamma_left, gamma_right):
            left_mask = x <= x0
            right_mask = ~left_mask  # Inverse of linear_mask

            #left_output = A / (1 + ((x[left_mask] - x0) / (gamma_left / 2)) ** 2)
            #right_output = A / (1 + ((x[right_mask] - x0) / (gamma_right / 2)) ** 2)

            left_output = _lorentzian(x[left_mask], A, x0, gamma_left)    
            right_output = _lorentzian(x[right_mask], A, x0, gamma_right)    

            output = np.empty_like(x)
            output[left_mask] = left_output
            output[right_mask] = right_output
            return output
        
        def _gaussian(x, amplitude, mu, sigma):
            return amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
        
        def _fit_singlet(XDATA, 
                        YDATA, 
                        method = 'lorentzian', 
                        window = tsp_window, 
                        expected_centre = tsp_ppm_location, 
                        show_fit = plot_fit, 
                        save = save_fig, 
                        fnm = plot_fnm,
                        title = title, 
                        return_uncertainty = return_uncertainty):

            window = window
            peak_loc = expected_centre
            x_mask = (XDATA >= peak_loc - window) * (XDATA <= peak_loc + window)
            x, y = XDATA[x_mask], YDATA[x_mask]
            max_y = np.nanmax(y)
            x_where_max_y = x[y == max_y][0]

            amp = max_y
            centre = x_where_max_y
            fwhm = 0.05
            
            if method == 'asymm_lorentzian':
                intital_guess = [amp, centre, fwhm, fwhm]
                popt, pcov = curve_fit(_asymm_lorentzian, x, y, p0=intital_guess)
                y_fit = _asymm_lorentzian(x, popt[0], popt[1], popt[2], popt[3])

            elif method == 'lorentzian':
                intital_guess = [amp, centre, fwhm]
                popt, pcov = curve_fit(_lorentzian, x, y, p0=intital_guess)
                y_fit = _lorentzian(x, popt[0], popt[1], popt[2])

            elif method == 'gaussian':
                intital_guess = [amp, centre, fwhm]
                popt, pcov = curve_fit(_gaussian, x, y, p0=intital_guess)
                y_fit = _gaussian(x, popt[0], popt[1], popt[2])

            integral_area = get_area(x, y_fit, avg = False)
            output = [x, y, y_fit, popt, integral_area]

            if return_uncertainty:
                # Calculate residuals
                residuals = y - y_fit

                # Calculate R-squared
                ss_total = np.sum((y - np.mean(y))**2)
                ss_residual = np.sum(residuals**2)
                r_squared = 1 - (ss_residual / ss_total)

                # Calculate RMSE
                rmse = np.sqrt(np.mean(residuals**2))

                #print(f"R-squared: {r_squared}")
                #print(f"RMSE: {rmse}")
                output = output + [r_squared, rmse]

            if show_fit:
                fig, axs = plt.subplots(1,2, figsize = [10,5])
                if title is None:
                    title = 'Method: %s, $R^{2} = %s$' % (method, round(r_squared, 4))
                axs[0].set_title(title)
                axs[0].plot(x, y, color = 'k', label = 'data')
                axs[0].plot(x, y_fit, color = 'red', label = 'fit: %s' % method)
                axs[0].invert_xaxis()
                axs[0].legend()
                lim1, lim2 = axs[0].get_ylim()
                axs[1].plot(x, y-y_fit, color = 'green', label = 'residual')
                axs[1].invert_xaxis()
                axs[1].legend()
                axs[1].set_ylim([lim1, lim2])
                if save:
                    if fnm is None:
                        fnm = '%s.png' % scan
                    plt.savefig(fnm)
                else:
                    plt.show()
                plt.close()

            return output

        out_tsp = _fit_singlet(self.initial_ppm, self.initial_amplitude)
        self.tsp_centre = out_tsp[3][1]
        self.tsp_amplitude = out_tsp[3][0]
        self.tsp_integral = out_tsp[4]
        self.tsp_fit = out_tsp[2]
        self.tsp_confidence = [out_tsp[5], out_tsp[6]]
        #intital_guess = [amp, centre, fwhm, fwhm]
        #output = [x, y, y_fit, popt, integral_area]

    def ScaleSpectra(self,
                     y_scaling = 1.,
                     x_shift = 0.0):

        self.processed_ppm = self.initial_ppm + x_shift
        self.processed_amplitude = self.initial_amplitude * y_scaling
  
    def QuickPlot(self,
                  plot_title = None):
        "A function to quickly plot ppm and amplitude to see what's going on"

        fig, axs = plt.subplots(1,1, figsize = [6,6])
        axs.set_title(plot_title)
        axs.plot(self.initial_ppm, self.initial_amplitude, label = 'Initial data')
        axs.plot(self.processed_ppm, self.processed_amplitude, label = 'Processed data')
        axs.set_xlabel('ppm')
        axs.set_ylabel('amplitude')
        axs.legend()
        axs.invert_xaxis()
        plt.show()
        plt.close()

    def SignalToNoise(self,
                      signal_bounds = [2.5, 2.7],
                      noise_bounds  = [-2.0, -1.0]):

        xdata = self.processed_ppm
        ydata = self.processed_amplitude

        noise_mask = (xdata > noise_bounds[0]) * (xdata < noise_bounds[1])
        sig_mask = (xdata > signal_bounds[0]) * (xdata < signal_bounds[1])

        noise = get_area(xdata[noise_mask], ydata[noise_mask], avg = True)
        sig = get_area(xdata[sig_mask], ydata[sig_mask], avg = True)
        
        self.spectra_signal = sig
        self.spectra_noise = noise
        self.snr  = sig/noise
        self.signal_bounds = sig_mask
        self.noise_bounds = noise_bounds

    def ClearSignalandNoise(self):
        self.spectra_signal = []
        self.spectra_noise  = []
        self.signal_bounds = []
        self.noise_bounds = []
