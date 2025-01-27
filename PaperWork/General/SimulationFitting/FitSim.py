import sys
import os 
import matplotlib as mpl
import csv
import pandas as pd 
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import h5py as h5
from scipy import stats
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '/Users/alexhill/Documents/GitHub/NMR_Analysis/Fourier90Runs')
parent_dir = os.path.abspath(os.path.join(os.getcwd(), '..'))
sys.path.append(parent_dir)
from NMRClassesH5 import *
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d



S = LoadSpectra()

concs = ["750uM", "1500uM", "3000uM", "6000uM"]

X_EXP = []
Y_EXP = []

X_SIM = []
Y_SIM = []

for conc in concs:

    S.ReadHDF5(path = "/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/General/DataSource",
                    filename = "Data.hdf5",
                    substance = 'Lactate',
                    concentration = conc,
                    field_strength = "80MHz",
                    pulse_seq = 'zg')

    x = S.initial_ppm[0][-1]
    y = S.initial_amplitude[0][-1]
    if x[0] > x[-1]:
        x = x[::-1]
        y = y[::-1]    





    A = AnalyseSpectra()
    A.InputData(x=x, y=y)
    A.FitTSP(plot_fit = False, fit_function = 'lorentzian')
    A.ScaleSpectra(x_shift = (-1) * A.tsp_centre)
    ppm_exp = A.processed_ppm
    amp_exp = A.processed_amplitude
    if ppm_exp[0] > ppm_exp[-1]:
        ppm_exp = ppm_exp[::-1]
        amp_exp = amp_exp[::-1]    


    path = '/Users/alexhill/Software/ccpnmr3.2.0/lactate31_nonoise_80MHz.hdf5'

    u = h5.File(path, 'r')
    ppm_sim = u['80.0/x'][()]
    amp_sim = u['80.0/toty'][()]
    u.close()


    def preprocessing(xy_sim, xy_exp):
        clip = [0.5, 2.5]
        x_sim, y_sim = xy_sim[0], xy_sim[1]
        x_exp, y_exp = xy_exp[0], xy_exp[1]
        
        x_sim, y_sim = x_sim[(x_sim >= clip[0]) * (x_sim <= clip[1])], y_sim[(x_sim >= clip[0]) * (x_sim <= clip[1])]
        x_exp, y_exp = x_exp[(x_exp>=clip[0]) * (x_exp<=clip[1])], y_exp[(x_exp>=clip[0]) * (x_exp<=clip[1])]

        yscale = np.nanmax(y_exp)/np.nanmax(y_sim)
        y_sim = y_sim #* yscale

        #fig, axs = plt.subplots(1,2)
        #axs[0].plot(x_sim, y_sim)
        #axs[1].plot(x_exp, y_exp)
        #plt.show()
        if x_sim[0] > x_sim[-1]:
            x_sim = x_sim[::-1]
            y_sim = y_sim[::-1]    
        if x_exp[0] > x_exp[-1]:
            x_exp = x_exp[::-1]
            y_exp = y_exp[::-1]

        y_sim_interp = np.interp(x_exp, x_sim, y_sim)

        return [x_sim, y_sim_interp], [x_exp, y_exp]

    xy_sim_pp, xy_exp_pp = preprocessing([ppm_sim, amp_sim], [ppm_exp, amp_exp])

    # Main execution
    x_data = xy_exp_pp[0]  # Your x-axis data
    simulated_data = xy_sim_pp[1]  # Your preprocessed simulated NMR data
    experimental_data = xy_exp_pp[1]  # Your preprocessed experimental NMR data

    def fitting_function(x, x_shift, y_scale, y_offset):
        """
        Fitting function that accounts for both x-axis and y-axis adjustments
        """
        return y_scale * np.interp(x + x_shift, x_data, simulated_data) + y_offset
    def fit_nmr_data(x_data, simulated_data, experimental_data):
        """
        Perform best fit of simulated NMR data to experimental data
        """
        try:
            # Initial guess for parameters (x_shift, y_scale, y_offset)
            initial_guess = [0.0, 1.0, 0.0] 
            
            # Perform the curve fit
            optimal_params, covariance = curve_fit(
                fitting_function, x_data, experimental_data, 
                p0=initial_guess
            )
            
            # Generate the fitted data
            fitted_data = fitting_function(x_data, *optimal_params)
            
            return optimal_params, fitted_data
        
        except Exception as e:
            print(f"An error occurred during fitting: {e}")
            return None, None

    # Example usage
    def plot_fit_results(x_data, simulated_data, experimental_data, fitted_data, optimal_params, filename):
        """
        Visualize the fitting results
        """
        plt.figure(figsize=(10, 6))
        plt.plot(x_data, experimental_data, label='Experimental Data', color='red')
        plt.plot(x_data, simulated_data, label='Original Simulated Data', color='blue', linestyle='--')
        plt.plot(x_data, fitted_data, label='Fitted Simulated Data', color='green')
        
        plt.title('NMR Data Fitting')
        plt.xlabel('X-axis')
        plt.ylabel('Intensity')
        plt.legend()
        
        # Add optimal parameters to the plot
        plt.text(0.05, 0.95, f'X-scale: {optimal_params[0]:.4f}\nX-shift: {optimal_params[0]:.4f}\nY-scale: {optimal_params[1]:.4f}\nY-offset: {optimal_params[2]:.4f}', 
                transform=plt.gca().transAxes, 
                verticalalignment='top')
        
        plt.tight_layout()
        plt.savefig(filename)
        plt.close()

    # Perform the fit
    optimal_params, fitted_data = fit_nmr_data(x_data, simulated_data, experimental_data)

    # Visualize results
    if optimal_params is not None:
        plot_fit_results(x_data, simulated_data, experimental_data, fitted_data, optimal_params, conc+'_fit.png')


    # Apply the fitting parameters to the original simulated data
    x_shift, y_scale, y_offset = optimal_params

    # Generate the fitted simulated data on the original x-range
    # Plot the results
    plt.figure(figsize=(10, 6))
    plt.plot(ppm_exp, amp_exp, label='Experimental Data', color='red')
    plt.plot(ppm_sim, amp_sim, label='Original Simulated Data', color='blue', linestyle='--')
    plt.plot(ppm_sim-x_shift, amp_sim * y_scale + y_offset, label='Fitted Simulated Data', color='green')
    #plt.plot(ppm_sim-x_shift, amp_sim * y_scale + y_offset, label='Fitted Simulated Data', color='limegreen')
    plt.title('NMR Data Fitting')
    plt.xlabel('ppm')
    plt.ylabel('Intensity')
    plt.legend()
    plt.tight_layout()
    plt.savefig(conc+'.png')
    plt.close()
    X_EXP.append(ppm_exp)
    Y_EXP.append(amp_exp)
    X_SIM.append(ppm_sim)
    Y_SIM.append(amp_sim * y_scale + y_offset)


fig, axs = plt.subplots(1,2)

for i in range(0, len(X_EXP)):
    axs[0].plot(X_EXP[i], Y_EXP[i])
    axs[1].plot(X_SIM[i], Y_SIM[i])

plt.show()
'''
x_scale, x_shift, y_scale, y_offset = optimal_params[0], optimal_params[1], optimal_params[2], optimal_params[3]


y_scale * np.interp(x_scale * ppm_sim + x_shift, ppm_exp, amp_sim) + y_offset



plt.plot(ppm_sim, amp_sim)
plt.plot(ppm_exp, amp_exp)
plt.plot(ppm_exp,    y_scale * np.interp(x_scale * ppm_sim + x_shift, ppm_exp, amp_sim) + y_offset)
    
    ppm_sim * x_scale + x_shift, amp_sim * y_scale + y_offset, lw = 0.5)
#plt.plot(ppm_sim - optimal_params[1], amp_sim * optimal_params[0], lw = 0.1)

plt.show()
#[ppm_sim, amp_sim], [ppm_exp, amp_exp]


'''


'''


intensity_theory_interp = xy_sim_pp[1]
ppm_exp, intensity_exp = xy_exp_pp[0], xy_exp_pp[1]

# Define the fitting model
def fit_model(ppm, scale, shift, fwhm):
    # Apply shift
    shifted_ppm = ppm - shift
    # Interpolate shifted theoretical data
    shifted_intensity = np.interp(ppm, ppm_exp - shift, intensity_theory_interp)
    # Apply Gaussian broadening
    broadened_intensity = gaussian_filter1d(shifted_intensity, sigma=fwhm / (2.355 * (ppm[1] - ppm[0])))
    return scale * broadened_intensity

# Perform the fitting
popt, pcov = curve_fit(fit_model, ppm_exp, intensity_exp, p0=[1, 0, 0.1])

# Extract fitted parameters
scale, shift, fwhm = popt

# Generate the fitted spectrum
fitted_spectrum = fit_model(ppm_exp, *popt)

# Plot the results
import matplotlib.pyplot as plt
plt.plot(ppm_exp, intensity_exp, label='Experimental')
plt.plot(ppm_exp, fitted_spectrum, label='Fitted', linestyle='--')
plt.plot(ppm_exp, intensity_theory_interp, label='Theoretical', linestyle=':')
plt.legend()
plt.gca().invert_xaxis()
plt.show()




plt.figure()
plt.scatter(xy_exp_pp[0], xy_exp_pp[1])
plt.scatter(xy_exp_pp[0], xy_sim_pp[1], alpha = 0.1)
#plt.scatter(xy_exp_pp[0], intensity_theory_interp)

plt.show()
'''


'''

#_,_,_ = getSNRs(snr_method = 'brk', plot_spectra = False, pulse_seq = 'zg30', substance = 'Glucose')
_,_,_ = getSNRs(snr_method = 'liv', plot_spectra = True, pulse_seq = 'zg', substance = 'Glucose', tsp_fit = 'lorentzian', plot_snrs = True, plot_name = 'glc_zg_liv.png')
_,_,_ = getSNRs(snr_method = 'liv', plot_spectra = True, pulse_seq = 'zg', substance = 'Citrate', tsp_fit = 'lorentzian', plot_snrs = True, plot_name = 'cit_zg_liv.png')
_,_,_ = getSNRs(snr_method = 'liv', plot_spectra = True, pulse_seq = 'zg', substance = 'Lactate', tsp_fit = 'lorentzian', plot_snrs = True, plot_name = 'lac_zg_liv.png')

_,_,_ = getSNRs(snr_method = 'liv', plot_spectra = False, pulse_seq = 'zg30', substance = 'Glucose', tsp_fit = 'lorentzian', plot_snrs = True, plot_name = 'glc_zg30_liv.png')
_,_,_ = getSNRs(snr_method = 'liv', plot_spectra = False, pulse_seq = 'zg30', substance = 'Citrate', tsp_fit = 'lorentzian', plot_snrs = True, plot_name = 'cit_zg30_liv.png')
_,_,_ = getSNRs(snr_method = 'liv', plot_spectra = False, pulse_seq = 'zg30', substance = 'Lactate', tsp_fit = 'lorentzian', plot_snrs = True, plot_name = 'lac_zg30_liv.png')

_,_,_ = getSNRs(snr_method = 'brk', plot_spectra = False, pulse_seq = 'zg', substance = 'Glucose', tsp_fit = 'lorentzian', plot_snrs = True, plot_name = 'glc_zg_brk.png')
_,_,_ = getSNRs(snr_method = 'brk', plot_spectra = False, pulse_seq = 'zg', substance = 'Citrate', tsp_fit = 'lorentzian', plot_snrs = True, plot_name = 'cit_zg_brk.png')
_,_,_ = getSNRs(snr_method = 'brk', plot_spectra = False, pulse_seq = 'zg', substance = 'Lactate', tsp_fit = 'lorentzian', plot_snrs = True, plot_name = 'lac_zg_brk.png')

_,_,_ = getSNRs(snr_method = 'brk', plot_spectra = False, pulse_seq = 'zg30', substance = 'Glucose', tsp_fit = 'lorentzian', plot_snrs = True, plot_name = 'glc_zg30_brk.png')
_,_,_ = getSNRs(snr_method = 'brk', plot_spectra = False, pulse_seq = 'zg30', substance = 'Citrate', tsp_fit = 'lorentzian', plot_snrs = True, plot_name = 'cit_zg30_brk.png')
_,_,_ = getSNRs(snr_method = 'brk', plot_spectra = False, pulse_seq = 'zg30', substance = 'Lactate', tsp_fit = 'lorentzian', plot_snrs = True, plot_name = 'lac_zg30_brk.png')

'''
'''
S=LoadSpectra()

S.ReadHDF5(path = "/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/General/DataSource",
                filename = "Data.hdf5",
                substance = 'Lactate',
                concentration = "1500uM",
                field_strength = "80MHz",
                pulse_seq = "zg30")

for I in range(0, len(np.array(S.initial_ppm)[0])):
    xa = np.array(S.initial_ppm)[0][I]
    xb = np.array(S.initial_ppm)[1][I]
    xc = np.array(S.initial_ppm)[2][I]

    ya = np.array(S.initial_amplitude)[0][I]
    yb = np.array(S.initial_amplitude)[1][I]
    yc = np.array(S.initial_amplitude)[2][I]


    plt.figure()
    plt.plot(xa, ya)
    plt.plot(xb, yb)
    plt.plot(xc, yc)
    
    plt.show()
    plt.close()
'''
def getLines(filename = None, sample = None, run = None):
    path = '/Users/alexhill/Documents/UOL/Research/Companies/ViBo/Metabolomics/Data_Analysis/80MHzLactate_Oct24'
    POSITIONS = [
        'POS1',
        'POS2',
        'POS3',
        'POS4',
        'POS5',
        'POS6',
        'POS7',
        'POS8',
        'POS9',
        'POS10',
        'POS11',
        'POS12'
        ]

    POSITIONS_SHORT = [
        '01',
        '02',
        '03',
        '04',
        '05',
        '06',
        '07',
        '08',
        '09',
        '10',
        '11',
        '12'
        ]

    SAMPLES = [
        'LAC750A',
        'LAC750B',
        'LAC750C',
        'LAC1500A',
        'LAC1500B',
        'LAC1500C',
        'LAC3000A',
        'LAC3000B',
        'LAC3000C',
        'LAC6000A',
        'LAC6000B',
        'LAC6000C'
        ]

    RUNS = [
        'zg30_1',
        'zg30_2',
        'zg30_4',
        'zg30_8',
        'zg30_16',
        'zg30_32',
        'zg30_64',
        'zg30_128',
        'zg30_256',
        'zg_1',
        'zg_2',
        'zg_4',
        'zg_8',
        'zg_16',
        'zg_32',
        'zg_64',
        'zg_128',
        'zg_256'
        ]

    LAB_ACTUAL = [
        '10-zg30',
        '20-zg30',
        '30-zg30',
        '40-zg30',
        '50-zg30',
        '60-zg30',
        '70-zg30',
        '80-zg30',
        '90-zg30',
        '100-zg',
        '110-zg',
        '120-zg',
        '130-zg',
        '140-zg',
        '150-zg',
        '160-zg',
        '170-zg',
        '180-zg'
        ]

    POSITIONS = np.array(POSITIONS).astype('<U19')
    POSITIONS_SHORT = np.array(POSITIONS_SHORT).astype('<U19')
    SAMPLES = np.array(SAMPLES).astype('<U19')
    LAB_ACTUAL = np.array(LAB_ACTUAL).astype('<U19')
    RUNS = np.array(RUNS).astype('<U19')

    if filename is None:
        prefix1 = POSITIONS[SAMPLES == sample]
        prefix2 = POSITIONS_SHORT[SAMPLES == sample]
        lact    = LAB_ACTUAL[RUNS == run]
        filename = prefix1[0]+'/2024-10-28-AH-LAC-'+prefix2[0]+'-'+lact[0]+'.txt'
    S = LoadSpectra()
    S.ReadTextFile(path = path, filename = filename)
    x = S.initial_ppm
    y = S.initial_amplitude    
    return x, y

def getData(filename = None, sample = None, run = None):
    path = '/Users/alexhill/Documents/UOL/Research/Companies/ViBo/Metabolomics/Data_Analysis/80MHzLactate_Oct24'
    POSITIONS = [
        'POS1',
        'POS2',
        'POS3',
        'POS4',
        'POS5',
        'POS6',
        'POS7',
        'POS8',
        'POS9',
        'POS10',
        'POS11',
        'POS12'
        ]

    POSITIONS_SHORT = [
        '01',
        '02',
        '03',
        '04',
        '05',
        '06',
        '07',
        '08',
        '09',
        '10',
        '11',
        '12'
        ]

    SAMPLES = [
        'LAC750A',
        'LAC750B',
        'LAC750C',
        'LAC1500A',
        'LAC1500B',
        'LAC1500C',
        'LAC3000A',
        'LAC3000B',
        'LAC3000C',
        'LAC6000A',
        'LAC6000B',
        'LAC6000C'
        ]

    RUNS = [
        'zg30_1',
        'zg30_2',
        'zg30_4',
        'zg30_8',
        'zg30_16',
        'zg30_32',
        'zg30_64',
        'zg30_128',
        'zg30_256',
        'zg_1',
        'zg_2',
        'zg_4',
        'zg_8',
        'zg_16',
        'zg_32',
        'zg_64',
        'zg_128',
        'zg_256'
        ]

    LAB_ACTUAL = [
        '10-zg30',
        '20-zg30',
        '30-zg30',
        '40-zg30',
        '50-zg30',
        '60-zg30',
        '70-zg30',
        '80-zg30',
        '90-zg30',
        '100-zg',
        '110-zg',
        '120-zg',
        '130-zg',
        '140-zg',
        '150-zg',
        '160-zg',
        '170-zg',
        '180-zg'
        ]

    POSITIONS = np.array(POSITIONS).astype('<U19')
    POSITIONS_SHORT = np.array(POSITIONS_SHORT).astype('<U19')
    SAMPLES = np.array(SAMPLES).astype('<U19')
    LAB_ACTUAL = np.array(LAB_ACTUAL).astype('<U19')
    RUNS = np.array(RUNS).astype('<U19')

    if filename is None:
        prefix1 = POSITIONS[SAMPLES == sample]
        prefix2 = POSITIONS_SHORT[SAMPLES == sample]
        lact    = LAB_ACTUAL[RUNS == run]
        filename = prefix1[0]+'/2024-10-28-AH-LAC-'+prefix2[0]+'-'+lact[0]+'.txt'
    S = LoadSpectra()
    S.ReadTextFile(path = path, filename = filename)
    A = AnalyseSpectra()
    A.InputData(x = S.initial_ppm, y = S.initial_amplitude)
    A.SignalToNoise(signal_bounds = lactic_acid_bounds, snr_choice = 'brk')
    sig = A.spectra_signal
    noise = A.spectra_noise
    snr = A.snr
    return sig, noise, snr

def sampleSNR(sample = 'LAC6000', exp = 'zg30'):
    sampA = sample+'A'
    sampB = sample+'B'
    sampC = sample+'C'
    A = []
    B = []
    C = []
    n_s = ['1','2','4','8','16','32','64','128','256']
    for N in n_s:
        run = exp+'_'+N
        print(run)
        sig, noise, snr = getData(sample = sampA, run = run)
        A.append([sig, noise, snr])

        sig, noise, snr = getData(sample = sampB, run = run)
        B.append([sig, noise, snr])

        sig, noise, snr = getData(sample = sampC, run = run)
        C.append([sig, noise, snr])
    return np.array(A), np.array(B), np.array(C)

'''
A6000, B6000, C6000 = sampleSNR(sample = 'LAC6000', exp = 'zg')
A3000, B3000, C3000 = sampleSNR(sample = 'LAC3000', exp = 'zg')
A1500, B1500, C1500 = sampleSNR(sample = 'LAC1500', exp = 'zg')
A750,  B750,  C750  = sampleSNR(sample = 'LAC750', exp = 'zg')


n_s = [1,2,4,8,16,32,64,128,256]


# Calculate means and standard errors for each group
def get_stats(A, B, C, data = 'sig'):
    idx = 0
    if data == 'noise':
        idx = 1
    if data == 'snr':
        idx = 2
    # Stack the arrays to get all measurements for each n_scan
    all_measurements = np.vstack((A[:,idx], B[:,idx], C[:,idx]))
    # Calculate mean and standard error
    means = np.mean(all_measurements, axis=0)
    std_err = np.std(all_measurements, axis=0) / np.sqrt(3)  # standard error for 3 measurements
    return means, std_err



# Define constants
data_types = ['sig', 'noise', 'snr']
scan_numbers = [750, 1500, 3000, 6000]
colors = ['orange', 'green', 'blue', 'red']
ylabels = ['Signal', 'Noise', 'SNR']

# Create arrays to store results
means = {dtype: np.zeros((4, len(n_s))) for dtype in data_types}
errors = {dtype: np.zeros((4, len(n_s))) for dtype in data_types}

# Calculate statistics for all combinations
sample_sets = [(A750, B750, C750), (A1500, B1500, C1500),
               (A3000, B3000, C3000), (A6000, B6000, C6000)]

for dtype in data_types:
    for i, (A, B, C) in enumerate(sample_sets):
        means[dtype][i], errors[dtype][i] = get_stats(A, B, C, data=dtype)

# Create and populate subplots
fig, axs = plt.subplots(1, 3, figsize=[18, 6])

for i, (dtype, ylabel) in enumerate(zip(data_types, ylabels)):
    for j, (scan_num, color) in enumerate(zip(scan_numbers, colors)):
        axs[i].errorbar(n_s, means[dtype][j], yerr=errors[dtype][j],
                       fmt='o', color=color, label=str(scan_num), capsize=5)
    
    axs[i].set_xlabel('N scans')
    axs[i].set_ylabel(ylabel)
    axs[i].set_xscale('log')
    axs[i].set_yscale('log')
    axs[i].legend()

plt.savefig('zg_1d.png')

# Stack SNR means if needed
snr_stack = np.vstack([means['snr'][i] for i in range(4)])

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm  # Add this import

import numpy as np
import matplotlib.pyplot as plt

# Create figure
plt.figure(figsize=(8, 6))

# Basic imshow plot
im = plt.imshow(np.log10(snr_stack), aspect='auto', origin='lower', vmin = 0, vmax = 2)

# Calculate centers of bins for y-axis (4 bins)
y_positions = np.arange(4)  # [0, 1, 2, 3]
y_labels = ['750', '1500', '3000', '6000']  # Custom labels
plt.yticks(y_positions, y_labels)

# Calculate centers of bins for x-axis (assuming 9 bins for the NScans values)
x_positions = np.arange(9)  # [0, 1, 2, 3, 4, 5, 6, 7, 8]
x_labels = ['1', '2', '4', '8', '16', '32', '64', '128', '256']  # Custom labels
plt.xticks(x_positions, x_labels)

# Add colorbar and labels
cbar = plt.colorbar(im)
cbar.set_label('log₁₀(SNR)')
plt.ylabel('Concentration (mM)')
plt.xlabel('NScans')

plt.tight_layout()
plt.savefig('zg_2d.png')



fnms = np.array([
'POS10/2024-10-28-AH-LAC-10-100-zg.txt',
'POS10/2024-10-28-AH-LAC-10-110-zg.txt',
'POS10/2024-10-28-AH-LAC-10-120-zg.txt',
'POS10/2024-10-28-AH-LAC-10-130-zg.txt',
'POS10/2024-10-28-AH-LAC-10-140-zg.txt',
'POS10/2024-10-28-AH-LAC-10-150-zg.txt',
'POS10/2024-10-28-AH-LAC-10-160-zg.txt',
'POS10/2024-10-28-AH-LAC-10-170-zg.txt',
'POS10/2024-10-28-AH-LAC-10-180-zg.txt'
])

plt.figure()

for i in range(0, len(fnms)):
    x, y = getLines(filename = fnms[i])
    plt.plot(x, y, label = fnms[i])

plt.legend()

plt.show()
'''
'''
# Calculate statistics for each group
means_6000_sig, err_6000_sig = get_stats(A6000, B6000, C6000, data = 'sig')
means_3000_sig, err_3000_sig = get_stats(A3000, B3000, C3000, data = 'sig')
means_1500_sig, err_1500_sig = get_stats(A1500, B1500, C1500, data = 'sig')
means_750_sig,  err_750_sig = get_stats(A750, B750, C750, data = 'sig')

means_6000_noise, err_6000_noise = get_stats(A6000, B6000, C6000, data = 'noise')
means_3000_noise, err_3000_noise = get_stats(A3000, B3000, C3000, data = 'noise')
means_1500_noise, err_1500_noise = get_stats(A1500, B1500, C1500, data = 'noise')
means_750_noise,  err_750_noise = get_stats(A750, B750, C750, data = 'noise')

means_6000_snr, err_6000_snr = get_stats(A6000, B6000, C6000, data = 'snr')
means_3000_snr, err_3000_snr = get_stats(A3000, B3000, C3000, data = 'snr')
means_1500_snr, err_1500_snr = get_stats(A1500, B1500, C1500, data = 'snr')
means_750_snr,  err_750_snr = get_stats(A750, B750, C750, data = 'snr')


# Create the scatter plot with error bars
fig, axs = plt.subplots(1,3, figsize = [18,6])
axs[0].errorbar(n_s, means_6000_sig, yerr=err_6000_sig, fmt='o', color='orange', label='6000', capsize=5)
axs[0].errorbar(n_s, means_3000_sig, yerr=err_3000_sig, fmt='o', color='green', label='3000', capsize=5)
axs[0].errorbar(n_s, means_1500_sig, yerr=err_1500_sig, fmt='o', color='blue', label='1500', capsize=5)
axs[0].errorbar(n_s, means_750_sig, yerr=err_750_sig, fmt='o', color='red', label='750', capsize=5)
axs[0].set_xlabel('N scans')
axs[0].set_ylabel('Signal')
axs[0].legend()

axs[1].errorbar(n_s, means_6000_noise, yerr=err_6000_noise, fmt='o', color='orange', label='6000', capsize=5)
axs[1].errorbar(n_s, means_3000_noise, yerr=err_3000_noise, fmt='o', color='green', label='3000', capsize=5)
axs[1].errorbar(n_s, means_1500_noise, yerr=err_1500_noise, fmt='o', color='blue', label='1500', capsize=5)
axs[1].errorbar(n_s, means_750_noise, yerr=err_750_noise, fmt='o', color='red', label='750', capsize=5)
axs[1].set_xlabel('N scans')
axs[1].set_ylabel('Noise')
axs[1].legend()

axs[2].errorbar(n_s, means_6000_snr, yerr=err_6000_snr, fmt='o', color='orange', label='6000', capsize=5)
axs[2].errorbar(n_s, means_3000_snr, yerr=err_3000_snr, fmt='o', color='green', label='3000', capsize=5)
axs[2].errorbar(n_s, means_1500_snr, yerr=err_1500_snr, fmt='o', color='blue', label='1500', capsize=5)
axs[2].errorbar(n_s, means_750_snr, yerr=err_750_snr, fmt='o', color='red', label='750', capsize=5)
axs[2].set_xlabel('N scans')
axs[2].set_ylabel('SNR')
axs[2].legend()

axs[0].set_yscale('log')
axs[0].set_xscale('log')

axs[1].set_yscale('log')
axs[1].set_xscale('log')

axs[2].set_yscale('log')
axs[2].set_xscale('log')

plt.show()

np.vstack((means_6000_snr,
            means_3000_snr,
            means_1500_snr,
            means_750_snr))
'''
'''

plt.close()
plt.figure()
plt.plot(n_s, A6000[:,0], label = 'sampA', color = 'orange')
plt.plot(n_s, B6000[:,0], label = 'sampB', color = 'orange')
plt.plot(n_s, C6000[:,0], label = 'sampC', color = 'orange')

plt.plot(n_s, A3000[:,0], label = 'sampA', color = 'green')
plt.plot(n_s, B3000[:,0], label = 'sampB', color = 'green')
plt.plot(n_s, C3000[:,0], label = 'sampC', color = 'green')


plt.plot(n_s, A1500[:,0], label = 'sampA', color = 'blue')
plt.plot(n_s, B1500[:,0], label = 'sampB', color = 'blue')
plt.plot(n_s, C1500[:,0], label = 'sampC', color = 'blue')


plt.plot(n_s, A750[:,0], label = 'sampA', color = 'red')
plt.plot(n_s, B750[:,0], label = 'sampB', color = 'red')
plt.plot(n_s, C750[:,0], label = 'sampC', color = 'red')


plt.xlabel('N scans')
plt.ylabel('Signal')
#plt.legend()
plt.show()
'''
'''
SIGS_0s, NOISES_0s, SNRS_0s = [], [], []
SIGS_1s, NOISES_1s, SNRS_1s = [], [], []

for j in range(1, 13):
    exp = int(j * 10)
    sig, noise, snr = getData(exp)
    SIGS_0s.append(sig)
    NOISES_0s.append(noise)
    SNRS_0s.append(snr)

    exp = int(j * 10 + 1)
    sig, noise, snr = getData(exp)
    SIGS_1s.append(sig)
    NOISES_1s.append(noise)
    SNRS_1s.append(snr)

# Calculate stats for each group
x_vals = [6, 3, 1.5, 0.75]

# Split into groups of 3 and calculate statistics
groups = [SIGS_0s[i:i+3] for i in range(0, len(SIGS_0s), 3)]
means = np.array([np.mean(group) for group in groups])
mins = np.array([np.min(group) for group in groups])
maxs = np.array([np.max(group) for group in groups])

# Calculate error bar lengths
yerr_lower = means - mins
yerr_upper = maxs - means
yerr = np.array([yerr_lower, yerr_upper])

# Linear fit
slope, intercept, r_value, p_value, std_err = stats.linregress(x_vals, means)
fit_line = slope * np.array(x_vals) + intercept

# Create the plot
plt.figure(figsize=(8, 6))
plt.errorbar(x_vals, means, yerr=yerr, fmt='o', capsize=5, label='Data')
plt.plot(x_vals, fit_line, 'r-', label=f'Linear fit (R² = {r_value**2:.3f})')

plt.xlabel('Concentration (mM)')
plt.ylabel('Signal')
plt.legend()
plt.grid(True)

# Optional: set y-axis to start at 0
plt.ylim(bottom=0)
plt.xlim([0, 7])

# Display the plot
plt.show()

# Print fit parameters
print(f"Slope: {slope:.2f}")
print(f"Intercept: {intercept:.2f}")
print(f"R²: {r_value**2:.3f}")

'''
