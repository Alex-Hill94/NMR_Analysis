import sys
import os 
import matplotlib as mpl
import csv
import pandas as pd 
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import matplotlib.patheffects as path_effects
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import h5py as h5
from scipy import stats
from scipy.stats import chi2
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '/Users/alexhill/Documents/GitHub/NMR_Analysis/Fourier90Runs')
from paul_col import *
from paul_tol_colours import *
parent_dir = os.path.abspath(os.path.join(os.getcwd(), '..'))
sys.path.append(parent_dir)
from NMRClassesH5 import *
sys.path.insert(1, '/Users/alexhill/Documents/GitHub/NMR_Analysis/Simulations/Fittings')
from FncLib import *
from scipy.optimize import curve_fit

lactic_acid_bounds = [1.17, 1.50]
glucose_bounds = [3.19, 3.98]
citric_acid_bounds = [2.37, 2.82]

lactic_acid_concs = [750, 1500, 3000, 6000]
glucose_concs = [1250, 2500, 5000, 10000]
citric_acid_concs = [50,  100,  200, 400]

scans = [1, 2, 4, 8, 16, 32, 64, 128, 256]
sample = 'D24'
snr_choice = 'liv'
times = []

def getSNRs(substance = 'Lactate', 
            pulse_seq = 'zg', 
            plot_spectra = True, 
            plot_snrs = True, 
            snr_method = 'liv', 
            tsp_fit = 'lorentzian', 
            plot_name = None):

    if substance == 'Lactate':
        concs = lactic_acid_concs
        bounds = lactic_acid_bounds
    if substance == 'Glucose':
        concs = glucose_concs
        bounds = glucose_bounds
    if substance == 'Citrate':
        concs = citric_acid_concs
        bounds = citric_acid_bounds

    conc_labs = [str(val)+"uM" for val in concs]

    S = LoadSpectra()
    if plot_spectra:
        fig, axs = plt.subplots(2,1, figsize = [12,6])
        for conc in conc_labs:
            S.ReadHDF5(path = "/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/General/DataSource",
                            filename = "Data.hdf5",
                            substance = substance,
                            concentration = conc,
                            field_strength = "80MHz",
                            pulse_seq = pulse_seq)

            x, y = S.initial_ppm[1][-1], S.initial_amplitude[1][-1]

            A = AnalyseSpectra()
            A.InputData(x=x, y=y)
            if tsp_fit is not None:
                A.FitTSP(plot_fit = False, fit_function = tsp_fit)
                A.ScaleSpectra(x_shift = (-1) * A.tsp_centre)
            axs[0].plot(x, y, label = conc)
            axs[1].plot(A.processed_ppm, A.processed_amplitude)
            axs[0].axvline(bounds[0], color='grey', linewidth=1.0, alpha=0.7, linestyle = '--')    
            axs[0].axvline(bounds[1], color='grey', linewidth=1.0, alpha=0.7, linestyle = '--')    
            axs[1].axvline(bounds[0], color='grey', linewidth=1.0, alpha=0.7, linestyle = '--')    
            axs[1].axvline(bounds[1], color='grey', linewidth=1.0, alpha=0.7, linestyle = '--')    
        axs[0].legend()
        plt.show()

    SIGS, NOISES, SNRs = [],[],[]

    for conc in conc_labs:
        S.ReadHDF5(path = "/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/General/DataSource",
                        filename = "Data.hdf5",
                        substance = substance,
                        concentration = conc,
                        field_strength = "80MHz",
                        pulse_seq = pulse_seq)
        for j in [0,1,2]:
            for k in range(0, len(S.initial_ppm[j])):
                #print(conc,j,k)
                x, y = S.initial_ppm[j][k], S.initial_amplitude[j][k]
                A = AnalyseSpectra()
                A.InputData(x=x, y=y)
                if tsp_fit is not None:
                    A.FitTSP(plot_fit = False, fit_function = tsp_fit)
                    A.ScaleSpectra(x_shift = (-1) * A.tsp_centre)
                A.SignalToNoise(signal_bounds = bounds, snr_choice = snr_method)
                SIGS.append(A.spectra_signal)
                NOISES.append(A.spectra_noise)
                SNRs.append(A.snr)

    SIGS = np.array(SIGS).reshape(4,3,9)
    NOISES = np.array(NOISES).reshape(4,3,9)
    SNRs = np.array(SNRs).reshape(4,3,9)
    
    if plot_snrs:
        fig, axs = plt.subplots(1,3, figsize = [12, 4])
        axs[0].scatter(S.nscans, np.mean(SIGS[0], axis = 0), color = 'orange')
        axs[0].scatter(S.nscans, np.mean(SIGS[1], axis = 0), color = 'green')
        axs[0].scatter(S.nscans, np.mean(SIGS[2], axis = 0), color = 'blue')
        axs[0].scatter(S.nscans, np.mean(SIGS[3], axis = 0), color = 'red')

        axs[1].scatter(S.nscans, np.mean(NOISES[0], axis = 0), color = 'orange')
        axs[1].scatter(S.nscans, np.mean(NOISES[1], axis = 0), color = 'green')
        axs[1].scatter(S.nscans, np.mean(NOISES[2], axis = 0), color = 'blue')
        axs[1].scatter(S.nscans, np.mean(NOISES[3], axis = 0), color = 'red')

        axs[2].scatter(S.nscans, np.mean(SNRs[0], axis = 0), color = 'orange')
        axs[2].scatter(S.nscans, np.mean(SNRs[1], axis = 0), color = 'green')
        axs[2].scatter(S.nscans, np.mean(SNRs[2], axis = 0), color = 'blue')
        axs[2].scatter(S.nscans, np.mean(SNRs[3], axis = 0), color = 'red')


        axs[0].set_xscale('log')
        axs[1].set_xscale('log')
        axs[2].set_xscale('log')

        axs[0].set_yscale('log')
        axs[1].set_yscale('log')
        axs[2].set_yscale('log')

        if plot_name is not None:
            plt.savefig(plot_name)
        else:    
            plt.show()  
        plt.close()
    return SIGS, NOISES, SNRs



sigs_glc, noises_glc, snrs_glc = getSNRs(substance = 'Glucose', pulse_seq = 'zg', plot_spectra = False, plot_snrs = False)
sigs_lac, noises_lac, snrs_lac = getSNRs(substance = 'Lactate', pulse_seq = 'zg', plot_spectra = False, plot_snrs = False)
sigs_cit, noises_cit, snrs_cit = getSNRs(substance = 'Citrate', pulse_seq = 'zg', plot_spectra = False, plot_snrs = False)

SIG_glc = np.mean(sigs_glc, axis = 1)
SIG_lac = np.mean(sigs_lac, axis = 1)
SIG_cit = np.mean(sigs_cit, axis = 1)



path = 'new_data.hdf5'
u = h5.File(path, 'r')
glc = u['glcnew'][()]
lac = u['lac'][()]
cit = u['cit'][()]
u.close()




fig, axs = plt.subplots(2, 3, figsize = [10, 8])

axs[0,0].scatter(scans, glc[0], label = '1250uM')
axs[0,0].scatter(scans, glc[1], label = '2500uM')
axs[0,0].scatter(scans, glc[2], label = '5000uM')
axs[0,0].scatter(scans, glc[3], label = '10000uM')

axs[0,1].scatter(scans, lac[0], label = '750uM')
axs[0,1].scatter(scans, lac[1], label = '1500uM')
axs[0,1].scatter(scans, lac[2], label = '3000uM')
axs[0,1].scatter(scans, lac[3], label = '6000uM')


axs[0,2].scatter(scans, cit[0], label = '50uM')
axs[0,2].scatter(scans, cit[1], label = '100uM')
axs[0,2].scatter(scans, cit[2], label = '200uM')
axs[0,2].scatter(scans, cit[3], label = '400uM')

axs[0,0].set_xlabel('Nscans')
axs[0,1].set_xlabel('Nscans')
axs[0,2].set_xlabel('Nscans')

axs[0,0].set_xscale('log')
axs[0,1].set_xscale('log')
axs[0,2].set_xscale('log')
axs[0,0].set_yscale('log')
axs[0,1].set_yscale('log')
axs[0,2].set_yscale('log')

axs[0,0].legend()
axs[0,1].legend()
axs[0,2].legend()



axs[1,0].scatter(scans, SIG_glc[0], label = '1250uM')
axs[1,0].scatter(scans, SIG_glc[1], label = '2500uM')
axs[1,0].scatter(scans, SIG_glc[2], label = '5000uM')
axs[1,0].scatter(scans, SIG_glc[3], label = '10000uM')

axs[1,1].scatter(scans, SIG_lac[0], label = '750uM')
axs[1,1].scatter(scans, SIG_lac[1], label = '1500uM')
axs[1,1].scatter(scans, SIG_lac[2], label = '3000uM')
axs[1,1].scatter(scans, SIG_lac[3], label = '6000uM')


axs[1,2].scatter(scans, SIG_cit[0], label = '50uM')
axs[1,2].scatter(scans, SIG_cit[1], label = '100uM')
axs[1,2].scatter(scans, SIG_cit[2], label = '200uM')
axs[1,2].scatter(scans, SIG_cit[3], label = '400uM')

axs[1,0].set_xlabel('Nscans')
axs[1,1].set_xlabel('Nscans')
axs[1,2].set_xlabel('Nscans')

axs[1,0].set_xscale('log')
axs[1,1].set_xscale('log')
axs[1,2].set_xscale('log')
axs[1,0].set_yscale('log')
axs[1,1].set_yscale('log')
axs[1,2].set_yscale('log')

plt.show()