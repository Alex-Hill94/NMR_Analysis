import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '/Users/alexhill/Documents/GitHub/NMR_Analysis/Fourier90Runs')
from NMRClasses_Lac import *
import matplotlib as mpl
import os 
import csv
import pandas as pd 
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import h5py as h5
from scipy import stats


path = '/Users/alexhill/Documents/UOL/Research/Companies/ViBo/Metabolomics/Data_Analysis/700MHzLactate'
lactic_acid_bounds = [1.17, 1.50]

suffixes = [10, 11, 20, 21, 30, 31, 40, 41, 50, 51, 60, 61, 70, 71, 80, 81, 90, 91, 100, 101, 110, 111, 120, 121]

labs = [6, 3, 1.5, 0.75]
labs_abc = ['A', 'B', 'C']
combined = [f"{conc}mM {letter}" for conc in labs for letter in labs_abc for _ in range(2)]

def getData(exp):
    filename = 'AHill_241025_7_exp%s.txt' % exp
    S = LoadSpectra()
    S.ReadTextFile(path = path, filename = filename)
    A = AnalyseSpectra()
    A.InputData(x = S.initial_ppm, y = S.initial_amplitude)
    A.SignalToNoise(signal_bounds = lactic_acid_bounds, snr_choice = 'brk')
    sig = A.spectra_signal
    noise = A.spectra_noise
    snr = A.snr
    return sig, noise, snr

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

