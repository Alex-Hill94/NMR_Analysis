import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '/Users/alexhill/Documents/GitHub/NMR_Analysis/Fourier90Runs')
from NMRClassesH5 import *
import matplotlib as mpl
import os 
import csv
import pandas as pd 
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import h5py as h5
from scipy import stats
import statsmodels.api as sm


path = '/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/General'

lactic_acid_bounds = [1.17, 1.50]
glucose_bounds = [3.19, 3.98]
citric_acid_bounds = [2.37, 2.82]
'''
lactic_acid_bounds = [1.31, 1.35]
glucose_bounds = [3.19, 3.97]
citric_acid_bounds = [2.52, 2.67]

lactic_acid_bounds = [1.21, 1.45]
glucose_bounds = [3.09, 4.07]
citric_acid_bounds = [2.42, 2.77]
'''

def getSNR(substance = 'Lactate'):

    if substance == 'Lactate':
        concentrations = ["750uM", "1500uM", "3000uM", "6000uM"]
        bounds = lactic_acid_bounds   
    
    elif substance == 'Citrate':
        concentrations = ["50uM", "100uM", "200uM", "400uM"]
        bounds = citric_acid_bounds

    elif substance == 'Glucose':
        concentrations = ["1250uM", "2500uM", "5000uM", "10000uM"]
        bounds = glucose_bounds
    

    SNRS = []
    NOISE = []
    SIGS = []
    Xs = []
    Ys = []
    for conc in concentrations:
        S = LoadSpectra()
        S.ReadHDF5(path = path,
                filename = 'Data.hdf5',
                substance = substance,
                concentration = conc,
                field_strength = '700MHz',
                pulse_seq = 'noesy')
        for i in [0, 1, 2]:
            x, y = S.initial_ppm[i], S.initial_amplitude[i]
            A = AnalyseSpectra()
            A.InputData(x = x, y = y)
            A.FitTSP()
            A.ScaleSpectra(x_shift = (-1.)*A.tsp_centre)
            #A.QuickPlot()
            A.SignalToNoise(signal_bounds = bounds, snr_choice = 'liv')
            sig = A.spectra_signal
            noise = A.spectra_noise
            snr = A.snr
            SNRS.append(snr)
            NOISE.append(noise)
            SIGS.append(sig)
            print(sig, conc, i)
            Ys.append(A.processed_amplitude)
            Xs.append(A.processed_ppm)
    return SIGS, NOISE, SNRS, Xs, Ys, bounds

for substance in ['Lactate', 'Glucose', 'Citrate']:

    if substance == 'Lactate':
        SIGS_0s, NOISES_0s, SNRS_0s, Xs, Ys, bounds = getSNR()
        x_vals = [750, 1500, 3000, 6000]
        xlims = [0, 7000]
    elif substance == 'Glucose':
        SIGS_0s, NOISES_0s, SNRS_0s, Xs, Ys, bounds = getSNR(substance = 'Glucose')
        x_vals = [1250, 2500, 5000, 10000]
        xlims = [0, 11000]
    elif substance == 'Citrate':
        SIGS_0s, NOISES_0s, SNRS_0s, Xs, Ys, bounds = getSNR(substance = 'Citrate')
        x_vals = [50,  100,  200, 400]
        xlims = [0, 500]

    plt.figure()
    for i in range(0, len(Xs)):
        plt.plot(Xs[i], Ys[i])
    plt.xlim(bounds)
    plt.savefig(substance+'_int_area.png')
    plt.close()

    scale = (x_vals[-1]/np.mean(SIGS_0s[-3:]))
    SIGS_0s = np.array(SIGS_0s)*scale

    x_vals = [val for val in x_vals for _ in range(3)]

    # Input data
    x = x_vals
    y = SIGS_0s

    # Calculate unique x values, means, and SEM
    unique_x = np.unique(x)
    means = np.array([np.mean(y[x == val]) for val in unique_x])
    std_devs = np.array([np.std(y[x == val], ddof=1) for val in unique_x])  # Sample standard deviation
    sems = std_devs / np.sqrt(3)  # Standard Error of the Mean

    # Weighted least squares (WLS) regression
    weights = 1 / (sems**2)  # Inverse variance weighting
    X = sm.add_constant(unique_x)  # Adds intercept term
    model = sm.WLS(means, X, weights=weights)
    results = model.fit()

    # Extract parameters
    slope = results.params[1]
    intercept = results.params[0]
    R2 = results.rsquared
    # Generate the fitted line
    x_fit = np.linspace(min(unique_x), max(unique_x), 100)
    y_fit = slope * x_fit + intercept

    # Plot the data and the fitted line
    plt.figure(figsize=(8, 6))
    plt.errorbar(unique_x, means, yerr=sems, fmt='o', label='Mean ± SEM', capsize=5, color='blue')
    plt.plot(x_fit, y_fit, label=f'Fitted Line: y = {slope:.2f}x + {intercept:.2f}\n R2 = {R2:.5f}', color='red', linewidth=2)
    plt.scatter(x, y, alpha=0.5, label='Original Data', color='gray')

    # Customize the plot
    plt.title(substance)
    plt.xlabel('Concentration [uM]')
    plt.ylabel('Intensity [Scaled Units]')
    plt.legend()
    plt.grid(True)
    plt.savefig(substance+'_80MHz_Bounds.png')
    plt.close()
    # Print the regression summary
    print(results.summary())


'''

S = LoadSpectra()
S.ReadHDF5(path = path,
           filename = 'Data.hdf5',
           substance = 'Lactate',
           concentration = '6000uM',
           field_strength = '700MHz',
           pulse_seq = 'noesy')

plt.figure()
plt.plot(S.initial_ppm[0], S.initial_amplitude[0], color = 'tab:blue')
plt.plot(S.initial_ppm[1], S.initial_amplitude[1], color = 'tab:blue')
plt.plot(S.initial_ppm[2], S.initial_amplitude[2], color = 'tab:blue')

S.ReadHDF5(path = path,
           filename = 'Data.hdf5',
           substance = 'Lactate',
           concentration = '3000uM',
           field_strength = '700MHz',
           pulse_seq = 'noesy')
plt.plot(S.initial_ppm[0], S.initial_amplitude[0], color = 'tab:red')
plt.plot(S.initial_ppm[1], S.initial_amplitude[1], color = 'tab:red')
plt.plot(S.initial_ppm[2], S.initial_amplitude[2], color = 'tab:red')

S.ReadHDF5(path = path,
           filename = 'Data.hdf5',
           substance = 'Lactate',
           concentration = '1500uM',
           field_strength = '700MHz',
           pulse_seq = 'noesy')

plt.plot(S.initial_ppm[0], S.initial_amplitude[0], color = 'tab:green')
plt.plot(S.initial_ppm[1], S.initial_amplitude[1], color = 'tab:green')
plt.plot(S.initial_ppm[2], S.initial_amplitude[2], color = 'tab:green')


S.ReadHDF5(path = path,
           filename = 'Data.hdf5',
           substance = 'Lactate',
           concentration = '750uM',
           field_strength = '700MHz',
           pulse_seq = 'noesy')

plt.plot(S.initial_ppm[0], S.initial_amplitude[0], color = 'tab:grey')
plt.plot(S.initial_ppm[1], S.initial_amplitude[1], color = 'tab:grey')
plt.plot(S.initial_ppm[2], S.initial_amplitude[2], color = 'tab:grey')


plt.show()
'''

'''

labs = [6, 3, 1.5, 0.75]
labs_abc = ['A', 'B', 'C']
combined = [f"{conc}mM {letter}" for conc in labs for letter in labs_abc for _ in range(2)]

def getData(exp):
    filename = 'AHill_241025_7_exp%s.txt' % exp
    S = LoadSpectra()
    S.ReadTextFile(path = path, filename = filename)
    A = AnalyseSpectra()
    A.InputData(x = S.initial_ppm, y = S.initial_amplitude)
    A.SignalToNoise(signal_bounds = lactic_acid_bounds, snr_choice = 'agi')
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


'''