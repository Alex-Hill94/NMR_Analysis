import sys
import os 
import matplotlib as mpl
import csv
import pandas as pd 
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import h5py as h5
from scipy import stats
from scipy.stats import chi2
from paul_col import *
from paul_tol_colours import *
# caution: path[0] is reserved for script path (or '' in REPL)

sys.path.insert(1, '/Users/alexhill/Documents/GitHub/NMR_Analysis/Fourier90Runs')
parent_dir = os.path.abspath(os.path.join(os.getcwd(), '..'))
sys.path.append(parent_dir)
from NMRClassesH5 import *

lactic_acid_bounds = [1.17, 1.50]
glucose_bounds = [3.19, 3.98]
citric_acid_bounds = [2.37, 2.82]

lactic_acid_concs = [750, 1500, 3000, 6000]
glucose_concs = [1250, 2500, 5000, 10000]
citric_acid_concs = [50,  100,  200, 400]

def goodness_of_fit_p_value(x, y, yerr, A, slope, scale_factor=1.0):
    """
    Perform a chi-squared goodness-of-fit test to assess how well the data fits 
    a power-law model of the form y = A * x^slope * scale_factor.

    Parameters:
        x (array-like): Independent variable (number of scans).
        y (array-like): Observed dependent variable (measured SNR).
        yerr (array-like): Standard error on the mean of y.
        A (float): Scaling coefficient from the best fit.
        slope (float): Slope from the best fit in log-log space.
        scale_factor (float): Factor to scale A (e.g. 1, 0.5, 0.25, 0.125).

    Returns:
        float: p-value indicating the probability that the data fits the model.
    """
    # Expected values based on the scaled model
    expected_y = A * scale_factor * x**slope

    # Chi-squared test statistic
    chi_sq = np.sum(((y - expected_y) / yerr) ** 2)

    # Degrees of freedom (number of data points - 1)
    dof = len(x) - 1

    # Compute the p-value from chi-squared distribution
    p_value = 1 - chi2.cdf(chi_sq, dof)

    return p_value

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

            x, y = S.initial_ppm[0][-1], S.initial_amplitude[0][-1]

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



def snr_fig(substance = 'Lactate', 
            pulse_seq = 'zg', 
            plot_spectra = False, 
            plot_snrs = False, 
            snr_method = 'liv', 
            tsp_fit = 'lorentzian', 
            plot_name = None):

    def _getylims(Xdata, Ydata, bounds):
        ys = Ydata[(Xdata <= bounds[1]) * (Xdata >= bounds[0])]
        ymax = np.nanmax(ys)
        ymin = np.nanmin(ys)
        ylims = [ymin, ymax]
        return ylims

    SIGS, NOISES, SNRs = getSNRs(snr_method = snr_method, plot_spectra = plot_spectra, pulse_seq = pulse_seq, substance = substance, tsp_fit = tsp_fit, plot_snrs = plot_snrs, plot_name = plot_name)

    nscans_long = np.arange(0, 9, 1)
    nscans_long_lab = np.array([2**val for val in nscans_long])


    MEANS = []
    SEMS  = []
    for i in range(0, len(SNRs)):
        y = SNRs[i]
        means = np.mean(y, axis = 0)
        std_devs = np.std(y, axis = 0)
        sems = std_devs / np.sqrt(3)  
        MEANS.append(means)
        SEMS.append(sems)

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
    concs_oth_lab = [(val/1000.) for val in concs]
    nscans_short = [4, 6, 8]
    nscans_long = np.arange(0, 9, 1)

    S = LoadSpectra()

    concs = [conc_labs[0], conc_labs[-1]]

    X = []
    Y = []
    for conc in concs:

        S.ReadHDF5(path = "/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/General/DataSource",
                        filename = "Data.hdf5",
                        substance = substance,
                        concentration = conc,
                        field_strength = "80MHz",
                        pulse_seq = pulse_seq)



        for n in nscans_short:
            x, y = S.initial_ppm[0][n], S.initial_amplitude[0][n]
            A = AnalyseSpectra()
            A.InputData(x=x, y=y)
            A.FitTSP(plot_fit = False, fit_function = tsp_fit)
            A.ScaleSpectra(x_shift = (-1) * A.tsp_centre, y_scaling = 1./A.tsp_integral)
            X.append(A.processed_ppm)
            Y.append(A.processed_amplitude)


    import matplotlib.pyplot as plt
    from matplotlib import gridspec

    cols = ['violet', 'green', 'red']
    matplotlib.rc('text', usetex=True)

    # Set the default font to Times New Roman
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.size'] = 14

    # Create figure
    fig = plt.figure(figsize=(10, 5))

    # Create GridSpec
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 2], height_ratios=[1, 1])

    # Create subplots
    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[1, 0])
    ax3 = plt.subplot(gs[:, 1])

    ax1.plot(X[0], Y[0], color = cols[0], label = '$n_{\mathrm{s}} =  2^{4}$')#, handlelength = 1)
    ax1.plot(X[1], Y[1], color = cols[1], label = '$n_{\mathrm{s}} = 2^{6}$')#, handlelength = 1)
    ax1.plot(X[2], Y[2], color = cols[2], label = '$n_{\mathrm{s}} = 2^{8}$')#, handlelength = 1)

    a = _getylims(X[0], Y[0], bounds)
    ax1.set_ylim([a[0], a[1] * 1.3])

    ax2.plot(X[3], Y[3], color = cols[0], label = '$n_{\mathrm{s}} =  2^{4}$')
    ax2.plot(X[4], Y[4], color = cols[1], label = '$n_{\mathrm{s}} = 2^{6}$')
    ax2.plot(X[5], Y[5], color = cols[2], label = '$n_{\mathrm{s}} = 2^{8}$')

    a = _getylims(X[3], Y[3], bounds)
    ax2.set_ylim([a[0], a[1] * 1.3])


    #ax1.set_ylim([-3.25*1e7, 2.5*1e8])
    #ax2.set_ylim([-8.5*1e7, 1.8*1e9])
    ax2.set_xlabel('$\delta$ [ppm]')

    ax1.set_ylabel('$I(\delta)$')
    ax2.set_ylabel('$I(\delta)$')

    my_axs = [ax1, ax2]
    for ax in my_axs:
        ax.set_xlim(bounds)
        ax.invert_xaxis()
        ax.set_yticks([])

    #ax2.set_ylim([-10, 355])

    my_axs = [ax1, ax2, ax3]

    for ax in my_axs:
        ax.tick_params(axis='x', direction='in', which='both')
        ax.tick_params(axis='y', direction='in', which='both')


    x_pos = sum(lactic_acid_bounds)/2.
    y_pos = 1.3*1e11

    ax1y = ax1.get_ylim()
    ax2y = ax2.get_ylim()

    ax1.annotate('$%s \mathrm{mmol/L}$' % concs_oth_lab[0], xy=(1.24, ax1y[1] * 0.87), xytext=(1.24, ax1y[1] * 0.87), horizontalalignment='center', color = 'k')
    ax2.annotate('$%s \mathrm{mmol/L}$' % concs_oth_lab[-1], xy=(1.24, ax2y[1] * 0.9), xytext=(1.24, ax2y[1] * 0.9), horizontalalignment='center', color = 'k')

    ax1.legend(loc = 'upper left', frameon = False, fontsize = 10)

    x = nscans_long_lab  # Independent variable (e.g., number of scans)
    y = MEANS[-1]        # Dependent variable (mean SNR)
    yerr = SEMS[-1]      # Standard error on the mean
    # Log-transform the data
    log_x = np.log(x)
    log_y = np.log(y)
    log_yerr = yerr / y  # Approximate error propagation in log space

    # Perform linear regression in log-log space
    coeffs = np.polyfit(log_x, log_y, 1, w=1/log_yerr)  # Weighted fit in log-log space
    slope, intercept = coeffs
    A = np.exp(intercept)  # Convert intercept back to original scale
    
    ax3.plot(nscans_long_lab, np.ones(len(nscans_long_lab)) * 3., color = 'grey', alpha = 0.5, ls = ':')
    ax3.plot(nscans_long_lab, np.ones(len(nscans_long_lab)) * 10., color = 'grey', alpha = 0.5, ls = ':')

    # Generate fitted values for plotting
    x_fit = np.linspace(min(x), max(x), 100)
    y_fit = A * x_fit**slope  # Convert fit back to original scale
    ax3.plot(x_fit, y_fit, '--', alpha = 0.5, color = 'tab:red')
    #p_value = goodness_of_fit_p_value(x, MEANS[-1], SEMS[-1], A, slope, scale_factor=1.0)
    #print(p_value)

    A_half = A / 2
    y_fit_half = A_half * x_fit**slope  # Half SNR power-law fit
    ax3.plot(x_fit, y_fit_half, '--', color = 'tab:green')
    #p_value = goodness_of_fit_p_value(x, MEANS[-2], SEMS[-2], A, slope, scale_factor=0.5)
    #print(p_value1)
    
    A_half = A / 4
    y_fit_half = A_half * x_fit**slope  # Half SNR power-law fit
    ax3.plot(x_fit, y_fit_half, '--', color = 'tab:orange')
    #p_value = goodness_of_fit_p_value(x, MEANS[-3], SEMS[-3], A, slope, scale_factor=0.25)
    #print(p_value2)

    A_half = A / 8
    y_fit_half = A_half * x_fit**slope  # Half SNR power-law fit
    ax3.plot(x_fit, y_fit_half, '--', color = 'tab:blue')
    #p_value = goodness_of_fit_p_value(x, MEANS[-4], SEMS[-4], A, slope, scale_factor=0.125)
    #print(p_value3)



    #Ps = [p_value, p_value1, p_value2, p_value3]
    cols = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
    for i in range(0, len(SNRs)):
        snr = SNRs[i]
        ax3.errorbar(nscans_long_lab, MEANS[i], yerr=SEMS[i], fmt='o', label='$%s \mathrm{mmol/L}$' % concs_oth_lab[i], capsize=5, color=cols[i])


    annotations = ['Best fit to data', '1/2 best fit', '1/4 best fit', '1/8 best fit']
    x_annotate = [15, 15, 15, 15]  # Example x-coordinates to place the text
    for i, (A_factor, text) in enumerate(zip([1, 0.5, 0.25, 0.125], annotations)):
        y_annotate = A * A_factor * np.array(x_annotate) ** slope  # Compute y position based on the line equation
        
        # Add annotation rotated to match the line's slope
        ax3.annotate(text, 
                    (x_annotate[i], y_annotate[i]*1.18),  # Position of annotation
                    xytext=(10, 0),  # Offset text a little to the right
                    textcoords='offset points',
                    fontsize=12,
                    color=np.flip(cols)[i],
                    rotation=np.degrees(np.arctan(slope)) - 4,  # Rotate parallel to line
                    ha='left')  # Align left

    #print(np.degrees(np.arctan(slope)))
    ax3.legend(frameon = False)
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_xlabel('Number of Scans')
    ax3.set_ylabel('$\mathrm{SNR}$')
    # Adjust layout
    plt.tight_layout()

    plt.title('Substance = %s, Pulse = %s, SNR = %s' % (substance, pulse_seq, snr_method))
    ## Show plot

    fnm = '%s_%s_%s.png' % (substance, pulse_seq, snr_method)
    plt.savefig(fnm, bbox_inches = 'tight', pad_inches = 0.05)
    plt.close()




#snr_fig(substance = 'Glucose', pulse_seq = 'zg')
#snr_fig(substance = 'Glucose', pulse_seq = 'zg30')
#
#snr_fig(substance = 'Lactate', pulse_seq = 'zg')
#snr_fig(substance = 'Lactate', pulse_seq = 'zg30')
#
#snr_fig(substance = 'Citrate', pulse_seq = 'zg')
#snr_fig(substance = 'Citrate', pulse_seq = 'zg30')


SIGS, NOISES, SNRs = getSNRs(substance = 'Glucose', pulse_seq = 'zg', plot_spectra = False, plot_snrs = False)
mat_g = np.nanmean(SNRs, axis = 1)

SIGS, NOISES, SNRs = getSNRs(substance = 'Lactate', pulse_seq = 'zg', plot_spectra = False, plot_snrs = False)
mat_l = np.nanmean(SNRs, axis = 1)

SIGS, NOISES, SNRs = getSNRs(substance = 'Citrate', pulse_seq = 'zg', plot_spectra = False, plot_snrs = False)
mat_c = np.nanmean(SNRs, axis = 1)

SIGS, NOISES, SNRs = getSNRs(substance = 'Glucose', pulse_seq = 'zg30', plot_spectra = False, plot_snrs = False)
mat_g30 = np.nanmean(SNRs, axis = 1)

SIGS, NOISES, SNRs = getSNRs(substance = 'Lactate', pulse_seq = 'zg30', plot_spectra = False, plot_snrs = False)
mat_l30 = np.nanmean(SNRs, axis = 1)

SIGS, NOISES, SNRs = getSNRs(substance = 'Citrate', pulse_seq = 'zg30', plot_spectra = False, plot_snrs = False)
mat_c30 = np.nanmean(SNRs, axis = 1)


x_vals = [0,1,2,3,4,5,6,7,8]
x_labs = [('$2^{%s}$' % val) for val in x_vals]

y_vals = [0, 1, 2, 3]
concs_glc_lab = [('%smM'% val) for val in np.array(glucose_concs)/1000.]
concs_lac_lab = [('%smM'% val) for val in np.array(lactic_acid_concs)/1000.]
concs_cit_lab = [('%smM'% val) for val in np.array(citric_acid_concs)/1000.]



# Example matrices (replace with your actual data)
matrices = [mat_g, mat_l, mat_c, mat_g30, mat_l30, mat_c30]

# Different y-axis labels for each row
y_labels_list = [
    concs_glc_lab,  # Row 1
    concs_lac_lab,  # Row 2
    concs_cit_lab  # Row 3
    ]

# X-axis labels (same for all)
x_vals = [0, 1, 2, 3, 4, 5, 6, 7, 8]
x_labels = [('$2^{%s}$' % val) for val in x_vals]

zg30_times  = [13, 16, 24, 37, 66, 123, 236, 465, 923]
zg_times    = [9,  17, 32, 62, 123, 244, 487, 971, 1940]

nscans_long = np.arange(0, 9, 1)
nscans_long_lab = np.array([2**val for val in nscans_long])

#fig, axs = plt.subplots(1,1)
#axs.plot(np.log10(zg_times), label = 'zg')
#axs.plot(np.log10(zg30_times), label = 'zg30')
#axs.set_xticks(nscans_long)
#axs.set_xticklabels(nscans_long_lab )
#plt.legend()
#plt.xlabel('scans')
#plt.ylabel('secs')
#plt.savefig('timings.png')


def plot_heatmap_grid(matrices, x_labels, y_labels_list):
    matplotlib.rc('text', usetex=True)
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.size'] = 14
    
    # Create figure with more right margin for colorbar
    fig = plt.figure(figsize=(8, 7))  # Slightly wider figure
    gs = fig.add_gridspec(3, 2, width_ratios=[1, 1])
    axes = np.array([[fig.add_subplot(gs[i, j]) for j in range(2)] for i in range(3)])
    
    # Find global color limits for consistency
    vmin = min(np.log10(mat).min() for mat in matrices)
    vmax = max(np.log10(mat).max() for mat in matrices)
    newcmp = tol_cmap('sunset')

    # Plot heatmaps
    im = axes[0,0].imshow(np.log10(mat_g), aspect='auto', cmap=newcmp, vmin=vmin, vmax=vmax, origin='lower', extent=[-0.5, len(x_labels)-0.5, -0.5, 3.5])
    im = axes[1,0].imshow(np.log10(mat_l), aspect='auto', cmap=newcmp, vmin=vmin, vmax=vmax, origin='lower', extent=[-0.5, len(x_labels)-0.5, -0.5, 3.5])
    im = axes[2,0].imshow(np.log10(mat_c), aspect='auto', cmap=newcmp, vmin=vmin, vmax=vmax, origin='lower', extent=[-0.5, len(x_labels)-0.5, -0.5, 3.5])
    im = axes[0,1].imshow(np.log10(mat_g30), aspect='auto', cmap=newcmp, vmin=vmin, vmax=vmax, origin='lower', extent=[-0.5, len(x_labels)-0.5, -0.5, 3.5])
    im = axes[1,1].imshow(np.log10(mat_l30), aspect='auto', cmap=newcmp, vmin=vmin, vmax=vmax, origin='lower', extent=[-0.5, len(x_labels)-0.5, -0.5, 3.5])
    im = axes[2,1].imshow(np.log10(mat_c30), aspect='auto', cmap=newcmp, vmin=vmin, vmax=vmax, origin='lower', extent=[-0.5, len(x_labels)-0.5, -0.5, 3.5])
    for ax in axes.flat:
        ax.set_xlim(-0.5, len(x_labels)-0.5)
        ax.set_ylim(-0.5, 3.5)
    #im = axes[0,0].imshow(np.log10(mat_g), aspect='auto', cmap=newcmp, vmin=vmin, vmax=vmax, 
    #                 origin='lower', extent=[-0.5, len(x_labels)-0.5, -0.5, 3.5])
    # Set axis labels and ticks
    for ax in axes.flat:
        ax.set_xticks(range(len(x_labels)))
        ax.set_xticklabels([''] * 9)

    
    axes[2,0].set_xticklabels(x_labels, rotation=45)
    axes[2,1].set_xticklabels(x_labels, rotation=45)
    
    # Set y-ticks for left column
    for i in range(3):
        axes[i,0].set_yticks([0, 1, 2, 3])
        axes[i,1].set_yticks([0, 1, 2, 3])
        axes[i,0].set_yticklabels(y_labels_list[i])
        axes[i,1].set_yticklabels([''] * 4)
        
    
    # Add labels
    axes[2,0].set_xlabel("Number of Scans (log scale)")
    axes[2,1].set_xlabel("Number of Scans (log scale)")
    axes[0,0].set_ylabel("Glucose Conc. (mM)")
    axes[1,0].set_ylabel("Lactate Conc. (mM)")
    axes[2,0].set_ylabel("Citrate Conc. (mM)")
    
    # Adjust subplot spacing
    # Modify the subplot spacing and margins
    plt.subplots_adjust(
        right=0.85,        # Reduce this slightly if colorbar is cut off
        hspace=0.07,       # Reduce vertical spacing between subplots
        wspace=0.05,       # Reduce horizontal spacing between subplots
        left=0.15,         # Add left margin for y-labels
        bottom=0.1,        # Add bottom margin for x-labels
        top=0.95          # Add top margin
    )    # Add shared colorbar
    top_pos = axes[0,1].get_position()

    bottom_pos = axes[2,1].get_position()

    colorbar_left = bottom_pos.x1 + 0.03
    colorbar_bottom = bottom_pos.y0
    colorbar_width = 0.03
    colorbar_height = top_pos.y1 - bottom_pos.y0
    colorbar_rect = [colorbar_left, colorbar_bottom, colorbar_width, colorbar_height]


    # Add colorbar
    cbar_ax = fig.add_axes(colorbar_rect)  # [left, bottom, width, height]
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label('log10(SNR)')
    cbar.ax.axhline(np.log10(3), c='orange', lw = 1.5)
    cbar.ax.axhline(1., c='green', lw = 1.5)
    # Inside your plotting loop, after the imshow commands, add:
    levels = [np.log10(3), np.log10(10)]  # The contour levels you want

    def create_stepped_boundary(data, threshold, grid_height, grid_width):
        """Create stepped boundary line segments at a given threshold value"""
        boundaries = []
        
        # Create a binary mask where True indicates values above threshold
        mask = data >= threshold
        
        # Scan each cell to find boundaries
        for i in range(grid_height):
            for j in range(grid_width):
                if mask[i,j]:
                    # Check cell above (if not at top edge)
                    if i > 0 and not mask[i-1,j]:
                        boundaries.append([(j, i), (j+1, i)])
                    # Check cell below (if not at bottom edge)
                    if i < grid_height-1 and not mask[i+1,j]:
                        boundaries.append([(j, i+1), (j+1, i+1)])
                    # Check cell left (if not at left edge)
                    if j > 0 and not mask[i,j-1]:
                        boundaries.append([(j, i), (j, i+1)])
                    # Check cell right (if not at right edge)
                    if j < grid_width-1 and not mask[i,j+1]:
                        boundaries.append([(j+1, i), (j+1, i+1)])
        
        return boundaries

    
    # In your plotting code, after imshow:

    axs = [axes[0,0], axes[1,0],axes[2,0],axes[0,1],axes[1,1],axes[2,1]]
    
    for I in range(0, len(axs)):
        data = np.log10(matrices[I])
        height, width = data.shape
        
        # Get boundaries for both thresholds
        boundaries_low = create_stepped_boundary(data, np.log10(3), height, width)
        boundaries_high = create_stepped_boundary(data, np.log10(10), height, width)
        
        # Plot the boundaries
        for line in boundaries_low:
            axs[I].plot([line[0][0]- 0.5, line[1][0]- 0.5], 
                        [line[0][1]- 0.5, line[1][1]- 0.5], 
                        color='orange', linewidth=1.3)
        
        for line in boundaries_high:
            axs[I].plot([line[0][0] - 0.5, line[1][0] - 0.5], 
                        [line[0][1] - 0.5, line[1][1] - 0.5], 
                        color='green', linewidth=1.3)

        axs[I].xaxis.set_tick_params(direction='in', which='both', right=True, top=True)
        axs[I].yaxis.set_tick_params(direction='in', which='both', right=True, top=True)


    ax_zg = axes[0,0].twiny()
    ax_zg.set_xlim([0., 9.])
    ax_zg.set_xticks(np.arange(.5, 9.5, 1.))
    ax_zg.set_xticklabels(zg_times, fontsize = 10)

    ax_zg30 = axes[0,1].twiny()
    ax_zg30.set_xlim([0., 9.])
    ax_zg30.set_xticks(np.arange(.5, 9.5, 1.))
    ax_zg30.set_xticklabels(zg30_times, fontsize = 10)

    ax_zg.xaxis.set_tick_params(direction='in', which='both', right=True, top=True)
    ax_zg30.xaxis.set_tick_params(direction='in', which='both', right=True, top=True)


    '''
    ax_zg = axes[0,0].twiny()
    ax_zg30 = axes[0,1].twiny()

    ax_zg.set_xticks(nscans_long)
    ax_zg.set_xticklabels(zg_times)


    ax_zg30.set_xticks(nscans_long)
    ax_zg30.set_xticklabels(zg30_times)
    '''
    #ax2.set_xlabel('Length of zg experiment (s)')  # we already handled the x-label with ax1
    #ax3.set_xlabel('Length of zg30 experiment (s)')  # we already handled the x-label with ax1

    axes[0,0].set_title('zg Experiment time (secs)', fontsize = 12)
    axes[0,1].set_title('zg30 Experiment time (secs)', fontsize = 12)


    #plt.show()
    plt.savefig('temp_hist2d.pdf', bbox_inches='tight', dpi=300)
    plt.close()

# Call the function
plot_heatmap_grid(matrices, x_labels, y_labels_list)

matrices  = matrices

x_labels = x_labels

y_labels_list = y_labels_list