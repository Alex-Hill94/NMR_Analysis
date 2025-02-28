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
from scipy.stats import norm

lactic_acid_bounds = [1.17, 1.50]
glucose_bounds = [3.19, 3.98]
citric_acid_old_bounds = [2.37, 2.82]
mid = np.mean(citric_acid_old_bounds)
rang = np.max(citric_acid_old_bounds) - np.min(citric_acid_old_bounds)
sep = rang/6.
citric_acid_bounds = [mid - sep, mid + sep]
citric_acid_bounds = citric_acid_old_bounds
noise_bounds       = [-2.0, -1.0]

lactic_acid_concs = [750, 1500, 3000, 6000]
glucose_concs = [1250, 2500, 5000, 10000]
citric_acid_concs = [50,  100,  200, 400]

scans = [1, 2, 4, 8, 16, 32, 64, 128, 256]
sample = 'D24'
snr_choice = 'liv'
times = []

def two_pan_exp_sim():
    S = LoadSpectra()
    nscans = [256, 128, 64, 32, 16, 8, 4, 2, 1]
    EXP = []
    SIM = []
    Y_tots = []
    Y_exps = []
    X_vals = []
    for nscan in nscans:
        S.ReadTextFile(nscan = nscan, sample =sample)
        #S.SubtractWater()
        #X_ORIG.append(S.water_sub_ppm)
        #Y_ORIG.append(S.water_sub_amplitude)
        A = AnalyseSpectra()
        A.InputData(x = S.initial_ppm, y = S.initial_amplitude)
        A.FitTSP()
        A.ScaleSpectra(x_shift = (-1) * A.tsp_centre)

        def glucose_fit(x, scale, x_shift, width):
            return Glucose(x_array=x, scale=scale, x_shift=x_shift, width=width)
        def lactate_fit(x, scale, x_shift, width):
            return Lactate(x_array=x, scale=scale, x_shift=x_shift, width=width)
        def citrate_fit(x, scale, x_shift, width):
            return Citrate(x_array=x, scale=scale, x_shift=x_shift, width=width)

        # Extract X and Y data
        X = A.processed_ppm
        Y = A.processed_amplitude
        glc_mask = (X >= min(glucose_bounds)) * (X <= max(glucose_bounds))
        lac_mask = (X >= min(lactic_acid_bounds)) * (X <= max(lactic_acid_bounds))
        cit_mask = (X >= min(citric_acid_bounds)) * (X <= max(citric_acid_bounds))

        def get_fit(substance = 'Glucose', plot = False):
            if substance == 'Glucose':
                mask = glc_mask
                fnc = Glucose
                fit_fnc = glucose_fit
            elif substance == 'Lactate':
                mask = lac_mask
                fnc = Lactate
                fit_fnc = lactate_fit
            elif substance == 'Citrate':
                mask = cit_mask
                fnc = Citrate
                fit_fnc = citrate_fit

            X_RED = X[mask]
            Y_RED = Y[mask]

            # Initial guesses
            initial_guesses = [
                np.nanmax(Y_RED),  # scale
                0,  # x_shift
                2 # linewidth
            ]

            # Perform curve fitting
            popt, pcov = curve_fit(fit_fnc, X_RED, Y_RED, p0=initial_guesses)

            # Get optimized parameters
            opt_scale, opt_x_shift, opt_width = popt

            Y_fit = fnc(x_array=X, scale=opt_scale, x_shift=opt_x_shift, width=opt_width)
            # Residuals
            # Predicted values
            Y_fit_red = fit_fnc(X_RED, *popt)
            residuals = Y_RED - Y_fit_red

            # R-squared
            SS_res = np.sum(residuals**2)
            SS_tot = np.sum((Y_RED - np.mean(Y_RED))**2)
            R_squared = 1 - (SS_res / SS_tot)

            # RMSE
            RMSE = np.sqrt(np.mean(residuals**2))

            # Reduced Chi-Square (if uncertainties are known)
            # Assuming uniform errors for demonstration, replace with actual errors if available
            errors = np.full_like(Y_RED, np.std(residuals))  # Example: constant error estimate
            chi_squared = np.sum((residuals / errors) ** 2)
            reduced_chi_squared = chi_squared / (len(Y_RED) - len(popt))  # DOF = N - num_params

            print(f"R-squared: {R_squared:.4f}")
            print(f"RMSE: {RMSE:.4f}")
            print(f"Reduced Chi-Square: {reduced_chi_squared:.4f}")

            if plot:
                # Plot results
                plt.figure()
                plt.plot(X, Y, label="Original Data")
                plt.plot(X, Y_fit, label="Fitted Curve", linestyle="--")
                plt.legend()
                plt.show()
            return Y_fit
        #### SIMS
        y_glc = get_fit(substance = 'Glucose')
        y_lac = get_fit(substance = 'Lactate')
        y_cit = get_fit(substance = 'Citrate')
        y_tot = y_glc + y_lac + y_cit
        s_glc = sum(y_glc[glc_mask])
        s_lac = sum(y_lac[lac_mask])
        s_cit = sum(y_cit[cit_mask])
        f_lg = s_lac/s_glc
        f_cg = s_cit/s_glc
        f_cl = s_cit/s_lac

        SIM.append([f_lg, f_cg, f_cl])

        #### EXPS
        s_glc_exp = sum(Y[glc_mask])
        s_lac_exp = sum(Y[lac_mask])
        s_cit_exp = sum(Y[cit_mask])
        f_lg_exp = s_lac_exp/s_glc_exp
        f_cg_exp = s_cit_exp/s_glc_exp
        f_cl_exp = s_cit_exp/s_lac_exp
        EXP.append([f_lg_exp, f_cg_exp, f_cl_exp])

        if nscan == 256 or nscan == 32:
            Y_tots.append(y_tot)
            Y_exps.append(Y)
            X_vals.append(X)

    return np.array(nscans), np.array(SIM), np.array(EXP), Y_tots, Y_exps, X_vals

def two_pan_exp_sim_with_uncertainties():
    S = LoadSpectra()
    nscans = [256, 128, 64, 32, 16, 8, 4, 2, 1]
    EXP = []
    SIM = []
    SIM_UNCERTAINTIES = []
    EXP_UNCERTAINTIES = []
    Y_tots = []
    Y_exps = []
    X_vals = []
    
    for nscan in nscans:
        print(nscan)
        S.ReadTextFile(nscan=nscan, sample=sample)
        A = AnalyseSpectra()
        A.InputData(x=S.initial_ppm, y=S.initial_amplitude)
        A.FitTSP()
        A.ScaleSpectra(x_shift=(-1) * A.tsp_centre)

        def glucose_fit(x, scale, x_shift, width):
            return Glucose(x_array=x, scale=scale, x_shift=x_shift, width=width)
        def lactate_fit(x, scale, x_shift, width):
            return Lactate(x_array=x, scale=scale, x_shift=x_shift, width=width)
        def citrate_fit(x, scale, x_shift, width):
            return Citrate(x_array=x, scale=scale, x_shift=x_shift, width=width)

        # Extract X and Y data
        X = A.processed_ppm
        Y = A.processed_amplitude
        glc_mask = (X >= min(glucose_bounds)) * (X <= max(glucose_bounds))
        lac_mask = (X >= min(lactic_acid_bounds)) * (X <= max(lactic_acid_bounds))
        cit_mask = (X >= min(citric_acid_bounds)) * (X <= max(citric_acid_bounds))
        noise_mask = (X >= min(noise_bounds)) * (X <= max(noise_bounds))
        noise_std = np.std(Y[noise_mask])
        n_points_glc = np.sum(glc_mask)
        n_points_lac = np.sum(lac_mask)
        n_points_cit = np.sum(cit_mask)

        # 1. Integration uncertainty
        int_uncrt_glc = noise_std * np.sqrt(n_points_glc)
        int_uncrt_lac = noise_std * np.sqrt(n_points_lac)
        int_uncrt_cit = noise_std * np.sqrt(n_points_cit)
        IU = [int_uncrt_glc, int_uncrt_lac, int_uncrt_cit]

        # 2. Baseline uncertainty
        # Estimate baseline offset using regions around peak
        BU = []
        for K in range(0, len([glc_mask, lac_mask, cit_mask])):
            mask = [glc_mask, lac_mask, cit_mask][K]
            n_points = [n_points_glc, n_points_lac, n_points_cit][K]
            peak_start = X[mask][0]
            peak_end = X[mask][-1]
            pre_peak = Y[(X > peak_start - 0.1) & (X < peak_start)]
            post_peak = Y[(X > peak_end) & (X < peak_end + 0.1)]
            baseline_level = np.mean(np.concatenate([pre_peak, post_peak]))
            #print(baseline_level)
            #print(n_points)
            baseline_uncertainty = np.abs(baseline_level) * n_points
            BU.append(baseline_uncertainty)

        total_uncertainty_exp = []
        for K in range(0, len([glc_mask, lac_mask, cit_mask])):
            TU = np.sqrt(IU[K] **2 + BU[K]**2)
            TU = np.sqrt(IU[K] **2)
            total_uncertainty_exp.append(TU)
        
        tu_glc = total_uncertainty_exp[0]
        tu_lac = total_uncertainty_exp[1]
        tu_cit = total_uncertainty_exp[2]
        

        
        def get_fit(substance='Glucose', plot=False):
            if substance == 'Glucose':
                mask = glc_mask
                fnc = Glucose
                fit_fnc = glucose_fit
            elif substance == 'Lactate':
                mask = lac_mask
                fnc = Lactate
                fit_fnc = lactate_fit
            elif substance == 'Citrate':
                mask = cit_mask
                fnc = Citrate
                fit_fnc = citrate_fit

            X_RED = X[mask]
            Y_RED = Y[mask]

            # Initial guesses
            initial_guesses = [
                np.nanmax(Y_RED),  # scale
                0,  # x_shift
                2  # linewidth
            ]

            # Perform curve fitting
            popt, pcov = curve_fit(fit_fnc, X_RED, Y_RED, p0=initial_guesses, bounds=([0, -np.inf, -np.inf], np.inf))

            # Get optimized parameters
            opt_scale, opt_x_shift, opt_width = popt

            Y_fit = fnc(x_array=X, scale=opt_scale, x_shift=opt_x_shift, width=opt_width)

            if plot:
                plt.figure()
                plt.plot(X, Y, label="Original Data")
                plt.plot(X, Y_fit, label="Fitted Curve", linestyle="--")
                plt.legend()
                plt.show()
                
            return Y_fit, popt, pcov

        # Get fits and fitting parameters
        y_glc, popt_glc, pcov_glc = get_fit(substance='Glucose')
        y_lac, popt_lac, pcov_lac = get_fit(substance='Lactate')
        y_cit, popt_cit, pcov_cit = get_fit(substance='Citrate')

        # Calculate ratios with uncertainties
        def calculate_ratio_uncertainties(X, Y, y_fits, masks, popt_dict, pcov_dict):
            def get_area_uncertainty(metabolite):
                mask = masks[metabolite]
                residuals = Y[mask] - y_fits[metabolite][mask]
                noise_std = np.std(residuals)
                
                n_points = np.sum(mask)
                params = popt_dict[metabolite]
                pcov = pcov_dict[metabolite]
                
                # Monte Carlo simulation
                n_samples = 1000
                param_samples = np.random.multivariate_normal(params, pcov, n_samples)
                areas = []
                
                for p in param_samples:
                    if metabolite == 'Glucose':
                        y = glucose_fit(X[mask], *p)
                    elif metabolite == 'Lactate':
                        y = lactate_fit(X[mask], *p)
                    else:  # Citrate
                        y = citrate_fit(X[mask], *p)
                    areas.append(np.sum(y))
                    
                fitting_std = np.std(areas)
                total_uncertainty = np.sqrt(noise_std**2 + fitting_std**2)
                
                return np.sum(y_fits[metabolite][mask]), total_uncertainty

            # Get areas and uncertainties
            areas = {}
            uncertainties = {}
            for met in ['Glucose', 'Lactate', 'Citrate']:
                areas[met], uncertainties[met] = get_area_uncertainty(met)

            def ratio_uncertainty(num, num_err, denom, denom_err):
                ratio = num / denom
                rel_error = np.sqrt((num_err/num)**2 + (denom_err/denom)**2)
                return ratio, abs(ratio) * rel_error

            # Calculate ratios and uncertainties
            ratios = {
                'Lactate/Glucose': ratio_uncertainty(
                    areas['Lactate'], uncertainties['Lactate'],
                    areas['Glucose'], uncertainties['Glucose']
                ),
                'Citrate/Glucose': ratio_uncertainty(
                    areas['Citrate'], uncertainties['Citrate'],
                    areas['Glucose'], uncertainties['Glucose']
                ),
                'Citrate/Lactate': ratio_uncertainty(
                    areas['Citrate'], uncertainties['Citrate'],
                    areas['Lactate'], uncertainties['Lactate']
                )
            }
            
            return ratios

        # Calculate uncertainties for this scan
        y_fits = {'Glucose': y_glc, 'Lactate': y_lac, 'Citrate': y_cit}
        masks = {'Glucose': glc_mask, 'Lactate': lac_mask, 'Citrate': cit_mask}
        popt_dict = {'Glucose': popt_glc, 'Lactate': popt_lac, 'Citrate': popt_cit}
        pcov_dict = {'Glucose': pcov_glc, 'Lactate': pcov_lac, 'Citrate': pcov_cit}
        
        ratios = calculate_ratio_uncertainties(
            X=X, 
            Y=Y, 
            y_fits=y_fits,
            masks=masks,
            popt_dict=popt_dict,
            pcov_dict=pcov_dict
        )

        # Extract ratio values and uncertainties
        sim_ratios = []
        uncertainties = []
        for ratio_name in ['Lactate/Glucose', 'Citrate/Glucose', 'Citrate/Lactate']:
            ratio, uncertainty = ratios[ratio_name]
            sim_ratios.append(ratio)
            uncertainties.append(uncertainty)

        SIM.append(sim_ratios)
        SIM_UNCERTAINTIES.append(uncertainties)

        # Calculate experimental ratios
        y_tot = y_glc + y_lac + y_cit
        s_glc_exp = sum(Y[glc_mask])
        s_lac_exp = sum(Y[lac_mask])
        s_cit_exp = sum(Y[cit_mask])
        f_lg_exp = s_lac_exp/s_glc_exp
        f_cg_exp = s_cit_exp/s_glc_exp
        f_cl_exp = s_cit_exp/s_lac_exp
        EXP.append([f_lg_exp, f_cg_exp, f_cl_exp])

        re_lg = np.sqrt(
            (tu_lac/s_lac_exp)**2 +
            (tu_glc/s_glc_exp)**2
        )
        U_lg = abs(f_lg_exp) * re_lg

        re_cg = np.sqrt(
            (tu_cit/s_cit_exp)**2 +
            (tu_glc/s_glc_exp)**2
        )
        U_cg = abs(f_cg_exp) * re_cg

        re_cl = np.sqrt(
            (tu_cit/s_cit_exp)**2 +
            (tu_lac/s_lac_exp)**2
        )
        U_cl = abs(f_cl_exp) * re_cl
        EXP_UNCERTAINTIES.append([U_lg,U_cg,U_cl])

        if nscan == 256 or nscan == 32:
            Y_tots.append(y_tot)
            Y_exps.append(Y)
            X_vals.append(X)

    return np.array(nscans), np.array(SIM), np.array(EXP), np.array(SIM_UNCERTAINTIES), np.array(EXP_UNCERTAINTIES), Y_tots, Y_exps, X_vals

n, s, e, u, ue, yt, y, x = two_pan_exp_sim_with_uncertainties()
'''
fig, axs = plt.subplots(1,1, figsize = [6, 6])

axs.plot(n, s[:,0], label = 'Lactate/Glucose Sim', color = 'tab:blue')
axs.plot(n, e[:,0], label = 'Lactate/Glucose Exp', color = 'tab:blue', ls = '--')
axs.plot(n, s[:,1], label = 'Citrate/Glucose Sim', color = 'tab:red')
axs.plot(n, e[:,1], label = 'Citrate/Glucose Exp', color = 'tab:red', ls = '--')

axs.set_xlabel('Nscans')
axs.set_ylabel('Ratio')
axs.legend()
axs.set_xscale('log')


plt.show()
'''

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

X = x[0]
glc_mask = (X >= min(glucose_bounds)) * (X <= max(glucose_bounds))
lac_mask = (X >= min(lactic_acid_bounds)) * (X <= max(lactic_acid_bounds))
cit_mask = (X >= min(citric_acid_bounds)) * (X <= max(citric_acid_bounds))


matplotlib.rc('text', usetex=True)

# Set the default font to Times New Roman
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 15

import matplotlib as mpl
cmap = mpl.colormaps['BrBG']

colours = cmap(np.linspace(0, 1, len(scans) + 2))

# Create figure
fig = plt.figure(figsize=(12, 5))

# Create GridSpec with 2 rows, 2 columns
# Set width_ratios to make left panels wider (e.g., 2:1 ratio)
gs = GridSpec(2, 2, figure=fig, width_ratios=[2, 1])

# Create the three subplots
ax1 = fig.add_subplot(gs[0, 0])  # Top left
ax2 = fig.add_subplot(gs[1, 0])  # Bottom left
ax3 = fig.add_subplot(gs[:, 1])  # Right (spans both rows)

glucose_colour = 'tab:blue'
citrate_colour = 'tab:green'
lactate_colour = 'tab:red'
noise_colour   = 'black'

y_max = np.nanmax(y[0][glc_mask])

ax1.plot(x[0], y[0]/y_max, color = 'darkslategrey', alpha = 1.0, ls = '-')
ax1.plot(x[0], yt[0]/y_max, color = 'tab:pink', lw = 1., alpha = 1.)

y_max = np.nanmax(y[1][glc_mask])

ax2.plot(x[1], y[1]/y_max, color = 'darkslategrey', alpha = 1.0, ls = '-')
ax2.plot(x[1], yt[1]/y_max, color = 'tab:pink', lw = 1., alpha = 1.)

#ax1.set_ylim([-1 * 1e8, 1.5*1e9])
#ax2.set_ylim([-1.27 * 1e7, 1.9*1e8])
ax1.set_ylim([-0.1, 1.3])
ax2.set_ylim([-0.1, 1.3])

ax1.set_xlim([0.9, 5.4])
ax2.set_xlim([0.9, 5.4])

ax1.invert_xaxis()
ax2.invert_xaxis()
custom_lines  = [(Line2D([0], [0], color='tab:pink',        linestyle = '-')),
                (Line2D([0], [0],          color='darkslategrey',         linestyle = '-'))]
ax1.legend(custom_lines, ['Simulation Fit', 'Experimental Data'], frameon = False, fontsize = 12)

#ax3.errorbar(n, s[:,0], yerr = u[:,0], label = 'Lactate/Glucose Sim', color = 'tab:purple')
#
#ax3.errorbar(n, s[:,1], yerr = u[:, 1], label = 'Citrate/Glucose Sim', color = 'tab:cyan')
#
#ax3.errorbar(n, e[:,0], yerr = ue[:,0], label = 'Lactate/Glucose Exp', color = 'tab:purple', ls = '--')
#
#ax3.errorbar(n, e[:,1], yerr = ue[:,1], label = 'Citrate/Glucose Exp', color = 'tab:cyan', ls = '--')
ax3.plot(n, s[:,0], label='Lactate/Glucose Sim', color='tab:purple', lw = 0.8)
ax3.fill_between(n, s[:,0] - u[:,0], s[:,0] + u[:,0], 
                 color='tab:purple', alpha=0.3)

ax3.plot(n, s[:,1], label='Citrate/Glucose Sim', color='tab:cyan', lw = 0.8)
ax3.fill_between(n, s[:,1] - u[:,1], s[:,1] + u[:,1], 
                 color='tab:cyan', alpha=0.3)

# Experimental data with uncertainty bands
ax3.plot(n, e[:,0], label='Lactate/Glucose Exp', color='tab:purple', ls='--', lw = 0.8)
ax3.fill_between(n, e[:,0] - ue[:,0], e[:,0] + ue[:,0], 
                 color='tab:purple', alpha=0.3, linestyle='--')

ax3.plot(n, e[:,1], label='Citrate/Glucose Exp', color='tab:cyan', ls='--', lw = 0.8)
ax3.fill_between(n, e[:,1] - ue[:,1], e[:,1] + ue[:,1], 
                 color='tab:cyan', alpha=0.3, linestyle='--')
ax3.set_xscale('log')
my_axs = [ax1, ax2, ax3]
fill_alpha = 0.07
for ax in my_axs:
    ax.tick_params(axis='x', direction='in', which='both', right=True, top=True)
    ax.tick_params(axis='y', direction='in', which='both', right=True, top=True)
for ax in my_axs[:2]:
    llim, ulim = ax.get_ylim()
    ax.fill_betweenx([llim, ulim], citric_acid_bounds[0], citric_acid_bounds[1], color= citrate_colour, alpha=fill_alpha)
    ax.fill_betweenx([llim, ulim], lactic_acid_bounds[0], lactic_acid_bounds[1], color= lactate_colour, alpha=fill_alpha)
    ax.fill_betweenx([llim, ulim], glucose_bounds[0], glucose_bounds[1], color=glucose_colour, alpha=fill_alpha)
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(['0.0','', '0.5', '', '1.0'])
    ax.grid(True, linestyle='--', alpha=0.2)  # Add gridlines for clarity


ax1.set_xticklabels([''])
#ax.text(x_ann, y_ann, lab, bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'), fontsize=10)

ax1.annotate('256 Scans', xy=(0.02, 0.85), xytext=(0.02, 0.85), horizontalalignment='left', xycoords = 'axes fraction', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'), fontsize=10)
ax2.annotate('32 Scans', xy=(0.02, 0.85), xytext=(0.02, 0.85), horizontalalignment='left', xycoords = 'axes fraction', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'), fontsize=10)
ax2.set_xlabel('$\delta$ [ppm]')
ax1.set_ylabel('$I(\delta)$')
ax2.set_ylabel('$I(\delta)$')

ax3.set_xlabel('Number of Scans')
ax3.set_ylabel('Signal Ratio')

custom_lines  = [(Line2D([0], [0.5], color='tab:purple',        linestyle = '-', linewidth = 3)),
        (Line2D([0], [0.5],          color='tab:cyan',         linestyle = '-', linewidth = 3)),
        (Line2D([0], [0],          color='k',         linestyle = '-')),
        (Line2D([0], [0],          color='k',         linestyle = '--'))]

ax3.legend(custom_lines, ['Lac/Glc', 'Cit/Glc', 'Sim', 'Exp'], frameon = False, fontsize = 12, loc = 'lower right')
ax3.grid(True, linestyle='--', alpha=0.6)  # Add gridlines for clarity

plt.tight_layout()



glc_label = '\\textbf{Glucose}'
cit_label = '\\textbf{Citrate}'
lac_label = '\\textbf{Lactate}'

FS = 9

x_pos_glc = sum(glucose_bounds)/2.
y_pos_main = 1.15
ax1.annotate(glc_label, xy=(x_pos_glc, y_pos_main), xytext=(x_pos_glc, y_pos_main), horizontalalignment='center', color = glucose_colour, fontsize = FS)


#ax_mid_left.annotate(glc_label, xy=(x_pos_glc, y_pos), xytext=(x_pos_glc, y_pos), horizontalalignment='center', color = glucose_colour)
x_pos_scaled = 4.12



x_pos = sum(citric_acid_bounds)/2.
ax1.annotate(cit_label, xy=(x_pos, y_pos_main), xytext=(x_pos, y_pos_main), horizontalalignment='center', color = citrate_colour, fontsize = FS)
x_pos = sum(citric_acid_bounds)/2. - 0.1
#ax_mid_right.annotate(cit_label, xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = citrate_colour)

y_pos_lac = 0.625
x_pos = sum(lactic_acid_bounds)/2.
ax1.annotate(lac_label, xy=(x_pos, y_pos_lac), xytext=(x_pos, y_pos_lac), horizontalalignment='center', color = lactate_colour, fontsize = FS)


#plt.show()
plt.savefig('sim_test_third_rang.pdf', bbox_inches = 'tight', pad_inches = 0.05)
plt.close()

def fit():
    S = LoadSpectra()
    nscans = [256, 128, 64, 32, 16, 8, 4, 2, 1]
    EXP = []
    SIM = []
    for nscan in nscans:
        S.ReadTextFile(nscan = nscan, sample =sample)
        #S.SubtractWater()
        #X_ORIG.append(S.water_sub_ppm)
        #Y_ORIG.append(S.water_sub_amplitude)
        A = AnalyseSpectra()
        A.InputData(x = S.initial_ppm, y = S.initial_amplitude)
        A.FitTSP()
        A.ScaleSpectra(x_shift = (-1) * A.tsp_centre)


        def glucose_fit(x, scale, x_shift, width):
            return Glucose(x_array=x, scale=scale, x_shift=x_shift, width=width)
        def lactate_fit(x, scale, x_shift, width):
            return Lactate(x_array=x, scale=scale, x_shift=x_shift, width=width)
        def citrate_fit(x, scale, x_shift, width):
            return Citrate(x_array=x, scale=scale, x_shift=x_shift, width=width)

        # Extract X and Y data
        X = A.processed_ppm
        Y = A.processed_amplitude
        glc_mask = (X >= min(glucose_bounds)) * (X <= max(glucose_bounds))
        lac_mask = (X >= min(lactic_acid_bounds)) * (X <= max(lactic_acid_bounds))
        cit_mask = (X >= min(citric_acid_bounds)) * (X <= max(citric_acid_bounds))

        def get_fit(substance = 'Glucose', plot = False):
            if substance == 'Glucose':
                mask = glc_mask
                fnc = Glucose
                fit_fnc = glucose_fit
            elif substance == 'Lactate':
                mask = lac_mask
                fnc = Lactate
                fit_fnc = lactate_fit
            elif substance == 'Citrate':
                mask = cit_mask
                fnc = Citrate
                fit_fnc = citrate_fit


            X_RED = X[mask]
            Y_RED = Y[mask]

            # Initial guesses
            initial_guesses = [
                np.nanmax(Y_RED),  # scale
                0,  # x_shift
                2 # linewidth
            ]

            # Perform curve fitting
            popt, pcov = curve_fit(fit_fnc, X_RED, Y_RED, p0=initial_guesses)

            # Get optimized parameters
            opt_scale, opt_x_shift, opt_width = popt

            # Compute best-fit curve
            Y_fit = fnc(x_array=X, scale=opt_scale, x_shift=opt_x_shift, width=opt_width)
            
            if plot:
                # Plot results
                plt.figure()
                plt.plot(X, Y, label="Original Data")
                plt.plot(X, Y_fit, label="Fitted Curve", linestyle="--")
                plt.legend()
                plt.show()
            return Y_fit

        y_glc = get_fit(substance = 'Glucose')
        y_lac = get_fit(substance = 'Lactate')
        y_cit = get_fit(substance = 'Citrate')
        Y_tot = y_glc + y_lac + y_cit

        y_exp_glc = sum(Y[glc_mask])
        y_exp_lac = sum(Y[lac_mask])
        y_exp_cit = sum(Y[cit_mask])
        Y_exp_tot = y_exp_glc + y_exp_lac# + y_exp_cit
        glc_exp_frac = y_exp_glc/Y_exp_tot
        lac_exp_frac = y_exp_lac/Y_exp_tot
        cit_exp_frac = y_exp_cit/Y_exp_tot

        ref_sum = sum(Y_tot)
        glc_frac = sum(y_glc)/ref_sum
        lac_frac = sum(y_lac)/ref_sum
        cit_frac    = sum(y_cit)/ref_sum
        #fig, axs = plt.subplots(2,1, figsize = [10, 5])
        #axs[0].set_title('nscan = %s, SIM =>  GLC = %s, LAC = %s, CIT = %s' % (nscan, glc_frac, lac_frac, cit_frac))
        #axs[0].plot(X, Y, label="Original Data")
        #axs[0].plot(X, Y_tot, label="Fitted Curve", linestyle="--")
        #axs[0].legend()
        #axs[0].set_xlim([-1, 6])
        #msk = (X > 0) *(X < 4)
        #axs[0].set_ylim([np.min(Y[msk]) * 10 , np.max(Y[msk]) * 1.3])
        #axs[1].set_title('EXP => GLC = %s, LAC = %s, CIT = %s' % (glc_exp_frac, lac_exp_frac, cit_exp_frac))
        #axs[1].plot(X, Y-Y_tot, label="Residual")
        #axs[1].legend()
        #axs[1].set_xlim([-1, 6])
        #axs[1].set_ylim([np.min(Y[msk]) * 10 , np.max(Y[msk]) * 1.3])
        #plt.savefig('%s.pdf' % nscan)
        #plt.close()
        SIM.append([glc_frac, lac_frac, cit_frac])
        EXP.append([glc_exp_frac, lac_exp_frac, cit_exp_frac])
        
    SIM = np.array(SIM)
    EXP = np.array(EXP)

    #return SIM, EXP, nscans


    SIM, EXP, NSCANS = fit()
    fig, axs = plt.subplots(1,2, figsize = [10, 5])
    axs[0].plot(NSCANS, SIM[:,0]/SIM[:,0], label = 'glc')
    axs[0].plot(NSCANS, SIM[:,1]/SIM[:,1], label = 'lac')
    axs[0].plot(NSCANS, SIM[:,2]/SIM[:,2], label = 'cit')

    axs[1].plot(NSCANS, EXP[:,0]/SIM[:,0], label = 'glc')
    axs[1].plot(NSCANS, EXP[:,1]/SIM[:,1], label = 'lac')
    axs[1].plot(NSCANS, EXP[:,2]/SIM[:,2], label = 'cit')

    axs[0].set_title('SIM')
    axs[1].set_title('EXP')

    axs[0].set_xlabel('NSCANS')
    axs[1].set_xlabel('NSCANS')

    axs[0].set_ylabel('Rel. contribution to signal')
    axs[1].set_ylabel('Rel. contribution to signal')
    axs[0].legend()

    axs[0].set_xscale('log')
    axs[1].set_xscale('log')
    #axs[0].set_ylim([0.7, 1.4])
    #axs[1].set_ylim([0.7, 1.4])

    plt.show()

def standard_fits():
    substance = 'Glucose'
    pulse_seq = 'zg'
    tsp_fit = 'lorentzian'
    S = LoadSpectra()
    def glucose_fit(x, scale, x_shift, width):
        return Glucose(x_array=x, scale=scale, x_shift=x_shift, width=width)

    concs = ['1250uM', '2500uM', '5000uM', '10000uM']

    SIG_MAIN = []
    #plt.figure()

    for conc in concs:
        S.ReadHDF5(path = "/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/General/DataSource",
                        filename = "Data.hdf5",
                        substance = substance,
                        concentration = conc,
                        field_strength = "80MHz",
                        pulse_seq = pulse_seq)

        SIG = []

        for j in [0,1,2]:
            s_temp = []
            for k in range(0, len(S.initial_ppm[j])):
                #print(conc,j,k)
                x, y = S.initial_ppm[j][k], S.initial_amplitude[j][k]
                A = AnalyseSpectra()
                A.InputData(x=x, y=y)
                A.FitTSP(plot_fit = False, fit_function = tsp_fit)
                A.ScaleSpectra(x_shift = (-1) * A.tsp_centre)
                X = A.processed_ppm
                Y = A.processed_amplitude
                mask = (X >= min(glucose_bounds)) * (X <= max(glucose_bounds))
                X_RED = X[mask]
                Y_RED = Y[mask]
                # Initial guesses
                initial_guesses = [
                    np.nanmax(Y_RED),  # scale
                    0,  # x_shift
                    2 # linewidth
                ]
                # Perform curve fitting
                try:
                    popt, pcov = curve_fit(glucose_fit, X_RED, Y_RED, p0=initial_guesses)

                    # Get optimized parameters
                    opt_scale, opt_x_shift, opt_width = popt

                    # Compute best-fit curve
                    Y_fit = Glucose(x_array=X, scale=opt_scale, x_shift=opt_x_shift, width=opt_width)
                    print(j, k)
                    s_temp.append(sum(Y_fit))
                    fig, axs = plt.subplots(2,1, figsize = [10, 5])
                    axs[0].plot(X, Y, label="Original Data")
                    axs[0].plot(X, Y_fit, label="Fitted Curve", linestyle="--")
                    axs[0].legend()
                    axs[0].set_xlim([-1, 6])
                    msk = (X > 0) *(X < 4)
                    axs[0].set_ylim([np.min(Y[msk]) * 10 , np.max(Y[msk]) * 1.3])
                    axs[1].plot(X, Y-Y_fit, label="Residual")
                    axs[1].legend()
                    axs[1].set_xlim([-1, 6])
                    axs[1].set_ylim([np.min(Y[msk]) * 10 , np.max(Y[msk]) * 1.3])
                    plt.savefig('%s_%s_%s.png' % (conc, j, k))
                    plt.close()

                except:
                    s_temp.append(np.nan)

            SIG.append(s_temp)


        #plt.scatter(scans, np.mean(SIG, axis = 0))
        SIG_MAIN.append(np.mean(SIG, axis = 0))

    path = 'new_data.hdf5'
    u = h5.File(path, 'r+')
    u.create_dataset('glcnew', data = np.array(SIG_MAIN))
    u.close()

    #plt.yscale('log')
    #plt.xscale('log')
    #plt.show()

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
    concs_oth = [concs_oth_lab[0], concs_oth_lab[-1]]
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
    matplotlib.rcParams['font.size'] = 15

    # Create figure
    fig = plt.figure(figsize=(10, 5))

    # Create GridSpec
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 2], height_ratios=[1, 1])

    # Create subplots
    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[1, 0])
    ax3 = plt.subplot(gs[:, 1])

    ax1.plot(X[0], np.zeros(len(X[0])), color = 'grey', alpha = 0.5, ls = '--', lw = 0.5)#, handlelength = 1)
    ax1.plot(X[0], Y[0], color = cols[0], label = '$n_{\mathrm{s}} =  2^{4}$')#, handlelength = 1)

    ax1.plot(X[1], Y[1], color = cols[1], label = '$n_{\mathrm{s}} = 2^{6}$')#, handlelength = 1)
    ax1.plot(X[2], Y[2], color = cols[2], label = '$n_{\mathrm{s}} = 2^{8}$')#, handlelength = 1)
    ax1.annotate('$%s \mathrm{mmol/L}$' % concs_oth[0], xycoords = 'axes fraction', xy=(0.83, 0.9), xytext=(0.83, 0.9), horizontalalignment='center', color = 'k', fontsize = 10)
    ax1.annotate('%s' % substance, xycoords = 'axes fraction', xy=(0.1, 0.9), xytext=(0.1, 0.9), horizontalalignment='center', color = 'k', fontsize = 10)

    a = _getylims(X[0], Y[0], bounds)
    ax1.set_ylim([a[0], a[1] * 1.3])
    ax2.annotate('$%s \mathrm{mmol/L}$' % concs_oth[1], xycoords = 'axes fraction', xy=(0.83, 0.9), xytext=(0.83, 0.9), horizontalalignment='center', color = 'k', fontsize = 10)

    ax2.plot(X[3], np.zeros(len(X[3])), color = 'grey', alpha = 0.5, ls = '--', lw = 0.5)
    ax2.plot(X[3], Y[3], color = cols[0], label = '$n_{\mathrm{s}} =  2^{4}$')
    ax2.plot(X[4], Y[4], color = cols[1], label = '$n_{\mathrm{s}} = 2^{6}$')
    ax2.plot(X[5], Y[5], color = cols[2], label = '$n_{\mathrm{s}} = 2^{8}$')

    a = _getylims(X[3], Y[3], bounds)
    if a[0]>0:
        a[0] = -10
    ax2.set_ylim([a[0], a[1] * 1.3])


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

    #ax1.annotate('$%s \mathrm{mmol/L}$' % concs_oth[0], xy=(1.24, ax1y[1] * 0.87), xytext=(1.24, ax1y[1] * 0.87), horizontalalignment='center', color = 'k')
    #ax2.annotate('$%s \mathrm{mmol/L}$' % concs_oth[1], xy=(1.24, ax2y[1] * 0.9), xytext=(1.24, ax2y[1] * 0.9), horizontalalignment='center', color = 'k')

    ax2.legend(frameon = False, fontsize = 11, handlelength = 1)

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
    
    ax3.plot(np.linspace(-1, 500, len(nscans_long_lab)), np.ones(len(nscans_long_lab)) * 3., color = 'grey', alpha = 0.5, ls = ':')
    ax3.plot(np.linspace(-1, 500, len(nscans_long_lab)), np.ones(len(nscans_long_lab)) * 10., color = 'grey', alpha = 0.5, ls = ':')

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

    ax3.annotate('LOD', 
                xy = (200, 2.6),  # Position of annotation
                xytext=(200, 2.6),  # Offset text a little to the right
                xycoords = 'data',
                fontsize=12,
                color='grey',
                ha='left') 

    ax3.annotate('LOQ', 
                xy = (200, 7.7),  # Position of annotation
                xytext=(200, 7.7),  # Offset text a little to the right
                xycoords = 'data',
                fontsize=12,
                color='grey',
                ha='left') 


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
                    (x_annotate[i], y_annotate[i]*1.3),  # Position of annotation
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
    ax3.set_xlim([0.75, 335])
    # Adjust layout
    plt.tight_layout()

    #plt.title('Substance = %s, Pulse = %s, SNR = %s' % (substance, pulse_seq, snr_method))
    ## Show plot

    fnm = '%s_%s_%s.pdf' % (substance, pulse_seq, snr_method)
    #plt.show()
    plt.savefig(fnm, bbox_inches = 'tight', pad_inches = 0.05)
    plt.close()

def snr_fig_threepan(pulse_seq = 'zg', 
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

    nscans_long = np.arange(0, 9, 1)
    nscans_long_lab = np.array([2**val for val in nscans_long])

    SUBSTANCES = ['Glucose', 'Lactate', 'Citrate']
    cols = ['violet', 'green', 'red']
    matplotlib.rc('text', usetex=True)

    # Set the default font to Times New Roman
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.size'] = 14

    fig, axes = plt.subplots(1, 3, figsize = [15, 5])


    for j in range(0, len(SUBSTANCES)):
        substance = SUBSTANCES[j]
        SIGS, NOISES, SNRs = getSNRs(snr_method = snr_method, plot_spectra = plot_spectra, pulse_seq = pulse_seq, substance = substance, tsp_fit = tsp_fit, plot_snrs = plot_snrs, plot_name = plot_name)

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
        concs_oth = [concs_oth_lab[0], concs_oth_lab[-1]]
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





        ax3 = axes[j]
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
        
        ax3.plot(np.linspace(-1, 500, len(nscans_long_lab)), np.ones(len(nscans_long_lab)) * 1., color = 'grey', alpha = 0.5, ls = ':')
        ax3.plot(np.linspace(-1, 500, len(nscans_long_lab)), np.ones(len(nscans_long_lab)) * 3., color = 'grey', alpha = 0.5, ls = ':')
        ax3.plot(np.linspace(-1, 500, len(nscans_long_lab)), np.ones(len(nscans_long_lab)) * 10., color = 'grey', alpha = 0.5, ls = ':')

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
        snr = SNRs[0]
        ax3.errorbar(nscans_long_lab, MEANS[0], yerr=SEMS[0], fmt='o', label=substance, capsize=5, color='white')
        for k in range(0, len(SNRs)):
            i = len(SNRs) - k -1 
            snr = SNRs[i]
            ax3.errorbar(nscans_long_lab, MEANS[i], yerr=SEMS[i], fmt='o', label='$%s \mathrm{mmol/L}$' % concs_oth_lab[i], capsize=5, color=cols[i], ms = 6)

        if j == 0:
            annotations = ['Best fit to data', '1/2 best fit', '1/4 best fit', '1/8 best fit']
            x_annotate = [15, 15, 15, 15]  # Example x-coordinates to place the text
            for i, (A_factor, text) in enumerate(zip([1, 0.5, 0.25, 0.125], annotations)):
                y_annotate = A * A_factor * np.array(x_annotate) ** slope  # Compute y position based on the line equation
                
                # Add annotation rotated to match the line's slope
                ax3.annotate(text, 
                            (x_annotate[i], y_annotate[i]*1.5),  # Position of annotation
                            xytext=(10, 0),  # Offset text a little to the right
                            textcoords='offset points',
                            fontsize=10,
                            color=np.flip(cols)[i],
                            rotation=np.degrees(np.arctan(slope)),  # Rotate parallel to line
                            ha='left')  # Align left

        #print(np.degrees(np.arctan(slope)))
        ax3.legend(frameon = False, loc = 'upper left', fontsize = 12)
        ax3.set_xscale('log')
        ax3.set_yscale('log')
        ax3.set_xlabel('Number of Scans')
        ax3.set_xlim([0.75, 335])
        ax3.set_ylim([0.8, 301])
        ax3.xaxis.set_tick_params(direction='in', which='both', right=True, top=True)
        ax3.yaxis.set_tick_params(direction='in', which='both', right=True, top=True)

        if j == 2:
            y_values = [3.0, 10.0]
            labels = ['LOD', 'LOQ']
            
            for y_val, label in zip(y_values, labels):
                # Plot the dotted line in two segments with a gap for text
                x_range = np.linspace(-1, 500, 100)
                gap_start = 0.9  # Position where to start the text gap
                gap_end = 2.0    # Position where to end the text gap
                
                ## First segment
                #axes[j].plot(x_range[x_range < gap_start], 
                #            np.ones(len(x_range[x_range < gap_start])) * y_val, 
                #            color='grey', alpha=0.5, ls=':')
                #
                ## Second segment
                #axes[j].plot(x_range[x_range > gap_end], 
                #            np.ones(len(x_range[x_range > gap_end])) * y_val, 
                #            color='grey', alpha=0.5, ls=':')
                
                # Add text with white background
                text = axes[j].text(1.3, y_val, label,
                                fontsize=12,
                                color='grey',
                                horizontalalignment='left',
                                verticalalignment='center',
                                bbox=dict(facecolor='white',
                                        edgecolor='none',
                                        alpha=1,
                                        pad=1))
                
                # Add a thin border around the text for better visibility
                text.set_path_effects([
                    path_effects.Stroke(linewidth=2, foreground='white'),
                    path_effects.Normal()
                ])


    axes[0].set_ylabel('$\mathrm{SNR}$')
    # Adjust layout

    #axes[2].annotate('LOD', 
    #            xy = (1, 3),  # Position of annotation
    #            xytext=(1, 3),  # Offset text a little to the right
    #            xycoords = 'data',
    #            fontsize=12,
    #            color='grey',
    #            ha='left') 


#getSNRs(substance = 'Glucose', pulse_seq = 'zg')
#getSNRs(substance = 'Citrate', pulse_seq = 'zg30')
#getSNRs(substance = 'Citrate', pulse_seq = 'wet')