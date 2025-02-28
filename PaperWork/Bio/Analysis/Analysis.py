import sys
import os 
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import h5py as h5
from scipy import stats
from scipy.stats import chi2
import argparse

# Import custom modules
sys.path.insert(1, '/Users/alexhill/Documents/GitHub/NMR_Analysis/Fourier90Runs')
from paul_col import *
from paul_tol_colours import *
parent_dir = os.path.abspath(os.path.join(os.getcwd(), '..'))
sys.path.append(parent_dir)
sys.path.insert(1, '/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/General')
from NMRClassesH5 import *

def main(bio):
    S = LoadSpectra()
    S.ReadTextFile(path="/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/Bio/Data", 
                   filename=f"2025-01-15-AH-BIO-0{bio}_180_zg.txt")
    x1, y1 = S.initial_ppm, S.initial_amplitude

    S.ReadTextFile(path="/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/Bio/Data", 
                   filename=f"AHill_Serum_240113_7_{bio}0_noesygppr1d.txt")
    x_noesy, y_noesy = S.initial_ppm, S.initial_amplitude

    S.ReadTextFile(path="/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/Bio/Data", 
                   filename=f"AHill_Serum_240113_7_{bio}1_cpmgpr1d.txt")
    x_cpmgp, y_cpmgp = S.initial_ppm, S.initial_amplitude

    S.ReadTextFile(path="/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/Bio/Data", 
                   filename=f"AHill_Serum_240113_7_{bio}2_ledbpgppr2s1d.txt")
    x_ledb, y_ledb = S.initial_ppm, S.initial_amplitude

    def plot_ax(axis, xdata, ydata, ylim_range=[-1, 3], Range=[-1, 5], leg_lab='80MHz'):
        ydata_ylims = ydata[(xdata >= ylim_range[0]) * (xdata <= ylim_range[1])]
        ydata_red = ydata[(xdata >= Range[0]) * (xdata <= Range[1])]
        xdata_red = xdata[(xdata >= Range[0]) * (xdata <= Range[1])]

        Min, Max = np.nanmin(ydata_ylims), np.nanmax(ydata_ylims)
        Span = Max - Min
        ulim = Max + Span * 0.1
        llim = Min - Span * 0.05
        ylims = [llim, ulim]
        axis.plot(xdata, ydata, color='darkslategrey', linewidth = 0.5)
        axis.set_ylim(ylims)
        axis.set_xlim(Range)
        axis.invert_xaxis()
        axis.tick_params(axis='x', direction='in', which='both')
        axis.tick_params(axis='y', direction='in', which='both')
        axis.set_xticks([0, 1, 2, 3, 4, 5])
        axis.set_xticklabels([''] * 6)
        yticks = np.linspace(ylims[0], ylims[1], 5)
        axis.set_yticks(yticks)
        axis.set_yticklabels([''] * len(yticks))

        axis.annotate(leg_lab, xy=(0.98, 0.8), xytext=(0.98, 0.8), 
                      horizontalalignment='right', color='k', xycoords='axes fraction')

    mpl.rc('text', usetex=True)
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['font.size'] = 14

    fig, axs = plt.subplots(4, 1, figsize=[8, 8])

    plot_ax(axs[0], x1, y1, Range=[0., 5.], leg_lab='80MHz: zg')
    plot_ax(axs[1], x_noesy, y_noesy, Range=[0., 5.], leg_lab='700MHz: noesygp')
    plot_ax(axs[2], x_cpmgp, y_cpmgp, Range=[0., 5.], leg_lab='700MHz: cpmg')
    plot_ax(axs[3], x_ledb, y_ledb, Range=[0., 5.], leg_lab='700MHz: ledbpgp')

    axs[3].set_xticklabels([0, 1, 2, 3, 4, 5])
    axs[3].set_xlabel('Chemical shift [ppm]')
    fig.supylabel('Relative intensity')
    fig.suptitle('Bio Sample \#%s' % bio)

    plt.savefig(f'bio_{bio}.pdf')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Specify the bio variable via the command line.")
    parser.add_argument("--bio", type=int, required=True, help="The bio sample number to analyze.")
    args = parser.parse_args()

    main(args.bio)

