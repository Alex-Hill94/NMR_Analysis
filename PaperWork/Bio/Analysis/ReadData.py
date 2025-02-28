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
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '/Users/alexhill/Documents/GitHub/NMR_Analysis/Fourier90Runs')
from paul_col import *
from paul_tol_colours import *
parent_dir = os.path.abspath(os.path.join(os.getcwd(), '..'))
sys.path.append(parent_dir)
sys.path.insert(1, '/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/General')
from NMRClassesH5 import *

bio = 1

S = LoadSpectra()
S.ReadTextFile(path = "/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/Bio/Data", filename = "2025-01-15-AH-BIO-0%s_180_zg.txt" % bio)
x1, y1 = S.initial_ppm, S.initial_amplitude

#S.ReadTextFile(path = "/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/Bio/Data", filename = "2025-01-15-AH-BIO-02_180_zg.txt")
#x2, y2 = S.initial_ppm, S.initial_amplitude
#
S.ReadTextFile(path = "/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/Bio/Data", filename = "2025-01-15-AH-BIO-03_180_zg.txt")
x3, y3 = S.initial_ppm, S.initial_amplitude

S.ReadTextFile(path = "/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/Bio/Data", filename = "AHill_Serum_240113_7_%s0_noesygppr1d.txt" % 1)
x_noesy1, y_noesy1 = S.initial_ppm, S.initial_amplitude

S.ReadTextFile(path = "/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/Bio/Data", filename = "AHill_Serum_240113_7_%s0_noesygppr1d.txt" % 3)
x_noesy3, y_noesy3 = S.initial_ppm, S.initial_amplitude


S.ReadTextFile(path = "/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/Bio/Data", filename = "AHill_Serum_240113_7_%s1_cpmgpr1d.txt" % 1)
x_cpmgp1, y_cpmgp1 = S.initial_ppm, S.initial_amplitude

S.ReadTextFile(path = "/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/Bio/Data", filename = "AHill_Serum_240113_7_%s1_cpmgpr1d.txt" % 3)
x_cpmgp3, y_cpmgp3 = S.initial_ppm, S.initial_amplitude


S.ReadTextFile(path = "/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/Bio/Data", filename = "AHill_Serum_240113_7_%s2_ledbpgppr2s1d.txt" % 1)
x_ledb1, y_ledb1 = S.initial_ppm, S.initial_amplitude

S.ReadTextFile(path = "/Users/alexhill/Documents/GitHub/NMR_Analysis/PaperWork/Bio/Data", filename = "AHill_Serum_240113_7_%s2_ledbpgppr2s1d.txt" % 3)
x_ledb3, y_ledb3 = S.initial_ppm, S.initial_amplitude


def plot_ax(axis, xdata, ydata, ylim_range = [-1, 3], Range = [-1, 5], leg_lab = '80MHz'):
    ydata_ylims = ydata[(xdata >= ylim_range[0]) * (xdata <= ylim_range[1])]
    ydata_red = ydata[(xdata >= Range[0]) * (xdata <= Range[1])]
    xdata_red = xdata[(xdata >= Range[0]) * (xdata <= Range[1])]

    Min, Max = np.nanmin(ydata_ylims), np.nanmax(ydata_ylims)
    Span = Max - Min
    ulim = Max + Span*0.1
    llim = Min - Span*0.05
    ylims = [llim, ulim]
    axis.plot(xdata, ydata, color = 'darkslategrey')
    axis.set_ylim(ylims)
    axis.set_xlim(Range)
    axis.invert_xaxis()
    axis.tick_params(axis='x', direction='in', which='both')
    axis.tick_params(axis='y', direction='in', which='both')
    axis.set_xticks([0,1,2,3,4,5])
    axis.set_xticklabels(['']*6)
    yticks = np.linspace(ylims[0], ylims[1], 5)
    axis.set_yticks(yticks)
    axis.set_yticklabels(['']*len(yticks))

    custom_lines = [
        Line2D([0], [0], color='white', marker='s', linestyle='None', markersize=10, label='Square'),
        Line2D([0], [0], color='white', marker='s', linestyle='None', markersize=10, label='Cross')
    ]
    #axis.legend(custom_lines, [leg_lab,''], loc = 'upper right', frameon = False)
    axis.annotate(leg_lab, xy=(0.98, 0.8), xytext=(0.98, 0.8), horizontalalignment='right', color = 'k', xycoords = 'axes fraction')

matplotlib.rc('text', usetex=True)

# Set the default font to Times New Roman
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 14
'''
fig, axs = plt.subplots(4,1,figsize = [8,8])

plot_ax(axs[0], x1, y1, Range = [0., 5.], leg_lab = '80MHz: zg')
plot_ax(axs[1], x_noesy, y_noesy, Range = [0., 5.], leg_lab = '700MHz: noesygp')
plot_ax(axs[2], x_cpmgp, y_cpmgp, Range = [0., 5.], leg_lab = '700MHz: cpmg')
plot_ax(axs[3], x_ledb, y_ledb, Range = [0., 5.], leg_lab = '700MHz: ledbpgp')

axs[3].set_xticklabels([0,1,2,3,4,5])

axs[3].set_xlabel('Chemical shift [ppm]')
fig.supylabel('Relative intensity')
fig.suptitle('Bio Sample \#%s' % bio)

#axs[0].plot(x1, y1, label = 'BIO1')
#axs.plot(x2, y2, label = 'BIO2')
#axs.plot(x3, y3, label = 'BIO3')
#axs[0].legend()
plt.savefig('bio_%s_old.pdf' % bio)
'''
def plot_ax2(axis, xdata, ydata, ylim_range = [-1, 3], Range = [-1, 5], leg_lab = '80MHz', color = 'darkslategrey', extra_yscale = None, extra_yshift = None):
    ydata_ylims = ydata[(xdata >= ylim_range[0]) * (xdata <= ylim_range[1])]
    ydata_red = ydata[(xdata >= Range[0]) * (xdata <= Range[1])]
    xdata_red = xdata[(xdata >= Range[0]) * (xdata <= Range[1])]

    Min, Max = np.nanmin(ydata_ylims), np.nanmax(ydata_ylims)
    ydata = (ydata - Min)/(Max-Min)

    if extra_yscale is not None:
        ydata = extra_yscale * ydata
    if extra_yshift is not None:
        ydata = ydata + extra_yshift

    Span = Max - Min
    ulim = Max + Span*0.1
    llim = Min - Span*0.05
    ylims = [-0.05, 1.1]
    axis.plot(xdata, ydata, color = color)
    axis.set_ylim(ylims)
    axis.set_xlim(Range)
    axis.invert_xaxis()
    axis.tick_params(axis='x', direction='in', which='both')
    axis.tick_params(axis='y', direction='in', which='both')
    axis.set_xticks([0,1,2,3,4,5])
    axis.set_xticklabels(['']*6)
    yticks = np.linspace(ylims[0], ylims[1], 5)
    axis.set_yticks(yticks)
    axis.set_yticklabels(['', '1', '2', '3', '4'])
    #axis.set_yticklabels(['']*len(yticks))

    custom_lines = [
        Line2D([0], [0], color='white', marker='s', linestyle='None', markersize=10, label='Square'),
        Line2D([0], [0], color='white', marker='s', linestyle='None', markersize=10, label='Cross')
    ]
    #axis.legend(custom_lines, [leg_lab,''], loc = 'upper right', frameon = False)


FS = 10

'''
fig, axs = plt.subplots(1,1,figsize = [9,4.5])

plot_ax2(axs, x1, y1, Range = [0., 5.], leg_lab = '80MHz: zg')
plot_ax2(axs, x_cpmgp, y_cpmgp, Range = [0., 5.], leg_lab = '700MHz: cpmg', color = 'tab:blue')
#plot_ax2(axs, x_cpmgp, y_cpmgp, Range = [0., 5.], leg_lab = '700MHz: cpmg', color = 'tab:red')
axs.annotate('80MHz: zg', xy=(0.98, 0.89), xytext=(0.98, 0.89), 
             horizontalalignment='right', color='darkslategrey', xycoords='axes fraction')

axs.annotate('700MHz: cpmg', xy=(0.98, 0.79), xytext=(0.98, 0.79), 
             horizontalalignment='right', color='tab:blue', xycoords='axes fraction')


axs.annotate('glucose', xy=(3.5, .18), xytext=(3.5, .18), 
             horizontalalignment='right', color='k', xycoords='data', fontsize = FS)


axs.annotate('water', xy=(4.2, .5), xytext=(4.2, .5), 
             horizontalalignment='right', color='k', xycoords='data', fontsize = FS)

axs.annotate('-N(CH$_{3}$)$_{3}$', xy=(2.97, .31), xytext=(2.97, .31), 
             horizontalalignment='right', color='k', xycoords='data', fontsize = FS)


axs.annotate('glycoprotein', xy=(1.8, .41), xytext=(1.8, .41), 
             horizontalalignment='right', color='k', xycoords='data', fontsize = FS)


axs.annotate('mobile\n lipid\n -CH$_{3}$', xy=(0.89, .25), xytext=(0.89, .25), 
             horizontalalignment='center', color='k', xycoords='data', fontsize = FS )

axs.annotate('lactate', xy=(1.55, .5), xytext=(1.55, .5), 
             horizontalalignment='center', color='k', xycoords='data', fontsize = FS )

axs.annotate('mobile\n lipid\n (-CH$_{2}$-)$_{n}$', xy=(.96, .57), xytext=(.96, .57), 
             horizontalalignment='center', color='k', xycoords='data', fontsize = FS )



axs.set_xticklabels([0,1,2,3,4,5])

axs.set_xlabel('Chemical shift [ppm]')
axs.set_ylabel('Relative intensity')
axs.set_title('Bio Sample \#%s' % bio)
plt.show()
'''
#plt.savefig('bio_%s_onepan.pdf' % bio)


fig, axes = plt.subplots(2,1,figsize = [7,7])

axs = axes[0]

#plot_ax2(axs, x1, y1, Range = [0., 5.], leg_lab = 'Bio1', color = 'darkslategrey')
plot_ax2(axs, x1, y1, Range = [0., 5.], leg_lab = 'Serum sample A', color = 'tab:blue', extra_yscale = 0.7, extra_yshift = -0.1)
plot_ax2(axs, x3, y3, Range = [0., 5.], leg_lab = 'Serum sample B', color = 'tab:orange')


plot_ax2(axes[1], x_ledb1, y_ledb1, Range = [0., 5.], ylim_range = [-1,2], leg_lab = 'Serum sample A', color = 'tab:blue', extra_yscale = .95)
plot_ax2(axes[1], x_ledb3, y_ledb3, Range = [0., 5.], ylim_range = [-1,2], leg_lab = 'Serum sample B', color = 'tab:orange', extra_yscale = 1.05)


axs.annotate('mobile\n lipid\n -CH$_{3}$', xy=(0.85, .37), xytext=(0.85, .37), 
             horizontalalignment='center', color='k', xycoords='data', fontsize = FS )

axs.annotate('mobile\n lipid\n (-CH$_{2}$-)$_{n}$', xy=(.96, .57), xytext=(.96, .57), 
             horizontalalignment='center', color='k', xycoords='data', fontsize = FS )
axs.annotate('-N(CH$_{3}$)$_{3}$', xy=(2.97, .17), xytext=(2.97, .17), 
             horizontalalignment='right', color='k', xycoords='data', fontsize = FS)


axs.annotate('glycoprotein', xy=(1.8, .17), xytext=(1.8, .17), 
             horizontalalignment='right', color='k', xycoords='data', fontsize = FS)

axs.annotate('water', xy=(3.9, .5), xytext=(3.9, .5), 
             horizontalalignment='right', color='k', xycoords='data', fontsize = FS)

low, high = axs.get_xlim()
yticks = axs.get_yticks()

for y in yticks:
    axes[0].plot([low, high], [y,y], color = 'grey', ls = ':', lw = 0.5)
    axes[1].plot([low, high], [y,y], color = 'grey', ls = ':', lw = 0.5)
axes[0].set_ylabel('Relative intensity')
axes[1].set_ylabel('Relative intensity')
axes[1].set_xlabel('Chemical shift [ppm]')
custom_lines = [
    Line2D([0], [0], color='tab:blue', label='Square'),
    Line2D([0], [0], color='tab:orange', label='Cross')
    ]

axes[0].legend(custom_lines, ['Serum sample A', 'Serum sample B'], frameon = False)
axes[0].set_title('80MHz NMR Experiment: zg')
axes[1].set_title('700MHz NMR Experiment: ledb')
#plt.savefig('showcase_fig_wt.pdf')
#axes[1].set_xticks([0,1,2,3,4,5])

axes[1].set_xticklabels([5, 4,3,2,1, 0])
plt.savefig('700_80_two_pan_ledb.png')


plt.close()
'''

plt.figure()

plt.plot(x_noesy1, y_noesy1)
plt.plot(x_noesy3, y_noesy3)
plt.show()
'''