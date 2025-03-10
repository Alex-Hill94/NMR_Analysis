import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '/Users/alexhill/Documents/GitHub/NMR_Analysis')
sys.path.insert(1, '/Users/alexhill/Documents/GitHub/NMR_Analysis/Fourier90Runs')
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from NMRClasses import *
from scipy.integrate import quad
import math

'''
def citrate():

    path = '/Users/alexhill/Software/ccpnmr3.2.0/citrate_nonoise_80MHz.hdf5'

    u = h5.File(path, 'r')
    xdata_sim = u['80.0/x'][()]
    ydata_sim = u['80.0/toty'][()]
    u.close()

    S = LoadSpectra()
    S.ReadTextFile(nscan = 256, sample = 'D24')
    xdata_data = S.initial_ppm
    ydata_data = S.initial_amplitude


    cdata = np.cumsum(ydata_sim)



    def find_bounds(xdata, ydata, sigma = 1.):
        def gaussian_pdf(u):
            """
            Gaussian probability density function.
            """
            return (1 / math.sqrt(2 * math.pi)) * math.exp(-u**2 / 2)

        # Specify the limits of integration
        lower_limit = -sigma  # -infinity
        upper_limit = sigma   # +infinity

        # Perform numerical integration
        result, error_estimate = quad(gaussian_pdf, lower_limit, upper_limit)
        print(result)
        llim = 0.5 - result/2
        ulim = result/2 + 0.5
        print(llim)
        print(ulim)
        cdata = np.cumsum(ydata)
        cdata = cdata/np.nanmax(cdata)
        x_low = np.min(xdata[cdata < llim])
        x_high = np.max(xdata[cdata > ulim])
        return llim, ulim, x_low, x_high


    # Enable LaTeX rendering
    matplotlib.rc('text', usetex=True)

    # Set the default font to Times New Roman
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.size'] = 12
    fig, axs = plt.subplots(1,1, figsize = [6, 6])
    #axs[1].plot(xdata_sim, cdata/np.nanmax(cdata), 'k')

    color = 'tab:blue'


    sigmas = [1,2,3,4,5]
    sigma_low = []
    sigma_high = []
    for sigma in sigmas:
        llim, ulim, x_low, x_high = find_bounds(xdata_sim, ydata_sim, sigma = sigma)
        axs.axvline(x_low, color='grey', linewidth=1.0, alpha=0.7, linestyle = '--')    
        axs.axvline(x_high, color='grey', linewidth=1.0, alpha=0.7, linestyle = '--')    
        sigma_low.append(x_low)
        sigma_high.append(x_high)


    axs.plot(xdata_sim, ydata_sim / (6.1678*1e8/(8.5*1e7)), color)
    #axs.plot(xdata_sim, np.nanmax(ydata_sim / (6.1678*1e8/(8.5*1e7))) * cdata/np.nanmax(cdata), 'red', lw = 1.0, alpha = 1.0)
    axs.plot(xdata_data, ydata_data, 'grey', lw = 0.5, alpha = 0.5)
    axs.tick_params(axis='y', labelcolor=color)
    axs.set_xlim([2.3,2.9])
    #axs[1].set_xlim([2.3,2.9])

    axs.set_ylim([-0.1*1e8, 1. * 1e8])
    #axs[1].set_ylim([0,1])

    axs.invert_xaxis()
    axs.set_xlabel('$\delta$ [ppm]')
    axs.set_ylabel('$I(\delta)$', color=color)
    ax2 = axs.twinx()


    color = 'tab:red'
    #ax2.set_ylabel('Normalised Cumulative Intensity', color=color)  # we already handled the x-label with ax1
    ax2.set_ylabel('$\sum_{\delta_{\mathrm{max}}}^{\delta_\mathrm{i}} I(\delta) d\delta$', color=color)  # we already handled the x-label with ax1
    #ax2.plot(, 'red', lw = 1.0, alpha = 1.0)
    ax2.plot(xdata_sim, cdata/np.nanmax(cdata), color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    my_axs = [axs, ax2]

    for ax in my_axs:
        ax.tick_params(axis='x', direction='in', which='both')
        ax.tick_params(axis='y', direction='in', which='both')

    llim, ulim = axs.get_ylim()

    yrange = ulim - llim
    y_pos = llim + 0.5* yrange

    x_pos = sigma_low[0] + 0.015
    axs.annotate('1$\sigma$', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = 'grey', rotation = 90)
    x_pos = sigma_high[0] - 0.015
    axs.annotate('1$\sigma$', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = 'grey', rotation = 270)

    x_pos = sigma_low[1] + 0.015
    axs.annotate('2$\sigma$', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = 'grey', rotation = 90)
    x_pos = sigma_high[1] - 0.015
    axs.annotate('2$\sigma$', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = 'grey', rotation = 270)

    #axs.set_title('2$\sigma$ $\in$ [%s, %s]' % (np.round(sigma_low[1], 2), np.round(sigma_high[1], 2)))


    custom_lines = [Line2D([0], [0], color='tab:blue', linestyle='-', lw = 2),
                    Line2D([0], [0], color='grey', linestyle='-',alpha=0.5, lw = 1)]

    axs.legend(custom_lines, ['CCPN Sim.', 'Fourier 80'], loc = 'upper left', frameon = False)
    #axs[1].invert_xaxis()

    plt.savefig('citrate_bounds.pdf')
    plt.close()

def glucose():

    path = '/Users/alexhill/Software/ccpnmr3.2.0/glucose_nonoise_80MHz.hdf5'

    u = h5.File(path, 'r')
    xdata_sim_tot = u['80.0/x'][()]
    ydata_sim_tot = u['80.0/toty'][()]
    u.close()

    S = LoadSpectra()
    S.ReadTextFile(nscan = 256, sample = 'D24')
    xdata_data = S.initial_ppm
    ydata_data = S.initial_amplitude

    def find_bounds(xdata, ydata, sigma = 1.):
        def gaussian_pdf(u):
            """
            Gaussian probability density function.
            """
            return (1 / math.sqrt(2 * math.pi)) * math.exp(-u**2 / 2)

        # Specify the limits of integration
        lower_limit = -sigma  # -infinity
        upper_limit = sigma   # +infinity

        # Perform numerical integration
        result, error_estimate = quad(gaussian_pdf, lower_limit, upper_limit)
        llim = 0.5 - result/2
        ulim = result/2 + 0.5
        cdata = np.cumsum(ydata)
        cdata = cdata/np.nanmax(cdata)
        x_low = np.min(xdata[cdata < llim])
        x_high = np.max(xdata[cdata > ulim])
        return llim, ulim, x_low, x_high

    xbounds_rough = [2.5, 4.3]
    xmask = (xdata_sim_tot >= np.min(xbounds_rough)) * (xdata_sim_tot <= np.max(xbounds_rough))
    xdata_sim = xdata_sim_tot[xmask]
    ydata_sim = ydata_sim_tot[xmask]
    cdata = np.cumsum(ydata_sim)
    xmask_data = (xdata_data >= np.min(xbounds_rough)) * (xdata_data <= np.max(xbounds_rough))

    xdata_data = xdata_data[xmask_data]
    ydata_data = ydata_data[xmask_data]

    # Enable LaTeX rendering
    matplotlib.rc('text', usetex=True)

    # Set the default font to Times New Roman
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.size'] = 12
    fig, axs = plt.subplots(1,1, figsize = [6, 6])
    #axs[1].plot(xdata_sim, cdata/np.nanmax(cdata), 'k')

    color = 'tab:blue'


    sigmas = [1,2,3]
    sigma_low = []
    sigma_high = []
    for sigma in sigmas:
        llim, ulim, x_low, x_high = find_bounds(xdata_sim, ydata_sim, sigma = sigma)
        axs.axvline(x_low, color='grey', linewidth=1.0, alpha=0.7, linestyle = '--')    
        axs.axvline(x_high, color='grey', linewidth=1.0, alpha=0.7, linestyle = '--')    
        sigma_low.append(x_low)
        sigma_high.append(x_high)


    axs.plot(xdata_sim, ydata_sim, color)
    #axs.plot(xdata_sim, np.nanmax(ydata_sim / (6.1678*1e8/(8.5*1e7))) * cdata/np.nanmax(cdata), 'red', lw = 1.0, alpha = 1.0)
    axs.plot(xdata_data, ydata_data/2., 'grey', lw = 0.5, alpha = 0.5)
    axs.tick_params(axis='y', labelcolor=color)
    axs.set_xlim([2.85,4.15])
    #axs[1].set_xlim([2.3,2.9])

    #axs.set_ylim([-0.1*1e8, 1. * 1e8])
    #axs[1].set_ylim([0,1])
    
    axs.invert_xaxis()
    axs.set_xlabel('$\delta$ [ppm]')
    axs.set_ylabel('$I(\delta)$', color=color)
    ax2 = axs.twinx()


    color = 'tab:red'
    #ax2.set_ylabel('Normalised Cumulative Intensity', color=color)  # we already handled the x-label with ax1
    ax2.set_ylabel('$\sum_{\delta_{\mathrm{max}}}^{\delta_\mathrm{i}} I(\delta) d\delta$', color=color)  # we already handled the x-label with ax1
    #ax2.plot(, 'red', lw = 1.0, alpha = 1.0)
    ax2.plot(xdata_sim, cdata/np.nanmax(cdata), color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    my_axs = [axs, ax2]

    for ax in my_axs:
        ax.tick_params(axis='x', direction='in', which='both')
        ax.tick_params(axis='y', direction='in', which='both')

    llim, ulim = axs.get_ylim()

    yrange = ulim - llim
    y_pos = llim + 0.5* yrange

    x_pos = sigma_low[0] + 0.055
    axs.annotate('1$\sigma$', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = 'grey', rotation = 90)
    x_pos = sigma_high[0] - 0.055
    axs.annotate('1$\sigma$', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = 'grey', rotation = 270)

    x_pos = sigma_low[1] + 0.055
    axs.annotate('2$\sigma$', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = 'grey', rotation = 90)
    x_pos = sigma_high[1] - 0.055
    axs.annotate('2$\sigma$', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = 'grey', rotation = 270)

    #axs.set_title('2$\sigma$ $\in$ [%s, %s]' % (np.round(sigma_low[1], 2), np.round(sigma_high[1], 2)))

    custom_lines = [Line2D([0], [0], color='tab:blue', linestyle='-', lw = 2),
                    Line2D([0], [0], color='grey', linestyle='-',alpha=0.5, lw = 1)]

    axs.legend(custom_lines, ['CCPN Sim.', 'Fourier 80'], loc = 'upper left', frameon = False)
    #axs[1].invert_xaxis()

    #plt.show()
    plt.savefig('glucose_bounds.pdf')
    plt.close()

def lactate():

    path = '/Users/alexhill/Software/ccpnmr3.2.0/lactate31_nonoise_80MHz.hdf5'

    u = h5.File(path, 'r')
    xdata_sim_tot = u['80.0/x'][()]
    ydata_sim_tot = u['80.0/toty'][()]
    u.close()

    S = LoadSpectra()
    S.ReadTextFile(nscan = 256, sample = 'D24')
    xdata_data = S.initial_ppm
    ydata_data = S.initial_amplitude



    def find_bounds(xdata, ydata, sigma = 1.):
        def gaussian_pdf(u):
            """
            Gaussian probability density function.
            """
            return (1 / math.sqrt(2 * math.pi)) * math.exp(-u**2 / 2)

        # Specify the limits of integration
        lower_limit = -sigma  # -infinity
        upper_limit = sigma   # +infinity

        # Perform numerical integration
        result, error_estimate = quad(gaussian_pdf, lower_limit, upper_limit)
        print(result)
        llim = 0.5 - result/2
        ulim = result/2 + 0.5
        print(llim)
        print(ulim)
        cdata = np.cumsum(ydata)
        cdata = cdata/np.nanmax(cdata)
        x_low = np.min(xdata[cdata < llim])
        x_high = np.max(xdata[cdata > ulim])
        return llim, ulim, x_low, x_high

    xbounds_rough = [0.5, 2.5]
    xmask = (xdata_sim_tot >= np.min(xbounds_rough)) * (xdata_sim_tot <= np.max(xbounds_rough))
    xmask_data = (xdata_data >= np.min(xbounds_rough)) * (xdata_data <= np.max(xbounds_rough))

    xdata_data = xdata_data[xmask_data]
    ydata_data = ydata_data[xmask_data]

    xdata_sim = xdata_sim_tot[xmask]
    ydata_sim = ydata_sim_tot[xmask]
    cdata = np.cumsum(ydata_sim)
    # Enable LaTeX rendering
    matplotlib.rc('text', usetex=True)

    # Set the default font to Times New Roman
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.size'] = 12
    fig, axs = plt.subplots(1,1, figsize = [6, 6])
    #axs[1].plot(xdata_sim, cdata/np.nanmax(cdata), 'k')

    color = 'tab:blue'


    sigmas = [1,2,3]
    sigma_low = []
    sigma_high = []
    for sigma in sigmas:
        llim, ulim, x_low, x_high = find_bounds(xdata_sim, ydata_sim, sigma = sigma)
        axs.axvline(x_low, color='grey', linewidth=1.0, alpha=0.7, linestyle = '--')    
        axs.axvline(x_high, color='grey', linewidth=1.0, alpha=0.7, linestyle = '--')    
        sigma_low.append(x_low)
        sigma_high.append(x_high)


    axs.plot(xdata_sim, ydata_sim*1.2, color)
    #axs.plot(xdata_sim, np.nanmax(ydata_sim / (6.1678*1e8/(8.5*1e7))) * cdata/np.nanmax(cdata), 'red', lw = 1.0, alpha = 1.0)
    axs.plot(xdata_data, ydata_data, 'grey', lw = 0.5, alpha = 0.5)
    axs.tick_params(axis='y', labelcolor=color)
    axs.set_xlim([1.1,1.6])
    #axs[1].set_xlim([2.3,2.9])

    #axs.set_ylim([-0.1*1e8, 1. * 1e8])
    #axs[1].set_ylim([0,1])
    
    axs.invert_xaxis()
    axs.set_xlabel('$\delta$ [ppm]')
    axs.set_ylabel('$I(\delta)$', color=color)
    ax2 = axs.twinx()


    color = 'tab:red'
    #ax2.set_ylabel('Normalised Cumulative Intensity', color=color)  # we already handled the x-label with ax1
    ax2.set_ylabel('$\sum_{\delta_{\mathrm{max}}}^{\delta_\mathrm{i}} I(\delta) d\delta$', color=color)  # we already handled the x-label with ax1
    #ax2.plot(, 'red', lw = 1.0, alpha = 1.0)
    ax2.plot(xdata_sim, cdata/np.nanmax(cdata), color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    my_axs = [axs, ax2]

    for ax in my_axs:
        ax.tick_params(axis='x', direction='in', which='both')
        ax.tick_params(axis='y', direction='in', which='both')

    llim, ulim = axs.get_ylim()

    yrange = ulim - llim
    y_pos = llim + 0.5* yrange

    x_pos = sigma_low[0] + 0.055
    axs.annotate('1$\sigma$', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = 'grey', rotation = 90)
    x_pos = sigma_high[0] - 0.055
    axs.annotate('1$\sigma$', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = 'grey', rotation = 270)

    x_pos = sigma_low[1] + 0.055
    axs.annotate('2$\sigma$', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = 'grey', rotation = 90)
    x_pos = sigma_high[1] - 0.055
    axs.annotate('2$\sigma$', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = 'grey', rotation = 270)

    #axs.set_title('2$\sigma$ $\in$ [%s, %s]' % (np.round(sigma_low[1], 2), np.round(sigma_high[1], 2)))

    custom_lines = [Line2D([0], [0], color='tab:blue', linestyle='-', lw = 2),
                    Line2D([0], [0], color='grey', linestyle='-',alpha=0.5, lw = 1)]

    axs.legend(custom_lines, ['CCPN Sim.', 'Fourier 80'], loc = 'upper left', frameon = False)
    #axs[1].invert_xaxis()

    plt.savefig('lactate_bounds.pdf')
    #plt.show()
    plt.close()

def all_spec():

    path = '/Users/alexhill/Software/ccpnmr3.2.0/lactate31_nonoise_80MHz.hdf5'

    u = h5.File(path, 'r')
    xdata_sim_lac = u['80.0/x'][()]
    ydata_sim_lac = u['80.0/toty'][()]
    u.close()

    path = '/Users/alexhill/Software/ccpnmr3.2.0/glucose_nonoise_80MHz.hdf5'

    u = h5.File(path, 'r')
    xdata_sim_glc = u['80.0/x'][()]
    ydata_sim_glc = u['80.0/toty'][()]
    u.close()

    path = '/Users/alexhill/Software/ccpnmr3.2.0/citrate_nonoise_80MHz.hdf5'

    u = h5.File(path, 'r')
    xdata_sim_cit = u['80.0/x'][()]
    ydata_sim_cit = u['80.0/toty'][()]
    u.close()

    S = LoadSpectra()
    S.ReadTextFile(nscan = 256, sample = 'D24')
    S.SubtractWater()
    xdata_data = S.initial_ppm
    ydata_data = S.initial_amplitude


    #cdata = np.cumsum(ydata_sim)

    custom_line = [Line2D([0], [0], color='white', linestyle='-', alpha = 0.0)]

    def find_bounds(xdata, ydata, sigma = 1.):
        def gaussian_pdf(u):
            """
            Gaussian probability density function.
            """
            return (1 / math.sqrt(2 * math.pi)) * math.exp(-u**2 / 2)

        # Specify the limits of integration
        lower_limit = -sigma  # -infinity
        upper_limit = sigma   # +infinity

        # Perform numerical integration
        result, error_estimate = quad(gaussian_pdf, lower_limit, upper_limit)
        print(result)
        llim = 0.5 - result/2
        ulim = result/2 + 0.5
        print(llim)
        print(ulim)
        cdata = np.cumsum(ydata)
        cdata = cdata/np.nanmax(cdata)
        x_low = np.min(xdata[cdata < llim])
        x_high = np.max(xdata[cdata > ulim])
        return llim, ulim, x_low, x_high

    matplotlib.rc('text', usetex=True)

    # Set the default font to Times New Roman
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.size'] = 12


    fig, axs = plt.subplots(5,1, figsize = [10, 6])
    axs[0].plot(xdata_data, ydata_data, 'grey')
    #axs[1].plot(xdata_sim_glc, ydata_sim_glc)
    axs[2].plot(xdata_sim_glc, ydata_sim_glc, color = 'teal')
    axs[3].plot(xdata_sim_lac, ydata_sim_lac, color = 'mediumspringgreen')
    axs[4].plot(xdata_sim_cit, ydata_sim_cit, color = 'limegreen')
    
    glc_height = 1.211*1e9
    lac_height = 6.29*1e8
    cit_height = 8.46*1e7

    scale_glc = glc_height/glc_height
    scale_lac = lac_height/glc_height
    scale_cit = cit_height/glc_height

    ydata_comb = (ydata_sim_glc * scale_glc) + (ydata_sim_lac * scale_lac) + (ydata_sim_cit * scale_cit)

    axs[1].plot(xdata_sim_glc, ydata_comb, color = 'dodgerblue')


    my_axs = [axs[0], axs[1], axs[2], axs[3], axs[4]]

    axs[0].set_ylim([0, 1.35*1e9])

    for ax in my_axs:
        ax.tick_params(axis='x', direction='in', which='both')
        ax.tick_params(axis='y', direction='in', which='both')
        ax.set_yticklabels([])        
        ax.set_xlim([0.3, 5.6])
        ax.invert_xaxis()

    for ax in my_axs[:-1]:
        ax.set_xticklabels([])

    axs[0].legend(custom_line, ['Experimental Data'], frameon = False, loc = 'upper right')    
    axs[1].legend(custom_line, ['Combined Sim.'], frameon = False, loc = 'upper right')
    axs[2].legend(custom_line, ['Glucose Sim.'], frameon = False, loc = 'upper right')
    axs[3].legend(custom_line, ['Lactate Sim.'], frameon = False, loc = 'upper right')
    axs[4].legend(custom_line, ['Citrate Sim.'], frameon = False, loc = 'upper right')

    axs[2].set_ylabel('$I(\delta)$')
    axs[4].set_xlabel('$\delta$ [ppm]')

    plt.subplots_adjust(hspace = 0)
    #plt.show()
    plt.savefig('comp_sims.pdf')
    plt.close()

    #fig, axs = plt.subplots(1,1, figsize = [4, 6])
    #axs.plot(xdata_data, ydata_data, 'grey')
    #plt.show()
'''

path = '/Users/alexhill/Documents/UOL/Research/Companies/ViBo/TemplateSpectra1.hdf5'
#metabs = ['Glucose', 'Lactate', 'Citrate']
freqs = np.linspace(80, 700.0, 63)
linewidth = 1.2
field = 80.0

ppms = []
amps = []

u = h5.File(path, 'r')

metab = 'Glucose'
X = u['%s/%s/ppm' % (metab, field)][()]
Y = u['%s/%s/amplitude' % (metab, field)][()]
ppms.append(X)
amps.append(Y)

X = u['%s/%s/ppm' % (metab, 700.0)][()]
Y = u['%s/%s/amplitude' % (metab, 700.0)][()]
ppms.append(X)
amps.append(Y)

u.close()

L = LoadSpectra()
L.ReadTextFile(path = './',filename = '10C_Glucose_zg_256_700MHz.txt')
ppm_700, amp_700 = L.initial_ppm, L.initial_amplitude

L.ReadTextFile(path = './',filename = '10C_Glucose_zg_256_80MHz.txt')
ppm_80, amp_80 = L.initial_ppm, L.initial_amplitude


glc_height = 1.211*1e9
lac_height = 6.29*1e8
cit_height = 8.46*1e7

glucose_colour = 'orange'
citrate_colour = 'green'
lactate_colour = 'blue'

colours = [glucose_colour, citrate_colour, lactate_colour]

scale_glc = glc_height/glc_height
scale_lac = lac_height/glc_height
scale_cit = cit_height/glc_height

scales = [scale_glc, scale_lac, scale_cit]


matplotlib.rc('text', usetex=True)

# Set the default font to Times New Roman
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 12
fig, axs = plt.subplots(2,1, figsize = [12, 6])

my_axs = [axs[0], axs[1]]

exp_val = abs(np.nanmin(amp_80))
sim_val = abs(np.nanmin(amps[0]))
scale = exp_val/sim_val

exp_X_80, exp_Y_80 = ppm_80, amp_80
sim_X_80, sim_Y_80 = ppms[0] + 0.13, 1.7*amps[0]*scale/1000
lims_80 = [np.nanmin(sim_Y_80[(sim_X_80>3) * (sim_X_80<4)]), np.nanmax(sim_Y_80[(sim_X_80>3) * (sim_X_80<4)])]

my_axs[0].plot(exp_X_80, exp_Y_80, color = 'grey', lw = linewidth, alpha = 0.8)
my_axs[0].plot(sim_X_80, sim_Y_80, color = colours[0], lw = linewidth)


m700_lims = [4.63, 4.67]
exp_val = np.nanmean(amp_700[(ppm_700>4.63) * (ppm_700<4.67)])
sim_val = np.nanmean(amps[1][(ppms[1]>4.63) * (ppms[1]<4.67)])
scale = exp_val/sim_val

exp_X_700, exp_Y_700 = ppm_700, amp_700
sim_X_700, sim_Y_700 = ppms[1], amps[1]*scale
lims_700 = [np.nanmin(sim_Y_700[(sim_X_700>3) * (sim_X_700<4)]), np.nanmax(sim_Y_700[(sim_X_700>3) * (sim_X_700<4)])]

my_axs[1].plot(exp_X_700, exp_Y_700, color = 'grey', lw = linewidth, alpha = 0.8)
my_axs[1].plot(sim_X_700, sim_Y_700, color = colours[0], lw = linewidth)


for my_ax in my_axs:
    my_ax.set_xlim([2.8, 5.6])
    my_ax.set_yticklabels([])
    my_ax.invert_xaxis()
    my_ax.tick_params(axis='x', direction='in', which='both')
    my_ax.tick_params(axis='y', direction='in', which='both')

for my_ax in my_axs[:-1]:
    my_ax.set_xticklabels([])

custom_line = [Line2D([0], [0], color='white', linestyle='-', alpha = 0.0),
                Line2D([0], [0], color='grey', linestyle='-', alpha = 1.0),
                Line2D([0], [0], color=colours[0], linestyle='-', alpha = 1.0)]

cl = Line2D([0], [0], color='white', linestyle='-', alpha = 0.0),

axs[0].legend(custom_line, ['80 MHz', 'Exp. (zg: ns = 256)', 'Sim. (CCPN)'], frameon = False, loc = 'upper right', labelcolor = 'k')
axs[1].legend(cl, ['700 MHz'], frameon = False, loc = 'upper right', labelcolor = 'k')
#axs[1].set_ylabel('$I(\delta)$')

axs[1].set_ylabel('Rel. Intensity')
axs[0].set_ylabel('Rel. Intensity')


plt.subplots_adjust(hspace = 0.)
my_axs[0].set_ylim([-1.*lims_80[0], lims_80[1] * 1.2])
my_axs[1].set_ylim([-1.*lims_700[0], lims_700[1] * 1.2])

#plt.show()
plt.savefig('demo.png')
plt.close()