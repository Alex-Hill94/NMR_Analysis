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
    ax2.set_ylabel('$\int I(\delta) d\delta$', color=color)  # we already handled the x-label with ax1
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


    custom_lines = [Line2D([0], [0], color='tab:blue', linestyle='-', lw = 2),
                    Line2D([0], [0], color='grey', linestyle='-',alpha=0.5, lw = 1)]

    axs.legend(custom_lines, ['CCPN Sim.', 'Fourier 80'], loc = 'upper left', frameon = False)
    #axs[1].invert_xaxis()

    plt.savefig('citrate_bounds.pdf')
    '''




glucose()
lactate()
citrate()

#all_spec()