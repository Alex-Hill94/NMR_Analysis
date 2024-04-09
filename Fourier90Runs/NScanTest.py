import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '/Users/alexhill/Documents/GitHub/NMR_Analysis')
from NMR_Tools import *



#x_values_1, data_1 = grab_data(path = "/Users/alexhill/Desktop/Metabolomics/Data_Analysis/BrukerRerun3Data", file_name = "D22_zg30_ns1.txt")

def comp_scan(  path = "/Users/alexhill/Desktop/Metabolomics/Data_Analysis/BrukerRerun3Data", 
                sample = 'D22', 
                pulse = 'zg30', 
                scans = [1, 2, 4],
                xlim = None,
                ylim = None,
                title = None,
                save = False,
                plot_name = 'dummy.pdf'):
    len_data = len(scans)

    XDATA = []
    YDATA = []

    for i in range(0, len_data):
        nscan = str(int(scans[i]))
        fnm = "%s_%s_ns%s.txt" % (sample, pulse, nscan)
        x, y = grab_data(path = path, file_name = fnm)
        XDATA.append(x)
        YDATA.append(y)

    fig, axs = plt.subplots(1,1, figsize = [12,6])

    for i in range(0, len_data):
        nscan = scans[i]
        x, y = XDATA[i], YDATA[i]
        axs.plot(x, y, label = 'nscan = %s' % nscan)
    axs.legend()
    if xlim is not None:
        axs.set_xlim(xlim)
    if ylim is not None:
        axs.set_ylim(ylim)
    if title is not None:
        axs.set_title(title)
    axs.invert_xaxis()
    axs.set_xlabel('ppm')
    axs.set_ylabel('amplitude')
    if save:
        plt.savefig(plot_name)
    else:
        plt.show()

    plt.close()

def comp_concs( path = "/Users/alexhill/Desktop/Metabolomics/Data_Analysis/BrukerRerun3Data", 
                samples = ['D22', 'D23', 'D24'], 
                pulse = 'zg30', 
                nscan = 64,
                xlim = None,
                ylim = None,
                title = None,
                save = False,
                plot_name = 'dummy.pdf'):

    len_data = len(samples)

    XDATA = []
    YDATA = []

    for i in range(0, len_data):
        sample = samples[i]
        fnm = "%s_%s_ns%s.txt" % (sample, pulse, nscan)
        x, y = grab_data(path = path, file_name = fnm)
        XDATA.append(x)
        YDATA.append(y)

    fig, axs = plt.subplots(1,1, figsize = [12,6])

    for i in range(0, len_data):
        sample = samples[i]
        x, y = XDATA[i], YDATA[i]
        axs.plot(x, y, label = 'Samp = %s' % sample)
    axs.legend()
    if xlim is not None:
        axs.set_xlim(xlim)
    if ylim is not None:
        axs.set_ylim(ylim)
    if title is not None:
        axs.set_title(title)
    axs.invert_xaxis()
    axs.set_xlabel('ppm')
    axs.set_ylabel('amplitude')
    if save:
        plt.savefig(plot_name)
    else:
        plt.show()

    plt.close()



def comp_area(  path = "/Users/alexhill/Desktop/Metabolomics/Data_Analysis/BrukerRerun3Data", 
                sample = 'D22', 
                pulse = 'zg30', 
                scans = [1, 2, 4],
                xrange = None,
                xlim = None,
                ylim = None,
                title = None,
                save = False,
                plot_name = 'dummy.pdf'):
    len_data = len(scans)

    XDATA = []
    YDATA = []

    for i in range(0, len_data):
        nscan = str(int(scans[i]))
        fnm = "%s_%s_ns%s.txt" % (sample, pulse, nscan)
        x, y = grab_data(path = path, file_name = fnm)
        XDATA.append(x)
        YDATA.append(y)
        print(trap(x,y, xrange = [4.65, 5])[-1])

    #x, y = XDATA[0], YDATA[0]


#comp_scan(scans = [256, 128, 64, 32], xlim = [-0.1, 5.1], ylim = [-10*1e6, 100*1e7], title = 'D22: 0.2mM Lac, 2mM Cit, 5mM Glc', save = True, plot_name = 'D22_zoom_view_hscans.pdf')
#comp_scan(scans = [16, 8, 4, 2, 1], xlim = [-0.1, 5.1], ylim = [-5*1e6, 5*1e7], title = 'D22: 0.2mM Lac, 2mM Cit, 5mM Glc', save = True, plot_name = 'D22_zoom_view_lscans.pdf')


#comp_scan(scans = [256, 128, 64, 32], xlim = [-0.1, 5.1], title = 'D22: 0.2mM Lac, 2mM Cit, 5mM Glc', save = True, plot_name = 'D22_tot_view_hscans.pdf')
comp_concs(path = "/Users/alexhill/Desktop/Metabolomics/Data_Analysis/BrukerRerun3Data", nscan = 64)
comp_concs(path = "/Users/alexhill/Desktop/Metabolomics/Data_Analysis/BrukerInPersonData", nscan = 64)

#comp_area(scans = [16, 8, 4, 2, 1], xlim = [-0.1, 5.1], title = 'D22: 0.2mM Lac, 2mM Cit, 5mM Glc', save = False, plot_name = 'D22_tot_view_lscans.pdf')


'''




fig, axs = plt.subplots(1,1, figsize = [12,6])
axs.plot(x_values_14, data_14, label = 'ID: 14. NS: 256')
axs.plot(x_values_140, data_140, label = 'ID: 140. NS: 64')
axs.plot(x_values_141, data_141, label = 'ID: 141. NS: 128')
axs.plot(x_values_142, data_142, label = 'ID: 142. NS: 192')
axs.legend()
axs.invert_xaxis()
axs.set_xlabel('ppm')
axs.set_ylabel('amplitude')

plt.show()

x14, y14, tot14 = trap(x_values_14, data_14, xrange = [4.4, 5.0])
x140, y140, tot140 = trap(x_values_140, data_140, xrange = [4.4, 5.0])
x141, y141, tot141 = trap(x_values_141, data_141, xrange = [4.4, 5.0])
x142, y142, tot142 = trap(x_values_142, data_142, xrange = [4.4, 5.0])

fig, axs = plt.subplots(1,2, figsize = [10, 5])

axs[0].plot(x140,  y140, label = 'nscan %s, tot = %s' % (64, np.round(np.log10(tot140), 3)), color = 'red')
axs[0].plot(x141,  y141, label = 'nscan %s, tot = %s' % (128, np.round(np.log10(tot141), 3)), color = 'k')
axs[0].plot(x142,  y142, label = 'nscan %s, tot = %s' % (192, np.round(np.log10(tot142), 3)), color = 'green')
axs[0].plot(x14, y14, label = 'nscan %s, tot = %s' % (256, np.round(np.log10(tot14), 3)), color = 'blue')
axs[0].legend()

axs[1].scatter([64], [tot140], marker = 'x', color = 'red')

axs[1].scatter([128], [tot141], marker = 'x', color = 'k')

axs[1].scatter([192], [tot142], marker = 'x', color = 'green')

axs[1].scatter([256], [tot14], marker = 'x', color = 'blue')

axs[0].set_xlabel('ppm')
axs[0].set_ylabel('amplitude')

axs[1].set_ylabel('Integrated area')
axs[1].set_xlabel('Number of scans')
axs[0].invert_xaxis()

plt.savefig('nscans_vs_amp.png')

'''