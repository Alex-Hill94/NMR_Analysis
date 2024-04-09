import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '/Users/alexhill/Documents/GitHub/NMR_Analysis')
from NMR_Tools import *

path = "/Users/alexhill/Desktop/Metabolomics/Data_Analysis/BrukerRerun3Data"
samples = ['D22', 'D23', 'D24']
pulse = 'zg30'
scans = [1, 2, 4, 8, 16, 32, 64, 128, 256]

xlim = None
ylim = None
title = None
save = False
method1 = 'gaussian'
method2 = 'lorentzian'
diff = 0.5
lw_stan =2
method3 = 'gpv'

def fit_water(XDATA, YDATA, method = 'lorentzian', plot = True, save_plot = False):
    def lorentzian(x, amplitude, center, fwhm):
        return amplitude / ((x - center) ** 2 + (fwhm / 2) ** 2)
    def gaussian(x, amplitude, mu, sigma):
        return amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
    def gpv(x, amplitude, center, fwhm, eta):
        return (1-eta) * gaussian(x, amplitude, center, fwhm / np.sqrt(2 * np.log(2))) + eta * lorentzian(x, amplitude, center, fwhm)
    x, y = XDATA, YDATA
    max_y = np.nanmax(y)
    x_where_max_y = x[y == max_y][0]


    amp = max_y
    mean = x_where_max_y
    stdev = 1
    initial_guess = [amp, mean, stdev]
    popt_gauss, pcov = curve_fit(gaussian, x, y, p0=initial_guess)
    output_gauss = gaussian(x, popt_gauss[0], popt_gauss[1], popt_gauss[2])


    amp = max_y
    cent = x_where_max_y
    fwhm = 1
    initial_guess = [amp, cent, fwhm]
    popt_lorentz, pcov = curve_fit(lorentzian, x, y, p0=initial_guess)
    output_lorentz = lorentzian(x, popt_lorentz[0], popt_lorentz[1], popt_lorentz[2])
        

    if plot:
        fig, axs = plt.subplots(1,2 , figsize = [10,5])
        axs[0].set_title('Samp: %s, Nscan = %s' % (sample, nscan))
        axs[0].plot(x, y, color = 'blue', ls = '-', lw = lw_stan, label = 'data_orig')
        axs[0].plot(x, output_gauss, color = 'red', ls = '--', lw = lw_stan/2., label = 'gauss')
        axs[0].plot(x, output_lorentz, color = 'green', ls = '--', lw = lw_stan/2., label = 'lorentz')
        #axs[0].plot(x, output_gpv, color = 'orange', ls = '--', lw = 0.5)
        axs[0].set_xlim([cent - diff, cent + diff])
        axs[0].invert_xaxis()
        axs[0].legend()

        axs[1].plot(x, y-y, color = 'blue', ls = '-', lw = lw_stan, label = 'data - data')
        axs[1].plot(x, y- output_gauss, color = 'red', ls = '--', lw = lw_stan/2.,  label ='data - gauss')
        axs[1].plot(x, y - output_lorentz, color = 'green', ls = '--', lw = lw_stan/2.,  label ='data - lorentz')
        #axs[1].plot(x, y - output_gpv, color = 'orange', ls = '--', lw = 0.5)
        axs[1].set_xlim([cent - diff, cent + diff])
        axs[1].invert_xaxis()
        axs[1].legend()
        if save_plot:
            plt.savefig('%s_%s.png' % (nscan,sample))
        else:
            plt.show()
    if method == 'gaussian':
        return x, output_gauss
    elif method == 'lorentzian':
        return x, output_lorentz


def fit_tsp(XDATA, YDATA, method = 'lorentzian', plot = True, save_plot = False, plot_name = 'dummy.png'):
    def lorentzian(x, amplitude, center, fwhm):
        return amplitude / ((x - center) ** 2 + (fwhm / 2) ** 2)
    def gaussian(x, amplitude, mu, sigma):
        return amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
    def gpv(x, amplitude, center, fwhm, eta):
        return (1-eta) * gaussian(x, amplitude, center, fwhm / np.sqrt(2 * np.log(2))) + eta * lorentzian(x, amplitude, center, fwhm)
    window = 0.2
    tsp_loc = 0.0
    x_mask = (XDATA >= tsp_loc - window) * (XDATA <= tsp_loc + window)
    x, y = XDATA[x_mask], YDATA[x_mask]
    max_y = np.nanmax(y)
    x_where_max_y = x[y == max_y][0]


    amp = max_y
    mean = x_where_max_y
    stdev = 1
    initial_guess = [amp, mean, stdev]
    popt_gauss, pcov = curve_fit(gaussian, x, y, p0=initial_guess, maxfev=5000)
    output_gauss = gaussian(x, popt_gauss[0], popt_gauss[1], popt_gauss[2])


    amp = max_y
    cent = x_where_max_y
    fwhm = 1
    initial_guess = [amp, cent, fwhm]
    popt_lorentz, pcov = curve_fit(lorentzian, x, y, p0=initial_guess, maxfev=5000)
    output_lorentz = lorentzian(x, popt_lorentz[0], popt_lorentz[1], popt_lorentz[2])
    DC = cent
    LC = popt_lorentz[1]
    DI = trapezoidal_rule(np.flip(x), np.flip(y))
    LI = trapezoidal_rule(np.flip(x), np.flip(output_lorentz))
    AMP = popt_lorentz[0]
    if plot:
        fig, axs = plt.subplots(1,2 , figsize = [10,5])
        axs[0].set_title('Samp: %s, Nscan = %s' % (sample, nscan))
        axs[0].plot(x, y, color = 'blue', ls = '-', lw = lw_stan, label = 'data_orig')
        axs[0].plot(x, output_gauss, color = 'red', ls = '--', lw = lw_stan/2., label = 'gauss')
        axs[0].plot(x, output_lorentz, color = 'green', ls = '--', lw = lw_stan/2., label = 'lorentz')
        axs[0].invert_xaxis()
        axs[0].legend()
        axs[1].set_title('DC = %s, LC = %s, DI = %s, LI = %s' % (np.round(DC, 4), np.round(LC, 4), np.round(DI, 0), np.round(LI, 0)))
        axs[1].plot(x, y-y, color = 'blue', ls = '-', lw = lw_stan, label = 'data - data')
        axs[1].plot(x, y- output_gauss, color = 'red', ls = '--', lw = lw_stan/2.,  label ='data - gauss')
        axs[1].plot(x, y - output_lorentz, color = 'green', ls = '--', lw = lw_stan/2.,  label ='data - lorentz')
        axs[1].invert_xaxis()
        axs[1].legend()
    if save_plot:
        plt.savefig(plot_name)
    else:
        plt.show()
    plt.close()

    print('Data centre = %s' %np.round(DC, 4))
    print('Data integral = %s' % np.round(DI, 4))
    print('Lorentzian centre = %s' % np.round(LC, 4))
    print('Lorentzian integral = %s' % np.round(LI, 4))
    return LC, LI, AMP

def get_spectra(nscan = 256):
    sample = 'D24'
    nscan  = nscan
    fnm = "%s_%s_ns%s.txt" % (sample, pulse, nscan)
    x, y = grab_data(path = path, file_name = fnm)
    _, y_lor = fit_water(x, y, plot = False)
    y_fit = y - y_lor
    return [x, y, y_fit, y_lor]

scans = [1, 2, 4, 8, 16, 32, 64, 128, 256]

#scan1 = get_spectra(scans[0])
#scan2 = get_spectra(scans[1])
#scan4 = get_spectra(scans[2])
#scan8 = get_spectra(scans[3])
#scan16 = get_spectra(scans[4])
#scan32 = get_spectra(scans[5])
#scan64 = get_spectra(scans[6])
#scan128 = get_spectra(scans[7])
#scan256 = get_spectra(scans[8])

#for i in range(0, len(scans)):
#    data = get_spectra(scans[i])
#    fig, axs = plt.subplots(1,2, figsize = [10,5])
#    axs[0].set_title('Scans = %s' % scans[i])
#    axs[0].plot(data[0], data[1], color = 'purple', label = 'data orig', linestyle = '-')
#    axs[0].plot(data[0], data[3], color = 'blue', label = 'water fit', linestyle = ':')
#    axs[0].plot(data[0], data[2], color = 'red', label = 'data water subbed', linestyle = '--')
#    axs[0].legend()
#    axs[0].set_ylim([-2.*1e7, 15*1e8])
#    axs[0].set_xlim([0.0, 5.5])
#    axs[0].invert_xaxis()
#    axs[1].set_title('Scans = %s' % scans[i])
#    axs[1].plot(data[0], data[1], color = 'purple', label = 'data orig', linestyle = '-')
#
#    axs[1].plot(data[0], data[2], color = 'blue', label = 'water fit', linestyle = ':')
#    axs[1].plot(data[0], data[3], color = 'red', label = 'data water subbed', linestyle = '--')
#    plt.savefig('scan_%s.png' % scans[i])
#    plt.close()

def get_area(X,Y, avg = True):
    if X[0] > X[-1]:
        X, Y = np.flip(X), abs(np.flip(Y))
    xrange = np.nanmax(X) - np.nanmin(X)
    integral = trapezoidal_rule(X, Y)
    if avg:
        integral = integral/xrange
    return integral

N = []
C = []
L = []

for i in range(0, len(scans)):
    data = get_spectra(scans[i])

    xdata, ydata = data[0], data[2]
    noise_mask = (xdata > -2.0) * (xdata < -1.0)
    cit_mask = (xdata > 2.5) * (xdata < 2.7)
    lac_mask = (xdata > 1.15) * (xdata < 1.5)
    

    noise = get_area(xdata[noise_mask], ydata[noise_mask])
    cit = get_area(xdata[cit_mask], ydata[cit_mask])
    lac = get_area(xdata[lac_mask], ydata[lac_mask])

    N.append(noise)
    C.append(cit)
    L.append(lac)
    fig, axs = plt.subplots(1,1, figsize = [5,5])
    axs.plot(xdata, ydata, lw = 0.5)
    #axs.plot(xdata[noise_mask], ydata[noise_mask], label = 'noise area')
    #axs.plot(xdata[cit_mask], ydata[cit_mask], label = 'citrate area')
    #axs.plot(xdata[lac_mask], ydata[lac_mask], label = 'lactate area')
    axs.invert_xaxis()
    axs.legend()
    plt.savefig('scan_fit_%s.png' % scans[i])
    plt.close()

N = np.array(N)
L = np.array(L)
C = np.array(C)


fig, axs = plt.subplots(1,2, figsize = [10, 5])
axs[0].set_title('Individual measures')
axs[0].plot(scans, N, color = 'k', label = 'Noise')
axs[0].plot(scans, C, color = 'red', label = 'Citrate')
axs[0].plot(scans, L, color = 'blue', label = 'Lactate')
axs[0].legend()
axs[1].set_title('Signal to noise')
axs[1].plot(scans, C/N, color = 'red', label = 'Citrate SNR')
axs[1].plot(scans, L/N, color = 'blue', label = 'Lactate SNR')
axs[0].set_ylabel('Integral')
axs[0].set_xlabel('Nscans')
axs[1].set_xlabel('Nscans')


plt.show()


def f():
    fig, axs = plt.subplots(1,2, figsize = [10, 5])

    axs[0].plot(x_D22, np.zeros(len(x_D22)), color = 'k', lw = 0.2)
    axs[0].plot(x_D22, y_D22, color = 'blue', label = 'D22')
    axs[0].plot(x_D23, y_D23, color = 'orange', label = 'D23')
    axs[0].plot(x_D24, y_D24, color = 'green', label = 'D24')

    axs[0].plot(x_D22, y_D22_lor, color = 'blue', ls = '--')
    axs[0].plot(x_D23, y_D23_lor, color = 'orange', ls = '--')
    axs[0].plot(x_D24, y_D24_lor, color = 'green', ls = '--')
    axs[0].set_xlim([3.0, 5.0])
    axs[0].set_ylim([-2.5 * 1e8, 2.1*1e9])
    axs[0].invert_xaxis()
    axs[0].legend()


    axs[1].plot(x_D22, np.zeros(len(x_D22)), color = 'k', lw = 0.2)
    axs[1].plot(x_D22, D22_fit)
    axs[1].plot(x_D23, D23_fit)
    axs[1].plot(x_D24, D24_fit)
    axs[1].set_xlim([3.0, 5.0])
    axs[1].set_ylim([-2.5 * 1e8, 2.1*1e9])
    axs[1].invert_xaxis()
    plt.show()


    c_22, n_22, amp_22 = fit_tsp(x_D22, D22_fit, plot = False)
    c_23, n_23, amp_23 = fit_tsp(x_D23, D23_fit, plot = False)
    c_24, n_24, amp_24 = fit_tsp(x_D24, D24_fit, plot = False)


    fig, axs = plt.subplots(1,2, figsize = [10,5])
    axs[0].plot(x_D22, D22_fit)
    axs[0].plot(x_D23, D23_fit)
    axs[0].plot(x_D24, D24_fit)

    axs[0].set_xlim([-0.05, 0.05])
    axs[0].invert_xaxis()
    axs[0].set_ylim([-2*1e7, 3*1e8])

    axs[1].plot(x_D22, D22_fit, color = 'blue')
    axs[1].plot(x_D22, D22_fit * (n_22/n_22), color = 'blue', ls = '--')
    axs[1].plot(x_D23 - (c_23 - c_22), D23_fit, ls = '-', color = 'green')
    axs[1].plot(x_D23 - (c_23 - c_22), D23_fit * (n_22/n_23), ls = '--', color = 'green')
    axs[1].plot(x_D24 - (c_24 - c_22), D24_fit, ls = '-', color = 'orange')
    axs[1].plot(x_D24 - (c_24 - c_22), D24_fit * (n_22/n_24), ls = '--', color = 'orange')

    axs[1].set_xlim([-0.05, 0.05])
    axs[1].invert_xaxis()
    axs[1].set_ylim([-2*1e7, 3*1e8])



    plt.show()


    '''

    x_lor, y_lor = fit_water(x, y, plot = False)

    fig, axs = plt.subplots(1,1, figsize = [10, 5])
    axs.plot(x, y )
    axs.plot(x, y - y_lor)
    axs.invert_xaxis()
    plt.show()
    '''
    '''

    for sample in samples:
        for nscan in scans:
            fnm = "%s_%s_ns%s.txt" % (sample, pulse, nscan)
            x, y = grab_data(path = path, file_name = fnm)
            try:
                fit_tsp(x, y, plot = True, save_plot = True, plot_name = 'tsp_%s_%s.png' % (sample, nscan))
            except:
                print('No fit found for %s %s' % (nscan,sample))
    '''
    #plt.close()

    '''

    fig, axs = plt.subplots(1,1, figsize = [5,5])
    axs.set_title('Samp: %s, Nscan = %s' % (sample, nscan))
    axs.plot(x, y, color = 'blue', ls = '-', lw = lw_stan, label = 'data_orig')
    axs.invert_xaxis()
    plt.show()
    '''



    import numpy as np

    def asymmetric_lorentzian(x, A, x0, gamma_left, gamma_right):
        """
        Compute the asymmetric Lorentzian function.

        Parameters:
        x : array-like
            The independent variable.
        A : float
            The amplitude of the peak.
        x0 : float
            The peak position.
        gamma_left : float
            The width on the left side of the peak.
        gamma_right : float
            The width on the right side of the peak.

        Returns:
        array-like
            The value of the asymmetric Lorentzian function at each point in x.
        """
        return A / (1 + ((x - x0) / (gamma_left / 2)) ** 2) / (1 + ((x - x0) / (gamma_right / 2)) ** 2)

    # Example usage
    x = np.linspace(-10, 10, 10000)
    A = 1.0
    x0 = 0.0
    gamma_left = 1.0
    gamma_right = 1.0

    y1 = asymmetric_lorentzian(x, A, x0, gamma_left, gamma_right)
    y2 = asymmetric_lorentzian(x, A, x0, gamma_left, gamma_right*2)
    y3 = asymmetric_lorentzian(x, A, x0, gamma_left, gamma_right*3)
    y4 = asymmetric_lorentzian(x, A, x0, gamma_left, gamma_right*4)

    plt.figure()
    plt.plot(x, y1)
    plt.plot(x, y2)
    plt.plot(x, y3)
    plt.plot(x, y4)
    plt.show()


    def lorentzian_single(x, A, x0, gamma):
        """
        Compute the Lorentzian function.

        Parameters:
        x : array-like
            The independent variable.
        A : float
            The amplitude of the peak.
        x0 : float
            The peak position.
        gamma : float
            The full width at half maximum (FWHM).

        Returns:
        array-like
            The value of the Lorentzian function at each point in x.
        """
        return A / (1 + ((x - x0) / (gamma / 2)) ** 2)

    #def lorentzian(x, amplitude, center, fwhm):
    #    return amplitude / ((x - center) ** 2 + (fwhm / 2) ** 2)



    A, x0, gamma = 1, 1, 1

    y_lorentzian = lorentzian_single(x, A, x0, gamma)
    y_lorentzian1 = lorentzian_single(x, A, x0, gamma)
    plt.figure()
    plt.plot(x, y_lorentzian)
    plt.plot(x, y_lorentzian1)
    plt.show()


    def lorentzian_mod(x, A, x0, gamma_left, gamma_right):
        left_mask = x <= x0
        right_mask = ~left_mask  # Inverse of linear_mask

        #left_output = A / (1 + ((x[left_mask] - x0) / (gamma_left / 2)) ** 2)
        #right_output = A / (1 + ((x[right_mask] - x0) / (gamma_right / 2)) ** 2)

        left_output = lorentzian_single(x[left_mask], A, x0, gamma_left)    
        right_output = lorentzian_single(x[right_mask], A, x0, gamma_right)    

        output = np.empty_like(x)
        output[left_mask] = left_output
        output[right_mask] = right_output
        return output

    x = np.linspace(-10, 10, 10000)
    A = 1.0
    x0 = 0.0
    gamma_left = 1.0
    gamma_right = 1.0


    y1 = lorentzian_mod(x, A, x0, gamma_left, gamma_right)
    y2 = lorentzian_mod(x, A, x0, gamma_left, gamma_right*2)
    y3 = lorentzian_mod(x, A, x0, gamma_left, gamma_right*3)
    y4 = lorentzian_mod(x, A, x0, gamma_left, gamma_right*4)

    plt.figure()
    plt.plot(x, y1)
    plt.plot(x, y2)
    plt.plot(x, y3)
    plt.plot(x, y4)
    plt.show()


    def fit_lor(XDATA, YDATA, fnc = 'asym', plot = True):
        x, y = np.flip(XDATA), np.flip(YDATA)
        max_y = np.nanmax(y)
        x_where_max_y = x[y == max_y][0]


        amp = max_y
        cent = x_where_max_y
        fwhm_left = 1
        fwhm_right = 1
        if fnc == 'asym':
            initial_guess = [amp, cent, fwhm_left, fwhm_right]
            popt, pcov = curve_fit(lorentzian_mod, x, y, p0=initial_guess)
            y_out = lorentzian_mod(x, popt[0], popt[1], popt[2], popt[3])
            
        elif fnc == 'sym':
            initial_guess = [amp, cent, fwhm_left]
            popt, pcov = curve_fit(lorentzian_single, x, y, p0=initial_guess)
            y_out = lorentzian_single(x, popt[0], popt[1], popt[2])
        print(fnc, popt)
        x_out = x
        return x_out, y, y_out



    samples = ['D22', 'D23', 'D24']

    nscan  = 256 

    sample = samples[0]
    fnm = "%s_%s_ns%s.txt" % (sample, pulse, nscan)
    x, y = grab_data(path = path, file_name = fnm)
    x_D22, ydata_D22, yfit_D22_asym = fit_lor(x, y, fnc = 'asym')
    _, _, yfit_D22_sym = fit_lor(x, y, fnc = 'sym')

    sample = samples[1]
    fnm = "%s_%s_ns%s.txt" % (sample, pulse, nscan)
    x, y = grab_data(path = path, file_name = fnm)
    x_D23, ydata_D23, yfit_D23_asym = fit_lor(x, y, fnc = 'asym')
    _, _, yfit_D23_sym = fit_lor(x, y, fnc = 'sym')

    sample = samples[2]
    fnm = "%s_%s_ns%s.txt" % (sample, pulse, nscan)
    x, y = grab_data(path = path, file_name = fnm)
    x_D24, ydata_D24, yfit_D24_asym = fit_lor(x, y, fnc = 'asym')
    _, _, yfit_D24_sym = fit_lor(x, y, fnc = 'sym')

    fig, axs = plt.subplots(1, 2, figsize = [10, 5])

    axs[0].set_title('Asymmetric')
    axs[0].plot(x_D22, ydata_D22, color = 'blue', ls = '-', label = 'D22')
    axs[0].plot(x_D22, yfit_D22_asym, color = 'blue', ls = '--')

    axs[0].plot(x_D23, ydata_D23, color = 'orange', ls = '-', label = 'D23')
    axs[0].plot(x_D23, yfit_D23_asym, color = 'orange', ls = '--')

    axs[0].plot(x_D24, ydata_D24, color = 'green', ls = '-', label = 'D24')
    axs[0].plot(x_D24, yfit_D24_asym, color = 'green', ls = '--')
    axs[0].legend()

    axs[1].set_title('Symmetric')
    axs[1].plot(x_D22, ydata_D22, color = 'blue', ls = '-', label = 'D22')
    axs[1].plot(x_D22, yfit_D22_sym, color = 'blue', ls = '--')

    axs[1].plot(x_D23, ydata_D23, color = 'orange', ls = '-', label = 'D23')
    axs[1].plot(x_D23, yfit_D23_sym, color = 'orange', ls = '--')

    axs[1].plot(x_D24, ydata_D24, color = 'green', ls = '-', label = 'D24')
    axs[1].plot(x_D24, yfit_D24_sym, color = 'green', ls = '--')

    axs[0].set_xlim([3.0, 5.0])
    axs[1].set_xlim([3.0, 5.0])

    axs[0].set_ylim([-0.25 * 1e9, 2.1*1e9])
    axs[1].set_ylim([-0.25 * 1e9, 2.1*1e9])

    axs[0].invert_xaxis()
    axs[1].invert_xaxis()
    plt.show()

    maskD22 = (x_D22>3.0) * (x_D22<5.0)
    maskD23 = (x_D23>3.0) * (x_D23<5.0)
    maskD24 = (x_D24>3.0) * (x_D24<5.0)

    fig, axs = plt.subplots(1, 3, figsize = [12, 4])


    axs[0].plot(x_D22[maskD22], ydata_D22[maskD22], color = 'blue', ls = '-', label = 'Data')
    axs[0].plot(x_D22[maskD22], yfit_D22_asym[maskD22], color = 'blue', ls = '--', label = 'Asym Fit')
    axs[0].plot(x_D22[maskD22], yfit_D22_sym[maskD22], color = 'blue', ls = ':', label = 'Sym Fit')
    axs[0].set_title('D22')
    axs[0].legend()

    axs[1].plot(x_D23[maskD23], ydata_D23[maskD23], color = 'green', ls = '-', label = 'D22')
    axs[1].plot(x_D23[maskD23], yfit_D23_asym[maskD23], color = 'green', ls = '--')
    axs[1].plot(x_D23[maskD23], yfit_D23_sym[maskD23], color = 'green', ls = ':')
    axs[1].set_title('D23')

    axs[2].plot(x_D24[maskD24], ydata_D24[maskD24], color = 'red', ls = '-', label = 'D22')
    axs[2].plot(x_D24[maskD24], yfit_D24_asym[maskD24], color = 'red', ls = '--')
    axs[2].plot(x_D24[maskD24], yfit_D24_sym[maskD24], color = 'red', ls = ':')
    axs[2].set_title('D24')

    my_ax = [axs[0], axs[1], axs[2]]
    for ax in my_ax:
        ax.set_ylim([-0.25 * 1e9, 2.1*1e9])
        ax.invert_xaxis()

    plt.savefig('fits_test.png')



    x, y = np.flip(XDATA), np.flip(YDATA)

    x,y  = x_D24, ydata_D24

    max_y = np.nanmax(y)
    x_where_max_y = x[y == max_y][0]


    amp = max_y
    cent = x_where_max_y
    fwhm_left = 1
    fwhm_right = 1

    initial_guess = [amp, cent, fwhm_left, fwhm_right]
    popt, pcov = curve_fit(lorentzian_mod, x, y, p0=initial_guess)

    ampl, centre, fwhm_l, fwhm_r = popt[0], popt[1], abs(popt[2]), abs(popt[3])
    llim = centre - 3*fwhm_l/2
    rlim = centre + 3*fwhm_r/2


    y1 = lorentzian_mod(x, popt[0], popt[1], popt[2], popt[3] * 1)
    y2 = lorentzian_mod(x, popt[0], popt[1], popt[2], popt[3] * 2)
    y3 = lorentzian_mod(x, popt[0], popt[1], popt[2], popt[3] * 3)


    mask = (x >= llim) * (x <= rlim)
    x_red, y_red = x[mask], y[mask]

    popt_new, pcov_new = curve_fit(lorentzian_mod, x_red, y_red, p0=[popt])
    y_new = lorentzian_mod(x, popt_new[0], popt_new[1], popt_new[2], popt_new[3])



    fig, axs = plt.subplots(1,2, figsize = [10,5])
    axs[0].set_title('D24 - water subtraction')
    axs[0].plot(x, y, label = 'data', color = 'blue')
    axs[0].plot(x, y1, label = 'old fit', color = 'green')
    axs[0].plot(x, y_new, label = 'new fit', color = 'orange')
    axs[0].plot([llim, rlim], [ampl/2, ampl/2], color = 'k')
    axs[0].invert_xaxis()
    axs[0].legend()

    axs[1].plot(x, y - y, color = 'blue', ls = ':')
    axs[1].plot(x, y - y1, color = 'green', ls = ':')
    axs[1].plot(x, y - y_new, color = 'orange', ls = ':')
    axs[1].invert_xaxis()
    plt.show()
    plt.savefig('d22_newfit.png')

    diff_y = np.diff(y)
    diff_x = np.diff(x)

    m = diff_y/diff_x

    print(fnc, popt)
    x_out = x




    '''
    sample = samples[0]
    fnm = "%s_%s_ns%s.txt" % (sample, pulse, nscan)
    XDATA, YDATA = grab_data(path = path, file_name = fnm)


    x, y = np.flip(XDATA), np.flip(YDATA)
    max_y = np.nanmax(y)
    x_where_max_y = x[y == max_y][0]


    amp = max_y
    cent = x_where_max_y
    fwhm_left = 1
    fwhm_right = 1

    initial_guess = [amp, cent, fwhm_left, fwhm_right]
    popt, pcov = curve_fit(lorentzian_mod, x, y, p0=initial_guess)
    y_out_mod = lorentzian_mod(x, popt[0], popt[1], popt[2], popt[3])
        
    initial_guess = [amp, cent, fwhm_left]
    popt, pcov = curve_fit(lorentzian_single, x, y, p0=initial_guess)
    y_out_single = lorentzian_single(x, popt[0], popt[1], popt[2])

    fig, axs = plt.subplots(1,1, figsize = [5,5])
    axs.plot(x, y, label = 'Data')
    axs.plot(x, y_out_mod, label = 'Sym. Lorentzian')
    axs.plot(x, y_out_single, label = 'Asym. Lorentzian')
    axs.legend()
    axs.invert_xaxis()
    plt.show()
    '''
