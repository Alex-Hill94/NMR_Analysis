import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '/Users/alexhill/Documents/GitHub/NMR_Analysis')
from NMR_Tools import *

# Code to compare concentration SNRs
# Fundamentally, we want to first fit to TSP, getting the area, peak height, and centre

def get_spectra(sample = 'D24', nscan = 256, pulse = "zg30", path = "/Users/alexhill/Desktop/Metabolomics/Data_Analysis/BrukerRerun3Data"):
    sample = sample
    nscan  = nscan
    fnm = "%s_%s_ns%s.txt" % (sample, pulse, nscan)
    print('Grabbing %s' % fnm)
    x, y = grab_data(path = path, file_name = fnm)
    #_, y_lor = fit_water(x, y, plot = False)
    #y_fit = y - y_lor
    return [x, y]

def get_area(X,Y, avg = True):
    if X[0] > X[-1]:
        X, Y = np.flip(X), abs(np.flip(Y))
    xrange = np.nanmax(X) - np.nanmin(X)
    integral = trapezoidal_rule(X, Y)
    if avg:
        integral = integral/xrange
    return integral

def _lorentzian_single(x, A, x0, gamma):
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

def _lorentzian_asymmetric(x, A, x0, gamma_left, gamma_right):
    left_mask = x <= x0
    right_mask = ~left_mask  # Inverse of linear_mask

    #left_output = A / (1 + ((x[left_mask] - x0) / (gamma_left / 2)) ** 2)
    #right_output = A / (1 + ((x[right_mask] - x0) / (gamma_right / 2)) ** 2)

    left_output = _lorentzian_single(x[left_mask], A, x0, gamma_left)    
    right_output = _lorentzian_single(x[right_mask], A, x0, gamma_right)    

    output = np.empty_like(x)
    output[left_mask] = left_output
    output[right_mask] = right_output
    return output

def _gaussian(x, amplitude, mu, sigma):
    return amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

def fit_singlet(XDATA, 
                YDATA, 
                method = 'lorentzian_asymmetric', 
                window = 0.2, 
                expected_centre = 0.0, 
                show_fit = True, 
                save = False, 
                fnm = None,
                title = None, 
                return_uncertainty = True):

    window = window
    peak_loc = expected_centre
    x_mask = (XDATA >= peak_loc - window) * (XDATA <= peak_loc + window)
    x, y = XDATA[x_mask], YDATA[x_mask]
    max_y = np.nanmax(y)
    x_where_max_y = x[y == max_y][0]

    amp = max_y
    centre = x_where_max_y
    fwhm = 0.5
    
    if method == 'lorentzian_asymmetric':
        intital_guess = [amp, centre, fwhm, fwhm]
        popt, pcov = curve_fit(_lorentzian_asymmetric, x, y, p0=intital_guess)
        y_fit = _lorentzian_asymmetric(x, popt[0], popt[1], popt[2], popt[3])

    elif method == 'lorentzian':
        intital_guess = [amp, centre, fwhm]
        popt, pcov = curve_fit(_lorentzian_single, x, y, p0=intital_guess)
        y_fit = _lorentzian_single(x, popt[0], popt[1], popt[2])

    elif method == 'gaussian':
        intital_guess = [amp, centre, fwhm]
        popt, pcov = curve_fit(_gaussian, x, y, p0=intital_guess)
        y_fit = _gaussian(x, popt[0], popt[1], popt[2])

    integral_area = get_area(x, y_fit)
    output = [x, y, y_fit, popt, integral_area]

    if return_uncertainty:
        # Calculate residuals
        residuals = y - y_fit

        # Calculate R-squared
        ss_total = np.sum((y - np.mean(y))**2)
        ss_residual = np.sum(residuals**2)
        r_squared = 1 - (ss_residual / ss_total)

        # Calculate RMSE
        rmse = np.sqrt(np.mean(residuals**2))

        print(f"R-squared: {r_squared}")
        print(f"RMSE: {rmse}")
        output = output + [r_squared, rmse]

    if show_fit:
        fig, axs = plt.subplots(1,2, figsize = [10,5])
        if title is not None:
            axs[0].set_title(title)
        axs[0].plot(x, y, color = 'k', label = 'data')
        axs[0].plot(x, y_fit, color = 'red', label = 'fit: %s' % method)
        axs[0].invert_xaxis()
        axs[0].legend()
        lim1, lim2 = axs[0].get_ylim()
        axs[1].plot(x, y-y_fit, color = 'green', label = 'residual')
        axs[1].invert_xaxis()
        axs[1].legend()
        axs[1].set_ylim([lim1, lim2])
        if save:
            if fnm is None:
                fnm = '%s.png' % scan
            plt.savefig(fnm)
        else:
            plt.show()
        plt.close()

    return output

def get_SignalNoise(xdata, ydata, noise_bounds = [-2.0, -1.0], signal_bounds = [2.5, 2.7]):
    noise_mask = (xdata > noise_bounds[0]) * (xdata < noise_bounds[1])
    sig_mask = (xdata > signal_bounds[0]) * (xdata < signal_bounds[1])
    noise = get_area(xdata[noise_mask], ydata[noise_mask], avg = True)
    sig = get_area(xdata[sig_mask], ydata[sig_mask], avg = True)
    return sig, noise

def fit_TSP(xdata, ydata):
    out_tsp = fit_singlet(xdata, ydata, method = 'lorentzian_asymmetric', expected_centre = 0.0, show_fit = True)
    out_nse = fit_singlet(xdata, ydata, method = 'lorentzian_asymmetric', expected_centre = -2.0, show_fit = True)
    out_tsp_rs      = out_tsp[-2]                
    out_tsp_rmse    = out_tsp[-1]                    
    out_nse_rs      = out_nse[-2]
    out_nse_rmse    = out_nse[-1]                    
    return out_tsp_rs, out_tsp_rmse, out_nse_rs, out_nse_rmse

scans = [1, 2, 4, 8, 16, 32, 64, 128, 256]

centres, areas = [], []
SIG_RS = []
SIG_RMSE = []
NOISE_RS = []
NOISE_RMSE = []


for scan in scans:
    spectra_a = get_spectra(sample = 'D24', nscan = scan)
    '''
    print('Nscan: %s' % scan)
    a = fit_singlet(spectra_a[0], spectra_a[1], method = 'lorentzian_asymmetric', fnm = 'sig_%s.png' % scan, title = 'nscan = %s' % scan)
    n = fit_singlet(spectra_a[0], spectra_a[1], method = 'lorentzian_asymmetric', expected_centre = -1.0, fnm = 'noise_%s.png' % scan, title = 'nscan = %s' % scan)
    area = a[4]
    centre = a[3][1]
    SIG_RS.append(a[5])
    SIG_RMSE.append(a[6])
    NOISE_RS.append(n[5])
    NOISE_RMSE.append(n[6])    
    areas.append(area)
    centres.append(centre)
    '''
    out_tsp_rs, out_tsp_rmse, out_nse_rs, out_nse_rmse = fit_TSP(spectra_a[0], spectra_a[1])



'''
rs_lab = '$\\mathrm{R^{2}} = \\sum{(y_{i} - \hat{y}_{i}})^{2}/\\sum{(y_{i} - \overline{y})^2}$'

fig, axs = plt.subplots(1,2, figsize = [10, 5])
axs[0].set_title(rs_lab)
axs[0].plot(scans, SIG_RS, label = 'TSP fit')
axs[0].plot(scans, NOISE_RS, label = 'Noise fit')
axs[0].set_xscale('log')
axs[0].legend()

axs[1].set_title('RMSE of Residuals')
axs[1].plot(scans, SIG_RMSE)
axs[1].plot(scans, NOISE_RMSE)

plt.show()

shifts = centres[-1] - centres
norms = areas[-1]/areas

fig, axs = plt.subplots(1,2, figsize = [10, 5])

for i in range(0, len(scans)):
    scan = scans[i]
    spectra_a = get_spectra(sample = 'D24', nscan = scan)
    xdata = spectra_a[0] + shifts[i]
    ydata = spectra_a[1] * norms[i]
    axs[0].plot(spectra_a[0], spectra_a[1])
    axs[1].plot(xdata, ydata)
    
axs[0].set_title('raw data')
axs[1].set_title('shifted and normed data')
plt.show()
'''
