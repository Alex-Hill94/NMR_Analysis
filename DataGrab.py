import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pylab import rc, Line2D
from scipy.optimize import curve_fit



def grab_data(path = "/Users/alexhill/Desktop/Metabolomics/Data_Analysis", file_name = "60MHz_standards_230607.txt"):
    # Define the file path
    file_path = path+"/"+file_name

    # Extract the limits and size
    left_limit = None
    right_limit = None
    size = None

    with open(file_path, 'r') as file:
        for line in file:
            if 'LEFT' in line:
                left_index = line.find('LEFT')
                right_index = line.find('RIGHT')
                if left_index != -1:
                    left_start = line.find('=', left_index) + 1
                    left_end = line.find('ppm', left_start)
                    left_limit = float(line[left_start:left_end].strip())
                if right_index != -1:
                    right_start = line.find('=', right_index) + 1
                    right_end = line.find('ppm', right_start)
                    right_limit = float(line[right_start:right_end].strip())
            if 'SIZE' in line:
                size = int(line.split('=')[1].split('(')[0].strip())
                break  # Stop reading the file after finding thesize

    # Generate the x-axis values
    x_values = np.linspace(left_limit, right_limit, size)

    # Read the file and extract the number list
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                data.append(float(line))
    return x_values, np.array(data)

def plot_lines(ax, file_name = "60MHz_standards_230607.txt", col = 'k', lw = 1.5, ls = '-', alpha = 0.7, label = None, zoom = 1.):
    x_values, data = grab_data(file_name = file_name)
    ax.plot(x_values, np.array(data)*zoom, color = col, lw = lw, ls = ls, alpha = alpha, label = label)

def align(spec1, spec2):
    #spec1 = np.array([0, 0, 1, 1, 1, 1, 0, 0, 0])
    #spec2 = np.array([0, 1, 1, 1, 0, 0, 0, 0, 0])

    # Define a function to calculate the Chi-Squared value between the spectra
    def chi_squared(shift, spec1, spec2):
        shifted_spec2 = np.roll(spec2, int(shift))
        chi_sq = np.sum(((spec1 - shifted_spec2) ** 2) / (spec1 + shifted_spec2))
        return chi_sq

    # Initial guess for the shift parameter
    initial_shift = 0

    # Use scipy's minimize function to find the optimal shift that minimizes Chi-Squared
    result = minimize(chi_squared, initial_shift, args=(spec1, spec2), method='Nelder-Mead')

    # Extract the optimal shift from the result
    optimal_shift = result.x[0]

    # Shift spec2 using the optimal shift
    aligned_spec2 = np.roll(spec2, int(optimal_shift))

    # Print the optimal shift and aligned spectrum
    print("Optimal Shift:", optimal_shift)
    print("Aligned Spectrum 2:", aligned_spec2)
    return optimal_shift

def cut_idx(xdata, centre = 1.17, window = 0.3):
    assert np.all(xdata == np.flip(np.sort(xdata))), 'Error!'
    llim, ulim = centre - window, centre + window
    idxs = np.arange(0, len(xdata), 1)
    valid_idx = idxs[(xdata>=llim) * (xdata<=ulim)]
    return valid_idx

def get_centre(x, y, method = 'max', plot = True):
    centre_guess = x[y == np.max(y)][0]
    if method == 'gaussian':
        def gaussian(x, amplitude, mu, sigma):
            return amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

        # Initial guesses for the parameters

        initial_guess = [1e10, centre_guess, 1]

        # Perform curve fitting using curve_fit
        popt, pcov = curve_fit(gaussian, x, y, p0=initial_guess)

        # Extracting fitted parameters
        amplitude, centre, sigma = popt
        if plot:
            # Plotting the data and the fitted Gaussian
            plt.figure(figsize=(8, 6))
            plt.plot(x, y, 'b.', label='Data')
            plt.plot(x, gaussian(x, amplitude, centre, sigma), 'r-', label='Fitted Gaussian')
            plt.legend()
            plt.xlabel('X-axis')
            plt.ylabel('Intensity')
            plt.title('Fitting Gaussian to Data')
            plt.show()

    elif method == 'max':
        centre = centre_guess
    return centre

class AlignSpectra:
    def __init__(self):
        self.reference_x = []
        self.reference_y = []
        self.data_x_input = []
        self.data_y_input = []
        self.data_x_output = []
        self.data_y_output = []
        self.shift = []
        self.normalise = []

    def Centre( self, 
                reference_x,
                reference_y,
                data_x_input,
                data_y_input,
                method = 'gaussian', 
                plot = True):
        
        self.reference_x = reference_x
        self.reference_y = reference_y
        self.data_x_input = data_x_input
        self.data_y_input = data_y_input
                
        
        def get_centre(x, y, method = 'max', plot = True):
            centre_guess = x[y == np.max(y)][0]
            if method == 'gaussian':
                def gaussian(x, amplitude, mu, sigma):
                    return amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

                # Initial guesses for the parameters

                initial_guess = [np.max(y), centre_guess, 1]

                # Perform curve fitting using curve_fit
                popt, pcov = curve_fit(gaussian, x, y, p0=initial_guess)

                # Extracting fitted parameters
                amplitude, centre, sigma = popt
                if plot:
                    # Plotting the data and the fitted Gaussian
                    plt.figure(figsize=(8, 6))
                    plt.plot(x, y, 'b.', label='Data')
                    plt.plot(x, gaussian(x, amplitude, centre, sigma), 'r-', label='Fitted Gaussian')
                    plt.legend()
                    plt.xlabel('X-axis')
                    plt.ylabel('Intensity')
                    plt.title('Fitting Gaussian to Data')
                    plt.show()

            elif method == 'max':
                centre = centre_guess
            return centre

        ref_cent = get_centre(self.reference_x, self.reference_y, method = method, plot = plot)
        data_cent = get_centre(self.data_x_input, self.data_y_input, method = method, plot = plot)
        shift = data_cent - ref_cent
        
        self.data_x_output = self.data_x_input - shift
        self.data_y_output = self.data_y_input
        self.shift = shift

    def Normalise(self, plot = False):

        def get_dbl_gauss(x, y, plot = True):
            def dbl_gaussian(x, a1, mu1, sigma1, a2, sigma2, offset):
                g1 = a1 * np.exp(-0.5 * ((x - mu1) / sigma1) ** 2)
                g2 = a2 * np.exp(-0.5 * ((x - (mu1+offset)) / sigma2) ** 2)
                return g1 + g2

            #a1_guess = np.max(y)
            #mu1_guess = x[y == np.max(y)][0]
            #sigma1_guess = 0.2
            #a2_guess = np.max(y)
            #sigma2_guess = 0.2
            #off_guess = 0.1
            a1 = 1e8
            a2 = 1e8
            sigma1 = 0.02
            sigma2 = 0.02 
            mu1 = 1.13
            offset = 0.1

            initial_guess = [a1, mu1, sigma1, a2, sigma2, offset]
            # Perform curve fitting using curve_fit
            popt, pcov = curve_fit(dbl_gaussian, x, y, p0=initial_guess)
            a1, mu1, sigma1, a2, sigma2, offset = popt

            if plot:
                plt.figure()
                plt.plot(x, y)
                plt.plot(x, dbl_gaussian(x, a1, mu1, sigma1, a2, sigma2, offset), 'r-', label='Fitted Gaussian')
                plt.show()
            return a1, a2


        def cut_idx(xdata, centre = 1.17, window = 0.3):
            assert np.all(xdata == np.flip(np.sort(xdata))), 'Error!'
            llim, ulim = centre - window, centre + window
            idxs = np.arange(0, len(xdata), 1)
            valid_idx = idxs[(xdata>=llim) * (xdata<=ulim)]
            return valid_idx

        ref_cut = cut_idx(self.reference_x, centre = 1.25, window = 0.35)
        data_cut = cut_idx(self.data_x_output, centre = 1.25, window = 0.35)

        ref_a1, ref_a2 = get_dbl_gauss(self.reference_x[ref_cut], self.reference_y[ref_cut], plot = plot)
        data_a1, data_a2 = get_dbl_gauss(self.data_x_output[data_cut], self.data_y_output[data_cut], plot = plot)

        norm = ((ref_a1/data_a1) + (ref_a2/data_a2))/2.
        self.data_y_output = self.data_y_output * norm
        self.normalise = norm

    def CompareSpectra(self):
        fig, axs = plt.subplots(1,1, figsize = [5, 5])
        axs.plot(self.reference_x, self.reference_y, color = 'k', label = 'Reference')
        axs.plot(self.data_x_input, self.data_y_input, color = 'blue', ls = '-', label = 'Data Orig.')
        axs.plot(self.data_x_output, self.data_y_output, color = 'blue', ls = '--', label = 'Data Shift and Norm.')
        axs.invert_xaxis()
        axs.set_xlabel('ppm')
        axs.set_ylabel('amplitude')
        axs.legend()
        plt.show()



def align(ref_x, ref_y, dat_x, dat_y):
    D = AlignSpectra()
    D.Centre(reference_x = ref_x,
                        reference_y = ref_y,
                        data_x_input = dat_x,
                        data_y_input = dat_y, 
                        plot = True)
    D.Normalise(plot = True)
    D.CompareSpectra()
    return D.data_x_output, D.data_y_output

if __name__ == "__main__":

    x_values_1, data_1 = grab_data(file_name = "Nov29-2023_D22_LC+NG.txt")
    x_values_2, data_2 = grab_data(file_name = "Nov29-2023_D23_LC+BG.txt")
    x_values_3, data_3 = grab_data(file_name = "Nov29-2023_D24_LC+AG.txt")

    x1, y1 = x_values_1, data_1
    x2, y2 = align(x1, y1, x_values_2, data_2)
    x3, y3 = align(x1, y1, x_values_3, data_3)
    '''
    fig, axs = plt.subplots(1,1, figsize = [5, 5])
    axs.plot(x1, y1, label = 'LC+NG')
    axs.plot(x2, y2, label = 'LC+BG')
    axs.plot(x3, y3, label = 'LC+AG')
    axs.invert_xaxis()
    axs.set_xlabel('ppm')
    axs.set_ylabel('amplitude')
    axs.legend()
    plt.show()

    fig, axs = plt.subplots(2, 2, figsize = [8, 8])

    my_axs = [axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]]
    axs[0,0].set_xlim([-2, 6])

    axs[0,1].set_xlim([1.06, 1.29])
    axs[0,1].set_ylim([-0.25*1e8, 1.5*1e8])

    axs[1,0].set_xlim([-0.4, 0.12])
    axs[1,0].set_ylim([-1*1e7, 4.6*1e7])

    axs[1,1].set_xlim([2.7, 5.8])
    axs[1,1].set_ylim([-0.26, 3.51*1e8])

    for ax in my_axs:
        ax.plot(x1, y1, label = 'LC+NG')
        ax.plot(x2, y2, label = 'LC+BG')
        ax.plot(x3, y3, label = 'LC+AG')
        ax.invert_xaxis()
        ax.set_xlabel('ppm')
        ax.set_ylabel('amplitude')

    axs[0,0].legend()


    plt.show()
    '''
    matplotlib.rcParams.update({'font.size': 15})
    fig, axs = plt.subplots(1, 1, figsize = [6, 6])
    axs.plot(x1, y1, label = '5 mM')
    axs.plot(x2, y2, label = '7 mM')
    axs.plot(x3, y3, label = '10 mM')
    axs.set_ylim([-0.26, 3.51*1e8])
    axs.invert_xaxis()
    axs.set_xlabel('ppm')
    axs.set_ylabel('amplitude')
    axs.set_title('Glucose')
    axs.set_xlim([3., 3.9])
    axs.legend()
    plt.savefig('glucose_exp_comp.png')
    plt.close()