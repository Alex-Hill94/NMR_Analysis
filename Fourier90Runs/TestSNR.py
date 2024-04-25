from NMRClasses import *

def find_peaks(x, y, threshold=0.1, min_distance=10):
    """
    Finds the peak positions in the given signal.
    
    Parameters:
    x (numpy.ndarray): The x-axis data (e.g., ppm).
    y (numpy.ndarray): The y-axis data (e.g., intensity).
    threshold (float): The minimum relative height of a peak (default is 0.1).
    min_distance (int): The minimum distance between peaks (default is 10).
    
    Returns:
    numpy.ndarray: The indices of the peak positions.
    """
    # Find all local maxima
    peaks = np.zeros_like(y, dtype=bool)
    peaks[1:-1] = (y[1:-1] > y[:-2]) & (y[1:-1] > y[2:])
    
    # Filter peaks based on threshold and minimum distance
    peak_heights = y[peaks]
    peak_positions = np.where(peaks)[0]
    
    # Filter peaks based on threshold
    mask = peak_heights > threshold * np.max(peak_heights)
    peak_positions = peak_positions[mask]
    
    # Filter peaks based on minimum distance
    if min_distance > 0:
        keep = np.ones_like(peak_positions, dtype=bool)
        for i in range(len(peak_positions)):
            if i > 0 and peak_positions[i] - peak_positions[i-1] < min_distance:
                keep[i] = False
            elif i < len(peak_positions) - 1 and peak_positions[i+1] - peak_positions[i] < min_distance:
                keep[i] = False
        peak_positions = peak_positions[keep]
    
    return peak_positions

def peak_snr(x, y, peak_indices):
    """
    Calculates the signal-to-noise ratio for each peak in the spectrum.
    
    Parameters:
    x (numpy.ndarray): The x-axis data (e.g., ppm).
    y (numpy.ndarray): The y-axis data (e.g., intensity).
    peak_indices (list): The indices of the peak maxima.
    
    Returns:
    numpy.ndarray: The signal-to-noise ratio for each peak.
    """
    snr_values = []
    
    for peak_idx in peak_indices:
        # Get the peak height
        peak_height = y[peak_idx]
        
        # Calculate the noise standard deviation in a window around the peak
        noise_std = np.std(y[peak_idx-10:peak_idx+10])
        
        # Calculate the SNR for this peak
        snr = peak_height / noise_std
        snr_values.append(snr)
    
    return np.array(snr_values)

def baseline_snr(x, y):
    """
    Calculates the signal-to-noise ratio using the baseline noise.
    
    Parameters:
    x (numpy.ndarray): The x-axis data (e.g., ppm).
    y (numpy.ndarray): The y-axis data (e.g., intensity).
    
    Returns:
    float: The signal-to-noise ratio.
    """
    # Find the maximum peak height
    max_peak = np.max(y)
    
    # Calculate the standard deviation of the baseline noise
    noise_std = np.std(y[np.abs(y) < 0.1 * max_peak])
    
    # Calculate the SNR
    snr = max_peak / noise_std
    
    return snr

def rms_snr(x, y):
    """
    Calculates the RMS signal-to-noise ratio.
    
    Parameters:
    x (numpy.ndarray): The x-axis data (e.g., ppm).
    y (numpy.ndarray): The y-axis data (e.g., intensity).
    
    Returns:
    float: The RMS signal-to-noise ratio.
    """
    # Calculate the RMS of the signal
    signal_rms = np.sqrt(np.mean(y**2))
    
    # Calculate the RMS of the noise
    noise_rms = np.sqrt(np.mean(y[np.abs(y) < signal_rms]**2))
    
    # Calculate the RMS SNR
    snr = signal_rms / noise_rms
    
    return snr

def peak_to_peak_snr(x, y):
    """
    Calculates the peak-to-peak signal-to-noise ratio.
    
    Parameters:
    x (numpy.ndarray): The x-axis data (e.g., ppm).
    y (numpy.ndarray): The y-axis data (e.g., intensity).
    
    Returns:
    float: The peak-to-peak signal-to-noise ratio.
    """
    # Find the maximum peak height
    max_peak = np.max(y)
    
    # Calculate the standard deviation of the noise
    noise_std = np.std(y)
    
    # Calculate the peak-to-peak SNR
    snr = max_peak / noise_std
    
    return snr

def plot_spectrum_with_peaks(x, y, peak_indices):
    """
    Plots the NMR spectrum with detected peak locations and data point separation.
    
    Parameters:
    x (numpy.ndarray): The x-axis data (e.g., ppm).
    y (numpy.ndarray): The y-axis data (e.g., intensity).
    peak_indices (numpy.ndarray): The indices of the detected peak positions.
    """
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot the spectrum
    ax.plot(x, y, color='black')
    
    # Plot the peak locations
    ax.scatter(x[peak_indices], y[peak_indices], color='red', marker='x', label='Peak Locations')
    
    # Add vertical lines for data point separation
    for i in range(0, len(x), 10):
        ax.axvline(x[i], color='lightgray', linewidth=0.5, alpha=0.5)    
    # Set labels and title
    ax.set_xlabel('ppm')
    ax.set_ylabel('Intensity')
    ax.set_title('NMR Spectrum with Peaks')
    ax.legend()
    
    # Show the plot
    plt.show()

def snr_agilent(x, y, noise_bounds, signal_bounds):
    noise_mask = (x >= np.min(noise_bounds)) * (x <= np.max(noise_bounds))
    noise_x, noise_y = x[noise_mask], y[noise_mask]
    noise = (sum(noise_y**2)/len(noise_y))**(0.5)

    signal_mask = (x >= np.min(signal_bounds)) * (x <= np.max(signal_bounds))
    signal_x, signal_y = x[signal_mask], y[signal_mask]
    signal = np.max(signal_y)
    SNR = signal/noise
    ppm_max = signal_x[signal_y == signal]
    print("AGILENT", signal, noise, SNR)
    return SNR, signal, ppm_max

def snr_liverpool(x, y, noise_bounds, signal_bounds):
    noise_mask = (x >= np.min(noise_bounds)) * (x <= np.max(noise_bounds))
    signal_mask = (x >= np.min(signal_bounds)) * (x <= np.max(signal_bounds))
    noise = get_area(x[noise_mask], y[noise_mask], avg = True)
    sig = get_area(x[signal_mask], y[signal_mask], avg = True)
    snr = sig/noise 
    #plt.figure()
    #plt.plot(x, y)
    #plt.plot(x[signal_mask], y[signal_mask])
    #plt.plot(x[noise_mask], y[noise_mask])
    #plt.show()
    return snr, sig, noise

def snr_bruker(x, y, noise_bounds, signal_bounds):
    def _is_odd(num):
        return num % 2 != 0


    def _bruker_noise(n_values):
    
        isodd = _is_odd(len(n_values))
        if isodd:
            n_values = n_values[:-1]
        N = len(n_values)
        n = int((N-1)/2)
        n_list = np.arange(-n, n + 1, 1)
        n_list2 = np.arange(1, n+1, 1)

        part_one = 0
        for i in range(0, len(n_list)):
            update = n_values[i]**2
            part_one+=update
        
        part_two = 0
        for i in range(0, len(n_list)):
            update = n_values[i]
            part_two+=update
        part_two = part_two**2

        part_three = 0
        for i in n_list2:
            update = i * (n_values[i] -  n_values[-i])
            part_three += update
        part_three = 3 * (part_three**2) / (N**2 - 1)


        numerator = part_one - 1/N * (part_two + part_three)
        denominator = N - 1
        noise = np.sqrt(numerator/denominator)

        return noise

    noise_mask = (x >= np.min(noise_bounds)) * (x <= np.max(noise_bounds))
    signal_mask = (x >= np.min(signal_bounds)) * (x <= np.max(signal_bounds))

    signal_values = y[signal_mask]
    SIGNAL = np.nanmax(signal_values)
    signal_loc = x[signal_mask][signal_values == SIGNAL]

    noise_values = y[noise_mask]
    NOISE  = _bruker_noise(noise_values)
    SNR = SIGNAL / (2.* NOISE)
    print('NOISF1: %s NOISF2: %s' % (np.max(noise_bounds), np.min(noise_bounds)))
    print('SIG F1: %s SIG F2: %s' % (np.max(x), np.min(x)))
    print('Singal (%s ppm) / Noise' % signal_loc)
    print('%s/(%s*2) SINO: %s' % (SIGNAL, NOISE, SNR))
    return SNR, SIGNAL, NOISE




if __name__ == '__main__':

    scans = [1, 2, 4, 8, 16, 32, 64, 128, 256]

    lactic_acid_bounds = [1.2, 1.45]
    glucose_bounds     = [3.0, 4.2]
    citric_acid_bounds = [2.5, 2.7]
    noise_bounds       = [-2.0, -1.0]

    bounds = [lactic_acid_bounds, glucose_bounds, citric_acid_bounds]
    bound_lab = ['lac', 'glc', 'cit']
    for i in range(0, len(bounds)):
        SNRs = []
        bound = bounds[i]
        for scan in scans:

            S = LoadSpectra()
            S.ReadTextFile(nscan = scan, 
                        sample = 'D24',
                        pulse = 'zg30')
            #S.SubtractWater()
            x = S.initial_ppm
            y = S.initial_amplitude


            snr_agi, sig_agi, ppm_loc_agi = snr_agilent(x, y, noise_bounds, bound)
            snr_liv, sig_liv, noise_liv = snr_liverpool(x, y, noise_bounds, bound)
            snr_brk, sig_brk, noise_brk = snr_bruker(x, y, noise_bounds, bound)

            snr_array = [snr_agi, snr_liv, snr_brk]

            SNRs.append(snr_array)

        SNRs = np.array(SNRs)

        SNR_AGI = SNRs[:,0]
        SNR_LIV = SNRs[:,1]
        SNR_BRK = SNRs[:,2]

        fig, ax = plt.subplots(1,1, figsize = [10, 5])
        ax.set_title('%s'  % bound_lab[i])
        ax.plot(scans, np.ones(len(scans)) * 3., color = 'k', alpha = 0.5, ls = '--')
        ax.plot(scans, SNR_AGI, color = 'orange', label = 'SNR Agilent')
        ax.plot(scans, SNR_LIV, color = 'blue', label = 'SNR Liverpool')
        ax.plot(scans, SNR_BRK, color = 'green', label = 'SNR Bruker')
        ax.set_ylabel('SNR')
        ax.set_xlabel('N Scan')
        ax.legend()
        ax.set_yscale('log')
        ax.set_ylim([1, 205])
        plt.show()
        #plt.savefig('Paper_Figs/snr_measure_noesy_%s.png' % bound_lab[i])
        #plt.close()


    #lactate_colour = 'blue'
    #noise_colour   = 'black'
    #fill_alpha         = 0.1
#
    #fig, ax = plt.subplots(1,1, figsize = [10, 5])
    #ax.plot(x, y)
    #llim, ulim = ax.get_ylim()
    #ax.fill_betweenx([llim, ulim], lactic_acid_bounds[0], lactic_acid_bounds[1], color= lactate_colour, alpha=fill_alpha)
    #ax.fill_betweenx([llim, ulim], noise_bounds[0], noise_bounds[1], color= noise_colour, alpha=fill_alpha)
    #ax.scatter(ppm_loc_agi, sig_agi)
    #ax.plot(noise_bounds, [sig_agi/snr_agi, sig_agi/snr_agi])
    #ax.invert_xaxis()
    #ax.set_xlabel('ppm')
    #plt.show()

    #x = S.water_sub_ppm
    #y = S.water_sub_amplitude
    #peak_positions = find_peaks(x, y, threshold = 0.0001)
    #plot_spectrum_with_peaks(x, y, peak_positions)
    #plt.close()
    #plt.figure()
    #plt.plot(x, y)
    #plt.scatter(x[peak_positions], y[peak_positions], color = 'red')
    #plt.show()
