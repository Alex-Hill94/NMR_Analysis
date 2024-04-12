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

if __name__ == '__main__':
    S = LoadSpectra()
    S.ReadTextFile(nscan = 256, 
                sample = 'D24',
                pulse = 'zg30')
    S.SubtractWater()
    x = S.initial_ppm
    y = S.initial_amplitude

    #x = S.water_sub_ppm
    #y = S.water_sub_amplitude
    peak_positions = find_peaks(x, y, threshold = 0.0001)
    plot_spectrum_with_peaks(x, y, peak_positions)
    #plt.close()
    #plt.figure()
    #plt.plot(x, y)
    #plt.scatter(x[peak_positions], y[peak_positions], color = 'red')
    #plt.show()
