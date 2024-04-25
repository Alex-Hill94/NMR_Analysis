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


def snr_agilent(x, y, noise_bounds, signal_bounds, verbose = False):
    noise_mask = (x >= np.min(noise_bounds)) * (x <= np.max(noise_bounds))
    noise_x, noise_y = x[noise_mask], y[noise_mask]
    noise = (sum(noise_y**2)/len(noise_y))**(0.5)

    signal_mask = (x >= np.min(signal_bounds)) * (x <= np.max(signal_bounds))
    signal_x, signal_y = x[signal_mask], y[signal_mask]
    signal = np.max(signal_y)
    SNR = signal/noise
    ppm_max = signal_x[signal_y == signal]
    if verbose:
        print("AGILENT", signal, noise, SNR)
    return SNR, signal, noise, ppm_max

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

def snr_bruker(x, y, noise_bounds, signal_bounds, verbose = False):
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
    if verbose:
        print('NOISF1: %s NOISF2: %s' % (np.max(noise_bounds), np.min(noise_bounds)))
        print('SIG F1: %s SIG F2: %s' % (np.max(x), np.min(x)))
        print('Singal (%s ppm) / Noise' % signal_loc)
        print('%s/(%s*2) SINO: %s' % (SIGNAL, NOISE, SNR))
    return SNR, SIGNAL, NOISE



def comp_time_snr(sample = 'D24', pulse = 'zg30', scan = 256, sig_bounds = [3.0, 4.2], noise_bounds = [-2.0, -1.0], snr_choice = 'brk'):
    if pulse == 'zg30':
        path = "/Users/alexhill/Desktop/Metabolomics/Rerun_Data/20240119_RERUN_3"
    else:
        path = "/Users/alexhill/Desktop/Metabolomics/Rerun_Data/20240119_RERUN_3"

    S = LoadSpectra()
    S.ReadTextFile(nscan = scan, 
                sample = 'D24',
                pulse = pulse)
    
    pulse1 = pulse
    if pulse == 'noesypr1dwv2':
        pulse1 = 'noesypr1dwv2.mjh'
    S.ReadRawData(path =  path,
        sample = 'D24',
        pulse = pulse1,
        nscan = scan)

    time_taken = S.time_taken
    x = S.initial_ppm
    y = S.initial_amplitude

    if snr_choice == 'brk':
        snr, sig, noise = snr_bruker(x, y, noise_bounds, sig_bounds)
    elif snr_choice == 'agi':
        snr, sig, noise, ppm_loc = snr_agilent(x, y, noise_bounds, sig_bounds)
    elif snr_choice == 'liv':
        snr, sig, noise = snr_liverpool(x, y, noise_bounds, sig_bounds)

    return snr, sig, noise, time_taken



def comp_snr():

    scans = [8, 16, 32, 64, 128, 256]

    lactic_acid_bounds = [1.2, 1.45]
    glucose_bounds     = [3.0, 4.2]
    citric_acid_bounds = [2.5, 2.7]
    noise_bounds       = [-2.0, -1.0]
    pulses = ['zg30', 'noesypr1dwv2']
    bounds = [lactic_acid_bounds, glucose_bounds, citric_acid_bounds]
    bound_lab = ['lac', 'glc', 'cit']
    snr_choices = ['brk', 'liv', 'agi']

    DATA = []

    for pulse in pulses:
        for i in range(0, len(bounds)):
            bound = bounds[i]
            blab = bound_lab[i]
            for scan in scans:
                for snr_choice in snr_choices:
                    SNR, SIG, NOISE, TIME = comp_time_snr(sample = 'D24', pulse = pulse, scan = scan, sig_bounds = bound, noise_bounds = [-2.0, -1.0], snr_choice = snr_choice)
                    data = [pulse, blab, scan, snr_choice, SNR, TIME]
                    DATA.append(data)

    glucose_colour = 'orange'
    citrate_colour = 'green'
    lactate_colour = 'blue'
    zg30_ls = '-'
    noes_ls = '--'
    zg30_mrk = 's'
    noes_mrk = 'o'
    
    DATA = np.array(DATA)
    pulse_array = DATA[:,0].astype('str')
    bound_array = DATA[:,1].astype('str')
    scan_array  = DATA[:,2].astype('int')
    snrchoice_array   = DATA[:,3].astype('str')
    snr_array = DATA[:,4].astype('float')
    time_array = DATA[:,5].astype('int')

    glc_mask = (bound_array == 'glc')
    lac_mask = (bound_array == 'lac')
    cit_mask = (bound_array == 'cit')

    agi_mask = (snrchoice_array == 'agi')
    liv_mask = (snrchoice_array == 'liv')
    brk_mask = (snrchoice_array == 'brk')

    zg3_mask = (pulse_array == 'zg30')
    noe_mask = (pulse_array == 'noesypr1dwv2')
    
    snrs = snr_array[glc_mask * brk_mask * zg3_mask]
    scans = scan_array[glc_mask * brk_mask * zg3_mask]

    plt.close()

    lw = 0.5

    matplotlib.rc('text', usetex=True)

    # Set the default font to Times New Roman
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.size'] = 12



    fig, axs = plt.subplots(1,2, figsize = [12, 5])

    axs[0].plot(scan_array[glc_mask * brk_mask * zg3_mask], snr_array[glc_mask * brk_mask * zg3_mask], color = glucose_colour, ls = zg30_ls, lw = lw)
    axs[0].plot(scan_array[glc_mask * brk_mask * noe_mask], snr_array[glc_mask * brk_mask * noe_mask], color = glucose_colour, ls = noes_ls, lw = lw)
    axs[0].scatter(scan_array[glc_mask * brk_mask * zg3_mask], snr_array[glc_mask * brk_mask * zg3_mask], edgecolor = 'k', color = glucose_colour, marker = zg30_mrk)
    axs[0].scatter(scan_array[glc_mask * brk_mask * noe_mask], snr_array[glc_mask * brk_mask * noe_mask], edgecolor = 'k', color = glucose_colour, marker = noes_mrk)

    axs[0].plot(scan_array[lac_mask * brk_mask * zg3_mask], snr_array[lac_mask * brk_mask * zg3_mask], color = lactate_colour, ls = zg30_ls, lw = lw)
    axs[0].plot(scan_array[lac_mask * brk_mask * noe_mask], snr_array[lac_mask * brk_mask * noe_mask], color = lactate_colour, ls = noes_ls, lw = lw)
    axs[0].scatter(scan_array[lac_mask * brk_mask * zg3_mask], snr_array[lac_mask * brk_mask * zg3_mask], edgecolor = 'k', color = lactate_colour, marker = zg30_mrk)
    axs[0].scatter(scan_array[lac_mask * brk_mask * noe_mask], snr_array[lac_mask * brk_mask * noe_mask], edgecolor = 'k', color = lactate_colour, marker = noes_mrk)

    axs[0].plot(scan_array[cit_mask * brk_mask * zg3_mask], snr_array[cit_mask * brk_mask * zg3_mask], color = citrate_colour, ls = zg30_ls, lw = lw)
    axs[0].plot(scan_array[cit_mask * brk_mask * noe_mask], snr_array[cit_mask * brk_mask * noe_mask], color = citrate_colour, ls = noes_ls, lw = lw)
    axs[0].scatter(scan_array[cit_mask * brk_mask * zg3_mask], snr_array[cit_mask * brk_mask * zg3_mask], edgecolor = 'k', color = citrate_colour, marker = zg30_mrk)
    axs[0].scatter(scan_array[cit_mask * brk_mask * noe_mask], snr_array[cit_mask * brk_mask * noe_mask], edgecolor = 'k', color = citrate_colour, marker = noes_mrk)

    axs[0].set_ylabel('$\mathrm{SNR}_{\mathrm{Brk.}}$')
    axs[0].set_xlabel('Number of scans ($n_{\mathrm{s}}$)')
    #axs[0].set_yscale('log')
    #axs[1].set_yscale('log')

    custom_lines = [Line2D([0], [0], color=glucose_colour, linestyle='-', lw = 2),
                    Line2D([0], [0], color=lactate_colour, linestyle='-', lw = 2),
                    Line2D([0], [0], color=citrate_colour, linestyle='-', lw = 2),
                    Line2D([0], [0], color='k', marker=zg30_mrk),
                    Line2D([0], [0], color='k', marker=noes_mrk, linestyle = '--')]
    
    axs[0].legend(custom_lines, ['Glucose [10 mM]', 'Lactate [2 mM]', 'Citrate [0.2 mM]', 'zg30', 'noesypr1dwv2'])


    axs[1].plot(time_array[glc_mask * brk_mask * zg3_mask], snr_array[glc_mask * brk_mask * zg3_mask], color = glucose_colour, ls = zg30_ls, lw = lw)
    axs[1].plot(time_array[glc_mask * brk_mask * noe_mask], snr_array[glc_mask * brk_mask * noe_mask], color = glucose_colour, ls = noes_ls, lw = lw)
    axs[1].scatter(time_array[glc_mask * brk_mask * zg3_mask], snr_array[glc_mask * brk_mask * zg3_mask], edgecolor = 'k', color = glucose_colour, marker = zg30_mrk)
    axs[1].scatter(time_array[glc_mask * brk_mask * noe_mask], snr_array[glc_mask * brk_mask * noe_mask], edgecolor = 'k', color = glucose_colour, marker = noes_mrk)

    axs[1].plot(time_array[lac_mask * brk_mask * zg3_mask], snr_array[lac_mask * brk_mask * zg3_mask], color = lactate_colour, ls = zg30_ls, lw = lw)
    axs[1].plot(time_array[lac_mask * brk_mask * noe_mask], snr_array[lac_mask * brk_mask * noe_mask], color = lactate_colour, ls = noes_ls, lw = lw)
    axs[1].scatter(time_array[lac_mask * brk_mask * zg3_mask], snr_array[lac_mask * brk_mask * zg3_mask], edgecolor = 'k', color = lactate_colour, marker = zg30_mrk)
    axs[1].scatter(time_array[lac_mask * brk_mask * noe_mask], snr_array[lac_mask * brk_mask * noe_mask], edgecolor = 'k', color = lactate_colour, marker = noes_mrk)

    axs[1].plot(time_array[cit_mask * brk_mask * zg3_mask], snr_array[cit_mask * brk_mask * zg3_mask], color = citrate_colour, ls = zg30_ls, lw = lw)
    axs[1].plot(time_array[cit_mask * brk_mask * noe_mask], snr_array[cit_mask * brk_mask * noe_mask], color = citrate_colour, ls = noes_ls, lw = lw)
    axs[1].scatter(time_array[cit_mask * brk_mask * zg3_mask], snr_array[cit_mask * brk_mask * zg3_mask], edgecolor = 'k', color = citrate_colour, marker = zg30_mrk)
    axs[1].scatter(time_array[cit_mask * brk_mask * noe_mask], snr_array[cit_mask * brk_mask * noe_mask], edgecolor = 'k', color = citrate_colour, marker = noes_mrk)

    axs[1].set_ylabel('$\mathrm{SNR}_{\mathrm{Brk.}}$')
    axs[1].set_xlabel('Experiment time [secs]')

    my_axs = [axs[0], axs[1]]

    for ax in my_axs:
        ax.tick_params(axis='x', direction='in', which='both')
        ax.tick_params(axis='y', direction='in', which='both')


    plt.show()
    #plt.savefig('Paper_Figs/pulses_snr.pdf')
    #plt.close()

   #plt.figure()
   #plt.plot(scan_array[cit_mask * brk_mask * zg3_mask], time_array[cit_mask * brk_mask * zg3_mask], label = 'zg30', ls = zg30_ls, color = 'k')
   #plt.plot(scan_array[cit_mask * brk_mask * noe_mask], time_array[cit_mask * brk_mask * noe_mask], label = 'zg30', ls = noes_ls, color = 'k')
   #plt.xlabel('Scans')
   #plt.ylabel('Time')
   #plt.savefig('Paper_Figs/scan_time_comp.png')



def plot_pulses():

    scans = [1, 2, 4, 8, 16, 32, 64, 128, 256]

    lactic_acid_bounds = [1.2, 1.45]
    glucose_bounds     = [3.0, 4.2]
    citric_acid_bounds = [2.5, 2.7]
    noise_bounds       = [-2.0, -1.0]
    pulses = ['zg30', 'noesypr1dwv2']
    bounds = [lactic_acid_bounds, glucose_bounds, citric_acid_bounds]
    bound_lab = ['lac', 'glc', 'cit']
    snr_choices = ['brk', 'liv', 'agi']
    glucose_colour = 'orange'
    citrate_colour = 'green'
    lactate_colour = 'blue'
    noise_colour   = 'black'
    fill_alpha         = 0.05

    import matplotlib as mpl
    cmap = mpl.colormaps['BrBG']
    colours = cmap(np.linspace(0, 1, len(scans) + 2))


    X_ZG = []
    Y_ZG = []
    X_NO = []
    Y_NO = []

    plt.close()

    # Enable LaTeX rendering
    matplotlib.rc('text', usetex=True)

    # Set the default font to Times New Roman
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.size'] = 12


    fig, axs = plt.subplots(2, 1, figsize = [12, 7])
    for i in range(3, len(scans)):
        scan = scans[i]
        S = LoadSpectra()
        S.ReadTextFile(nscan = scan, sample = 'D24')
        X_ZG.append(S.initial_ppm)
        Y_ZG.append(S.initial_amplitude)
        axs[0].plot(S.initial_ppm, S.initial_amplitude, color = colours[i+1], lw = 1)
        
        S = LoadSpectra()
        S.ReadTextFile(nscan = scan, sample = 'D24', pulse = 'noesypr1dwv2')
        X_NO.append(S.initial_ppm)
        Y_NO.append(S.initial_amplitude)
        axs[1].plot(S.initial_ppm, S.initial_amplitude, color = colours[i+1], lw = 1)


    custom_line = [Line2D([0], [0], color='w', linestyle='-')]
    
    label = 'Solvent peak\n $\\times$ 180'

    x_pos_glc = 0.7 + sum(glucose_bounds)/2. 
    y_pos_main = 1.2*1e9
    axs[0].annotate(label, xy=(x_pos_glc - 0.12, y_pos_main), xytext=(x_pos_glc - 0.12, y_pos_main), horizontalalignment='center', color = 'k')



    glc_label = '\\textbf{Glucose}\n \\textbf{10mM}'
    cit_label = '\\textbf{Citrate}\n \\textbf{0.2mM}'
    lac_label = '\\textbf{Lactate}\n \\textbf{2mM}'
    noise_label = '\\textbf{Noise}'

    #axs[0].annotate('   ', xy=(x_pos_glc + 0.35, y_pos_main*1.12), xytext=(x_pos_glc + 0.35, y_pos_main*1.12 - (0.1*1e9)), horizontalalignment='center', color = 'k',
     #               arrowprops=dict(facecolor='black', shrink=0.05))
    axs[0].arrow(x = x_pos_glc + 0.3, y = y_pos_main*1.12 - (0.1*1e9), dx = 0, dy = 0.1*1e9, head_length = 0.03*1e9, width = 0.025, facecolor= 'k')
    x_pos = 0.03+sum(glucose_bounds)/2.
    y_pos = 1.2*1e9
    axs[0].annotate(glc_label, xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = glucose_colour, fontsize = 10)

    x_pos = sum(citric_acid_bounds)/2.
    y_pos = 1.2*1e9
    axs[0].annotate(cit_label, xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = citrate_colour, fontsize = 10)


    x_pos = sum(lactic_acid_bounds)/2.
    y_pos = 1.2*1e9
    axs[0].annotate(lac_label, xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = lactate_colour, fontsize = 10)


    my_axs = [axs[0], axs[1]]

    axs[0].set_ylim([-2. * 1e7, 1.41 *1e9])
    axs[1].set_ylim([-0.75 * 1e9, 2.304 *1e9])
    for ax in my_axs:
        llim, ulim = ax.get_ylim()
        ax.fill_betweenx([llim, ulim], glucose_bounds[0], glucose_bounds[1], color= glucose_colour, alpha=fill_alpha)
        ax.fill_betweenx([llim, ulim], lactic_acid_bounds[0], lactic_acid_bounds[1], color= lactate_colour, alpha=fill_alpha)
        ax.fill_betweenx([llim, ulim], citric_acid_bounds[0], citric_acid_bounds[1], color= citrate_colour, alpha=fill_alpha)
        ax.tick_params(axis='x', direction='in', which='both')
        ax.tick_params(axis='y', direction='in', which='both')

        ax.set_xlim([-0.5, 5.5])
        ax.invert_xaxis()

    axs[0].set_xticklabels([])
    axs[1].set_xlabel('ppm')
    zg_lab = '\\textbf{zg30}'
    noes_lab = '\\textbf{noesypr1dwv2}'
    
    axs[0].legend(custom_line, [zg_lab], frameon = False, loc = 'upper right')
    axs[1].legend(custom_line, [noes_lab], frameon = False, loc = 'upper right')

    import matplotlib.colorbar as colorbar
    import matplotlib.colors as clr

    cb_colors = []
    xtick_labs = []

    for i in range(0, len(scans)):
        cb_colors.append(colours[i + 1])
        xtick_labs.append("%s" % scans[i])

    num_colors = len(cb_colors)
    cmap_ = clr.ListedColormap(cb_colors)

    top_pos = axs[0].get_position()
    bottom_pos = axs[1].get_position()

    # Calculate the position for the colorbar
    colorbar_left = 0.92
    colorbar_bottom = bottom_pos.y0
    colorbar_width = 0.02
    colorbar_height = top_pos.y1 - bottom_pos.y0
    colorbar_rect = [colorbar_left, colorbar_bottom, colorbar_width, colorbar_height]


    ax = fig.add_axes(colorbar_rect)

    cb = colorbar.ColorbarBase(ax, orientation='vertical',
                            cmap=cmap_, norm=plt.Normalize(-0.5, num_colors - 0.5))

    #print(range(num_colors))
    #print(xtick_labs)
#
    cb.set_ticks(range(num_colors))
    cb.ax.set_yticklabels(xtick_labs)
    cb.set_label('Number of scans ($n_{\mathrm{s}}$)', fontsize = 13)
    


    plt.savefig('Paper_Figs/pulse_comparison.pdf')
    #plt.show()

    plt.close()

def get_snr(sample = 'D24', pulse = 'noesypr1dwv2', scan = 256, sig_bounds = [3.0, 4.2], noise_bounds = [-2.0, -1.0], snr_choice = 'brk'): 

    SNR, SIG, NOISE, TIME = comp_time_snr(sample = sample, pulse = pulse, scan = scan, sig_bounds = sig_bounds, noise_bounds =noise_bounds, snr_choice = snr_choice)
    data = [pulse, scan, snr_choice, SNR, SIG, NOISE, TIME]
    return data

#plot_pulses()
comp_snr()
#data = get_snr()

#if __name__ == '__main__':
#    plot_pulses()