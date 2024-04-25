from NMRClasses import *
import matplotlib

def big_plot():
    scans = [1, 2, 4, 8, 16, 32, 64, 128, 256]
    sample = 'D24'
    times = []

    for scan in scans:
        S = LoadSpectra()
        S.ReadRawData(nscan = scan, 
                    sample = sample,
                    pulse = 'zg30')
        tt = S.time_taken
        times.append(tt)
    print (scans)
    print(times)

    S = LoadSpectra()
    S.ReadTextFile(sample = sample, nscan = 256)
    S.ReadRawData(sample=sample)
    S.SubtractWater()

    A = AnalyseSpectra()
    A.InputData(x = S.water_sub_ppm, y = S.water_sub_amplitude)
    A.FitTSP()

    ref_x, ref_y = A.tsp_centre, A.tsp_integral

    citric_acid_bounds = [2.37, 2.82]
    lactic_acid_bounds = [1.17, 1.50]
    glucose_bounds     = [3.19, 3.98]
    noise_bounds       = [-2.0, -1.0]
    import matplotlib as mpl

    cmap = mpl.colormaps['BrBG']

    colours = cmap(np.linspace(0, 1, len(scans) + 2))

    X_ORIG = []
    Y_ORIG = []

    X_SCALED = []
    Y_SCALED = []

    SNR_ORIG = []
    SNR_SCALED = []


    for scan in scans:
        print(scan)
        S = LoadSpectra()
        S.ReadTextFile(nscan = scan, sample =sample)
        S.SubtractWater()
        X_ORIG.append(S.water_sub_ppm)
        Y_ORIG.append(S.water_sub_amplitude)

        A = AnalyseSpectra()
        A.InputData(x = S.water_sub_ppm, y = S.water_sub_amplitude)
        
        A.SignalToNoise(signal_bounds = citric_acid_bounds, snr_choice = 'liv')
        citric_acid_snr = A.snr

        A.SignalToNoise(signal_bounds = lactic_acid_bounds, snr_choice = 'liv')
        lactic_acid_snr = A.snr

        A.SignalToNoise(signal_bounds = glucose_bounds, snr_choice = 'liv')
        glucose_snr = A.snr

        snrs = [citric_acid_snr, lactic_acid_snr, glucose_snr]

        SNR_ORIG.append(snrs)


        A.FitTSP()
        tsp_x, tsp_y = A.tsp_centre, A.tsp_integral
        y_scale = ref_y/tsp_y
        x_shift = ref_x - tsp_x
        A.ScaleSpectra(y_scaling = y_scale, x_shift = x_shift)

        A.SignalToNoise(signal_bounds = citric_acid_bounds, snr_choice = 'liv')
        citric_acid_snr = A.snr

        A.SignalToNoise(signal_bounds = lactic_acid_bounds, snr_choice = 'liv')
        lactic_acid_snr = A.snr

        A.SignalToNoise(signal_bounds = glucose_bounds, snr_choice = 'liv')
        glucose_snr = A.snr

        snrs = [citric_acid_snr, lactic_acid_snr, glucose_snr]


        X_SCALED.append(A.processed_ppm)
        Y_SCALED.append(A.processed_amplitude)
        SNR_SCALED.append(snrs)

    snrs_orig = np.array(SNR_ORIG)
    snrs_scale = np.array(SNR_SCALED)


    glucose_colour = 'orange'
    citrate_colour = 'green'
    lactate_colour = 'blue'
    noise_colour   = 'black'
    fill_alpha         = 0.1


    custom_lines = []
    legend_labels = []

    # Create figure and gridspec
    fig = plt.figure(figsize=(10, 6))
    gs = fig.add_gridspec(2, 2)

    # Top panel spanning two columns
    ax_top = fig.add_subplot(gs[0, :])
    for i in range(0, len(scans)):
        ax_top.plot(X_ORIG[i], Y_ORIG[i], label = 'nscan = %s' % scans[i], lw = 0.5, color = colours[i + 1])
        custom_lines.append(Line2D([0], [0], color=colours[i+1], linestyle='-'))
        legend_labels.append('ns = %s' % scans[i])

    llim, ulim = ax_top.get_ylim()
    ax_top.fill_betweenx([llim, ulim], citric_acid_bounds[0], citric_acid_bounds[1], color= citrate_colour, alpha=fill_alpha)
    ax_top.fill_betweenx([llim, ulim], lactic_acid_bounds[0], lactic_acid_bounds[1], color= lactate_colour, alpha=fill_alpha)
    ax_top.fill_betweenx([llim, ulim], glucose_bounds[0], glucose_bounds[1], color=glucose_colour, alpha=fill_alpha)
    ax_top.fill_betweenx([llim, ulim], noise_bounds[0], noise_bounds[1], color= noise_colour, alpha=fill_alpha)

    ax_top.set_xlim([-2.2, 5.5])
    ax_top.set_ylim([-0.5 * 1e8, 1.6 *1e9])
    ax_top.set_xlabel('ppm')
    ax_top.invert_xaxis()
    ax_top.legend(custom_lines, legend_labels, frameon = False)

    # Bottom left panel
    ax_bottom_left = fig.add_subplot(gs[1, 0])

    ax_bottom_left.plot(scans, np.ones(len(scans))*3., color = 'k', lw = 0.7, alpha = 0.5, ls = '--')
    ax_bottom_left.plot(scans, snrs_orig[:,0], color = citrate_colour, alpha = fill_alpha*8)
    ax_bottom_left.plot(scans, snrs_orig[:,1], color = lactate_colour, alpha = fill_alpha*8)
    ax_bottom_left.plot(scans, snrs_orig[:,2], color = glucose_colour, alpha = fill_alpha*8)
    ax_bottom_left.set_xlabel('Scans')
    ax_bottom_left.set_ylabel('SNR')




    # Bottom right panel
    ax_bottom_right = fig.add_subplot(gs[1, 1])

    ax_bottom_right.plot(times, np.ones(len(times))*3., color = 'k', lw = 0.7, alpha = 0.5, ls = '--')
    ax_bottom_right.plot(times, snrs_orig[:,0], color = citrate_colour, alpha = fill_alpha*8)
    ax_bottom_right.plot(times, snrs_orig[:,1], color = lactate_colour, alpha = fill_alpha*8)
    ax_bottom_right.plot(times, snrs_orig[:,2], color = glucose_colour, alpha = fill_alpha*8)
    ax_bottom_right.set_xlabel('Time (secs)')
    ax_bottom_right.set_ylabel('SNR')
    ax_bottom_right.set_yscale('log')
    ax_bottom_right.set_xscale('log')

    x_pos = sum(glucose_bounds)/2.
    y_pos = 1.3*1e9
    ax_top.annotate('Glucose\n10mM', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = glucose_colour)

    x_pos = sum(citric_acid_bounds)/2.
    y_pos = 1.3*1e9
    ax_top.annotate('Citrate\n0.2mM', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = citrate_colour)


    x_pos = sum(lactic_acid_bounds)/2.
    y_pos = 1.3*1e9
    ax_top.annotate('Lactate\n2mM', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = lactate_colour)

    x_pos = sum(noise_bounds)/2.
    y_pos = 0.075*1e9
    ax_top.annotate('Noise', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = noise_colour)


    my_axs = [ax_top, ax_bottom_left, ax_bottom_right]

    for ax in my_axs:
        ax.tick_params(axis='x', direction='in', which='both')
        ax.tick_params(axis='y', direction='in', which='both')

    # Adjust layout
    fig.tight_layout()

    # Show plot

    plt.savefig('classes_d24.pdf')
    #plt.show()
    plt.close()

    #'matplotlib-label-line'

def just_spectra():


    fig, axs = plt.subplots(1,1, figsize = [12, 4])
    for i in range(0, len(scans)):
        axs.plot(X_ORIG[i], Y_ORIG[i], label = 'nscan = %s' % scans[i], lw = 0.5, color = colours[i + 1])

    axs.set_xlim([-2.5, 5.5])
    axs.set_ylim([-2. * 1e8, 0.7 *1e9])
    axs.invert_xaxis()
    axs.legend(ncol = 2)
    axs.set_title('Water subbed spectra')
    axs.set_xlabel('ppm')
    llim00, ulim00 = axs.get_ylim()
    axs.fill_betweenx([llim00, ulim00], citric_acid_bounds[0], citric_acid_bounds[1], color= citrate_colour, alpha=fill_alpha)
    axs.fill_betweenx([llim00, ulim00], lactic_acid_bounds[0], lactic_acid_bounds[1], color= lactate_colour, alpha=fill_alpha)
    axs.fill_betweenx([llim00, ulim00], glucose_bounds[0], glucose_bounds[1], color=glucose_colour, alpha=fill_alpha)
    axs.fill_betweenx([llim00, ulim00], noise_bounds[0], noise_bounds[1], color= noise_colour, alpha=fill_alpha)

    plt.show()

def scaled_vs_not_scaled():


    fig, axs = plt.subplots(2,2, figsize = [8, 8])

    for i in range(0, len(scans)):
        axs[0,0].plot(X_ORIG[i], Y_ORIG[i], label = 'nscan = %s' % scans[i], lw = 0.5, color = colours[i + 1])
        axs[0,1].plot(X_SCALED[i], Y_SCALED[i], lw = 0.5, color = colours[i + 1])

    axs[0,0].set_xlim([-2.5, 5.5])
    axs[0,1].set_xlim([-2.5, 5.5])
    axs[0,0].set_ylim([-2. * 1e8, 0.7 *1e9])
    #axs[1].set_ylim([-2. * 1e8, 0.7 *1e9])



    axs[0, 0].invert_xaxis()
    axs[0, 1].invert_xaxis()
    axs[0, 0].legend(ncol = 2)

    axs[1,0].plot(scans, snrs_orig[:,0], color = citrate_colour, alpha = fill_alpha*3)
    axs[1,0].plot(scans, snrs_orig[:,1], color = lactate_colour, alpha = fill_alpha*3)
    axs[1,0].plot(scans, snrs_orig[:,2], color = glucose_colour, alpha = fill_alpha*3)

    axs[1,1].plot(scans, snrs_scale[:,0], label = 'Cit', color = citrate_colour, alpha = fill_alpha*3)
    axs[1,1].plot(scans, snrs_scale[:,1], label = 'Lac', color = lactate_colour, alpha = fill_alpha*3)
    axs[1,1].plot(scans, snrs_scale[:,2], label = 'Glc', color = glucose_colour, alpha = fill_alpha*3)

    axs[1,1].legend()

    axs[0,0].set_title('Water subbed spectra')

    axs[0,1].set_title('Scaled + shifted spectra')

    axs[0,0].set_xlabel('ppm')
    axs[0,1].set_xlabel('ppm')
    axs[1,0].set_xlabel('nscans')
    axs[1,1].set_xlabel('nscans')

    axs[1,0].set_ylabel('SNR')
    axs[1,1].set_ylabel('SNR')


    llim00, ulim00 = axs[0,0].get_ylim()
    llim01, ulim01 = axs[0,1].get_ylim()

    axs[0,0].fill_betweenx([llim00, ulim00], citric_acid_bounds[0], citric_acid_bounds[1], color= citrate_colour, alpha=fill_alpha)
    axs[0,1].fill_betweenx([llim01, ulim01], citric_acid_bounds[0], citric_acid_bounds[1], color= citrate_colour, alpha=fill_alpha)

    axs[0,0].fill_betweenx([llim00, ulim00], lactic_acid_bounds[0], lactic_acid_bounds[1], color= lactate_colour, alpha=fill_alpha)
    axs[0,1].fill_betweenx([llim01, ulim01], lactic_acid_bounds[0], lactic_acid_bounds[1], color= lactate_colour, alpha=fill_alpha)

    axs[0,0].fill_betweenx([llim00, ulim00], glucose_bounds[0], glucose_bounds[1], color=glucose_colour, alpha=fill_alpha)
    axs[0,1].fill_betweenx([llim01, ulim01], glucose_bounds[0], glucose_bounds[1], color=glucose_colour, alpha=fill_alpha)

    axs[0,0].fill_betweenx([llim00, ulim00], noise_bounds[0], noise_bounds[1], color= noise_colour, alpha=fill_alpha)
    axs[0,1].fill_betweenx([llim01, ulim01], noise_bounds[0], noise_bounds[1], color= noise_colour, alpha=fill_alpha)


    plt.show()

def big_plot_scaled(scan = 256):
    scans = [1, 2, 4, 8, 16, 32, 64, 128, 256]
    sample = 'D22'
    times = []

    for scan in scans:
        S = LoadSpectra()
        S.ReadRawData(nscan = scan, 
                    sample = sample,
                    pulse = 'zg30')
        tt = S.time_taken
        times.append(tt)

    S = LoadSpectra()
    S.ReadTextFile(nscan = 256, sample = sample)
    S.ReadRawData(sample = sample)
    S.SubtractWater()

    A = AnalyseSpectra()
    A.InputData(x = S.water_sub_ppm, y = S.water_sub_amplitude)
    A.FitTSP()

    ref_x, ref_y = A.tsp_centre, A.tsp_integral
    
    citric_acid_bounds = [2.37, 2.82]
    lactic_acid_bounds = [1.17, 1.50]
    glucose_bounds     = [3.19, 3.98]
    noise_bounds       = [-2.0, -1.0]
    import matplotlib as mpl

    cmap = mpl.colormaps['BrBG']

    colours = cmap(np.linspace(0, 1, len(scans) + 2))

    X_ORIG = []
    Y_ORIG = []

    X_SCALED = []
    Y_SCALED = []

    SNR_ORIG = []
    SNR_SCALED = []

    SIG_ORIG = []
    SIG_SCALED = []


    for scan in scans:
        print(scan)
        S = LoadSpectra()
        S.ReadTextFile(nscan = scan, sample = sample)
        S.SubtractWater()
        X_ORIG.append(S.water_sub_ppm)
        Y_ORIG.append(S.water_sub_amplitude)

        A = AnalyseSpectra()
        A.InputData(x = S.water_sub_ppm, y = S.water_sub_amplitude)
        
        A.SignalToNoise(signal_bounds = citric_acid_bounds)
        citric_acid_snr = A.snr
        citric_acid_signal = A.spectra_signal

        A.SignalToNoise(signal_bounds = lactic_acid_bounds)
        lactic_acid_snr = A.snr
        lactic_acid_signal = A.spectra_signal

        A.SignalToNoise(signal_bounds = glucose_bounds)
        glucose_snr = A.snr
        glucose_signal = A.spectra_signal
        noise_signal = A.spectra_noise
        
        snrs = [citric_acid_snr, lactic_acid_snr, glucose_snr]
        signals = [citric_acid_signal, lactic_acid_signal, glucose_signal, noise_signal]

        SNR_ORIG.append(snrs)
        SIG_ORIG.append(signals)

        A.FitTSP(plot_fit = True,  save_fig = True, plot_fnm = 'tsp_%s.png' % scan)

        tsp_x, tsp_y = A.tsp_centre, A.tsp_integral

        y_scale = ref_y/tsp_y
        x_shift = ref_x - tsp_x
        A.ScaleSpectra(y_scaling = y_scale, x_shift = x_shift)

        A.SignalToNoise(signal_bounds = citric_acid_bounds)
        citric_acid_snr = A.snr
        citric_acid_signal = A.spectra_signal


        A.SignalToNoise(signal_bounds = lactic_acid_bounds)
        lactic_acid_snr = A.snr
        lactic_acid_signal = A.spectra_signal


        A.SignalToNoise(signal_bounds = glucose_bounds)
        glucose_snr = A.snr
        glucose_signal = A.spectra_signal
        noise_signal = A.spectra_noise


        snrs = [citric_acid_snr, lactic_acid_snr, glucose_snr]
        signals = [citric_acid_signal, lactic_acid_signal, glucose_signal, noise_signal]


        X_SCALED.append(A.processed_ppm)
        Y_SCALED.append(A.processed_amplitude)
        SNR_SCALED.append(snrs)
        SIG_SCALED.append(signals)

    snrs_orig = np.array(SNR_ORIG)
    snrs_scale = np.array(SNR_SCALED)

    sigs_orig = np.array(SIG_ORIG)
    sigs_scale = np.array(SIG_SCALED)



    glucose_colour = 'orange'
    citrate_colour = 'green'
    lactate_colour = 'blue'
    noise_colour   = 'black'
    fill_alpha         = 0.1

    custom_lines = []
    legend_labels = []

    # Create figure and gridspec
    fig = plt.figure(figsize=(10, 6))
    gs = fig.add_gridspec(2, 2)

    # Top panel spanning two columns
    ax_top = fig.add_subplot(gs[0, :])
    for i in range(0, len(scans)):
        ax_top.plot(X_SCALED[i], Y_SCALED[i], label = 'nscan = %s' % scans[i], lw = 0.5, color = colours[i + 1])
        custom_lines.append(Line2D([0], [0], color=colours[i+1], linestyle='-'))
        legend_labels.append('ns = %s' % scans[i])

    llim, ulim = ax_top.get_ylim()
    ax_top.fill_betweenx([llim, ulim], citric_acid_bounds[0], citric_acid_bounds[1], color= citrate_colour, alpha=fill_alpha)
    ax_top.fill_betweenx([llim, ulim], lactic_acid_bounds[0], lactic_acid_bounds[1], color= lactate_colour, alpha=fill_alpha)
    ax_top.fill_betweenx([llim, ulim], glucose_bounds[0], glucose_bounds[1], color=glucose_colour, alpha=fill_alpha)
    ax_top.fill_betweenx([llim, ulim], noise_bounds[0], noise_bounds[1], color= noise_colour, alpha=fill_alpha)

    ax_top.set_xlim([-0.25, 4.4])
    ax_top.set_ylim([-0.25 * 1e9, 1.56 *1e9])
    ax_top.set_xlabel('ppm')
    ax_top.invert_xaxis()
    ax_top.legend(custom_lines, legend_labels, frameon = False, loc = 'upper right', fontsize = 7)

    # Bottom left panel
    ax_bottom_left = fig.add_subplot(gs[1, 0])

    #ax_bottom_left.plot(scans, np.ones(len(scans))*3., color = 'k', lw = 0.7, alpha = 0.5, ls = '--')
    #ax_bottom_left.plot(scans, sigs_orig[:,0]/sigs_orig[:,0][-1], color = citrate_colour, alpha = fill_alpha*8, linestyle = '-')
    #ax_bottom_left.plot(scans, sigs_orig[:,1]/sigs_orig[:,1][-1], color = lactate_colour, alpha = fill_alpha*8, linestyle = '-')
    #ax_bottom_left.plot(scans, sigs_orig[:,2]/sigs_orig[:,2][-1], color = glucose_colour, alpha = fill_alpha*8, linestyle = '-')

    ax_bottom_left.scatter(scans, sigs_scale[:,0], color = citrate_colour, alpha = fill_alpha*8, marker = 'x', label = 'Cit')
    ax_bottom_left.scatter(scans, sigs_scale[:,1], color = lactate_colour, alpha = fill_alpha*8, marker = 'x', label = 'Lac')
    ax_bottom_left.scatter(scans, sigs_scale[:,2], color = glucose_colour, alpha = fill_alpha*8, marker = 'x', label = 'Glc')
    ax_bottom_left.scatter(scans, sigs_scale[:,3], color = 'k',            alpha = fill_alpha*8, marker = 'x', label = 'Noise')
    ax_bottom_left.set_xscale('log')

    custom_lines2 = [Line2D([0], [0], color=citrate_colour, linestyle='-'),
                    Line2D([0], [0], color=lactate_colour, linestyle='-'),
                    Line2D([0], [0], color=glucose_colour, linestyle='-'),
                    Line2D([0], [0], color='k', linestyle='-')]

    ax_bottom_left.legend(custom_lines2, ['Cit', 'Lac', 'Glc', 'Noise'], frameon = False, fontsize = 7)

    ax_bottom_left.set_xlabel('Scans')
    ax_bottom_left.set_ylabel('Re-scaled spectra integral')




    # Bottom right panel
    ax_bottom_right = fig.add_subplot(gs[1, 1])

    ax_bottom_right.plot(times, np.ones(len(times))*3., color = 'k', lw = 0.7, alpha = 0.5, ls = '--')
    ax_bottom_right.plot(times, snrs_scale[:,0], color = citrate_colour, alpha = fill_alpha*8)
    ax_bottom_right.plot(times, snrs_scale[:,1], color = lactate_colour, alpha = fill_alpha*8)
    ax_bottom_right.plot(times, snrs_scale[:,2], color = glucose_colour, alpha = fill_alpha*8)
    ax_bottom_right.set_xlabel('Time (secs)')
    ax_bottom_right.set_ylabel('SNR')
    ax_bottom_right.set_yscale('log')
    ax_bottom_right.set_xscale('log')

    x_pos = sum(glucose_bounds)/2.
    y_pos = 1.3*1e9
    ax_top.annotate('Glucose\n5mM', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = glucose_colour)

    x_pos = sum(citric_acid_bounds)/2.
    y_pos = 1.3*1e9
    ax_top.annotate('Citrate\n0.2mM', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = citrate_colour)


    x_pos = sum(lactic_acid_bounds)/2.
    y_pos = 1.3*1e9
    ax_top.annotate('Lactate\n2mM', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = lactate_colour)

    x_pos = sum(noise_bounds)/2.
    y_pos = 0.075*1e9
    ax_top.annotate('Noise', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = noise_colour)


    my_axs = [ax_top, ax_bottom_left, ax_bottom_right]

    for ax in my_axs:
        ax.tick_params(axis='x', direction='in', which='both')
        ax.tick_params(axis='y', direction='in', which='both')

    # Adjust layout
    fig.tight_layout()

    # Show plot

    #plt.savefig('big_fig_scaled_integral.pdf')
    plt.show()
    plt.close()

    #'matplotlib-label-line'

#big_plot_scaled()

def d22_vs_d4(scan = '256'):
    sample = 'D22'
    pulse = 'zg30'


    D22 = LoadSpectra()
    D22.ReadTextFile(nscan = scan, 
                sample = 'D22',
                pulse = pulse)
    D22.SubtractWater()
    A22 = AnalyseSpectra()
    A22.InputData(x = D22.water_sub_ppm, y = D22.water_sub_amplitude)
    A22.FitTSP()



    D24 = LoadSpectra()
    D24.ReadTextFile(nscan = scan, 
                sample = 'D24',
                pulse = pulse)
    D24.SubtractWater()
    A24 = AnalyseSpectra()
    A24.InputData(x = D24.water_sub_ppm, y = D24.water_sub_amplitude)
    A24.FitTSP()

    glucose_bounds     = [3.0, 4.2]
    lactic_acid_bounds = [1.2, 1.45]



    glucose_mask = (D24.initial_ppm > glucose_bounds[0]) * (D24.initial_ppm < glucose_bounds[1]) 
    lactic_acid_mask = (D24.initial_ppm > lactic_acid_bounds[0]) * (D24.initial_ppm < lactic_acid_bounds[1]) 
    
    glc_ylims = [np.nanmin(D24.initial_amplitude[glucose_mask])  ,  np.nanmax(D24.initial_amplitude[glucose_mask])*1.3]
    lac_ylims = [np.nanmin(D24.initial_amplitude[lactic_acid_mask])  ,  np.nanmax(D24.initial_amplitude[lactic_acid_mask])*1.3]

    D24.initial_amplitude


    fig, axs = plt.subplots(1,2, figsize = [12, 5])

    axs[0].plot(D22.water_sub_ppm, D22.water_sub_amplitude, label = 'D22', color = 'black')
    axs[0].plot(D24.water_sub_ppm, D24.water_sub_amplitude, label = 'D24', color = 'red')
    axs[0].set_xlabel('ppm')
    axs[0].set_ylabel('amplitude')
    axs[0].legend()
    axs[0].invert_xaxis()
    glucose_xlims     = [2.3, 4.5]
    axs[0].set_xlim(glucose_xlims)
    axs[0].set_ylim(glc_ylims)

    axs[1].plot(D22.water_sub_ppm, D22.water_sub_amplitude, label = 'water_sub data', color = 'black')
    axs[1].plot(D24.water_sub_ppm, D24.water_sub_amplitude, label = 'water_sub data', color = 'red')
    axs[1].set_xlabel('ppm')
    axs[1].set_ylabel('amplitude')
    axs[1].invert_xaxis()
    lactic_acid_xlims = [1., 1.65]

    axs[1].set_xlim(lactic_acid_xlims)
    axs[1].set_ylim(lac_ylims)

    plt.savefig('')

    '''
    fig, axs = plt.subplots(1,2, figsize = [12, 5])
    axs[0].plot(D22.initial_ppm, D22.initial_amplitude, label = 'Initial data')
    axs[0].plot(D22.water_sub_ppm, D22.water_sub_amplitude, label = 'water_sub data')
    if len(D22.fit != 0):
        axs[0].plot(D22.water_sub_ppm, D22.fit, label = 'Fit', ls = '--')
    axs[0].set_xlabel('ppm')
    axs[0].set_ylabel('amplitude')
    axs[0].legend()
    axs[0].invert_xaxis()

    axs[1].plot(D24.initial_ppm, D24.initial_amplitude, label = 'Initial data')
    axs[1].plot(D24.water_sub_ppm, D24.water_sub_amplitude, label = 'water_sub data')
    if len(D24.fit != 0):
        axs[1].plot(D24.water_sub_ppm, D24.fit, label = 'Fit', ls = '--')
    axs[1].set_xlabel('ppm')
    axs[1].set_ylabel('amplitude')
    axs[1].legend()
    axs[1].invert_xaxis()
    plt.show()
    '''

scans = [1, 2, 4, 8, 16, 32, 64, 128, 256]
citric_acid_bounds = [2.37, 2.82]
lactic_acid_bounds = [1.17, 1.50]
glucose_bounds     = [3.19, 3.98]
noise_bounds       = [-2.0, -1.0]

def getim(scan, sample, pulse):
    S = LoadSpectra()
    S.ReadTextFile(nscan = scan, 
                sample = sample,
                pulse = pulse)
    S.SubtractWater()
    A = AnalyseSpectra()
    A.InputData(x = S.water_sub_ppm, y = S.water_sub_amplitude)

    A.SignalToNoise(signal_bounds = glucose_bounds)
    glc_snr = A.snr
    glc_signal = A.spectra_signal
    noise_signal = A.spectra_noise

    A.SignalToNoise(signal_bounds = lactic_acid_bounds)
    lac_snr = A.snr
    lac_signal = A.spectra_signal

    A.SignalToNoise(signal_bounds = citric_acid_bounds)
    cit_snr = A.snr
    cit_signal = A.spectra_signal

    return glc_snr, glc_signal, lac_snr, lac_signal, cit_snr, cit_signal,  noise_signal

def new():

    D22 = []
    D24 = []
    pulse = 'zg30'

    for scan in scans:
        d22_data = getim(scan, 'D22', pulse)
        d24_data = getim(scan, 'D24', pulse)
        D22.append(d22_data)
        D24.append(d24_data)

    D22 = np.array(D22)
    D24 = np.array(D24)


    glucose_colour = 'orange'
    citrate_colour = 'green'
    lactate_colour = 'blue'
    noise_colour   = 'black'

    fill_alpha         = 0.1

    custom_lines = [Line2D([0], [0], color=glucose_colour, linestyle='-'),
                    Line2D([0], [0], color=citrate_colour, linestyle='-'),
                    Line2D([0], [0], color=lactate_colour, linestyle='-'),
                    Line2D([0], [0], color='k', linestyle='-'),
                    Line2D([0], [0], color='k', linestyle='--')]



    fig, axs = plt.subplots(1,1, figsize = [6,6])

    axs.plot(scans, D22[:,0], color = glucose_colour, alpha = fill_alpha*8, linestyle = '-')
    axs.plot(scans, D24[:,0], color = glucose_colour, alpha = fill_alpha*8, linestyle = '--')

    axs.plot(scans, D22[:,2], color = lactate_colour, alpha = fill_alpha*8, linestyle = '-')
    axs.plot(scans, D24[:,2], color = lactate_colour, alpha = fill_alpha*8, linestyle = '--')

    axs.plot(scans, D22[:,4], color = citrate_colour, alpha = fill_alpha*8, linestyle = '-')
    axs.plot(scans, D24[:,4], color = citrate_colour, alpha = fill_alpha*8, linestyle = '--')

    axs.legend(custom_lines, ['Glc', 'Lac', 'Cit', 'D22', 'D24'])

    axs.set_xlabel('Scans')
    axs.set_ylabel('SNR')
    plt.show()

def combo():
    scans = [1, 2, 4, 8, 16, 32, 64, 128, 256]
    sample = 'D24'
    times = []

    for scan in scans:
        S = LoadSpectra()
        S.ReadRawData(nscan = scan, 
                    sample = sample,
                    pulse = 'zg30')
        tt = S.time_taken
        times.append(tt)

    S = LoadSpectra()
    S.ReadTextFile(sample = sample, nscan = 256)
    S.ReadRawData(sample=sample)
    S.SubtractWater()

    A = AnalyseSpectra()
    A.InputData(x = S.water_sub_ppm, y = S.water_sub_amplitude)
    A.FitTSP()

    ref_x, ref_y = A.tsp_centre, A.tsp_integral

    citric_acid_bounds = [2.37, 2.82]
    lactic_acid_bounds = [1.17, 1.50]
    glucose_bounds     = [3.19, 3.98]
    noise_bounds       = [-2.0, -1.0]
    import matplotlib as mpl

    cmap = mpl.colormaps['BrBG']

    colours = cmap(np.linspace(0, 1, len(scans) + 2))

    X_ORIG = []
    Y_ORIG = []

    X_SCALED = []
    Y_SCALED = []

    SNR_ORIG = []
    SNR_SCALED = []


    for scan in scans:
        print(scan)
        S = LoadSpectra()
        S.ReadTextFile(nscan = scan, sample =sample)
        S.SubtractWater()
        X_ORIG.append(S.water_sub_ppm)
        Y_ORIG.append(S.water_sub_amplitude)

        A = AnalyseSpectra()
        A.InputData(x = S.water_sub_ppm, y = S.water_sub_amplitude)
        
        A.SignalToNoise(signal_bounds = citric_acid_bounds)
        citric_acid_snr = A.snr

        A.SignalToNoise(signal_bounds = lactic_acid_bounds)
        lactic_acid_snr = A.snr

        A.SignalToNoise(signal_bounds = glucose_bounds)
        glucose_snr = A.snr

        snrs = [citric_acid_snr, lactic_acid_snr, glucose_snr]

        SNR_ORIG.append(snrs)


        A.FitTSP()
        tsp_x, tsp_y = A.tsp_centre, A.tsp_integral
        y_scale = ref_y/tsp_y
        x_shift = ref_x - tsp_x
        A.ScaleSpectra(y_scaling = y_scale, x_shift = x_shift)

        A.SignalToNoise(signal_bounds = citric_acid_bounds)
        citric_acid_snr = A.snr

        A.SignalToNoise(signal_bounds = lactic_acid_bounds)
        lactic_acid_snr = A.snr

        A.SignalToNoise(signal_bounds = glucose_bounds)
        glucose_snr = A.snr

        snrs = [citric_acid_snr, lactic_acid_snr, glucose_snr]


        X_SCALED.append(A.processed_ppm)
        Y_SCALED.append(A.processed_amplitude)
        SNR_SCALED.append(snrs)

    snrs_orig = np.array(SNR_ORIG)
    snrs_scale = np.array(SNR_SCALED)


    glucose_colour = 'orange'
    citrate_colour = 'green'
    lactate_colour = 'blue'
    noise_colour   = 'black'
    fill_alpha         = 0.1


    custom_lines = []
    legend_labels = []

    # Create figure and gridspec
    fig = plt.figure(figsize=(10, 6))
    gs = fig.add_gridspec(2, 2)

    # Top panel spanning two columns
    ax_top = fig.add_subplot(gs[0, :])
    for i in range(0, len(scans)):
        ax_top.plot(X_ORIG[i], Y_ORIG[i], label = 'nscan = %s' % scans[i], lw = 0.5, color = colours[i + 1])
        custom_lines.append(Line2D([0], [0], color=colours[i+1], linestyle='-'))
        legend_labels.append('ns = %s' % scans[i])

    llim, ulim = ax_top.get_ylim()
    ax_top.fill_betweenx([llim, ulim], citric_acid_bounds[0], citric_acid_bounds[1], color= citrate_colour, alpha=fill_alpha)
    ax_top.fill_betweenx([llim, ulim], lactic_acid_bounds[0], lactic_acid_bounds[1], color= lactate_colour, alpha=fill_alpha)
    ax_top.fill_betweenx([llim, ulim], glucose_bounds[0], glucose_bounds[1], color=glucose_colour, alpha=fill_alpha)
    ax_top.fill_betweenx([llim, ulim], noise_bounds[0], noise_bounds[1], color= noise_colour, alpha=fill_alpha)

    ax_top.set_xlim([-2.1, 5.5])
    ax_top.set_ylim([-0.5 * 1e8, 1.6 *1e9])
    ax_top.set_xlabel('ppm')
    ax_top.invert_xaxis()
    ax_top.legend(custom_lines, legend_labels, frameon = False)

    # Bottom left panel
    ax_bottom_right = fig.add_subplot(gs[1, 1])
    ax_bottom_right.plot(scans, np.ones(len(scans))*3., color = 'k', lw = 0.7, alpha = 0.5, ls = '--')
    ax_bottom_right.scatter(scans, snrs_orig[:,2], color = glucose_colour, alpha = fill_alpha*8, marker = 'x', label = 'Glucose')
    ax_bottom_right.scatter(scans, snrs_orig[:,1], color = lactate_colour, alpha = fill_alpha*8, marker = 'x', label = 'Lactate')
    ax_bottom_right.scatter(scans, snrs_orig[:,0], color = citrate_colour, alpha = fill_alpha*8, marker = 'x', label = 'Citrate')
    ax_bottom_right.plot(scans, snrs_orig[:,2], color = glucose_colour, alpha = fill_alpha*8, lw = 0.3)
    ax_bottom_right.plot(scans, snrs_orig[:,1], color = lactate_colour, alpha = fill_alpha*8, lw = 0.3)
    ax_bottom_right.plot(scans, snrs_orig[:,0], color = citrate_colour, alpha = fill_alpha*8, lw = 0.3)
    ax_bottom_right.set_xscale('log')
    ax_bottom_right.set_xlabel('ns')
    ax_bottom_right.set_ylabel('SNR')
    ax_bottom_right.legend(frameon = False, loc = 'upper left')



    # Bottom right panel
    ax_bottom_left = fig.add_subplot(gs[1, 0])
    x = X_SCALED[0]
    for i in range(0, len(x), 10):
        ax_bottom_left.axvline(x[i], color='lightgray', linewidth=0.5, alpha=0.5)    

    for i in range(0, len(scans)):
        ax_bottom_left.plot(X_SCALED[i], Y_SCALED[i], lw = 1.0, color = colours[i + 1])

    ax_bottom_left.set_xlabel('ppm')
    ax_bottom_left.set_xlim([lactic_acid_bounds[0] - 0.1, lactic_acid_bounds[1] + 0.1])
    ax_bottom_left.fill_betweenx([llim, ulim], lactic_acid_bounds[0], lactic_acid_bounds[1], color= lactate_colour, alpha=fill_alpha)
    ax_bottom_left.set_ylim([-0.9 * 1e8, 0.8 *1e9])
    ax_bottom_left.invert_xaxis()
    #ax_bottom_right.set_ylabel('I')

    #ax_bottom_right.plot(times, np.ones(len(times))*3., color = 'k', lw = 0.7, alpha = 0.5, ls = '--')
    #ax_bottom_right.plot(times, snrs_orig[:,0], color = citrate_colour, alpha = fill_alpha*8)
    #ax_bottom_right.plot(times, snrs_orig[:,1], color = lactate_colour, alpha = fill_alpha*8)
    #ax_bottom_right.plot(times, snrs_orig[:,2], color = glucose_colour, alpha = fill_alpha*8)
    #ax_bottom_right.set_xlabel('Time (secs)')
    #ax_bottom_right.set_ylabel('SNR')
    #ax_bottom_right.set_yscale('log')
    #ax_bottom_right.set_xscale('log')

    x_pos = sum(glucose_bounds)/2.
    y_pos = 1.3*1e9
    ax_top.annotate('Glucose\n10mM', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = glucose_colour)

    x_pos = sum(citric_acid_bounds)/2.
    y_pos = 1.3*1e9
    ax_top.annotate('Citrate\n0.2mM', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = citrate_colour)


    x_pos = sum(lactic_acid_bounds)/2.
    y_pos = 1.3*1e9
    ax_top.annotate('Lactate\n2mM', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = lactate_colour)

    x_pos = sum(noise_bounds)/2.
    y_pos = 0.075*1e9
    ax_top.annotate('Noise', xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = noise_colour)


    my_axs = [ax_top, ax_bottom_left, ax_bottom_right]

    for ax in my_axs:
        ax.tick_params(axis='x', direction='in', which='both')
        ax.tick_params(axis='y', direction='in', which='both')

    # Adjust layout
    fig.tight_layout()

    # Show plot

    plt.savefig('Paper_Figs/fig1_lactate.pdf')
    plt.close()

    #'matplotlib-label-line'

def mega_combo():
    scans = [1, 2, 4, 8, 16, 32, 64, 128, 256]
    sample = 'D24'
    snr_choice = 'brk'
    times = []

    for scan in scans:
        S = LoadSpectra()
        S.ReadRawData(nscan = scan, 
                    sample = sample,
                    pulse = 'zg30')
        tt = S.time_taken
        times.append(tt)

    S = LoadSpectra()
    S.ReadTextFile(sample = sample, nscan = 256)
    S.ReadRawData(sample=sample)
    S.SubtractWater()

    A = AnalyseSpectra()
    A.InputData(x = S.water_sub_ppm, y = S.water_sub_amplitude)
    A.FitTSP()

    ref_x, ref_y = A.tsp_centre, A.tsp_integral
    
    citric_acid_bounds = [2.37, 2.82]
    lactic_acid_bounds = [1.17, 1.50]
    glucose_bounds     = [3.19, 3.98]
    #citric_acid_bounds = [2.5, 2.7]
    #lactic_acid_bounds = [1.2, 1.45]
    #glucose_bounds     = [3.0, 4.2]
    noise_bounds       = [-2.0, -1.0]
    import matplotlib as mpl

    cmap = mpl.colormaps['BrBG']

    colours = cmap(np.linspace(0, 1, len(scans) + 2))

    X_ORIG = []
    Y_ORIG = []

    X_SCALED = []
    Y_SCALED = []

    SNR_ORIG = []
    SNR_SCALED = []


    for scan in scans:
        print(scan)
        S = LoadSpectra()
        S.ReadTextFile(nscan = scan, sample =sample)
        S.SubtractWater()
        X_ORIG.append(S.water_sub_ppm)
        Y_ORIG.append(S.water_sub_amplitude)

        A = AnalyseSpectra()
        A.InputData(x = S.water_sub_ppm, y = S.water_sub_amplitude)
        
        A.SignalToNoise(signal_bounds = citric_acid_bounds, snr_choice = snr_choice)
        citric_acid_snr = A.snr

        A.SignalToNoise(signal_bounds = lactic_acid_bounds, snr_choice = snr_choice)
        lactic_acid_snr = A.snr

        A.SignalToNoise(signal_bounds = glucose_bounds, snr_choice = snr_choice)
        glucose_snr = A.snr

        snrs = [citric_acid_snr, lactic_acid_snr, glucose_snr]

        SNR_ORIG.append(snrs)


        A.FitTSP()
        tsp_x, tsp_y = A.tsp_centre, A.tsp_integral
        y_scale = ref_y/tsp_y
        x_shift = ref_x - tsp_x
        A.ScaleSpectra(y_scaling = y_scale, x_shift = x_shift)

        A.SignalToNoise(signal_bounds = citric_acid_bounds, snr_choice = snr_choice)
        citric_acid_snr = A.snr

        A.SignalToNoise(signal_bounds = lactic_acid_bounds, snr_choice = snr_choice)
        lactic_acid_snr = A.snr

        A.SignalToNoise(signal_bounds = glucose_bounds, snr_choice = snr_choice)
        glucose_snr = A.snr

        snrs = [citric_acid_snr, lactic_acid_snr, glucose_snr]


        X_SCALED.append(A.processed_ppm)
        Y_SCALED.append(A.processed_amplitude)
        SNR_SCALED.append(snrs)

    snrs_orig = np.array(SNR_ORIG)
    snrs_scale = np.array(SNR_SCALED)

    glucose_colour = 'orange'
    citrate_colour = 'green'
    lactate_colour = 'blue'
    noise_colour   = 'black'
    fill_alpha         = 0.07


    custom_lines = []
    legend_labels = []

    # Enable LaTeX rendering
    matplotlib.rc('text', usetex=True)

    # Set the default font to Times New Roman
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.size'] = 10


    # Create figure and gridspec
    fig = plt.figure(figsize=(10, 8))
    gs = fig.add_gridspec(3, 2, width_ratios = [1.0, 1.0])

    # Top panel spanning two columns
    ax_top = fig.add_subplot(gs[0, :])
    for i in range(0, len(scans)):
        ax_top.plot(X_ORIG[i], Y_ORIG[i], label = 'nscan = %s' % scans[i], lw = 0.5, color = colours[i + 1])
        custom_lines.append(Line2D([0], [0], color=colours[i+1], linestyle='-'))
        legend_labels.append('ns = %s' % scans[i])

    llim, ulim = ax_top.get_ylim()
    ax_top.fill_betweenx([llim, ulim], citric_acid_bounds[0], citric_acid_bounds[1], color= citrate_colour, alpha=fill_alpha)
    ax_top.fill_betweenx([llim, ulim], lactic_acid_bounds[0], lactic_acid_bounds[1], color= lactate_colour, alpha=fill_alpha)
    ax_top.fill_betweenx([llim, ulim], glucose_bounds[0], glucose_bounds[1], color=glucose_colour, alpha=fill_alpha)
    ax_top.fill_betweenx([llim, ulim], noise_bounds[0], noise_bounds[1], color= noise_colour, alpha=fill_alpha)
    ax_top.axvline(-0.05, color='lightgray', linewidth=0.5, alpha=0.5)    
    ax_top.axvline(0.05, color='lightgray', linewidth=0.5, alpha=0.5)    

    ax_top.set_xlim([-2.1, 5.5])
    ax_top.set_ylim([-0.5 * 1e8, 1.6 *1e9])
    ax_top.set_xlabel('$\delta$ [ppm]')
    ax_top.invert_xaxis()
    #ax_top.legend(custom_lines, legend_labels, frameon = False, fontsize = 8, ncol = 3, loc = 7)


    # MID LEFT panel
    ax_mid_left = fig.add_subplot(gs[1, 0])
    x = X_SCALED[0]
    for i in range(0, len(x), 100):
        ax_mid_left.axvline(x[i], color='lightgray', linewidth=0.5, alpha=0.5)    

    for i in range(0, len(scans)):
        ax_mid_left.plot(X_SCALED[i], Y_SCALED[i], lw = 1.0, color = colours[i + 1])

    ax_mid_left.set_xlabel('$\delta$ [ppm]')
    ax_mid_left.set_xlim([glucose_bounds[0] - 0.1, glucose_bounds[1] + 0.1])
    ax_mid_left.fill_betweenx([llim, ulim], glucose_bounds[0], glucose_bounds[1], color= glucose_colour, alpha=fill_alpha)
    ax_mid_left.set_ylim([-0.9 * 1e8, 1.4*1e9])
    ax_mid_left.invert_xaxis()


    # MID RIGHT panel
    ax_mid_right = fig.add_subplot(gs[1, 1])
    x = X_SCALED[0]
    for i in range(0, len(x), 10):
        ax_mid_right.axvline(x[i], color='lightgray', linewidth=0.5, alpha=0.5)    

    for i in range(0, len(scans)):
        ax_mid_right.plot(X_SCALED[i], Y_SCALED[i], lw = 1.0, color = colours[i + 1])

    ax_mid_right.set_xlabel('$\delta$ [ppm]')
    ax_mid_right.set_xlim([citric_acid_bounds[0] - 0.1, citric_acid_bounds[1] + 0.1])
    ax_mid_right.fill_betweenx([llim, ulim], citric_acid_bounds[0], citric_acid_bounds[1], color= citrate_colour, alpha=fill_alpha)
    ax_mid_right.set_ylim([-0.9 * 1e8, 0.15*1e9])
    ax_mid_right.invert_xaxis()
    

    # Bottom left panel
    ax_bottom_left = fig.add_subplot(gs[2, 0])
    x = X_SCALED[0]
    for i in range(0, len(x), 10):
        ax_bottom_left.axvline(x[i], color='lightgray', linewidth=0.5, alpha=0.5)    

    for i in range(0, len(scans)):
        ax_bottom_left.plot(X_SCALED[i], Y_SCALED[i], lw = 1.0, color = colours[i + 1])

    ax_bottom_left.set_xlabel('$\delta$ [ppm]')
    ax_bottom_left.set_xlim([lactic_acid_bounds[0] - 0.1, lactic_acid_bounds[1] + 0.1])
    ax_bottom_left.fill_betweenx([llim, ulim], lactic_acid_bounds[0], lactic_acid_bounds[1], color= lactate_colour, alpha=fill_alpha)
    ax_bottom_left.set_ylim([-0.5 * 1e8, 0.7*1e9])
    ax_bottom_left.invert_xaxis()


    # Bottom right panel
    ax_bottom_right = fig.add_subplot(gs[2, 1])
    #ax_bottom_right.plot(scans, np.ones(len(scans))*3., color = 'k', lw = 0.5, alpha = 0.5, ls = '-')
    #ax_bottom_right.plot(scans, np.ones(len(scans))*10., color = 'k', lw = 0.5, alpha = 0.5, ls = '-')

    ax_bottom_right.scatter(scans, snrs_orig[:,2], color = glucose_colour, marker = 's', edgecolor = 'k', label = 'Glucose')#, alpha = fill_alpha*8)
    ax_bottom_right.scatter(scans, snrs_orig[:,1], color = lactate_colour, marker = 's', edgecolor = 'k', label = 'Lactate')#, alpha = fill_alpha*8)
    ax_bottom_right.scatter(scans, snrs_orig[:,0], color = citrate_colour, marker = 's', edgecolor = 'k', label = 'Citrate')#, alpha = fill_alpha*8)
    ax_bottom_right.plot(scans, snrs_orig[:,2], color = glucose_colour, lw = 0.3)#, alpha = fill_alpha*8)
    ax_bottom_right.plot(scans, snrs_orig[:,1], color = lactate_colour, lw = 0.3)#, alpha = fill_alpha*8)
    ax_bottom_right.plot(scans, snrs_orig[:,0], color = citrate_colour, lw = 0.3)#, alpha = fill_alpha*8)
    ax_bottom_right.set_xscale('log')
    ax_bottom_right.set_xlabel('$n_{\mathrm{s}}$')

    if snr_choice == 'brk':
        snr_lab = '$\mathrm{SNR}_{\mathrm{Brk.}}$'
    elif snr_choice == 'agi':
        snr_lab = '$\mathrm{SNR}_{\mathrm{Agi.}}$'
    elif snr_choice == 'liv':
        snr_lab = '$\mathrm{SNR}_{\mathrm{Int.}}$'

    ax_bottom_right.set_ylabel(snr_lab)
    ax_bottom_right.set_xlim([0.7, 335])



    custom_lines  = [(Line2D([0], [0], color=glucose_colour,markeredgecolor = 'k',  marker = 's', linestyle='-')),
                    (Line2D([0], [0], color=lactate_colour, markeredgecolor = 'k', marker = 's', linestyle='-')),
                    (Line2D([0], [0], color=citrate_colour, markeredgecolor = 'k', marker = 's', linestyle='-'))]

    ax_bottom_right.legend(custom_lines, ['Glucose', 'Lactate', 'Citrate'], frameon = False, loc = 'upper left')


    x_min, x_max = ax_bottom_right.get_xlim()
    y_min, y_max = ax_bottom_right.get_ylim()

    

    ax_bottom_right.fill_between([x_min, x_max], 0, 3, color= 'lightsalmon', alpha=0.2)
    ax_bottom_right.fill_between([x_min, x_max], 3, 10, color= 'navajowhite', alpha=0.2)
    ax_bottom_right.fill_between([x_min, x_max], 10, y_max*10, color= 'lightgreen', alpha=0.2)
    ax_bottom_right.set_ylim([0.9, y_max*1.2])
    ax_bottom_right.set_yscale('log')

    glc_label = '\\textbf{Glucose}\n \\textbf{10mM}'
    cit_label = '\\textbf{Citrate}\n \\textbf{0.2mM}'
    lac_label = '\\textbf{Lactate}\n \\textbf{2mM}'
    tsp_label = '\\textbf{TSP}\n \\textbf{XmM}'
    noise_label = '\\textbf{Noise}'

    x_pos_glc = sum(glucose_bounds)/2.
    y_pos_main = 1.3*1e9
    ax_top.annotate(glc_label, xy=(x_pos_glc, y_pos_main), xytext=(x_pos_glc, y_pos_main), horizontalalignment='center', color = glucose_colour)
    ax_top.annotate(tsp_label, xy=(0.0, y_pos_main), xytext=(0.0, y_pos_main), horizontalalignment='center', color = 'k')
    
    
    y_pos = 1.*1e9
    ax_mid_left.annotate(glc_label, xy=(x_pos_glc, y_pos), xytext=(x_pos_glc, y_pos), horizontalalignment='center', color = glucose_colour)
    x_pos_scaled = 4.12
    y_pos_scaled = y_pos + 0.25*1e9



    x_pos = sum(citric_acid_bounds)/2.
    y_pos = 1.3*1e9
    ax_top.annotate(cit_label, xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = citrate_colour)
    y_pos = 0.1*1e9
    x_pos = sum(citric_acid_bounds)/2. - 0.1
    ax_mid_right.annotate(cit_label, xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = citrate_colour)


    x_pos = sum(lactic_acid_bounds)/2.
    y_pos = 1.3*1e9
    ax_top.annotate(lac_label, xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = lactate_colour)
    y_pos = 0.56 *1e9
    ax_bottom_left.annotate(lac_label, xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = lactate_colour)

    x_pos = sum(noise_bounds)/2.
    y_pos = 0.075*1e9
    y_pos = 1.3*1e9
    ax_top.annotate(noise_label, xy=(x_pos, y_pos), xytext=(x_pos, y_pos), horizontalalignment='center', color = 'grey')


    my_axs = [ax_top, ax_mid_left, ax_mid_right, ax_bottom_left]
    scale_10 = 'Scaled\n 10 d.p. per $\\vert\\vert$'
    scale_100 = 'Scaled\n 100 d.p. per $\\vert\\vert$'
    inps =   ['Raw', scale_100, scale_10, scale_10]

    # Add a boxed annotation to each subplot
    for i in range(0, len(my_axs)):
        ax = my_axs[i]
        lab = inps[i]
        # Get the current axis limits
        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()
        
        # Calculate the x and y coordinates for the annotation
        if i == 0:
            x_ann = x_min + 0.02 * (x_max - x_min)
            y_ann = y_max - 0.11 * (y_max - y_min)
        else:
            x_ann = x_min + 0.05 * (x_max - x_min)
            y_ann = y_max - 0.2 * (y_max - y_min)
        
        # Add the boxed annotation
        ax.text(x_ann, y_ann, lab, bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'), fontsize=10)





    my_axs = [ax_top, ax_mid_left, ax_mid_right, ax_bottom_left, ax_bottom_right]

    for ax in my_axs:
        ax.tick_params(axis='x', direction='in', which='both')
        ax.tick_params(axis='y', direction='in', which='both')

    my_axs = [ax_top, ax_mid_left, ax_mid_right, ax_bottom_left]
    for ax in my_axs:
        ax.set_ylabel('$I(\delta)$')
    #plt.subplots_adjust(left=0.05, right=20)

    #fig.subplots_adjust(right=0.8)
    #cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])    # Adjust layout
    #fig.colorbar(im, cax=cbar_ax)
    #'matplotlib-label-line'
    import matplotlib.colorbar as colorbar
    import matplotlib.colors as clr

    cb_colors = []
    xtick_labs = []

    for i in range(0, len(scans)):
        cb_colors.append(colours[i + 1])
        xtick_labs.append("%s" % scans[i])

    num_colors = len(cb_colors)
    cmap_ = clr.ListedColormap(cb_colors)

    top_pos = ax_top.get_position()
    bottom_pos = ax_bottom_right.get_position()

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

    cb.set_ticks(range(num_colors))
    cb.ax.set_yticklabels(xtick_labs)
    cb.set_label('$n_{\mathrm{s}}$', fontsize = 13)
    
    plt.subplots_adjust(hspace = 0.3)
    
    
    #fig.tight_layout()

    # Show plot
    #plt.show()
    plt.savefig('Paper_Figs/test_%s.pdf' % snr_choice, bbox_inches = 'tight', pad_inches = 0.2 )
    plt.close()


mega_combo()
#big_plot()


