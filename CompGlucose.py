import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from pylab import Line2D
from DataGrab import *
path = 'mag_strength_test.hdf5'

freqs = np.linspace(10.0, 800.0, 80)
idxs  = np.arange(0, len(freqs), 1)

u = h5.File(path, 'r')

x = []
y = []

i = idxs[freqs == 80.0][0]

custom_lines = [(Line2D([0], [0], color='red', lw=4)),
                (Line2D([0], [0], color='blue', lw=4))]

x_sim = u['%s/x' % freqs[i]][()]
y_sim = u['%s/toty' % freqs[i]][()]

u.close()

x_exp, y_exp = grab_data(file_name = "Nov29-2023_D24_LC+AG.txt")

fig, axs = plt.subplots(1,1, figsize = [5,5])
axs.plot(x_sim, y_sim, color = 'red', lw = 1., ls = '-')
axs.plot(x_exp - 0.04, y_exp * 5.77/2.37, color = 'blue', lw = 1., ls = '-')
axs.legend(custom_lines, ['Simulation - %s MHz' % freqs[i], 'Experiment - %s MHz' % freqs[i]], loc = 'upper left')
axs.set_xlabel('ppm')
axs.set_ylabel('amplitude')
axs.set_xlim([3.2, 4.1])
axs.set_ylim([0.0, 7.1*1e8])
axs.set_title('Glucose')
axs.invert_xaxis()
plt.savefig('glucose.png')
plt.close()


