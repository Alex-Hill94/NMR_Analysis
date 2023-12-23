import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from pylab import Line2D

path = 'mag_strength_test.hdf5'

freqs = np.linspace(10.0, 800.0, 80)
idxs  = np.arange(0, len(freqs), 1)

u = h5.File(path, 'r')

x = []
y = []

custom_lines = [Line2D([0], [0], color='white', lw=4)]

for i in range(0, len(idxs)):
    idx = idxs[i]
    x = u['%s/x' % freqs[i]][()]
    y = u['%s/toty' % freqs[i]][()]

    fig, axs = plt.subplots(1,1, figsize = [5,5])
    axs.plot(x, y, color = 'red', lw = 1., ls = '-')
    axs.legend(custom_lines, ['%s MHz' % freqs[idx]], loc = 'upper left')
    axs.set_xlabel('ppm')
    axs.set_ylabel('amplitude')
    axs.set_xlim([2.2, 3.])
    axs.invert_xaxis()
    plt.savefig('Ims/%s.png' % i)
    plt.close()

u.close()
'''
path = 'mag_strength_test.hdf5'
u = h5.File(path, 'r')

x_80 = u['80.0/x'][()]
y_80 = u['80.0/toty'][()]

x_800 = u['800.0/x'][()]
y_800 = u['800.0/toty'][()]
u.close()

fig, axs = plt.subplots(1,1, figsize = [5,5])
axs.plot(x_80, y_80, color = 'red', lw = 1., ls = '-', label = '80 MHz')
axs.plot(x_800, y_800, color = 'blue', lw = 1., ls = '-', label = '800 MHz')
axs.set_xlabel('ppm')
axs.set_ylabel('amplitude')
axs.set_xlim([2.2, 3.])
axs.invert_xaxis()
plt.show()
'''