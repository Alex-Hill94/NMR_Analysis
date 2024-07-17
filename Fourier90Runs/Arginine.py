from NMRClasses import *
import matplotlib
from bounds import *
import h5py as h5

S = LoadSpectra()
S.ReadTextFile(nscan = 256, 
            sample = 'A2000',
            pulse = 'zg30', 
            path = '/Users/alexhill/Desktop/Metabolomics/Data_Analysis/LiverpoolFilipData')


path = '/Users/alexhill/Software/ccpnmr3.2.0/bin/arginine_nonoise_80MHz.hdf5'

u = h5.File(path, 'r')
Y = u['80.0/toty'][()]
X = u['80.0/x'][()]
u.close()

fig, axs = plt.subplots(3,1, figsize = [10, 6])

axs[0].plot(S.initial_ppm, S.initial_amplitude, label = 'Initial data', color = 'blue')
axs[0].plot(X, Y, label = 'Simulations', color = 'red')

axs[1].plot(S.initial_ppm, S.initial_amplitude, label = 'Initial data', color = 'blue')

axs[2].plot(X, Y, label = 'Simulations', color = 'red')

axs[0].set_xlim([-0.1, 5.6])
axs[1].set_xlim([-0.1, 5.6])
axs[2].set_xlim([-0.1, 5.6])
axs[0].invert_xaxis()
axs[1].invert_xaxis()
axs[2].invert_xaxis()
axs[0].legend()
plt.show()