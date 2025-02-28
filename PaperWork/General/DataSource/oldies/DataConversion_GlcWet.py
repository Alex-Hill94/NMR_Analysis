import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '/Users/alexhill/Documents/GitHub/NMR_Analysis/Fourier90Runs')
from NMRClasses import *
import matplotlib as mpl
import os 
import csv
import pandas as pd 
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import h5py as h5
from scipy import stats

def getData(path, filename):
    #filename = 'AHill_241025_7_exp%s.txt' % exp
    S = LoadSpectra()
    S.ReadTextFile(path = path, filename = filename)
    return S.initial_ppm, S.initial_amplitude

root_path = '/Users/alexhill/Documents/UOL/Research/Companies/ViBo/Metabolomics/Data_Analysis/80MHzGlucoseWET_Feb25/'

subdir = 'POS'
#filename_base = '2024-10-28-AH-GLC-01-10-zg30.txt'
exps = np.linspace(10,90,9 ).astype('int')
wet_exps = exps

positions = np.arange(1, 13, 1)
pos = [f"{num:02}" for num in positions]
ns = [1, 2, 4, 8, 16, 32, 64, 128, 256]


#WET = np.zeros((12, 9, 2, 65536))

WET = []

for i in range(0, len(positions)):
    path = root_path+subdir+str(positions[i])
    print(path)
    wet_temp = []
    for j in range(0, len(wet_exps)):
        w_e = wet_exps[j]
        filename = '2025-02-05-AH-GLC-K_WET_f-'+pos[i]+'-'+str(w_e)+'.txt'
        print(filename)
        x, y = getData(path, filename)
        wet_temp.append([x, y])
    WET.append(wet_temp)
WET = np.array(WET)

data_file = 'DataWet.hdf5'
concs = ['1250uM', '2500uM', '5000uM', '10000uM']



u = h5.File(data_file, 'r+')


for i in range(0, len(concs)):
    start = int(i*3)
    conc = concs[i]
    
    root = 'Glucose/80MHz/'+conc+'/'

    XA = WET[start + 0][:, 0, :]
    YA = WET[start + 0][:, 1, :]

    XB = WET[start + 1][:, 0, :]
    YB = WET[start + 1][:, 1, :]

    XC = WET[start + 2][:, 0, :]
    YC = WET[start + 2][:, 1, :]

    u.create_dataset(root +'wet/ppm_a', data =XA)
    u.create_dataset(root +'wet/amp_a', data = YA)

    u.create_dataset(root + 'wet/ppm_b', data = XB)
    u.create_dataset(root + 'wet/amp_b', data = YB)

    u.create_dataset(root + 'wet/ppm_c', data = XC)
    u.create_dataset(root + 'wet/amp_c', data = YC)

    u.create_dataset(root + 'wet/n_scans', data = ns)


u.close()

u = h5.File(data_file, 'r')

for conc in concs:
    amp_a = u['Glucose/80MHz/'+conc+'/wet/amp_a'][()]
    amp_b = u['Glucose/80MHz/'+conc+'/wet/amp_b'][()]
    amp_c = u['Glucose/80MHz/'+conc+'/wet/amp_c'][()]
    n_scans = u['Glucose/80MHz/'+conc+'/wet/n_scans'][()]
    ppm_a = u['Glucose/80MHz/'+conc+'/wet/ppm_a'][()]
    ppm_b = u['Glucose/80MHz/'+conc+'/wet/ppm_b'][()]
    ppm_c = u['Glucose/80MHz/'+conc+'/wet/ppm_c'][()]

    fig, axs = plt.subplots(1, 3)
    axs[0].set_title(conc)
    for i in range(0, len(amp_a)):
        axs[0].plot(ppm_a[i], amp_a[i])
        axs[1].plot(ppm_b[i], amp_b[i])
        axs[2].plot(ppm_c[i], amp_c[i])
        axs[0].set_xlim([0,5])
        axs[1].set_xlim([0,5])
        axs[2].set_xlim([0,5])
    plt.show()
    plt.close()

u.close()
