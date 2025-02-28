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

files = ['80MHzLactate_Oct24', '80MHzGlucose_Nov24', '80MHzCitrate_Nov24', '700MHzLactate', '700MHzGlucose', '700MHzCitrate']
root_path = '/Users/alexhill/Documents/UOL/Research/Companies/ViBo/Metabolomics/Data_Analysis/80MHzGlucose_Nov24/'

subdir = 'POS'
#filename_base = '2024-10-28-AH-LAC-01-10-zg30.txt'
exps = np.linspace(10,180,18 ).astype('int')
zg30_exps = exps[:9]
zg_exps   = exps[9:]
positions = np.arange(1, 13, 1)
pos = [f"{num:02}" for num in positions]
ns = [1, 2, 4, 8, 16, 32, 64, 128, 256]

ZG30 = []
ZG   = []

for i in range(0, len(positions)):
    path = root_path+subdir+str(positions[i])
    print(path)
    zg30_temp = []
    zg_temp = []
    for j in range(0, len(zg30_exps)):
        z3_e = zg30_exps[j]
        filename = '2024-11-26-AH-GLC-'+pos[i]+'_'+str(z3_e)+'_zg30.txt'
        print(filename)
        x, y = getData(path, filename)
        zg30_temp.append([x, y])
    for j in range(0, len(zg_exps)):
        z_e = zg_exps[j]
        filename = '2024-11-26-AH-GLC-'+pos[i]+'_'+str(z_e)+'_zg.txt'
        print(filename)
        x, y = getData(path, filename)
        zg_temp.append([x, y])
    ZG30.append(zg30_temp)
    ZG.append(zg_temp)

ZG30 = np.array(ZG30)
ZG = np.array(ZG)

data_file = 'Data.hdf5'
concs = ['1250uM', '2500uM', '5000uM', '10000uM']
'''
u = h5.File(data_file, 'r+')


for i in range(0, len(concs)):
    start = int(i*3)
    conc = concs[i]
    
    root = 'Glucose/80MHz/'+conc+'/'

    XA = ZG[start + 0][:, 0, :]
    YA = ZG[start + 0][:, 1, :]

    XB = ZG[start + 1][:, 0, :]
    YB = ZG[start + 1][:, 1, :]

    XC = ZG[start + 2][:, 0, :]
    YC = ZG[start + 2][:, 1, :]

    u.create_dataset(root +'zg/ppm_a', data =XA)
    u.create_dataset(root +'zg/amp_a', data = YA)

    u.create_dataset(root + 'zg/ppm_b', data = XB)
    u.create_dataset(root + 'zg/amp_b', data = YB)

    u.create_dataset(root + 'zg/ppm_c', data = XC)
    u.create_dataset(root + 'zg/amp_c', data = YC)

    u.create_dataset(root + 'zg/n_scans', data = ns)

    XA = ZG30[start + 0][:, 0, :]
    YA = ZG30[start + 0][:, 1, :]

    XB = ZG30[start + 1][:, 0, :]
    YB = ZG30[start + 1][:, 1, :]

    XC = ZG30[start + 2][:, 0, :]
    YC = ZG30[start + 2][:, 1, :]

    u.create_dataset(root +'zg30/ppm_a', data =XA)
    u.create_dataset(root +'zg30/amp_a', data = YA)

    u.create_dataset(root + 'zg30/ppm_b', data = XB)
    u.create_dataset(root + 'zg30/amp_b', data = YB)

    u.create_dataset(root + 'zg30/ppm_c', data = XC)
    u.create_dataset(root + 'zg30/amp_c', data = YC)

    u.create_dataset(root + 'zg30/n_scans', data = ns)

u.close()

'''
u = h5.File(data_file, 'r')

for conc in concs:
    amp_a = u['Glucose/80MHz/'+conc+'/zg/amp_a'][()]
    amp_b = u['Glucose/80MHz/'+conc+'/zg/amp_b'][()]
    amp_c = u['Glucose/80MHz/'+conc+'/zg/amp_c'][()]
    n_scans = u['Glucose/80MHz/'+conc+'/zg/n_scans'][()]
    ppm_a = u['Glucose/80MHz/'+conc+'/zg/ppm_a'][()]
    ppm_b = u['Glucose/80MHz/'+conc+'/zg/ppm_b'][()]
    ppm_c = u['Glucose/80MHz/'+conc+'/zg/ppm_c'][()]

    fig, axs = plt.subplots(1, 3)
    axs[0].set_title(conc)
    for i in range(0, len(amp_a)):
        axs[0].plot(ppm_a[i], amp_a[i])
        axs[1].plot(ppm_b[i], amp_b[i])
        axs[2].plot(ppm_c[i], amp_c[i])
    plt.show()
    plt.close()

u.close()
