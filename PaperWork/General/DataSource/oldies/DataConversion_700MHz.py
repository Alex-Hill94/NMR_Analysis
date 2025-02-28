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
lacate = True
glucose = False
citrate = False
if glucose:
    '''

    root_path = '/Users/alexhill/Documents/UOL/Research/Companies/ViBo/Metabolomics/Data_Analysis/700MHzGlucose/'
    pref = 'Ahill_241122_7_exp'
    suff = '0012.txt'
    ray = np.arange(1, 13, 1)
    path = root_path
    print(path)

    NOES = []
    for i in ray:
        filename = pref+str(i)+suff
        print(filename)
        x, y = getData(path, filename)
        NOES.append([x,y])


    NOES = np.array(NOES)

    data_file = 'Data.hdf5'
    concs = ['1250uM', '2500uM', '5000uM', '10000uM']

    u = h5.File(data_file, 'r+')

    for i in range(0, len(concs)):
        start = int(i*3)
        conc = concs[i]
        
        root = 'Glucose/700MHz/'+conc+'/'

        XA = NOES[start + 0][0]
        YA = NOES[start + 0][1]

        XB = NOES[start + 1][0]
        YB = NOES[start + 1][1]

        XC = NOES[start + 2][0]
        YC = NOES[start + 2][1]

        u.create_dataset(root +'noesy/ppm_a', data =XA)
        u.create_dataset(root +'noesy/amp_a', data = YA)

        u.create_dataset(root + 'noesy/ppm_b', data = XB)
        u.create_dataset(root + 'noesy/amp_b', data = YB)

        u.create_dataset(root + 'noesy/ppm_c', data = XC)
        u.create_dataset(root + 'noesy/amp_c', data = YC)

        u.create_dataset(root + 'noesy/n_scans', data = 32)

    u.close()
    '''

if lacate:
    root_path = '/Users/alexhill/Documents/UOL/Research/Companies/ViBo/Metabolomics/Data_Analysis/700MHzLactate/'
    pref = 'AHill_241025_7_exp'
    suff = '.txt'
    ray = np.arange(10, 130, 10)
    path = root_path
    print(path)

    NOES = []
    for i in ray:
        filename = pref+str(i)+suff
        print(filename)
        x, y = getData(path, filename)
        NOES.append([x,y])


    NOES = np.array(NOES)

    data_file = 'Data.hdf5'
    concs = ['750uM', '1500uM', '3000uM', '6000uM']
    concs = ['6000uM', '3000uM', '1500uM', '750uM']

    u = h5.File(data_file, 'r+')

    for i in range(0, len(concs)):
        start = int(i*3)
        conc = concs[i]
        
        root = 'Lactate/700MHz/'+conc+'/'

        XA = NOES[start + 0][0]
        YA = NOES[start + 0][1]

        XB = NOES[start + 1][0]
        YB = NOES[start + 1][1]

        XC = NOES[start + 2][0]
        YC = NOES[start + 2][1]

        u.create_dataset(root +'noesy/ppm_a', data =XA)
        u.create_dataset(root +'noesy/amp_a', data = YA)

        u.create_dataset(root + 'noesy/ppm_b', data = XB)
        u.create_dataset(root + 'noesy/amp_b', data = YB)

        u.create_dataset(root + 'noesy/ppm_c', data = XC)
        u.create_dataset(root + 'noesy/amp_c', data = YC)

        u.create_dataset(root + 'noesy/n_scans', data = 32)

    u.close()


if citrate:
    root_path = '/Users/alexhill/Documents/UOL/Research/Companies/ViBo/Metabolomics/Data_Analysis/700MHzCitrate/Dummy/'
    pref = 'Ahill_241122_7_exp'
    suff = '.txt'
    ray = np.arange(1, 13, 1)
    path = root_path
    print(path)

    NOES = []
    for i in ray:
        filename = pref+str(i)+suff
        print(filename)
        x, y = getData(path, filename)
        NOES.append([x,y])


    NOES = np.array(NOES)

    data_file = 'Data.hdf5'
    concs = ['50uM', '100uM', '200uM', '400uM']

    u = h5.File(data_file, 'r+')

    for i in range(0, len(concs)):
        start = int(i*3)
        conc = concs[i]
        
        root = 'Citrate/700MHz/'+conc+'/'

        XA = NOES[start + 0][0]
        YA = NOES[start + 0][1]

        XB = NOES[start + 1][0]
        YB = NOES[start + 1][1]

        XC = NOES[start + 2][0]
        YC = NOES[start + 2][1]

        u.create_dataset(root +'noesy/ppm_a', data =XA)
        u.create_dataset(root +'noesy/amp_a', data = YA)

        u.create_dataset(root + 'noesy/ppm_b', data = XB)
        u.create_dataset(root + 'noesy/amp_b', data = YB)

        u.create_dataset(root + 'noesy/ppm_c', data = XC)
        u.create_dataset(root + 'noesy/amp_c', data = YC)

        u.create_dataset(root + 'noesy/n_scans', data = 32)

    u.close()



'''
u = h5.File(data_file, 'r')
plt.figure()

for conc in concs:
    amp_a = u['Citrate/700MHz/'+conc+'/noesy/amp_a'][()]
    amp_b = u['Citrate/700MHz/'+conc+'/noesy/amp_b'][()]
    amp_c = u['Citrate/700MHz/'+conc+'/noesy/amp_c'][()]
    n_scans = u['Citrate/700MHz/'+conc+'/noesy/n_scans'][()]
    ppm_a = u['Citrate/700MHz/'+conc+'/noesy/ppm_a'][()]
    ppm_b = u['Citrate/700MHz/'+conc+'/noesy/ppm_b'][()]
    ppm_c = u['Citrate/700MHz/'+conc+'/noesy/ppm_c'][()]

    #plt.title(conc)
    plt.plot(ppm_a, amp_a)
    plt.plot(ppm_b, amp_b)
    plt.plot(ppm_c, amp_c)
plt.show()
plt.close()

u.close()
'''