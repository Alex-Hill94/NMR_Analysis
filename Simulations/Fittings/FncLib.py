import numpy as np
import random
import pandas as pd
import os
import h5py as h5
from LineshapeCreator import createLineshape
import matplotlib.pyplot as plt

#points = 65536
#limits = (12, -1)
#width = 2.0
#scale = 3*1e8

def GenData(x_array, scale, x_shift = 0, substance = 'Glucose', width = 2.):
    #width = 2.0
    points = len(x_array)
    limits = (np.nanmax(x_array), np.nanmin(x_array))
    frequency = 80.0
    signals = pd.read_csv("/Users/alexhill/Documents/GitHub/NMR_Analysis/Simulations/Fittings/%s_signals.csv" % substance)
    signalList = [(row['chemical_shift'] + x_shift, row['height'], width / frequency) for index, row in signals.iterrows()]
    x, y = createLineshape(signalList, points=points, limits=limits)
    y = y * scale
    return y

def Glucose(x_array = np.linspace(12, -1, 65536), scale = 1., x_shift = 0.0, width = 2.0):
    y = GenData(x_array, scale, x_shift = x_shift, width = width, substance = 'Glucose')
    return y

def Lactate(x_array = np.linspace(12, -1, 65536), scale = 1., x_shift = 0.0, width = 2.0):
    y = GenData(x_array, scale, x_shift = x_shift, width = width, substance = 'Lactate')
    return y
    
def Citrate(x_array = np.linspace(12, -1, 65536), scale = 1., x_shift = 0.0, width = 2.0):
    y = GenData(x_array, scale, x_shift = x_shift, width = width, substance = 'Citrate')
    return y
    
if __name__ == '__main__':

    x = np.linspace(5, 0, 65536)
    y = Glucose(x_array = x)

    x1 = np.linspace(7, -1, 655360)
    y1 = Glucose(x_array = x1, scale = 2., x_shift = 1.0)

    plt.figure()
    plt.plot(x1, y1)
    plt.plot(x, y)
    plt.show()