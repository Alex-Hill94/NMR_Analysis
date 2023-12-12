import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pylab import rc, Line2D

custom_lines = [Line2D([0], [0], color= 'white',ls = '-', lw=2, alpha = 0),
				Line2D([0], [0], color= 'white',ls = '-', lw=2, alpha = 0)]


def grab_data(path = "/Users/alexhill/Desktop/Metabolomics/Data_Analysis", file_name = "60MHz_standards_230607.txt"):
    # Define the file path
    file_path = path+"/"+file_name

    # Extract the limits and size
    left_limit = None
    right_limit = None
    size = None

    with open(file_path, 'r') as file:
        for line in file:
            if 'LEFT' in line:
                left_index = line.find('LEFT')
                right_index = line.find('RIGHT')
                if left_index != -1:
                    left_start = line.find('=', left_index) + 1
                    left_end = line.find('ppm', left_start)
                    left_limit = float(line[left_start:left_end].strip())
                if right_index != -1:
                    right_start = line.find('=', right_index) + 1
                    right_end = line.find('ppm', right_start)
                    right_limit = float(line[right_start:right_end].strip())
            if 'SIZE' in line:
                size = int(line.split('=')[1].split('(')[0].strip())
                break  # Stop reading the file after finding thesize

    # Generate the x-axis values
    x_values = np.linspace(left_limit, right_limit, size)

    # Read the file and extract the number list
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                data.append(float(line))
    return x_values, data

def plot_lines(ax, file_name = "60MHz_standards_230607.txt", col = 'k', lw = 1.5, ls = '-', alpha = 0.7, label = None, zoom = 1.):
    x_values, data = grab_data(file_name = file_name)
    ax.plot(x_values, np.array(data)*zoom, color = col, lw = lw, ls = ls, alpha = alpha, label = label)

x,y = grab_data(file_name = "Nov29-2023_D22_LC+NG.txt")
