from __future__ import division
from numpy import matrix
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.colorbar
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker   
import matplotlib.gridspec as gridspec 
from matplotlib.patches import Ellipse
from matplotlib.colors import LogNorm
import numpy as np

paul_blue = ( 0./255., 119./255., 187./255.)
paul_teal = (0./255., 153./255., 136./255.)
paul_red  = (204./255., 51./255., 17./255.)
paul_magenta  = (238./255., 51./255., 119./255.)
paul_cyan  = (51./255., 187./255., 238./255.)

paul_rainbow = ((25./255., 101./255., 176./255.), 
				(123./255., 175./255., 222./255.),
				(78./255., 178./255., 101./255.),
				(202./255., 224./255., 171./255.),
				(247./255., 240./255., 86./255.),
				(220./255., 5./255., 12./255.))

paul_option  = ((111./255., 76./255., 155./255.), 
				(85./255., 104./255., 184./255.),
				(77./255., 138./255., 198./255.),
				(84./255., 158./255., 179./255.),
				(96./255., 171./255., 158./255.),
				(119./255., 183./255., 125./255.))

paul_p_blue =  (187./255., 204./255., 238./255.)
paul_p_cyan =  (204./255., 238./255., 255./255.)
paul_p_green = (204./255., 221./255., 170./255.)
paul_p_red =   (255./255., 204./255., 204./255.)



paul_dark_blue = (54./255. , 75./255., 154./255.)
paul_med_blue = (74./255. , 123./255., 183./255.)
paul_light_blue = (146./255. , 197./255., 222./255.)


paul_dark_green = (27./255. , 120./255., 55./255.)
paul_med_green = (90./255. , 174./255., 97./255.)
paul_light_green = (172./255. , 211./255., 158./255.)


paul_dark_red = (165./255. , 0./255., 38./255.)
paul_med_red = (221./255. , 61./255., 45./255.)
paul_light_red = (244./255. , 165./255., 130./255.)

paul_muted =  ((136./255., 204./255., 238./255.), 
				(68./255., 170./255., 153./255.),
				(221./255., 204./255., 119./255.),
				(204./255., 102./255., 119./255.),
				(136./255., 34./255., 85./255.))
