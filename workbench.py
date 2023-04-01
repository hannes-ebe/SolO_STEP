'''Workbench zum Herumprobieren'''

import matplotlib as mpl
import pylab as plt
import numpy as np
import datetime as dt
import matplotlib.dates as mdates
import load_nom_II as ld
import plot_nom_II as pt
from scipy.optimize import curve_fit
from analyze_nom_II import STEP

ebins = np.array([  0.98 ,   2.144,   2.336,   2.544,   2.784,   3.04 ,   3.312,
         3.6  ,   3.92 ,   4.288,   4.672,   5.088,   5.568,   6.08 ,
         6.624,   7.2  ,   7.84 ,   8.576,   9.344,  10.176,  11.136,
        12.16 ,  13.248,  14.4  ,  15.68 ,  17.152,  18.688,  20.352,
        22.272,  24.32 ,  26.496,  28.8  ,  31.36 ,  34.304,  37.376,
        40.704,  44.544,  48.64 ,  52.992,  57.6  ,  62.72 ,  68.608,
        74.752,  81.408,  89.088,  97.28 , 105.984, 115.2  , 125.44 ,
       137.216, 149.504, 162.816, 178.176, 194.56 , 211.968, 230.4  ,
       372.736])

# dat = STEP('data/',2021,12,4)

# pldat, pltime, vmax = dat.data_prep(ebins,res = '1min', head = 0, period = [dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)], norm = 'tmax', overflow = True, esquare = False)

# pdat = pldat[4]
# ptime = np.append(pltime[4],pltime[4][-1]+dt.timedelta(seconds=60))

# print(pdat)
# print(len(pdat))
# print(pdat.T)
# print(len(pdat.T))

# print(ptime)
# print(len(ptime))
# print(ptime.T)
# print(len(ptime.T))

# print(len(ebins))

a = np.array([0.0,0.0,2.0,3.4,0.0])
print(np.mean(a))