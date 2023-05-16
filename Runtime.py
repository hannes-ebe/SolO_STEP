import matplotlib as mpl
import pylab as plt
import numpy as np
import datetime as dt
import matplotlib.dates as mdates
import load_nom_II as ld
import plot_nom_II as pt
import math
import mag
import scipy.constants as const

from PIL import Image
import glob
from STEP import STEP

from scipy.optimize import curve_fit

ebins = np.array([  0.98 ,   2.144,   2.336,   2.544,   2.784,   3.04 ,   3.312,
         3.6  ,   3.92 ,   4.288,   4.672,   5.088,   5.568,   6.08 ,
         6.624,   7.2  ,   7.84 ,   8.576,   9.344,  10.176,  11.136,
        12.16 ,  13.248,  14.4  ,  15.68 ,  17.152,  18.688,  20.352,
        22.272,  24.32 ,  26.496,  28.8  ,  31.36 ,  34.304,  37.376,
        40.704,  44.544,  48.64 ,  52.992,  57.6  ,  62.72 ,  68.608,
        74.752,  81.408,  89.088,  97.28 , 105.984, 115.2  , 125.44 ,
       137.216, 149.504, 162.816, 178.176, 194.56 , 211.968, 230.4  ,
       372.736])
cmap = mpl.cm.jet
cmap.set_under('w')
hmap = mpl.cm.seismic



class STEP_Runtime(STEP):
    def __init__(self,year,month,day,rpath = '/data/projects/solo/step_v0008/',rpath_mag = '/data/projects/solo/mag/l2_soar/rtn_1minute',magnet=False,lastofmonth=False):
        super().__init__(year,month,day,rpath=rpath,rpath_mag=rpath_mag,magnet=magnet,lastofmonth=lastofmonth)
        
    def rel_speed(self,E,m):
        '''Berechnet die relativistische Geschwindigkeit zur übergebenen kinetischen Energie E und Masse m.'''
        ayuda = E/m/const.c/const.c + 1.0
        return const.c*np.sqrt(1.0-1.0/ayuda/ayuda)
    



# dat = STEP_Runtime(2021, 12, 4)
dat = STEP_Runtime(2021, 12, 4, rpath='data/')

box_list = [[[15,35],[30,38]],[[20,45],[25,30]],[[26,80],[20,25]],[[30,120],[15,20]],[[40,155],[10,15]],[[50,170],[3,10]]]
period = (dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30))

pixel_means = dat.calc_energy_means(ebins=ebins,head=-1, period=period, box_list=box_list)
pixel_means2 = dat.calc_energy_means(ebins=ebins,head=-1, period=period, box_list=box_list, window_width=1)


# Geschwindigkeit:
# Korrektur, da Energie-Mittelwerte in keV angegeben sind.
v = dat.rel_speed(pixel_means[3]*1.6022e-16, const.m_e)

# Filter die nans raus
divided_v = 1/v
mask = np.isfinite(divided_v)
divided_v = divided_v[mask]

v2 = dat.rel_speed(pixel_means2[3]*1.6022e-16, const.m_e)
divided_v2 = 1/v2
mask2 = np.isfinite(divided_v2)
divided_v2 = divided_v2[mask2]

def lin(x,m,b):
    return m*x+b

# Unix time stamp (Sekunden seit 1. Januar 1970) für curve_fit()
ts = np.array([t.timestamp() for t in pixel_means[0]])[mask]
popt, pcov = curve_fit(lin,ts,divided_v)

ts2 = np.array([t.timestamp() for t in pixel_means2[0]])[mask2]
popt2, pcov2 = curve_fit(lin,ts2,divided_v2)

plt.plot(ts,lin(ts,popt[0],popt[1]),label='fit one minute data')
plt.plot(ts2,lin(ts2,popt2[0],popt2[1]),label='fit five minute data')

    
# Counts über fünf Minuten oder eine Minute summiert.
plt.scatter(ts,divided_v,marker='x',label='five minutes')
plt.scatter(ts2,divided_v2,marker='x',label='one minute')
plt.title('pixel 3; 2021-21-04')
plt.xlabel('unix time stamp [s]')
plt.ylabel(r'$v^{-1}~[\frac{\mathrm{m}}{\mathrm{s}}]$')
plt.legend()
plt.savefig('runtime/test_runtime_pixel3_2021_12_04.png')
plt.show()