######################################################
### Creating some extra plots for the SolO meeting ###
######################################################

from STEP import STEP
from mag import MAGdata
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from cdflib import CDF


ebins = np.array([  0.98 ,   2.144,   2.336,   2.544,   2.784,   3.04 ,   3.312,
         3.6  ,   3.92 ,   4.288,   4.672,   5.088,   5.568,   6.08 ,
         6.624,   7.2  ,   7.84 ,   8.576,   9.344,  10.176,  11.136,
        12.16 ,  13.248,  14.4  ,  15.68 ,  17.152,  18.688,  20.352,
        22.272,  24.32 ,  26.496,  28.8  ,  31.36 ,  34.304,  37.376,
        40.704,  44.544,  48.64 ,  52.992,  57.6  ,  62.72 ,  68.608,
        74.752,  81.408,  89.088,  97.28 , 105.984, 115.2  , 125.44 ,
       137.216, 149.504, 162.816, 178.176, 194.56 , 211.968, 230.4  ,
       372.736])


### 2021-12-4 ###

# dat = STEP(2021, 12, 4, magnet_default_path=True)
dat = STEP(2021, 12, 4, rpath='data/STEP/', rpath_mag='data/mag/')
period =[dt.datetime(2021,12,4,13,50),dt.datetime(2021,12,4,14,30)]
def grenz(t):
    return -0.5*t + 20


pixel_means5 = dat.calc_energy_means(ebins=ebins,head=-1, period=period, grenzfunktion=grenz)



### example plot for energy means pixel 5 ###

dat.plot_energy_means(pixel_means5, 'SolO_meeting/energy_means_2021_12_04_pixel_4.png', pixel_list=[4])

### time series of pitch angles ###

dat.mag.pw_ts('SolO_meeting/pw_2021_12_04.png',[4],period=period,window_width=5)

### time series of magnetic field ###

dat.mag.mag_ts('SolO_meeting/magnetic_field_2021_12_04.png',period=period)


plt.close('all')