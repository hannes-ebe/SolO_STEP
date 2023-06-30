#################################################################################################
# Neuer Plausibilitätscheck für Pitchwinkeleinflüsse. Die Idee ist, dass Pixel, die             #
# annähernd die gleichen Pitchwinkel sehen, auch die gleichen Energiemittelwerte haben sollten. #
# Versuche einen ersten Vergleich durchhzuführen.                                               #
#################################################################################################

from STEP import STEP
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt


ebins = np.array([  0.98 ,   2.144,   2.336,   2.544,   2.784,   3.04 ,   3.312,
         3.6  ,   3.92 ,   4.288,   4.672,   5.088,   5.568,   6.08 ,
         6.624,   7.2  ,   7.84 ,   8.576,   9.344,  10.176,  11.136,
        12.16 ,  13.248,  14.4  ,  15.68 ,  17.152,  18.688,  20.352,
        22.272,  24.32 ,  26.496,  28.8  ,  31.36 ,  34.304,  37.376,
        40.704,  44.544,  48.64 ,  52.992,  57.6  ,  62.72 ,  68.608,
        74.752,  81.408,  89.088,  97.28 , 105.984, 115.2  , 125.44 ,
       137.216, 149.504, 162.816, 178.176, 194.56 , 211.968, 230.4  ,
       372.736])


dat = STEP(2021, 12, 4, magnet_default_path=True)
# STEP(2021, 12, 4, rpath='data/STEP/')

period =(dt.datetime(2021,12,4,13,50),dt.datetime(2021,12,4,14,30))

def grenz(t):
    return -0.5*t + 20

pixel_means5 = dat.calc_energy_means(ebins=ebins,head=-1, period=period, grenzfunktion=grenz)
# pixel_means1 = dat.calc_energy_means(ebins=ebins,head=-1, period=period, grenzfunktion=grenz, window_width=1)

pw5 = dat.calc_pw(period, window_width=5)
# pw1 = dat.calc_pw(period, window_width=1)


### Plots für Pixel mit ähnlichen Pitchwinkeln: ###

dat.pixel_comparison(pixel_means5, pw5, 'plausibility_pw/comparison_5_minutes_pixel_1_5.png', 1, 5)
dat.pixel_comparison(pixel_means5, pw5, 'plausibility_pw/comparison_5_minutes_pixel_2_4.png', 2, 4)
dat.pixel_comparison(pixel_means5, pw5, 'plausibility_pw/comparison_5_minutes_pixel_6_10.png', 6, 10)
dat.pixel_comparison(pixel_means5, pw5, 'plausibility_pw/comparison_5_minutes_pixel_7_9.png', 7, 9)
dat.pixel_comparison(pixel_means5, pw5, 'plausibility_pw/comparison_5_minutes_pixel_11_15.png', 11, 15)
dat.pixel_comparison(pixel_means5, pw5, 'plausibility_pw/comparison_5_minutes_pixel_12_14.png', 12, 14)
    
# dat.pixel_comparison(pixel_means1, pw1, 'plausibility_pw/comparison_1_minutes_pixel_1_5.png', 1, 5)
# dat.pixel_comparison(pixel_means1, pw1, 'plausibility_pw/comparison_1_minutes_pixel_2_4.png', 2, 4)
# dat.pixel_comparison(pixel_means1, pw1, 'plausibility_pw/comparison_1_minutes_pixel_6_10.png', 6, 10)
# dat.pixel_comparison(pixel_means1, pw1, 'plausibility_pw/comparison_1_minutes_pixel_7_9.png', 7, 9)
# dat.pixel_comparison(pixel_means1, pw1, 'plausibility_pw/comparison_1_minutes_pixel_11_15.png', 11, 15)
# dat.pixel_comparison(pixel_means1, pw1, 'plausibility_pw/comparison_1_minutes_pixel_12_14.png', 12, 14)


plt.close("all")