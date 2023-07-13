#################################################################################################
# Neuer Plausibilitätscheck für Pitchwinkeleinflüsse. Die Idee ist, dass Pixel, die             #
# annähernd die gleichen Pitchwinkel sehen, auch die gleichen Energiemittelwerte haben sollten. #
# Versuche einen ersten Vergleich durchhzuführen.                                               #
#################################################################################################

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
# # STEP(2021, 12, 4, rpath='data/STEP/')
# period =(dt.datetime(2021,12,4,13,50),dt.datetime(2021,12,4,14,30))
# def grenz(t):
#     return -0.5*t + 20


### 2022-2-7 ###

# dat = STEP(2022,2,7, magnet_default_path=True)
# period=[dt.datetime(2022,2,7,20,15),dt.datetime(2022,2,7,21,15)]
# grenz = None


### 2022-11-12 ###

# dat = STEP(2022,11,12, magnet_default_path=True)
# period = [dt.datetime(2022,11,12,2,30),dt.datetime(2022,11,12,4)]
# grenz = None


### 2022-11-19 ###

dat = STEP(2022,11,19, magnet_default_path=True)
period = [dt.datetime(2022,11,19,13,30),dt.datetime(2022,11,19,16,30)]
grenz = None





# pixel_means5 = dat.calc_energy_means(ebins=ebins,head=-1, period=period, grenzfunktion=grenz)
# pixel_means1 = dat.calc_energy_means(ebins=ebins,head=-1, period=period, grenzfunktion=grenz, window_width=1)

### Achtung hier mal andere Normierung!!! ###
pixel_means5 = dat.calc_energy_means(ebins=ebins,head=-1, period=period, grenzfunktion=grenz, norm='ptmax')

pw5, pw5_time = dat.calc_pw(period, window_width=5)
# pw1 pw1_time = dat.calc_pw(period, window_width=1)





### Plot der Zeitreihe als STEP-Plot ###

# dat.plot_ts(period=period, head=-1, save='plausibility_pw/2021_12_04/', norm='tmax',grenzfunktion=grenz)
# dat.plot_ts(period=period, head=-1, save='plausibility_pw/2022_02_07/', norm='ptmax',grenzfunktion=grenz)
# dat.plot_ts(period=period, head=-1, save='plausibility_pw/2022_11_12/', norm='tmax',grenzfunktion=grenz)
dat.plot_ts(period=period, head=-1, save='plausibility_pw/2022_11_19/', norm='ptmax',grenzfunktion=grenz)



### Plot der Verteilung von Energiemittelwert und Pitchwinkel ###

# dat.distribution_ring('time_series_energy_means_pw_2021_12_04','mean of energy (time series)',head=-1,window_width=1,norm='tmax',save='plausibility_pw/2021_12_04/',period=period,grenzfunktion=grenz,below=True,close=True)
# dat.distribution_ring('time_series_energy_means_pw_2022_02_07','mean of energy (time series)',head=-1,window_width=1,norm='ptmax',save='plausibility_pw/2022_02_07/',period=period,grenzfunktion=grenz,below=True,close=True)
# dat.distribution_ring('time_series_energy_means_pw_2022_11_12','mean of energy (time series)',head=-1,window_width=1,norm='tmax',save='plausibility_pw/2022_11_12/',period=period,grenzfunktion=grenz,below=True,close=True)
dat.distribution_ring('time_series_energy_means_pw_2022_11_19','mean of energy (time series)',head=-1,window_width=1,norm='ptmax',save='plausibility_pw/2022_11_19/',period=period,grenzfunktion=grenz,below=True,close=True)



### Plots für Pixel mit ähnlichen Pitchwinkeln: ###

def compare_opposing_pixels(dir):
    dat.pixel_comparison(pixel_means5, pw5, pw5_time, 'plausibility_pw/' + dir + '/comparison_5_minutes_pixel_1_5.png', 1, 5)
    dat.pixel_comparison(pixel_means5, pw5, pw5_time, 'plausibility_pw/' + dir + '/comparison_5_minutes_pixel_2_4.png', 2, 4)
    dat.pixel_comparison(pixel_means5, pw5, pw5_time, 'plausibility_pw/' + dir + '/comparison_5_minutes_pixel_6_10.png', 6, 10)
    dat.pixel_comparison(pixel_means5, pw5, pw5_time, 'plausibility_pw/' + dir + '/comparison_5_minutes_pixel_7_9.png', 7, 9)
    dat.pixel_comparison(pixel_means5, pw5, pw5_time, 'plausibility_pw/' + dir + '/comparison_5_minutes_pixel_11_15.png', 11, 15)
    dat.pixel_comparison(pixel_means5, pw5, pw5_time, 'plausibility_pw/' + dir + '/comparison_5_minutes_pixel_12_14.png', 12, 14)

    dat.pixel_comparison_corrected(pixel_means5, pw5, 'plausibility_pw/' + dir + '/comparison_corrected_5_minutes_pixel_1_5.png', 1, 5)
    dat.pixel_comparison_corrected(pixel_means5, pw5, 'plausibility_pw/' + dir + '/comparison_corrected_5_minutes_pixel_2_4.png', 2, 4)
    dat.pixel_comparison_corrected(pixel_means5, pw5, 'plausibility_pw/' + dir + '/comparison_corrected_5_minutes_pixel_6_10.png', 6, 10)
    dat.pixel_comparison_corrected(pixel_means5, pw5, 'plausibility_pw/' + dir + '/comparison_corrected_5_minutes_pixel_7_9.png', 7, 9)
    dat.pixel_comparison_corrected(pixel_means5, pw5, 'plausibility_pw/' + dir + '/comparison_corrected_5_minutes_pixel_11_15.png', 11, 15)
    dat.pixel_comparison_corrected(pixel_means5, pw5, 'plausibility_pw/' + dir + '/comparison_corrected_5_minutes_pixel_12_14.png', 12, 14)
    

# compare_opposing_pixels('2021_12_04')
# compare_opposing_pixels('2022_02_07')
# compare_opposing_pixels('2022_11_12')
compare_opposing_pixels('2022_11_19')



# dat.pixel_comparison(pixel_means1, pw1, pw1_time, 'plausibility_pw/comparison_1_minutes_pixel_1_5.png', 1, 5)
# dat.pixel_comparison(pixel_means1, pw1, pw1_time, 'plausibility_pw/comparison_1_minutes_pixel_2_4.png', 2, 4)
# dat.pixel_comparison(pixel_means1, pw1, pw1_time, 'plausibility_pw/comparison_1_minutes_pixel_6_10.png', 6, 10)
# dat.pixel_comparison(pixel_means1, pw1, pw1_time, 'plausibility_pw/comparison_1_minutes_pixel_7_9.png', 7, 9)
# dat.pixel_comparison(pixel_means1, pw1, pw1_time, 'plausibility_pw/comparison_1_minutes_pixel_11_15.png', 11, 15)
# dat.pixel_comparison(pixel_means1, pw1, pw1_time, 'plausibility_pw/comparison_1_minutes_pixel_12_14.png', 12, 14)


plt.close("all")


### Calculating errors in pw necessary to do proper correction. ###

# # Assuming that both pw have same error
# # 4th data point of pixel 12 and 14

# energy = 19.5
# energy_corrected_err = 4
# pw_err = np.radians(10)
# pw1 = np.radians(46.15)
# pw2 = np.radians(46.05)
# factor = 2.0*energy*np.cos(pw2)/np.cos(pw1)/np.cos(pw1) * (np.sin(pw1)*np.cos(pw2)/np.cos(pw1) + np.sin(pw2))
# pw_err = energy_corrected_err/factor
# print(pw_err)
# print(np.degrees(pw_err))


# # Errors of pitchangles

# dat = MAGdata(period=period)

# print(dat.pitchangles_err)

# dat.pw_ts(err=True)






### Checking out the cdf-files of SolO mag data ###

# dat = CDF('/data/projects/solo/mag/l2_soar/rtn_1minute/2021/solo_L2_mag-rtn-normal-1-minute_20211204_V01.cdf')
# mag = dat.varget('B_RTN')
# lbl = dat.varget('LBL1_B_RTN')
# rep = dat.varget('REP1_B_RTN')
# print(mag)
# print(lbl)
# print(rep)
# print(dat.cdf_info())