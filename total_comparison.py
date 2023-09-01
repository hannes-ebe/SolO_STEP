#####################################################################
# Comparing each pixel with each other. Correcting the energies for #
# every combination (in both directions).                           #
##################################################################### 

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

# dat = STEP(2021, 12, 4, mag_path = 'default', mag_frame = 'srf')
# # dat = STEP(2021, 12, 4, rpath='data/STEP/', mag_path='data/mag/srf', mag_frame = 'srf')
# period =(dt.datetime(2021,12,4,13,50),dt.datetime(2021,12,4,14,30))
# def grenz(t):
#     return -0.5*t + 20


### 2022-11-12 ###

dat = STEP(2022, 11, 12, mag_path='default', mag_frame='srf')
period = (dt.datetime(2022,11,12,2,40),dt.datetime(2022,11,12,3,25))
grenz = None


### 2022-11-19 ###

# dat = STEP(2022, 11, 19, mag_path = 'default', mag_frame = 'srf')
# period = (dt.datetime(2022,11,19,14),dt.datetime(2022,11,19,16,15))
# grenz = None



def func_total_comparison(dat=dat, period=period, grenz=grenz):
    pixel_means5, pixel_var5 = dat.calc_energy_means(ebins=ebins,head=-1, period=period, grenzfunktion=grenz, norm='ptmax')
    pw5, pw5_time = dat.calc_pw(period, window_width=5)
    year = str(period[0].year - 2000)
    if period[0].month < 10:
        month = '0' + str(period[0].month)
    else:
        month = str(period[0].month)
    if period[0].day < 10:
        day = '0' + str(period[0].day)
    else:
        day = str(period[0].day)
    
    
    dat.plot_ts(period=period, head=-1, save='total_comparison/', norm='ptmax',grenzfunktion=grenz)
    dat.distribution_ring(f'time_series_energy_means_pw_{year}_{month}_{day}','mean of energy (time series)',head=-1,window_width=1,norm='ptmax',save='total_comparison/',period=period,grenzfunktion=grenz,below=True,close=True)
    print(f'{year}_{month}_{day}')
    
    # for i in range(1,16):
    #     for j in range(1,16):
    #         # Bei gleichen Pixeln dürfte sich nichts ändern...
    #         dat.pixel_comparison_corrected(pixel_means5, pixel_var5, pw5, f'total_comparison/corrected_energies_5_minutes_pixel_{i}_{j}.png', i, j)
    #         plt.close('all')
    
    
    ### Untersuche Pixel 4 etwas genauer ###
    
    for j in range(1,16):
        dat.pixel_comparison_total(pixel_means5, pixel_var5, pw5, pw5_time, f'total_comparison/pixel4_{year}_{month}_{day}/corrected_energies_total_info_pixel_i=4_{j}.png', 4, j)
        plt.close('all')
    
    for i in range(1,16):
        dat.pixel_comparison_total(pixel_means5, pixel_var5, pw5, pw5_time, f'total_comparison/pixel4_{year}_{month}_{day}/corrected_energies_total_info_pixel_{i}_j=4.png', i, 4)
        plt.close('all')
    
    
    ### Teste mal Korrektur auf Pixel 11 ###
    
    for j in range(1,16):
        dat.pixel_comparison_total(pixel_means5, pixel_var5, pw5, pw5_time, f'total_comparison/pixel11_{year}_{month}_{day}/corrected_energies_total_info_pixel_i=11_{j}.png', 11, j)
        plt.close('all')
    
    for i in range(1,16):
        dat.pixel_comparison_total(pixel_means5, pixel_var5, pw5, pw5_time, f'total_comparison/pixel11_{year}_{month}_{day}/corrected_energies_total_info_pixel_{i}_j=11.png', i, 11)
        plt.close('all')
        
        
def step_plot_correction(dat=dat, period=period, grenz=grenz):
    pixel_means, pixel_var = dat.calc_energy_means(ebins=ebins,head=-1, period=period, grenzfunktion=grenz, norm='ptmax')
    pw, pw_time = dat.calc_pw(period, window_width=5)
    year = str(period[0].year - 2000)
    if period[0].month < 10:
        month = '0' + str(period[0].month)
    else:
        month = str(period[0].month)
    if period[0].day < 10:
        day = '0' + str(period[0].day)
    else:
        day = str(period[0].day)
    
    fig, ax = dat.step_plot('time', 'mean of energy [keV]', 'energy means with pitch angle correction')
    
    pixel1 = 4
    for pixel2 in range(1,16):
        ax[pixel2].errorbar(pixel_means[0],pixel_means[pixel1],yerr=np.sqrt(pixel_var[pixel1]),marker='x',label=f'mean pixel {pixel1}')
            
        # Übergebe willkürliche Fehler, da ich diese eh nicht brauche.
        energy2_corrected = dat.energy_correction(pixel_means[pixel2],pw[pixel1-1],pw[pixel2-1],2,2)[0]
            
        ax[pixel2].errorbar(pixel_means[0],energy2_corrected,yerr=np.sqrt(pixel_var[pixel2]),marker='x',label=f'mean pixel {pixel2}')
        ax[pixel2].tick_params(axis='x',labelrotation=45)
        ax[pixel2].legend()
    plt.savefig(f'total_comparison/step_plot_total_correction_energy_means_pixel{pixel1}_{year}_{month}_{day}.png')
    plt.close('all')
    

    fig, ax = dat.step_plot('time', 'difference of energy means [keV]', 'difference of energy means with correction')
        
    for pixel2 in range(1,16):
        # Übergebe willkürliche Fehler, da ich diese eh nicht brauche.
        energy2_corrected = dat.energy_correction(pixel_means[pixel2],pw[pixel1-1],pw[pixel2-1],2,2)[0]
        diff_corrected = energy2_corrected - pixel_means[pixel1]
        
        ax[pixel2].plot(pixel_means[0],diff_corrected,marker='x')
        ax[pixel2].axhline(0,color='tab:red')
            
        ax[pixel2].tick_params(axis='x',labelrotation=45)
    plt.savefig(f'total_comparison/step_plot_total_correction_differences_energy_pixel{pixel1}_{year}_{month}_{day}.png')
    plt.close('all')


func_total_comparison(dat,period,grenz)
step_plot_correction(dat,period,grenz)