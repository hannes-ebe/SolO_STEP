# Suche optimalen Magnetfeld-Offset für jeden Zeitpunkt einer Zeitreihe

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import cma
from STEP import STEP
import pickle

ebins = np.array([  0.98 ,   2.144,   2.336,   2.544,   2.784,   3.04 ,   3.312,
         3.6  ,   3.92 ,   4.288,   4.672,   5.088,   5.568,   6.08 ,
         6.624,   7.2  ,   7.84 ,   8.576,   9.344,  10.176,  11.136,
        12.16 ,  13.248,  14.4  ,  15.68 ,  17.152,  18.688,  20.352,
        22.272,  24.32 ,  26.496,  28.8  ,  31.36 ,  34.304,  37.376,
        40.704,  44.544,  48.64 ,  52.992,  57.6  ,  62.72 ,  68.608,
        74.752,  81.408,  89.088,  97.28 , 105.984, 115.2  , 125.44 ,
       137.216, 149.504, 162.816, 178.176, 194.56 , 211.968, 230.4  ,
       372.736])

# ### 2021-12-04 ###

# def grenz(t):
#     return -0.5*t + 20

# try:
#     dat = STEP(2021, 12, 4, rpath='data/STEP/', mag_path='data/mag/srf', mag_frame = 'srf')
# except:
#     dat = STEP(2021,12,4,mag_path='default',mag_frame='srf')
# period =(dt.datetime(2021,12,4,13,50),dt.datetime(2021,12,4,14,30))

# ### 2022-11-12 ###

# dat = STEP(2022, 11, 12, mag_path='default', mag_frame='srf')
# period = (dt.datetime(2022,11,12,2,40),dt.datetime(2022,11,12,3,25))
# grenz = None

### 2022-11-19 ###

# dat = STEP(2022, 11, 19, mag_path = 'default', mag_frame = 'srf')
# period = (dt.datetime(2022,11,19,14),dt.datetime(2022,11,19,16,15))
# grenz = None


# ### 2022-12-01 ###

dat = STEP(2022, 12, 1, mag_path = 'default', mag_frame = 'srf')
period = [dt.datetime(2022,12,1,6,30),dt.datetime(2022,12,1,9,30)]
grenz = None



def pw(flow,B,B_offset):
    '''Übergebe den particle flow-Vektor als Geschwindigkeit und den Magnetfeldvektor (am besten in SRF) und berechne die Pitchwinkel über das Skalarprodukt.
    Zusätzlich kann für die Magnetfeldkomponenten ein konstanter Offset übergeben werden.'''
    len_flow = np.sqrt(flow[0]**2 + flow[1]**2 + flow[2]**2)
    len_B = np.sqrt((B[0]+B_offset[0])**2 + (B[1]+B_offset[1])**2 + (B[2]+B_offset[2])**2)
    argument = (flow[0]*(B[0]+B_offset[0]) + flow[1]*(B[1]+B_offset[1]) + flow[2]*(B[2]+B_offset[2]))/len_flow/len_B
    result = np.arccos(argument)
    return result
        
def calc_pw(dat,B_offset):
    '''Berechne die Pitchwinkel für die Elektronen, welche auf STEP treffen in erster Näherung.
    Dafür wird der Winkel zwischen dem particle flow vector der Pixel und dem Magnetfeld herangezogen.
    Kann wieder einen Offset für das Magnetfeld übergeben.'''
    pws =  []
    for i in range(15):
        pws.append(pw(dat.flow_vector[i],[dat.B_R,dat.B_T,dat.B_N],B_offset))
    return np.array(pws)

def average_pw(dat,period,pitchangles,window_width=5):
    '''Berechnung der gemittelten Pitchwinkel'''
    # Maske, da ich nur die Magnetfelddaten innerhalb von period brauche:
    mask = (dat.time > period[0]) * (dat.time <= period[1])
    pw = [[] for i in range(15)]
    pw_time = []
        
    i = 0
    while (period[0] + dt.timedelta(minutes=(i+1)*window_width)) <= period[1]:
        pw_time.append(period[0] + dt.timedelta(minutes=(i+0.5)*window_width))
            
        for k in [i for i in range(1,16)]:
            # Mittelung der Pitchwinkel (k-1, da ich keine Zeit im array stehen habe)
            pw_data = pitchangles[k-1][mask]
            new_pw = np.sum(pw_data[i*window_width:(i+1)*window_width])/window_width
            pw[k-1].append(new_pw)
        i +=1
    return pw, pw_time



pixel_means = dat.calc_energy_means(ebins=ebins,head=-1, period=period, grenzfunktion=grenz, norm='ptmax')[0]



def func_to_min_each_point_in_time(B_offset,pixel_means,time_index,reference_pixel=3):
    '''Berechne gemittelten Abstand der korrigierten Energie zum Referenzpixel für jeden einzelnen Zeitpunkt.
    Minimiere den mittleren quadratischen Abstand (root mean square).'''
    deviation_of_pixels = []

    for pixel in range(1,16):
            pitchangles = calc_pw(dat.mag,B_offset)
            pw, pw_time = average_pw(dat.mag,period,pitchangles)

            energy2_corrected = dat.energy_correction(pixel_means[pixel],pw[reference_pixel-1],pw[pixel-1],2,2)[0]
            diff_corrected = energy2_corrected[time_index] - pixel_means[reference_pixel][time_index]
            deviation_of_pixels.append(diff_corrected)

    total_deviation = np.sqrt(np.sum(np.array(deviation_of_pixels)**2)/len(deviation_of_pixels))
    print('calculated function to minimize')
    return total_deviation



B_offset0 = [0,0,0]
B_offsets_ts = []
es_pixel = []

len_ts = len(pixel_means[0])

datei = open('B_offsets_ts.txt','x')  # x erzeugt neue Datei und gibt Fehler, falls sie bereits existiert.
datei.write("Optimales Magnetfeld fuer die einzelnen Zeitpunkte:")

for point_in_time in range(0,len_ts):
    print(f'Started calculating point in time {point_in_time}')
    x_opt, es = cma.fmin2(objective_function = func_to_min_each_point_in_time, x0 = [B_offset0], sigma0 = 0.5, args = [pixel_means,point_in_time], options = {'maxfevals':500})
    B_offsets_ts.append(x_opt)
    es_pixel.append(es)
    datei.write(f"\npot{point_in_time}: {x_opt}")
    print(f'Finished calculating point in time {point_in_time}')
    
datei.write("\nZeitreihe der optimalen Magnetfelder:")
datei.write(f"\n{B_offsets_ts}")
datei.close()

print('Idealer Magnetfeld-Offset (time series):')
for t in range(0,len_ts):
    print(f'Point in time {t}: {B_offsets_ts[t]}')
    
def step_plot_ideal_offsets_each_ts(dat, period, grenz, Offsets_ts, title=None):
    '''Übergebe Zeitreihe von idealen Magnetfeld-Offsets, um damit ideale Korrektur zu berechnen.'''

    pixel_means, pixel_var = dat.calc_energy_means(ebins=ebins,head=-1, period=period, grenzfunktion=grenz, norm='ptmax')
    
    len_ts = len(pixel_means[0])
    
    # Offsets = np.array([np.zeros(11),np.zeros(11),np.array([i for i in range(-5,6)])]).T

    pixel1 = 3
    
    fig, ax = dat.step_plot('time', 'difference of energy means [keV]', f'difference of corrected energy means to pixel {pixel1}')

    year = str(period[0].year - 2000)
    if period[0].month < 10:
         month = '0' + str(period[0].month)
    else:
        month = str(period[0].month)
    if period[0].day < 10:
        day = '0' + str(period[0].day)
    else:
        day = str(period[0].day)

    for pixel2 in range(1,16):
        diff_of_pixel = []
        for pot in range(len_ts):
            # Gehen die folgenden zwei zeilen einfacher?
            pitchangles = calc_pw(dat.mag,Offsets_ts[pot])
            pw, pw_time = average_pw(dat.mag,period,pitchangles)
            
            # Übergebe willkürliche Fehler, da ich diese eh nicht brauche.
            energy2_corrected = dat.energy_correction(pixel_means[pixel2][pot],pw[pixel1-1][pot],pw[pixel2-1][pot],2,2)[0]
            diff_corrected = energy2_corrected - pixel_means[pixel1][pot]
            diff_of_pixel.append(diff_corrected)
            
        ax[pixel2].plot(pixel_means[0],diff_of_pixel,marker='x')
        ax[pixel2].axhline(0,color='tab:red')
                
        ax[pixel2].tick_params(axis='x',labelrotation=45)
    if type(title) == str:
        plt.savefig(f'mag_variation_cmaes/' + title + '.png')
    else:
        plt.savefig(f'mag_variation_cmaes/step_plot_total_correction_differences_energy_pixel{pixel1}_{year}_{month}_{day}_multiple_offsets.png')

    # Noch einen Plot der Zeitreihe der idealen Abweichung:
    plt.figure(figsize=(10,8))
    offs = np.array(Offsets_ts).T
    labels = [r'$B_X$',r'$B_Y$',r'$B_T$']
    for i,off in enumerate(offs):
        plt.plot(pixel_means[0],off,label=labels[i])
    plt.tick_params(axis='x',labelrotation=45)
    plt.title('Ideal Offsets for Magnetic Field')
    plt.legend()
    plt.xlabel('time')
    plt.ylabel('magnetic field component in nT')

    if type(title) == str:
        plt.savefig(f'mag_variation_cmaes/' + title + '_mag.png')
    else:
        plt.savefig(f'mag_variation_cmaes/step_plot_total_correction_differences_energy_pixel{pixel1}_{year}_{month}_{day}_multiple_offsets_mag.png')
    plt.close('all')

# step_plot_ideal_offsets_each_ts(dat,period,grenz,B_offsets_ts,title='ideal_mag_offsets_ts_2021_12_04_rms')
# step_plot_ideal_offsets_each_ts(dat,period,grenz,B_offsets_ts,title='ideal_mag_offsets_ts_2022_11_12_rms')
# step_plot_ideal_offsets_each_ts(dat,period,grenz,B_offsets_ts,title='ideal_mag_offsets_ts_2022_11_19_rms')
step_plot_ideal_offsets_each_ts(dat,period,grenz,B_offsets_ts,title='ideal_mag_offsets_ts_2022_12_01_rms')