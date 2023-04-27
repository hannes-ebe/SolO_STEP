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
import mag

ebins = np.array([  0.98 ,   2.144,   2.336,   2.544,   2.784,   3.04 ,   3.312,
         3.6  ,   3.92 ,   4.288,   4.672,   5.088,   5.568,   6.08 ,
         6.624,   7.2  ,   7.84 ,   8.576,   9.344,  10.176,  11.136,
        12.16 ,  13.248,  14.4  ,  15.68 ,  17.152,  18.688,  20.352,
        22.272,  24.32 ,  26.496,  28.8  ,  31.36 ,  34.304,  37.376,
        40.704,  44.544,  48.64 ,  52.992,  57.6  ,  62.72 ,  68.608,
        74.752,  81.408,  89.088,  97.28 , 105.984, 115.2  , 125.44 ,
       137.216, 149.504, 162.816, 178.176, 194.56 , 211.968, 230.4  ,
       372.736])

# Field of View of STEP (Phi, Theta):
fov = np.array([[-25,24],[-25,12],[-25,0],[-25,-12],[-25,-24],[-35,24],[-35,12],[-35,0],[-35,-12],[-35,-24],[-45,24],[-45,12],[-45,0],[-45,-12],[-45,-24]])



# Versuch Magnetfelddaten aus cdf zu laden

# dat = mag.MAGdata(period=[dt.datetime(2021,12,4),dt.datetime(2021,12,5)])

# dat.pw_ts()




### Plausibilitätscheck Pitchwinkel-Berechnung ###

def pw1(v_theta,v_phi,B_theta,B_phi):
    return (np.arccos(np.sin(v_theta)*np.cos(v_phi)*np.sin(B_theta)*np.cos(B_phi) + np.sin(v_theta)*np.sin(v_phi)*np.sin(B_theta)*np.sin(B_phi) + np.cos(v_theta)*np.cos(B_theta)))/np.pi*180

def pw2(v_theta,v_phi,B_theta,B_phi):
    return (np.arccos(np.sin(v_phi)*np.cos(v_theta)*np.sin(B_phi)*np.cos(B_theta) + np.sin(v_phi)*np.sin(v_theta)*np.sin(B_phi)*np.sin(B_theta) + np.cos(v_phi)*np.cos(B_phi)))/np.pi*180

B_phi = -35
B_theta = 0

plausibilitätscheck_pw = []

for i in range(len(fov)):
    plausibilitätscheck_pw.append(pw2(fov[i][1],fov[i][0],B_theta,B_phi))
    
print(plausibilitätscheck_pw)

fig, ax = plt.subplots(1,1,figsize=(8,7))

x_corners = [0,1,2,3,4,5]
y_corners = [0,1,2,3]
vmin = min(plausibilitätscheck_pw)
vmax = max(plausibilitätscheck_pw)
     
pw_list = np.array([plausibilitätscheck_pw[10:15],plausibilitätscheck_pw[5:10],plausibilitätscheck_pw[0:5]])
            
tmp = ax.pcolormesh(x_corners,y_corners,pw_list,vmin=vmin,vmax=vmax)
plt.colorbar(tmp,label='pitch angle [°]')
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
ax.set_title('plausibility check of pitch angle calculation')

plt.show()





# dat = STEP('data/',2021,12,4)
# box_list = [[[20,40],[25,35]],[[40,60],[15,25]]]

# pldat, pltime, vmax = dat.data_prep(ebins,res = '1min', head = 0, period = [dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)], norm = 'tmax', overflow = True, esquare = False)
# box_data = dat.data_boxes(pldat,box_list)

# pdat = pldat[4]
# ptime = np.append(pltime[4],pltime[4][-1]+dt.timedelta(seconds=60))

# print(len(pldat[4]))
# print(len(box_data[4]))

# for i in range(19,62):
#        print(box_data[4][i])

# print(pldat)
# print(len(pldat))

# print(pdat)
# print(len(pdat))
# print(pdat.T)
# print(len(pdat.T))

# print(ptime)
# print(len(ptime))
# print(ptime.T)
# print(len(ptime.T))

# print(len(ebins))

# a = np.zeros((181,56),dtype='float')
# print(len(a.T))