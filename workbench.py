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
# import mag

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





### Plausibilitätscheck Pitchwinkel-Berechnung ###
'''Der Winkel \phi gibt die Drehung in der RT-Ebene und \theta in der RN-Ebene. Für die Umrechnung müsste mit v als Länge des Vektors gelten:
    R = v*cos(phi)*cos(theta)
    T = v*sin(phi)*cos(theta)
    N = v*sin(theta)'''

def pw(v_theta,v_phi,B_theta,B_phi):
    '''Wichtig ist die Normierung zu beachten. Wenn ich die RTN-Koordinaten über die Winkel beschreibe, kann ich die Vorfaktoren 1 setzen. Numpy nimmt die
    Winkel in Radians. Muss den Spezialfall von Werten größer 1 im arccos abfangen.'''
    v_theta = v_theta/360*2*np.pi
    v_phi = v_phi/360*2*np.pi
    B_theta = B_theta/360*2*np.pi
    B_phi = B_phi/360*2*np.pi
    argument = np.cos(v_theta)*np.cos(v_phi)*np.cos(B_theta)*np.cos(B_phi) + np.cos(v_theta)*np.sin(v_phi)*np.cos(B_theta)*np.sin(B_phi) + np.sin(v_theta)*np.sin(B_theta)
    if argument > 1.0 and v_theta == B_theta and v_phi == B_phi:
        argument = 1.0
    pitchangle = np.arccos(argument)/np.pi*180
    return pitchangle

B_phi = -35
B_theta = 0

def plausibility_check(B_phi,B_theta):

    plausibilitätscheck_pw = []

    for i in range(len(fov)):
        plausibilitätscheck_pw.append(pw(fov[i][1],fov[i][0],B_theta,B_phi))
        
    print(plausibilitätscheck_pw)

    fig, ax = plt.subplots(1,1,figsize=(10,6))

    x_corners = [0,1,2,3,4,5]
    y_corners = [0,1,2,3]
    vmin = min(plausibilitätscheck_pw)
    vmax = max(plausibilitätscheck_pw)
        
    pw_list = np.array([plausibilitätscheck_pw[10:15],plausibilitätscheck_pw[5:10],plausibilitätscheck_pw[0:5]])
                
    tmp = ax.pcolormesh(x_corners,y_corners,pw_list,vmin=vmin,vmax=vmax)
    plt.colorbar(tmp,label='pitch angle [°]')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_title(r'plausibility check of pitch angle calculation ($\phi_B=' + str(B_phi) + r'°,~\theta_B=' + str(B_theta) + r'°$)')

    plt.savefig(f'plausibility_check/plausibility_check_phi_{B_phi}_theta_{B_theta}.png')

# for angles in fov:
#     plausibility_check(angles[0],angles[1])



### Funktion als Maske benutzen ###

a = np.ones((30,20))
# print(a)

for i in range(20):
    a[:,i]*=i
# print(a)

def g(t):
    return t+1

# m = a < g(np.arange(30))[:,np.newaxis]
# print(m)
# print(g(np.arange(30))[:,np.newaxis])
# print(g(np.arange(30)))
# print(a[m])
# b = 1.0*a
# a[~m] = 0    # Tilde invertiert den boolschen Wert???
# print(a[~m])
# tmp = plt.pcolormesh(a.T)
# plt.colorbar(tmp)
# plt.show()


ebins = np.array([0.98 ,2.144,2.336,2.544,2.784,3.04,3.312])
data = np.array([[200,200,50,50,50,300],[200,50,50,50,300,50],[200,50,50,300,50,50],[50,50,300,50,50,50],[50,50,300,50,50,50],[50,50,300,50,50,50],[50,300,50,50,50,50],[50,300,50,50,50,50],[50,300,50,50,50,50],[300,50,50,50,50,50],[300,50,50,50,50,50],[300,50,50,50,50,50]]).T
print(data.shape)
print(data)

# Wie genau funktioniert die  Grenzfunktion???

def create_masks(grenzfunktion,shape):
    '''Übergebe eine Grenzfunktion als Funktion der Indizes der zweiten array-Dimension (Zeitreihe). Die Grenzfunktion berechnet dann den Index in der ersten array-Dimension (Energie-Bins),
    ab dem die Daten beginnen sollen. Beachte: Die Energie-Bins werden im STEP-Plot von unten nach oben gezählt.'''
    len_ebins = shape[0]
    len_time = shape[1]
    a = np.ones((len_ebins,len_time))
    for i in range(len_ebins):
        a[i,:]*=i    # Setze in erster array-Dimension die Werte auf den entsprechenden Index
    mask = a < grenzfunktion(np.arange(len_time))[np.newaxis,:]
    return mask

def grenz(i):
    return -i**0.5 + 2

mask = create_masks(grenz,data.shape)
print(mask)
print(~mask)
data[mask] = 0
print(data)


def func():
    pass
print(type(func)==func)



tmp = plt.pcolormesh(data)
plt.colorbar(tmp)
# plt.yscale('log')
# plt.show()


# b = 1.*a
# b[~m] = 0 

# Energies = array([1,3,6,...,100])
# Earr = zeros(Datenarry.shape)
# Earr[:] = Energies
# MittlereEnergien = (Datenarray * Earr).sum(axis = 1)/Datenarray.sum(axis = 1)