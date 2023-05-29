import matplotlib as mpl
import pylab as plt
import numpy as np
import datetime as dt
import matplotlib.dates as mdates
import load_nom_II as ld
import plot_nom_II as pt
import math
import mag

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



class STEP_Events(STEP):
    def __init__(self,year,month,day,rpath = '/data/projects/solo/step_v0008/',rpath_mag = '/data/projects/solo/mag/l2_soar/rtn_1minute',magnet_default_path=False,lastofmonth=False):
        super().__init__(year,month,day,rpath=rpath,rpath_mag=rpath_mag,magnet_default_path=magnet_default_path,lastofmonth=lastofmonth)
    
    def sum_pixel(self,ebins,head,period,norm):
        '''Summiere erst über alle Pixel.'''
        
        pldat, pltime, vmax = self.data_prep(ebins=ebins,head=head,period=period,norm=norm)
        
        '''Was hat Lars sich bei der nächsten Zeile gedacht???'''
        # ptime = np.append(pltime[1],pltime[1][-1]+dt.timedelta(seconds=60))
        ptime = pltime[1]
        
        # Hintergrund-Pixel wird weggelassen
        summe = pldat[1]

        for i in range(2,16):
            summe += pldat[i]
            
        return summe, pldat, ptime
            
    def sum_energy(self,summe, pldat):
        '''Summiere über ebins.'''
        total_sum = []
        for i in summe:
            total_sum.append(sum(i))
   
        total_sum = np.array(total_sum)
        
        return total_sum
    
    def plot_total_sum(self,time,total_sum,period):
        '''Plotte die Summe über alle Pixel und ebins gegen die Zeit.'''
        fig, ax = plt.subplots(figsize=(10,8))
        
        ax.plot(time,total_sum)
        plt.title(f'Time series of sum of all pixels and ebins from {period[0]} to {period[1]}')
        ax.set_ylabel('sum in a.u.')
        ax.set_xlabel('time')
        ax.tick_params(axis='x',labelrotation=45)
        
        plt.show()
        
    def total_sum_wrapper(self,ebins=ebins[1:-1],head=-1,period=(dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)),norm='ptmax'):
        '''Wrapper für totale Summe der Energie'''
        summe, pldat, pltime = self.sum_pixel(ebins=ebins,head=head,period=period,norm=norm)
        total_sum = self.sum_energy(summe,pldat)
        self.plot_total_sum(pltime,total_sum,period)
            
            
dat = STEP_Events(2021,12,4)
dat.total_sum_wrapper(ebins=ebins[18:-1])        


    # def plot(self):
    #     fig, ax = plt.subplots()
        
    #     #plt.xlabel()
    #     #plt.ylabel()
        
    #     tmp = ax.pcolormesh(ptime,ebins,sum.T,cmap=cmap,vmin=0.1)
        
    #     plt.show()
        