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



class STEP_Events():
    def __init__(self,year,month,day,rpath = '/data/projects/solo/step_v0008/',rpath_mag = '/data/projects/solo/mag/l2_soar/rtn_1minute',lastofmonth=False):
        '''Magnetfeld wird gleich mitgeladen.'''
        # Loading data
        if lastofmonth:
            if month!=12:
                self.itime, self.idata = ld.load_nom(rpath=rpath,period=(dt.datetime(year,month,day),dt.datetime(year,month+1,1)), products=('M','A'))
            else:
                self.itime, self.idata = ld.load_nom(rpath=rpath,period=(dt.datetime(year,month,day),dt.datetime(year+1,1,1)), products=('M','A'))
        else:
            self.itime, self.idata = ld.load_nom(rpath=rpath,period=(dt.datetime(year,month,day),dt.datetime(year,month,day+1)), products=('M','A'))
        print('STEP-Data loaded successfully.')
        
        # Combining data (Main and Auxiliary Product)
        ld.combine_data(self.itime, self.idata)
        print('STEP-Data combined successfully.')
        
        # Loading MAG-Data
        # if lastofmonth:
        #     if month!=12:
        #         self.mag = mag.MAGdata(path = rpath_mag, period = (dt.datetime(year,month,day),dt.datetime(year,month+1,1)))
        #     else:
        #         self.mag = mag.MAGdata(path = rpath_mag, period = (dt.datetime(year,month,day),dt.datetime(year+1,1,1)))
        # else:
        #     self.mag = mag.MAGdata(path = rpath_mag, period=(dt.datetime(year,month,day),dt.datetime(year,month,day+1)))

    def cut_data(self,t0,t1):
        cdat = {}
        ctime = {}
        for p in self.idata.keys():
            m = (self.itime[p]>=t0) * (self.itime[p]<t1)
            ctime[p] = self.itime[p][m] 
            cdat[p] = self.idata[p][m] 
        return ctime, cdat

    def data_prep(self,ebins=ebins,res = '1min', head = 0, period = None, norm = False, overflow = True, esquare = False):
        '''Vorbereitung der STEP-Daten basierend auf Lars Skripten.'''
        vmax = None # Falls keine Normierung angegben wird, gebe ich None zurück. Was soll vmax eigentlich darstellen?

        if period:
            time,dat = self.cut_data(period[0]-dt.timedelta(seconds=59),period[1])
        else:
            time,dat = self.itime, self.idata
        pldat = []
        pltime = []
        if type(res) == dt.timedelta:
            for i in range(16):
                if head in [0,1]:
                    tmptime = [time['STEP_C_%i_%.2i'%(head,i)][0]]
                    while tmptime[-1]<time['STEP_C_%i_%.2i'%(head,i)][-1]:
                        tmptime.append(tmptime[-1]+res)
                if head in [0,1]:
                    tdat = 1.*dat['STEP_C_%i_%.2i'%(head,i)]
                    ttime = time['STEP_C_%i_%.2i'%(head,i)]
                    tmpdat = []
                    for en in range(tdat.shape[1]):
                        H,t = np.histogram(ttime,bins = tmptime,weights = tdat[:,en])
                        tmpdat.append(H)
                    tmpdat = np.array(tmpdat).T
                    pldat.append(tmpdat)
                    pltime.append(np.array(tmptime[:-1]))
                else:
                    pldat.append(1.*dat['STEP_C_0_%.2i'%(i)]-dat['STEP_C_1_%.2i'%(i)])
        elif res == '1min':
            for i in range(16):
                if head in [0,1]:
                    pldat.append(1.*dat['STEP_C_%i_%.2i'%(head,i)])
                    pltime.append(time['STEP_C_%i_%.2i'%(head,i)])
                else:
                    pldat.append(1.*dat['STEP_C_0_%.2i'%(i)]-dat['STEP_C_1_%.2i'%(i)])
                    pltime.append(time['STEP_C_0_%.2i'%(i)])
                if esquare:
                    pldat[-1]= pldat[-1]/(np.diff(ebins)**2)[np.newaxis,:]
                    #pldat[-1]= pldat[-1]/ np.diff(ebins)[np.newaxis,:]
        elif res == '1s':
            for i in range(16):
                if head in [0,1]:
                    pldat.append(dat['STEP_C_%i_%.2i'%(head,i)]/60.)
                    pltime.append(time['STEP_C_%i_%.2i'%(head,i)])
                else:
                    pldat.append(dat['STEP_C_0_%.2i'%(i)]/60.-dat['STEP_C_1_%.2i'%(i)]/60.)
                    pltime.append(time['STEP_C_0_%.2i'%(i)])
                if esquare:
                    pldat[-1]= pldat[-1]/(np.diff(ebins)**2)[np.newaxis,:]
            for i in range(16):
                if head in [0,1]:
                    pldat.append(dat['STEP_M_%i_%.2i'%(head,i)])
                    pltime.append(time['STEP_M_%i_%.2i'%(head,i)])
                else:
                    pldat.append(dat['STEP_M_0_%.2i'%(i)]-dat['STEP_M_1_%.2i'%(i)])
                    pltime.append(time['STEP_M_0_%.2i'%(i)])
                if esquare:
                    pldat[-1]= pldat[-1]/(np.diff(ebins)[8:40]**2)[np.newaxis,:]
        if not overflow:
            for i in range(len(pldat)):
                pldat[i][:,-1] = 0.
        if norm == 'tmax':
            vmax = np.zeros(pltime[0].shape)
            for i in range(len(pldat)):
                vpmax = np.amax(pldat[i],axis=1)
                vmax[vpmax>vmax] = vpmax[vpmax>vmax]
        elif norm == 'ptmax' or norm == 'ptemax':
            vmax = np.zeros((16,pltime[0].shape[0]))
            for i in range(16):
                vpmax = np.amax(pldat[i],axis=1)
                vmax[i][vpmax>vmax[i]] = vpmax[vpmax>vmax[i]]
        elif norm == 'max' or norm == 'logmax':
            vmax = 0.
            for i in range(len(pldat)):
                vpmax = np.amax(pldat[i])
                vmax = max(vpmax,vmax)
        elif norm == 'pemax':
            vmax = np.zeros((16,pldat[i].shape[1]))
            for i in range(16):
                vpmax = np.amax(pldat[i],axis=0)
                vmax[i] = vpmax
        elif norm == 'emax':
            vmax = np.zeros(pldat[i].shape[1])
            for i in range(16):
                vpmax = np.amax(pldat[i],axis=0)
                vmax[vpmax>vmax] = vpmax[vpmax>vmax]

        return pldat, pltime, vmax
    
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
        