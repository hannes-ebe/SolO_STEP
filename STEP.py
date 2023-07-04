import matplotlib as mpl
import pylab as plt
import numpy as np
import datetime as dt
import matplotlib.dates as mdates
import load_nom_II as ld
import plot_nom_II as pt
import math
import mag
import matplotlib.gridspec as gridspec

from PIL import Image
import glob
import os

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


class STEP():
    def __init__(self,year,month,day,rpath = '/data/projects/solo/step_v0008/',rpath_mag = None,magnet_default_path=False,lastofmonth=False):
        '''Magnetfeld wird gleich mitgeladen, wenn rpath_mag übergeben wird. Wird magnet_default_path auf True gesetzt wird automatisch der 
        korrekte Dateipfad für die Uni-Rechner genutzt.'''
        self.ebins = ebins = np.array([  0.98 ,   2.144,   2.336,   2.544,   2.784,   3.04 ,   3.312,
         3.6  ,   3.92 ,   4.288,   4.672,   5.088,   5.568,   6.08 ,
         6.624,   7.2  ,   7.84 ,   8.576,   9.344,  10.176,  11.136,
        12.16 ,  13.248,  14.4  ,  15.68 ,  17.152,  18.688,  20.352,
        22.272,  24.32 ,  26.496,  28.8  ,  31.36 ,  34.304,  37.376,
        40.704,  44.544,  48.64 ,  52.992,  57.6  ,  62.72 ,  68.608,
        74.752,  81.408,  89.088,  97.28 , 105.984, 115.2  , 125.44 ,
       137.216, 149.504, 162.816, 178.176, 194.56 , 211.968, 230.4  ,
       372.736])
        
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
        if magnet_default_path == True:
            self.magnet_default_path = '/data/projects/solo/mag/l2_soar/rtn_1minute'
            rpath_mag = self.magnet_default_path

        if type(rpath_mag) == str:
            if lastofmonth:
                if month!=12:
                    self.mag = mag.MAGdata(path = rpath_mag, period = (dt.datetime(year,month,day),dt.datetime(year,month+1,1)))
                else:
                    self.mag = mag.MAGdata(path = rpath_mag, period = (dt.datetime(year,month,day),dt.datetime(year+1,1,1)))
            else:
                self.mag = mag.MAGdata(path = rpath_mag, period=(dt.datetime(year,month,day),dt.datetime(year,month,day+1)))

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
    
    def data_boxes(self,pldat,box_list):
        '''Nehme pldat und, um die Daten auf die Boxen einzuschränken. Dabei werden die Werte, die ignoriert werden sollen durch 0.0 ersetzt.
        box_data wird zurückgegeben.'''
        # Erster Index von pdat müsste die Zeitreihe sein, der zweite der Energie-Bin.
        # [[[timelow,timeup],[energylow,energyup]],[[timelow,timeup],[energylow,energyup]]]
        box_data = []
        for pdat in pldat:
            ayuda = np.zeros(shape=(len(pdat),len(pdat.T)),dtype='float')
            for box in box_list:
                tlow = box[0][0]
                tup = box[0][1]
                elow = box[1][0]
                eup = box[1][1]

                # Beschleunigung durch Nutzung von Masken???
                # Genauer ansehen, wenn ich in Laufzeitprobleme laufe...
                # Code auf jeden Fall schon so angepasst, dass data_boxes time bekommt
                # mask = ()

                # Loope erst durch die entsprechende Zeitreihe und dann jeweils durch die Energiebins
                for t in range(tlow,tup+1):
                    for e in range(elow,eup+1):
                        ayuda[t][e] = pdat[t][e]
            box_data.append(ayuda)
        return box_data
    
    def create_masks(self,grenzfunktion,shape):
        '''Übergebe eine Grenzfunktion als Funktion der Indizes der ersten array-Dimension (Zeitreihe). Die Grenzfunktion berechnet dann den Index in der nullten array-Dimension (Energie-Bins),
        ab dem die Daten beginnen sollen. Beachte: Die Energie-Bins werden im STEP-Plot von unten nach oben gezählt.'''
        len_ebins = shape[0]
        len_time = shape[1]
        a = np.ones((len_ebins,len_time))
        for i in range(len_ebins):
            a[i,:]*=i    # Setze in erster array-Dimension die Werte auf den entsprechenden Index
        mask = a < grenzfunktion(np.arange(len_time))[np.newaxis,:]
        return mask
    
    def set_zero(self,pldat,mask):
        '''Übergebe eine Maske und die Daten aller Pixel. Werte, die laut Maske True sind, werden auf 0 gesetzt.'''
        for pdat in pldat:
            pdat[mask] = 0
        return pldat


    def step_plot(self,xlabel,ylabel,title):
        '''Erzeugt den bekannten STEP-Plot mit allen 16 Pixeln. Gibt Liste mit axes zurück.'''
        fig = plt.figure(figsize = (15,10))
        fig.subplots_adjust(wspace = 0, hspace = 0)
        ax = []        
        for i in range(16):
            if i == 0:
                ax.append(fig.add_subplot(4,5,3))
            else:
                if i == 1:
                    ax.append(fig.add_subplot(4,5,5+i))
                else:
                    ax.append(fig.add_subplot(4,5,5+i,sharex = ax[1], sharey = ax[1]))              
                if i == 6:
                    ax[i].set_ylabel(ylabel)
                if i == 13:
                    ax[i].set_xlabel(xlabel)
            if i not in [0,1,6,11]:
                for t in ax[i].get_yticklabels():
                    t.set_visible(False)
            if i < 11:
                for t in ax[i].get_xticklabels():
                    t.set_visible(False)
        ax[0].set_title(title)
        return fig, ax
    
    def plot_ts(self,ebins=ebins,res = '1min', head = 0, period = None, save = False, norm = False, overflow = True, esquare = False, grenzfunktion = None, box_list= False):
        '''Plottet die Zeitreihen der Energien. Übergebe verschachtelte Liste mit Grenzen für Boxen, die eingezeichnet werden. Grenzen als Indizes der Listen übergeben.
        Wenn sowohl box_list als auch grenzfunktion gegeben sind, wird grenzfunktion genutzt und eingezeichnet.'''
        pldat, pltime, vmax = self.data_prep(ebins,res,head,period,norm,overflow,esquare)
        fig, ax = self.step_plot('Date', 'Energy [keV]', 'Head %i'%head)
        for i in range(16):
            pdat = pldat[i]
            ptime = np.append(pltime[i],pltime[i][-1]+dt.timedelta(seconds=60))

            ax[i].set_yscale('log')
            if type(res) == dt.timedelta:
                if not norm:
                    tmp = ax[i].pcolormesh(ptime,ebins,pdat.T, cmap = cmap, vmin = 0.1)
                elif norm == 'tmax':
                    tmp = ax[i].pcolormesh(ptime,ebins,pdat.T/vmax, cmap = cmap, vmin = np.amin(1/vmax)*0.95,vmax = 1.)
                elif norm == 'ptmax':
                    tmp = ax[i].pcolormesh(ptime,ebins,pdat.T/vmax[i], cmap = cmap, vmin = np.amin(1/vmax)*0.95,vmax = 1.)
                elif norm == 'ptemax':
                    pdat2 = (pdat.T/vmax[i]).T
                    vemax = np.amax(pdat2,axis=0)
                    print(pdat.shape)
                    print(vemax.shape)
                    pdat2 = (pdat2/vemax)
                    pdat2[np.isinf(pdat2)] = np.NaN
                    tmp = ax[i].pcolormesh(ptime,ebins,pdat2.T, cmap = cmap, vmin = np.nanmin(pdat2),vmax = 1.)
                elif norm == 'max':
                    tmp = ax[i].pcolormesh(ptime,ebins,pdat.T, cmap = cmap, vmin = 0.,vmax = vmax)
                elif norm == 'logmax':
                    tmp = ax[i].pcolormesh(ptime,ebins,np.log10(pdat).T, cmap = cmap, vmin = 0.,vmax = np.log10(vmax))
                elif norm == 'pemax':
                    tmp = ax[i].pcolormesh(ptime,ebins,(pdat/vmax[i]).T, cmap = cmap, vmin = np.amin(1/vmax[i])*0.99,vmax = 1.)
                elif norm == 'emax':
                    tmp = ax[i].pcolormesh(ptime,ebins,(pdat/vmax).T, cmap = cmap, vmin = np.amin(1/vmax)*0.99,vmax = 1.)

            elif res == '1min':
                if not norm:
                    tmp = ax[i].pcolormesh(ptime,ebins,pdat.T, cmap = cmap, vmin = 0.1)
                elif norm == 'tmax':
                    tmp = ax[i].pcolormesh(ptime,ebins,pdat.T/vmax, cmap = cmap, vmin = np.amin(1/vmax)*0.95,vmax = 1.)
                elif norm == 'ptmax':
                    tmp = ax[i].pcolormesh(ptime,ebins,pdat.T/vmax[i], cmap = cmap, vmin = np.amin(1/vmax)*0.95,vmax = 1.)
                elif norm == 'ptemax':
                    pdat2 = (pdat.T/vmax[i]).T
                    vemax = np.amax(pdat2,axis=0)
                    print(pdat.shape)
                    print(vemax.shape)
                    pdat2 = (pdat2/vemax)
                    pdat2[np.isinf(pdat2)] = np.NaN
                    tmp = ax[i].pcolormesh(ptime,ebins,pdat2.T, cmap = cmap, vmin = np.nanmin(pdat2),vmax = 1.)
                elif norm == 'max':
                    tmp = ax[i].pcolormesh(ptime,ebins,pdat.T, cmap = cmap, vmin = 0.,vmax = vmax)
                elif norm == 'logmax':
                    tmp = ax[i].pcolormesh(ptime,ebins,np.log10(pdat).T, cmap = cmap, vmin = 0.,vmax = np.log10(vmax))
                elif norm == 'pemax':
                    tmp = ax[i].pcolormesh(ptime,ebins,(pdat/vmax[i]).T, cmap = cmap, vmin = np.amin(1/vmax[i])*0.99,vmax = 1.)
                elif norm == 'emax':
                    tmp = ax[i].pcolormesh(ptime,ebins,(pdat/vmax).T, cmap = cmap, vmin = np.amin(1/vmax)*0.99,vmax = 1.)

            elif res == '1s':
                if not norm:
                    tmp = ax[i].pcolormesh(ptime,ebins[:9],pdat[:,:8].T, cmap = cmap, vmin = 0.1)
                    tmp = ax[i].pcolormesh(ptime,ebins[40:],pdat[:,40:].T, cmap = cmap, vmin = 0.1)
                elif norm == 'tmax':
                    tmp = ax[i].pcolormesh(ptime,ebins,pdat.T/vmax, cmap = cmap, vmin = np.amin(1/vmax)*0.99,vmax = 1.)
                elif norm == 'ptmax':
                    tmp = ax[i].pcolormesh(ptime,ebins,pdat.T/vmax[i], cmap = cmap, vmin = np.amin(1/vmax)*0.99,vmax = 1.)
                elif norm == 'max':
                    tmp = ax[i].pcolormesh(ptime,ebins,pdat.T, cmap = cmap, vmin = 0.,vmax = vmax)
                elif norm == 'logmax':
                    tmp = ax[i].pcolormesh(ptime,ebins[:9],np.log10(pdat[:,:8]).T, cmap = cmap, vmin = np.log10(1./60.),vmax = np.log10(vmax))
                    tmp = ax[i].pcolormesh(ptime,ebins[40:],np.log10(pdat[:,40:]).T, cmap = cmap, vmin = np.log10(1./60.),vmax = np.log10(vmax))
                elif norm == 'pemax':
                    tmp = ax[i].pcolormesh(ptime,ebins,(pdat/vmax[i]).T, cmap = cmap, vmin = np.amin(1/vmax[i])*0.99,vmax = 1.)
                elif norm == 'emax':
                    tmp = ax[i].pcolormesh(ptime,ebins,(pdat/vmax).T, cmap = cmap, vmin = np.amin(1/vmax)*0.99,vmax = 1.)
                pdat = pldat[i+16]
                ptime = np.append(pltime[i+16],pltime[i+16][-1]+dt.timedelta(seconds=1))
                if not norm:
                    tmp = ax[i].pcolormesh(ptime,ebins[8:41],pdat.T, cmap = cmap, vmin = 0.1)
                elif norm == 'tmax':
                    tmp = ax[i].pcolormesh(ptime,ebins,pdat.T/vmax, cmap = cmap, vmin = np.amin(1/vmax)*0.99,vmax = 1.)
                elif norm == 'ptmax':
                    tmp = ax[i].pcolormesh(ptime,ebins,pdat.T/vmax[i], cmap = cmap, vmin = np.amin(1/vmax)*0.99,vmax = 1.)
                elif norm == 'max':
                    tmp = ax[i].pcolormesh(ptime,ebins,pdat.T, cmap = cmap, vmin = 0.,vmax = vmax)
                elif norm == 'logmax':
                    tmp = ax[i].pcolormesh(ptime,ebins[8:41],np.log10(pdat).T, cmap = cmap, vmin = np.log10(1./60.),vmax = np.log10(vmax))
                elif norm == 'pemax':
                    tmp = ax[i].pcolormesh(ptime,ebins,(pdat/vmax[i]).T, cmap = cmap, vmin = np.amin(1/vmax[i])*0.99,vmax = 1.)
                elif norm == 'emax':
                    tmp = ax[i].pcolormesh(ptime,ebins,(pdat/vmax).T, cmap = cmap, vmin = np.amin(1/vmax)*0.99,vmax = 1.)
                
            if callable(grenzfunktion):
                x_ind = np.arange(len(ptime))
                y_ind = (np.rint(grenzfunktion(x_ind))).astype(int)
                ayuda_mask = (y_ind < len(ebins)-1) * (y_ind >= 0)
                x = ptime[ayuda_mask]
                y = ebins[y_ind[ayuda_mask]]
                ax[i].step(x,y,where='post',color='tab:red')
            
            # PLots der Boxen, Grenzen als  Indizes übergeben.
            # Erster Index von pdat müsste die Zeitreihe sein, der zweite der Energie-Bin.
            # [[[timelow,timeup],[energylow,energyup]],[[timelow,timeup],[energylow,energyup]]]
            if type(box_list) == list:
                for box in box_list:
                    tlow = ptime[box[0][0]]
                    tup = ptime[box[0][1]]
                    elow = ebins[box[1][0]]
                    eup = ebins[box[1][1]]
                    ax[i].vlines([tlow,tup],elow,eup,color='firebrick')
                    ax[i].hlines([elow,eup],tlow,tup,color='firebrick')
                
            if (norm and i == 0): # or not norm:
                tax = fig.add_subplot(4,50,31)
                if norm == 'tmax':
                    plt.colorbar(tmp,cax = tax, label = 'Counts(t)/max(Counts(t)')
                elif norm == 'ptmax':
                    plt.colorbar(tmp,cax = tax, label = 'Counts(t)/max(Counts(t)')
                elif norm == 'max':
                    plt.colorbar(tmp,cax = tax, label = 'Counts')
                elif norm == 'logmax':
                    plt.colorbar(tmp,cax = tax, label = 'Log10(C)')
                elif norm == 'pemax':
                    plt.colorbar(tmp,cax = tax, label = 'C(E,pixel)/max(C(E,pixel)')
                elif norm == 'emax':
                    plt.colorbar(tmp,cax = tax, label = 'C(E)/max(C(E)')

            ax[i].hlines(ebins[8],ptime[0],ptime[-1])
            ax[i].hlines(ebins[40],ptime[0],ptime[-1])            
            if i > 10:
                ax[i].tick_params(labelrotation=90)
            if period:
                ax[0].set_xlim(period[0],period[1])
        if save:
            if type(save) == str:
                plt.savefig(save + 'TS_%.4i_%.2i_%.2i_%.2i-%.2i-%.2i-%i_%.2i-%.2i-%.2i_H%i_%s_%s.png'%(ptime[0].year,ptime[0].month,ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
                #plt.savefig(save + 'head%i/'%head + 'TS_%i_%.2i:%.2i:%.2i-%i_%.2i:%.2i:%.2i_H%i_%s_%s.pdf'%(ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
            else:
                plt.savefig('TS_%.4i_%.2i_%.2i_%.2i-%.2i-%.2i-%i_%.2i-%.2i-%.2i_H%i_%s_%s.png'%(ptime[0].year,ptime[0].month,ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
                #plt.savefig('TS_%i_%.2i:%.2i:%.2i-%i_%.2i:%.2i:%.2i_H%i_%s_%s.pdf'%(ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
        print('Plotted TS_%.4i_%.2i_%.2i_%.2i-%.2i-%.2i-%i_%.2i-%.2i-%.2i_H%i_%s_%s.png'%(ptime[0].year,ptime[0].month,ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res) + ' successfully.')
        
        
    def calc_pw(self,period,window_width=5):
        '''Berechnung der Pitchwinkel'''
        # Maske, da ich nur die Magnetfelddaten innerhalb von period brauche:
        mask = (self.mag.time > period[0]) * (self.mag.time <= period[1])
        pw = [[] for i in range(15)]
        pw_time = []
        
        i = 0
        while (period[0] + dt.timedelta(minutes=(i+1)*window_width)) <= period[1]:
            pw_time.append(period[0] + dt.timedelta(minutes=(i+0.5)*window_width))
            
            for k in [i for i in range(1,16)]:
                # Mittelung der Pitchwinkel (k-1, da ich keine Zeit im array stehen habe)
                pw_data = self.mag.pitchangles[k-1][mask]
                new_pw = np.sum(pw_data[i*window_width:(i+1)*window_width])/window_width
                pw[k-1].append(new_pw)
            i +=1
        return pw, pw_time

    
    
    def calc_energy_means(self,ebins=ebins,res = '1min', head = -1, period = (dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)), grenzfunktion=None, below=True, box_list=None, window_width=5, pixel_list=[i for i in range(1,16)], norm = 'tmax', overflow = True, esquare = False):
        '''Wichtig: Die Energie ist in keV angegeben!!!
        Falls andere Zeitauflösung als Minuten gewählt wird, kann es Probleme mit window_width und der Berechnung der Zeit geben.
        below gehört zu grenzfunktion und gibt an, ob die Werte über oder unter der Funktion auf Null gesetzt werden.
        Wenn sowohl box_list als auch  grenzfunktion gegeben sind, wird grenzfunktion genutzt.'''
        i = 0
        pixel_means = [[] for i in range(16)]     # Liste mit Listen der Mittelwerte der einzelnen Pixel. Die erste Liste enthält die Zeitstempel (jeweils Mitte der Zeitfenster)
        
        if callable(grenzfunktion):
            '''Bereite zunächst little_helper_dat vor und erzeuge dann eine Maske, die ja auf alle Pixel angewendet werden kann.
            Anschließend wird set_zero die eine Maske und ein array mit allen Pixeln übergeben, um die entsprechenden Werte auf 
            Null zu setzen.'''
            little_helper_dat = np.array(self.data_prep(ebins,res,head,period,norm,overflow,esquare)[0])
            # Transponieren scheint hier notwendig zu sein, damit Zeit und Energie den korrekten array-Dimensionen entsprechen...
            # Nutze ersten Eintrag in little_helper_dat, shape sollte für alle Pixel gleich sein.
            ayuda_shape = (little_helper_dat[0].shape[1],little_helper_dat[0].shape[0])
            mask = self.create_masks(grenzfunktion=grenzfunktion,shape=ayuda_shape)
            if below == False:
                mask = ~mask  # Tilde invertiert (logical-not)
            pldat = self.set_zero(little_helper_dat,mask.T)
            
        elif type(box_list) == list:
            little_helper_dat = self.data_prep(ebins,res,head,period,norm,overflow,esquare)[0]
            pldat = self.data_boxes(little_helper_dat,box_list)
        else:
            pldat = self.data_prep(ebins,res,head,period,norm,overflow,esquare)[0]

        # Mittelwerte der einzelnen Energie-Bins
        mittelwerte_energie_bins = []
        for k in range(len(ebins)-1):
            mittelwerte_energie_bins.append(0.5*(ebins[k+1]+ebins[k]))
        mittelwerte_energie_bins = np.array(mittelwerte_energie_bins)

        while (period[0] + dt.timedelta(minutes=(i+1)*window_width)) <= period[1]:
            pixel_means[0].append(period[0] + dt.timedelta(minutes=(i+0.5)*window_width))
            
            # Alte Berechnung der Mittelwerte:

            # for k in pixel_list:
            #     pdat = pldat[k][i*window_width:(i+1)*window_width]
            #     integral = np.sum(pdat,axis=0)
            #     # calculating mean
            #     mean = 0.0
            #     for j in range(len(integral)):
            #         # Für Bestimmung des Mittelwertes der Verteilung wird Mitte der Bins gewählt
            #         mean += integral[j]*0.5*(ebins[j+1]+ebins[j])
            #     mean = mean/np.sum(integral)
            #     pixel_means[k].append(mean)
            # i +=1


            # Neue Berechnung der Mittelwerte 
            
            for k in pixel_list:
                pdat = pldat[k][i*window_width:(i+1)*window_width]
                integral = np.sum(pdat,axis=0)
                # calculating mean
                mean = np.sum(mittelwerte_energie_bins*integral)/np.sum(integral)
                pixel_means[k].append(mean)
            i +=1
            
        
        for i in range(len(pixel_means)):
            # Array erstellen
            pixel_means[i] = np.array(pixel_means[i])
            
        return pixel_means
    
    def calc_sliding_energy_means(self,ebins=ebins,res = '1min', head = -1, period = (dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)), grenzfunktion=None, below=True, box_list=None, window_width=5, sliding_width=4, pixel_list=[i for i in range(1,16)], norm = 'tmax', overflow = True, esquare = False):
        '''Wichtig: Die Energie ist in keV angegeben!!!
        Falls andere Zeitauflösung als Minuten gewählt wird, kann es Probleme mit window_width und der Berechnung der Zeit geben.
        Nehme immer vier benachbarte energy-bins und berechne den Mittelwert (Sowohl der normierten counts und der Energie). Dieses Fenster gleitet über alle bins. 
        Nehme dann aus diesen Werten den maximalen Wert. Versuche so die Effekte der Daten-Boxen zu minimieren, weil die Tail-Länge an der Seite den Mittelwert beeinflusst.
        Wenn sowohl box_list als auch  grenzfunktion gegeben sind, wird grenzfunktion genutzt.'''
        i = 0
        pixel_means = [[] for i in range(16)]     # Liste mit Listen der Mittelwerte der einzelnen Pixel. Die erste Liste enthält die Zeitstempel (jeweils Mitte der Zeitfenster)
        
        if callable(grenzfunktion):
            '''Bereite zunächst little_helper_dat vor und erzeuge dann eine Maske, die ja auf alle Pixel angewendet werden kann.
            Anschließend wird set_zero die eine Maske und ein array mit allen Pixeln übergeben, um die entsprechenden Werte auf 
            Null zu setzen.'''
            little_helper_dat = np.array(self.data_prep(ebins,res,head,period,norm,overflow,esquare)[0])
            # Transponieren scheint hier notwendig zu sein, damit Zeit und Energie den korrekten array-Dimensionen entsprechen...
            # Nutze ersten Eintrag in little_helper_dat, shape sollte für alle Pixel gleich sein.
            ayuda_shape = (little_helper_dat[0].shape[1],little_helper_dat[0].shape[0])
            mask = self.create_masks(grenzfunktion=grenzfunktion,shape=ayuda_shape)
            if below == False:
                mask = ~mask  # Tilde invertiert (logical-not)
            pldat = self.set_zero(little_helper_dat,mask.T)
            
        elif type(box_list) == list:
            little_helper_dat = self.data_prep(ebins,res,head,period,norm,overflow,esquare)[0]
            pldat = self.data_boxes(little_helper_dat,box_list)
        else:
            pldat = self.data_prep(ebins,res,head,period,norm,overflow,esquare)[0]

        while (period[0] + dt.timedelta(minutes=(i+1)*window_width)) <= period[1]:
            pixel_means[0].append(period[0] + dt.timedelta(minutes=(i+0.5)*window_width))
            
            # Berechnung der Mittelwerte:
            for k in pixel_list:
                pdat = pldat[k][i*window_width:(i+1)*window_width]
                integral = np.sum(pdat,axis=0)
                # calculating mean
                # Breite des Sliding-windows sollte natürlicher Teiler der 56 energy-bins sein und wird in Zahl der Indizes angegeben.
                j_indizes=[i*sliding_width for i in range(int(len(integral)/sliding_width))]
                j_indizes=[i for i in range(0,len(integral)-sliding_width)]
                means_energy = [] # mittlere Energien
                means_normedcounts = [] # Mittelung der normierten Counts im Window
                for j in j_indizes:
                    # Für Bestimmung des Mittelwertes der Verteilung wird Mitte der Bins gewählt
                    mean = 0.0
                    count =  0.0
                    for l in range(sliding_width):
                        mean += integral[j+l]*0.5*(ebins[j+l+1]+ebins[j+l])
                        count += integral[j+l]*(ebins[j+l+1]-ebins[j+l])      # Multipliziere hier nochmal mit der Bin-Breite; Keine Normierung auf eins, damit ich nachher max, count wählen kann
                    mean = mean/np.sum(integral[j:j+sliding_width])
                    means_energy.append(mean)
                    means_normedcounts.append(count)
                total_mean = means_energy[means_normedcounts.index(max(means_normedcounts))]  # liefert mir Energiemittelwert für maximalen normierten Count beim Sliding-Window
                pixel_means[k].append(total_mean)
            i +=1
        
        for i in range(len(pixel_means)):
            # Array erstellen
            pixel_means[i] = np.array(pixel_means[i])
            
        return pixel_means
    
    def plot_energy_means(self,pixel_means, rpath, pixel_list=[i for i in range(1,16)]):
        '''Plottet die übergebenen Mittelwerte für alle Pixel.'''
        fig, ax = plt.subplots(figsize=(10,6))
        for i in pixel_list:
            ax.scatter(pixel_means[0],pixel_means[i],marker='x',label=f'pixel {i}')
        ax.set_ylabel('mean of energy [keV]')
        ax.set_xlabel('time')
        ax.tick_params(axis='x',labelrotation=45)
        plt.title('energy means')
        plt.legend()
        plt.savefig(rpath)
        
    def pixel_comparison(self, pixel_means, pw, pw_time, rpath, pixel1, pixel2):
        '''Vergleich zweier Pixel mit Energie-Mittelwerten und Pitchwinkel.
        Die Zeit der Energiemittelwerte steckt in pixel_means. Für die Pitchwinkel
        wird sie extra übergeben.'''
        plt.figure(figsize=(10,12))
        gs = gridspec.GridSpec(2,1,height_ratios=[3,1])
        ax = plt.subplot(gs[0])
        ax_pw = ax.twinx()
        ax_quot = plt.subplot(gs[1],sharex=ax)
        ax_quot_pw = ax_quot.twinx()
        for i in [pixel1, pixel2]:
            ax.plot(pixel_means[0],pixel_means[i],marker='x',label=f'energy mean of pixel {i}')
            ax_pw.plot(pixel_means[0],pw[i-1],marker="^",label=f'pitch angle of pixel {i}')
        ax.set_ylabel('mean of energy [keV]')
        ax_pw.set_ylabel('pitch angle [°]')


        # Arrays für die Betrachtung der Quotienten
        # Mal sehen, ob ich eine brauchbare Korrektur finden kann...
        
        pixel_means1_ayuda = np.sqrt(np.array(pixel_means[pixel1]))
        pixel_means2_ayuda = np.sqrt(np.array(pixel_means[pixel2]))

        pw1_ayuda = np.cos(np.radians(np.array(pw[pixel1-1])))
        pw2_ayuda = np.cos(np.radians(np.array(pw[pixel2-1])))

        ratio_means = pixel_means1_ayuda/pixel_means2_ayuda
        ratio_pw = pw2_ayuda/pw1_ayuda


        ax_quot.plot(pixel_means[0],ratio_means,marker='x',label='quotient of energy means',c='tab:red')
        ax_quot_pw.plot(pw_time,ratio_pw,marker='^',label='quotient of cosine of pitch angles',c='tab:green')
        ax_quot.set_ylabel('quotient of energy means')
        ax_quot_pw.set_ylabel('quotient of cosine of pitch angles')
        ax_quot.set_xlabel('time')
        ax_quot.tick_params(axis='x',labelrotation=45)
        ax.set_title('Comparison of energy means and pitch angles')
        ax.legend(loc='upper left')
        ax_pw.legend(loc='upper right')
        ax_quot.legend(loc='upper left')
        ax_quot_pw.legend(loc='upper right')
        plt.subplots_adjust(hspace=0.001)
        plt.savefig(rpath)
        plt.close('all')
        
        
                
    def energy_correction(self,energy,pw1_degree,pw2_degree):
        '''Korrigiert die übergebene Energie über die übergebenen Pitchwinkel und die implementierte Korrektur.
        energy entspricht der Energie des Pixels, der pw2 gesehen hat.'''
        
        pw1 = np.radians(pw1_degree)
        pw2 = np.radians(pw2_degree)
        corr = np.cos(pw2)*np.cos(pw2)/np.cos(pw1)/np.cos(pw1)
        return corr*energy
        
        
        
    def pixel_comparison_corrected(self, pixel_means, pw, rpath, pixel1, pixel2):
        '''Plot der originalen Energiemittelwerte sowie der Korrektur.'''
        
        plt.figure(figsize=(10,9))
        gs = gridspec.GridSpec(2,1,height_ratios=[2,1])
        
        ax = plt.subplot(gs[0])
        
        ax.plot(pixel_means[0],pixel_means[pixel1],marker='x',label=f'energy mean of pixel {pixel1}')
        ax.plot(pixel_means[0],pixel_means[pixel2],marker='x',label=f'energy mean of pixel {pixel2}')
        
        energy2_corrected = self.energy_correction(pixel_means[pixel2],pw[pixel1-1],pw[pixel2-1])
        
        ax.plot(pixel_means[0],energy2_corrected,marker='x',label=f'corrected energy mean of pixel {pixel2}')
        
        ax.set_ylabel('mean of energy [keV]')
        ax.legend()
        ax.set_title('Correction of energy means')
        
        # Plots der Differenzen
        ax2 = plt.subplot(gs[1],sharex=ax)
        
        diff_original = pixel_means[pixel2] - pixel_means[pixel1]
        diff_corrected = energy2_corrected - pixel_means[pixel1]
        
        ax2.plot(pixel_means[0],diff_original,marker='x',label='difference of energy means')
        ax2.plot(pixel_means[0],diff_corrected,marker='x',label='difference of energy means with correction')
        
        ax2.set_ylabel('difference of energy means [keV]')
        ax2.set_xlabel('time')
        ax2.tick_params(axis='x',labelrotation=45)
        ax2.legend()
        
        plt.subplots_adjust(hspace=0.001)
        plt.savefig(rpath)
        plt.close('all')

        
        
    def distribution_ring(self, filename, title, head, norm, period, grenzfunktion=None,below=True,box_list=None, norm_pixel=3, correction = False, save='gif/', res = '1min', overflow = True, esquare = False,window_width = 5, close=True, sorted_by_energy=False):
        '''Darstellung von means der einzelnen Pixel und Pitchwinkel als GIF. Es soll die ringförmige Verteilung deutlich werden. Code basiert auf minütlichen Daten!!!'''
        
        # Berechnung der Energie-Mittelwerte und Mittelung der Pitchwinkel über Intervalle der Länge window_width
        
        # Zunächst alle Bilddateien aus gif/-Ordner löschen, um Probleme beim erstellen der Gifs zu vermeiden.
        dir = "gif_images"
        filelist = glob.glob(os.path.join(dir,"*"))
        for f in filelist:
            os.remove(f)
        
        
        # Berechnung der Energiemittelwerte
        # pixel_list kann bei Aufruf von distribution_ring nicht verändert werden, da ich eh alle Pixel brauche...
        pixel_means = self.calc_energy_means(ebins=ebins,res=res,head=head,period=period,grenzfunktion=grenzfunktion,below=below,box_list=box_list,window_width=window_width,pixel_list=[i for i in range(1,16)],norm=norm,overflow=overflow,esquare=esquare)
        
        # Berechnung der Pitchwinkel
        pw, pw_time = self.calc_pw(period=period,window_width=window_width)

            
        if sorted_by_energy:
            ind = np.argsort(pixel_means[norm_pixel])
        else:
            ind = np.array([i for i in range(len(pixel_means[norm_pixel]))])
            
        norm_factor = np.array(pixel_means[norm_pixel]) 
        vmin = 1.0
        vmax = 1.0
        vmin_pw = 1.0
        vmax_pw = 1.0
        
        for k in range(0,16):
            if k == 0:
                pixel_means[k] = np.array(pixel_means[k])
            else:
                pixel_means[k] = np.array(pixel_means[k])
                pw[k-1] = np.array(pw[k-1])
        
        for k in range(1,16):
            if correction:
                corr = self.correction_pw(pw[norm_pixel-1],pw[k-1])
            else:
                corr = 1.0
            pixel_means[k] = pixel_means[k]/norm_factor*corr
            pw[k-1] = pw[k-1]*corr
            if np.nanmax(pixel_means[k]) > vmax:
                vmax = np.nanmax(pixel_means[k])
            if np.nanmin(pixel_means[k]) < vmin:
                vmin = np.nanmin(pixel_means[k])
            if np.nanmax(pw[k-1]) > vmax_pw:
                vmax_pw = np.nanmax(pw[k-1])
            if np.nanmin(pw[k-1]) < vmin_pw:
                vmin_pw = np.nanmin(pw[k-1])
                                  
                    
        
        # Plotting
        x_corners = [0,1,2,3,4,5]
        y_corners = [0,1,2,3]
        for k in range(len(pixel_means[1])):
            
            a = [pixel_means[j][ind][k] for j in range(1,6)]
            b = [pixel_means[j][ind][k] for j in range(6,11)]
            c = [pixel_means[j][ind][k] for j in range(11,16)]        
            data_means = np.array([c,b,a])
            
            d = [pw[j-1][ind][k] for j in range(1,6)]
            e = [pw[j-1][ind][k] for j in range(6,11)]
            f = [pw[j-1][ind][k] for j in range(11,16)]  
            pw_list = np.array([f,e,d])
            
            
            fig, axes = plt.subplots(2,1,figsize=(10,12))
            
            ax = axes[0]
            
            # Nehme betraglich größte Differenz für colorbar
            v_ayuda = max(abs(vmax-1.0),abs(1.0-vmin))
            
            tmp = ax.pcolormesh(x_corners,y_corners,data_means,vmin=1.0-v_ayuda,vmax=1.0+v_ayuda,cmap=hmap)
            plt.colorbar(tmp,label=f'mean of energy/(mean of energy of pixel {norm_pixel})')
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            
            if sorted_by_energy == True:
                ax.set_title(title + f'\n(head {head}, {window_width} minute steps, normed to pixel {norm_pixel})\nfrom ' + str(period[0]) + ' to ' + str(period[1]) + f'\n Mean Energy of Pixel {norm_pixel}: ' + str(round(norm_factor[ind][k],2)) + ' [keV]')
            else:
                ax.set_title(title + f'\n(head {head}, {window_width} minute steps, normed to pixel {norm_pixel})\nfrom ' + str(period[0]) + ' to ' + str(period[1]) + '\n Time: ' + str(pixel_means[0][ind][k]))
            
            
            ax = axes[1]
            tmp = ax.pcolormesh(x_corners,y_corners,pw_list,vmin=vmin_pw,vmax=vmax_pw)
            plt.colorbar(tmp,label=r'pitch angle $\varphi$ [°]')
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.set_title(r'Corresponding pitch angles $\varphi$ for the pixels')
    
    
    
            # Quick and dirty führende Nullen für korrekte Sortierung
            if k < 10:
                plt.savefig('gif_images/image0' + str(k) + '.png')
            else:
                plt.savefig('gif_images/image' + str(k) + '.png')
            
            if close:
                plt.close('all')
        
        # Erstellen des GIF's
        frames = []
        imgs = glob.glob("gif_images/image*.png")
        imgs.sort()
        for i in imgs:
            new_frame = Image.open(i)
            frames.append(new_frame)

        # Save the png images into a GIF file that loops forever
        frames[0].save(save + f'{filename}.gif', format='GIF', append_images=frames[1:], save_all=True, duration=500, loop=0)
        print('Created GIF successfully!')