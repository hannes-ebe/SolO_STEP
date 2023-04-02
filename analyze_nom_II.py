'''Funktionen um STEP-Daten zu analysieren.'''

import matplotlib as mpl
import pylab as plt
import numpy as np
import datetime as dt
import matplotlib.dates as mdates
import load_nom_II as ld
import plot_nom_II as pt

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


class STEP:
    '''Lädt STEP-Daten zur anschließenden Analyse'''

    def __init__(self,rpath,year,month,day,lastofmonth=False):
        # Loading data
        if lastofmonth:
            if month!=12:
                self.itime, self.idata = ld.load_nom(rpath=rpath,period=(dt.datetime(year,month,day),dt.datetime(year,month+1,1)), products=('M','A'))
            else:
                self.itime, self.idata = ld.load_nom(rpath=rpath,period=(dt.datetime(year,month,day),dt.datetime(year+1,1,1)), products=('M','A'))
        else:
            self.itime, self.idata = ld.load_nom(rpath=rpath,period=(dt.datetime(year,month,day),dt.datetime(year,month,day+1)), products=('M','A'))
        print('Data loaded successfully.')
        
        # Combining data (Main and Auxiliary Product)
        ld.combine_data(self.itime, self.idata)
        print('Data combined successfully.')

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
                # Loope erst durch die entsprechende Zeitreihe und dann jeweils durch die Energiebins
                for t in range(tlow,tup+1):
                    for e in range(elow,eup+1):
                        ayuda[t][e] = pdat[t][e]
            box_data.append(ayuda)
        return box_data
            

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
    
    def plot_ts(self,ebins=ebins,res = '1min', head = 0, period = None, save = False, norm = False, overflow = True, esquare = False, box_list= False):
        '''Plottet die Zeitreihen der Energien. Übergebe verschachtelte Liste mit Grenzen für Boxen, die eingezeichnet werden. Grenzen als Indizes der Listen übergeben.'''
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
    
    def pixel_integral_window(self, filename, ebins=ebins,res = '1min', head = 0, period = None, save = False, norm = False, overflow = True, esquare = False):
        '''Funktion, die Daten für alle 16 Pixel bekommt und diese darstellt. Übergebe ein Zeitfenster, über welches die jeweiligen Energien integriert werden.
        save gibt den directory-Pfad an unter dem gespeichert wird.'''
        # Vorbereitung Daten
        pldat, pltime, vmax = self.data_prep(ebins,res,head,period,norm,overflow,esquare)
        
        # Plotting
        fig = plt.figure(figsize = (15,10))
        fig.subplots_adjust(wspace = 0, hspace = 0)
        ax = []
        means = []
        for i in range(16):
            pdat = pldat[i]
            if i == 0:
                ax.append(fig.add_subplot(4,5,3))
            else:
                if i == 1:
                    ax.append(fig.add_subplot(4,5,5+i,sharex = ax[0]))
                else:
                    ax.append(fig.add_subplot(4,5,5+i,sharex = ax[0], sharey = ax[1]))
                ax[-1].set_xscale('log')
                integral = np.sum(pdat,axis=0)
                ax[-1].step(ebins[1:],integral,where='pre')
                # calculating mean
                mean = 0.0
                for j in range(len(integral)):
                    # Für Bestimmung des Mittelwertes der Verteilung wird Mitte der Bins gewählt
                    mean += integral[j]*0.5*(ebins[j+1]+ebins[j])
                mean = mean/np.sum(integral)
                means.append(mean)
                ax[i].axvline(mean,color='blue') #,label='Mean')
                
                if i == 6:
                    ax[i].set_ylabel('integral along time')
                if i == 13:
                    ax[i].set_xlabel('Energy [keV]')
                ax[i].axvline(ebins[8],color='firebrick')
                ax[i].axvline(ebins[40],color='firebrick') #,label='Energy range of STEP')
                # ax[i].legend()
            if i not in [0,1,6,11]:
                for t in ax[i].get_yticklabels():
                    t.set_visible(False)
            if i < 11:
                for t in ax[i].get_xticklabels():
                    t.set_visible(False)
        # Plot der Means:
        ax[0].set_xscale('log')
        ax[0].set_ylabel('pixel')
        ax[0].scatter(means,[i for i in range(1,16)],marker='x',label='Mean')
        ax[0].legend()
        ax[0].text(1.5, 0.5,'Red Lines: Energy range of STEP\nBlue Line: Mean of energy distribution',transform=ax[0].transAxes)
        ax[0].set_title('Energy distribution and mean for head ' + str(head) + ' from ' + str(period[0]) + ' to ' + str(period[1]))
        if save:
            if type(save) == str:
                plt.savefig(save + filename)
            else:
                plt.savefig(filename)
        print('Plotted ' + filename + ' successfully.')

    def evolution_energy_means(self, filename, ebins=ebins,res = '1min', head = 0, period = None, window_width = 5, save = False, norm = False, overflow = True, esquare = False, box_list = False):
        '''Übergebe period und Zeitfensterbreite. Intergriere die jeweiligen Energien in den Zeitfenstern und stelle die Entwicklung der means für die einzelnen Pixel da.'''
        i = 0
        pixel_means = [[] for i in range(16)]     # Liste mit Listen der Mittelwerte der einzelnen Pixel. Die erste Liste enthält die Zeitstempel (jeweils Mitte der Zeitfenster)
        
        while (period[0] + dt.timedelta(minutes=(i+1)*window_width)) <= period[1]:
            window = [period[0] + dt.timedelta(minutes=i*window_width), period[0] + dt.timedelta(minutes=(i+1)*window_width)]

            if type(box_list) == list:
                little_helper = self.data_prep(ebins,res,head,window,norm,overflow,esquare)[0]
                pldat = self.data_boxes(little_helper,box_list)
            else:
                pldat = self.data_prep(ebins,res,head,window,norm,overflow,esquare)[0]

            pixel_means[0].append(period[0] + dt.timedelta(minutes=(i+0.5)*window_width))
            
            # Berechnung der Mittelwerte:
            for k in range(1,16):
                pdat = pldat[k]
                integral = np.sum(pdat,axis=0)
                # calculating mean
                mean = 0.0
                for j in range(len(integral)):
                    # Für Bestimmung des Mittelwertes der Verteilung wird Mitte der Bins gewählt
                    mean += integral[j]*0.5*(ebins[j+1]+ebins[j])
                mean = mean/np.sum(integral)
                pixel_means[k].append(mean)
            i +=1
        
        fig, ax = self.step_plot('Time', 'Mean of energy distribution [keV]', 'Evolution of mean of energy distribution for head ' + str(head) + ' from ' + str(period[0]) + ' to ' + str(period[1]) + ' (' + str(window_width) + ' minute steps)')
        for i in range(1,16):
            ax[i].plot(pixel_means[0],pixel_means[i])
            ax[i].set_yscale('log')
            if i > 10:
                ax[i].tick_params(axis='x',labelrotation=90)
            ax[i].axhline(ebins[8],color='firebrick')
            ax[i].axhline(ebins[40],color='firebrick') #,label='Energy range of STEP')
            
        if save:
            if type(save) == str:
                plt.savefig(save + filename)
            else:
                plt.savefig(filename)
        print('Plotted ' + filename + ' successfully.')            



    def landau(self,x,A,B,C,D):
        return A/np.sqrt(2*np.pi)*np.exp(-B*0.5*((x+C) + np.exp(-(x+C)))) + D

    def landau_fit(self,xdata,ydata,p0):
        popt, pcov = curve_fit(self.landau,xdata,ydata,p0)
        return popt,pcov

    def marginal_distribution(self,ebins=ebins,res = '1min', head = 0, pixel = 0, period = None, save = False, norm = False, overflow = True, esquare = False, fit = False):
        '''Plot des Histogramms eines einzelnen Pixels und der Projektionen auf Zeit- und Energieachse'''

        pldat, pltime, vmax = self.data_prep(ebins,res,head,period,norm,overflow,esquare)

        fig = plt.figure(figsize = (8,15))
        ax = []
        for i in range(3):
            if i == 0:
                pdat = pldat[pixel]
                ptime = np.append(pltime[pixel],pltime[pixel][-1]+dt.timedelta(seconds=60))
                ax.append(fig.add_subplot(3,1,1))

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

            if i == 1:  
                # Projektion auf Energie-Achse  
                pdat = pldat[pixel]    
                ax.append(fig.add_subplot(3,1,i+1))
                ax[-1].set_xscale('log')
                ydata = np.sum(pdat,axis=0)
                ax[-1].step(ebins[1:],ydata,where='pre')
                if fit:
                    xdata = ebins[1:] - np.diff(ebins)
                    popt,pcov = self.landau_fit(xdata,ydata)
                    print('Landau-Fit:')
                    print('Parameter: ', popt)
                    print('Kovarianz: ', pcov)
            
            if i == 2:
                # Projektion auf Zeit-Achse
                pdat = pldat[pixel]    
                ax.append(fig.add_subplot(3,1,i+1))
                ax[-1].step(ptime[1:],np.sum(pdat,axis=1),where='pre')
                
                    
            if (norm and i == 0): # or not norm:
                tax = fig.add_subplot(4,50,50)
                if norm == 'tmax':
                    plt.colorbar(tmp,cax = tax, label = 'Counts(t)/max(Counts(t))')
                elif norm == 'ptmax':
                    plt.colorbar(tmp,cax = tax, label = 'Counts(t)/max(Counts(t))')
                elif norm == 'max':
                    plt.colorbar(tmp,cax = tax, label = 'Counts')
                elif norm == 'logmax':
                    plt.colorbar(tmp,cax = tax, label = 'Log10(C)')
                elif norm == 'pemax':
                    plt.colorbar(tmp,cax = tax, label = 'C(E,pixel)/max(C(E,pixel))')
                elif norm == 'emax':
                    plt.colorbar(tmp,cax = tax, label = 'C(E)/max(C(E))')

            if i == 0:
                tl = [l.get_text() for l in ax[i].get_xticklabels()]
                ax[i].set_title('Head %i'%head+'/pixel %i'%pixel)
            # ax[i].text(0.05,0.05,'%i'%(i),transform=ax[i].transAxes,color = 'k', backgroundcolor = 'w', fontsize = 6)
            if i == 0:
                ax[i].hlines(ebins[8],ptime[0],ptime[-1])
                ax[i].hlines(ebins[40],ptime[0],ptime[-1])
                ax[i].tick_params(axis='x',labelrotation=0)
                ax[i].set_xlabel('Date')
                ax[i].set_ylabel('Energy [keV]')

            if i == 1:
                ax[i].set_ylabel('sum along date axis')
                ax[i].set_xlabel('Energy [keV]')
                ax[i].axvline(ebins[8],color='firebrick')
                ax[i].axvline(ebins[40],color='firebrick',label='energy range of STEP')
                ax[i].legend()
            if i == 2:
                ax[i].set_ylabel('sum along energy axis')
                ax[i].set_xlabel('Date')
                # if norm == 'tmax':
                #     ax[i].set_ylabel('Counts(t)/max(Counts(t)')
                # elif norm == 'ptmax':
                #     ax[i].set_ylabel('Counts(t)/max(Counts(t)')
                # elif norm == 'max':
                #     ax[i].set_ylabel('Counts')
                # elif norm == 'logmax':
                #     ax[i].set_ylabel('Log10(C)')
                # elif norm == 'pemax':
                #     ax[i].set_ylabel('C(E,pixel)/max(C(E,pixel)')
                # elif norm == 'emax':
                #     ax[i].set_ylabel('C(E)/max(C(E)')

        if period:
            ax[0].set_xlim(period[0],period[1])
        if save:
            if type(save) == str:
                plt.savefig(save + 'marginal_pixel%i_'%pixel + 'head%i_'%head + 'TS_%.4i_%.2i_%.2i_%.2i-%.2i-%.2i-%i_%.2i-%.2i-%.2i_H%i_%s_%s.png'%(ptime[0].year,ptime[0].month,ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
                #plt.savefig(save + 'head%i/'%head + 'TS_%i_%.2i:%.2i:%.2i-%i_%.2i:%.2i:%.2i_H%i_%s_%s.pdf'%(ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
            else:
                plt.savefig('marginal_pixel%i_'%pixel + 'head%i_'%head + 'TS_%.4i_%.2i_%.2i_%.2i-%.2i-%.2i-%i_%.2i-%.2i-%.2i_H%i_%s_%s.png'%(ptime[0].year,ptime[0].month,ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
                #plt.savefig('TS_%i_%.2i:%.2i:%.2i-%i_%.2i:%.2i:%.2i_H%i_%s_%s.pdf'%(ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
        print('Analyzed marginal distribution successfully.')

    def fit_energy(self,ebins=ebins,res = '1min', head = 0, pixel = 0, period = None, save = False, norm = False, overflow = True, esquare = False, fit = False, p0 = None):
        pldat, pltime, vmax = self.data_prep(ebins,res,head,period,norm,overflow,esquare)
        fig = plt.figure(figsize = (8,5))
        pdat = pldat[pixel]
        ptime = np.append(pltime[pixel],pltime[pixel][-1]+dt.timedelta(seconds=60))

        if res == '1s':
            # Eigentlich Index i+16... Erstmal nicht relevant...
            pdat = pldat[0+16]
            ptime = np.append(pltime[0+16],pltime[0+16][-1]+dt.timedelta(seconds=1))

        # Projektion auf Energie-Achse  
        pdat = pldat[pixel]    
        ax = fig.add_subplot(1,1,1)
        # ax[-1].set_xscale('log')
        ydata = np.sum(pdat,axis=0)
        ax.step(ebins[1:],ydata,where='pre')

        ax.set_ylabel('sum of counts along date axis')
        ax.set_xlabel('Energy [keV]')
        ax.axvline(ebins[8],color='firebrick')
        ax.axvline(ebins[40],color='firebrick',label='energy range of STEP')
        ax.legend()
        ax.set_xlim(ebins[0],ebins[-13])
    
        if fit:
            xdata = ebins[1:] - np.diff(ebins)
            popt,pcov = self.landau_fit(xdata,ydata,p0)
            print('Landau-Fit:')
            print('Parameter: ', popt)
            print('Kovarianz: ', pcov)
            A = round(popt[0],2)
            B = round(popt[1],2)
            C = round(popt[2],2)
            D = round(popt[3],2)
            xlin = np.linspace(ebins[0],ebins[-1],1000)
            ax.plot(xlin,self.landau(xlin,A,B,C,D),color='orange')
            ax.set_title(r'Landau-Fit: $C(E)\approx\frac{%s'%A+r'}{\sqrt{2\pi}}\exp{(-\frac{%s'%B+r'}{2}((E+%s'%C+r')+\exp{(-(E+%s'%C+r')))})}+%s'%D+r'$', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        else:
            ax.set_title('Landau-Fit')
               
        if save:
            if type(save) == str:
                plt.savefig(save + 'energy_fit_pixel%i'%pixel + 'head%i_'%head + 'TS_%.4i_%.2i_%.2i_%.2i-%.2i-%.2i-%i_%.2i-%.2i-%.2i_H%i_%s_%s.png'%(ptime[0].year,ptime[0].month,ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
                #plt.savefig(save + 'head%i/'%head + 'TS_%i_%.2i:%.2i:%.2i-%i_%.2i:%.2i:%.2i_H%i_%s_%s.pdf'%(ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
            else:
                plt.savefig('energy_fit_pixel%i'%pixel + 'head%i_'%head + 'TS_%.4i_%.2i_%.2i_%.2i-%.2i-%.2i-%i_%.2i-%.2i-%.2i_H%i_%s_%s.png'%(ptime[0].year,ptime[0].month,ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
                #plt.savefig('TS_%i_%.2i:%.2i:%.2i-%i_%.2i:%.2i:%.2i_H%i_%s_%s.pdf'%(ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
        print('Fitted energy successfully.')
        if fit:
            return popt,pcov