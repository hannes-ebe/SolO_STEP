'''Funktionen um STEP-Daten zu analysieren.'''

import matplotlib as mpl
import pylab as plt
import numpy as np
import datetime as dt
import matplotlib.dates as mdates
import load_nom_II as ld
import plot_nom_II as pt

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

    def __init__(self,year,month,day,lastofmonth=False):
        # Loading data
        if lastofmonth:
            if month!=12:
                self.itime, self.idata = ld.load_nom(period=(dt.datetime(year,month,day),dt.datetime(year,month+1,1)), products=('M','A'))
            else:
                self.itime, self.idata = ld.load_nom(period=(dt.datetime(year,month,day),dt.datetime(year+1,1,1)), products=('M','A'))
        else:
            self.itime, self.idata = ld.load_nom(period=(dt.datetime(year,month,day),dt.datetime(year,month,day+1)), products=('M','A'))
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

        if period:
            time,dat = self.cut_data(self.idata,self.itime,period[0]-dt.timedelta(seconds=59),period[1])
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

    def marginal_distribution(self,ebins=ebins,res = '1min', head = 0, pixel = 0, period = None, save = False, norm = False, overflow = True, esquare = False):
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
                    
            # To Do...
            # Wie sieht Struktur von pdat aus??? Wie genau mache ich die Histogramme???
            # pdat == pldat[pixel]
            print(pldat[pixel])
            print(len(pldat[pixel][0]))
            print(ebins)
            print(len(pldat[pixel].T[0]))
            print(ptime)

            if i == 1:  
                # Projektion auf Energie-Achse  
                pdat = pldat[pixel]    
                ax.append(fig.add_subplot(3,1,i+1))
                # ax[-1].set_xscale('log')
                print(np.sum(pdat,axis=0))
                ax[-1].hist(np.sum(pdat,axis=0),bins=ebins)
            
            if i == 2:
                # Projektion auf Zeit-Achse
                pdat = pldat[pixel]    
                ax.append(fig.add_subplot(3,1,i+1))
                print(np.sum(pdat,axis=1))
                ax[-1].hist(np.sum(pdat,axis=1),bins=ptime)
                
                    
            if (norm and i == 0): # or not norm:
                tax = fig.add_subplot(4,50,50)
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

            if i == 0:
                tl = [l.get_text() for l in ax[i].get_xticklabels()]
                ax[i].set_title('Head %i'%head)
            # ax[i].text(0.05,0.05,'%i'%(i),transform=ax[i].transAxes,color = 'k', backgroundcolor = 'w', fontsize = 6)
            if i == 0:
                ax[i].hlines(ebins[8],ptime[0],ptime[-1])
                ax[i].hlines(ebins[40],ptime[0],ptime[-1])
                ax[i].tick_params(axis='x',labelrotation=45)
                ax[i].set_xlabel('Date')
                ax[i].set_ylabel('Energy [keV]')

            if i == 1 or i == 2:
                if norm == 'tmax':
                    ax[i].set_ylabel('Counts(t)/max(Counts(t)')
                elif norm == 'ptmax':
                    ax[i].set_ylabel('Counts(t)/max(Counts(t)')
                elif norm == 'max':
                    ax[i].set_ylabel('Counts')
                elif norm == 'logmax':
                    ax[i].set_ylabel('Log10(C)')
                elif norm == 'pemax':
                    ax[i].set_ylabel('C(E,pixel)/max(C(E,pixel)')
                elif norm == 'emax':
                    ax[i].set_ylabel('C(E)/max(C(E)')

        if period:
            ax[0].set_xlim(period[0],period[1])
        if save:
            if type(save) == str:
                plt.savefig(save + 'marginal_pixel%i_'%pixel + 'head%i_'%head + 'TS_%.4i_%.2i_%.2i_%.2i-%.2i-%.2i-%i_%.2i-%.2i-%.2i_H%i_%s_%s.png'%(ptime[0].year,ptime[0].month,ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
                #plt.savefig(save + 'head%i/'%head + 'TS_%i_%.2i:%.2i:%.2i-%i_%.2i:%.2i:%.2i_H%i_%s_%s.pdf'%(ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
            else:
                plt.savefig('TS_%.4i_%.2i_%.2i_%.2i-%.2i-%.2i-%i_%.2i-%.2i-%.2i_H%i_%s_%s.png'%(ptime[0].year,ptime[0].month,ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
                #plt.savefig('TS_%i_%.2i:%.2i:%.2i-%i_%.2i:%.2i:%.2i_H%i_%s_%s.pdf'%(ptime[0].day,ptime[0].hour,ptime[0].minute,ptime[0].second,ptime[-2].day,ptime[-2].hour,ptime[-2].minute,ptime[-2].second,head,norm,res))
        print('Analyzed marginal distribution successfully.')