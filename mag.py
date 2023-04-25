######################################################################
# Author : Lars Berger                                               #
# Date : 23.09.2020                                                  #
#                                                                    #
# Load SOLO/MAG data currently from CDF files provided by Tim H.     #
# directly, later CSV format might/should follow to allow working    #
# without use of a cdf interface                                     #
######################################################################

import numpy as np
import datetime as dt
import pylab as pl
import os
from cdflib import CDF,cdfepoch
import matplotlib.pyplot as plt

class MAGdata(object):
    def __init__(self, path = '/data/projects/solo/mag/l2_soar/rtn_1minute', period = (dt.datetime(2022,1,1),dt.datetime(2022,1,2)), res = 'min'):
        '''Nutze automatisch Level 2 Daten. Wenn ich einen Tag laden will, muss ich bei period den nächsten Tag als obere Grenze angeben.
        Lade die minütlichen Daten. Ich habe auch nicht ewig Zeit. Für 8 Hz-Daten muss res = None.'''
    
        self.path = path
        self.period = period
        self.res = res           # res = None -> 8 Hz data
        self.B_R = np.zeros((0))
        self.B_T = np.zeros((0))
        self.B_N = np.zeros((0))
        self.time = []
        # Blickrichtung der einzelnen Pixel (Phi, Theta):
        self.fov = np.array([[-25,24],[-25,12],[-25,0],[-25,-12],[-25,-24],[-35,24],[-35,12],[-35,0],[-35,-12],[-35,-24],[-45,24],[-45,12],[-45,0],[-45,-12],[-45,-24]])
        
        self._load()
        self._calc_angles()
        self._calc_pw()
        
    def _load(self):
        d = self.period[0]
        while d < self.period[1]:
            y = d.year
            m = d.month
            dd = d.day
            if not self.res:
                l = os.listdir(self.path + '/%.4i/'%y)
                l.sort()
                ll = ['solo_L2_mag-rtn-normal_%.4i%.2i%.2i'%(y,m,dd) in s for s in l]
                fname =''
                for i,n in enumerate(ll):
                    if n:
                        fname = self.path + '/%.4i/'%(y) + l[i]

                print(d,fname)
            if self.res == 'min':
                l = os.listdir(self.path + '/%.4i/'%y)
                l.sort()
                ll = ['solo_L2_mag-rtn-normal-1-minute_%.4i%.2i%.2i'%(y,m,dd) in s for s in l]
                fname =''
                for i,n in enumerate(ll):
                    if n:
                        fname = self.path + '/%.4i/'%(y) + l[i]

                print(d,fname)
            try:
                dat = CDF(fname)
                t = dat.varget('EPOCH')
                t = cdfepoch.encode_tt2000(t)
                for tt in t:
                    self.time.append(dt.datetime.fromisoformat(str(tt)[:26]))
                    self.time[-1] = self.time[-1].replace(tzinfo = None)
                b = dat.varget('B_RTN')
                self.B_R = np.append(self.B_R,np.array(b[:,0]))
                self.B_T = np.append(self.B_T,np.array(b[:,1]))
                self.B_N = np.append(self.B_N,np.array(b[:,2]))
            except:
                print('No MAG data : ', d)
            d = d + dt.timedelta(days = 1)
        self.time = np.array(self.time).flatten()
        self.B_R = np.array(self.B_R).flatten()
        self.B_T = np.array(self.B_T).flatten()
        self.B_N = np.array(self.B_N).flatten()

    def _calc_angles(self):
        '''Winkel der Magnetfeldvektoren im Spacecraft-Frame'''
        self.phi = np.arctan2(self.B_T,self.B_R)
        self.theta = np.arctan2(self.B_N,np.sqrt(self.B_T**2+self.B_R**2))
        
    def _calc_pw(self):
        '''Berechne die Pitchwinkel für die Elektronen, welche auf STEP treffen in erster Näherung.
        Dafür wird der Winkel zwischen den Blickrichtungen der Pixel und
        dem Magnetfeld herangezogen.'''
        # Lars hat etwas von Normierung geschrieben... Warum??? Müsste es nicht reichen die Winkel im Spacecraft-Frame zu nehmen???
        # Nutze theta und phi für Kugelkoordinaten wie bei Wiki...
        # Umrechnung in Grad nicht vergessen...
        B_phi = self.phi
        B_theta = self.theta
        self.pw = []
        
        for i in range(15):
            v_theta = np.radians(self.fov[i][1])   # Umrechnung in radians für numpy
            v_phi = np.radians(self.fov[i][0])
            # Umrechnung in Grad ist bedenken
            self.pw.append((np.arccos(np.sin(v_theta)*np.cos(v_phi)*np.sin(B_theta)*np.cos(B_phi) + np.sin(v_theta)*np.sin(v_phi)*np.sin(B_theta)*np.sin(B_phi) + np.cos(v_theta)*np.cos(B_theta)))/np.pi*180)
        
        self.pw = np.array(self.pw)
        
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
    
    def pw_ts(self):

        fig, ax = self.step_plot(r'time', r'$\varphi$ [°]', 'Time series of pitch angle from ' + str(self.period[0]) + ' to ' + str(self.period[1]))

        for i in range(1,16):
            ax[i].scatter(self.time,self.pw[i-1],marker='x')
            ax[i].tick_params(axis='x', labelrotation=90)
        
        plt.show()
        
