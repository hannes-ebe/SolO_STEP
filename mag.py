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
        self.pitchangles = []
        self.pitchangles_err = []
        # Blickrichtung der einzelnen Pixel (Phi, Theta):
        self.fov = np.array([[-25,24],[-25,12],[-25,0],[-25,-12],[-25,-24],[-35,24],[-35,12],[-35,0],[-35,-12],[-35,-24],[-45,24],[-45,12],[-45,0],[-45,-12],[-45,-24]])
        
        self._load()
        self._calc_angles()
        self._calc_pw()
        # self._calc_pw_err()
        
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
        
    def pw(self,v_phi,v_theta):
        '''Es ist wichtig die Normierung zu beachten. Wenn ich die RTN-Koordinaten über die Winkel beschreibe, kann ich die Vorfaktoren 1 setzen.
        Numpy nimmt die Winkel in Radians. Muss den Spezialfall von Werten größer 1 im arccos abfangen Dieser Sonderfall ist hier erstmal ignoriert. 
        Für das Magnetfeld, werden die Daten des Objekts genommen.
        Der Winkel \phi gibt die Drehung in der RT-Ebene und \theta in der RN-Ebene.
        Für die Umrechnung müsste mit v als Länge des Vektors gelten:
            R = v*cos(phi)*cos(theta)
            T = v*sin(phi)*cos(theta)
            N = v*sin(theta)
        Die Attribute phi und theta beziehen sich bei Lars auf das Magnetfeld.'''
        v_theta = v_theta/360*2*np.pi
        v_phi = v_phi/360*2*np.pi
        B_theta = self.theta/360*2*np.pi
        B_phi = self.phi/360*2*np.pi
        argument = np.cos(v_theta)*np.cos(v_phi)*np.cos(B_theta)*np.cos(B_phi) + np.cos(v_theta)*np.sin(v_phi)*np.cos(B_theta)*np.sin(B_phi) + np.sin(v_theta)*np.sin(B_theta)
        # if argument > 1.0 and v_theta == B_theta and v_phi == B_phi:
        #     argument = 1.0
        pitchangle = np.arccos(argument)/np.pi*180
        return pitchangle
        
    def _calc_pw(self):
        '''Berechne die Pitchwinkel für die Elektronen, welche auf STEP treffen in erster Näherung.
        Dafür wird der Winkel zwischen den Blickrichtungen der Pixel und dem Magnetfeld herangezogen.'''
        for i in range(15):
            v_phi = self.fov[i][0]
            v_theta = self.fov[i][1]
            self.pitchangles.append(self.pw(v_phi,v_theta))
        self.pitchangles = np.array(self.pitchangles)
        
        
    def pw_err(self,v_phi,v_theta):
        ''' Implementierung der Fehler der Pitchwinkel. Nehme erstmal an, dass der Fehler des absoluten Feldes von 0.2 nT
        nur in in B_N-Richtung auftritt. Die Magnetfelder sind in nT gegeben.'''
        B_N_err = 0.2
        B_theta_err = 1/(1+self.B_N*self.B_N/(self.B_T*self.B_T+self.B_R*self.B_R)) * 1/np.sqrt(self.B_T*self.B_T+self.B_R*self.B_R) * B_N_err
        v_theta = np.radians(v_theta)
        v_phi = np.radians(v_phi)
        B_theta = np.radians(self.theta)
        B_phi = np.radians(self.phi)
        
        pw_err = np.abs((np.sin(B_phi)*np.sin(v_phi)+np.cos(B_phi)*np.cos(v_phi))*np.cos(v_theta)*np.sin(B_theta)-np.sin(v_theta)*np.cos(B_theta)/np.sqrt(1-(np.sin(v_theta)*np.sin(B_theta)+np.sin(B_phi)*np.sin(v_phi)*np.cos(v_theta)*np.cos(B_theta)+np.cos(B_phi)*np.cos(v_phi)*np.cos(v_theta)*np.cos(B_theta))**2))*B_theta_err
        return pw_err
        
        
    def _calc_pw_err(self):
        '''Berechne die Fehler der Pitchwinkel mittels der Implementierung unter pw_err.''' 
        for i in range(15):
            v_phi = self.fov[i][0]
            v_theta = self.fov[i][1]
            self.pitchangles_err.append(self.pw_err(v_phi,v_theta))
        self.pitchangles_err = np.array(self.pitchangles_err)
        
        
        
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
    
    
    def pw_ts(self,err=False):
        '''Bei err == True werden Fehler mit dargestellt.'''

        fig, ax = self.step_plot(r'time', r'$\varphi$ [°]', 'Time series of pitch angle from ' + str(self.period[0]) + ' to ' + str(self.period[1]))

        for i in range(1,16):
            if err == True:
                ax[i].errorbar(self.time,self.pitchangles[i-1],yerr=self.pitchangles_err[-1],marker='x')
            else:
                ax[i].scatter(self.time,self.pitchangles[i-1],marker='x')
            ax[i].tick_params(axis='x', labelrotation=90)
        
        plt.show()
        
