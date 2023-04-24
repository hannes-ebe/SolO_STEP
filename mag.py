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

class MAGdata(object):
    def __init__(self, path = '/data/projects/solo/mag/l2_soar/rtn_high_time_res/', lvl='L2', period = (dt.datetime(2022,1,1),dt.datetime(2022,1,2)), res = None):
        self.path = path
        self.lvl = lvl
        self.period = period
        self.res = res           # res = None -> 8 Hz data
        self.B_R = np.zeros((0))
        self.B_T = np.zeros((0))
        self.B_N = np.zeros((0))
        self.time = []
        self._load()
        self._calc_angles()
        
    def _load(self):
        path = self.path + self.lvl + '/'
        d = self.period[0]
        while d < self.period[1]:
            y = d.year
            m = d.month
            dd = d.day
            if not self.res:
                l = os.listdir(self.path + self.lvl + '/%.4i/'%y)
                l.sort()
                ll = ['solo_L2_mag-rtn-normal_%.4i%.2i%.2i'%(y,m,dd) in s for s in l]
                fname =''
                for i,n in enumerate(ll):
                    if n:
                        fname = self.path + self.lvl + '/%.4i/'%(y) + l[i]

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
        self.phi = np.arctan2(self.B_T,self.B_R)/np.pi*180
        self.theta = np.arctan2(self.B_N,np.sqrt(self.B_T**2+self.B_R**2))/np.pi*180
