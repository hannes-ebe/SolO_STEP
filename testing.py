###################################################################################################
### Testing and checking of MAG data and pitch angle calculation in space craft reference frame ###
###################################################################################################

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import load_nom_II as ld
import plot_nom_II as pt
from cdflib import CDF
from mag import MAGdata


period =(dt.datetime(2021,12,4,13,50),dt.datetime(2021,12,4,14,30))

dat = CDF(path='data/mag/2021/srf/solo_L2_mag-srf-normal_20211204_V01.cdf')
mag = dat.varget('B_SRF')
# lbl = dat.varget('LBL1_B_RTN')
# rep = dat.varget('REP1_B_RTN')
print(mag)
# print(lbl)
# print(rep)
print(dat.cdf_info())

