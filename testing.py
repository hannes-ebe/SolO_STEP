###################################################################################################
### Testing and checking of MAG data and pitch angle calculation in space craft reference frame ###
###################################################################################################

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import load_nom_II as ld
import plot_nom_II as pt
import time
from cdflib import CDF, cdfepoch
from mag import MAGdata

### EPOCH is given in ns since beginning of 2000 ###

period =(dt.datetime(2021,12,4,13,50),dt.datetime(2021,12,4,14,30))

dat_srf = CDF(path='data/mag/2021/srf/solo_L2_mag-srf-normal_20211204_V01.cdf')
time_srf = dat_srf.varget('EPOCH')
mag_srf = dat_srf.varget('B_SRF').T

t = cdfepoch.encode_tt2000(time_srf)
time_srf = []
for tt in t:
    time_srf.append(dt.datetime.fromisoformat(str(tt)[:26]))

dat_rtn = CDF('data/mag/2021/rtn/solo_L2_mag-rtn-normal-1-minute_20211204_V01.cdf')
time_rtn = dat_rtn.varget('EPOCH')
mag_rtn = dat_rtn.varget('B_RTN').T

t = cdfepoch.encode_tt2000(time_rtn)
time_rtn = []
for tt in t:
    time_rtn.append(dt.datetime.fromisoformat(str(tt)[:26]))

# print(mag_srf)
# print(mag_rtn)
# print(time_srf)
# print(time_rtn)

fig, (ax1,ax2) = plt.subplots(2,1,figsize=(10,6))

for mag_component in mag_srf:
    ax1.plot(time_srf,mag_component)

for mag_component in mag_rtn:
    ax2.plot(time_rtn,mag_component)

ax1.set_title('SRF')
ax2.set_title('RTN')

plt.savefig('pw_srf/comparison_mag_data.png')
plt.show()

# print(dat_srf.cdf_info())
# print(dat_rtn.cdf_info())

