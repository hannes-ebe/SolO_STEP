'''Erste Analyse von Elektronenevents in STEP-Daten'''

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt

from analyze_nom_II import STEP

#dat = STEP(2021,12,4,rpath='data/')
dat = STEP(2021,12,4)

box_list = [[[15,35],[30,38]],[[20,45],[25,30]],[[26,80],[20,25]],[[30,120],[15,20]],[[40,155],[10,15]],[[50,170],[3,10]]]
pixel_list = [4,9,14]

# dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/')
# dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/',period=[dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)],fit=True)
# dat.fit_energy(pixel=4,norm='tmax',save='analyze/',period=[dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)],fit=True,p0=[17500,1,-10,0])

# dat.pixel_integral_window(filename='test.png',norm='tmax',save='etracks/',period=[dt.datetime(2021,12,4,14,30),dt.datetime(2021,12,4,14,35)])
# dat.evolution_energy_means(filename='test_evolution.png',norm='tmax',save='etracks/',period=[dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)],box_list=box_list)
# dat.plot_ts(period=[dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)], head=0, save='etracks/', norm='tmax',box_list=box_list)

dat.evolution_energy_means_ts(filename='test_evolution.png',head=0,norm='tmax',save='etracks/',period=[dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)],box_list=box_list,pixel_list=pixel_list)


# dat = STEP(2021,12,5)
# dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/')

# dat = STEP(2022,1,14)
# dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/')

# dat = STEP(2022,1,18)
# dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/')

# Da Funktioniert etwas nicht...
# dat = STEP(2021,2,7)
# dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/')