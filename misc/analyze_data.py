'''Erste Analyse von Elektronenevents in STEP-Daten'''

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt

from analyze_nom_II import STEP

# dat = STEP(2021,12,4,rpath='data/')
# dat = STEP(2021,12,4)

# box_list = [[[15,35],[30,38]],[[20,45],[25,30]],[[26,80],[20,25]],[[30,120],[15,20]],[[40,155],[10,15]],[[50,170],[3,10]]]
# box_list_2 = [[[15,35],[30,35]],[[20,45],[25,30]],[[30,60],[20,25]],[[30,80],[15,20]],[[50,100],[10,15]],[[100,150],[3,10]]]
pixel_list = [i for i in range(1,16)] # [4,9,14] 

# dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/')
# dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/',period=[dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)],fit=True)
# dat.fit_energy(pixel=4,norm='tmax',save='analyze/',period=[dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)],fit=True,p0=[17500,1,-10,0])

# dat.pixel_integral_window(filename='test.png',norm='tmax',save='etracks/',period=[dt.datetime(2021,12,4,14,30),dt.datetime(2021,12,4,14,35)])
# dat.evolution_energy_means(filename='test_evolution.png',norm='tmax',save='etracks/',period=[dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)],box_list=box_list)
# dat.plot_ts(period=[dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)], head=0, save='etracks/', norm='tmax',box_list=box_list)

# dat.evolution_energy_means_ts(filename='test_evolution.png',head=0,norm='tmax',save='etracks/',period=[dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)],box_list=box_list,pixel_list=pixel_list,norm_pixel=3,close=True)
# dat.evolution_energy_means_ts(filename='test_evolution_2.png',head=0,norm='tmax',save='etracks/v2',period=[dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)],box_list=box_list_2,pixel_list=pixel_list,norm_pixel=4,close=True)

# dat.distribution_ring(filename='Test', title='Mean of energy sorted by energy of pixel 3', head=0, norm='tmax', save='etracks/', period=[dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)],box_list=box_list, norm_pixel = 3, res = '1min', overflow = True, esquare = False,window_width = 5, close=True, sorted_by_energy=True)
# dat.distribution_ring(filename='Test_time', title='Mean of energy (time series)', head=0, norm='tmax', save='etracks/', period=[dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)],box_list=box_list, norm_pixel = 3, res = '1min', overflow = True, esquare = False,window_width = 5, close=True, sorted_by_energy=False)

# dat.wrapper_distribution_ring('21_12_04', head=0, norm='tmax', period=[dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)], box_list=box_list, norm_pixel=3)

# dat = STEP(2021,12,5)
# dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/')

# dat = STEP(2022,1,14)
# dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/')

# dat = STEP(2022,1,18)
# dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/')

# Da Funktioniert etwas nicht...
# dat = STEP(2021,2,7)
# dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/')



'''Untersuchung von Events zu denen Magnetfeldaten vorliegen'''

'''Ereignisse vor Mitte April 2022'''

dat = STEP(2022,2,7)

box_list=[[[35,65],[12,23]],[[40,85],[8,12]],[[55,100],[2,8]]]
period=[dt.datetime(2022,2,7,19,30),dt.datetime(2022,2,7,22)]

dat.evolution_energy_means_ts(filename='22_02_07.png',head=-1,norm='tmax',save='etracks_test/',period=period,box_list=box_list,pixel_list = pixel_list,norm_pixel=3)
dat.wrapper_distribution_ring('22_02_07_electrons', head=-1, norm='tmax', period=period, box_list=box_list, norm_pixel=3)


# dat = STEP(2021,12,4)

# box_list = [[[15,35],[30,38]],[[20,45],[25,30]],[[26,80],[20,25]],[[30,120],[15,20]],[[40,155],[10,15]],[[50,170],[3,10]]]
# period = [dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)]

# dat.evolution_energy_means_ts(filename='21_12_04.png',head=-1,norm='tmax',save='etracks_test/',period=period,box_list=box_list,pixel_list = pixel_list,norm_pixel=3)
# dat.wrapper_distribution_ring('21_12_04_electrons', head=-1, norm='tmax', period=period, box_list=box_list, norm_pixel=3)


# dat = STEP(2022,11,19)

# box_list=[[[25,50],[28,35]],[[30,60],[23,28]],[[35,70],[18,23]],[[40,110],[12,18]],[[60,165],[2,12]]]
# period = [dt.datetime(2022,11,19,13,30),dt.datetime(2022,11,19,17,30)]

# dat.evolution_energy_means_ts(filename='22_11_19.png',head=-1,norm='tmax',save='etracks_test/',period=period,box_list=box_list,pixel_list = pixel_list,norm_pixel=3)
# dat.wrapper_distribution_ring('22_11_19_electrons', head=-1, norm='tmax', period=period, box_list=box_list, norm_pixel=3)


# dat = STEP(2022,12,1)

# box_list=[[[58,67],[35,40]],[[60,70],[32,35]],[[65,75],[28,32]],[[70,82],[23,28]],[[74,90],[18,23]],[[74,100],[14,18]],[[77,132],[2,14]]]
# period = [dt.datetime(2022,12,1,6,30),dt.datetime(2022,12,1,9,30)]

# dat.evolution_energy_means_ts(filename='22_12_01.png',head=-1,norm='tmax',save='etracks_test/',period=period,box_list=box_list,pixel_list = pixel_list,norm_pixel=3)
# dat.wrapper_distribution_ring('22_12_01_electrons', head=-1, norm='tmax', period=period, box_list=box_list, norm_pixel=3)


# dat = STEP(2022,12,12)

# box_list=[[[25,50],[23,30]],[[30,70],[17,23]],[[40,80],[12,17]],[[50,85],[8,12]],[[60,90],[2,8]]]
# period = [dt.datetime(2022,12,12,7,30),dt.datetime(2022,12,12,10)]
        
# dat.evolution_energy_means_ts(filename='22_12_12.png',head=-1,norm='tmax',save='etracks_test/',period=period,box_list=box_list,pixel_list = pixel_list,norm_pixel=3)
# dat.wrapper_distribution_ring('22_12_12_electrons', head=-1, norm='tmax', period=period, box_list=box_list, norm_pixel=3)


# dat = STEP(2022,12,24)

# box_list=[[[20,35],[35,40]],[[25,40],[30,35]],[[28,50],[23,30]],[[33,58],[17,23]],[[40,64],[12,17]],[[45,75],[8,12]],[[50,80],[2,8]]]
# period = [dt.datetime(2022,12,24,4),dt.datetime(2022,12,24,6)]

# dat.evolution_energy_means_ts(filename='22_12_24.png',head=-1,norm='tmax',save='etracks_test/',period=period,box_list=box_list,pixel_list = pixel_list,norm_pixel=3)
# dat.wrapper_distribution_ring('22_12_24_electrons', head=-1, norm='tmax', period=period, box_list=box_list, norm_pixel=3)














# dat1 = STEP(2021,12,4)
# dat2 = STEP(2022,11,19)
# dat3 = STEP(2022,12,1)
# dat4 = STEP(2022,12,12)
# dat5 = STEP(2022,12,24)


# # Before Mid-April 2022:

# # 2021-12-04:
# box_list = [[[15,35],[30,38]],[[20,45],[25,30]],[[26,80],[20,25]],[[30,120],[15,20]],[[40,155],[10,15]],[[50,170],[3,10]]]
# period = [dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)]

# dat1.evolution_energy_means_ts(filename='21_12_04.png',head=0,norm='tmax',save='etracks/',period=period,box_list=box_list,pixel_list = pixel_list,norm_pixel=3)
# dat1.wrapper_distribution_ring('21_12_04', head=0, norm='tmax', period=period, box_list=box_list, norm_pixel=3)


# # After Mid-April 2022:

# # 2022-11-19:
# box_list=[[[25,50],[28,35]],[[30,60],[23,28]],[[35,70],[18,23]],[[40,110],[12,18]],[[60,165],[2,12]]]
# period = [dt.datetime(2022,11,19,13,30),dt.datetime(2022,11,19,17,30)]

# dat2.evolution_energy_means_ts(filename='22_11_19.png',head=0,norm='tmax',save='etracks/',period=period,box_list=box_list,pixel_list = pixel_list,norm_pixel=3)
# dat2.wrapper_distribution_ring('22_11_19', head=0, norm='tmax', period=period, box_list=box_list, norm_pixel=3)

# # 2022-12-01:
# box_list=[[[58,67],[35,40]],[[60,70],[32,35]],[[65,75],[28,32]],[[70,82],[23,28]],[[74,90],[18,23]],[[74,100],[14,18]],[[77,132],[2,14]]]
# period = [dt.datetime(2022,12,1,6,30),dt.datetime(2022,12,1,9,30)]

# dat3.evolution_energy_means_ts(filename='22_12_01.png',head=0,norm='tmax',save='etracks/',period=period,box_list=box_list,pixel_list = pixel_list,norm_pixel=3)
# dat3.wrapper_distribution_ring('22_12_01', head=0, norm='tmax', period=period, box_list=box_list, norm_pixel=3)

# # 2022-12-12:
# box_list=[[[25,50],[23,30]],[[30,70],[17,23]],[[40,80],[12,17]],[[50,85],[8,12]],[[60,90],[2,8]]]
# period = [dt.datetime(2022,12,12,7,30),dt.datetime(2022,12,12,10)]
        
# dat4.evolution_energy_means_ts(filename='22_12_12.png',head=0,norm='tmax',save='etracks/',period=period,box_list=box_list,pixel_list = pixel_list,norm_pixel=3)
# dat4.wrapper_distribution_ring('22_12_12', head=0, norm='tmax', period=period, box_list=box_list, norm_pixel=3)

# # 2022-12-24:
# box_list=[[[20,35],[35,40]],[[25,40],[30,35]],[[28,50],[23,30]],[[33,58],[17,23]],[[40,64],[12,17]],[[45,75],[8,12]],[[50,80],[2,8]]]
# period = [dt.datetime(2022,12,24,4),dt.datetime(2022,12,24,6)]

# dat5.evolution_energy_means_ts(filename='22_12_24.png',head=0,norm='tmax',save='etracks/',period=period,box_list=box_list,pixel_list = pixel_list,norm_pixel=3)
# dat5.wrapper_distribution_ring('22_12_24', head=0, norm='tmax', period=period, box_list=box_list, norm_pixel=3)






