'''Erste Analyse von Elektronenevents in STEP-Daten'''

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt

from analyze_nom_II import STEP

dat = STEP(2021,12,4)
# dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/')
# dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/',period=[dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)],fit=True)
# dat.fit_energy(pixel=4,norm='tmax',save='analyze/',period=[dt.datetime(2021,12,4,13,30),dt.datetime(2021,12,4,16,30)],fit=True,p0=[17500,1,-10,0])
dat.pixel_integral_window(filename='test.png',norm='tmax',save='etracks/',period=[dt.datetime(2021,12,4,14,30),dt.datetime(2021,12,4,14,35)])



# dat = STEP(2021,12,5)
# dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/')

# dat = STEP(2022,1,14)
# dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/')

# dat = STEP(2022,1,18)
# dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/')

# Da Funktioniert etwas nicht...
# dat = STEP(2021,2,7)
# dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/')