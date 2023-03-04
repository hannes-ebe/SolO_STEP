'''Erste Analyse von Elektronenevents in STEP-Daten'''

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt

from analyze_nom_II import STEP

dat = STEP(2021,12,4)
dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/')

dat = STEP(2021,12,5)
dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/')

dat = STEP(2022,1,14)
dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/')

dat = STEP(2022,1,18)
dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/')

dat = STEP(2021,2,7)
dat.marginal_distribution(pixel=4,norm='tmax',save='analyze/')