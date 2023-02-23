''' Testing and working with STEP data'''

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import load_nom_II as ld
import plot_nom_II as pt

# Loading data
time, data = ld.load_nom(period=(dt.datetime(2022,2,27),dt.datetime(2022,2,28)), products=('M','A'))
print('Data loaded successfully.')

# Combining data (Main and Auxiliary Product)
ld.combine_data(time, data)
print('Data combined successfully.')

# Plotting data (Aufruf zum Speichern in Plotfunktion)
pt.plot_ts(data, time, head=0, save='plots/')
pt.plot_ts(data, time, head=1, save='plots/')
pt.plot_ts_diff(data,time, save='plots/', single=False)
print('Data plotted successfully.')