''' Testing and working with STEP data'''

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import load_nom_II as ld
import plot_nom_II as pt

# Loading data
# Muss noch gemacht werden

# Plotting data
pt.plot_ts(data, time, head=0)
pt.plot_ts(data, time, head=1)
pt.plot_ts_diff(data,time)