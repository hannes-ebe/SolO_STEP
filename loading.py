'''Loading of STEP data in extra file'''

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import load_nom_II as ld
import plot_nom_II as pt


# Loading data
time, data = ld.load_nom(period=(dt.datetime(2022,2,27),dt.datetime(2022,2,28)), products=('M','A'))

# Combining data (Main and Auxiliary Product)
ld.combine_data(time, data)
