'''Workbench zum Herumprobieren'''

import matplotlib as mpl
import pylab as plt
import numpy as np
import datetime as dt
import matplotlib.dates as mdates
import load_nom_II as ld
import plot_nom_II as pt

from scipy.optimize import curve_fit

print([i for i in range(1,16)])

# print(str(dt.datetime(2023,3,30,8,30)))

# save = 'etracks/'
# filename = 'test.png'

# fig = plt.figure(figsize = (15,10))
# fig.subplots_adjust(wspace = 0, hspace = 0)
# ax = []
# for i in range(16):
#     if i == 0:
#         ax.append(fig.add_subplot(4,5,3))
#     else:
#         ax.append(fig.add_subplot(4,5,5+i,sharex = ax[0], sharey = ax[0]))
# if save:
#     if type(save) == str:
#         plt.savefig(save + filename)
#     else:
#         plt.savefig(filename)
# print('Plotted ' + filename + ' successfully!')