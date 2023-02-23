''' Plotting of head0, head1 and diff.'''

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import load_nom_II as ld
import plot_nom_II as pt


def plot(year,month,day,lastofmonth=False):
    # Loading data
    if lastofmonth:
        if month!=12:
            time, data = ld.load_nom(period=(dt.datetime(year,month,day),dt.datetime(year,month+1,1)), products=('M','A'))
        else:
            time, data = ld.load_nom(period=(dt.datetime(year,month,day),dt.datetime(year,1,1)), products=('M','A'))
    else:
        time, data = ld.load_nom(period=(dt.datetime(year,month,day),dt.datetime(year,month,day+1)), products=('M','A'))
    print('Data loaded successfully.')
    
    # Combining data (Main and Auxiliary Product)
    ld.combine_data(time, data)
    print('Data combined successfully.')
    
    # Plotting data (Aufruf zum Speichern in Plotfunktion)
    pt.plot_ts(data, time, head=0, save='plots/')
    pt.plot_ts(data, time, head=1, save='plots/')
    pt.plot_ts_diff(data,time, save='plots/', single=False)
    print('Data plotted successfully.')
    
    
# Oktober 2021 
# for i in range(22,32):
#     mon = 10
#     y = 2021
#     if i != 31:
#         plot(y,mon,i)
#     else:
#         plot(y,mon,i,lastofmonth=True)
        
# November 2021
for i in range(1,31):
    mon = 11
    y = 2021
    if i != 30:
        plot(y,mon,i)
    else:
        plot(y,mon,i,lastofmonth=True)

# Dezember 2021
for i in range(1,32):
    mon = 12
    y = 2021
    if i != 31:
        plot(y,mon,i)
    else:
        plot(y,mon,i,lastofmonth=True)

# Januar 2022
for i in range(1,32):
    mon = 1
    y = 2022
    if i != 31:
        plot(y,mon,i)
    else:
        plot(y,mon,i,lastofmonth=True)

# Februar 2022
for i in range(1,29):
    mon = 2
    y = 2022
    if i != 28:
        plot(y,mon,i)
    else:
        plot(y,mon,i,lastofmonth=True)
        
# März 2022
for i in range(1,9):
    mon = 3
    y =2022
    if i != 31:
        plot(y,mon,i)
    else:
        plot(y,mon,i,lastofmonth=True)
