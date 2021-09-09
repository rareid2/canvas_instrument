import numpy as np
import matplotlib.pyplot as plt
import datetime as dt

from readDEMETER import find_files, get_data_from_parsed
from plots import plot_spectrogram, plot_TD, plot_map

import cartopy.crs as ccrs

# set data path to hard drive
# need to do a file tree walk
datapath = '/home/rileyannereid/workspace/canvas_demeter/data/all_data'
parsed_files = find_files(datapath)
dunits = ['mV/m','nT']

for fo in parsed_files:
#    if '20070415_04' in fo and '1131' in fo:
#        if '1131' in fo:
#            d_unit = dunits[0]
#        else:
#            d_unit = dunits[1]

        data = get_data_from_parsed(fo)
        st_date = dt.datetime(2007,4,15,5,0)
        en_date = st_date + dt.timedelta(minutes=5)

        fig = plt.figure(figsize=(6.25, 8))
        E_TD = plt.subplot(1, 1, 1)

        #projection=ccrs.PlateCarree()
        #E_map = plt.subplot(2, 1, 2, projection=projection)

        #plot_map(E_map, projection, data, st_date, en_date)
        plot_TD(E_TD, data, st_date, en_date,d_unit)

        fig.suptitle('DEMETER Burst Data at ' + dt.datetime.strftime(st_date, '%Y-%m-%d %H:%M:%S'))
        plt.savefig('/home/rileyannereid/workspace/canvas_demeter/data/bursts_pngs/' + dt.datetime.strftime(st_date, '%Y%m%d_%H%M') + 'burstdata')
        plt.show()
        plt.close()
        #plot_spectrogram(data, st_date, en_date, d_unit)

