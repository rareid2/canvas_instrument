import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
import datetime as dt

import cartopy.crs as ccrs
import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.feature.nightshade import Nightshade

# add back in the convertsion of u to m in the TD 

def plot_spectrogram(fig, ax, data, st_date, en_date, d_unit):

    #fig = plt.figure()

    #combine_packets(data)
    fs = 40000.
    tds = int(1e6/fs) # microseconds
    T_array = []
    E_array = []

    for i in range(int(len(data))):
        packetn = 'packet ' + str(i)
        time = data[packetn]['time']
        
        for m in range(int(len(data[packetn]['Edata']))):
            T_array.append(time + dt.timedelta(microseconds=tds)) 
            E_array.append(data[packetn]['Edata'][m])

    T_array = np.array(T_array)
    T_stamp = [(t - dt.datetime(1970,1,1)) for t in T_array]
    T_stamp = np.array(T_stamp)

    st_date1 = (st_date - dt.datetime(1970, 1, 1))
    en_date1 = (en_date - dt.datetime(1970, 1, 1))

    st_idx = (np.abs(T_stamp - st_date1)).argmin()
    en_idx = (np.abs(T_stamp - en_date1)).argmin()
    T_array = np.array(T_array[st_idx:en_idx+1])
    E_array = np.array(E_array[st_idx:en_idx+1])

    overlap = 0.5
    nfft = 8192
    window = 'hanning'

    ff, tt, FE = scipy.signal.spectrogram(E_array, fs=fs, window=window,
                                        nperseg=nfft, mode='psd', scaling='density')
    E_mag = np.log10(FE)
    if d_unit == 'mV/m':
        pcm = ax.pcolormesh(tt, ff/1e3, E_mag, cmap='jet',vmin=-2, vmax=2)
    else:
        pcm = ax.pcolormesh(tt, ff/1e3, E_mag, cmap='jet',vmin=-10, vmax=-3)
    
    #cbar = fig.colorbar(pcm, ax=ax)
    #cbar.set_label('log[uV^2/Hz]')
    current_lt = data['packet 0']["b'Local_Time "]
    ax.set_title('Spectrogram Local Time = '+str(current_lt))
    #ax.set_xlabel('seconds after ' + str(T_array[0]))
    ax.set_ylabel('Frequency [kHz]')

def plot_TD(maxes, data, st_date, en_date, d_unit):

    #combine_packets(data)
    fs = 40000.
    tds = int(1e6/fs) # microseconds
    T_array = []
    E_array = []

    for i in range(int(len(data))):
        packetn = 'packet ' + str(i)
        time = data[packetn]['time']
        
        for m in range(int(len(data[packetn]['Edata']))):
            T_array.append(time + dt.timedelta(microseconds=tds)) 
            E_array.append(data[packetn]['Edata'][m])

    T_array = np.array(T_array)
    T_stamp = [(t - dt.datetime(1970,1,1)) for t in T_array]
    T_stamp = np.array(T_stamp)

    st_date1 = (st_date - dt.datetime(1970, 1, 1))
    en_date1 = (en_date - dt.datetime(1970, 1, 1))

    st_idx = (np.abs(T_stamp - st_date1)).argmin()
    en_idx = (np.abs(T_stamp - en_date1)).argmin()
    T_array = np.array(T_array[st_idx:en_idx+1])
    E_array = np.array(E_array[st_idx:en_idx+1])

    newt = []
    for t in T_array:
        sec_after = (t - T_array[0]).total_seconds()
        newt.append(sec_after)

    #E_array = np.abs(E_array)
    #E_avg = np.mean(E_array)
    #Esort = E_array.tolist()
    #Esort.sort()
    
    if d_unit == 'mV/m':
        #print(max(E_array))
        maxes.plot(newt, E_array/1e3,'-') 
    else:
        maxes.plot(newt, E_array,'-')
        #print(T_array[0],E_avg,Esort[-1])  
    maxes.set_ylabel(d_unit)
    maxes.set_xlabel('seconds after ' + str(T_array[0]))
    maxes.set_title('time domain data')

def plot_map(ax, projection, data, st_date, en_date):
    ax.coastlines()
    ax.set_xticks(np.linspace(-180, 180, 5), crs=projection)
    ax.set_yticks(np.linspace(-90, 90, 5), crs=projection)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.stock_img()
    ax.gridlines()


    lats = []
    lons = []
    times = []

    for i in range(int(len(data))):
        packetn = 'packet ' + str(i)
        lats.append(data[packetn]["b'Geocentric_lat "])
        lons.append(data[packetn]["b'Geocentric_long "])
        times.append(data[packetn]["time"])
    #print(lats)
    T_array = np.array(times)
    T_stamp = [(t - dt.datetime(1970,1,1)) for t in T_array]
    T_stamp = np.array(T_stamp)

    st_date1 = (st_date - dt.datetime(1970, 1, 1))
    en_date1 = (en_date - dt.datetime(1970, 1, 1))

    st_idx = (np.abs(T_stamp - st_date1)).argmin()
    en_idx = (np.abs(T_stamp - en_date1)).argmin()
    T_array = np.array(T_array[st_idx:en_idx+1])
    lats = np.array(lats[st_idx:en_idx+1])
    lons = np.array(lons[st_idx:en_idx+1])
    nlons = []
    for ll in lons:
        if ll > 180:
            dll = ll-180
            newll = -180 + dll 
        else:
            newll = ll
        nlons.append(newll)

    ax.plot(nlons, lats, c='r', linewidth=3)
    ax.set_ylim([-90,90])
    ax.set_xlim([-180,180])


    ax.add_feature(Nightshade(times[0]))
