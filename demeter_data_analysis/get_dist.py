import numpy as np
import matplotlib.pyplot as plt
import datetime as dt

from readDEMETER import find_files, get_data_from_parsed
from plots import plot_spectrogram, plot_TD, plot_map

import cartopy.crs as ccrs
import os
from os import listdir
from os.path import isfile, join


# gets maximum amplite from all data

def get_max(datapath,d_unit):

    file1 = open(datapath+"/edge_cases.txt","a")

    parsed_files = find_files(datapath)

    fcount = 0
    for fo in parsed_files:

        # next find the correct file names

        if fo[-4:] == 'DATp':
            print(fo)
            data = get_data_from_parsed(fo)
            
            #combine_packets(data)
            fs = 40000.
            tds = int(1e6/fs) # microseconds
            T_array = []
            E_array = []
            LTarray = []

            for i in range(int(len(data))):
                packetn = 'packet ' + str(i)
                time = data[packetn]['time']
                
                for m in range(int(len(data[packetn]['Edata']))):
                    T_array.append(time + dt.timedelta(microseconds=tds)) 
                    E_array.append(data[packetn]['Edata'][m])

                LTarray.append(data[packetn]["b'Local_Time "])

            T_array = np.array(T_array)
            E_array = np.array(E_array)
            LTarray = np.array(LTarray)

            max_list = []
            localtimelist = []
            
            savetimef = T_array[0].strftime('%Y_%m_%d_%H_%M')

            # chunk by every second ish 
            for p in range(0,len(E_array),32768):
                # max value of abs val of array in uV/m
                max_list.append(max(np.abs(E_array[p:p+32768])))
                #print(max(np.abs(E_array[p:p+32768])))
                if d_unit == 'mV/m':
                    if max(np.abs(E_array[p:p+32768])) > 10e3:
                        print(max(np.abs(E_array[p:p+32768])),'over')
                        file1.write(str(T_array[p])+'\n')
                else:
                    if max(np.abs(E_array[p:p+32768])) > 1:
                        print(max(np.abs(E_array[p:p+32768])),'over')
                        file1.write(str(T_array[p])+'\n')

            # save LT
            # save every 4th
            for k in range(0,len(LTarray),4):
                localtimelist.append(LTarray[k])

            # timing
            savetime = T_array[0].strftime('%Y_%m_%d_%H_%M')
            # save maximums
            np.savetxt(datapath+'/max_'+savetime+'.txt',np.array(max_list))
            np.savetxt(datapath+'/lt_'+savetime+'.txt',np.array(localtimelist))
        fcount+=1
        print('file ' + str(fcount) + ' of ' + str(len(parsed_files)))


def is_outlier(points, thresh=3.7):

    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh



def plot_dist(datapath,d_unit):

    rawfiles = [f for f in listdir(datapath) if isfile(join(datapath, f))]
    lt_day = []
    lt_night = []

    for r in rawfiles:
        if 'lt' in r and '.txt' in r:
            save_lt = r
            for r in rawfiles:
                if 'max' in r and save_lt[3:] in r and '.txt' in r:
                    save_max = r
            local_t = np.loadtxt(datapath+'/'+save_lt)
            if local_t[0] > 15:
                mm = np.loadtxt(datapath+'/'+save_max)
                lt_night.append(mm)

            else:
                mm = np.loadtxt(datapath+'/'+save_max)
                lt_day.append(mm)

    flat_day = [item for sublist in lt_day for item in sublist]
    flat_night = [item for sublist in lt_night for item in sublist]

    if d_unit == 'mV/m':
        flat_day = np.array(flat_day)/1e3
        flat_night = np.array(flat_night)/1e3
    else:
        flat_day = np.array(flat_day)
        flat_night = np.array(flat_night)

    filtered_day = flat_day[~is_outlier(flat_day)]
    filtered_night = flat_night[~is_outlier(flat_night)]

    import seaborn as sns

    sns.set_theme(style="darkgrid")

    ltimes = ['day', 'night']

    for ltime in ltimes:
        if ltime == 'day':
            filtered_data = filtered_day
            unfiltered_data = flat_day
            hist, bins, _ = plt.hist(unfiltered_data, bins=50)
            plt.close()
            
        else:
            filtered_data = filtered_night
            unfiltered_data = flat_night
            hist, bins, _ = plt.hist(unfiltered_data, bins=50)
            plt.close()

        fig, axs = plt.subplots(3, 1)
        axs[0].hist(unfiltered_data,bins=50)
        axs[0].set_title(ltime + ' max amplitude (unfiltered)')

        axs[1].hist(filtered_data,bins=50)
        axs[1].set_title(ltime + ' max amplitude (filtered)')

        # histogram on log scale. 
        # Use non-equal bin sizes, such that they look equal on log scale.
        logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
        axs[2].hist(unfiltered_data, bins=logbins)
        axs[2].set_xscale('log')
        axs[2].set_title(ltime + ' max amplitude (log)')
        
        fig.text(0.55, 0.02, d_unit, ha='center', va='center')
        fig.text(0.03, 0.5, 'counts', ha='center', va='center', rotation='vertical')
        fig.tight_layout()
        plt.savefig(datapath+'/' + ltime + '_distribution.png')
        plt.close()

        # print ltime
        overcount = 0
        for uf in unfiltered_data:
            if d_unit == 'mV/m':
                if uf > 20:
                    overcount+=1
            elif d_unit == 'nT':
                if uf > 10:
                    overcount+=1
        
        percent_over = overcount / len(unfiltered_data)
        print(ltime, 'percent over = ', percent_over)


def plot_edge_cases(datapath,d_unit):
    from readDEMETER import get_data_from_parsed
    from plots import plot_spectrogram

    f = open(datapath+'/'+'edge_cases.txt','r')
    edge_cases = []
    edge_case_start = 0
    last_dc = dt.datetime(2008,5,1,0,0,0)
    start_edge_case = last_dc

    alldates = []
    for line in f:
        #print(line)
        line = line.rstrip('\n')
        dc = dt.datetime.strptime(line,"%Y-%m-%d %H:%M:%S.%f")
        alldates.append(dc)
        time_change = dc - last_dc
        if np.abs(time_change.total_seconds()) < 3*60:
            if edge_case_start == 0:
                start_edge_case = last_dc
                edge_case_start = 1
                edge_cases.append(start_edge_case)
        elif np.abs(time_change.total_seconds()) > 120 and np.abs((start_edge_case - last_dc).total_seconds()) < 6*60:
            end_edge_case = last_dc
            edge_cases.append(end_edge_case)
        else:
            edge_case_start = 0
        last_dc = dc
    f.close()
    # also check if it is alone
    lone_edge_cases = []
    for tdate in alldates:
        hit = 0
        for pdate in alldates:
            checkt = pdate - tdate
            if np.abs(checkt.total_seconds()) > 1 and np.abs(checkt.total_seconds()) < 5*60: 
                hit = 1
        if hit == 0:
            lone_edge_cases.append(tdate)
    
    # final clean up
    edge_cases_clean = []
    for ei,edgec in enumerate(edge_cases[2:]):
        if ei == 0:
            continue
        if np.abs((edgec - edge_cases[ei-1+2]).total_seconds()) < 240:
            edge_cases_clean.append(edge_cases[ei-1+2])
            edge_cases_clean.append(edgec)
    
    rawfiles = [f for f in listdir(datapath) if isfile(join(datapath, f))]    
    for i in range(0,len(edge_cases_clean),2):
        current_case = edge_cases_clean[i]
        print(current_case)
        for r in rawfiles:
            if 'DATp' in r:
                rday = int(r[25:27])
                rhour = int(r[28:30])
                rmin = int(r[30:32])
                if rday == current_case.day and np.abs(rhour - current_case.hour) < 2 and np.abs(rmin - current_case.minute) < 5:
                    print(r)
                    data = get_data_from_parsed(datapath+'/'+r)
                       
                    st_date = current_case - dt.timedelta(seconds=10)
                    en_date = edge_cases_clean[i+1] + dt.timedelta(seconds=10)
                    print(st_date,en_date)

                    fig, axs = plt.subplots(2,sharex=True)
                    
                    plot_spectrogram(fig, axs[0],data, st_date, en_date, d_unit)
                    plot_TD(axs[1], data, st_date, en_date, d_unit)
                    plt.savefig(datapath +'/plots/' + dt.datetime.strftime(st_date, '%Y%m%d_%H%M') + '_plot')
    
    for i in range(0,len(lone_edge_cases)):
        current_case = lone_edge_cases[i]
        print(current_case)
        for r in rawfiles:
            if 'DATp' in r:
                rday = int(r[25:27])
                rhour = int(r[28:30])
                rmin = int(r[30:32])
                if rday == current_case.day and np.abs(rhour - current_case.hour) < 2 and np.abs(rmin - current_case.minute) < 5:
                    print(r)
                    data = get_data_from_parsed(datapath+'/'+r)
                        
                    st_date = current_case - dt.timedelta(seconds=10)
                    en_date = current_case + dt.timedelta(seconds=10)
                    print(st_date,en_date)

                    fig, axs = plt.subplots(2,sharex=True)
                    
                    plot_spectrogram(fig, axs[0],data, st_date, en_date, d_unit)
                    plot_TD(axs[1], data, st_date, en_date, d_unit)
                    
                    plt.savefig(datapath +'/plots/' + dt.datetime.strftime(st_date, '%Y%m%d_%H%M') + '_plot')


# whooo

datapath = '/media/rileyannereid/DEMETER/2008/01/1136'
d_unit = 'nT'
# find the maximums - processes all data, takes ~3hrs
#get_max(datapath,d_unit) 

# plot distribution and edge cases 
# distribuition is pretty quick, edge cases ~2hrs
#plot_dist(datapath,d_unit)
plot_edge_cases(datapath,d_unit)

# or here is code to plot one individually 
""""
r = 'DMT_N1_1136_206371_20080512_111032_20080512_111730.DATp'
current_case = dt.datetime(2008,5,12,11,10,32)

data = get_data_from_parsed(datapath+'/'+r)

st_date = current_case + dt.timedelta(seconds=0)
en_date = current_case + dt.timedelta(seconds=80)
print(st_date,en_date)

fig, axs = plt.subplots(2,sharex=True)

plot_spectrogram(fig, axs[0], data, st_date, en_date, d_unit)
plot_TD(axs[1], data, st_date, en_date, d_unit)

plt.savefig(datapath +'/plots/' + dt.datetime.strftime(st_date, '%Y%m%d_%H%M') + '_plot')
"""