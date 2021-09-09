# this script reads the demeter data files
import os
import datetime as dt
from os import listdir
from os.path import isfile, join
import numpy as np

def get_data_from_parsed(filename): # send in full path to parsed data file (from IDL output), return dictionary with packets

    # open the file, in read binary mode
    f = open(filename, 'rb')
    datalines = [line for line in f]

    # data we might want to keep
    inds = np.arange(15, 33, 1) # block 2 info
    time_ind = 6 # block 4 info

    # create a dictionary
    data = {}

    # first clean up the lines
    for p in range(0, len(datalines), 1429):
        n = int(p//1429)
        packetn = 'packet ' + str(n)
        data[packetn] = {}
        for ind in inds:
            line = datalines[ind+p]
            key_str = ""
            line = str(line)
            for ci, char in enumerate(line):
                if char == '=':
                    # everything after here is the float we want
                    key_val = line[ci+1:-3]
                    if ind == 30 or ind == 32:
                        kv = key_val.strip()
                        kv = kv.split()
                        kv = np.array((float(kv[0]), float(kv[1]), float(kv[2])))
                    else:
                        kv = key_val.strip()
                        kv = float(kv)
                    break
                    
                key_str += str(char)
            data[packetn][key_str] = kv
            
        # finally, grab the time
        time = datalines[6+p].strip()
        time = time.split()
        # note that dateime objs take in microseconds as last arg. 
        time_obj = dt.datetime(int(time[1]), int(time[2]), int(time[3]), int(time[4]), int(time[5]), int(time[6]), int(time[7])*1000)  
        data[packetn]['time'] = time_obj
        data[packetn]['fs'] = 40e3 # Hz 
        data[packetn]['unit'] = 'uV/m'

        # now we've got all the header info we want, let's get the actual data
        packet_data = datalines[p+63:p+1429]
        if packet_data == []: # werid fix
            continue
        packet_data[0] = packet_data[0][8:] # remove the data intro

        packet_strip = [pd.strip() for pd in packet_data]
        packet_split = [pd.split() for pd in packet_strip]
        
        # add in final data cleaned up
        pdata = []
        for ps in packet_split:
            for pd in ps:
                #print(pd)
                pdata.append(float(pd))

        data[packetn]['Edata'] = np.array(pdata) # add data to dict. for this packet
    print('found ' + str(n+1) + ' packets')
    f.close()
    return data

def find_files(datapath):
    rawfiles = [f for f in listdir(datapath) if isfile(join(datapath, f))]
    rawfiles_final = []
    for r in rawfiles:
        if 'DATp' in r:
            pass 
        elif r + 'p' in rawfiles:
            pass
        else:
            rawfiles_final.append(r)

    if len(rawfiles_final) != 0:
        for f in rawfiles_final:

            fn = datapath + '/' + f
            if '.txt' in fn or '.png' in fn:
                continue
            cwd = os.getcwd()
            os.chdir('/home/rileyannereid/workspace/canvas_demeter/IDL_edited')
            os.system("idl -e 'rd_dmt_n1' -arg " + fn)
            os.chdir(cwd)

    parsefiles = [f for f in listdir(datapath) if isfile(join(datapath, f))]
    final_parsefiles = []
    for p in parsefiles:
        if 'DATp' in p:
            p = datapath + '/' + p
            final_parsefiles.append(p)
    return final_parsefiles