import numpy as np
import pandas as pd

from spacepy import pycdf
from datetime import datetime, timedelta
import bisect
from scipy import signal
from scipy import interpolate
import spacepy.datamodel as dm

from scipy.fftpack import fft
from scipy import stats

import matplotlib.pyplot as plt
import matplotlib.colors as colors
plt.style.use('seaborn')
import seaborn as sns
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 200

# ----------------------------------------------------------

def getCdfData(cdf, start, stop, epoch_str, data_str):
    '''Finds indicies for given time interval to grab corresponding cdf data'''
    
    # determine indicies for start and end times
    start_ind = bisect.bisect_left(cdf[epoch_str], start)
    stop_ind = bisect.bisect_right(cdf[epoch_str], stop)
    
    # grab data for specified time window only
    time = cdf[epoch_str][start_ind:stop_ind]
    data = cdf[data_str][start_ind:stop_ind] #,:]
    
    # return epoch and data files
    return(time, data)

# ----------------------------------------------------------

def adcUnitConversion(data, data_tag):
    '''input array of potential or search coil measurements and convert from ADC units to mV and nT'''
    
    n_bits = 16
    
    if data_tag == 'vb2':
        meas_range = 12.5e3 #mV
        ad_res = 0.4 #mV
        
    elif data_tag == 'mscb2':
        meas_range = 12 #nT #had as 11 for some reason?
        ad_res = 0.36/1e3 #nT
    
    converted_data = data*(2*meas_range/(2**n_bits))
    return(converted_data)

# ----------------------------------------------------------

def nearest_ind(items, pivot):
    '''finds data index nearest to pivot'''
    
    time_diff = np.abs([date - pivot for date in items])
    return int(time_diff.argmin(0))

# ----------------------------------------------------------

def newBurstFiles(vb2_fn, mscb2_fn, burst_start, burst_stop, file_time): 
    '''given data files and time interval saves one burst of data in a new file with units converted'''

    ## SPECIFY CDF FILES ##
    vb2_cdf = pycdf.CDF(vb2_fn)
    mscb2_cdf = pycdf.CDF(mscb2_fn)

    ## GRAB DATA FOR BURST WINDOW ##
    vb2_time, vb2_data = getCdfData(vb2_cdf, burst_start, burst_stop, 'epoch', 'vb2')
    mscb2_time, mscb2_data = getCdfData(mscb2_cdf, burst_start, burst_stop, 'epoch', 'mscb2')

    ## CONVERT UNITS ##
    vb2_converted = adcUnitConversion(vb2_data, 'vb2')
    mscb2_converted = adcUnitConversion(mscb2_data, 'mscb2')

    V1 = vb2_converted[:,0]               
    V2 = vb2_converted[:,1]
    V3 = vb2_converted[:,2]
    V4 = vb2_converted[:,3]
    V5 = vb2_converted[:,4]
    V6 = vb2_converted[:,5]
    Bu = mscb2_converted[:,0]
    Bv = mscb2_converted[:,1]
    Bw = mscb2_converted[:,2]

    ## SAVE DATA ##
    data = dm.SpaceData()

    data['Vt'] = vb2_time
    data['V1'] = V1
    data['V2'] = V2
    data['V3'] = V3
    data['V4'] = V4
    data['V5'] = V5
    data['V6'] = V6

    data['Bt'] = mscb2_time
    data['Bu'] = Bu
    data['Bv'] = Bv
    data['Bw'] = Bw

    save_fn = 'b2_data_'+file_time
    dm.toJSONheadedASCII(save_fn, data)

# ----------------------------------------------------------

def str2datetime(time_array):
    '''converts various arrays of time strings into datetime objects'''
    
    for i in range(0, len(time_array)):
        
        if len(time_array[i]) == len('2012-01-22T00:14:24') and 'T' in time_array[i]:
            time_array[i] = (datetime.strptime(time_array[i][:], "%Y-%m-%dT%H:%M:%S"))
            
        elif len(time_array[i]) == len('2012-01-22T00:14:24') and 'T' not in time_array[i]:
            time_array[i] = (datetime.strptime(time_array[i][:], "%Y-%m-%d %H:%M:%S"))
            
        elif len(time_array[i]) != len('2012-01-22T00:14:24') and 'T' not in time_array[i]:
            time_array[i] = (datetime.strptime(time_array[i][:], "%Y-%m-%d %H:%M:%S.%f"))
            
        else:
            time_array[i] = (datetime.strptime(time_array[i][:], "%Y-%m-%dT%H:%M:%S.%f"))

    return(time_array)

# ----------------------------------------------------------


# ADAPTED FROM DAVID MALASPINA
def rotation_matrices(B_in):
    '''creates rotation matrices from magnetometer date to rotate to a coordinate system with z-axis aligned with the B-field'''
    
    tlen = np.shape(B_in)[0]
    
    M = np.zeros((tlen, 3, 3))
    
    for i in range(0, tlen):

        Bx = B_in[i,0]
        By = B_in[i,1]
        Bz = B_in[i,2]

        # ******************************
        # Determine angles for rotation matricies
        phi = np.arctan(np.abs(Bx)/np.abs(By))
        theta = np.arctan(np.sqrt(Bx**2 + By**2) / np.abs(Bz))
        # this is arctan2! remove abs() when switching to arctan2 - and then I dont need quantrant corrections?
        # np.sqrt(Bx**2 + By**2) --> np.norm2(Bx, By)

        # ******************************
        # Add in quantrant specific changes to get signs right

        if Bx > 0 and By > 0: 
            phi = phi
        if Bx > 0 and By < 0: 
            phi = np.pi - phi
        if Bx < 0 and By > 0: 
            phi = 2*np.pi - phi 
        if Bx < 0 and By < 0: 
            phi = phi + np.pi

        if Bz > 0: 
            theta = theta
        if Bz < 0: 
            theta = np.pi - theta


        # ******************************
        # Build Phi rotation array

        M2 = np.zeros((3,3))

        M2[0,0] = np.cos(phi)
        M2[1,0] = -np.sin(phi)
        M2[2,0] = 0

        M2[0,1] = np.sin(phi)
        M2[1,1] = np.cos(phi)
        M2[2,1] = 0

        M2[0,2] = 0
        M2[1,2] = 0
        M2[2,2] = 1

        # ******************************
        # Build theta rotation array

        M3 = np.zeros((3,3))

        M3[0,0] = 1
        M3[1,0] = 0
        M3[2,0] = 0

        M3[0,1] = 0
        M3[1,1] = np.cos(theta)
        M3[2,1] = -np.sin(theta)

        M3[0,2] = 0
        M3[1,2] = np.sin(theta)
        M3[2,2] = np.cos(theta)

        # ******************************
        # Combine array transformations
        # numpy or scipy has rotation matrices! rotx, roty, rotz... can replace all the M2 and M3 stuff 
        M4 = np.matmul(M2,M3)
        M[i,:,:] = M4

    return (M)

# ----------------------------------------------------------

def intepolateM(t_efw, t_emfisis, M):
    # need to rename inputs
# t_efw =  Bt from mscb2   
# t_emfisis = "Epoch"
    '''interpolates rotation matrices from magnetometer time axis to search coil time axis'''
    
    # create array to save matrices into
    M_interp = np.zeros((len(t_efw), 3, 3))

    for row in range(0,3):
        for col in range(0,3):

            x = np.asarray([t.timestamp() for t in t_emfisis])
            y = M[:, row, col]
            tck = interpolate.splrep(x, y) #, s=0)

            xnew = [t.timestamp() for t in t_efw]
            shift = x[0] - xnew[0] 
            xnew += shift
            ynew = interpolate.splev(xnew, tck) #, der=0)
            xnew_dt = np.asarray([datetime.fromtimestamp(t) for t in xnew])

            M_interp[:, row, col] = ynew
            
    return M_interp, x, xnew

# ----------------------------------------------------------

def newAlignedFiles(M_interp_B, M_interp_E, Bu, Bv, Bw, Eu, Ev, Ew):
    '''rotates B and E data into field aligned coordinate system using rotation matrices'''

    B_aligned = np.zeros((np.shape(M_interp_B)[0], 3))
    E_aligned = np.zeros((np.shape(M_interp_E)[0], 3))

    for i in range(0,len(Bu)):
        Bx = Bu[i]
        By = Bv[i]
        Bz = Bw[i]
        M_B = M_interp_B[i,:,:]
        
        vin = [Bx, By, Bz]
        B_out = np.matmul(vin,M_B)
        B_aligned[i,:] = B_out[:]
    
    for i in range(0,len(Eu)):
        Ex = Eu[i]
        Ey = Ev[i]
        Ez = Ew[i]
        M_E = M_interp_E[i,:,:]
        
        vin = [Ex, Ey, Ez]
        E_out = np.matmul(vin,M_E)
        E_aligned[i,:] = E_out[:]
        
    return B_aligned, E_aligned

# ----------------------------------------------------------

def get_fft(sig, fs):
    '''scipy.signal.spectrogram for complex and PSD output using hann window'''
    # see if mode=['complex', 'psd'] gives both outputs

    f, t, Sxx = signal.spectrogram(sig, fs, window=('hann'),
                                   nperseg=512, mode='complex')
    
    f_psd, t_psd, PSD = signal.spectrogram(sig, fs, window=('hann'),
                                               nperseg=512, mode='psd')
    
    
    return(f, t, Sxx, PSD)


# ----------------------------------------------------------

def bin_centers(bin_edges, mean):
    # use aritmetic if bins are linearly spaced
    # use geometric if bins are log-spaced
    
    if mean is 'arithmetic':
        bin_centers = (bin_edges[1:] + bin_edges[:-1])/2
        
    elif mean is 'geometric':
        bin_centers = np.sqrt(bin_edges[1:] * bin_edges[:-1])
        
    # length check
    if len(bin_edges) == len(bin_centers)+1:
        return(bin_centers)

# ----------------------------------------------------------

def histFit(PSD):

    hist, bin_edges = np.histogram(np.log10(PSD), bins=20, range=None, normed=None, weights=None, density=1)
    bin_cent = bin_centers(bin_edges, mean='arithmetic') #check that this is actually right

    # Gaussian Fit
    mu, sigma = stats.norm.fit(np.log10(PSD))
    median = np.median(np.log10(PSD))

    threshold_2sigma = mu + (2*sigma)
    threshold_1sigma = mu + (1*sigma)
    threshold_8median = np.log10(8*(10**median))

    fit_stats = [mu, sigma, median, threshold_2sigma, threshold_1sigma, threshold_8median]

    # Skewnorm Fit
    a, loc, scale = stats.skewnorm.fit(np.log10(PSD))
    mean, var, skew, kurt = stats.skewnorm.stats(a, loc=loc, scale=scale, moments='mvsk')

    fit_stats = [a, loc, scale, mean, var, skew, kurt]
        
    return(fit_stats)
    
# ----------------------------------------------------------

def str2complex(array):
    '''converts str to np.complex'''
    new_array = array
    rows = np.shape(array)[0]
    cols = np.shape(array)[1]
    
    for i in range(0, rows):
        for j in range(0, cols):
            if type(array[i,j]) is str:
                new_array[i,j] = np.complex(array[i,j])
    
    return new_array

# ----------------------------------------------------------

def santolik_B(freqs, S_Bx, S_By, S_Bz, fs, tidx):
    ''' The Santolik method, for time-domain signals E and B '''
    
    FBx = S_Bx[:,tidx]
    FBy = S_By[:,tidx]
    FBz = S_Bz[:,tidx]

    # output space
    # wave vector stuff:
    theta_vec = np.zeros_like(freqs)
    phi_vec = np.zeros_like(freqs)
    k1_vec = np.zeros_like(freqs)
    k2_vec = np.zeros_like(freqs)
    k3_vec = np.zeros_like(freqs)

    # polarization stuff:
    F_vec = np.zeros_like(freqs)
    Lp_vec = np.zeros_like(freqs)

    # power?
    Smag_vec = np.zeros_like(freqs, dtype=complex)
    
    # The algorithm!
    for fi in np.arange(1,len(freqs)):
#         print(fi)
        
        # B components
        B1 = FBx[fi]
        B2 = FBy[fi]
        B3 = FBz[fi]
        
        # S components - simplify to np.outer(z,np.conj(z))
        # i = 1, j = 1,2,3
        S11 = B1*np.conj(B1)
        S12 = B1*np.conj(B2)
        S13 = B1*np.conj(B3)
        # i = 2, j = 1,2,3
        S21 = B2*np.conj(B1)
        S22 = B2*np.conj(B2)
        S23 = B2*np.conj(B3)
        # i = 3, j = 1,2,3
        S31 = B3*np.conj(B1)
        S32 = B3*np.conj(B2)
        S33 = B3*np.conj(B3)
        
        Smag = S11 + S22 + S33
        Smag_vec[fi] = Smag
        
        # Build A -- try something like np.concatenate([np.real(A),np.imag(A)],axis=0)
        A = [[np.real(S11),  np.real(S12),  np.real(S13)],
             [np.real(S12),  np.real(S22),  np.real(S23)],
             [np.real(S13),  np.real(S23),  np.real(S33)],
             [0,            -np.imag(S12), -np.imag(S13)],
             [np.imag(S12),  0,            -np.imag(S23)],
             [np.imag(S13),  np.imag(S23),  0           ]]

        A = np.nan_to_num(A) # check why I would be expecting a nan - need to add a check if its something that shouldnt happen
        # SVD
        U,s,VT=np.linalg.svd(A,full_matrices=False)
        
        # define singular vales in ascending order (and corresponding rows of VT)
        w1 = s[2] 
        w2 = s[1]
        w3 = s[0]
        
        VT1 = VT[2,:]
        VT2 = VT[1,:]
        VT3 = VT[0,:]

        # solve for polarization parameters
        if w3!=0:
            F = 1 - np.sqrt(w1/w3)
            Lp = w2/w3
        if w3==0:
            F = np.nan
            Lp = np.nan
        
        # save in output vectors
        F_vec[fi] = F
        Lp_vec[fi] = Lp
        
        # define unit wave vector components
        k1 = VT1[0]
        k2 = VT1[1]
        k3 = VT1[2]
            
        # solve for wave vector angles -- might be able to use arctan2? check that np.abs should be there
        theta = np.abs(np.arctan(np.sqrt(k1**2 + k2**2)/k3))
        if k1>=0:
            phi = np.arctan(k2/k1)
        if k1<0 and k2<0:
            phi = np.arctan(k2/k1) - np.pi
        if k1<0 and k2>=0:
            phi = np.arctan(k2/k1) + np.pi
        
        # save in output vectors
        theta_vec[fi] = theta
        phi_vec[fi] = phi
        k1_vec[fi] = k1
        k2_vec[fi] = k2
        k3_vec[fi] = k3
    
        

    return freqs, theta_vec, phi_vec, F_vec, Lp_vec, Smag_vec, k1_vec, k2_vec, k3_vec

    