import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import os
import datetime
import scipy.signal
import scipy.linalg
import matplotlib.gridspec as gridspec
from scipy import interpolate
import bisect
from scipy.fftpack import fft
from scipy.signal.windows import hann
from scipy.signal import blackman

# canvas science packets
# APID 0 spectra (real)
# APID 1 x-spectra (real)
# APID 2 x-spectra (imag)


def get_data(packet_file, APID):
    # parses the packets as they come in
    
    if APID == 0:
        # parse into array of 5 values
        # with 67 arrays (one for each fbin) 
        packet_data = np.array([1,2,3,4,5])

    return packet_data


def get_Ez(packets):
    # takes all packets for one fbin
    # resolve the 6th channel
    return

def 

    




def santolik(Q):

    A = []
    B = []
    for l in range(6):
        A_tmp = np.zeros([3,3],'complex')
        B_tmp = np.zeros([3],'complex')

        # Valid entries in Levi-Civita symbol (minus 1, because python is zero-indexed)
        for e,i,j,k in ([1,0,1,2],[1,1,2,0],[1,2,0,1],[-1,0,2,1],[-1,1,0,2],[-1,2,1,0]):
            A_tmp[i,j] = e*Q[k+3,l]
            B_tmp[i]= Q[i,l]
        A.extend(A_tmp)
        B.extend(B_tmp)


    A2 = np.concatenate([np.real(A),np.imag(A)],axis=0)
    B2 = np.concatenate([np.real(B),np.imag(B)],axis=0)

    U,s,VT=np.linalg.svd(A2,full_matrices=False)
    V = VT.T
    W_inv = np.linalg.inv(np.diag(s))
    n = np.linalg.multi_dot([V,W_inv,U.T,B2])
    


        # angles
        theta_calc = np.arctan(np.sqrt(n[0]**2 + n[1]**2)/n[2])
#         print("theta is:",theta*R2D,"deg")
        phi_calc   = np.arctan(n[1]/n[0])
        if (n[0]<0 and n[2] < 0):
            phi_calc -= np.pi
        if (n[0]<0 and n[2] >= 0):
            phi_calc += np.pi
            
            
        # wrap phi
        if phi_calc > np.pi:
            phi_calc -= 2.*np.pi
        if phi_calc <= -np.pi:
            phi_calc+= 2.*np.pi

        theta_vec[fi] = theta_calc
        phi_vec[fi] = phi_calc
        n_vec[fi,:] = n
        
        # Planarity:
        beta = A2.dot(n) 
        big_N = np.sum(pow(beta - B2,2))
        big_D = np.sum(pow(np.abs(beta) + np.abs(B2),2))
        planarity[fi] = 1 - np.sqrt(big_N/big_D)
        
        # Uh... spectral power or something? I'm making this up
        Qmag[fi] = np.sum(np.trace(Q))
        
        # Poyinting flux
        poynting_flux[fi] = np.linalg.norm(np.cross(np.array([FBx[fi], FBy[fi], FBz[fi]],'complex'),
                        np.array([FEx[fi], FEy[fi], FEz[fi]],'complex'))/mu)
        
        # Poynting flux again, but computed from the elements of Q rather than E(f) and B(f):
        pflux_from_Q[fi] = np.sqrt(
                            np.abs(Q[1,3] + Q[4,0])**2 +
                            np.abs(Q[0,5] + Q[3,2])**2 +
                            np.abs(Q[2,4] + Q[5,1])**2)/c/mu


    return freqs, theta_vec, phi_vec, n_vec, planarity, Qmag, poynting_flux, pflux_from_Q