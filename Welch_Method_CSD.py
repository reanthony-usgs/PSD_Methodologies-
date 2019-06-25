#!/usr/bin/env python

from obspy.core import read, UTCDateTime
from matplotlib.mlab import csd
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
import numpy as np
import matplotlib.pyplot as plt
from obspy.signal.invsim import evalresp
import glob
import sys

debug = False

net = "TA"

# get list of all stations in the network
stas = glob.glob("/msd/" + net + "*")
#print(stas)
stas = [sta.split("_")[1] for sta in stas]
print (stas)

stas = ['O18A']

#sys.exit()
print(stas)

loc = "--"



# Number of points of windows and overlap for PSD calculation 
# Easiest to split up in powers of two
nfft = 2**15
windlap = 0.5

# Set the start time - we do this since ASL uses Jdays 
stime = UTCDateTime("2008-227T12:00:00")
etime = stime + 3600
print(stime)

for sta in stas:

    # Set the path to the data
    path = "/home/reanthony/IRIS_mseed_Data/" + net + "_" + sta + "/" + str(stime.year) + "/" 
    path+=  str(stime.year) + "_" + str(stime.julday).zfill(3) + "_" + net + "_" + sta + "/" + "*.seed"
    
    if debug:
        print(path)
    
    st = read(path)
    
    # trim to the selected data
    st.trim(stime,etime)
    
    if debug:
        print(st)
    
    
    
    fig = plt.figure(1)
    
    for tr in st:
    
    
        power, freq = csd(tr.data, tr.data, NFFT=nfft,
                        noverlap=nfft*windlap, Fs = 1./tr.stats.delta, 
                        scale_by_freq=True)
                        
        if debug:
            print(power)
            print(freq)
        
        
        # Remove first Point     
        freq = freq [1:]
        power = power[1:] 
        
        # make power a real quantity 
        power = np.abs(power)
        
        # path to where the responses are
        resppath = "/home/reanthony/IRIS_mseed_Data/" + net + "_" + sta + "/" + str(stime.year) + "/" 
        resppath+= str(stime.year) + "_" + str(stime.julday).zfill(3) + "_" + net + "_" + sta + "/" + "RESP."
    
        # Add the usual station information 
        resppath += tr.id
        
        resp = evalresp(t_samp = tr.stats.delta, nfft=nfft, 
                filename=resppath, date = tr.stats.starttime,
                station=tr.stats.station, channel=tr.stats.channel,
                locid=tr.stats.location, network=tr.stats.network,
                units="ACC")
        
        
        # Plot the response 
        # Remove first value of Response as we did above 
        resp = resp[1:]  
        
        #fig = plt.figure(1)
        #plt.semilogx(1./freq, 20*np.log10(resp))
        #plt.xlabel('Period (s)')
        #plt.ylabel('Power (dB)')
        #plt.show()
        
        # Convert to dB and remove the Response 
        poweR = 10.*np.log10(power/(np.abs(resp)**2))
        powerNR = 10.*np.log10(((2*np.pi*freq)**2)*power/((1490.*(2.**26)/40.)**2))
        
        NLNMper,NLNMpower = get_nlnm()
        NHNMper,NHNMpower = get_nhnm()
        
        
        
        plt.semilogx(1./freq, poweR, label=tr.stats.station)
        #plt.semilogx(1./freq, powerNR, 'b', label='Sensitivity Only')
plt.semilogx(NLNMper, NLNMpower, linewidth=4., color='k')
plt.semilogx(NHNMper, NHNMpower, linewidth=4., color='k',label='NLNM/NHNM')
plt.xlabel('Period (s)')
plt.ylabel('Power (dB)')
plt.xlim((0.05,1000.))
plt.ylim((-200., -80.))
plt.legend()
plt.show()
        
import csv
f=open('Adam_O18A_PSD_5Seg.txt','w')

writer = csv.writer(f, delimiter='\t')
writer.writerows(zip(1./freq,poweR))

f.close()    
