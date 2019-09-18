#!/usr/bin/env python


# This script calculates 30 s PSDs from the RidgeCrest Aftershock using 10
# s data windows. Plots up the median and 95% Observation interval statistics


from obspy.core import Stream, read, UTCDateTime
import matplotlib.pyplot as plt
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
import numpy as np
from matplotlib.mlab import csd
from obspy.signal.invsim import evalresp
from math import pi
import sys
from itertools import combinations
import math, copy
import random
import csv


debug = True

net = "GS"
sta = "CA04"
loc = "00"
chan = "HHZ"

# Set the start day - we will set various start dates in the window below
sday = UTCDateTime("2019-193T00:00:00")


# Set the number of windows you want to compute PSDs over 
Nwin = 200

# Set the path to the data
path = "/msd/" + net + "_" + sta + "/" + str(sday.year) + "/" + str(sday.julday).zfill(3) 
path += "/" + loc + "_" + chan + "*.seed"
       


# Path to Response 

# path to where the responses are

resppath = "/APPS/metadata/RESPS/RESP." + net + "." + sta + "." + loc + "." + chan

if debug:
    print(resppath) 


# Number of points of windows and overlap for PSD calculation 
# Easiest to split up in powers of two
nfft = 2**10
windlap = 0.5
    

st = read(path)
st.detrend('constant')
st.merge(fill_value=0.)

# Start the loop of PSD calculations 

for i in range(0, Nwin):
    
    st2 = st.copy()
    
    # get a random start time in seconds
    # Note this will always be an integer second
    
    ss = random.randint(0,86349)
    
    if "H_day" not in vars():
        H_day = ss/(60.*60.)

    else:
        H_day = np.vstack((H_day,ss/(60.*60.)))
    
    stime = sday + ss
    
    # Trim a 20 s window
    etime = stime + 30
    
    st2.trim(stime,etime)
    
    # Now get the PSD
    print(i)
    tr = st2[0]
    power, freq = csd(tr.data, tr.data, NFFT=nfft,
                    noverlap=nfft*windlap, Fs = 1./tr.stats.delta, 
                    scale_by_freq=True)
    
    freq = freq [1:]
    power = power[1:] 
    
    power = np.abs(power)
    
    # Get the Response
    resp = evalresp(t_samp = tr.stats.delta, nfft=nfft, 
                filename=resppath, date = tr.stats.starttime,
                station=tr.stats.station, channel=tr.stats.channel,
                locid=tr.stats.location, network=tr.stats.network,
                units="ACC")
    resp = resp[1:]
    
    # Convert to dB and remove the response 
    
    power = 10.*np.log10(power/(np.abs(resp)**2))
    
    # Assemble the Matrix of PSDs
    
    if "Power_M" not in vars():
        Power_M = power
        
    else:
        Power_M = np.vstack((Power_M,power))
        
# get the mean and STD

Power_Avg = np.mean(Power_M,axis=0)
Power_STD = np.std(Power_M,axis=0)

print(Power_Avg.shape)        


# Plot up the mean PSD and 95% confidence intervals (2*STD)

NLNMper,NLNMpower = get_nlnm()
NHNMper,NHNMpower = get_nhnm()

freq = np.matrix.transpose(freq)


print(type(H_day[0]))
print(type(freq[0]))

f=open('RidgeCrest_HHZ_freqs.txt','w')

for fw in freq:
    f.write(str(fw) + '\n')
f.close()

H_day = np.asarray(H_day)

f2=open('RidgeCrest_Hours.txt','w')

for HD in H_day:
    f2.write(str(float(HD)) + '\n')
    
f2.close()



f3=open('RidgeCrest_PSDs.txt','w')

writer = csv.writer(f3, delimiter='\t')
writer.writerows(Power_M)
        
        
fig = plt.figure(1)

plt.semilogx(1./freq, Power_Avg, linewidth=5., color='k', label='Mean PSD')
plt.semilogx(1./freq, Power_Avg - (2*Power_STD), linewidth=2., color='r')
plt.semilogx(1./freq, Power_Avg + (2*Power_STD), linewidth=2., color='r', label='95th percentile')
plt.semilogx(NLNMper, NLNMpower, linewidth=2., color='k')
plt.semilogx(NHNMper, NHNMpower, linewidth=2., color='k',label='NLNM/NHNM')
plt.xlabel('Period (s)')
plt.ylabel('Power (dB)')
plt.xlim((0.05,1000.))
plt.ylim((-200., -80.))
plt.legend()
plt.show()        
        
        
