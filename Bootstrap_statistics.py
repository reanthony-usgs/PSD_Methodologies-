#!/usr/bin/env python


# This script loads in a day-long PSD and calculates the peridogram
# Estimate for several randomly selected windows. 

# Confidence intervals are then assigned based on the distribution of 
# these peridogram estimates in various bins.



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


net = "IU"
sta = "ANMO"
loc = "00"
chan = "LHZ"

# Set the start day - we will set various start dates in the window below
sday = UTCDateTime("2018-200T00:00:00")


# Set the number of windows you want to compute PSDs over 
Nwin = 100000


# read in several months of data

st = Stream()
Daynum = 0.

for day in range(200,366):
    st += read("/msd/" + net + "_" + sta + "/" + str(sday.year) + "/" + str(day).zfill(3) + "/" + loc + "_" + chan + "*.seed")
    Daynum=Daynum+1.
    
# Add in 2019 data    
    
for day in range(1,131):
    st += read("/msd/" + net + "_" + sta + "/2019/" + str(day).zfill(3) + "/" + loc + "_" + chan + "*.seed")
    Daynum=Daynum+1.
    

print(Daynum)

# Path to Response 

# path to where the responses are

resppath = "/APPS/metadata/RESPS/RESP." + net + "." + sta + "." + loc + "." + chan
print(resppath) 


# Number of points of windows and overlap for PSD calculation 
# Easiest to split up in powers of two
nfft = 2**10
windlap = 0.5
    
total_seconds = ((60.*60.*24.*Daynum)-(3601.*3.))

st.detrend('constant')
st.merge(fill_value=0.)



# Start the loop of PSD calculations 


for i in range(0, Nwin):
    
    st2 = st.copy()
    
    # get a random start time in seconds
    
    # Note this will always be an integer second
    
    
    
    ss = random.randint(0,total_seconds)
    
    
    if "H_day" not in vars():
        H_day = ss/(60.*60.*24)
        
    else:
        H_day = np.vstack((H_day,ss/(60.*60.*24)))
    
    
    
    stime = sday + ss
    
    # Trim a 3-hour window
    etime = stime + (3600.*3.)
    
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


# Write to file
# np.savez('ANMO_2018_309_LHZ.npz',Power_M,freq)

#np.savez('Freqs_BHZ.npz', freq, fmt='%.4f')


freq = np.matrix.transpose(freq)


print(type(H_day[0]))
print(type(freq[0]))

#f=open('LHZ_freqs.txt','w')

#for fw in freq:
#    f.write(str(fw) + '\n')
    
#f.close()

H_day = np.asarray(H_day)

f2=open('Hours_3Hr_ANMO_2019_STS6_3.txt','w')

for HD in H_day:
    f2.write(str(float(HD)) + '\n')
    
f2.close()



f3=open('ANMO_3Hr_2019_STS6_LHZ_3.txt','w')

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
        
        
    
    
                
    
    
    
    
                        
    
    
    
    
    
        
    
    






    



