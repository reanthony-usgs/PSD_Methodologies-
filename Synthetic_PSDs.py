#!/usr/bin/env python
from obspy.core import read
from obspy.signal.spectral_estimation import get_nhnm, get_nlnm
import matplotlib.pyplot as plt
from scipy.signal import periodogram, welch, gaussian, tukey, get_window
import numpy as np
import math
import matplotlib as mpl
import sys
#mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=16)

debug = True


# Window size defined in samples
def peterson_algorithm(st, debug = False):
    
    for st_window in st.slide(2**13, 2**13):
        f, p2 = welch(st_window[0].data, fs = 1., nperseg=int(np.floor(0.75*st_window[0].stats.npts/12.)),nfft=2**(math.ceil(np.log2(np.floor(1.75*st_window[0].stats.npts/12.))))) 
        # This is our 5 pt smoothing algorithm
        p2 = np.convolve(p2, np.ones((5,))/5)[(5-1):]
        if 'p' in vars():
            p = np.vstack((p, p2))
        else:
            p = p2
    if debug:
        print(p)
    p2 = np.percentile(p, 1., axis=0)
    p = np.median(p, axis=0)
    return f, p, p2

# Window size defined in seconds
def mcnamara_algorithm(data, debug = False):
    for st_window in st.slide(60*60, 60*30):
        NFFT = 2**(math.ceil(np.log2((float(st_window[0].stats.npts)/13.)))-1)
        
        if debug:
            print('Here is NFFT:' + str(NFFT))

        nlen = int(np.floor(.75*float(st_window[0].stats.npts)/13.))
        if debug:
            print(nlen)
        wind = get_window(('tukey',10), nlen)
        if debug:
            print(wind)
        f, p2 = welch(st_window[0].data, fs =1., window = wind , nfft=NFFT, nperseg = nlen)
        #f =0., p2 = 0.
        if 'p' in vars():
            p = np.vstack((p, p2))
        else:
            p = p2
    p2 = np.percentile(p, 20., axis=0)
    p = np.median(p, axis=0)
    
    return f, p, p2
 
def gausswin(N, alpha=2.5, debug = False): 
    if N == 1:
        return 1.
    n = np.arange(N) - float(N-1)/2. 
    if debug:
        print(n)


    gw = np.exp(-.5*((alpha*n)/(float(N-1)/2.))**2)
    if debug:
        print('Here is gausswin')
        print(gw)
    return gw



    
def gausssmooth(f, p, debug = False):
    np2, nf = [],[]
    
    # the calculation of setting the window width seems strange to me
    bw = 1./14.
    df = f[1]- f[0]
    for ind in range(len(f)):
        fc = f[ind]
        # Lower limit
        f1 = fc*(10**(bw*-1))
        # upper limit
        f2 = fc*(10**(bw))
        
        nf1 = int(np.floor(f1/f[0]))
        nf2 = int(np.ceil(f2/f[0]))
        #nf1 = int(np.floor(len(f)*f1/f[-1]))
        #nf2 = int(np.ceil(len(f)*f2/f[-1]))
        
        
        # make sure we stay within frequency bounds of FFT
        nf1 = max([1, nf1])
        nf2 = min([len(f)-1, nf2-1])
        
        if debug:
            print('Here is nf1 and nf2: ' + str(nf1) + ' ' + str(nf2))
            print('Here is center frequency:' + str(fc))
            print('Here are are freq limits:' + str(f[nf1]) + ' ' + str(f[nf2]))
        
        sp = p[nf1:nf2]
        gw = gausswin(len(sp))
        #print(gw)
        #gw = gaussian(len(sp), 1.)
        #if debug:
            #print('Here is gw')
            #print(gw)
        np2.append(np.dot(sp,gw)/np.sum(gw))
        nf.append(fc)

    np2 = np.asarray(np2)
    nf = np.asarray(nf)

# This is likely incorrect, frequency should not change    
    return f, np2    
    
# Use 2 Hour Windows with 50% overlap 
def berger_algorithm(data, debug = True):
    for st_window in st.slide(2*60*60, 2*30*60):
        st2 = st_window.copy()
        st2[0].detrend('linear')
        st2[0].detrend('constant')
        NFFT = int(2**(math.ceil(np.log2(st_window[0].stats.npts))))
        # Check periodogram estimate
        f, p = periodogram(st_window[0].data, fs=1., window='hann',  nfft=NFFT)
        # print(p)
        f = f[1:]
        p=p[1:]
        f, p = gausssmooth(f,p)

        if 'pB' in vars():
            pB = np.vstack((pB, p))
        else:
            pB = p
    pB2 = np.percentile(pB, 1., axis=0)     
    print(pB2)  
    pB = np.median(pB, axis = 0)
    
    return f, pB, pB2


def octavesmooth(x, octaveFraction = 1.0/3.0 ):
    y = []
    rightFactor = math.pow(2.0, octaveFraction / 2.0 )
    leftFactor =  1.0 / rightFactor
    for n in range(len(x)):
        left = long(n * leftFactor)
        right = long(n * rightFactor)
        y.append( sum( x[left : right + 1] ) / (right - left + 1) )
		
    return y





st = read('XX.XXXX.00.BNZ.msd')
#st[0].data = np.append(st[0].data, st[0].data) 
st[0].data=np.asarray(st[0].data, dtype=np.float64)/(2.**32)
replace =[]

f = open('Test_synth_020sps_17280000.txt','r')
for line in f:
    replace.append(float(line))
# We're replacing all the data in the minised file    
st[0].data = np.asarray(replace)
st[0].stats.sampling_rate = 20

# What is decimated.........
st[0].decimate(5)
st[0].decimate(2)
st[0].decimate(2)
f,p = periodogram(st[0].data, window='hann', fs=1.)
#po = octavesmooth(p, 1.)







    #plt.semilogx(1./f, 10.*np.log10(po), label='Octave Smoothing', color='C4')
fig = plt.figure(2, figsize=(12,12))
ax1 = plt.subplot(2,1,1)
ax2=plt.subplot(2,1,2)
for idx in range(2):
    plt.subplot(2,1,idx+1)
    # plot Peterson Models
    per_nlnm, pow_nlnm = get_nlnm()
    plt.semilogx(per_nlnm,pow_nlnm, linewidth=2, color='k')
    per_nhnm, pow_nhnm = get_nhnm()
    plt.semilogx(per_nhnm,pow_nhnm, linewidth=2, color='k', label='NLNM/NHNM')
    
    # straight periodogram of full 24 hours
    f,p = periodogram(st[0].data, window='hann', fs=1.)
    plt.semilogx(1./f, 10.*np.log10(p), label='Periodogram', alpha=.3)
    
    
    # Plot the Berger Model
    f = open('GSN_noisemodel.txt','r')
    freqs, amps =[], []
    for line in f:
        line = ' '.join(line.split())
        line = line.split(' ')
        freqs.append(1./float(line[0]))
        amps.append(float(line[1]))
    
    freqs = np.asarray(freqs)
    plt.semilogx(1./freqs, amps, label='GSN Vertical Noise Model', color='.5', linewidth=2)
    f.close()



# Calculate PSDs using Peterson method
f,p, p2 = peterson_algorithm(st)

f2=open('Peterson_Statistics.txt','w')

for triple in zip(f,p,p2):
    f2.write(str(triple[0])+ ", " + str(triple[1]) + ", " + str(triple[2]) + '\n')
  
f2.close  

plt.subplot(2,1,1)
plt.semilogx(1./f, 10.*np.log10(p), label='Peterson (1993)', color='C1')
plt.subplot(2,1,2)
plt.semilogx(1./f, 10.*np.log10(p2), label='Peterson (1993) Minimum', color='C1')


# Calculate PSDs using Berger Method without Gaussian Smoother
f, p, p2 = berger_algorithm(st)


f3=open('Berger_Statistics.txt','w')

for triple in zip(f,p,p2):
    f3.write(str(triple[0])+ ", " + str(triple[1]) + ", " + str(triple[2]) + '\n')
  
f3.close  
        
if debug:
    print('Here are the Berger models')
    print(f)
    print(p)
    print(p2)
#sys.exit()

# Lets write the Berger model out to a flat file




plt.subplot(2,1,1)
plt.semilogx(1./f, 10.*np.log10(p), label='Berger et al. (2004)', color='C2')
plt.subplot(2,1,2)
plt.semilogx(1./f, 10.*np.log10(p2), label='Berger et al. (2004) Minimum', color='C2')


# McNamara and Buland Method 
f, p, p2 = mcnamara_algorithm(st)
plt.subplot(2,1,1)
plt.semilogx(1./f, 10.*np.log10(p), label='McNamara and Buland (2004)', color='C3')

plt.subplot(2,1,2)
plt.semilogx(1./f, 10.*np.log10(p2), label='McNamara and Buland (2004)', color='C3')
po = octavesmooth(p, 1.)
plt.subplot(2,1,1)

# McNamara and Buland after 1 octave smoothing 
plt.semilogx(1./f, 10.*np.log10(po), label='Octave Smoothing', color='C4')
plt.subplot(2,1,2)
po2 = octavesmooth(p2, 1.)
plt.semilogx(1./f, 10.*np.log10(po2), label='Octave Smoothing', color='C4')


for idx in range(2):
    plt.subplot(2,1,idx+1)
    #plt.legend(loc=1)
    plt.ylim((-200., -120.))
    plt.xlabel('Period (s)')
    plt.ylabel('Power (dB rel. 1 $(m/s^2)^2/Hz$)')
    plt.xlim((2., 2000))

ax1.text(-0.08, 1., '(a)', transform=ax1.transAxes,
      fontsize=26, fontweight='bold', va='top', ha='right')
ax2.text(-0.08, 1., '(b)', transform=ax2.transAxes,
      fontsize=26, fontweight='bold', va='top', ha='right')

handles, labels = ax1.get_legend_handles_labels()
ax = plt.gca()
leg = fig.legend(handles, labels, loc = 'lower center', ncol = 3, fontsize = 15)
plt.subplots_adjust(bottom = 0.14)


plt.show()
#plt.savefig('SyntheticDataexample.jpg', format='JPEG', dpi=400)
#plt.savefig('SyntheticDataexample.pdf', format='PDF', dpi=400)


