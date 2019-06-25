""" synthnoise.py
Create synthetic time-series from noise models.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from obspy.core import *


class Noise():
  """ Implementation of different earth noise models
    Noise(type='NLMN',units='ACC')
    Currently defaults are the only ones implemented
    NLNM - New low-noise model (Peterson 1993)  
    
    Noise at a particular period or a an array of periods is determined by the function
    Noise.calculate(period) 
    This then applies the requested noise model by mapping private functions to the calculate
    function.
    The function calculate returns either a single amplitude or an array of amplitudes
    depending on input.
  """
  _nlnm_coeff=np.array([[.1,-162.36,5.64],
    [.17,-166.7,0.00],
    [.4,-170.,-8.3],
    [.8,-166.4,28.9],
    [1.24,-168.60,52.48],
    [2.4,-159.98,29.81],
    [4.3,-141.1,0.],
    [5.,-71.36,-99.77],
    [6.,-97.26,-66.49],
    [10.,-132.18,-31.57],
    [12.,-205.27,36.16],
    [15.6,-37.65,-104.33],
    [21.9,-114.37,-47.10],
    [31.6,-160.58,-16.28],
    [45.,-187.5,0.],
    [70.,-216.47,15.70],
    [101.,-185.,0.],
    [154.,-168.34,-7.61],
    [328.,-217.43,11.9],
    [600.,-258.28,26.6],
    [10000.,-346.88,48.75],
    [100000.,np.nan,np.nan]])

  _nhnm_coeff=np.array([[0.1,-108.73,-17.23],
    [.22,-150.34,-80.5],
    [.32,-122.31,-23.87],
    [.8,-116.85,32.51],
    [3.8,-108.48,18.08],
    [4.6,-74.66,-32.95],
    [6.3,0.66,-127.18],
    [7.9,-93.37,-22.42],
    [15.4,73.54,-162.98],
    [20.,-151.52,10.01],
    [354.8,-206.66,31.63],
    [100000.,np.nan,np.nan]])
    
  def __init__(self,type='NLNM',units='ACC'):
    if type=='NLNM':
      self.units='ACC'
      self.type=type
      #self.calculate=_perterson
  
  def calculate(self,p):
    """ Implement it the slow way at the moment this could be linearized, but for now it 
    will work.  
    """
    p=np.asarray(p,dtype=np.float64)
    out=np.zeros(np.shape(p))
    p_out=np.zeros(np.shape(p))
    for j, period in np.ndenumerate(p):

      for c in self._nlnm_coeff:
        if period >= c[0]:
          if self.units=='ACC':
            out[j]=c[1]+c[2]*np.log10(period)
            p_out[j]=period
    return p_out,out
  
class SyntheticTS():
  """ Create a synthetic time series given an appropriate PSD of a noise model"""
  def __init__(self,psd=[],period=[],nseg=10,sps=20,durseg=2*3600,sigma_whitenoise=1.):
    self.nseg=nseg
    self.sps=sps
    self.durseg=durseg
    self.sigma_whitenoise=sigma_whitenoise
    self.N=sps*durseg
    self.psd=np.asarray(psd)  #PSD to create synthetic for
    self.period=np.asarray(period) #Periods matching PSD values
    
  def get_synthetic_frequencies(self):
    return np.fft.rfftfreq(self.sps*self.durseg,d=1./self.sps)   
    
  def create(self,period=[],psd=[]):
    """Create a synthetic time series from a PSD with dB (units**2/Hz).  It is up to the
    user to maintain proper units.
    """
    if len(self.psd)==0 and len(psd)>0:
      self.psd=np.asarray(psd)
    if len(self.period)==0 and len(period)>0:
      self.period=np.asarray(period)
    # Convert from PSD to spectral amplitude
    # Note no conversion for /Hz unit??????
    noise_specamp=np.sqrt(np.power(10,self.psd/10.))
    noise_specamp[0]=0. #Force DC component to zero acceleration may be better to go really small    

    ts_seg=self.sigma_whitenoise*np.random.randn(self.nseg,self.N)
    for i in np.arange(0,self.nseg):
      init_fft=np.fft.rfft(ts_seg[i,:],norm='ortho')#Apply the real fft to the randomly generated data
      # Scale our FFT appropriately
      # In this case we scale our amplitude by 2 twice once for each RFFT and inverse??????
      ts_seg[i,:]=np.fft.irfft(noise_specamp*init_fft*2.,norm='ortho')*(2.)

    full_ts=ts_seg[0,:]
    # print(np.shape(full_ts))
    for i in np.arange(1,self.nseg):
       full_ts=np.append(full_ts,ts_seg[i,:])
    self.ts=full_ts    
  
  def savetxt(self,fpath): 
    """ wrapper for numpy savetxt where fpath is the file path and name to save to"""
    np.savetxt(fpath,self.ts,fmt='%g')
    
  def to_obspy(self,net='XX',sta='XXXX',loc='00',chan='BNZ',starttime='',countspervolt=2**32):
    tr=Trace()
    tr.stats.sampling_rate=self.sps
    tr.stats.delta=1./self.sps
    tr.stats.npts=len(self.ts)
    tr.stats.network=net
    tr.stats.location=loc
    tr.stats.station=sta
    tr.stats.channel=chan
    if starttime=='':
      starttime=UTCDateTime.now()
    else: 
      starttime=UTCDateTime(starttime)
    tr.stats.starttime=starttime
    tr.data=np.asarray(self.ts*countspervolt,dtype=np.int32)
    return tr
    
  def plot(self,type='PSD'):
    """ type can be PSD or TS for timeseries or type ALL for all options"""
    if type=='PSD' or type=='ALL':
      f,P=scipy.signal.welch(self.ts,nperseg=len(self.ts), 
        window=scipy.signal.get_window('boxcar',len(self.ts)),scaling='density',fs=self.sps)
      per=1./(f)
      per[0]=0.

      plt.figure()
      plt.semilogx(per,decibels(P,1.),label='Welch',alpha=.5) #*len(self.ts)
      plt.semilogx(self.period,self.psd,color='k',label='NLNM',linewidth=2)
      plt.legend()
      plt.xlabel('Period (s)')
      plt.ylabel('PSD Acceleration (dB rel. 1 (m/s**2)**2/Hz)')
    if type=='TS' or type=='ALL':
      plt.figure()
      plt.plot()
    plt.show()        

def decibels(x,ref):
  return 10.*np.log10(x/ref)
  


# Create our synthetic Timeseries object
synth=SyntheticTS(sps=20,durseg=7200,nseg=100,sigma_whitenoise=1.)
# Create our Noise model object
nlnm=Noise()
# Determine the frequencies we will be using in our fft's
per_fft=1./synth.get_synthetic_frequencies()
# Calculate our noise model at the correct frequencies
period,psd=nlnm.calculate(per_fft)
# Create our synthetic by convolving white noise with our psd we want
synth.create(period=period,psd=psd)
synth.savetxt("Test2_synth_%03dsps_%d.txt" %(synth.sps,len(synth.ts)))
synth.plot(type='PSD')
tr=synth.to_obspy()
tr.write("%s.msd" % (tr.id),format='MSEED')
