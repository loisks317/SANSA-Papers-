#SDalgorithm.py
#
# take the median of all SD values between range gates 25-40 at each time
# compare this to EISCAT values at all heights
# take the closest value
# this is the 'height' the SD is seeing
#
# LKS, SANSA, Feb 2016. Listening to Selena Gomez and questioning life choices.
# but it's so bad that I can focus. Seriously, google 'Slow Down', it's terrible
#
# love my comments, Feb 2016, Svalbard (bitches)
#
import numpy as np
import matplotlib as mpl
import os
import readProcessedEISCAT
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt # best package ever
import matplotlib.dates as dt
from matplotlib.dates import HourLocator, DayLocator, DateFormatter, MinuteLocator
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter,FuncFormatter
import readProcessedEISCAT
from numpy import ma
import pandas as pd
#
# functions
def nearestDate(dates, pivot):
  return dates.index(min(dates, key=lambda x: abs(x - pivot)))

def runningMedian(seq, M):
    """
     Purpose: Find the median for the points in a sliding window (odd number in size) 
              as it is moved from left to right by one point at a time.
      Inputs:
            seq -- list containing items for which a running median (in a sliding window) 
                   is to be calculated
              M -- number of items in window (window size) -- must be an integer > 1
      Otputs:
         medians -- list of medians with size N - M + 1
       Note:
         1. The median of a finite list of numbers is the "center" value when this list
            is sorted in ascending order. 
         2. If M is an even number the two elements in the window that
            are close to the center are averaged to give the median (this
            is not by definition)
    """
    from itertools import islice
    from collections import deque,Counter
    from bisect import insort, bisect_left

    seq = iter(seq)
    s = []   
    m = M // 2

    # Set up list s (to be sorted) and load deque with first window of seq
    s = [item for item in islice(seq,M)]    
    d = deque(s)

    # Simple lambda function to handle even/odd window sizes    
    median = lambda : s[m] if bool(M&1) else (s[m-1]+s[m])*0.5

    # Sort it in increasing order and extract the median ("center" of the sorted window)
    s.sort()    
    medians = [median()]   

    # Now slide the window by one point to the right for each new position (each pass through 
    # the loop). Stop when the item in the right end of the deque contains the last item in seq
    for item in seq:
        old = d.popleft()          # pop oldest from left
        d.append(item)
        # push newest in from right
        try:
          del s[bisect_left(s, old)] # locate insertion point and then remove old 
              # insert newest such that new sort is not required
          insort(s, item)       
          medians.append(np.nanmedian(s))  
        except:
          medians.append(np.nan)
          pass

    return medians

#
def sciNot(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def searchUpper(EISCATnes, uBound, EISCATalts, upperBound, upperBoundH, iTime, indexSD): 
    for j in range(indexSD,len(EISCATnes)):
        try:
          if EISCATnes[j]>= uBound:
            # this is the upper bound
            upperBound[iTime]=EISCATnes[j-1]
            upperBoundH[iTime]=EISCATalts[j-1]
            return(upperBound, upperBoundH)
        except:
            pass
    upperBound[iTime]=np.nan
    upperBoundH[iTime]=np.nan
    return(upperBound, upperBoundH)

def searchLower(EISCATnes, lBound, EISCATalts, lowerBound, lowerBoundH, iTime): 
    for j in range(len(EISCATnes)):
        try:
          if EISCATnes[j]>= lBound:
            # this is the lower bound
            lowerBound[iTime]=EISCATnes[j]
            lowerBoundH[iTime]=EISCATalts[j]
            return(lowerBound, lowerBoundH)
        except:
            pass
    lowerBound[iTime]=np.nan
    lowerBoundH[iTime]=np.nan
    return(lowerBound, lowerBoundH)
#
# from processed guisdap data
dataEISCAT=readProcessedEISCAT.ReadData()
alt=dataEISCAT['alt'] # 42 altitudes between 76 km to 613 km
Times=dataEISCAT['stim']
ne=dataEISCAT['ne']
os.chdir('ProcessedEISCAT')
B=np.genfromtxt('EISCATB.txt')
os.chdir('..')
#
# load in SuperDARN
os.chdir('CUTLASS')
import readVelSTEREO2 as SD
#
dataSD=SD.readSUPERDARN()
#AllneSD=np.nansum(np.array(dataSD[0]),axis=0)
AllneSD=np.array(dataSD[0])
AllFreqsSD=np.array(dataSD[1])
AllfpSD=dataSD[2]
AllTimesSD=np.array(dataSD[3])
AllNs=np.array(dataSD[4])

os.chdir('..')
#
# now for each time average between gates 25 - 40
medNeSD15=np.zeros(len(AllTimesSD))


#
#temp=np.swapaxes(AllneSD,1,0)[25:35] # gates 25 45
#tempMed=np.nanmedian(temp, axis=0)
tempMed=np.abs(AllneSD)#*(AllNs**2)
#
# now get the closest eiscat value
# first get nearest time
NearestNeEIS=np.ones(len(AllTimesSD))*np.nan
NearestNeSD=np.ones(len(AllTimesSD))*np.nan
Neartesth=np.ones(len(AllTimesSD))*np.nan
NearestHeat=np.ones(len(AllTimesSD))*np.nan
Neartesth=np.ones(len(AllTimesSD))*np.nan
upperBound=np.ones(len(AllTimesSD))*np.nan
lowerBound=np.ones(len(AllTimesSD))*np.nan
upperBoundH=np.ones(len(AllTimesSD))*np.nan
lowerBoundH=np.ones(len(AllTimesSD))*np.nan

#
# resample the SuperDARN to 30 second increments between 10 to 12 UT
pdTime=pd.DatetimeIndex(AllTimesSD)
df1=pd.DataFrame({'data':tempMed}, index=pdTime[:-1])
rng2=pd.date_range('20150312', periods=2880, freq='30S')#[1200:1440] # 10 to 12 UT
dataR=np.array(df1['data'].resample('30S', how='mean').reindex(index=rng2, fill_value=np.nan))[1200:1440]

#
# heater density
heaterDens=(((6*2*np.pi*1e6)**2) * (8.85*1e-12) * (9.11*1e-31))/(1.6*1e-19)**2
for iTime in range(240):

        # array choices at time index
        EISCATalts=np.swapaxes(alt,1,0)[0]
        # ones in our proper ray tracing frame
        EISCATnes=np.swapaxes(ne,1,0)[iTime]
        #

        try:
            indexSD=  np.nanargmin(np.abs(EISCATnes-dataR[iTime]))
            indexHeater=np.nanargmin(np.abs(EISCATnes-heaterDens))
            
            Neartesth[iTime]=EISCATalts[indexSD]
            NearestHeat[iTime]=EISCATalts[indexHeater]
                        
            NearestNeEIS[iTime]=EISCATnes[indexSD]
            NearestNeSD[iTime]=dataR[iTime]
            
        except(ValueError):  # past EISCAT range
            NearestNeEIS[iTime]=np.nan
            NearestNeSD[iTime]=tempMed[iTime]
            Neartesth[iTime]=np.nan
            upperBound[iTime]=np.nan
            lowerBound[iTime]=np.nan
            upperBoundH[iTime]=np.nan
            lowerBoundH[iTime]=np.nan
            NearestHeat[iTime]=np.nan
        #
    #
os.chdir('/Users/loisks/Desktop/Functions/')
import colormaps as cmaps
os.chdir('/Users/loisks/Documents/ResearchProjects/SANSAProject/')    
# plot the results
fig=plt.figure()
ax=fig.add_subplot(111)
plt.subplots_adjust(right=0.75, top=0.92, bottom=0.11, left=0.11)
days = DayLocator(interval=1) 
hours = MinuteLocator(interval=30) 
hours2 = MinuteLocator(interval=10) 
daysFmt = DateFormatter('%H:%M')
fig.gca().xaxis.set_major_locator(hours)
fig.gca().xaxis.set_major_formatter(daysFmt)
fig.gca().xaxis.set_minor_locator(hours2)
font = {'family' : 'normal',
              'weight' : 'bold',
              'size'   : 22}
plt.rc('font', **font)
#
# first plot the EISCAT ne as background
ax.set_ylabel('Altitude [km]', fontsize=22, fontweight='bold')
time=dt.date2num(Times)
Altitudes=np.nanmean(alt, axis=1)
ax.set_ylim(100,300)
ax.set_xlim(time[0], time[-1])
X,Y=np.meshgrid(time, Altitudes)
#data=ma.masked_invalid(np.log10(ne))#.transpose()
data=ma.masked_invalid(ne)
ax.set_xlabel("Time", fontsize=20, fontweight='bold')
vmin=8*1e10
vmax=8*1e11
col=ax.pcolormesh(X,Y,data, cmap='viridis', vmin=vmin, vmax=vmax)
font = {'family' : 'normal',
    'weight' : 'bold',
    'size'   : 22}
plt.rc('font', **font)
cbaxes = fig.add_axes([0.8, 0.2, 0.03, 0.7])
cb = plt.colorbar(col, cax = cbaxes,format=FuncFormatter(sciNot),ticks=np.linspace(vmin,vmax,5)) 
cb.set_label('m$^{-3}$', fontsize=25)
fig.set_size_inches(13,9)
plt.draw()
#
# overplot
#    
time2=dt.date2num(AllTimesSD)
#points = np.array([time2,Neartesth]).T.reshape(-1, 1, 2)
#segments = np.concatenate([points[:-1], points[1:]], axis=1)
#lc = LineCollection(segments, cmap=plt.get_cmap('RdPu'),
#                    norm=plt.Normalize(10,12))
#cbaxes2= fig.add_axes([0.8, 0.55, 0.03, 0.37])
#lowerBoundH[lowerBoundH<=120]=np.nan
#upperBoundH[upperBoundH>=170]=np.nan
#ax.errorbar(time2,Neartesth my own median function

#


#newx=dt.date2num(Times)
newx=rng2.to_pydatetime()[1200:1440]
newy=np.array(Neartesth[0:240])
newy2=np.array(NearestHeat[0:240])
newc=np.array(NearestNeSD[0:240])
cc=np.ones(len(newy2))*heaterDens

col2=ax.scatter(newx, newy,c=newc, s=90, cmap= 'viridis', vmin=vmin, vmax=vmax, edgecolor='violet', lw=2)

col2=ax.scatter(newx, newy2,c=cc, s=40, cmap= 'viridis', vmin=vmin, vmax=vmax, edgecolor='red', lw=2)
#
# plot where we expect the heater to be

# find the nearest EISCAT density



ax.plot(newx, np.ones(len(newx))*120, lw=5, ls='--',c='magenta')
ax.plot(newx, np.ones(len(newx))*220, lw=5, ls='--',c='magenta')
#ax.axhspan(120, 220, alpha=0.3, color='k')
#cb2 = plt.colorbar(col2,cax = cbaxes2,format=FuncFormatter(sciNot),ticks=np.linspace(vmin,vmax,5)) 
#cb2.set_label('m$^{-3}$', fontsize=25)    
subdir_name='EISCAT_Spectrogram'

if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name) 
os.chdir(subdir_name)#
plt.savefig('All_20150312_EISCATSD_STEREO_CompareSpect.png')
os.chdir('..')
plt.close()


        
        
