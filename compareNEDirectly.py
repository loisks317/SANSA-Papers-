# compareNEDirectly.py
#
# from approximation of range gates from Ray Tracing
# look at the CUTLASS ne versus the EISCAT ne
#
# LKS, January 2016 SANSA
#
#
import numpy as np
import matplotlib.pyplot as plt
import os
import readProcessedEISCAT
import matplotlib.dates as dt
from matplotlib.dates import HourLocator, DayLocator, DateFormatter, MinuteLocator
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
#
# load in the IRI data
os.chdir('ModelData')
IRI=np.genfromtxt('IRIdata_range=880km.txt')
IRIalt=IRI[:,0]
IRIne=IRI[:,1]
os.chdir('..')

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
# load in EISCAT
dataEISCAT=readProcessedEISCAT.ReadData()
altE=np.nanmean(dataEISCAT['alt'], axis=1) # 42 altitudes between 76 km to 613 km
rightALTs=np.where((altE >= 100) & (altE <= 500))[0] # I think, need to confirm
TimesE=dt.date2num(dataEISCAT['stim'])
neE=dataEISCAT['ne']
#
# load in SuperDARN
os.chdir('CUTLASS')
import readCUTLASSvel as SD
#
dataSD=SD.readSUPERDARN()
# divide ne by 28 because scatter range of EISCAT = 8 km2
# SuperDARN = 15 * 15 = 225 km2
# 225 / 8 = 28.1
# so the scattering area of the SuperDARN Ne is 28 times larger than EISCAT
# we apply this later in the plotting section
# 
AllneSD=np.array(dataSD[0])

AllFreqsSD=np.array(dataSD[1])
AllfpSD=dataSD[2]
AllTimesSD=np.array(dataSD[3])
AllNs=np.array(dataSD[4])

os.chdir('..')
#
# so we need ne at certain freqs
#cFreqs=np.where((FreqsSD > 14500000) & (FreqsSD < 16000000))[0][:-1]
#cFreqs=np.where(FreqsSD > 0)[0][:-1] # just say all 
#
labels=['15_MHz', '16_MHz', '17_MHz', '18_MHz']
for iFREQ in range(len(AllFreqsSD)):
     # good SuperDARN data
     cNE=np.swapaxes(AllneSD[iFREQ],1,0) # to make it gates x Time
     cNE[cNE==0]=np.nan
     cTimes=dt.date2num(AllTimesSD[iFREQ])
     gates=np.linspace(28,40,13)
     #
     # side by side compare
     # plot for different alts
     cFactor=[]
     for iRNG in range(len(gates)):
       for iALT in range(len(rightALTs)):
          fig=plt.figure()
          ax=fig.add_subplot(111)
          plt.subplots_adjust(right=0.7, top=0.92, bottom=0.11, left=0.11)
          days = DayLocator(interval=1) 
          hours = MinuteLocator(interval=30) 
          hours2 = MinuteLocator(interval=60) 
          daysFmt = DateFormatter('%H:%M')
          fig.gca().xaxis.set_major_locator(hours)
          fig.gca().xaxis.set_major_formatter(daysFmt)
          fig.gca().xaxis.set_minor_locator(hours2)
          font = {'family' : 'normal',
                  'weight' : 'bold',
                  'size'   : 22}
          plt.rc('font', **font)
          t1=max(cTimes[0], TimesE[0])
          t2=min(cTimes[-1], TimesE[-1])
          # find where closest
          cIRIalt=min(range(len(IRIalt)), key=lambda i: abs(IRIalt[i]-altE[rightALTs[iALT]]))
          cNE[gates[iRNG]][cNE[gates[iRNG]]<0]=np.nan
          s1mask=np.isfinite(cNE[gates[iRNG]])
          medWin=8
          
           # for EISCAT data
          tempDat=cNE[gates[iRNG]][s1mask]
          #tempDat[tempDat>1e12]=np.nan
          #tempDat[tempDat<1.0e10]=np.nan
     
          ax.plot(runningMedian(cTimes[s1mask],medWin),runningMedian(tempDat,medWin), lw=3, marker='o',ls='-', color='blue')
          ax.plot(runningMedian(TimesE,medWin), runningMedian(neE[rightALTs[iALT]], medWin), lw=3, marker='o',ls='-', color='gold')
          ax.plot(TimesE, [IRIne[cIRIalt]]*len(TimesE), lw=3, color='red')
          #ax.set_yscale('log')
          ax.set_xlabel('Time')
#          ax.set_ylim(1e10.8,1e11.2)
          #ax.set_xlim(t1, t2)
          ax.set_ylabel('Electron Density [m$^{-3}$]')
          plt.legend(['SuperDARN', 'EISCAT', 'IRI'], bbox_to_anchor=(1.52, 0.7))
          fig.set_size_inches(13,9)
          subdir_name='Compare_EISCAT_CUTLASS'
          if not os.path.exists(subdir_name):
              os.umask(0) # unmask if necessary
              os.makedirs(subdir_name) 
          os.chdir(subdir_name)#
          plt.savefig(labels[iFREQ]+'_20150312Ne_alt='+str(altE[rightALTs[iALT]])+'_gate='+str(gates[iRNG])+'.png')
          os.chdir('..')
          plt.close()
