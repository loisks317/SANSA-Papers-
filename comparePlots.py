# comparePlots.py
#
# module for plotting routine for the compare all
# needs to be really flexible because of all the crazy data we're given
#
# LKS, March 2016. Tromsoooooo
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import datetime
import readProcessedEISCAT
import matplotlib.dates as dt
from matplotlib.dates import HourLocator, DayLocator, DateFormatter, MinuteLocator
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import comparePlots as CP # for plotting routine
import pandas as pd
#
# contains:
#
# runningMedian(seq, M) * just needed a place to put this function
# correctedCompare(data,times,Edata,Etime,NeIRI,legendLabels,dtStart, dtEnd, color=['crimson', 'mediumseagreen', 'orange', 'blue', 'gold'], constant=8)
# differencePlot(data,times,Etime,legendLabels,dtStart, dtEnd, color=['crimson', 'mediumseagreen', 'orange', 'blue', 'gold'])

############################################################################

def correctedCompare(data,maxData,times,Edata,Etime,NeIRI,IRItime,legendLabels,dtStart, dtEnd,date, ylimlow,ylimhigh,color=['crimson', 'mediumseagreen', 'orange', 'blue', 'gold'], constant=1):
  # terms
  # data =  an array of data arrays
  # maxData= max Ne expected based on plasma frequency
  # times = array of time arrays matching data
  # Edata = EISCAT data
  # Etime = EISCAT time
  # NeIRI = IRI ne
  # IRItime = IRI time in datetime format
  # legendLabels= exactly what it sounds like
  # dtStart,dtEnd = x limits on plot
  # date = date of experiment
  # color = array of color matching data, default has 5 elements
  # constant = constant to correct the data, default is 8
  #
  fig=plt.figure()
  days = DayLocator(interval=1) 
  hours = MinuteLocator(interval=60) 
  hours2 = MinuteLocator(interval=10) 
  daysFmt = DateFormatter('%H:%M')
  fig.gca().xaxis.set_major_locator(hours)
  fig.gca().xaxis.set_major_formatter(daysFmt)
  fig.gca().xaxis.set_minor_locator(hours2)
  font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
  plt.rc('font', **font)
  ax=fig.add_subplot(111)
  plt.subplots_adjust(right=0.7, top=0.92, bottom=0.28, left=0.11)
  ax.set_ylabel('$N_e$ [m$^{-3}$]', fontsize=22, fontweight='bold')
  ax.set_ylim(ylimlow,ylimhigh)
  #y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
  #ax.yaxis.set_major_formatter(y_formatter)

  #ax.set_yscale('log')
  # plot the SuperDARN data

  for iplot in range(len(data)):
      ax.scatter(times[iplot],np.array(data[iplot])/constant, s=130, marker='>', c=color[iplot],edgecolors='None')
      
  #
  # plot the EISCAT data 
  ax.scatter(Etime, Edata, s=100, c='darkgray', marker='d', edgecolors='None')
  #
  # plot the IRI data 
  ax.plot(dt.date2num(IRItime), NeIRI, lw=3,ls='--', c='navy')
  #
  
  myyfmt = matplotlib.ticker.ScalarFormatter(useOffset=True)
  myyfmt._set_offset(1e12)
  ax.yaxis.set_major_formatter(myyfmt)
  plt.legend(legendLabels, bbox_to_anchor=[1.4, 0.9],scatterpoints=1, fontsize=20)
  # plot max lines
  for iplot in range(len(data)):
      temp=np.array(maxData[iplot])
      ax.plot(times[iplot],(temp/temp)*np.nanmedian(maxData[iplot]), lw=3, ls='--', c=color[iplot], label=None)
      
  ax.set_xlim(dtStart, dtEnd)
  ax.set_xlabel("UT Time", fontsize=20, fontweight='bold')
  fig.set_size_inches(13,9)
  subdir_name='Ne_Densities'
  if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name) 
  os.chdir(subdir_name)#
  plt.savefig(date+'_Comparison_All.png')
  os.chdir('..')
  plt.close()    

############################################################################

def differencePlot(data,times,Etime,legendLabels,dtStart, dtEnd,date, color=['crimson', 'mediumseagreen', 'orange', 'blue', 'gold']):
  # terms
  # data =  an array of data arrays
  # times = array of time arrays matching data
  # Etime = EISCAT time
  # legendLabels= exactly what it sounds like
  # dtStart,dtEnd = x limits on plot
  # date = date of experiment
  # color = array of color matching data, default has 5 elements
  #
  fig=plt.figure()
  days = DayLocator(interval=1) 
  hours = MinuteLocator(interval=60) 
  hours2 = MinuteLocator(interval=10) 
  daysFmt = DateFormatter('%H:%M')
  fig.gca().xaxis.set_major_locator(hours)
  fig.gca().xaxis.set_major_formatter(daysFmt)
  fig.gca().xaxis.set_minor_locator(hours2)
  font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
  plt.rc('font', **font)
  ax=fig.add_subplot(111)
  plt.subplots_adjust(right=0.7, top=0.92, bottom=0.28, left=0.11)
  ax.set_ylabel('CUTLASS / EISCAT N$_e$', fontsize=22, fontweight='bold')
  ax.plot(Etime, np.ones(len(Etime)), lw=2, color='k',ls='--')  
  for iplot in range(len(data)):
      ax.scatter(times[iplot], data[iplot], s=120, marker='>', c=color[iplot], label=legendLabels[iplot],edgecolors='None')
  plt.legend( bbox_to_anchor=[1.4, 0.7], fontsize=20, scatterpoints=1)
  ax.set_xlim(dtStart, dtEnd)
  ax.set_xlabel("UT Time", fontsize=20, fontweight='bold')
  fig.set_size_inches(13,9)
  ax.set_ylim(0, 20)
  subdir_name='Ne_Densities'
  if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name) 
  os.chdir(subdir_name)#
  plt.savefig(date+'_Diff_Comparison.png')
  os.chdir('..')
  plt.close()
  
############################################################################

def shiftStereo(dataStereo,timesStereo, dataShifting, timesShifting,legendLabels,dtStart, dtEnd, date,ylimlow,ylimhigh, color=['crimson', 'mediumseagreen', 'orange', 'blue', 'gold'], color2=['black', 'gray', 'silver', 'white', 'cream']):
  # terms
  # dataStereo =  an array of Stereo arrays
  # timesStereo = array of Stereo time arrays matching data
  # dataShifting = array of small frequency shifts
  # timesShifting = array of small times shifts
  # legendLabels= exactly what it sounds like
  # dtStart,dtEnd = x limits on plot
  # date = date of experiment
  # color = array of color matching data, default has 5 elements
  # color2 = array of greyscale, default has 5 elements
  #
  fig=plt.figure()
  days = DayLocator(interval=1) 
  hours = HourLocator(interval=1) 
  hours2 = MinuteLocator(interval=10) 
  daysFmt = DateFormatter('%H:%M')
  fig.gca().xaxis.set_major_locator(hours)
  fig.gca().xaxis.set_major_formatter(daysFmt)
  fig.gca().xaxis.set_minor_locator(hours2)
  font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
  plt.rc('font', **font)
  ax=fig.add_subplot(111)
  ax.set_ylim(ylimlow,ylimhigh)
  plt.subplots_adjust(right=0.7, top=0.92, bottom=0.28, left=0.11)
  ax.set_ylabel('$N_e$ [m$^{-3}$]', fontsize=22, fontweight='bold')
  #
  for iplot in range(len(dataStereo)):
      ax.scatter(timesStereo[iplot], dataStereo[iplot], s=120, marker='>', c=color[iplot],edgecolors='None')
  for iplot in range(len(dataShifting)):
         ax.scatter(timesShifting[iplot], dataShifting[iplot], s=80, marker='h', c=color2[iplot])

  plt.legend(legendLabels, bbox_to_anchor=[1.4, 0.9],scatterpoints=1, fontsize=20)
  ax.set_xlim(dtStart, dtEnd)
  class FixedOrderFormatter(ScalarFormatter):
    """Formats axis ticks using scientific notation with a constant order of 
    magnitude"""
    def __init__(self, order_of_mag=0, useOffset=True, useMathText=False):
        self._order_of_mag = order_of_mag
        ScalarFormatter.__init__(self, useOffset=useOffset, 
                                 useMathText=useMathText)
    def _set_orderOfMagnitude(self, range):
        """Over-riding this to avoid having orderOfMagnitude reset elsewhere"""
        self.orderOfMagnitude = self._order_of_mag

# Force the y-axis ticks to use 1e-9 as a base exponent 
  ax.yaxis.set_major_formatter(FixedOrderFormatter(12))

# Make the x-axis ticks formatted to 0 decimal places
 # ax.taxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.0f'))
  ax.set_xlabel("UT Time", fontsize=20, fontweight='bold')
  fig.set_size_inches(13,9)
  subdir_name='Ne_Densities'
  if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name) 
  os.chdir(subdir_name)#
  plt.savefig(date+'_shifting_stereo-compare.png')
  os.chdir('..')
  plt.close()





############################################################################
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
