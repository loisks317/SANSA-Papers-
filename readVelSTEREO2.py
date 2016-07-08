# readVelSTEREO2.py
#
# This method is different in that it takes a median range gate and then
# compares it with the next step
#
#
# LKS February 2016, Svalbard 
#
#
# First we load specific modules
import numpy as np # always should load numpy
import datetime # nice for keeping track of time
import matplotlib.pyplot as plt # best package ever
import os # to switch directories
# need these for time intervals later
import matplotlib.dates as dt
from matplotlib.dates import HourLocator, DayLocator, DateFormatter, MinuteLocator
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import pandas as pd
#
boltzmann=1.38*1e-23
me=9.11*1e-31
#
def readSUPERDARN():
     # open the file and read in the lines
     # the '2:' means we start after the header on the
     # third line and we go to the end of the array
     file=open('20150312f_vel.ascii', 'rb')
     lines=file.readlines()[2:]
     #
     # stuff to hold the data
     # [] = empty list, nothing fancy
     # {} = dictionary, this is fancy. we're going to use it
     #      to hold the data arrays at each time step
     # 
     Times=[] # n points, where n = 75 here 
     Freqs=[] # n points
     tLen=9980
     Data=np.zeros((tLen,75))   # n x altitude , in the dictionary
     #
     # loop through the data file
     # pattern is
     # Ln 1: HHMM SSs (071) frequency
     # Ln 2: data
     # Ln 3: \n
     #
     # set a counter
     count=0
     dataCount=0
     for iLine in lines:
         if count == 0: # get the header info
             temp=iLine.split()
             HrMn=str(temp[0])
             Hr=int(HrMn[0:2])
             Mn=int(HrMn[2:])
             Ss=int(str(temp[1])[:-1])
             dT=datetime.datetime(2015,3,12,Hr,Mn,Ss)
             # append on to array
             Times.append(dT)
             #
             # get the freq
             Freqs.append(int(temp[3]))
             count+=1
         elif count ==1:
             temp=iLine.split()
             for iTemp in range(len(temp)):
                 Data[dataCount][iTemp]=int(temp[iTemp])
             dataCount+=1
             count+=1
         elif count ==2:
             # skip this and reset the count
             count=0          
     #
     # bad data
     #
     # swap the axes for next loop
     #Data=np.swapaxes(Data,1,0)
     AllNe=np.ones((len(Times)-1))*np.nan
     Allfp=np.ones((len(Times)-1))*np.nan
     AllNs=np.ones((len(Times)-1))*np.nan

     #
     # resample
     pdTime=pd.DatetimeIndex(Times)
     

     #
     # iterate through by height
     # find pairs where we have v measured at t1 and then t2
     # calculate the plasma frequency
     #
     # fp = sqrt[ (f1**2 (1- vm1**2/vm2**2))/(1- vm1**2 f1**2 /(vm2**2 f2**2))]
     # ns = sqrt[1 - fp**2/f**2]
     # Ne = 0.0124 fp**2
     for iTime in range(len(Times)-1):
         print np.shape(Data)
         dataR=Data[iTime][25:38] # should be a list in time
         dataR2=Data[iTime+1][25:38]
         dataR[dataR==10000]=np.nan
         dataR2[dataR2==10000]=np.nan
         #
         # get the nans
         mData=np.nanmean(dataR)
         mData2=np.nanmean(dataR2)
    
         f1=Freqs[iTime]*1000
         f2=Freqs[iTime+1]*1000

         v1= mData
         v2= mData2

         ff1=f1**2        
         ff2=f2**2
         vm1= (v1)**2
         vm2= (v2)**2
         print vm1
         print vm2
         Allfp[iTime]=np.abs((ff1)*(1-(vm1/vm2)) / (1 - (vm1*ff1)/(vm2*ff2)))
         AllNe[iTime]=((2*np.pi)**2)*Allfp[iTime]/(3175.2)
         AllNs[iTime]=np.sqrt(1-(np.abs(Allfp[iTime])/((ff1+ff2)/2.0)))


     return AllNe, Freqs, Allfp, Times ,AllNs, Data
