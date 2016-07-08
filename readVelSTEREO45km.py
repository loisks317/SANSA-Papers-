# readCUTLASSvel45km.py
#
# our goal is to read in the CUTLASS-Hankasalmi SuperDARN data
# for 3-12-2015 to start with. This data was provided by Mike Kosch.
# File format is a little weird, but the data columns are height
# not sure of the exact height profile right now. Working on that...
#
# updated to look at 45 km range gates instead of 15km
#
#
# LKS January 2016, SANSA
# February 2016, Svalbard update for 45 km range gates
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
     Data=np.zeros((tLen,25))   # n x altitude , in the dictionary
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
             gate=0
             for iH in range(25):
                 t1=int(temp[gate])
                 t2=int(temp[gate+1])
                 t3=int(temp[gate+2])
                 if t1==10000:
                      t1=np.nan
                 if t2==10000:
                      t2=np.nan
                 if t3==10000:
                      t3=np.nan
                 Data[dataCount][iH]=np.nanmedian([t1,t2,t3])
                 gate=gate+3
             dataCount+=1
             count+=1
         elif count ==2:
             # skip this and reset the count
             count=0          
     #
     # bad data
     #
     # swap the axes for next loop
    # Data=np.swapaxes(Data,1,0)
     # do 20 second averages to start
     window=40
     AllNe1=np.ones((len(Freqs)/window,25))*np.nan
     Allfp1=np.ones((len(Freqs)/window,25))*np.nan
     AllNs1=np.ones((len(Freqs)/window,25))*np.nan
     AllNe2=np.ones((len(Freqs)/window,25))*np.nan
     Allfp2=np.ones((len(Freqs)/window,25))*np.nan
     AllNs2=np.ones((len(Freqs)/window,25))*np.nan
     AllNe3=np.ones((len(Freqs)/window,25))*np.nan
     Allfp3=np.ones((len(Freqs)/window,25))*np.nan
     AllNs3=np.ones((len(Freqs)/window,25))*np.nan
     AllVel1=np.ones((len(Freqs),25))*np.nan
     AllVel2=np.ones((len(Freqs),25))*np.nan
     AllVel3=np.ones((len(Freqs),25))*np.nan
     AllVel4=np.ones((len(Freqs),25))*np.nan
     AllFreq1=np.ones((len(Freqs)))*np.nan
     AllFreq2=np.ones((len(Freqs)))*np.nan
     AllFreq3=np.ones((len(Freqs)))*np.nan
     AllFreq4=np.ones((len(Freqs)))*np.nan

     AllNe=[AllNe1, AllNe2, AllNe3]
     Allfp=[Allfp1, Allfp2, Allfp3]
     AllNs=[AllNs1,AllNs2,AllNs3]

     #
     # resample
     pdTime=pd.DatetimeIndex(Times)
     #


     #
     # iterate through by height
     # find pairs where we have v measured at t1 and then t2
     # calculate the plasma frequency
     #
     # fp = sqrt[ (f1**2 (1- vm1**2/vm2**2))/(1- vm1**2 f1**2 /(vm2**2 f2**2))]
     # ns = sqrt[1 - fp**2/f**2]
     # Ne = 0.0124 fp**2
     #
     lbound=[14800,15800,16800,17800]
     ubound=[15800,16800,17800,18800]
     # need to sort by Frequencies
     for iH in range(25):
       for iSort in range(len(Freqs)):
         if ((Freqs[iSort] < 15800) & (Freqs[iSort] >= 14800)):
             AllVel1[iSort][iH]=Data[iSort][iH]
             AllFreq1[iSort]=Freqs[iSort]
         elif ((Freqs[iSort] < 16800) & (Freqs[iSort] >= 15800)):
             AllVel2[iSort][iH]=Data[iSort][iH]
             AllFreq2[iSort]=Freqs[iSort]
         elif ((Freqs[iSort] < 17800) & (Freqs[iSort] >= 16800)):
             AllVel3[iSort][iH]=Data[iSort][iH]
             AllFreq3[iSort]=Freqs[iSort]
         elif ((Freqs[iSort] < 18800) & (Freqs[iSort] >= 17800)):
             AllVel4[iSort][iH]=Data[iSort][iH]
             AllFreq4[iSort]=Freqs[iSort]
     AllVel=[AllVel1,AllVel2,AllVel3, AllVel4]
     AllFreq=[AllFreq1,AllFreq2,AllFreq3, AllFreq4]


     for iFreq in range(3):
      for iH in range(25):
       for iTime in range((len(Freqs)/window)-1):

         v1=np.nanmedian(AllVel[iFreq][iTime*window:iTime*window + window,iH])
         v2=np.nanmedian(AllVel[iFreq+1][iTime*window:iTime*window + window,iH])
         
         #

         f1=np.nanmedian(AllFreq[iFreq][iTime*window:iTime*window + window]*1000)
         f2=np.nanmedian(AllFreq[iFreq+1][iTime*window:iTime*window + window]*1000)
         # look at the 18 MHz
       #  if (np.isnan(f2)==True) & (iFreq==1):
       #      v2=np.nanmedian(AllVel[3][iTime*window:iTime*window + window,iH])
       #      f2=np.nanmedian(AllFreq[3][iTime*window:iTime*window + window]*1000)

         ff1=f1**2        
         ff2=f2**2
         vm1=(v1)**2
         vm2=(v2)**2        

         Allfp[iFreq][iTime][iH]=ff1*(1-(vm1/vm2)) / (1 - (vm1*ff1)/(vm2*ff2))
         AllNe[iFreq][iTime][iH]=(0.0124)*Allfp[iFreq][iTime][iH]
         AllNs[iFreq][iTime][iH]=np.sqrt(1-(Allfp[iFreq][iTime][iH]/((ff1+ff2)/2.0)))



     return AllNe, AllFreq, Allfp,  Times ,AllNs, AllVel
