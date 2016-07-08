# readCUTLASSvel.py
#
# our goal is to read in the CUTLASS-Hankasalmi SuperDARN data
# for 3-12-2015 to start with. This data was provided by Mike Kosch.
# File format is a little weird, but the data columns are height
# not sure of the exact height profile right now. Working on that... 
#
#
# LKS January 2016, SANSA 
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
     Data=np.swapaxes(Data,1,0)
     AllNe1=np.ones((240,75))*np.nan
     Allfp1=np.ones((240,75))*np.nan
     AllNs1=np.ones((240,75))*np.nan
     AllNe2=np.ones((240,75))*np.nan
     Allfp2=np.ones((240,75))*np.nan
     AllNs2=np.ones((240,75))*np.nan
     AllNe3=np.ones((240,75))*np.nan
     Allfp3=np.ones((240,75))*np.nan
     AllNs3=np.ones((240,75))*np.nan

     AllNe=[AllNe1, AllNe2, AllNe3]
     Allfp=[Allfp1, Allfp2, Allfp3]
     AllNs=[AllNs1,AllNs2,AllNs3]
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
     #
     lbound=[14800,15800,16800,17800]
     ubound=[15800,16800,17800,18800]
     Freqs=np.array(Freqs)
     for iFreq in range(3):
       for iH in range(75):
         dataR=Data[iH] # should be a list in time
         dataR[dataR==10000]=np.nan
 
         overlap=np.where((Freqs>=lbound[iFreq])&(Freqs<=ubound[iFreq]))[0]
         overlap2=np.where((Freqs>=lbound[iFreq+1])&(Freqs<=ubound[iFreq+1]))[0]
         df1=pd.DataFrame({'data':dataR[overlap]}, index=pdTime[overlap])
         df2=pd.DataFrame({'data':dataR[overlap2]}, index=pdTime[overlap2])
         dff1=pd.DataFrame({'freqs':Freqs[overlap]}, index=pdTime[overlap])
         dff2=pd.DataFrame({'freqs':Freqs[overlap2]}, index=pdTime[overlap2])
         rng2=pd.date_range('20150312', periods=2880, freq='30S')#[1200:1440] # 10 to 12 UT
         dataR=np.array(df1['data'].resample('30S', how='mean').reindex(index=rng2, fill_value=np.nan))
         dataR2=np.array(df2['data'].resample('30S', how='mean').reindex(index=rng2, fill_value=np.nan))
         frf1=np.array(dff1['freqs'].resample('30S', how='mean').reindex(index=rng2, fill_value=np.nan))
         frf2=np.array(dff2['freqs'].resample('30S', how='mean').reindex(index=rng2, fill_value=np.nan))
         #
         # get the right time out of it 
         for iR in range(240):
              # check to make sure the frequencies are different
            
              f1=frf1[1200:1440][iR]*1000
              f2=frf2[1200:1440][iR]*1000
              v1= dataR[1200:1440][iR]
              v2= dataR2[1200:1440][iR]

              ff1=f1**2        
              ff2=f2**2
              vm1= (v1)**2
              vm2= (v2)**2                 
              Allfp[iFreq][iR][iH]=(ff1)*(1-(vm1/vm2)) / (1 - (vm1*ff1)/(vm2*ff2))
              AllNe[iFreq][iR][iH]=((2*np.pi)**2)*Allfp[iFreq][iR][iH]/(3175.2)
              AllNs[iFreq][iR][iH]=np.sqrt(1-(Allfp[iFreq][iR][iH]/((ff1+ff2)/2.0)))


     return AllNe, Freqs, Allfp,  rng2.to_pydatetime()[1200:1440] ,AllNs, Data
