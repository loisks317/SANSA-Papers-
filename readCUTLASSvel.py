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
     #noDats=np.where(Data==10000)[0]
     #print(len(noDats))
     
     #Data[noDats]=np.nan
     
     #
     # swap the axes for next loop
     Data=np.swapaxes(Data,1,0)
     #
     # iterate through by height
     # find pairs where we have v measured at t1 and then t2
     # calculate the plasma frequency
     #
     # fp = sqrt[ (f1**2 (1- vm1**2/vm2**2))/(1- vm1**2 f1**2 /(vm2**2 f2**2))]
     # ns = sqrt[1 - fp**2/f**2]
     # Ne = 0.0124 fp**2
     #
     Freqs=np.array(Freqs) 

     # separate into proper cats
     freq15i = np.where((Freqs>= 15000) & (Freqs < 16000))[0]
     freq16i = np.where((Freqs>= 16000) & (Freqs < 17000))[0]
     freq17i = np.where((Freqs>= 17000) & (Freqs < 18000))[0]
     freq18i = np.where((Freqs>= 18000) & (Freqs < 19000))[0]
     freqs15=Freqs[freq15i]*1000
     times15=np.array(Times)[freq15i]
     data15=Data[:,freq15i]
     freqs16=Freqs[freq16i]*1000
     times16=np.array(Times)[freq16i]
     data16=Data[:,freq16i]
     freqs17=Freqs[freq17i]*1000
     times17=np.array(Times)[freq17i]
     data17=Data[:,freq17i]
     freqs18=Freqs[freq18i]*1000
     times18=np.array(Times)[freq18i]
     data18=Data[:,freq18i]

     fp15=np.ones((len(freqs15),75))*np.nan
     Ne15=np.ones((len(freqs15),75))*np.nan
     Ns15=np.ones((len(freqs15),75))*np.nan
     fp16=np.ones((len(freqs16),75))*np.nan
     Ne16=np.ones((len(freqs16),75))*np.nan
     Ns16=np.ones((len(freqs16),75))*np.nan
     fp17=np.ones((len(freqs17),75))*np.nan
     Ne17=np.ones((len(freqs17),75))*np.nan
     Ns17=np.ones((len(freqs17),75))*np.nan
     fp18=np.ones((len(freqs18),75))*np.nan
     Ne18=np.ones((len(freqs18),75))*np.nan
     Ns18=np.ones((len(freqs18),75))*np.nan

     AllFreqs=[freqs15,freqs16,freqs17,freqs18]
     AllData=[data15,data16,data17,data18]
     AllTime=[times15,times16,times17,times18]
     AllNe=[Ne15,Ne16,Ne17,Ne18]
     AllNs=[Ns15,Ns16,Ns17,Ns18]
     Allfp=[fp15,fp16,fp17,fp18]

     for iFREQ in range(len(AllFreqs)):
       for iH in range(75):
         dataR=AllData[iFREQ][iH]
         iR=0
         tLen=len(AllFreqs[iFREQ])
         while iR < (tLen-1):
             # check to see if we have a two adjacent velocity measures
              tfre=AllFreqs[iFREQ][iR]
              try:
                while dataR[iR]==10000:
                   iR+=1
              except(IndexError):
                   break
              # now we have a first value
              v1=dataR[iR]
              m1=iR # memory IR
              try:
                while dataR[iR+1]==10000: # bad data 
                   iR+=1
              except(IndexError):
                   break
              if iR < tLen-2: # break out of loop
                 v2= dataR[iR+1]                   
                 ff1=AllFreqs[iFREQ][m1]**2        
                 ff2=AllFreqs[iFREQ][iR+1]**2
                 vm1= (v1)**2
                 vm2= (v2)**2
                 
                 Allfp[iFREQ][m1][iH]=(ff1)*(1-(vm1/vm2)) / (1 - (vm1*ff1)/(vm2*ff2))
                 AllNe[iFREQ][m1][iH]=(1/(4*np.pi**2))*(1.0e6)*Allfp[iFREQ][m1][iH]/(8980**2)
                 AllNs[iFREQ][m1][iH]=np.sqrt(1-(Allfp[iFREQ][m1][iH]/((ff1+ff2)/2.0)))

              # count up
              iR+=1

     
     return AllNe, AllFreqs, Allfp,  AllTime ,AllNs, AllData
