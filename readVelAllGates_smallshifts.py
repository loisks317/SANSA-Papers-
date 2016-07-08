# readVelAllGates_smallshifts.py
#
# our goal is to read in the CUTLASS-Hankasalmi SuperDARN data
# for 3-12-2015 to start with. This data was provided by Mike Kosch.
# File format is a little weird, but the data columns are height
# not sure of the exact height profile right now. Working on that...
#
# the updated code only looks at small shifts for comparison
#
# updated to look at 20-35 instead of just 45 km
# also will do 1 minute windows instead of 40 second windows
#
#
# LKS January 2016, SANSA
# March 2016, Tromso update based on Gillies 2011 paper
# 
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
def readSUPERDARN(date):
     # open the file and read in the lines
     # the '2:' means we start after the header on the
     # third line and we go to the end of the array
     os.chdir('Data')
     if date[0:4]=='2016':
          file1=open(date+'f_vel_chA.ascii','rb')
          file2=open(date+'f_vel_chB.ascii','rb')
          lines1=file1.readlines()[2:]
          lines2=file2.readlines()[2:]         
          flag1=np.ones(len(lines1))
          flag0=np.zeros(len(lines2))
          lines=list(lines1)+list(lines2)
          flag=list(flag1)+list(flag0)
     elif date[0:4]=='2015':
          file=open(date+'f_vel.ascii','rb')
          lines=file.readlines()[2:]
          flag=np.zeros(len(lines))
     os.chdir('..')
     #
     # stuff to hold the data
     # [] = empty list, nothing fancy
     # {} = dictionary, this is fancy. we're going to use it
     #      to hold the data arrays at each time step
     # 
     Times=[] # n points, where n = 75 here 
     Freqs=[] # n points
     tLen=len(lines)/3 +1
     Data=np.zeros((tLen))   # n x altitude , in the dictionary
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
     totalfp={}
     totalne={}
     totaltime={}
     for iLine in range(len(lines)):
          if ((flag[iLine]==0) and (flag[iLine-1]==1)):
               print ('switch to B')
               count=0
          if count == 0: # get the header info
               temp=lines[iLine].split()
               HrMn=str(temp[0])
               try:
                         Hr=int(HrMn[0:2])
                         Mn=int(HrMn[2:])
                         Ss=int(str(temp[1])[:-1])
                         if Ss > 59:
                              Ss=59
                         dT=datetime.datetime(int(date[0:4]),int(date[4:6]),int(date[6:8]),Hr,Mn,Ss)               
                         Times.append(dT)

                         Freqs.append(int(temp[3]))
               except(ValueError):
                         Times.append(np.nan)
                         Freqs.append(np.nan)
                         print 'Error in reading in'
               count+=1
          elif count ==1:
                    temp=lines[iLine].split()
                    # check the flag
                    # right now they all seem to be 20:40 may change later
                    if flag[iLine]==0: # 15 km resolution
                         t1=[float(i) for i in temp[29:34]]
                    if flag[iLine]==1: # 45 km resolution
                         t1=[float(i) for i in temp[29:34]]
                    t1=np.array(t1)
                    a=np.where(t1==10000)[0]
                    t1[a]=np.nan
                    Data[dataCount]=np.nanmean(t1)
                    dataCount+=1
                    count+=1
          elif count ==2:
             # skip this and reset the count
             count=0          

     #
     # adjust the data
     Freqs=np.array(Freqs)
     Data=np.array(Data)
     f13band=np.ones(len(Data))*np.nan
     f15band=np.ones(len(Data))*np.nan
     f16band=np.ones(len(Data))*np.nan
     f17band=np.ones(len(Data))*np.nan
     f18band=np.ones(len(Data))*np.nan
     f19band=np.ones(len(Data))*np.nan
     #
     f13F=np.ones(len(Data))*np.nan
     f15F=np.ones(len(Data))*np.nan
     f16F=np.ones(len(Data))*np.nan
     f17F=np.ones(len(Data))*np.nan
     f18F=np.ones(len(Data))*np.nan
     f19F=np.ones(len(Data))*np.nan
     #
     # set up to do simultaneous comparison
     f13band[Freqs<14800 ]=Data[Freqs< 14800]
     f15band[(Freqs<15800) & (Freqs>= 14800)]=Data[(Freqs<15800) & (Freqs>= 14800)]
     f16band[(Freqs<16800) & (Freqs>= 15800)]=Data[(Freqs<16800) & (Freqs>= 15800)]
     f17band[(Freqs<18000) & (Freqs>= 16800)]=Data[(Freqs<18000) & (Freqs>= 16800)]
     f18band[(Freqs<18800) & (Freqs>= 18000)]=Data[(Freqs<18800) & (Freqs>= 18000)]
     f19band[Freqs > 18800]=Data[Freqs > 18800]
     #
     f13F[Freqs < 14800] = Freqs[Freqs<14800]
     f15F[(Freqs<15800) & (Freqs>= 14800)]=Freqs[(Freqs<15800) & (Freqs>= 14800)]
     f16F[(Freqs<16800) & (Freqs>= 15800)]=Freqs[(Freqs<16800) & (Freqs>= 15800)]
     f17F[(Freqs<18000) & (Freqs>= 16800)]=Freqs[(Freqs<18000) & (Freqs>= 16800)]
     f18F[(Freqs<18800) & (Freqs>= 18000)]=Freqs[(Freqs<18800) & (Freqs>= 18000)]
     f19F[Freqs>18800]=Freqs[Freqs>18800]
     #
     f13Times=Times
     f15Times=Times
     f16Times=Times
     f17Times=Times
     f18Times=Times
     f19Times=Times

     # do not both with putting into a data frame
     bands=[f13band, f15band, f16band, f17band, f18band,f19band]
     timeT=[f13Times, f15Times, f16Times, f17Times, f18Times, f19Times]
     freqs=[f13F, f15F, f16F, f17F, f18F,f19F]
     
     for i in range(6):
         tempBand=np.array(bands[i])
         tempTime=np.array(timeT[i])
         tt=tempTime[~np.isnan(tempBand)]
         ff=np.array(freqs[i])[~np.isnan(tempBand)]
         bb=tempBand[~np.isnan(tempBand)]

         bands[i]=bb
         timeT[i]=tt
         #
         # now calculate the ne
         if len(bb) !=0:
           Allfp=np.ones(len(bb)-1)*np.nan
           AllNe=np.ones(len(bb)-1)*np.nan
           for iPoint in range(len(bb)-1):
             v1=bb[iPoint]
             v2=bb[iPoint+1]
             f1=ff[iPoint]*1000
             f2=ff[iPoint+1]*1000
             ff1=f1**2
             ff2=f2**2
             vm1=v1**2
             vm2=v2**2
             Allfp[iPoint]=(np.abs(ff1*(1-(vm1/vm2)) / (1 - (vm1*ff1)/(vm2*ff2))))
             AllNe[iPoint]=0.0124*Allfp[iPoint]
           totalfp[i]=Allfp
           totalne[i]=AllNe
           totaltime[i]=tt
         else:
           totalfp[i]=[np.nan]
           totalne[i]=[np.nan]
           totaltime[i]=[np.nan]
           
    
     return totalne,  totalfp, totaltime
