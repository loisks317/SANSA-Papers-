# readCUTLASSvel45km.py
#
# our goal is to read in the CUTLASS-Hankasalmi SuperDARN data
# for 3-12-2015 to start with. This data was provided by Mike Kosch.
# File format is a little weird, but the data columns are height
# not sure of the exact height profile right now. Working on that...
#
# updated to look at 20-35 instead of just 45 km
# also will do 1 minute windows instead of 40 second windows
#
#
# LKS January 2016, SANSA
# March 2016, Tromso update based on Gillies 2011 paper
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
     Times=[] # n points, where n = 75 here 
     Freqs=[] # n points
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

          
     tLen=len(lines)/3 + 1
     Data=np.zeros((tLen))   # n x altitude , in the dictionary
     # set a counter
     count=0
     dataCount=0
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
     # now i think the best way to do this is to create
     # four pandas arrays  based on frequency
     #
     pdTime=pd.DatetimeIndex(Times)
     Freqs=np.array(Freqs)
     rng2=pd.date_range(date, periods=720, freq='2T') # 10 to 12 UT
     # adjust the data
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
     f13band[Freqs<14800 ]=Data[Freqs< 14800]
     f15band[(Freqs<15800) & (Freqs>= 14800)]=Data[(Freqs<15800) & (Freqs>= 14800)]
     f16band[(Freqs<16800) & (Freqs>= 15800)]=Data[(Freqs<16800) & (Freqs>= 15800)]
     f17band[(Freqs<18000) & (Freqs>= 16800)]=Data[(Freqs<18000) & (Freqs>= 16800)]
     f18band[(Freqs<18500) & (Freqs>= 18000)]=Data[(Freqs<18500) & (Freqs>= 18000)]
     f19band[Freqs > 18500]=Data[Freqs>18500]
     #
     f13F[Freqs < 14800] = Freqs[Freqs<14800]
     f15F[(Freqs<15800) & (Freqs>= 14800)]=Freqs[(Freqs<15800) & (Freqs>= 14800)]
     f16F[(Freqs<16800) & (Freqs>= 15800)]=Freqs[(Freqs<16800) & (Freqs>= 15800)]
     f17F[(Freqs<18000) & (Freqs>= 16800)]=Freqs[(Freqs<18000) & (Freqs>= 16800)]
     f18F[(Freqs<18500) & (Freqs>= 18000)]=Freqs[(Freqs<18500) & (Freqs>= 18000)]
     f19F[Freqs>18500]=Freqs[Freqs>18500]
    
     #
     # now put into data frames
     df13=pd.DataFrame({'13':f13band}, index=pdTime)
     df15=pd.DataFrame({'15':f15band}, index=pdTime)
     df16=pd.DataFrame({'16':f16band}, index=pdTime)
     df17=pd.DataFrame({'17':f17band}, index=pdTime)
     df18=pd.DataFrame({'18':f18band}, index=pdTime)
     df19=pd.DataFrame({'19':f19band}, index=pdTime)

     # print out all of this
     print '13: ' + str(len(np.where(np.isnan(np.array(df13)) == False)[0]))
     print '15: ' + str(len(np.where(np.isnan(np.array(df15)) == False)[0]))
     print '16: ' + str(len(np.where(np.isnan(np.array(df16)) == False)[0]))
     print '17: ' + str(len(np.where(np.isnan(np.array(df17)) == False)[0]))
     print '18: ' + str(len(np.where(np.isnan(np.array(df18)) == False)[0]))
     print '19: ' + str(len(np.where(np.isnan(np.array(df19)) == False)[0]))


     dfF13=pd.DataFrame({'13':f13F}, index=pdTime)
     dfF15=pd.DataFrame({'15':f15F}, index=pdTime)
     dfF16=pd.DataFrame({'16':f16F}, index=pdTime)
     dfF17=pd.DataFrame({'17':f17F}, index=pdTime)
     dfF18=pd.DataFrame({'18':f18F}, index=pdTime)
     dfF19=pd.DataFrame({'19':f19F}, index=pdTime)

     #
     # and resample to 2 minutes
     data13=np.array(df13['13'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))
     data15=np.array(df15['15'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))
     data16=np.array(df16['16'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))
     data17=np.array(df17['17'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))
     data18=np.array(df18['18'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))
     data19=np.array(df19['19'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))

     freqs13=np.array(dfF13['13'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))
     freqs15=np.array(dfF15['15'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))
     freqs16=np.array(dfF16['16'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))
     freqs17=np.array(dfF17['17'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))
     freqs18=np.array(dfF18['18'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))
     freqs19=np.array(dfF19['19'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))

     datas=[data13,data15,data16,data17,data18,data19]
     freqies=[freqs13,freqs15,freqs16,freqs17,freqs18,freqs19]
     labels=['13','15','16','17','18','19']
     #
     # empty arrays
     Allfp=[[[ np.nan for i in range(len(rng2))] for j in range(6)] for k in range(6)]
     Maxfp=[[[ np.nan for i in range(len(rng2))] for j in range(6)] for k in range(6)]
     AllNe=[[[ np.nan for i in range(len(rng2))] for j in range(6)] for k in range(6)]
     MaxNe=[[[ np.nan for i in range(len(rng2))] for j in range(6)] for k in range(6)]

     # compare each possible combination
     for iFreq in range(6):
      for iFreq2 in range(6):
       for iTime in range(len(rng2)):

         v1=datas[iFreq][iTime]
         v2=datas[iFreq2][iTime]

         # this is bad data according to reviewer 2
         #
         f1=(freqies[iFreq][iTime])*1000
         f2=(freqies[iFreq2][iTime])*1000
         ff1=f1**2        
         ff2=f2**2
         if ff2 > ff1:
           Maxfp[iFreq][iFreq2][iTime]=f1
           if np.abs(v1) > np.abs(v2):
             # print v1
             # print v2
              v1 = np.nan
              v2 = np.nan
         else:
           Maxfp[iFreq][iFreq2][iTime]=f2
           if np.abs(v2) > np.abs(v1):
             # print v1
             # print v2
              v1 = np.nan
              v2 = np.nan              

         
         vm1=(v1)**2
         vm2=(v2)**2
         #mf=np.nanmin([ff1,ff2])

         Allfp[iFreq][iFreq2][iTime]=(np.abs(ff1*(1-(vm1/vm2)) / (1 - (vm1*ff1)/(vm2*ff2))))
         AllNe[iFreq][iFreq2][iTime]=0.0124*Allfp[iFreq][iFreq2][iTime]
         MaxNe[iFreq][iFreq2][iTime]=0.0124*(Maxfp[iFreq][iFreq2][iTime]**2)

         #if np.isnan(Allfp[iFreq][iFreq2][iTime]) == False:
         #     print "density parameters: "
         #     print ff1
         #     print ff2
         #     print (1-(vm1/vm2))
         #     print (1 - (vm1*ff1)/(vm2*ff2))
         #     print Allfp[iFreq][iFreq2][iTime]
         #     print AllNe[iFreq][iFreq2][iTime]
     
     return AllNe,  Allfp, rng2, datas, freqies, Maxfp, MaxNe
