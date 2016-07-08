# compareAll.py
#
# compare IRI, EISCAT, Ne frequency shifting method, and Ne_Stereo method
# then compare Ne_Stereo Method at 1 MHz interval vs 2 MHz interval
#
# LKS February 2016, Svalbard. Doing the 1-7 AM work shift again ;)
#
import numpy as np
import matplotlib.pyplot as plt
import os
import datetime
import readProcessedEISCAT
import matplotlib.dates as dt
from matplotlib.dates import HourLocator, DayLocator, DateFormatter, MinuteLocator
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import comparePlots as CP # for plotting routine
import pandas as pd
from scipy import stats
#
# globals
me=9.11*1e-31
q=1.6*1e-19
eps0=8.85*1e-12
date='20150312'
#date='20160303'
ylimlow=4*1e10
ylimhigh=5*1e12
ylimlow2=1*1e11
ylimhigh2=1*1e13
if date[0:4]=='2016':
    dtStart=datetime.datetime.strptime(date, '%Y%m%d')
    HrStart=14
    HrEnd=18
    dtStart+=datetime.timedelta(seconds=3600*HrStart)
    dtEnd=dtStart+datetime.timedelta(seconds=3600*4)
    iSt2=HrStart*30
    eSt2=HrEnd*30
    if date=='20160303':
        iSt=1681
        eSt=2160
        alt=220
    elif date=='20160302':
        iSt=1680
        eSt=2166        
    elif date=='20160304':
        iSt=1696
        eSt=2160
elif date[0:4]=='2015':
    dtStart=datetime.datetime(2015,3,12,10,0,0)
    dtEnd=datetime.datetime(2015,3,12,12,0,0)
    iSt=1200
    eSt=1440
    iSt2=10*30
    eSt2=12*30
    alt=240
IRIhours=5
IRItime=[ dtStart+datetime.timedelta(seconds=3600*i) for i in range(IRIhours)]        
#
# load in the IRI data
os.chdir('ModelData')
IRI=np.genfromtxt('IRI_'+str(alt)+'_'+date+'.txt')
#IRIalt=IRI[:,0]
NeIRI=IRI[:,1]
#
# IRI 2
IRI2=np.genfromtxt('IRIdata_'+date+'.txt')
IRIprof=IRI2[:,1]
IRIprofalts=IRI2[:,0]
os.chdir('..')
#
# load in EISCAT
dataEISCAT=readProcessedEISCAT.ReadData(date)
altE=np.nanmean(dataEISCAT['alt'], axis=1) # 42 altitudes between 76 km to 613 km
rightALTs=np.where((altE >= 100) & (altE <= 500))[0] # I think, need to confirm
Etime=dt.date2num(dataEISCAT['stim'])
neE=dataEISCAT['ne']
#

os.chdir('..')
plasmaHeater=0.0124*((6.2) * 1e6)**2 
#
# load in SuperDARN
os.chdir('CUTLASS')
import readVelAllGates as Stereo
import readVelAllGates_smallshifts as ss
dataStereo=Stereo.readSUPERDARN(date)
datass=ss.readSUPERDARN(date)
os.chdir('..')
NeEiscat=neE[min(range(len(altE)), key=lambda i: abs(altE[i]-alt))]
#NeIRI=IRIne[min(range(len(IRIalt)), key=lambda i: abs(IRIalt[i]-alt))]
# adjust?
if date=='20150312':
    NeEiscat=NeEiscat[np.argsort(Etime)]
    Etime=np.sort(Etime)
#dfE=pd.DataFrame({'E':NeEiscat}, index=dt.num2date(Etime))
rng2=pd.date_range(date, periods=720, freq='2T') 
cTime=dt.date2num(rng2.to_pydatetime())
df2=pd.DataFrame(data=NeEiscat, index=Etime)


#
# now make it a choice
if date=='20150312':
    # first experiment
    #
    # assume range gate 11 for SuperDARN
    NeStereo1516=np.abs(dataStereo[0][1][2]) # pick the 15 to 16 MHz comparison
    NeStereo1617=np.abs(dataStereo[0][2][3])
    NeStereo1618=np.abs(dataStereo[0][2][4])
    NeStereo1718=np.abs(dataStereo[0][3][4])

    # print the values
    print '15-16: ' + str(len(np.where(np.isnan(NeStereo1516) == False)[0]))
    print '16-17: ' + str(len(np.where(np.isnan(NeStereo1617) == False)[0]))
    print '16-18: ' + str(len(np.where(np.isnan(NeStereo1618) == False)[0]))
    print '17-18: ' + str(len(np.where(np.isnan(NeStereo1718) == False)[0]))   

    #
    maxNeStereo1516=np.abs(dataStereo[6][1][2])
    maxNeStereo1617=np.abs(dataStereo[6][2][3])
    maxNeStereo1618=np.abs(dataStereo[6][2][4])
    maxNeStereo1718=np.abs(dataStereo[6][3][4])
    
    Ness15=np.abs(datass[0][1])
    Ness16=np.abs(datass[0][2])
    Ness17=np.abs(datass[0][3])
    Ness18=np.abs(datass[0][4])

    TimesStereo=dt.date2num(dataStereo[2].to_pydatetime())
    Times15ss=datass[2][1]
    Times16ss=datass[2][2]
    Times17ss=datass[2][3]
    Times18ss=datass[2][4]
    #
    # resample the superDARN small shift data
    df15=pd.DataFrame({'15':Ness15}, index=pd.DatetimeIndex(Times15ss[:-1]))
    df16=pd.DataFrame({'16':Ness16}, index=pd.DatetimeIndex(Times16ss[:-1]))
    df17=pd.DataFrame({'17':Ness17}, index=pd.DatetimeIndex(Times17ss[:-1]))
    df18=pd.DataFrame({'18':Ness18}, index=pd.DatetimeIndex(Times18ss[:-1]))

    # resample to 2 minutes
    mNess15=np.array(df15['15'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))
    mNess16=np.array(df16['16'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))
    mNess17=np.array(df17['17'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))
    mNess18=np.array(df18['18'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))

    p1=dt.date2num(dtStart)
    p2=dt.date2num(dtEnd)    
    mTimesStereo=np.array(TimesStereo)
    
    # get the right superDARN and EISCAT times
    SDtime=mTimesStereo[np.where(mTimesStereo>=p1)[0][0]:np.where(mTimesStereo<=p2)[0][-1]]
    #Eiscattime=mTimesEiscat[np.where(mTimesEiscat>=p1)[0][0]:np.where(mTimesEiscat<=p2)[0][-1]]
    SDdata1516=np.array(NeStereo1516[np.where(mTimesStereo>=p1)[0][0]:np.where(mTimesStereo<=p2)[0][-1]])
    SDdata1617=np.array(NeStereo1617[np.where(mTimesStereo>=p1)[0][0]:np.where(mTimesStereo<=p2)[0][-1]])
    SDdata1618=np.array(NeStereo1618[np.where(mTimesStereo>=p1)[0][0]:np.where(mTimesStereo<=p2)[0][-1]])
    SDdata1718=np.array(NeStereo1718[np.where(mTimesStereo>=p1)[0][0]:np.where(mTimesStereo<=p2)[0][-1]])
    #
    # and for the max value
    maxSDdata1516=np.array(maxNeStereo1516[np.where(mTimesStereo>=p1)[0][0]:np.where(mTimesStereo<=p2)[0][-1]])
    maxSDdata1617=np.array(maxNeStereo1617[np.where(mTimesStereo>=p1)[0][0]:np.where(mTimesStereo<=p2)[0][-1]])
    maxSDdata1618=np.array(maxNeStereo1618[np.where(mTimesStereo>=p1)[0][0]:np.where(mTimesStereo<=p2)[0][-1]])
    maxSDdata1718=np.array(maxNeStereo1718[np.where(mTimesStereo>=p1)[0][0]:np.where(mTimesStereo<=p2)[0][-1]])

    df1516=pd.DataFrame(data=SDdata1516, index=SDtime)
    df1617=pd.DataFrame(data=SDdata1617, index=SDtime)
    df1618=pd.DataFrame(data=SDdata1618, index=SDtime)
    df1718=pd.DataFrame(data=SDdata1718, index=SDtime)

    redf1516=np.array(df2.reindex(df1516.index, method='nearest'))[:,0]
    redf1617=np.array(df2.reindex(df1617.index, method='nearest'))[:,0]
    redf1618=np.array(df2.reindex(df1618.index, method='nearest'))[:,0]
    redf1718=np.array(df2.reindex(df1718.index, method='nearest'))[:,0]

    normalize1516=np.array((SDdata1516)/(1.0*redf1516))
    normalize1617=np.array((SDdata1617)/(1.0*redf1617))
    normalize1618=np.array((SDdata1618)/(1.0*redf1618))
    normalize1718=np.array((SDdata1718)/(1.0*redf1718))

    # final arrays
    cDataStereo=[NeStereo1516, NeStereo1617, NeStereo1618, NeStereo1718]
    maxDataStereo=[maxNeStereo1516, maxNeStereo1617, maxNeStereo1618, maxNeStereo1718]
    cTimeStereo=[mTimesStereo]*4
    dDataStereo=[normalize1516,normalize1617,normalize1618,normalize1718]
    dTimeStereo=[SDtime]*4
    ssData=[mNess15,mNess16,mNess17,mNess18]
    ssTime=[rng2.to_pydatetime()]*4
    legendLabels=['15 MHz', '16 MHz', '17 MHz', '18 MHz']
    legendLabels2 =['15 - 16 MHz', '16 - 17 MHz', '16 - 18 MHz', '17 - 18 MHz']
    legendLabels3=legendLabels2 + legendLabels

#
# for next experiment
if date=='20160302':
    # first experiment
    #
    # assume range gate 11 for SuperDARN
    NeStereo1516=np.abs(dataStereo[0][0][1]) # pick the 15 to 16 MHz comparison

    Ness15=np.abs(datass[0][1])
    Ness16=np.abs(datass[0][2])


    TimesStereo=dt.date2num(dataStereo[2].to_pydatetime())
    Times15ss=datass[2][1]
    Times16ss=datass[2][2]
    #
    # resample the superDARN small shift data
    df15=pd.DataFrame({'15':Ness15}, index=pd.DatetimeIndex(Times15ss[:-1]))
    df16=pd.DataFrame({'16':Ness16}, index=pd.DatetimeIndex(Times16ss[:-1]))

    # resample to 2 minutes
    mNess15=np.array(df15['15'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))
    mNess16=np.array(df16['16'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))

    p1=dt.date2num(dtStart)
    p2=dt.date2num(dtEnd)    
    mTimesStereo=np.array(TimesStereo)
    # get the right superDARN and EISCAT times
    SDtime=mTimesStereo[np.where(mTimesStereo>=p1)[0][0]:np.where(mTimesStereo<=p2)[0][-1]]
    #Eiscattime=mTimesEiscat[np.where(mTimesEiscat>=p1)[0][0]:np.where(mTimesEiscat<=p2)[0][-1]]
    SDdata1516=np.array(NeStereo1516[np.where(mTimesStereo>=p1)[0][0]:np.where(mTimesStereo<=p2)[0][-1]])
   
    df1516=pd.DataFrame(data=SDdata1516, index=SDtime)
    #Edata=np.array(mNeEiscat[np.where(mTimesEiscat>=p1)[0][0]:np.where(mTimesEiscat<=p2)[0][-1]])
    #df2=pd.DataFrame(data=NeEiscat,index=Etime)
    redf1516=np.array(df2.reindex(df1516.index, method='nearest'))[:,0]

    normalize1516=np.array((SDdata1516)/(1.0*redf1516))

    # final arrays
    cDataStereo=[NeStereo1516]
    cTimeStereo=[mTimesStereo]
    dDataStereo=[normalize1516]
    dTimeStereo=[SDtime]
    ssData=[mNess15,mNess16]
    ssTime=[rng2.to_pydatetime()]*2
    legendLabels=['15 MHz', '16 MHz']
    legendLabels2=['15-16 MHz']
    legendLabels3=legendLabels2 + legendLabels
#
#
# for next experiment
if date=='20160303':
    # first experiment
    #
    # assume range gate 11 for SuperDARN
    
    NeStereo1315=np.abs(dataStereo[0][0][1])
    NeStereo1316=np.abs(dataStereo[0][0][2])
    NeStereo1516=np.abs(dataStereo[0][1][2]) # pick the 15 to 16 MHz comparison

    # print the values
    print '13-15: ' + str(len(np.where(np.isnan(NeStereo1315) == False)[0]))
    print '13-16: ' + str(len(np.where(np.isnan(NeStereo1316) == False)[0]))
    print '15-16: ' + str(len(np.where(np.isnan(NeStereo1516) == False)[0]))  


    
    maxNeStereo1315=np.abs(dataStereo[6][0][1])
    maxNeStereo1316=np.abs(dataStereo[6][0][2])
    maxNeStereo1516=np.abs(dataStereo[6][1][2]) # pick the 15 to 16 MHz comparison

    Ness13=np.abs(datass[0][0])
    Ness15=np.abs(datass[0][1])
    Ness16=np.abs(datass[0][2])

    TimesStereo=dt.date2num(dataStereo[2].to_pydatetime())
    Times13ss=datass[2][0]
    Times15ss=datass[2][1]
    Times16ss=datass[2][2]
    #
    # resample the superDARN small shift data
    df13=pd.DataFrame({'13':Ness13}, index=pd.DatetimeIndex(Times13ss[:-1]))
    df15=pd.DataFrame({'15':Ness15}, index=pd.DatetimeIndex(Times15ss[:-1]))
    df16=pd.DataFrame({'16':Ness16}, index=pd.DatetimeIndex(Times16ss[:-1]))

    # resample to 2 minutes
    mNess13=np.array(df13['13'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))
    mNess15=np.array(df15['15'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))
    mNess16=np.array(df16['16'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))

    # NEED TO MAKE SURE INDEXING IS RIGHT.... RIGHT NOW it's NOT

    p1=dt.date2num(dtStart)
    p2=dt.date2num(dtEnd)    
    mTimesStereo=np.array(TimesStereo)
    # get the right superDARN and EISCAT times
    SDtime=mTimesStereo[np.where(mTimesStereo>=p1)[0][0]:np.where(mTimesStereo<=p2)[0][-1]]

    SDdata1315=np.array(NeStereo1315[np.where(mTimesStereo>=p1)[0][0]:np.where(mTimesStereo<=p2)[0][-1]])
    SDdata1316=np.array(NeStereo1316[np.where(mTimesStereo>=p1)[0][0]:np.where(mTimesStereo<=p2)[0][-1]])
    SDdata1516=np.array(NeStereo1516[np.where(mTimesStereo>=p1)[0][0]:np.where(mTimesStereo<=p2)[0][-1]])

    df1315=pd.DataFrame(data=SDdata1315, index=SDtime)
    df1316=pd.DataFrame(data=SDdata1316, index=SDtime)
    df1516=pd.DataFrame(data=SDdata1516, index=SDtime)
    
    redf1315=np.array(df2.reindex(df1315.index, method='nearest'))[:,0]
    redf1316=np.array(df2.reindex(df1316.index, method='nearest'))[:,0]
    redf1516=np.array(df2.reindex(df1516.index, method='nearest'))[:,0]

    normalize1315=np.array((SDdata1315)/(1.0*redf1315))
    normalize1316=np.array((SDdata1316)/(1.0*redf1316))
    normalize1516=np.array((SDdata1516)/(1.0*redf1516))

    # final arrays
    cDataStereo=[NeStereo1315, NeStereo1316, NeStereo1516]
    maxDataStereo=[maxNeStereo1315, maxNeStereo1316, maxNeStereo1516]
    cTimeStereo=[mTimesStereo]*3
    dDataStereo=[normalize1315,normalize1316,normalize1516]
    dTimeStereo=[SDtime]*3
    ssData=[mNess13,mNess15,mNess16]
    ssTime=[rng2.to_pydatetime()]*3
    legendLabels=['13 MHz','15 MHz', '16 MHz']
    legendLabels2=['13 - 15 MHz',  '13 - 16 MHz', '15 - 16 MHz']
    legendLabels3=legendLabels2 + legendLabels
#
#
if date=='20160304':
    # first experiment
    #
    # assume range gate 11 for SuperDARN
    NeStereo1819=np.abs(dataStereo[0][4][5]) # pick the 15 to 16 MHz comparison
 
    Ness18=np.abs(datass[0][4])
    Ness19=np.abs(datass[0][5])

    TimesStereo=dt.date2num(dataStereo[2].to_pydatetime())
    Times18ss=datass[2][4]
    Times19ss=datass[2][5]
    #
    # resample the superDARN small shift data
    df18=pd.DataFrame({'18':Ness18}, index=pd.DatetimeIndex(Times18ss[:-1]))
    df19=pd.DataFrame({'19':Ness19}, index=pd.DatetimeIndex(Times19ss[:-1]))

    # resample to 2 minutes
    mNess18=np.array(df18['18'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))
    mNess19=np.array(df19['19'].resample('2T', how='mean').reindex(index=rng2, fill_value=np.nan))

    p1=dt.date2num(dtStart)
    p2=dt.date2num(dtEnd)    
    mTimesStereo=np.array(TimesStereo)
    # get the right superDARN and EISCAT times
    SDtime=mTimesStereo[np.where(mTimesStereo>=p1)[0][0]:np.where(mTimesStereo<=p2)[0][-1]]
    #Eiscattime=mTimesEiscat[np.where(mTimesEiscat>=p1)[0][0]:np.where(mTimesEiscat<=p2)[0][-1]]
    SDdata1819=np.array(NeStereo1819[np.where(mTimesStereo>=p1)[0][0]:np.where(mTimesStereo<=p2)[0][-1]])


    df1819=pd.DataFrame(data=SDdata1819, index=SDtime)
    
    #df2=pd.DataFrame(data=NeEiscat,index=Etime)
    redf1819=np.array(df2.reindex(df1819.index, method='nearest'))[:,0]

    normalize1819=np.array((SDdata1819)/(1.0*redf1819))


    # final arrays
    cDataStereo=[NeStereo1819]
    cTimeStereo=[mTimesStereo]
    dDataStereo=[normalize1819]
    dTimeStereo=[SDtime]
    ssData=[mNess18,mNess19]
    ssTime=[rng2.to_pydatetime()]*2
    legendLabels=['18 MHz', '19 MHz']
    legendLabels2=['18 - 19 MHz']
    legendLabels3=legendLabels2 + legendLabels
#
# resample EISCAT here
dfei=pd.DataFrame({'data':NeEiscat}, index=dt.num2date(Etime))
eiscatdat=dfei.resample('2T', how='mean')
eiscatdatIRI=np.array(dfei.resample('50T', how='mean'))
eiscattime=dt.date2num(eiscatdat.index.to_pydatetime())

#
# get correlation coefficient
#for icor in range(len(legendLabels2)):
#    index= cDataStereo[icor][iSt2:eSt2]
#    print 'Correlation coefficient for ' + str(legendLabels2[icor]) +' : ' +str( stats.linregress(np.array(eiscatdat)[:,0][~np.isnan(index)], index[~np.isnan(index)]))
#try:
#    print 'Correlation coefficient for IRI ' +' : ' +str( stats.linregress(eiscatdatIRI[:,0], NeIRI))
#except:
#    print 'Correlation coefficient for IRI ' +' : ' +str( stats.linregress(eiscatdatIRI[:,0], NeIRI[1:-1]))    

# NOW PLOT
CP.correctedCompare(cDataStereo,maxDataStereo,cTimeStereo,eiscatdat,eiscattime,NeIRI,IRItime,[ 'IRI']+legendLabels2 + ['EISCAT'],dtStart, dtEnd, date, ylimlow,ylimhigh)
CP.differencePlot(dDataStereo,dTimeStereo,Etime,legendLabels2,dtStart, dtEnd,date)
CP.shiftStereo(cDataStereo,cTimeStereo, ssData, ssTime,legendLabels2+legendLabels,dtStart, dtEnd, date, ylimlow2,ylimhigh2)

#
# plot EISCAT
import matplotlib as mpl
avgNe=np.nanmean(neE, axis=1)
fig=plt.figure()
ax2=fig.add_subplot(111)
plt.subplots_adjust(right=0.80, top=0.92, bottom=0.28, left=0.11)
ax2.set_ylabel('Altitude [km]', fontsize=22, fontweight='bold')
ax2.set_xlabel('$N_e$ [m$^{-3}$]', fontsize=22, fontweight='bold')
ax2.set_xscale('log')
ax2.set_ylim(80,500)
# resample over a 1 minute period
ax2.plot(avgNe,altE,lw=2, c='DarkOrange')
ax2.plot(IRIprof, IRIprofalts, lw=2, c='b')
ax2.set_xlim(1e10,1e12)
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['xtick.minor.width'] = 3
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['ytick.minor.size'] = 5
mpl.rcParams['ytick.minor.width'] = 3
plt.legend(['EISCAT', 'IRI'], bbox_to_anchor=[1.2, 0.9])
font = {'family' : 'normal',
                 'weight' : 'bold',
                 'size'   : 22}
plt.rc('font', **font)
fig.set_size_inches(13,9)
plt.savefig(date+'_EISCATprofile.png')
plt.close()

# number of points
#  2016-03-03
#  Down Sampled to 2 minutes 
#  13-15: 39
#  13-16: 22
#  15-16: 32

#  All of the Data 
#  13: 1089
#  15: 1562
#  16: 1710
#  
#  2015-03-12
#  Down Sampled to 2 minutes 
#  15-16: 23
#  16-17: 7
#  16-18: 14
#  17-18: 6

#  All of the Data 
#  15: 1024
#  16: 1365
#  17: 185
#  18: 168
