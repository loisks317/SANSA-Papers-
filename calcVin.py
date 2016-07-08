# calcVin.py
#
# our goal is to read in the CUTLASS-Hankasalmi SuperDARN data
# for 3-12-2015 to start with. This data was provided by Mike Kosch.
# File format is a little weird, but the data columns are height
# not sure of the exact height profile right now. Working on that...
# calculate the ion neutral collision frequency
#
#
# LKS February 2016, SANSA
# modified in February 2015, Svalbard to invlude the gradient term
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
import readCUTLASSvel 
#
           
#
# plot this
# initialize the figure
data = readCUTLASSvel.readSUPERDARN()
AllNe=data[0]
AllFreqs=data[1]
Allfp=data[2]
AllTimes=data[3]
Velocity=data[5]

labels=['15_MHz', '16_MHz', '17_MHz', '18_MHz']
vin15=np.ones((len(AllTimes[0]),75))*np.nan
vin16=np.ones((len(AllTimes[1]),75))*np.nan
vin17=np.ones((len(AllTimes[2]),75))*np.nan
vin18=np.ones((len(AllTimes[3]),75))*np.nan
nu15=np.ones((len(AllTimes[0]),75))*np.nan
nu16=np.ones((len(AllTimes[1]),75))*np.nan
nu17=np.ones((len(AllTimes[2]),75))*np.nan
nu18=np.ones((len(AllTimes[3]),75))*np.nan
nu=[nu15,nu16,nu17,nu18]
vin=[vin15,vin16,vin17,vin18]


for i in range(4):
 for iH in range(75-1):
     vel=Velocity[i][iH]
     vel2=Velocity[i][iH+1]
     Times=AllTimes[i]
     Freqs=AllFreqs[i]
     #
     # calculate the velocity derivative
     tLen=len(Times)
     m1=iR=0
     while iR in range(tLen):
         try:
             while vel[iR]==10000:
                   iR+=1
         except(IndexError):
             break
              # now we have a first value
         v1=vel[iR]
         t1=Times[iR]
         m1=iR # memory IR
         
         try:
             while vel[iR+1]==10000: # bad data 
                   iR+=1
         except(IndexError):
             break
         if iR < tLen-2:
             # break out of loop
            v2=vel[iR+1]
            vR2=vel2[iR]
            t2=Times[iR+1]
            diffv=v2-v1
            difft=(t2-t1)
            vin[i][m1][iH]=np.abs((diffv)/(((v2+v1)/2.0)*1.0*difft.seconds))
            #
            # have to add in the spatial gradient term
            if( vR2 == 10000) or (v1==10000):
               diffv2=np.nan
            else:
               diffv2=vR2-v1
            alt=15.0*1000
            #print(vin[i][m1][iH]/(np.abs(diffv2/(1.0*alt))))
            vin[i][m1][iH]=vin[i][m1][iH]+np.abs(diffv2/(1.0*alt))
            Ti=300
            CN2Oplus=6.82*1e10
            CO2Oplus=8.10*1e10
            CN2O2plus=4.13*1e10
            CO2N2plus=4.49*1e10
            O2O2plus=(2.59*1e11)*(1.0/np.sqrt(Ti))*(1/(1-0.073*np.log10(Ti))**2)
            N2N2plus=(5.14*1e11)*(1.0/np.sqrt(Ti))*(1/(1-0.069*np.log10(Ti))**2)
            nu[i][m1][iH]=(1.0e6)*vin[i][m1][iH]*0.5*(O2O2plus+N2N2plus+CN2Oplus+CO2Oplus+CN2O2plus+CO2N2plus)
            
            
            #nuN2[i][m1][iH]=(1.0e6)*(vin[i][m1][iH]/2.0)/ ((5.14*1e-11)*np.sqrt(Ti)*(1-0.069*np.log10(Ti))**2)
            iR+=1
         else:
          iR+=1

 #
 # calculated neutral density
 #
 

 
 fig=plt.figure()
 #
 from numpy import ma
 os.chdir('/Users/loisks/Desktop/Functions/')
 import colormaps as cmaps
 os.chdir('/Users/loisks/Documents/ResearchProjects/SANSAProject/CUTLASS/')
 # this is a colormap I like a lot but I have to load it extra special :)
 plt.register_cmap(name='viridis', cmap=cmaps.viridis) 
 fig=plt.figure()
 days = DayLocator(interval=1) 
 hours = MinuteLocator(interval=30) 
 hours2 = HourLocator(interval=1) 
 daysFmt = DateFormatter('%H:%M')
 fig.gca().xaxis.set_major_locator(hours)
 fig.gca().xaxis.set_major_formatter(daysFmt)
 fig.gca().xaxis.set_minor_locator(hours2)
 font = {'family' : 'normal',
         'weight' : 'bold',
         'size'   : 22}
 plt.rc('font', **font)
 ax=fig.add_subplot(111)
 plt.subplots_adjust(right=0.80, top=0.92, bottom=0.28, left=0.11)
 ax.set_ylabel('Approximate Distance Away [km]', fontsize=22, fontweight='bold')
 time=dt.date2num(Times)
 Altitudes=480+np.array(range(75))*15
 X,Y=np.meshgrid(time, Altitudes)
 data=ma.masked_invalid(np.log10(vin[i])).transpose()
 vmin=-3
 vmax=0
 ax.set_ylim(800,1100)
 #ax.set_xlim(time[int(3*len(time)/5)], time[-1])
 ax.set_xlabel("Time", fontsize=20, fontweight='bold')
 #ax.scatter(hope_epoch[start_hope:end_hope], PA_labels,    , marker='|', s=40
 col=ax.pcolormesh(X,Y,data, cmap='viridis', vmin=vmin, vmax=vmax)
 font = {'family' : 'normal',
         'weight' : 'bold',
         'size'   : 22}
 plt.rc('font', **font)
 cbaxes = fig.add_axes([0.84, 0.27, 0.03, 0.65])
 cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(vmin, vmax+1)) 
 cb.set_label('Collision Frequency [Hz]', fontsize=25)
 fig.set_size_inches(13,9)
 subdir_name='CUTLASS_Spectrogram'
 if not os.path.exists(subdir_name):
     os.umask(0) # unmask if necessary
     os.makedirs(subdir_name) 
 os.chdir(subdir_name)#
 plt.savefig(labels[i]+'20150312_CollisionFreq.png')
 os.chdir('..')
 plt.close()



#
#
# neutral density
 fig=plt.figure()
 #
 from numpy import ma
 os.chdir('/Users/loisks/Desktop/Functions/')
 import colormaps as cmaps
 os.chdir('/Users/loisks/Documents/ResearchProjects/SANSAProject/CUTLASS/')
 # this is a colormap I like a lot but I have to load it extra special :)
 plt.register_cmap(name='viridis', cmap=cmaps.viridis) 
 fig=plt.figure()
 days = DayLocator(interval=1) 
 hours = MinuteLocator(interval=30) 
 hours2 = HourLocator(interval=1) 
 daysFmt = DateFormatter('%H:%M')
 fig.gca().xaxis.set_major_locator(hours)
 fig.gca().xaxis.set_major_formatter(daysFmt)
 fig.gca().xaxis.set_minor_locator(hours2)
 font = {'family' : 'normal',
         'weight' : 'bold',
         'size'   : 22}
 plt.rc('font', **font)
 ax=fig.add_subplot(111)
 plt.subplots_adjust(right=0.80, top=0.92, bottom=0.28, left=0.11)
 ax.set_ylabel('Approximate Distance Away [km]', fontsize=22, fontweight='bold')
 time=dt.date2num(Times)
 Altitudes=480+np.array(range(75))*15
 X,Y=np.meshgrid(time, Altitudes)
 data=ma.masked_invalid(np.log10(nu[i])).transpose()
 vmin=15
 vmax=17
 ax.set_ylim(800,1100)
 #ax.set_xlim(time[int(3*len(time)/5)], time[-1])
 ax.set_xlabel("Time", fontsize=20, fontweight='bold')
 #ax.scatter(hope_epoch[start_hope:end_hope], PA_labels,    , marker='|', s=40
 col=ax.pcolormesh(X,Y,data, cmap='viridis', vmin=vmin, vmax=vmax)
 font = {'family' : 'normal',
         'weight' : 'bold',
         'size'   : 22}
 plt.rc('font', **font)
 cbaxes = fig.add_axes([0.84, 0.27, 0.03, 0.65])
 cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(vmin, vmax+1)) 
 cb.set_label('m$^{-3}$', fontsize=25)
 fig.set_size_inches(13,9)
 subdir_name='CUTLASS_Spectrogram'
 if not os.path.exists(subdir_name):
     os.umask(0) # unmask if necessary
     os.makedirs(subdir_name) 
 os.chdir(subdir_name)#
 plt.savefig(labels[i]+'20150312_NeutralDensity.png')
 os.chdir('..')
 plt.close()


fig=plt.figure()
#
from numpy import ma
os.chdir('/Users/loisks/Desktop/Functions/')
import colormaps as cmaps
os.chdir('/Users/loisks/Documents/ResearchProjects/SANSAProject/CUTLASS/')
# this is a colormap I like a lot but I have to load it extra special :)
plt.register_cmap(name='viridis', cmap=cmaps.viridis) 
fig=plt.figure()
days = DayLocator(interval=1) 
hours = MinuteLocator(interval=30) 
hours2 = HourLocator(interval=1) 
daysFmt = DateFormatter('%H:%M')
fig.gca().xaxis.set_major_locator(hours)
fig.gca().xaxis.set_major_formatter(daysFmt)
fig.gca().xaxis.set_minor_locator(hours2)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
ax=fig.add_subplot(111)
plt.subplots_adjust(right=0.80, top=0.92, bottom=0.28, left=0.11)
ax.set_ylabel('Range Gates', fontsize=22, fontweight='bold')
time0=dt.date2num(AllTimes[0])
time1=dt.date2num(AllTimes[1])
time2=dt.date2num(AllTimes[2])
time3=dt.date2num(AllTimes[3])

#Altitudes=480+np.array(range(75))*15
Altitudes=np.array(range(75))
X0,Y0=np.meshgrid(time0, Altitudes)
data0=ma.masked_invalid(np.log10(nu[0])).transpose()
X1,Y1=np.meshgrid(time1, Altitudes)
data1=ma.masked_invalid(np.log10(nu[1])).transpose()
X2,Y2=np.meshgrid(time2, Altitudes)
data2=ma.masked_invalid(np.log10(nu[2])).transpose()
X3,Y3=np.meshgrid(time3, Altitudes)
data3=ma.masked_invalid(np.log10(nu[3])).transpose()
vmin=15
vmax=17
ax.set_ylim(25,35)
#ax.set_xlim(time[int(3*len(time)/5)], time[-1])
ax.set_xlabel("Time", fontsize=20, fontweight='bold')
#ax.scatter(hope_epoch[start_hope:end_hope], PA_labels,    , marker='|', s=40
col0=ax.pcolormesh(X0,Y0,data0, cmap='viridis', vmin=vmin, vmax=vmax)
col1=ax.pcolormesh(X1,Y1,data1, cmap='viridis', vmin=vmin, vmax=vmax)
col2=ax.pcolormesh(X2,Y2,data2, cmap='viridis', vmin=vmin, vmax=vmax)
col3=ax.pcolormesh(X3,Y3,data3, cmap='viridis', vmin=vmin, vmax=vmax)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
cbaxes = fig.add_axes([0.84, 0.27, 0.03, 0.65])
cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(vmin, vmax+1) )
cb.set_label('m$^{-3}$', fontsize=25)
fig.set_size_inches(13,9)
subdir_name='CUTLASS_Spectrogram'
if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name) 
os.chdir(subdir_name)#
plt.savefig('AllFreqs_20150312_NeutralDensity.png')
os.chdir('..')
plt.close()

fig=plt.figure()
#
from numpy import ma
os.chdir('/Users/loisks/Desktop/Functions/')
import colormaps as cmaps
os.chdir('/Users/loisks/Documents/ResearchProjects/SANSAProject/CUTLASS/')
# this is a colormap I like a lot but I have to load it extra special :)
plt.register_cmap(name='viridis', cmap=cmaps.viridis) 
fig=plt.figure()
days = DayLocator(interval=1) 
hours = MinuteLocator(interval=30) 
hours2 = HourLocator(interval=1) 
daysFmt = DateFormatter('%H:%M')
fig.gca().xaxis.set_major_locator(hours)
fig.gca().xaxis.set_major_formatter(daysFmt)
fig.gca().xaxis.set_minor_locator(hours2)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
ax=fig.add_subplot(111)
plt.subplots_adjust(right=0.80, top=0.92, bottom=0.28, left=0.11)
ax.set_ylabel('Approximate Distance Away [km]', fontsize=22, fontweight='bold')
time0=dt.date2num(AllTimes[0])
time1=dt.date2num(AllTimes[1])
time2=dt.date2num(AllTimes[2])
time3=dt.date2num(AllTimes[3])

Altitudes=480+np.array(range(75))*15
X0,Y0=np.meshgrid(time0, Altitudes)
data0=ma.masked_invalid(np.log10(vin[0])).transpose()
X1,Y1=np.meshgrid(time1, Altitudes)
data1=ma.masked_invalid(np.log10(vin[1])).transpose()
X2,Y2=np.meshgrid(time2, Altitudes)
data2=ma.masked_invalid(np.log10(vin[2])).transpose()
X3,Y3=np.meshgrid(time3, Altitudes)
data3=ma.masked_invalid(np.log10(vin[3])).transpose()
vmin=-3
vmax=0
ax.set_ylim(800,1100)
#ax.set_xlim(time[int(3*len(time)/5)], time[-1])
ax.set_xlabel("Time", fontsize=20, fontweight='bold')
#ax.scatter(hope_epoch[start_hope:end_hope], PA_labels,    , marker='|', s=40
col0=ax.pcolormesh(X0,Y0,data0, cmap='viridis', vmin=vmin, vmax=vmax)
col1=ax.pcolormesh(X1,Y1,data1, cmap='viridis', vmin=vmin, vmax=vmax)
col2=ax.pcolormesh(X2,Y2,data2, cmap='viridis', vmin=vmin, vmax=vmax)
col3=ax.pcolormesh(X3,Y3,data3, cmap='viridis', vmin=vmin, vmax=vmax)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
cbaxes = fig.add_axes([0.84, 0.27, 0.03, 0.65])
cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(vmin, vmax+1) )
cb.set_label('Collision Frequency [Hz]', fontsize=25)
fig.set_size_inches(13,9)
subdir_name='CUTLASS_Spectrogram'
if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name) 
os.chdir(subdir_name)#
plt.savefig('AllFreqs_20150312_CollFreq.png')
os.chdir('..')
plt.close()
