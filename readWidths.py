# readwidths.py
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
import readCUTLASSwidths
import readCUTLASSvel
#
           
#
# plot this
# initialize the figure
data = readCUTLASSwidths.readSUPERDARN()
velocity=readCUTLASSvel.readSUPERDARN()
AllData=data[0]
AllFreqs=data[1]
AllTimes=data[2]
AllVel=velocity[-1] # should be velocity

labels=['15_MHz', '16_MHz', '17_MHz', '18_MHz']


for i in range(len(data[0])):
     Widths=AllData[i]
     Vel=AllVel[i]
     Freqs=AllFreqs[i]
     Times=AllTimes[i]
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
     Widths[Widths<0]=np.nan
     data=ma.masked_invalid(Widths)#.transpose()
     vmin=0
     vmax=50
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
     cb = plt.colorbar(col, cax = cbaxes) 
     cb.set_label('Spectral Width [m/s]', fontsize=25)
     fig.set_size_inches(13,9)
     subdir_name='CUTLASS_Spectrogram'
     if not os.path.exists(subdir_name):
         os.umask(0) # unmask if necessary
         os.makedirs(subdir_name) 
     os.chdir(subdir_name)#
     plt.savefig(labels[i]+'20150312_SpectralWidths.png')
     os.chdir('..')
     plt.close()



Widths0=AllData[0]
Freqs0=AllFreqs[0]
Widths1=AllData[1]
Freqs1=AllFreqs[1]
Widths2=AllData[2]
Freqs2=AllFreqs[2]
Widths3=AllData[3]
Freqs3=AllFreqs[3]
Times0=AllTimes[0]
Times1=AllTimes[1]
Times2=AllTimes[2]
Times3=AllTimes[3]
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
ax.set_ylabel('Range Gates', fontsize=22, fontweight='bold')
time0=dt.date2num(Times0)
time1=dt.date2num(Times1)
time2=dt.date2num(Times2)
time3=dt.date2num(Times3)
Altitudes=range(75)
X0,Y0=np.meshgrid(time0, Altitudes)
X1,Y1=np.meshgrid(time1, Altitudes)
X2,Y2=np.meshgrid(time2, Altitudes)
X3,Y3=np.meshgrid(time3, Altitudes)
Widths0[Widths0<0]=np.nan
Widths1[Widths1<0]=np.nan
Widths2[Widths2<0]=np.nan
Widths3[Widths3<0]=np.nan

data0=ma.masked_invalid(Widths0)#.transpose()
data1=ma.masked_invalid(Widths1)
data2=ma.masked_invalid(Widths2)
data3=ma.masked_invalid(Widths3)
vmin=0
vmax=50
ax.set_ylim(25,38)
ax.set_xlim(dt.date2num(datetime.datetime(2015,3,12,10,0,0)), datetime.datetime(2015,3,12,12,0,0))
#ax.set_xlim(time[int(3*len(time)/5)], time[-1])
ax.set_xlabel("Time", fontsize=20, fontweight='bold')
#ax.scatter(hope_epoch[start_hope:end_hope], PA_labels,    , marker='|', s=40
col=ax.pcolormesh(X0,Y0,data0, cmap='viridis', vmin=vmin, vmax=vmax)
col=ax.pcolormesh(X1,Y1,data1, cmap='viridis', vmin=vmin, vmax=vmax)
col=ax.pcolormesh(X2,Y2,data2, cmap='viridis', vmin=vmin, vmax=vmax)
col=ax.pcolormesh(X3,Y3,data3, cmap='viridis', vmin=vmin, vmax=vmax)

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
cbaxes = fig.add_axes([0.84, 0.27, 0.03, 0.65])
cb = plt.colorbar(col, cax = cbaxes) 
cb.set_label('Spectral Width [m/s]', fontsize=25)
fig.set_size_inches(13,9)
subdir_name='CUTLASS_Spectrogram'
if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name) 
os.chdir(subdir_name)#
plt.savefig('20150312_AllSpectralWidths.png')
os.chdir('..')
plt.close()
#
#
#
# spectrogram of velocities 
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
time0=dt.date2num(Times0)
time1=dt.date2num(Times1)
time2=dt.date2num(Times2)
time3=dt.date2num(Times3)
Altitudes=range(75)
X0,Y0=np.meshgrid(time0, Altitudes)
X1,Y1=np.meshgrid(time1, Altitudes)
X2,Y2=np.meshgrid(time2, Altitudes)
X3,Y3=np.meshgrid(time3, Altitudes)
#Widths0[Widths0<0]=np.nan
#Widths1[Widths1<0]=np.nan
#Widths2[Widths2<0]=np.nan
#Widths3[Widths3<0]=np.nan
AllVel[0][AllVel[0]==10000]=np.nan
AllVel[1][AllVel[1]==10000]=np.nan
AllVel[2][AllVel[2]==10000]=np.nan
AllVel[3][AllVel[3]==10000]=np.nan

data0=ma.masked_invalid(AllVel[0])#.transpose()
data1=ma.masked_invalid(AllVel[1])
data2=ma.masked_invalid(AllVel[2])
data3=ma.masked_invalid(AllVel[3])
vmin=-200
vmax=50
ax.set_ylim(25,38)
ax.set_xlim(dt.date2num(datetime.datetime(2015,3,12,10,0,0)), datetime.datetime(2015,3,12,12,0,0))
#ax.set_xlim(time[int(3*len(time)/5)], time[-1])
ax.set_xlabel("Time", fontsize=20, fontweight='bold')
#ax.scatter(hope_epoch[start_hope:end_hope], PA_labels,    , marker='|', s=40
col=ax.pcolormesh(X0,Y0,data0, cmap='viridis', vmin=vmin, vmax=vmax)
col=ax.pcolormesh(X1,Y1,data1, cmap='viridis', vmin=vmin, vmax=vmax)
col=ax.pcolormesh(X2,Y2,data2, cmap='viridis', vmin=vmin, vmax=vmax)
col=ax.pcolormesh(X3,Y3,data3, cmap='viridis', vmin=vmin, vmax=vmax)

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
cbaxes = fig.add_axes([0.84, 0.27, 0.03, 0.65])
cb = plt.colorbar(col, cax = cbaxes) 
cb.set_label('Line of Sight Velocities [m/s]', fontsize=25)
fig.set_size_inches(13,9)
subdir_name='CUTLASS_Spectrogram'
if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name) 
os.chdir(subdir_name)#
plt.savefig('20150312_AllLOSVel.png')
os.chdir('..')
plt.close()
