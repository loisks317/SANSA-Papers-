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
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter,FuncFormatter
import readCUTLASSvel45km

def sciNot(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)
#
# may a good plan is to average all range gates together in 25-38 range or
# where we think the heater is and then do small freq shifts
# based on time (kHz scale shifts)
           
#
# plot this
# initialize the figure
data = readCUTLASSvel45km.readSUPERDARN()
AllNe=data[0]
AllFreqs=data[1]
Allfp=data[2]
AllTimes=data[3]

labels=['15_MHz', '16_MHz', '17_MHz', '18_MHz']


for i in range(len(data[0])):
     Ne=AllNe[i]
     Freqs=AllFreqs[i]
     fp=Allfp[i]
     Times=AllTimes
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
     plt.subplots_adjust(right=0.70, top=0.92, bottom=0.28, left=0.11)
     ax.set_ylabel('Range Gate Number', fontsize=22, fontweight='bold')
     time=dt.date2num(Times)[::40][1:-1]
     #Altitudes=480+np.array(range(75))*15
     Altitudes=range(25)
     X,Y=np.meshgrid(time, Altitudes)
     data=ma.masked_invalid(Ne[:-1]).transpose()
     vmin=3*1e12
     vmax=3.5*1e12
     ax.set_ylim(8,13)
     #ax.set_xlim(time[int(3*len(time)/5)], time[-1])
     ax.set_xlabel("Time", fontsize=20, fontweight='bold')
     #ax.scatter(hope_epoch[start_hope:end_hope], PA_labels,    , marker='|', s=40
     col=ax.pcolormesh(X,Y,data, cmap='viridis', vmin=vmin, vmax=vmax)
     font = {'family' : 'normal',
             'weight' : 'bold',
             'size'   : 22}
     plt.rc('font', **font)
     cbaxes = fig.add_axes([0.84, 0.27, 0.03, 0.65])
     cb = plt.colorbar(col, cax = cbaxes,format=FuncFormatter(sciNot),ticks=np.linspace(vmin,vmax,5))
     cb.set_label('Electron Number Density [m$^{-3}$]', fontsize=25)
     fig.set_size_inches(13,9)
     subdir_name='CUTLASS_Spectrogram'
     if not os.path.exists(subdir_name):
         os.umask(0) # unmask if necessary
         os.makedirs(subdir_name) 
     os.chdir(subdir_name)#
     plt.savefig(labels[i]+'20150312Ne_freqshifts.png')
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
plt.subplots_adjust(right=0.7, top=0.92, bottom=0.28, left=0.11)
ax.set_ylabel('Range Gate Number', fontsize=22, fontweight='bold')
#time0=dt.date2num(AllTimes[0])
#time1=dt.date2num(AllTimes[1])
#time2=dt.date2num(AllTimes[2])
#time3=dt.date2num(AllTimes[3])

Altitudes=range(25)
#X0,Y0=np.meshgrid(time, Altitudes)
data0=ma.masked_invalid(AllNe[0][:-1]).transpose()
data1=ma.masked_invalid(AllNe[1][:-1]).transpose()
data2=ma.masked_invalid(AllNe[2][:-1]).transpose()
data3=ma.masked_invalid(AllNe[3][:-1]).transpose()
vmin=3*1e12
vmax=4*1e12
ax.set_ylim(8,13)
ax.set_xlim(dt.date2num(datetime.datetime(2015,3,12,10,0,0)), dt.date2num(datetime.datetime(2015,3,12,12,0,0)))
#ax.set_xlim(time[int(3*len(time)/5)], time[-1])
ax.set_xlabel("Time", fontsize=20, fontweight='bold')
#ax.scatter(hope_epoch[start_hope:end_hope], PA_labels,    , marker='|', s=40
col0=ax.pcolormesh(X,Y,data0, cmap='viridis', vmin=vmin, vmax=vmax)
col1=ax.pcolormesh(X,Y,data1, cmap='viridis', vmin=vmin, vmax=vmax)
col2=ax.pcolormesh(X,Y,data2, cmap='viridis', vmin=vmin, vmax=vmax)
col3=ax.pcolormesh(X,Y,data3, cmap='viridis', vmin=vmin, vmax=vmax)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
cbaxes = fig.add_axes([0.8, 0.27, 0.03, 0.65])
cb = plt.colorbar(col0, cax = cbaxes,format=FuncFormatter(sciNot),ticks=np.linspace(vmin,vmax,5))
cb.set_label('Electron Number Density [m$^{-3}$]', fontsize=25)
fig.set_size_inches(13,9)
subdir_name='CUTLASS_Spectrogram'
if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name) 
os.chdir(subdir_name)#
plt.savefig('AllFreqs_20150312Ne_freqshifts.png')
os.chdir('..')
plt.close()
