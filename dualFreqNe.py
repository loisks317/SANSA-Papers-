# dualFreqNe.py
#
# calculate Ne using the STEREO mode of the CUTLASS Radar
# compare beteween 15 and 16 MHz and then 16,17,and 18 MHz at the
# end of the experiment
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
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter,FuncFormatter
import readCUTLASSvel

def sciNot(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)
#
           
#
# plot this
# initialize the figure
data = readCUTLASSvel.readSUPERDARN()
AllNe=data[0]
AllFreqs=data[1]
Allfp=data[2]
AllTimes=data[3]

labels=['15_MHz', '16_MHz', '17_MHz', '18_MHz']


for i in range(len(data[0])):
     Ne=AllNe[i]
     Freqs=AllFreqs[i]
     fp=Allfp[i]
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
     ax.set_ylabel('Range Gate Number', fontsize=22, fontweight='bold')
     time=dt.date2num(Times)
     #Altitudes=480+np.array(range(75))*15
     Altitudes=range(75)
     X,Y=np.meshgrid(time, Altitudes)
     data=ma.masked_invalid(np.log10(Ne[:-1])).transpose()
     vmin=10.8
     vmax=11.1
     ax.set_ylim(25,35)
     #ax.set_xlim(time[int(3*len(time)/5)], time[-1])
     ax.set_xlabel("Time", fontsize=20, fontweight='bold')
     #ax.scatter(hope_epoch[start_hope:end_hope], PA_labels,    , marker='|', s=40
     col=ax.pcolormesh(X,Y,data, cmap='viridis', vmin=vmin, vmax=vmax)
     font = {'family' : 'normal',
             'weight' : 'bold',
             'size'   : 22}
     plt.rc('font', **font)
     cbaxes = fig.add_axes([0.84, 0.27, 0.03, 0.65])
     cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%0.2f}$'),ticks=np.linspace(vmin-.2, vmax+.1,5)) 
     cb.set_label('Electron Number Density [m$^{-3}$]', fontsize=25)
     fig.set_size_inches(13,9)
     subdir_name='CUTLASS_Spectrogram'
     if not os.path.exists(subdir_name):
         os.umask(0) # unmask if necessary
         os.makedirs(subdir_name) 
     os.chdir(subdir_name)#
     plt.savefig(labels[i]+'20150312Ne.png')
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
plt.subplots_adjust(right=0.75, top=0.92, bottom=0.28, left=0.11)
ax.set_ylabel('Range Gate Number', fontsize=22, fontweight='bold')
time0=dt.date2num(AllTimes[0])
time1=dt.date2num(AllTimes[1])
time2=dt.date2num(AllTimes[2])
time3=dt.date2num(AllTimes[3])

Altitudes=range(75)
X0,Y0=np.meshgrid(time0, Altitudes)
data0=ma.masked_invalid(AllNe[0][:-1]).transpose()
X1,Y1=np.meshgrid(time1, Altitudes)
data1=ma.masked_invalid(AllNe[1][:-1]).transpose()
X2,Y2=np.meshgrid(time2, Altitudes)
data2=ma.masked_invalid(AllNe[2][:-1]).transpose()
X3,Y3=np.meshgrid(time3, Altitudes)
data3=ma.masked_invalid(AllNe[3][:-1]).transpose()
vmin=7*1e10
vmax=1*1e11
ax.set_ylim(25,35)
ax.set_xlim(dt.date2num(datetime.datetime(2015,3,12,10,0,0)), dt.date2num(datetime.datetime(2015,3,12,12,0,0)))
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
cbaxes = fig.add_axes([0.8, 0.27, 0.03, 0.65])
cb = plt.colorbar(col0, cax = cbaxes,format=FuncFormatter(sciNot),ticks=np.linspace(vmin,vmax,5))
cb.set_label('Electron Number Density [m$^{-3}$]', fontsize=25)
fig.set_size_inches(13,9)
subdir_name='CUTLASS_Spectrogram'
if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name) 
os.chdir(subdir_name)#
plt.savefig('AllFreqs_20150312Ne.png')
os.chdir('..')
plt.close()
