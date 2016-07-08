# readCUTLASSvelSTEREO.py
#
# our goal is to read in the CUTLASS-Hankasalmi SuperDARN data
# for 3-12-2015 to start with. This data was provided by Mike Kosch.
# File format is a little weird, but the data columns are height
# not sure of the exact height profile right now. Working on that...
# STEREO feature accounts for the 1 MHz dual freq mode
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
import readVelSTEREO45km

def sciNot(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)
#
           
#
# plot this
# initialize the figure
data = readVelSTEREO45km.readSUPERDARN()
Ne=data[0]
Freqs=data[1]
fp=data[2]
Times=data[3]

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
hours = HourLocator(interval=1) 
hours2 = MinuteLocator(interval=10) 
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
time=dt.date2num(Times)[::40][1:-1]

#Altitudes=480+np.array(range(75))*15
Altitudes=range(25)
X,Y=np.meshgrid(time, Altitudes)
data1=ma.masked_invalid(Ne[0][:-1]).transpose()
data2=ma.masked_invalid(Ne[1][:-1]).transpose()
data3=ma.masked_invalid(Ne[2][:-1]).transpose()

vmin=1*1e12
vmax=4*1e12
#ax.set_ylim(8,13)
#ax.set_xlim(time[int(3*len(time)/5)], time[-1])
ax.set_xlabel("Time", fontsize=20, fontweight='bold')
#ax.scatter(hope_epoch[start_hope:end_hope], PA_labels,    , marker='|', s=40
col=ax.pcolormesh(X,Y,data1, cmap='viridis', vmin=vmin, vmax=vmax)
col=ax.pcolormesh(X,Y,data2, cmap='viridis', vmin=vmin, vmax=vmax)
col=ax.pcolormesh(X,Y,data3, cmap='viridis', vmin=vmin, vmax=vmax)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
cbaxes = fig.add_axes([0.8, 0.27, 0.03, 0.65])
cb = plt.colorbar(col, cax = cbaxes,format=FuncFormatter(sciNot),ticks=np.linspace(vmin-.2, vmax+.1,6)) 
cb.set_label('Electron Number Density [m$^{-3}$]', fontsize=25)
fig.set_size_inches(13,9)
subdir_name='CUTLASS_Spectrogram'
if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name) 
os.chdir(subdir_name)#
plt.savefig('Stereo_Allgates_20150312Ne.png')
os.chdir('..')
plt.close()


