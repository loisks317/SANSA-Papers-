# calculateNs.py
#
# calculate the number density as a function of frequency
# make a guess at height
# may take a few times to iterate to convergence
#
# LKS, SANSA February 2016
#
#
import numpy as np # always should load numpy
import datetime # nice for keeping track of time
import matplotlib.pyplot as plt # best package ever
import os # to switch directories
# need these for time intervals later
import matplotlib.dates as dt
from matplotlib.dates import HourLocator, DayLocator, DateFormatter, MinuteLocator
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import readCUTLASSvel as readVel


#
# read in the data
data=readVel.readSUPERDARN()

AllTime=data[3]
AllNs=data[4]
#
# calculate elevation angle now
#
h=150 # km, this is a guess
Re=6374 # km, earth radius
inclination=78.3887 # for 150 km

for j in range(28,36):
  elvAngle={}
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
  ax.set_ylabel('SuperDarn Elevation Angle', fontsize=22, fontweight='bold')
  ax.set_xlabel('Time', fontsize=22, fontweight='bold')
  ax.set_ylim(0,90)
  ax.set_xlim(dt.date2num(np.min(AllTime[0])), dt.date2num(np.max(AllTime[-1])))
  colors=['red', 'green', 'blue', 'orange']
  p1=[]
  for i in range(4):
    NsTemp=AllNs[i][:,j]
    NsTemp[NsTemp>1]=np.nan
    #NsT=np.nanmedian(NsTemp,axis=1)
    #
    elvAngle[i]=90-np.arccos(NsTemp*np.sin(inclination*np.pi/180.)*(Re+h*1.0)/(1.0*Re))*180/(1.0*np.pi)
    p1.append(ax.scatter(dt.date2num(AllTime[i]),elvAngle[i], color=colors[i]))
  plt.legend(p1,['15 MHz', '16 MHz', '17 MHz', '18 MHz'], bbox_to_anchor=[1.2, 0.7])
  fig.set_size_inches(13,9)
  subdir_name='ElevationAngleEstimate'
  if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name) 
  os.chdir(subdir_name)#
  plt.savefig('ElevationAngleEstimates150km_rangeGate='+str(j)+'.png')
  plt.close()
  os.chdir('..')
