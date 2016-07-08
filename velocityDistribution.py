# velocityDistribution.py
#
# show the distribution of velocities for each frequency band to check for
# randomoness
#
# LKS, February 2016 Svalbard!
#
#
import numpy as np # always should load numpy
import datetime # nice for keeping track of time
import matplotlib.pyplot as plt # best package ever
import os # to switch directories
# need these for time intervals later
import matplotlib.dates as dt
from matplotlib.dates import HourLocator, DayLocator, DateFormatter, MinuteLocator
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter,FuncFormatter
import readCUTLASSvel
#
# 
data = readCUTLASSvel.readSUPERDARN()
AllNe=data[0]
AllFreqs=data[1]
Allfp=data[2]
AllTimes=data[3]
AllVel=data[5]
#
labels=['15', '16', '17', '18']
#
#
for ilabel in range(len(labels)):
    tempV=AllVel[ilabel]
    tempV[tempV==10000]=np.nan
    cV=tempV[~np.isnan(tempV)]
    fig=plt.figure()
    ax=fig.add_subplot(111)
    plt.subplots_adjust(right=0.88, top=0.92, bottom=0.28, left=0.18)
    font = {'family' : 'normal',
        'weight' : 'bold',
            'size'   : 20}
    plt.rc('font', **font)
    binwidth=50
    ax.set_xlim(-1000,1000)

    ax.hist(cV, bins=np.arange(min(cV), max(cV) + binwidth, binwidth))
    ax.set_ylabel('Number of Points') 
    ax.set_xlabel('Velocities')
    subdir_name='VelocityDistribution'
    if not os.path.exists(subdir_name):
        os.umask(0) # unmask if necessary
        os.makedirs(subdir_name) 
    os.chdir(subdir_name)#
    plt.savefig(labels[ilabel]+'_velDist.png')
    os.chdir('..')
    plt.close()    
