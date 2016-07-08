# calculateEISCATConductivity.py
#
# calculates the conductivities out of the EISCAT data at each height
# in the data files
# uses equations from schunk and nagy 2009 along with equations
# from kyoto dst site for hall and pedersen conductivity formulas.
#
# LKS, January 2016
#
# imports
import numpy as np
import os
import readProcessedEISCAT
import matplotlib.pyplot as plt # best package ever
# need these for time intervals later
import matplotlib.dates as dt
from matplotlib.dates import HourLocator, DayLocator, DateFormatter, MinuteLocator
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter 
#
#
# constants
e= 1.6*1e-19
me=9.11*1e-31
mi=1.67*1e-27
#
# now read in electron temperature, heights, electron ne, and ion collisions
# from processed guisdap data
dataEISCAT=readProcessedEISCAT.ReadData()
alt=np.nanmean(dataEISCAT['alt'], axis=1) # 42 altitudes between 76 km to 613 km
Times=dataEISCAT['stim']
ne=np.swapaxes(dataEISCAT['ne'],1,0)
ti=np.swapaxes(dataEISCAT['ti'],1,0)
te=np.swapaxes(dataEISCAT['te'],1,0)
vin=np.swapaxes(dataEISCAT['vin'],1,0)
vi=np.swapaxes(dataEISCAT['vel'],1,0)
os.chdir('ProcessedEISCAT')
B=np.genfromtxt('EISCATB.txt')
os.chdir('..')

os.chdir('ModelData')
#
# The order of MSIS columns
# height between 100 km to 800 km with 5 km steps
# output params: height, O, N2, O2, T, He, H, N
# where molecular species are in cm^-3 and temperature in K
#
MSISData=np.genfromtxt('MSISneutralDensity.txt')
os.chdir('..')
#
# now find nearest MSIS neutral densities at EISCAT heights
MSISheights=MSISData[:,0]
def find_nearest(array,value):
    idx = min(enumerate(array), key=lambda x: np.abs(x[1]-value))[0]
    return idx
#
#
# now combine this with shunk and nagy equations to get electron neutral
# collision frequencies

# weird, it doesn't depend on ion density at all... that's odd

#
# now create an altitude loop
vei=np.zeros((len(Times), len(alt)))
vee=np.zeros((len(Times), len(alt)))
sigP=np.zeros((len(Times), len(alt)))
sigH=np.zeros((len(Times), len(alt)))
cowling=np.zeros((len(Times), len(alt)))
for h in range(len(alt)):
  for itime in range(len(ne)):
    # equations for election-ion collisions and electron-electron collisions
    # are eqns 4.114 and 4.115 from schunk and nagy
    # assume quasi-neutrality, that ni = ne
    # assume charge state is 1

    vei[itime][h]=54.5 *1.0e-6* ne[itime][h] / (te[itime][h]**(3/2.))
    vee[itime][h]= 54.5 * 1.0e-6 *ne[itime][h]/(np.sqrt(2)*(te[itime][h]**(3/2.)))
    #ve=vei[itime,h]+vee[itime,h]
    #
    # calculate electron neutral collision frequencies
    # first get nearest height in MSIS
    idx=find_nearest(MSISheights, alt[h])
    Ven_N2=(2.33*1e-11)*(MSISData[idx, 2])*(1-1.21*1e-4 * te[itime][h])*te[itime][h]
    Ven_O2=(1.82*1e-10)*(MSISData[idx,3])*(1+3.6*1e-2 * np.sqrt(te[itime][h]))*np.sqrt(te[itime][h])
    Ven_O=(8.9*1e-11)*(MSISData[idx,1])*(1+5.7*1e-4 *te[itime][h])*np.sqrt(te[itime][h])
    Ven_He=(4.6*1e-10)*(MSISData[idx,5])*(np.sqrt(te[itime][h]))
    Ven_H=(4.5*1e-9)*(MSISData[idx,6])*(1-1.35*1e-4 *te[itime][h])*np.sqrt(te[itime][h])
    Ven=Ven_N2+Ven_O2+Ven_O+Ven_He+Ven_H
    ve=vei[itime][h]+Ven
    #
    # now calculate conductivity
    # assume vei >> ven, not sure if this is fair
    sig0= ne[itime,h]*(e**2)/(me*ve)
    # get B from IGRF
    # assuming Tromso location... not sure if this is correct
    # IGRF wants it in degrees
    # from https://www.eiscat.se/groups/Documentation/BasicInfo/locations.html
    #phi=19
    #theta =69.5
    iB=B[h]*1e-9
    # cyclotron freq
    cycE=e*iB/me
    # ASK MIKE ABOUT CYCLOTRON FREQ
    cycI=e*iB/(mi)
    # short hand for bigger eqs
    X= (cycE*cycI)/(vin[itime,h] * ve)
    #
    # pedersen
    sigP[itime,h]=sig0*((1+X)*(ve**2))/(((1+X)**2) * ve**2 + cycE**2)
    sigH[itime,h]=sig0*cycE*ve/((1+X)**2 *(ve**2) + cycE**2)
    cowling[itime,h]=sigP[itime,h] + sigH[itime,h]**2 / sigP[itime,h]

# now plot these quantities
dataC=[sigP,sigH,cowling]
labels=['Pedersen Conductivity', 'Hall Conductivity', 'Cowling Conductivity']
vminC=[-7,-9,-7]
vmaxC=[-5,-7,-5]
#
# units are approximatley 10-3 Siemans/meter for conductivity
for iDat in range(len(dataC)):
# initialize the figure
   fig=plt.figure()
   from numpy import ma
   os.chdir('/Users/loisks/Desktop/Functions/')
   import colormaps as cmaps
   os.chdir('/Users/loisks/Documents/ResearchProjects/SANSAProject/')
   plt.register_cmap(name='viridis', cmap=cmaps.viridis) 
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
   plt.subplots_adjust(right=0.83, top=0.92, bottom=0.28, left=0.11)
   ax.set_ylabel('Altitude [km]', fontsize=22, fontweight='bold')
   time=dt.date2num(Times)
   Altitudes=alt
   ax.set_ylim(alt[0], alt[-1])
   #ax.set_xlim(time[0], time)# last 17 missing
   X,Y=np.meshgrid(time,Altitudes)
   data=ma.masked_invalid(np.log10(dataC[iDat])).transpose()
   vmin=vminC[iDat]
   vmax=vmaxC[iDat]
   ax.set_xlabel("Time", fontsize=20, fontweight='bold')
   col=ax.pcolormesh(X,Y,data, cmap='viridis', vmin=vmin, vmax=vmax)
   font = {'family' : 'normal',
           'weight' : 'bold',
           'size'   : 22}
   plt.rc('font', **font)
   cbaxes = fig.add_axes([0.87, 0.27, 0.03, 0.65]) 
   cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%0.2f}$'),ticks=np.linspace(vmin-.1, vmax+.1,6))
   cb.set_label(labels[iDat]+' [S/m]', fontsize=25)
   fig.set_size_inches(13,9)
   subdir_name='EISCAT_Spectrogram'
   if not os.path.exists(subdir_name):
       os.umask(0) # unmask if necessary
       os.makedirs(subdir_name) 
   os.chdir(subdir_name)#
   plt.savefig('20150312_'+labels[iDat]+'.png')
   os.chdir('..')
   plt.close()    
#
# plot ne
fig=plt.figure()
from numpy import ma
os.chdir('/Users/loisks/Desktop/Functions/')
import colormaps as cmaps
os.chdir('/Users/loisks/Documents/ResearchProjects/SANSAProject/')
plt.register_cmap(name='viridis', cmap=cmaps.viridis) 
days = DayLocator(interval=1) 
#hours = MinuteLocator(interval=1) 
#hours2 = MinuteLocator(interval=2)
hours = MinuteLocator(interval=30) 
hours2 = MinuteLocator(interval=60) 
daysFmt = DateFormatter('%H:%M')
fig.gca().xaxis.set_major_locator(hours)
fig.gca().xaxis.set_major_formatter(daysFmt)
fig.gca().xaxis.set_minor_locator(hours2)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
ax=fig.add_subplot(111)
plt.subplots_adjust(right=0.83, top=0.92, bottom=0.28, left=0.11)
ax.set_ylabel('Altitude [km]', fontsize=22, fontweight='bold')
time=dt.date2num(Times)
Altitudes=alt
#ax.set_ylim(200,500)
#ax.set_xlim(time[0], time[10])# last 17 missing
X,Y=np.meshgrid(time, Altitudes)
data=ma.masked_invalid(np.log10(ne)).transpose()
vmin=11
vmax=12
ax.set_xlabel("Time", fontsize=20, fontweight='bold')
   # ax.scatter(hope_epoch[start_hope:end_hope], PA_labels,    , marker='|', s=40
col=ax.pcolormesh(X,Y,data, cmap='viridis', vmin=vmin, vmax=vmax)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
cbaxes = fig.add_axes([0.87, 0.27, 0.03, 0.65]) 
cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%0.2f}$'),ticks=np.linspace(vmin-.1, vmax+.1,6)) 
cb.set_label('m$^{-3}$', fontsize=25)
fig.set_size_inches(13,9)
subdir_name='EISCAT_Spectrogram'
if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name) 
os.chdir(subdir_name)#
plt.savefig('20150312_ne.png')
os.chdir('..')
plt.close()    


#
# plot te, ti
#
# plot ne
fig=plt.figure()
from numpy import ma
os.chdir('/Users/loisks/Desktop/Functions/')
import colormaps as cmaps
os.chdir('/Users/loisks/Documents/ResearchProjects/SANSAProject/')
plt.register_cmap(name='viridis', cmap=cmaps.viridis) 
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
plt.subplots_adjust(right=0.83, top=0.92, bottom=0.28, left=0.11)
ax.set_ylabel('Altitude [km]', fontsize=22, fontweight='bold')
time=dt.date2num(Times)
Altitudes=alt
ax.set_ylim(100,500)
ax.set_xlim(time[0], time[-1])# last 17 missing
X,Y=np.meshgrid(time, Altitudes)
data=ma.masked_invalid(np.log10(te)).transpose()
vmin=3
vmax=4
ax.set_xlabel("Time", fontsize=20, fontweight='bold')
   # ax.scatter(hope_epoch[start_hope:end_hope], PA_labels,    , marker='|', s=40
col=ax.pcolormesh(X,Y,data, cmap='viridis', vmin=vmin, vmax=vmax)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
cbaxes = fig.add_axes([0.87, 0.27, 0.03, 0.65]) 
cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(vmin, vmax+1)) 
cb.set_label('Electron Temperature [K]', fontsize=25)
fig.set_size_inches(13,9)
subdir_name='EISCAT_Spectrogram'
if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name) 
os.chdir(subdir_name)#
plt.savefig('20150312_te.png')
os.chdir('..')
plt.close()

# plot ti

fig=plt.figure()
from numpy import ma
os.chdir('/Users/loisks/Desktop/Functions/')
import colormaps as cmaps
os.chdir('/Users/loisks/Documents/ResearchProjects/SANSAProject/')
plt.register_cmap(name='viridis', cmap=cmaps.viridis) 
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
plt.subplots_adjust(right=0.83, top=0.92, bottom=0.28, left=0.11)
ax.set_ylabel('Altitude [km]', fontsize=22, fontweight='bold')
time=dt.date2num(Times)
Altitudes=alt
ax.set_ylim(100,500)
ax.set_xlim(time[0], time[-1])# last 17 missing
X,Y=np.meshgrid(time, Altitudes)
data=ma.masked_invalid(np.log10(ti)).transpose()
vmin=2
vmax=4
ax.set_xlabel("Time", fontsize=20, fontweight='bold')
   # ax.scatter(hope_epoch[start_hope:end_hope], PA_labels,    , marker='|', s=40
col=ax.pcolormesh(X,Y,data, cmap='viridis', vmin=vmin, vmax=vmax)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
cbaxes = fig.add_axes([0.87, 0.27, 0.03, 0.65]) 
cb = plt.colorbar(col, cax = cbaxes,format=FormatStrFormatter('$10^{%d}$'),ticks=range(vmin, vmax+1)) 
cb.set_label('Ion Temperature [K]', fontsize=25)
fig.set_size_inches(13,9)
subdir_name='EISCAT_Spectrogram'
if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name) 
os.chdir(subdir_name)#
plt.savefig('20150312_ti.png')
os.chdir('..')
plt.close()

# plot vi

fig=plt.figure()
from numpy import ma
os.chdir('/Users/loisks/Desktop/Functions/')
import colormaps as cmaps
os.chdir('/Users/loisks/Documents/ResearchProjects/SANSAProject/')
plt.register_cmap(name='viridis', cmap=cmaps.viridis) 
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
plt.subplots_adjust(right=0.83, top=0.92, bottom=0.28, left=0.11)
ax.set_ylabel('Altitude [km]', fontsize=22, fontweight='bold')
time=dt.date2num(Times)
Altitudes=alt
ax.set_ylim(100,500)
ax.set_xlim(time[0], time[-1])# last 17 missing
X,Y=np.meshgrid(time, Altitudes)
data=ma.masked_invalid(vi).transpose()
vmin=-100
vmax=100
ax.set_xlabel("Time", fontsize=20, fontweight='bold')
   # ax.scatter(hope_epoch[start_hope:end_hope], PA_labels,    , marker='|', s=40
col=ax.pcolormesh(X,Y,data, cmap='viridis', vmin=vmin, vmax=vmax)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)
cbaxes = fig.add_axes([0.87, 0.27, 0.03, 0.65]) 
cb = plt.colorbar(col, cax = cbaxes) 
cb.set_label('Ion Velocity [m/s]', fontsize=25)
fig.set_size_inches(13,9)
subdir_name='EISCAT_Spectrogram'
if not os.path.exists(subdir_name):
    os.umask(0) # unmask if necessary
    os.makedirs(subdir_name) 
os.chdir(subdir_name)#
plt.savefig('20150312_vi.png')
os.chdir('..')
plt.close() 




