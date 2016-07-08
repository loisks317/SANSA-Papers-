# calculateSuperDARNconductivity.py
#
# calculate the SuperDARN conductivity using height estimates from the EISCAT
# superDARn comparison plot (guess at 140 km for now based on ray tracing)
# use MSIS for neutral density collision frequencies
# use the EISCAT data for ion temperature. It would be stupid to use
# IRI when EISCAT is available for 140-150 km
#
# LKS, plane from CPT -> DOH,  February 20, 2016
#
import numpy as np
import os
import datetime
import matplotlib.pyplot as plt # best package ever
# need these for time intervals later
import matplotlib.dates as dt
from matplotlib.dates import HourLocator, DayLocator, DateFormatter, MinuteLocator
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter,FuncFormatter
import readCUTLASSvel
os.chdir('..')
import readProcessedEISCAT
#
# defined functions
#
def sciNot(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)
#
# get the EISCAT data
dataEISCAT=readProcessedEISCAT.ReadData()
alt=np.nanmean(dataEISCAT['alt'], axis=1) # 42 altitudes between 76 km to 613 km
# find the nearest altitude index
desiredAlt=140
hindex=np.nanargmin(np.abs(alt-desiredAlt))
print(alt[hindex])
Times=dataEISCAT['stim']
ti=np.swapaxes(dataEISCAT['ti'],1,0)[:,hindex]
te=np.swapaxes(dataEISCAT['te'],1,0)[:,hindex]
vin=np.swapaxes(dataEISCAT['vin'],1,0)[:,hindex]
vi=np.swapaxes(dataEISCAT['vel'],1,0)[:,hindex]
os.chdir('CUTLASS')
#
#
# constants
e= 1.6*1e-19
me=9.11*1e-31
mi=1.67*1e-27
#
#
# plot this
# initialize the figure
# just do 16 MHz
data = readCUTLASSvel.readSUPERDARN()
AllNe=data[0][1]
AllFreqs=data[1]
Allfp=data[2][1]
AllTimes=data[3][1]
gates=75 # just how it is
labels=['15_MHz', '16_MHz', '17_MHz', '18_MHz']
#
# put the EISCAT data on SuperDARN time scales
# already at the right heights

EISti=np.zeros(len(AllTimes))
EISte=np.zeros(len(AllTimes))
EISvin=np.zeros(len(AllTimes))
EISvi=np.zeros(len(AllTimes))
#
# Tnos
TnosSD=dt.date2num(AllTimes)
TnosEIS=dt.date2num(Times)
for iT in range(len(AllTimes)):
    # get the SuperDARN time
    Rindex=np.nanargmin(np.abs(TnosSD[iT]-TnosEIS))
    EISti[iT]=ti[Rindex]
    EISte[iT]=te[Rindex]
    EISvin[iT]=vin[Rindex]
    EISvi[iT]=vi[Rindex]    
#
# Get the B Data then MSIS
#
os.chdir('..')
os.chdir('ProcessedEISCAT')
B=np.genfromtxt('EISCATB.txt')[hindex]
os.chdir('..')
#
# The order of MSIS columns
# height between 100 km to 800 km with 5 km steps
# output params: height, O, N2, O2, T, He, H, N
# where molecular species are in cm^-3 and temperature in K
#
os.chdir('ModelData')
MSISData=np.genfromtxt('MSISneutralDensity.txt')
os.chdir('..')
#
# now find nearest MSIS neutral densities at EISCAT heights
MSISheights=MSISData[:,0]
vei=np.zeros((len(AllTimes), gates))
vee=np.zeros((len(AllTimes), gates))
sigP=np.zeros((len(AllTimes),gates))
sigH=np.zeros((len(AllTimes), gates))
for iH in range(gates): # number of range gates 
    for itime in range(len(AllTimes)):
    # equations for election-ion collisions and electron-electron collisions
    # are eqns 4.114 and 4.115 from schunk and nagy
    # assume quasi-neutrality, that ni = ne
    # assume charge state is 1
            vei[itime][iH]=54.5 *1.0e-6* AllNe[itime][iH] / (EISte[itime]**(3/2.))
           # vee[itime][iH]= 54.5 * 1.0e-6 *AllNe[itime][iH]/(np.sqrt(2)*(EISte[itime]**(3/2.)))
            # calculate electron neutral collision frequencies
            # first get nearest height in MSIS
            idx=np.nanargmin(np.abs(MSISheights-desiredAlt) )
            Ven_N2=(2.33*1e-11)*(MSISData[idx, 2])*(1-1.21*1e-4 * EISte[itime])*EISte[itime]
            Ven_O2=(1.82*1e-10)*(MSISData[idx,3])*(1+3.6*1e-2 * np.sqrt(EISte[itime]))*np.sqrt(EISte[itime])
            Ven_O=(8.9*1e-11)*(MSISData[idx,1])*(1+5.7*1e-4 *EISte[itime])*np.sqrt(EISte[itime])
            Ven_He=(4.6*1e-10)*(MSISData[idx,5])*(np.sqrt(EISte[itime]))
            Ven_H=(4.5*1e-9)*(MSISData[idx,6])*(1-1.35*1e-4 *EISte[itime])*np.sqrt(EISte[itime])
            Ven=Ven_N2+Ven_O2+Ven_O+Ven_He+Ven_H
            ve=vei[itime][iH]+Ven
            # now calculate conductivity
            # assume vei >> ven, not sure if this is fair
            sig0=AllNe[itime,iH]*(e**2)/(me*ve)
            # get B from IGRF
            # assuming Tromso location... not sure if this is correct
            # IGRF wants it in degrees
            # from https://www.eiscat.se/groups/Documentation/BasicInfo/locations.html
            #phi=19
            #theta =69.5
            iB=B*1e-9
            # cyclotron freq
            cycE=e*iB/me
            # ASK MIKE ABOUT CYCLOTRON FREQ
            cycI=e*iB/(mi)
            # short hand for bigger eqs
            X= (cycE*cycI/((2*np.pi)**2))/(EISvin[itime] * ve)
            #
            # pedersen
            sigP[itime,iH]=sig0*((1+X)*(ve**2))/(((1+X)**2) * ve**2 + cycE**2)
            sigH[itime,iH]=sig0*cycE*ve/((1+X)**2 *(ve**2) + cycE**2)
#
# whew that should give us the superDARN conductivities
# now plot these quantities
#
dataC=[sigP,sigH]
labels=['Pedersen Conductivity', 'Hall Conductivity']
vminC=[2*1e-5,1.7*1e-6]
vmaxC=[2.3*1e-5,1.9*1e-6]
#
# go to CUTLASS directory
os.chdir('CUTLASS')
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
   plt.subplots_adjust(right=0.7, top=0.92, bottom=0.28, left=0.11)
   ax.set_ylabel('Range Gates', fontsize=22, fontweight='bold')
   Altitudes=range(75)
   ax.set_ylim(25, 35)
   ax.set_xlim(datetime.datetime(2015,3,12,10,0,0), datetime.datetime(2015,3,12,12,0,0))
   #ax.set_xlim(time[0], time)# last 17 missing
   X,Y=np.meshgrid(TnosSD,Altitudes)
   data=ma.masked_invalid(dataC[iDat]).transpose()
   vmin=vminC[iDat]
   vmax=vmaxC[iDat]
   ax.set_xlabel("Time", fontsize=20, fontweight='bold')
   col=ax.pcolormesh(X,Y,data, cmap='viridis', vmin=vmin, vmax=vmax)
   font = {'family' : 'normal',
           'weight' : 'bold',
           'size'   : 22}
   plt.rc('font', **font)
   cbaxes = fig.add_axes([0.75, 0.28, 0.03, 0.63])
   cb = plt.colorbar(col, cax = cbaxes,format=FuncFormatter(sciNot),ticks=np.linspace(vmin,vmax,5)) 
   cb.set_label(labels[iDat]+' [S/m]', fontsize=25)
   fig.set_size_inches(13,9)
   subdir_name='SuperDARN_Conductivity'
   if not os.path.exists(subdir_name):
       os.umask(0) # unmask if necessary
       os.makedirs(subdir_name) 
   os.chdir(subdir_name)#
   plt.savefig('20150312_'+labels[iDat]+'.png')
   os.chdir('..')
   plt.close()    

os.chdir('CUTLASS')
