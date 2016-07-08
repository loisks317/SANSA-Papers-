# readProcessedEISCAT.py
#
# the EISCAT data is in individual ascii files labelled with the parameter
# contained within. The data is processed via the guisdap file Amore'
# provided... there are probably things to worry about it
#
# LKS, January 2016. SANSA.
#
import numpy as np
import os
import datetime as dt
#
#
def ReadData(date):
    
     #
     # read the EISCAT data processed files
     # return a data array with a whole of lot parameters
     #
     
     #
     # change to data directory
     os.chdir('/Users/loisks/Documents/ResearchProjects/SANSAProject/Processed30secondEISCAT/'+date)
     #
     # create python dictionary to hold the files
     data={}
     fileParams=['alt', 'az', 'el', 'ne', 'ti', 'te', 'vin', 'vel', 'comp1', 'evel', 'ene', 'evel', 'resfit', 'resfitsd', 'Pt',  'status', 'pp', 'pprange']
     #
     # each one of these requires it's own special tweaking
     # easy one
     for iParam in fileParams:
         data[iParam]=np.genfromtxt(iParam+'.txt', delimiter=',')
     
     
     # gross datetimes
     # needs to be put into datetimes
     with open('stim.txt','rb') as f:
         contentB=f.readlines()
     with open('etim.txt','rb') as f:
         contentE=f.readlines()
     datetimesBEG=[]
     datetimesEND=[]
     for iCon in range(len(contentB)):
         dateConvB=str(contentB[iCon])
         dateConvE=str(contentE[iCon])
         try:
             periodB=dateConvB.index('.')
             tempB=dateConvB[:periodB]
         except(ValueError):
              tempB=dateConvB[:periodB]+'0'
         try:
             periodE=dateConvE.index('.')
             tempE=dateConvE[:periodE]
         except(ValueError):
              tempE=dateConvE[:periodE]+'0'
         try:
             datetimesBEG.append(dt.datetime.strptime(tempB, '%Y,%m,%d,%H,%M,%S'))
         except(ValueError):
             # not the greatest correction
             tempErr=tempB.index('60')
             tempL=list(tempB)
             tempL[tempErr-3:tempErr-1]=str(int(tempB[tempErr-3:tempErr-1]))
             tempL[tempErr:tempErr+2]='00'
             tempB=''.join(tempL)
             datetimesBEG.append(dt.datetime.strptime(tempB, '%Y,%m,%d,%H,%M,%S'))
         try:
             datetimesEND.append(dt.datetime.strptime(tempE, '%Y,%m,%d,%H,%M,%S'))
         except(ValueError):
             # not the greatest correction
             tempErr=tempE.index('60')
             tempL=list(tempE)
             tempL[tempErr-3:tempErr-1]=str(int(tempE[tempErr-3:tempErr-1]))
             tempL[tempErr:tempErr+2]='0'
             tempE=''.join(tempL)
             datetimesEND.append(dt.datetime.strptime(tempE, '%Y,%m,%d,%H,%M,%S'))
     data['stim']=datetimesBEG
     data['etim']=datetimesEND
     
     os.chdir('..')
     return data
