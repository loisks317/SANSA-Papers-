# getSANAEB.py
#
# This REQUIRES python 2.7!!!!
#
#
import numpy as np
import os

#
#


os.chdir('/Users/loisks/Documents/ResearchProjects/SANSAProject/ProcessedEISCAT/')
alt=np.genfromtxt('alt.txt')
phi=19
theta=69.5
B=[]
os.chdir('/Users/loisks/Desktop/Functions')
import IGRF
for h in range(len(alt)):
    res=IGRF.IGRF_fortran( phi, theta , (alt[h])+6374 ,int(2015))
    B.append(np.sqrt(res[0]**2 + res[1]**2 + res[2]**2))

os.chdir('/Users/loisks/Documents/ResearchProjects/SANSAProject/ProcessedEISCAT/')
file=open('EISCATB.txt', 'wb')
for item in B:
    print >> file, item
        
