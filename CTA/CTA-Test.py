#!/usr/bin/python
#
# Read Effective area for CTA
#
#
import re
import numpy as np
import math
import EBL

eblfilename = 'EBL_Gilmore.txt'
ebl = EBL.EBL(eblfilename)

#bckname = "CTA-Performance-North-5h-Background.txt"
#bckg = np.loadtxt(bckname)
#background = 0
#for i in range(0,len(bckg)):
#    background = background + bckg[i][2] * 3600 * 5
#print background
# This results in 7885.17828094 counts

background = 7885.17828094
    
# I asssume I know I have 80 energy bins

filename = "CTA-Performance-North-5h-EffArea.txt"
EffArea = np.loadtxt(filename)

# Generate spectrum as a numpy array
# Begin with 1e-12 (E/TEV)^(-2.1) 1/(TeV.cm^2.s)
# This is a very bright spectrum approximately (but not exactly) that of a
# single point source responsible for the IceCube diffuse flux

spectrum = np.empty([80,1])
energy = np.empty([80,1])
area = np.empty([80,1])
delta = np.empty([80,1])

redshift = 0.3

signal = 0
for i in range(0,77):
    spectrum[i] = 1e-12* math.pow(EffArea[i,0],-2.1) * math.exp(-ebl.TAU(EffArea[i,0],redshift))
    energy[i] = EffArea[i][0]
    area[i] = EffArea[i][1]
    lowlog = math.log10(energy[i])-0.025
    hilog = math.log10(energy[i])+0.025
    delta[i] = math.pow(10,hilog)-math.pow(10,lowlog)    
    # 10,000 to convert to cm^2
    # Delta is binwidth
    # 3600 * 5 seconds of signal
    signal = signal + spectrum[i] * area[i] * 10000 * delta[i] * 3600 * 5
print "Signal", signal, " Background", background

print "Significance", signal/math.sqrt(background)
