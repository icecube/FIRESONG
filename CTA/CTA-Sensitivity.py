#!/usr/bin/python
#
import re
import numpy as np
import math
import EBL
import argparse

#
# Process command line options
#
parser = argparse.ArgumentParser()
parser.add_argument('--ObsTime', action='store', dest='hours', default="50",
                    help='Observation time in hours')
parser.add_argument("--index", action="store", dest='index', type=float, default=-2.1,
                    help="Spectral Index")
parser.add_argument("--fluxnorm", action="store", dest='fluxnorm', type=float, default=1e-12,
                    help="Flux normalization")
parser.add_argument("--redshift", action="store", dest='redshift', type=float, default=0.,
                    help="Redshift")
options = parser.parse_args()

bckname = "Performance/CTA-Performance-North-" + options.hours + "h-Background.txt"
areaname = "Performance/CTA-Performance-North-" + options.hours + "h-EffArea.txt"
eblname = 'EBL_Gilmore.txt'

ebl = EBL.EBL(eblname)

bckg = np.loadtxt(bckname)
background = 0
for i in range(0,len(bckg)):
    background = background + bckg[i][2] * 3600 * float(options.hours)
# ->>> This results in  7885.17828094 counts in the 5 hour search
# ->>> This results in 44411.1408094 counts in the 50 hour search

EffArea = np.loadtxt(areaname)

spectrum = np.empty([80,1])
energy = np.empty([80,1])
area = np.empty([80,1])
delta = np.empty([80,1])

signal = 0
# HACK ALERT!!!!! I integrate only the first 77 entries, to ensure the energy is 100 TeV or lower
# The EBL attenuation is not available above 100 TeV
for i in range(0,77):
#   Generate spectrum as a numpy array
#   Units: 1/(TeV.cm^2.s)
    spectrum[i] = options.fluxnorm * math.pow(EffArea[i,0],options.index) * math.exp(-ebl.TAU(EffArea[i,0],options.redshift))
    energy[i] = EffArea[i][0]
    area[i] = EffArea[i][1]
    lowlog = math.log10(energy[i])-0.025
    hilog = math.log10(energy[i])+0.025
    delta[i] = math.pow(10,hilog)-math.pow(10,lowlog)    
    # 10,000 to convert to cm^2
    # Delta is binwidth
    # 3600 * 5 seconds of signal
    signal = signal + spectrum[i] * area[i] * 10000 * delta[i] * 3600 * float(options.hours)
#print "Signal", signal, " Background", background
print signal/math.sqrt(background)
