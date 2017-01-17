#!/usr/bin/python
#
from __future__ import division
import re
import numpy as np
import math
import EBL
import argparse

def Significance(fluxnorm,index,redshift,obsTime):
# units: fluxnor = GeV/(cm^2.s), obsTime = h
    bckname = "Performance/CTA-Performance-North-" + obsTime + "h-Background.txt"
    areaname = "Performance/CTA-Performance-North-" + obsTime + "h-EffArea.txt"
    eblname = 'EBL_Gilmore.txt'

    ebl = EBL.EBL(eblname)

    bckg = np.loadtxt(bckname) # units: Hz
    background = 0
    for i in range(0,len(bckg)):
        background = background + bckg[i][2] * 3600 * float(obsTime)
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
        energy[i] = EffArea[i][0]
        area[i] = EffArea[i][1] # m^2
        spectrum[i] = fluxnorm/1000. * math.pow(energy[i],index) * math.pow(100., -index-2.0) * math.exp(-ebl.TAU(energy[i],redshift))
        # 1/(TeV.cm^2.s) = TeV/cm^2.s * TeV^-2.13 * TeV^0.13 *ebl
        lowlog = math.log10(energy[i])-0.025
        hilog = math.log10(energy[i])+0.025
        delta[i] = math.pow(10,hilog)-math.pow(10,lowlog)  # TeV  
        # 10,000 to convert to cm^2
        # Delta is binwidth
        # 3600 * 5 seconds of signal
        signal = signal + spectrum[i] * area[i] * 10000 * delta[i] * 3600 * float(obsTime)
        # N_sig = N_sig + 1/(TeV.cm^2.s) * m^2 * cm^2/m^2 * TeV * s
    return (signal/math.sqrt(background))[0]

