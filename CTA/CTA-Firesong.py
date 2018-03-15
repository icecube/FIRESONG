#!/usr/bin/python
#
#

from __future__ import division
import glob
import subprocess
import re
import argparse
import CTAsensitivity_LiMa

#### HARDCODE ALERT
parser = argparse.ArgumentParser()
parser.add_argument('input', metavar='input', type=str, 
                    help='Input file')
parser.add_argument('-o', action='store', dest='output',default= 'CTA-Firesong.out',
                    help='Output filename')
parser.add_argument('--index', dest='index', action='store', type=float, default=2.13,
                    help='index of neutrino spectrum')
parser.add_argument("--transient", action='store_true',
                    dest='Transient', default=False,
                    help='Simulate transient sources, NOT TESTED YET!')
parser.add_argument("--timescale", action='store', dest='timescale', type=float,
                    default=1000., help='time scale of transient sources, default is 1000sec.')
parser.add_argument("--array", action='store', dest='array', type=str,
                    default='North', help='Choose between the CTA North or South, and between the BL and TS arrays: North, North-TS, South, South-TS')
parser.add_argument("--zenith", action='store', dest='zenith', type=str,
                    default='20deg', help='Choose between low zenith 20deg and high zenith 40deg observations with CTA')
options = parser.parse_args()

try:
    alerts = open(options.input)
except:
    print("Couldn't open " + options.input)
    quit()
try:
    output = open(options.output,"w")
except:
    print("Couldn't open" + options.output)
    quit()
        
Observed05 = 0
Observed5 = 0
Observed50 = 0
Neutrinos = 0
# CTA prod3b IRFs for 30 minutes observation
bckname_30m = "Performance/Prod3b/CTA-Performance-" +options.array+ "-" + options.zenith + "-average-30m_20170627_Background.txt"
areaname_30m = "Performance/Prod3b/CTA-Performance-" +options.array+ "-" + options.zenith + "-average-30m_20170627_EffArea.txt"
# CTA prod3b IRFs for 5 hours observation
bckname_5h = "Performance/Prod3b/CTA-Performance-" +options.array+ "-" + options.zenith + "-average-05h_20170627_Background.txt"
areaname_5h = "Performance/Prod3b/CTA-Performance-" +options.array+ "-" + options.zenith + "-average-05h_20170627_EffArea.txt"
# CTA prod3b IRFs for 50 hours observation
bckname_50h = "Performance/Prod3b/CTA-Performance-" +options.array+ "-" + options.zenith + "-average-50h_20170627_Background.txt"
areaname_50h = "Performance/Prod3b/CTA-Performance-" +options.array+ "-" + options.zenith + "-average-50h_20170627_EffArea.txt"

for neutrino in alerts:
    if re.search("#",neutrino):
        continue
    (dec,z,f) = neutrino.split()
    dec = float(dec)
    z = float(z)
    f = float(f)
# Gilmore table only allows to check up to z = 9
# These should be undetectable anyway
# It'd be better to make CTA-Sensisitivity in to a library
    if (z<9): 
        significance = CTAsensitivity_LiMa.Significance(f, z, 0.5, bckname_30m, areaname_30m, options)
        if float(significance)>5:
            Observed05 = Observed05+1
            #output.write('{0:.3f} {1:.4f} {2:.6e} {3:.2f}\n'.format(dec, z, f, significance))
        significance = CTAsensitivity_LiMa.Significance(f, z, 5., bckname_5h, areaname_5h, options)    
        if significance>5:
            Observed5 = Observed5+1
            #output.write('{0:.3f} {1:.4f} {2:.6e} {3:.2f}\n'.format(dec, z, f, significance))
        significance = CTAsensitivity_LiMa.Significance(f, z, 50., bckname_50h, areaname_50h, options)    
        if significance>5:
            Observed50 = Observed50+1
            #output.write('{0:.3f} {1:.4f} {2:.6e} {3:.2f}\n'.format(dec, z, f, significance))
    Neutrinos = Neutrinos +1
    if Neutrinos%50==0:
        print("Processed Neutrinos:", Neutrinos)
Prob05 = Observed05/Neutrinos
Prob5 = Observed5/Neutrinos
Prob50 = Observed50/Neutrinos
output.write('{0:4f} {1:4f} {2:4f} \n'.format(Prob05,Prob5,Prob50))

alerts.close()
output.close()
