#!/usr/bin/python
#
#

from __future__ import division
import glob
import subprocess
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('input', metavar='input', type=str, 
                    help='Input file')
parser.add_argument('-o', action='store', dest='output',default= 'CTA-Firesong.out',
                    help='Output filename')
options = parser.parse_args()

try:
    alerts = open(options.input)
except:
    print "Couldn't open " + options.input
    quit()
try:
    output = open(options.output,"w")
except:
    print "Couldn't open" + options.output
    quit()
        
Observed05 = 0
Observed5 = 0
Observed50 = 0
Neutrinos = 0

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
        out = subprocess.check_output("./CTA-Sensitivity.py --ObsTime 0.5 --index -2.13 --fluxnorm " + str(f) + " --redshift " + str(z), shell=True)
        significance = float(out[2:-2])
        if float(significance)>5:
            Observed05 = Observed05+1
            output.write('{:.3f} {:.4f} {:.6e}\n'.format(dec, z, f))
        out = subprocess.check_output("./CTA-Sensitivity.py --ObsTime 5 --index -2.13 --fluxnorm " + str(f) + " --redshift " + str(z), shell=True)
        significance = float(out[2:-2])
        if significance>5:
            Observed5 = Observed5+1
            output.write('{:.3f} {:.4f} {:.6e}\n'.format(dec, z, f))
        out = subprocess.check_output("./CTA-Sensitivity.py --ObsTime 50 --index -2.13 --fluxnorm " + str(f) + " --redshift " + str(z), shell=True)
        significance = float(out[2:-2])
        if significance>5:
            Observed50 = Observed50+1
            output.write('{:.3f} {:.4f} {:.6e}\n'.format(dec, z, f))
    Neutrinos = Neutrinos +1
    if Neutrinos%50==0:
        print "Processed Neutrinos:", Neutrinos
Prob05 = Observed05/Neutrinos
Prob5 = Observed5/Neutrinos
Prob50 = Observed50/Neutrinos
output.write('{:4f} {:4f} {:4f} \n'.format(Prob05,Prob5,Prob50))

alerts.close()
output.close()
