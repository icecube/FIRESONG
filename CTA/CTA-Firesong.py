#!/usr/bin/python
#
#

import glob
import subprocess
import re

list = glob.glob("../Results/CTA/cta*")

TotalObserved = 0
TotalNeutrinos = 0

runsummary = open("RunSummary.txt","w")
for filename in list:
    file = open(filename,"r")
    Observed = 0
    Neutrinos = 0
    for line in file:
        Neutrinos = Neutrinos +1
        (dec,z,f,obs) = line.split()
        out = subprocess.check_output("./CTA-Sensitivity.py --ObsTime 5 --index -2.13 --fluxnorm " + str(f) + " --redshift " + z, shell=True)        
        significance = out[2:-2]
        if float(significance)>5:
            Observed = Observed+1
            print dec,z,f,obs,significance
    runsummary.write(str(Observed)+ " " + str(Neutrinos) + "\n")
    file.close()
    TotalObserved = TotalObserved + Observed
    TotalNeutrinos = TotalNeutrinos + Neutrinos
runsummary.write("#" + str(TotalObserved) + " " + str(TotalNeutrinos) + "\n")
runsummary.close()
