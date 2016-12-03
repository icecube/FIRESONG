#!/usr/bin/python
#
#

import glob
import subprocess
import re

list = glob.glob("../Results/CTA/cta*")

TotalObserved05 = 0
TotalObserved5 = 0
TotalObserved50 = 0
TotalNeutrinos = 0

runsummary = open("RunSummary.txt","w")
for filename in list:
    file = open(filename,"r")
    Observed05 = 0
    Observed5 = 0
    Observed50 = 0   
    Neutrinos = 0
    for line in file:
        Neutrinos = Neutrinos +1
        (dec,z,f,obs) = line.split()
        out = subprocess.check_output("./CTA-Sensitivity.py --ObsTime 0.5 --index -2.13 --fluxnorm " + str(f) + " --redshift " + z, shell=True)     
        significance = out[2:-2]
        if float(significance)>5:
            Observed05 = Observed05+1
            print dec,z,f,obs,significance
        out = subprocess.check_output("./CTA-Sensitivity.py --ObsTime 5 --index -2.13 --fluxnorm " + str(f) + " --redshift " + z, shell=True)       
        significance = out[2:-2]
        if float(significance)>5:
            Observed5 = Observed5+1
            print dec,z,f,obs,significance
        out = subprocess.check_output("./CTA-Sensitivity.py --ObsTime 50 --index -2.13 --fluxnorm " + str(f) + " --redshift " + z, shell=True)      
        significance = out[2:-2]
        if float(significance)>5:
            Observed50 = Observed50+1
            print dec,z,f,obs,significance
    runsummary.write(str(Observed05)+ " " + str(Observed5) + " " + str(Observed50) + " " + str(Neutrinos) + "\n")
    file.close()
    TotalObserved05 = TotalObserved05 + Observed05
    TotalObserved5 = TotalObserved5 + Observed5
    TotalObserved50 = TotalObserved50 + Observed50
    TotalNeutrinos = TotalNeutrinos + Neutrinos
runsummary.write("#" + str(TotalObserved05) + " " + str(TotalObserved5) + " " + str(TotalObserved50) + " " + str(TotalNeutrinos) + "\n")
runsummary.close()
