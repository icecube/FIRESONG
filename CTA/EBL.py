#!/usr/bin/python
#

import math

class EBL:
    def __init__(self,filename):
        self.filename = filename
        file = open(self.filename)
        while 1:
            line = file.readline()
            splitline = line.split()
            if splitline[0] == 'BINZ':
                self.zbins = len(splitline)-1
                self.z = splitline[1:len(splitline)]
            if splitline[0] == 'BINE':
                self.Ebins = len(splitline)-1
                self.E = splitline[1:len(splitline)]
            if splitline[0] == 'TAU':
                break
        for i in range(self.zbins):
            self.z[i] = float(self.z[i])
        for i in range(self.Ebins):
            self.E[i] = float(self.E[i]) # convert to GeV. Hardcoded, I know.
            
        # I assume we have read LOG10E and THETA
        self.tau = [[0]*self.zbins for x in range(self.Ebins)]
        for i in range(self.Ebins):
            line = file.readline()
            splitline = line.split()
            for j in range(self.zbins):
                self.tau[i][j] = float(splitline[j])
        file.close()
    def TAU(self,energy,redshift):
        if (redshift > 0):
        # We want to find the bin that is just below the value
            i = 0 # bin for energy
            while (energy > self.E[i+1]):
                i=i+1
            j = 0 # bin for redshift
            while (redshift > self.z[j+1]):
                j=j+1
            tau1 = self.tau[i][j] + ((self.tau[i][j+1] - self.tau[i][j]) * (redshift - self.z[j]) / (self.z[j+1] - self.z[j]));
            tau2 = self.tau[i+1][j] + ((self.tau[i+1][j+1] - self.tau[i+1][j]) * (redshift - self.z[j]) / (self.z[j+1]-self.z[j]));
            TAU = tau2 + (tau2-tau1) * (energy-self.E[i])/(self.E[i+1]-self.E[i]);
            return TAU
        else:
            return 0.
