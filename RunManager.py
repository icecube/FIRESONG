#!/usr/bin/python

import os

density = 1e-8

for i in range(1,100):
    os.system("./StandardCandleSFR.py -o Results/Firesong_Density_" + str(density) + "_" + str(i) + ".txt -d " + str(density) + " -p")
