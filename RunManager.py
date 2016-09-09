#!/usr/bin/python

import os

density = 1e-8

for i in range(1,100):
    filename = "Results/Firesong_Density_" + str(density) + "_" + str(i).zfill(3) + ".txt"
    os.system("./StandardCandleSFR.py -o " + filename + " -d " + str(density) + " -p")
    os.system("gzip " + filename)
