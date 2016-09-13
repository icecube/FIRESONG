#!/usr/bin/python

import os

density = 1e-8

for i in range(1,20):
    filename = "Results/Firesong_SFH_Density_" + str(density) + "_" + str(i).zfill(3) + ".txt"
    os.system("./StandardCandle.py -o " + filename + " -d " + str(density) + " -p")
    os.system("gzip " + filename)
