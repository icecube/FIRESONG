#!/usr/bin/python

import os

density = 3e-8

for i in range(1,10000):
    filename = "Firesong_SFH_Density_" + str(density) + "_" + str(i).zfill(4) + ".txt"
    os.system("./StandardCandle.py -o " + filename + " -d " + str(density) + " --hawc")
    os.system("gzip " + filename)
