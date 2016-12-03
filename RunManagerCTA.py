#!/usr/bin/python

import os
import shutil

density = 1e-6

for i in range(1,100):
    filename = "Firesong_SFH_Density_" + str(density) + "_" + str(i).zfill(4) + ".txt"
    os.system("./StandardCandle.py -o " + filename + " -d " + str(density) + " --cta")
    os.system("gzip " + filename)
    shutil.move(filename + ".gz", "Results/CTA_North/1e-06")
