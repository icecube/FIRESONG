#!/usr/bin/python

#	Calculate the number of observed events by carrying out numerical integration:
#	N_evt = Sum{Flux * Exposure * dE * dSteradian} 

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--fluxnorm", action="store",
                    dest="fluxnorm", default=0.9e-8,
                    help="Flux normalization in A(E/100TeV)^-index units GeV.cm^-2.s^-1.sr^-1")
parser.add_argument("--index", action="store",
                    dest="index", default=2.13,
                    help="Spectral index (positive)")
options = parser.parse_args()
index = float(options.index)
fluxnorm = float(options.fluxnorm)

print ("##############################################################################")
print ("Diffuse flux event number calculation")
print ("Using Aachen's Diffuse flux exposure above 200 TeV")
print ("Spectrum assumed: " + str(options.fluxnorm) + "(E/100TeV)^-" + str(options.index) + " GeV.cm^-2.s^-1.sr^-1")



#Get file name for all truncated data files
filelist = open('filelist.txt', 'r')

#This is the total event
Total = 0

#sinDec list
sinDec = np.arange(-0.075, 1.025, 0.05)
i = 0

print('# sin(Declin)  Observed_Events')
#loop through all files
for fn in filelist :
	fn = fn.rstrip('\n')
	input = open(fn, 'r')

	#Set array for storing the energy and the exposure
	e_list = np.array([])
	expose_list = np.array([])

	#Read each line
	for line in input :
		line = line.split()
		if line[0] != '#' :
			e_list = np.append(e_list, float(line[0]))
			expose_list = np.append(expose_list, float(line[1]))

	input.close()

	#number of event per energy (i.e. flux * exposure * delta_steradian)
	Devt = 1e4*0.5*fluxnorm*np.power(e_list/1e5, -(index-2.))*np.power(e_list, -2)*expose_list*2*np.pi*0.05

	#get mid point of the number of event, too get more accurate result
	Devt_mid = Devt[:-1] + np.diff(Devt)/2.

	#multiply with the bin width of energy, i.e. delta_E
	evt = Devt_mid * np.diff(e_list)

	#print the total number of event in each steradian bin
	print (str(sinDec[i])+ " " + str(np.sum(evt)))

	#increase i
	i += 1

	#print the total number of event for all steradian bin
	Total += np.sum(evt)

print ('# The total observed events above 200TeV is:' + str(Total))

filelist.close()
