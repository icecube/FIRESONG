#!/usr/bin/python

#	Calculate the number of observed events by carrying out numerical integration:
#	N_evt = Sum{Flux * Exposure * dE * dSteradian} 
import os
import numpy as np

try:
    firesongdir = os.environ['FIRESONG']
except:
    print "Enviromental variable FIRESONG not set"
    quit()
icecubedir = firesongdir + "/IceCube/"

def IceCubeEvt(norm, index, notruncate):

	#Get file name for all truncated data files
	if notruncate == False:
		filelist = open(icecubedir + 'filelist.txt', 'r')
	else:
		filelist = open(icecubedir + 'filelist_2.txt', 'r')

	#sinDec list
	sinDec = np.arange(-0.075, 1.025, 0.05)
	i = 0

	evtlist = [0]*len(sinDec)

	#loop through all files
	for fn in filelist :
		fn = icecubedir + fn.rstrip('\n')
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
		Devt = 1e4*0.5*norm*np.power(e_list/1e5, -(index-2.))*np.power(e_list, -2)*expose_list*2*np.pi*0.05

		#get mid point of the number of event, too get more accurate result
		Devt_mid = Devt[:-1] + np.diff(Devt)/2.

		#multiply with the bin width of energy, i.e. delta_E
		evt = Devt_mid * np.diff(e_list)

		evtlist[i] = np.sum(evt)

		#increase i
		i += 1


	filelist.close()
	
	return evtlist
