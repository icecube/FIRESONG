#!/usr/bin/python

#	Calculate the number of observed events by carrying out numerical integration:
#	N_evt = Sum{Flux * Exposure * dE * dSteradian} 


import numpy as np

#Diffuse flux is expressed as E^2dN/dE = fit * (E/100TeV)^index GeV/cm^2/sr/s
def IceCubeEvt(fit, index):
	#Get file name for all truncated data files
	filelist = ["multi-year-diffuse-exposure-cosZenMin0.05-cosZenMax0.10-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin0.00-cosZenMax0.05-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin-0.05-cosZenMax0.00-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin-0.10-cosZenMax-0.05-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin-0.15-cosZenMax-0.10-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin-0.20-cosZenMax-0.15-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin-0.25-cosZenMax-0.20-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin-0.30-cosZenMax-0.25-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin-0.35-cosZenMax-0.30-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin-0.40-cosZenMax-0.35-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin-0.45-cosZenMax-0.40-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin-0.50-cosZenMax-0.45-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin-0.55-cosZenMax-0.50-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin-0.60-cosZenMax-0.55-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin-0.65-cosZenMax-0.60-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin-0.70-cosZenMax-0.65-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin-0.75-cosZenMax-0.70-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin-0.80-cosZenMax-0.75-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin-0.85-cosZenMax-0.80-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin-0.90-cosZenMax-0.85-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin-0.95-cosZenMax-0.90-truncated-200TeV.txt",
				"multi-year-diffuse-exposure-cosZenMin-1.00-cosZenMax-0.95-truncated-200TeV.txt"]

	#This is the total event
	Total = 0

	#sinDec list
	sinDec = np.arange(-0.075, 1.025, 0.05)
	i = 0

	evtlist = [0]*len(sinDec)

	#print('# sin(Declin)  Observed_Events')
	#loop through all files
	for fn in filelist :
		fn = "IceCubeEffectiveArea/" + fn
		infile = open(fn, 'r')

		#Set array for storing the energy and the exposure
		e_list = np.array([])
		expose_list = np.array([])

		#Read each line
		for line in infile :
			line = line.split()
			if line[0] != '#' :
				e_list = np.append(e_list, float(line[0]))
				expose_list = np.append(expose_list, float(line[1]))

		infile.close()

		#number of event per energy (i.e. flux * exposure * delta_steradian)
		Devt = 1e4*0.5*fit*np.power(e_list/1e5, index)*np.power(e_list, -2)*expose_list*2*np.pi*0.05

		#get mid point of the number of event, too get more accurate result
		Devt_mid = Devt[:-1] + np.diff(Devt)/2.

		#multiply with the bin width of energy, i.e. delta_E
		evt = Devt_mid * np.diff(e_list)

		#print the total number of event in each steradian bin
		#print (str(sinDec[i])+ " " + str(np.sum(evt)))
		evtlist[i] = np.sum(evt)
		#increase i
		i += 1

		#print the total number of event for all steradian bin
		#Total += np.sum(evt)

	#print ('# The total observed events above 200TeV is:' + str(Total))

	return evtlist
