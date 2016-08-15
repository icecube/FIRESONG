#!/usr/bin/python
#
# Sumulation plan:
# The flux is known
# a) Find the local source density to match IceCube's flux assumming a given luminosity function and assuming SFR distribution
# b) Sample from SFR and Luminosity function. 
#
#

import random
import math
import cosmolopy
import numpy
from scipy.interpolate import UnivariateSpline
import argparse

#cosmology = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
# Plack 2015 parameters
cosmology = {'omega_M_0' : 0.308, 'omega_lambda_0' : 0.692, 'h' : 0.678}
cosmology = cosmolopy.distance.set_omega_k_0(cosmology)

def StarFormationHistory(x):
    if x<0.30963:
        return math.pow(10,3.28*x-1.82)
    if (x>=0.30963)and(x<0.73878):
        return math.pow(10,-0.26*x-0.724)
    if x>=0.73878:
        return math.pow(10,-8.0*x+4.99)

def ObservedRate(x):
    return 4*3.141592*StarFormationHistory(x)*cosmolopy.distance.diff_comoving_volume(math.pow(10,x)-1, **cosmology)

# This is the number of events as a function of declination 
def IceCubeResponse(sinDec):
    sinDecbin = numpy.array([-0.075,-0.025,0.025,0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775,0.825,0.875,0.925,0.975])
    response = numpy.array([105.896, 185.11, 177.775, 155.445, 133.492, 113.77, 99.8513, 87.4213, 76.7264, 63.8566, 57.6465, 48.2527, 44.1256, 34.2095, 30.4975, 25.8617, 22.4174, 19.0691, 15.683, 9.20381, 5.12658, 3.40891])
    spline = UnivariateSpline(sinDecbin,response)
    if sinDec>=-0.1 and sinDec<=1.:
        return spline(sinDec)
    else:
        return 0

parser = argparse.ArgumentParser()
parser.add_argument('-f', action='store', dest='filename',
                    help='Output filename')
options = parser.parse_args()

output = open(options.filename,"w")



    
dL1 = cosmolopy.distance.luminosity_distance(1.0, **cosmology)

TotalFlux = 0

# I think I need this number to be 18,570,000 to match 10^-6 /Mpc^3


for i in range(185700):
    x = random.random() # This is log10(1+z) 
    z = math.pow(10,x)-1 # redshift
    test = random.random()*7e10 # SFH * dVc * 4pi
#    test = random.random()*0.2
    if test < ObservedRate(x):        
#    if test < StarFormationHistory(x):
        dL = cosmolopy.distance.luminosity_distance(z, **cosmology)
        #This number should be 2.15e-15 for 18,570,000 trials above (rate 10^-6)
        #This number should be 2.15e-14 for 1,857,000 trials (rate 10^-7)
        #This number should be 2.15e-13 for 185,700 trials (rate 10^-7)
        flux = 2.15e-13*(dL1*dL1)/(dL*dL) # This numbers magically gives a flux normalization of 10^-8 GeV/cm^2.s.sr. It's otherwise meaningless.
        # It needs to be matched to the number of trials in range() above.

        sinDec = random.random()*2. -1.
        events = IceCubeResponse(sinDec)*flux/1e-8

        # Calculate the number of observed events
        obs = numpy.random.poisson(events,1)
        output.write(str(180*math.asin(sinDec)/3.141592) + " " + str(z) + " " + str(flux) + " " + str(obs[0]) + "\n")
        TotalFlux = TotalFlux + flux
output.write("# " + str(TotalFlux) + "\n")

output.close()
