from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import random
import math
import cosmolopy
from scipy.interpolate import UnivariateSpline
import argparse

# Plack 2015 parameters
cosmology = {'omega_M_0' : 0.308, 'omega_lambda_0' : 0.692, 'h' : 0.678}
cosmology = cosmolopy.distance.set_omega_k_0(cosmology) #Flat universe

#StarFormationHistory (SFR), from Hopkins and Beacom, unit = M_sun/yr/Mpc^3 
def StarFormationHistory(x):
    if x<0.30963:
        return math.pow(10,3.28*x-1.82)
    if (x>=0.30963)and(x<0.73878):
        return math.pow(10,-0.26*x-0.724)
    if x>=0.73878:
        return math.pow(10,-8.0*x+4.99)

def ObservedRate(x):
    return 4*np.pi*StarFormationHistory(x)*cosmolopy.distance.diff_comoving_volume(math.pow(10,x)-1, **cosmology)

#Generate the bins
numbin = 10000                                     # this one is arbitrary right now
bins = range(0, numbin)
for i in range(0,numbin):
	bins[i] = bins[i]/numbin

#Generate the histogram
binmid = bins[:-1] + np.diff(bins)/2.
hist = [ObservedRate(binmid[i]) for i in range(0,len(binmid))]

#Generate the random number of log(1+z)
cdf = np.cumsum(hist)
cdf = cdf / cdf[-1]
values = np.random.rand(100000)                     # We will have 100000 trials
value_bins = np.searchsorted(cdf, values)
random_from_cdf = binmid[value_bins]

# This is the number of events as a function of declination (Effective Area?)
def IceCubeResponse(sinDec):
    sinDecbin = np.array([-0.075,-0.025,0.025,0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775,0.825,0.875,0.925,0.975])
    response = np.array([105.896, 185.11, 177.775, 155.445, 133.492, 113.77, 99.8513, 87.4213, 76.7264, 63.8566, 57.6465, 48.2527, 44.1256, 34.2095, 30.4975, 25.8617, 22.4174, 19.0691, 15.683, 9.20381, 5.12658, 3.40891])
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

#Standard Canle, all source has same lumionsity at z = 1.0
dL1 = cosmolopy.distance.luminosity_distance(1.0, **cosmology)


#The idea is to put all calculation outside for loop
#Generate the declination
sinDec = np.random.rand(100000)
sinDec = sinDec*2. -1.
#Get z
z = np.power(10, random_from_cdf)-1
#get luminosity distance at z
dL = cosmolopy.distance.luminosity_distance(z, **cosmology)
#get flux at z
flux = 1.86e-13 * (dL1*dL1)/(dL*dL) 

#total flux
TotalFlux = np.sum(flux)

#Calculate the effectivearea at icecube for each source
EA = [0]*len(sinDec)
for i in range(0, len(sinDec)):
	EA[i] = IceCubeResponse(sinDec[i])

#Get mean no. of event due to each source
events = EA*flux/1e-8
#Get the number of event from poisson distribution
obs = np.random.poisson(events, (1, len(events)))

#Calculate the declination
declin = 180*np.arcsin(sinDec)/np.pi

### The following two lines may be a faster way to output the array, will save it for later
#finalset = zip(declin, z, flux, obs[0])
#np.savetxt(options.filename, finalset, delimiter=" ", fmt='%f, %f, %f, %i')

for i in range(0, len(random_from_cdf)):
	output.write(str(declin[i]) + " " + str(z[i]) + " " + str(flux[i]) + " " + str(obs[0][i]) + "\n")

print ('Total Flux is ' + str(TotalFlux))
output.write("# " + str(TotalFlux) + "\n")

output.close()








