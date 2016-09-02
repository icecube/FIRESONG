from __future__ import division
import numpy as np
import random
import math
import cosmolopy
from scipy.interpolate import UnivariateSpline
import argparse
import scipy

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

def ObservedRate(z):
    return 4*np.pi*StarFormationHistory(np.log10(1+z))*cosmolopy.distance.diff_comoving_volume(z, **cosmology)

#Total area under the Observed Rate function
norm = scipy.integrate.quad(lambda z: ObservedRate(z), 0, 10)[0]

#area under the Observed Rate, from z=0 to z=0.1
area = scipy.integrate.quad(lambda z: ObservedRate(z), 0, 0.1)[0]

#comoving volume of the local universe
vlocal = cosmolopy.distance.comoving_volume(0.1, **cosmology)

#guess of neutrino source density
rho0 = 1e-6

Ntotal = rho0 * vlocal / (area/norm)

print Ntotal

#Max value of Observed Rate Function is 7.412E10
Ntotalmc = rho0 * vlocal / (area / (7.412E10 * 10))

print Ntotalmc

dL1 = cosmolopy.distance.luminosity_distance(1.0, **cosmology)
def dL(z):
	return cosmolopy.distance.luminosity_distance(z, **cosmology)

Flux_norm = 1e-8 / scipy.integrate.quad(lambda z: Ntotal*dL1*dL1/(dL(z)*dL(z))*ObservedRate(z)/norm, 0, 10)[0]

print Flux_norm












