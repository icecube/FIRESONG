#!/usr/bin/python
#
# x = log10(1+z)

import numpy as np
import scipy
import cosmolopy
cosmology = {'omega_M_0' : 0.308, 'omega_lambda_0' : 0.692, 'h' : 0.678}
cosmology = cosmolopy.distance.set_omega_k_0(cosmology) #Flat universe

# This is the redshift distribution for arbitrary evolution for either transient or steady sources 
def RedshiftDistribution(z,options):
    if options.Transient == False:
        return 4*np.pi* Evolution(np.log10(1+z),options.Evolution) * cosmolopy.distance.diff_comoving_volume(z,**cosmology)
    else:
        return 4*np.pi* Evolution(np.log10(1+z),options.Evolution) * 1/(1+z) * cosmolopy.distance.diff_comoving_volume(z,**cosmology)
      
def Evolution(x,evol):
    if (evol=="HB2006SFR"):
        return HopkinsBeacom2006StarFormationRate(x)
    elif (evol=="NoEvolution"):
        return NoEvolution(x)
    elif (evol=="CC2015SNR"):
        return CandelsClash2015SNRate(x)
    else:
        print "Source evolution " +  evol +  " not recognized"
        quit()

#StarFormationHistory (SFR), from Hopkins and Beacom 2006, unit = M_sun/yr/Mpc^3 
def HopkinsBeacom2006StarFormationRate(x):
    if x<0.30963:
        return np.power(10,3.28*x-1.82)
    if (x>=0.30963)and(x<0.73878):
        return np.power(10,-0.26*x-0.724)
    if x>=0.73878:
        return np.power(10,-8.0*x+4.99)

def NoEvolution(x):
    return 1.

def CandelsClash2015SNRate(x):
    a = 0.015
    b = 1.5
    c = 5.0
    d = 6.1
    density = a*(10.**x)**c / ((10.**x / b)**d+1.)
    return density

def StandardCandleSources(options):
  norm = scipy.integrate.quad(lambda z: RedshiftDistribution(z, options), 0, options.zmax)[0]
  area = scipy.integrate.quad(lambda z: RedshiftDistribution(z, options), 0, 0.01)[0]

  vlocal = cosmolopy.distance.comoving_volume(0.01, **cosmology)
  Ntotal = options.density * vlocal / (area/norm)
  dL1 = dL1 = cosmolopy.distance.luminosity_distance(1.0, **cosmology)
  Fluxnorm = 4*np.pi*options.fluxnorm / scipy.integrate.quad(lambda z: Ntotal*dL1*dL1/np.power(cosmolopy.distance.luminosity_distance(z, **cosmology), 2)*RedshiftDistribution(z, options)/norm*((1+z)/2.)**(-options.index+1), 0, options.zmax)[0]
  return [int(Ntotal), Fluxnorm] 

# Wrapper fucntion - so that cosmolopy is only imported here.
def LuminosityDistance(z):
    return cosmolopy.distance.luminosity_distance(z, **cosmology)
