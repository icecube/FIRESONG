#!/usr/bin/python
# Authors: Chris Tung
#          Igncio Taboada
# Date: December 5, 2016
      

from __future__ import division
import os
import numpy as np
import random
import math
import cosmolopy
import scipy.integrate
from scipy.interpolate import UnivariateSpline, InterpolatedUnivariateSpline
import argparse
from evt_calculation import IceCubeEvt
from fluxmodel import PowerLawFlux, LognormalFlux

# Plack 2015 parameters
cosmology = {'omega_M_0' : 0.308, 'omega_lambda_0' : 0.692, 'h' : 0.678}
cosmology = cosmolopy.distance.set_omega_k_0(cosmology) #Flat universe

#
# There are two options for Evolution: 
# a)Star Formation History
# b)No evolution

#StarFormationHistory (SFR), from Hopkins and Beacom, unit = M_sun/yr/Mpc^3 
def StarFormationHistory(x):
    if x<0.30963:
        return math.pow(10,3.28*x-1.82)
    if (x>=0.30963)and(x<0.73878):
        return math.pow(10,-0.26*x-0.724)
    if x>=0.73878:
        return math.pow(10,-8.0*x+4.99)

def NoEvolution(x):
    return 1.

def sourcefluxdistribution(fluxoption, mean, width, size):
    if fluxoption == 0:
        return mean
    if fluxoption == 1:
        return LognormalFlux(mean, size, width)
    if fluxoption == 2:
        return PowerLawFlux(mean, -2, size, width)

try:
    firesongdir = os.environ['FIRESONG']
except:
    print "Enviromental variable FIRESONG not set"
    quit()
outputdir = firesongdir + "/Results/"
#
# Process command line options
#
parser = argparse.ArgumentParser()
parser.add_argument('-o', action='store', dest='filename',default= 'Firesong.out',
                    help='Output filename')
parser.add_argument('-d', action='store', dest='density', type=float, default = 1e-9,
                    help='Local neutrino source density [1/Mpc^3]')
parser.add_argument("-p", action="store_false",
                    dest="NoPSComparison", default=True,
                    help="Calculate detectable point sources")
parser.add_argument("--hawc", action="store_false",
                    dest="NoHAWC", default=True,
                    help="Calculate detectable point sources for HAWC")
parser.add_argument("--cta", action="store_false",
                    dest="NoCTA", default=True,
                    help="Calculate detectable point sources for CTA")
parser.add_argument("--noevolution", action="store_false",
                    dest="NoEvolution", default=True,
                    help="Disable Star Formation History Evolution")
parser.add_argument("--transient", action='store_true',
                    dest='transient', default=False,
                    help='Simulate transient sources, NOT FULLY WORKING YET!')
parser.add_argument("--zmax", action="store",
                    dest="zmax", default=10.,
                    help="Highest redshift")
#The following two options normal and index, are for tweaking the IceCubeResponse function
parser.add_argument("--fluxnorm", action="store", dest='fluxnorm', type=float, default=0.9e-8,
                    help="flux normalization of the diffuse flux for IceCubeResponse")
parser.add_argument("--index", action="store", dest='index', type=float, default=2.13,
                    help="Index of the diffuse flux for IceCubeResponse")
parser.add_argument("--notruncate", action="store_true", dest='notruncate', default=False,
                    help="Switch to use non-truncated exposure data")
parser.add_argument("--fluxmodel", action="store", 
                    dest="fluxmodel", type=int, default=0,
                    help="Type of fluxmodel, 0 for standard model, 1 for lognormal, 2 for powerlaw")
parser.add_argument("--sigma", action="store",
                    dest="sigma", type=float, default=0.3,
                    help="The width of log normal distribution, in unit of dex, default is 0.3")
parser.add_argument("--nohistogram", action="store_false",
                    dest="histogram", default=True,
                    help="Output the statistic of events, default is True, recommended False for notruncate")

options = parser.parse_args()
output = open(outputdir+str(options.filename),"w")

if (options.NoEvolution==False):
    Evolution=NoEvolution
else:
    Evolution=StarFormationHistory

icecubeevtlist = IceCubeEvt(options.fluxnorm, options.index, options.notruncate)

zmax = float(options.zmax)

# This is the redshift distribution for arbitrary evolution    
def Redshift_distribution(z, transient):
    if transient == False:
        return 4*np.pi*Evolution(np.log10(1+z))*cosmolopy.distance.diff_comoving_volume(z, **cosmology)
    if transient == True:
        return 4*np.pi*Evolution(np.log10(1+z))*1/(1+z)*cosmolopy.distance.diff_comoving_volume(z, **cosmology)

def NumberOfSourcesStandardCandle(rho0, norm, transient):
  norm = scipy.integrate.quad(lambda z: Redshift_distribution(z, transient), 0, zmax)[0]
  area = scipy.integrate.quad(lambda z: Redshift_distribution(z, transient), 0, 0.01)[0]
  vlocal = cosmolopy.distance.comoving_volume(0.01, **cosmology)
  Ntotal = rho0 * vlocal / (area/norm)
  dL1 = dL1 = cosmolopy.distance.luminosity_distance(1.0, **cosmology)
  Fluxnorm = 4*np.pi*options.fluxnorm / scipy.integrate.quad(lambda z: Ntotal*dL1*dL1/np.power(cosmolopy.distance.luminosity_distance(z, **cosmology), 2)*Redshift_distribution(z, transient)/norm, 0, zmax)[0]
  return [int(Ntotal), Fluxnorm] 

N_sample, candleflux = NumberOfSourcesStandardCandle(options.density, options.fluxnorm, options.transient)

flux_z1 = sourcefluxdistribution(options.fluxmodel, candleflux, options.sigma, N_sample)

print ("##############################################################################")
print ("FIRESONG initializing")
if options.fluxmodel == 0:
    print ("Standard candle sources")
if options.fluxmodel == 1:
    print ("Lognormal distributed sources")
if options.fluxmodel == 2:
    print ("PowerLaw distributed sources")
print ("Star formation evolution? " + str(options.NoEvolution))
print ("Number of neutrinos sources in the Universe: " + str(N_sample))
print ("Uses neutrino diffuse flux: E^2 dN/dE = " + str(options.fluxnorm) + " (E/100 TeV)^(" + str(-(options.index-2.)) + ") GeV/cm^2.s.sr")
print ("Local density of neutrino sources: " + str(options.density) + "/Mpc^3")
print ("Use 200TeV truncated exposure data? " + str(not options.notruncate))
print ("Redshift range: 0 - " + str(options.zmax)) 
print ("FIRESONG initialization done")

#Generate the bins of z
redshift_bins = np.arange(0, zmax, 0.001)

#Generate the histogram
redshift_binmid = redshift_bins[:-1] + np.diff(redshift_bins)/2.
sfh = [Redshift_distribution(redshift_binmid[i], options.transient) for i in range(0,len(redshift_binmid))]

#Generate the random number of z
cdf = np.cumsum(sfh)
cdf = cdf / cdf[-1]
test = np.random.rand(N_sample)                     
bin_index = np.searchsorted(cdf, test)
redshift_list = redshift_binmid[bin_index]

# This is the number of events as a function of declination (Effective Area?)
def IceCubeResponse(sinDec):
    sinDecbin = np.array([-0.075,-0.025,0.025,0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775,0.825,0.875,0.925,0.975])
    response = np.array(icecubeevtlist)
    spline = UnivariateSpline(sinDecbin,response)
    SignalExpectation = spline(sinDec)
    if sinDec>=-0.1 and sinDec<=1. and SignalExpectation>0:
        return SignalExpectation
    else:
        return 0

#Standard Canle, all source has same lumionsity at z = 1.0, so we calculate the luminosity distance at z=1
dL1 = cosmolopy.distance.luminosity_distance(1.0, **cosmology)

#The idea is to put all calculations outside for loop
#Generate the declination
sinDec = np.random.rand(N_sample)
sinDec = sinDec*2. -1.
#Get z
z = redshift_list
#get luminosity distance at z
dL = cosmolopy.distance.luminosity_distance(z, **cosmology)
#get point source flux at z, the coefficient here should match the Total flux = 0.9e-8 Gev / cm^2 / s / sr
#for N_sample = 14921752, facctor is 2.20e-14
flux = flux_z1 * (dL1*dL1)/(dL*dL) 

#total flux
TotalFlux = np.sum(flux)


#Calculate the number of events at icecube for each source's declination angle
EA = [IceCubeResponse(sinDec[i]) for i in range(0, len(sinDec))]

#Get mean no. of event due to each source
events = EA*flux/(options.fluxnorm*2*np.pi*0.05)
#Get the number of event from poisson distribution
obs = np.random.poisson(events, (1, len(events)))

#Calculate the declination
declin = 180*np.arcsin(sinDec)/np.pi

### The following two lines may be a faster way to output the array, will save it for later
#finalset = zip(declin, z, flux, obs[0])
#np.savetxt(options.filename, finalset, delimiter=" ", fmt='%f, %f, %f, %i')

output.write("# FIRESONG Output description\n")
output.write("# Declination: degrees\n")
output.write("# Redshift\n")
output.write("# flux: E^2 dN/dE assuming " + str(options.fluxnorm) + " (E/100 TeV)^(" + str(-(options.index-2.)) + ") GeV/cm^2.s.sr\n")
output.write("#     Note that as of 2016, IceCube can detect point sources of ~10^-9 in the\n")
output.write("#     qunits used here.\n")
output.write("# Observed: Number of >200 TeV neutrino events detected, using 6 year Diffuse effective area by Sebastian+Leif\n")
output.write("# declination     z      flux       observed" + "\n")

for i in range(0, len(redshift_list)):
    output.write(str(declin[i]) + " " + str(z[i]) + " " + str(flux[i]) + " " + str(obs[0][i]) + "\n")

output.write("# E^2 dNdE = " + str(TotalFlux/(4*np.pi)) + "\n")

if options.histogram == True:
    histofreq, histobin = np.histogram(obs[0], bins=int(obs[0].max())+1, range=(0, obs[0].max()+1)) 

print ("RESULTS")
print ('E^2 dNdE = ' + str(TotalFlux/(4*np.pi)))
if options.histogram == True:
    print ('Event Distribution')
    for i, j in zip(histobin, histofreq):
        print(str(i)+"  "+str(j))
    print('-END-')
else:
    print ('Total no. of event is: '+str(np.sum(obs)))

################################
#Find out if there are sources that exceed IceCube's limit
#
if (options.NoPSComparison==False):
    sens = np.load('sens.npy')
    sensspline = InterpolatedUnivariateSpline(sens['dec'], sens['2'], k=5, ext=3)
    fluxToPointSourceSentivity = flux/(sensspline(np.arcsin(sinDec))*1e3)  #1e3 to change from TeV to GeV
    detectable  = [[i, j] for i, j in zip(fluxToPointSourceSentivity, declin) if i >= 1.]
    print ("Detectable sources are: ")
    print detectable
    output.write("# Fluxes exceeding Point Source limits " + str(detectable) + "\n")

if (options.NoHAWC==False):
    hawc_output = open(outputdir + "/HAWC/hawc_" + options.filename,"w")
    detectable = ([[i, j, k] for i, j, k in zip(flux, declin, redshift_list) if j>-26. and j < 64. and k<0.1])
    for i in range(0,len(detectable)):
        hawc_output.write('%.3e %.3f %.3f\n' % (detectable[i][0], detectable[i][1], detectable[i][2]))
    hawc_output.close()

if (options.NoCTA==False):
    cta_output = open(outputdir + "/CTA/cta_" + options.filename,"w")
    for i in range(0, len(redshift_list)):
        if obs[0][i]>0:
            cta_output.write( str(declin[i]) + " " + str(z[i]) + " " + str(flux[i]) + " " + str(obs[0][i]) + "\n")
    cta_output.close()
    
output.close()

