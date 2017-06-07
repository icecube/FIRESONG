#!/usr/bin/python
# Authors: Chris Tung
#          Igncio Taboada
#

# General imports
from __future__ import division
import os
import gzip
import re
import argparse
# Numpy / Scipy
import numpy as np
import scipy.integrate
# Firesong code
from Evolution import RedshiftDistribution, StandardCandleSources, LuminosityDistance, LtoFlux
from Luminosity import LuminosityFunction, LuminosityPDF


#
# Check that the Firesong environmental variable is set
# This is needed for output and to read exposures, effective areas, etc
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
parser.add_argument('-N', action='store', dest='AlertNumber',type=int,default= 1,
                    help='Number of neutrinos to generate')
parser.add_argument('-o', action='store', dest='filename',default= 'Firesong.out',
                    help='Output filename')
parser.add_argument('-d', action='store', dest='density', type=float, default = 1e-9,
                    help='Local neutrino source density [1/Mpc^3]')
parser.add_argument("--evolution", action="store",
                    dest="Evolution", default='HB2006SFR',
                    help="Source evolution options:  HB2006SFR (default),  NoEvolution")
parser.add_argument("--transient", action='store_true',
                    dest='Transient', default=False,
                    help='Simulate transient sources, NOT TESTED YET!')
parser.add_argument("--timescale", action='store', dest='timescale', type=float,
                    default=1000., help='time scale of transient sources, default is 1000sec.')
parser.add_argument("--zmax", action="store", type=float,
                    dest="zmax", default=10.,
                    help="Highest redshift to be simulated")
parser.add_argument("--fluxnorm", action="store", dest='fluxnorm', type=float, default=0.9e-8,
                    help="Astrophysical neutrino flux normalization A on E^2 dN/dE = A (E/100 TeV)^(-index+2) GeV/cm^2.s.sr")
parser.add_argument("--index", action="store", dest='index', type=float, default=2.13,
                    help="Astrophysical neutrino spectral index on E^2 dN/dE = A (E/100 TeV)^(-index+2) GeV/cm^2.s.sr")
parser.add_argument("--LF",action="store", dest="LF",default="SC",                    
                    help="Luminosity function, SC for standard candles, LG for lognormal, PL for powerlaw")
parser.add_argument("--sigma", action="store",
                    dest="sigma", type=float, default=1.0,
                    help="Width of a log normal Luminosity function in dex, default: 1.0")
parser.add_argument("--L", action="store",
                    dest="luminosity", type=float, default=0.0,
                    help="Set luminosity for each source, will reset fluxnorm option, unit erg/yr")

options = parser.parse_args()
if re.search('.gz$', options.filename):
    output = gzip.open(outputdir+str(options.filename), 'wb')
else:
    output = open(outputdir+str(options.filename),"w")

N_sample, candleflux = StandardCandleSources(options)
## Integrate[EdN/dE, {E, 10TeV, 10PeV}] * 4*Pi * dL1^2 * unit conversion
luminosity = candleflux * (1.e-5) * scipy.integrate.quad(lambda E: 2.**(-options.index+2)*(E/1.e5)**(-options.index+1), 1.e4, 1.e7)[0] * 4*np.pi * (LuminosityDistance(1.)*3.086e24)**2. *50526
## If luminosity of the sources is specified, re-calculate candleflux
if options.luminosity != 0.0:
    candleflux = LtoFlux(options)
    luminosity = options.luminosity
flux_z1 = LuminosityFunction(options, N_sample, candleflux)

print ("##############################################################################")
print ("##### FIRESONG initializing - Calculating Neutrino CDFs #####")
if options.LF == "SC":
    print ("Standard candle sources")
if options.LF == "LG":
    print ("Lognormal distributed sources")
if options.LF == "PL":
    print ("PowerLaw distributed sources")
print ("Source evolution assumed: " + str(options.Evolution))
print ("Local density of neutrino sources: " + str(options.density) + "/Mpc^3")
print ("Total number of neutrinos sources in the Universe: " + str(N_sample))
print ("Desired neutrino diffuse flux: E^2 dN/dE = " + str(options.fluxnorm) + " (E/100 TeV)^(" + str(-(options.index-2.)) + ") GeV/cm^2.s.sr")
print ("Redshift range: 0 - " + str(options.zmax)) 
print ("Standard Candle Luminosity: {:.4e} erg/yr".format(luminosity))
print ("##### FIRESONG initialization done #####")

##################################################
#        Simulation starts here
##################################################

output.write("# FIRESONG Output description\n")
output.write("# Desired neutrino diffuse flux:\n")
output.write("#      E^2 dN_{diffuse}/dE = " + str(options.fluxnorm) + " (E/100 TeV)^(" + str(-(options.index-2.)) + ") [GeV/cm^2.s.sr]\n") 
output.write("# Neutrino point source fluxes listed below are of \'A\' where the flux is:\n")
output.write("#      E^2 dN_{PS}/dE = A * (E/100 TeV)^(" + str(-(options.index-2.)) + ") [GeV/cm^2.s.sr]\n") 
output.write("# Standard Candle Luminosity: {:.4e} erg/yr \n".format(luminosity))
output.write("# Note that using 4 years, IceCube sensitivity in the northern hemisphere\n")
output.write("# is approximately 10^-9 in the units used for A\n")
output.write("# Dec(deg) Redshift A\n")

# Luminosity distace for z=1. Internally, fluxes are scaled to this distance.
dL1 = LuminosityDistance(1.)

# Generate a histogram to store redshifts. Starts at z = 0.0005 and increases in steps of 0.001
redshift_bins = np.arange(0.0005,options.zmax, 0.001)


# Calculate the redshift z PDF for neutrino events
if options.Transient == False:
    NeutrinoPDF_z = [RedshiftDistribution(z, options)*((1+z)/2.)**(-options.index+2)/(LuminosityDistance(z)**2.) for z in redshift_bins]
else:
    NeutrinoPDF = [RedshiftDistribution(redshift_bins[i], options)*((1+redshift_bins[i])/2.)**(-options.index+3)*flux_z1*(dL1**2.)/(LuminosityDistance(redshift_bins[i])**2.) for i in range(0,len(redshift_bins))]
NeutrinoCDF_z = np.cumsum(NeutrinoPDF_z)
NeutrinoCDF_z = NeutrinoCDF_z / NeutrinoCDF_z[-1]

# Obtain the flux_z1 PDF for neutrino event
if options.LF != "SC":
    f1_list, NeutrinoCDF_f = LuminosityPDF(options, candleflux)

for i in range(0,options.AlertNumber):
    # Random variates from the above constructed PDFs
    test = np.random.rand()
    bin_index_z = np.searchsorted(NeutrinoCDF_z, test)
    z = redshift_bins[bin_index_z]
    if options.LF != 'SC':
        test = np.random.rand()
        bin_index_f = np.searchsorted(NeutrinoCDF_f, test)
        flux_z1 = f1_list[bin_index_f]
    # Random declination over the entire sky
    sinDec = 2*np.random.rand() -1
    declin = 180*np.arcsin(sinDec)/np.pi
    dL = LuminosityDistance(z)
    flux = flux_z1 * (dL1 / dL)**2 * ((1+z)/2.)**(-options.index+2)
    if options.Transient == True:
        flux = flux/(options.timescale)
    output.write('{:.3f} {:.4f} {:.6e}\n'.format(declin, z, flux))
output.close()
