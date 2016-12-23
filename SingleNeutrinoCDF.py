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
from Evolution import RedshiftDistribution, StandardCandleSources, LuminosityDistance
from Luminosity import LuminosityFunction
import hawc

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

options = parser.parse_args()
if re.search('.gz$', options.filename):
    output = gzip.open(outputdir+str(options.filename), 'wb')
else:
    output = open(outputdir+str(options.filename),"w")

N_sample, candleflux = StandardCandleSources(options)
flux_z1 = LuminosityFunction(options,N_sample,candleflux)

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
print ("##### FIRESONG initialization done #####")

##################################################
#        Simulation starts here
##################################################

output.write("# FIRESONG Output description\n")
output.write("# Desired neutrino diffuse flux:\n")
output.write("#      E^2 dN_{diffuse}/dE = " + str(options.fluxnorm) + " (E/100 TeV)^(" + str(-(options.index-2.)) + ") [GeV/cm^2.s.sr]\n") 
output.write("# Neutrino point source fluxes listed below are of \'A\' where the flux is:\n")
output.write("#      E^2 dN_{PS}/dE = A * (E/100 TeV)^(" + str(-(options.index-2.)) + ") [GeV/cm^2.s.sr]\n") 
output.write("# Note that using 4 years, IceCube sensitivity in the northern hemisphere\n")
output.write("# is approximately 10^-9 in the units used for A\n")
output.write("# Dec(deg) Redshift A\n")

# Luminosity distace for z=1. Internally, fluxes are scaled to this distance.
dL1 = LuminosityDistance(1.)

# Generate a histogram to store redshifts. Starts at z = 0.0005 and increases in steps of 0.001
redshift_bins = np.arange(0.0005,options.zmax, 0.001)

NeutrinoPDF = [RedshiftDistribution(redshift_bins[i], options)*flux_z1*(dL1*dL1)/(LuminosityDistance(redshift_bins[i])*LuminosityDistance(redshift_bins[i])) for i in range(0,len(redshift_bins))]
NeutrinoCDF = np.cumsum(NeutrinoPDF)
NeutrinoCDF = NeutrinoCDF / NeutrinoCDF[-1]
    
options = parser.parse_args()
if re.search('.gz$', options.filename):
    output = gzip.open(outputdir+str(options.filename), 'wb')
else:
    output = open(outputdir+str(options.filename),"w")

for i in range(0,len(NeutrinoCDF)):
    z = redshift_bins[i] 
    dL = LuminosityDistance(z)
    flux = flux_z1 * (dL1*dL1)/(dL*dL)
    output.write('{:.4f} {:.6e} {:.6e}\n'.format(z,flux,NeutrinoCDF[i]))
output.close()
