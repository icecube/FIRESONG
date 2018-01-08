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
from Luminosity import LuminosityFunction

def firesong_simulation(options, outputdir):
    if re.search('.gz$', options.filename):
        output = gzip.open(outputdir+str(options.filename), 'wb')
    else:
        output = open(outputdir+str(options.filename),"w")

    if (options.zNEAR>0):
        if re.search('.gz$', options.filename):
            near_output = gzip.open(outputdir + "Near_" + options.filename,"w")
        else:
            near_output = open(outputdir + "Near_" + options.filename,"w")

    #Calculate total number of sources in the universe, and the flux from each source
    N_sample, candleflux = StandardCandleSources(options)
    ## Integrate[EdN/dE, {E, 10TeV, 10PeV}] * 4*Pi * dL1^2 * unit conversion
    luminosity = candleflux * (1.e-5) * scipy.integrate.quad(lambda E: 2.**(options.index-2)*(E/1.e5)**(-options.index+1), 1.e4, 1.e7)[0] * 4*np.pi * (LuminosityDistance(1.)*3.086e24)**2. *50526
    ## the candle flux of Transient mode is the candle fluence, so we need to convert it to flux first to obtain the luminosity
    if options.Transient == True:
        luminosity = luminosity / options.timescale
    ## If luminosity of the sources is specified, re-calculate candleflux
    if options.luminosity != 0.0:
        candleflux = LtoFlux(options)
        ## Need to change from flux to fluence for transient mode
        if options.Transient == True:
            candleflux = candleflux*options.timescale
        luminosity = options.luminosity
    flux_z1 = LuminosityFunction(options,N_sample,candleflux)


    print ("##############################################################################")
    print ("##### FIRESONG initializing #####")
    if options.LF == "SC":
        print ("Standard candle sources")
    if options.LF == "LG":
        print ("Lognormal distributed sources")
    if options.LF == "PL":
        print ("PowerLaw distributed sources")
    if options.Transient == True:
        print ("Transient Sources, time scale = " + str(options.timescale) + "s")
    print ("Source evolution assumed: " + str(options.Evolution))
    print ("Local density of neutrino sources: " + str(options.density) + "/Mpc^3")
    print ("Total number of neutrinos sources in the Universe: " + str(N_sample))
    if options.luminosity == 0.0:
        print ("Desired neutrino diffuse flux: E^2 dN/dE = " + str(options.fluxnorm) + " (E/100 TeV)^(" + str(-(options.index-2.)) + ") GeV/cm^2.s.sr")
    print ("Redshift range: 0 - " + str(options.zmax)) 
    print ("Standard Candle Luminosity: {:.4e} erg/yr".format(luminosity))
    print ("##### FIRESONG initialization done #####")

    ##################################################
    #        Simulation starts here
    ##################################################

    output.write("# FIRESONG Output description\n")
    if options.LF == "SC":
        output.write("# Standard candle sources\n")
    if options.LF == "LG":
        output.write("# Lognormal distributed sources\n")
    if options.LF == "PL":
        output.write("# PowerLaw distributed sources\n")
    if options.Transient == True:
        output.write("# Transient Sources, time scale = " + str(options.timescale) + "s \n")
    output.write("# Desired neutrino diffuse flux:\n")
    output.write("#      E^2 dN_{diffuse}/dE = " + str(options.fluxnorm) + " (E/100 TeV)^(" + str(-(options.index-2.)) + ") [GeV/cm^2.s.sr]\n") 
    output.write("# Neutrino point source fluxes listed below are of \'A\' where the flux is:\n")
    output.write("#      E^2 dN_{PS}/dE = A * (E/100 TeV)^(" + str(-(options.index-2.)) + ") [GeV/cm^2.s]\n") 
    output.write("# Standard Candle Luminosity: {:.4e} erg/yr \n".format(luminosity))
    output.write("# Note that using 4 years, IceCube sensitivity in the northern hemisphere\n")
    output.write("# is approximately 10^-9 in the units used for A\n")
    output.write("# Dec(deg) Redshift A\n")


    # Luminosity distace for z=1. Internally, fluxes are scaled to this distance.
    dL1 = LuminosityDistance(1.)
    # Generate a histogram to store redshifts. Starts at z = 0.0005 and increases in steps of 0.001
    redshift_bins = np.arange(0.0005,options.zmax, options.zmax/10000.)

    # RedshiftCDF is used for inverse transform sampling
    RedshiftPDF = [RedshiftDistribution(redshift_bins[i], options) for i in range(0,len(redshift_bins))]
    RedshiftCDF = np.cumsum(RedshiftPDF)
    RedshiftCDF = RedshiftCDF / RedshiftCDF[-1]

    TotalFlux = 0

    for i in range(0,N_sample):
        # Generate a random redshift using inverse transform sampling
        test = np.random.rand()
        bin_index = np.searchsorted(RedshiftCDF, test)
        z = redshift_bins[bin_index]
        # Random declination over the entire sky
        sinDec = 2*np.random.rand() -1
        declin = 180*np.arcsin(sinDec)/np.pi
        dL = LuminosityDistance(z)
        ## IMPORTANT notice, in the following "flux" means fluence in Transient mode, but flux in steady source mode, until TotalFlux(TotalFluence)
        ## is calculated
        if options.LF != 'SC':
            flux = flux_z1[i] * (dL1*dL1)/(dL*dL) * ((1.+z)/2.)**(-options.index+2)
        else:
            flux = flux_z1 * (dL1*dL1)/(dL*dL) * ((1.+z)/2.)**(-options.index+2)
        if options.Transient == True:
            flux = flux*(1.+z)/2.
        TotalFlux = TotalFlux + flux
        # For transient sources, the flux measured on Earth will be red-shifted-fluence/{(1+z)*burst duration} 
        if options.Transient == True:
            flux = flux / ((1.+z)*options.timescale)
        output.write('{:.4f} {:.4f} {:.4e}\n'.format(declin, z, flux))
        if i%100000==0 and i>0:
            print "Generated ", i, " neutrino sources"
    ############# This is the place to plug in Detector output modules ############
        if (z<options.zNEAR):
            near_output.write('{:.4e} {:.4f} {:.4f}\n'.format(flux,declin, z))

    #For transient source, we calculate the total fluence from all sources, then obtain the diffuse flux by doing a time average over a year
    if options.Transient == True:
        TotalFlux = TotalFlux / (86400*365)
    output.write("# E^2 dNdE = " + str(TotalFlux/(4*np.pi)) + "\n")
    print "Actual diffuse flux simulated :  E^2 dNdE = " + str(TotalFlux/(4*np.pi)) + " (E/100 TeV)^(" + str(-(options.index-2.)) + ") [GeV/cm^2.s.sr]" 
    output.close()
    if (options.zNEAR>0):
        near_output.close()

def firesong_simulation_by_hand(outputdir,
                                filename   ='Firesong.out',
                                density    = 1e-9,
                                Evolution  = "HB2006SFR", 
                                Transient  = False,
                                timescale  = 1000.,
                                zmax       = 10.,
                                fluxnorm   = 0.9e-8,
                                index      = 2.13,
                                LF         = "SC",
                                sigma      = 1.0,
                                zNEAR      = -1,
                                luminosity = 0.0):
    
    options = argparse.Namespace( filename   = filename,
                                  density    = density, 
                                  Evolution  = Evolution, 
                                  Transient  = Transient, 
                                  timescale  = timescale, 
                                  zmax       = zmax, 
                                  fluxnorm   = fluxnorm, 
                                  index      = index, 
                                  LF         = LF, 
                                  sigma      = sigma, 
                                  zNEAR      = zNEAR, 
                                  luminosity = luminosity)
    firesong_simulation(options, outputdir)

if __name__ == "__main__":

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
    parser.add_argument("--zNEAR", action="store",dest="zNEAR", type=float,
                        default=-1, help="Write down a sepaarate file for sources closer than specified redshift. If nothing is specfied, no file is written.")
    parser.add_argument("--L", action="store",
                        dest="luminosity", type=float, default=0.0,
                        help="Set luminosity for each source, will reset fluxnorm option, unit erg/yr")

    options = parser.parse_args()
    
    firesong_simulation(options, outputdir)
