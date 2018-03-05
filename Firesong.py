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
# Firesong code
from Evolution import get_evolution, SourcePopulation
from Evolution import TransientSourcePopulation, cosmology, Simulation
from Luminosity import get_LuminosityFunction
from input_output import output_writer, print_str

def firesong_simulation(options, outputdir):
    emin = 1e4
    emax = 1e7

    if options.Transient is True:
        population = TransientSourcePopulation(cosmology,
                                               get_evolution(options.Evolution),
                                               timescale=options.timescale)
    else:
        population = SourcePopulation(cosmology,
                                      get_evolution(options.Evolution))

    N_sample = int(population.Nsources(options.density, options.zmax))

    if options.luminosity == 0.0:
        ## If luminosity not specified calculate luminosity from diffuse flux
        luminosity = population.StandardCandleLuminosity(options.fluxnorm,
                                                         options.density,
                                                         options.zmax,
                                                         options.index,
                                                         emin=emin,
                                                         emax=emax)
    else:
        luminosity = options.luminosity

    delta_gamma = 2-options.index
    print_str(options.LF, options.Transient, options.timescale,
              options.Evolution, options.density, N_sample,
              options.luminosity, options.fluxnorm, delta_gamma,
              options.zmax, luminosity)

    ##################################################
    #        Simulation starts here
    ##################################################

    simulation = Simulation(population,
                            get_LuminosityFunction(options, luminosity),
                            index=options.index,
                            zmax=options.zmax,
                            emin=emin,
                            emax=emax)

    out = output_writer(outputdir, options.filename)
    out.write_header(options.LF, options.Transient, options.timescale,
                     options.fluxnorm, delta_gamma, options.luminosity)

    # CHECK loop can be more efficient
    TotalFlux = 0
    for i in range(0, N_sample):
        # IMPORTANT notice, in the following "flux" means fluence in
        # Transient mode, but flux in steady source mode,
        # until TotalFlux(TotalFluence) is calculated
        flux, z, declin = simulation.sample_flux()
        TotalFlux = TotalFlux + flux

        # For transient sources, the flux measured on Earth will be
        # red-shifted-fluence/{(1+z)*burst duration}
        if options.Transient:
            flux = population.fluence2flux(flux, z)

        out.write(declin, z, flux)
        if i % 100000 == 0 and i > 0:
            print "Generated ", i, " neutrino sources"

    # For transient source, we calculate the total fluence from all sources,
    # then obtain the diffuse flux by doing a time average over a year
    if options.Transient:
        TotalFlux = TotalFlux / population.yr2sec
    TotalFlux /= 4*np.pi  # give in per sr

    out.finish(TotalFlux)
    print "Actual diffuse flux simulated :  E^2 dNdE = {TotalFlux} (E/100 TeV)^({delta_gamma}) [GeV/cm^2.s.sr]".format(**locals())


def firesong_simulation_by_hand(outputdir,
                                filename='Firesong.out',
                                density=1e-9,
                                Evolution="HB2006SFR",
                                Transient=False,
                                timescale=1000.,
                                zmax=10.,
                                fluxnorm=0.9e-8,
                                index=2.13,
                                LF="SC",
                                sigma=1.0,
                                zNEAR=-1,
                                luminosity=0.0):

    options = argparse.Namespace(filename=filename,
                                 density=density,
                                 Evolution=Evolution,
                                 Transient=Transient,
                                 timescale=timescale,
                                 zmax=zmax,
                                 fluxnorm=fluxnorm,
                                 index=index,
                                 LF=LF,
                                 sigma=sigma,
                                 zNEAR=zNEAR,
                                 luminosity=luminosity)
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
    parser.add_argument('-o', action='store',
                        dest='filename',
                        default='Firesong.out',
                        help='Output filename')
    parser.add_argument('-d', action='store', dest='density',
                        type=float, default=1e-9,
                        help='Local neutrino source density [1/Mpc^3]')
    parser.add_argument("--evolution", action="store",
                        dest="Evolution", default='HB2006SFR',
                        help="Source evolution options:  HB2006SFR (default),  NoEvolution")
    parser.add_argument("--transient", action='store_true',
                        dest='Transient', default=False,
                        help='Simulate transient sources, NOT TESTED YET!')
    parser.add_argument("--timescale", action='store',
                        dest='timescale', type=float,
                        default=1000., help='time scale of transient sources, default is 1000sec.')
    parser.add_argument("--zmax", action="store", type=float,
                        dest="zmax", default=10.,
                        help="Highest redshift to be simulated")
    parser.add_argument("--fluxnorm", action="store", dest='fluxnorm',
                        type=float, default=0.9e-8,
                        help="Astrophysical neutrino flux normalization A on E^2 dN/dE = A (E/100 TeV)^(-index+2) GeV/cm^2.s.sr")
    parser.add_argument("--index", action="store", dest='index',
                        type=float, default=2.13,
                        help="Astrophysical neutrino spectral index on E^2 dN/dE = A (E/100 TeV)^(-index+2) GeV/cm^2.s.sr")
    parser.add_argument("--LF", action="store", dest="LF", default="SC",
                        help="Luminosity function, SC for standard candles, LG for lognormal, PL for powerlaw")
    parser.add_argument("--sigma", action="store",
                        dest="sigma", type=float, default=1.0,
                        help="Width of a log normal Luminosity function in dex, default: 1.0")
    parser.add_argument("--zNEAR", action="store", dest="zNEAR",
                        type=float, default=-1,
                        help="Write down a sepaarate file for sources closer than specified redshift. If nothing is specfied, no file is written.")
    parser.add_argument("--L", action="store",
                        dest="luminosity", type=float, default=0.0,
                        help="Set luminosity for each source, will reset fluxnorm option, unit erg/yr")

    options = parser.parse_args()

    firesong_simulation(options, outputdir)
