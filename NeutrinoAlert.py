#!/usr/bin/python
# Authors: Chris Tung
#          Igncio Taboada
#

# General imports
from __future__ import division
import argparse
# Numpy / Scipy
import numpy as np
# Firesong code
from Evolution import get_evolution, SourcePopulation
from Evolution import TransientSourcePopulation, cosmology
from Luminosity import get_LuminosityFunction
from input_output import output_writer, print_config, get_outputdir
from sampling import InverseCDF


def calc_NeutrinoAlert(outputdir,
                       filename='Firesong.out',
                       density=1e-9,
                       Evolution="HB2006SFR",
                       Transient=False,
                       timescale=1000.,
                       zmin=0.0005,
                       zmax=10.,
                       bins=10000,
                       fluxnorm=0.9e-8,
                       index=2.13,
                       LF="SC",
                       sigma=1.0,
                       luminosity=0.0,
                       emin=1e4,
                       emax=1e7,
                       AlertNumber=1):
    if Transient:
        population = TransientSourcePopulation(cosmology,
                                               get_evolution(Evolution),
                                               timescale=timescale)
    else:
        population = SourcePopulation(cosmology, get_evolution(Evolution))

    N_sample = int(population.Nsources(density, zmax))

    if luminosity == 0.0:
        ## If luminosity not specified calculate luminosity from diffuse flux
        luminosity = population.StandardCandleLuminosity(fluxnorm,
                                                         density,
                                                         zmax,
                                                         index,
                                                         emin=emin,
                                                         emax=emax)

    luminosity_function = get_LuminosityFunction(luminosity, LF=LF,
                                                 sigma=sigma)

    delta_gamma = 2-index
    print_config(LF, Transient, timescale, Evolution, density, N_sample,
                 luminosity, fluxnorm, delta_gamma, zmax, luminosity,
                 mode=" - Calculating Neutrino CDFs ")

    ##################################################
    #        Simulation starts here
    ##################################################
    out = output_writer(outputdir, filename)
    out.write_header(LF, Transient, timescale, fluxnorm, delta_gamma, luminosity)

    # Generate a histogram to store redshifts.
    # Default: Starts at z = 0.0005 and increases in steps of 0.001
    redshift_bins = np.arange(zmin, zmax, zmax/float(bins))

    z0 = 1
    # Calculate the redshift z PDF for neutrino events
    if Transient:
        NeutrinoPDF_z = [population.RedshiftDistribution(z)*((1+z)/(1+z0))**(-index+3)/(population.LuminosityDistance(z)**2.) for z in redshift_bins]
    else:
        NeutrinoPDF_z = [population.RedshiftDistribution(z)*((1+z)/(1+z0))**(-index+2)/(population.LuminosityDistance(z)**2.) for z in redshift_bins]
    invCDF = InverseCDF(redshift_bins, NeutrinoPDF_z)

    for i in range(0, options.AlertNumber):
        z = invCDF.sample()
        # Random declination over the entire sky
        sinDec = 2 * np.random.rand() - 1
        declin = np.degrees(np.arcsin(sinDec))
        lumi = luminosity_function.sample_distribution()
        flux = population.Lumi2Flux(lumi, index, emin, emax, z)
        out.write(declin, z, flux)
    out.finish()

if __name__ == "__main__":
    outputdir = get_outputdir()

    # Process command line options
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', action='store', dest='filename',
                        default='Firesong.out', help='Output filename')
    parser.add_argument('-d', action='store', dest='density',
                        type=float, default=1e-9,
                        help='Local neutrino source density [1/Mpc^3]')
    parser.add_argument("--evolution", action="store",
                        dest="Evolution", default='HB2006SFR',
                        help="Source evolution options:  HB2006SFR (default), NoEvolution")
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
    parser.add_argument("--L", action="store",
                        dest="luminosity", type=float, default=0.0,
                        help="Set luminosity for each source, will reset fluxnorm option, unit erg/yr")
    parser.add_argument('-N', action='store', dest='AlertNumber',
                        type=int, default=1,
                        help='Number of neutrinos to generate')
    options = parser.parse_args()

    calc_NeutrinoAlert(outputdir,
                       filename=options.filename,
                       density=options.density,
                       Evolution=options.Evolution,
                       Transient=options.Transient,
                       timescale=options.timescale,
                       zmax=options.zmax,
                       fluxnorm=options.fluxnorm,
                       index=options.index,
                       LF=options.LF,
                       sigma=options.sigma,
                       luminosity=options.luminosity,
                       AlertNumber=options.AlertNumber)
