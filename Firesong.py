#!/usr/bin/python
# Authors: Chris Tung
#          Ignacio Taboada
#

# General imports
#from __future__ import division
import argparse
# Numpy / Scipy
import numpy as np
# Firesong code
from Evolution import get_evolution, SourcePopulation
from Evolution import TransientSourcePopulation, cosmology
from Luminosity import get_LuminosityFunction
from input_output import output_writer, print_config, get_outputdir
from sampling import InverseCDF

def firesong_simulation(outputdir,
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
                        seed=None,
                        zNEAR=-1):

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

    rng = np.random.RandomState(seed)

    redshift_bins = np.arange(zmin, zmax, zmax/float(bins))
    RedshiftPDF = population.RedshiftDistribution(redshift_bins)
    invCDF = InverseCDF(redshift_bins, RedshiftPDF)

    out = output_writer(outputdir, filename)
    out.write_header(LF, Transient, timescale, fluxnorm,
                     delta_gamma, luminosity)

    # IMPORTANT notice, in the following "flux" means fluence in
    # Transient mode, but flux in steady source mode,
    # until TotalFlux(TotalFluence) is calculated
    
    # sample source
    zs = invCDF(rng.uniform(low=0.0, high=1.0, size = N_sample))
    lumis = luminosity_function.sample_distribution(nsources=N_sample, rng=rng)
    if np.ndim(lumis) < 1:
        lumis = np.array([lumis]*N_sample)
    fluxes = population.Lumi2Flux(lumis, index, emin, emax, zs)
    # Random declination over the entire sky
    sinDecs = rng.uniform(-1, 1, size=N_sample)
    declins = np.degrees(np.arcsin(sinDecs))

    TotalFlux = np.sum(fluxes)

    # For transient sources, the flux measured on Earth will be
    # red-shifted-fluence/{(1+z)*burst duration}
    if Transient:
        fluxes = population.fluence2flux(fluxes, zs)

    out.write(declins, zs, fluxes)

    # For transient source, we calculate the total fluence from all sources,
    # then obtain the diffuse flux by doing a time average over a year
    if Transient:
        TotalFlux = TotalFlux / population.yr2sec
    TotalFlux /= 4*np.pi  # give in per sr

    out.finish(TotalFlux)
    print("Actual diffuse flux simulated :")
    log = "E^2 dNdE = {TotalFlux} (E/100 TeV)^({delta_gamma}) [GeV/cm^2.s.sr]"
    print(log.format(**locals()))

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
                        help="Source evolution options:  HB2006SFR (default), YMKBH2008SFR, CC2015SNR , MD2014SFR, NoEvolution")
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
    parser.add_argument("--zNEAR", action="store", dest="zNEAR",
                        type=float, default=-1,
                        help="Write down a sepaarate file for sources closer than specified redshift. If nothing is specfied, no file is written.")
    options = parser.parse_args()

    firesong_simulation(outputdir,
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
                        zNEAR=options.zNEAR)
