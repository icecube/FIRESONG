#!/usr/bin/python
# Authors: Theo Glauch
#

# General imports
from __future__ import division
import argparse
# Numpy / Scipy
import numpy as np
# Firesong code
from Evolution import get_evolution, SourcePopulation, cosmology
from Luminosity import get_LuminosityFunction
from input_output import get_outputdir


def flux_pdf(outputdir,
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
             LumMin=1e45, LumMax=1e54, nLbins=120,
             logFMin=-10, logFMax=6, nFluxBins=200,
             with_dFdz=False):

    if Transient:
        population = TransientSourcePopulation(cosmology,
                                               get_evolution(Evolution),
                                               timescale=timescale)
    else:
        population = SourcePopulation(cosmology, get_evolution(Evolution))

    N_sample = population.Nsources(density, zmax)

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
                 mode=" - Calculating Flux PDF")

    # Setup Arrays
    int_norm = population.RedshiftIntegral(zmax)

    zs = np.linspace(zmin, zmax, bins)
    deltaz = float(zmax-zmin)/bins

    Ls = np.linspace(np.log10(LumMin), np.log10(LumMax), nLbins)
    deltaL = (np.log10(LumMax)-np.log10(LumMin))/nLbins

    logFlux_array = np.linspace(logFMin, logFMax, nFluxBins)
    deltaLogFlux = float(logFMax-logFMin) / nFluxBins

    Count_array = np.zeros(nFluxBins)
    fluxOutOfBounds = []
    if with_dFdz:
        Flux_from_fixed_z = np.zeros(bins)

    # Integration
    # Loop over redshift bins
    for i, z in enumerate(zs):
        tot_flux_from_z = 0.
        # Loop over Luminosity bins
        for lum in Ls:
            # Number of Sources in
            dN = N_sample * luminosity_function.pdf(10**lum) * deltaL * \
                (population.RedshiftDistribution(z)/int_norm) * deltaz

            # Flux to Source Strength
            logF = np.log10(population.Lumi2Flux(10**lum,
                                                 index=index,
                                                 emin=emin,
                                                 emax=emax,
                                                 z=z))

            # Add dN to Histogram
            if logF < logFMax and logF > logFMin:
                idx = int((logF-logFMin) / deltaLogFlux)
                Count_array[idx] += dN
                if with_dFdz:
                    tot_flux_from_z += dN*10**logF
            else:
                fluxOutOfBounds.append(logF)
        if with_dFdz:
            Flux_from_fixed_z[i] = tot_flux_from_z

    if with_dFdz:
        return logFlux_array, Count_array, zs, Flux_from_fixed_z
    return logFlux_array, Count_array


if __name__ == "__main__":

    outputdir = get_outputdir()

    # Process command line options
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
                        help="Source evolution options:  HB2006SFR (default), NoEvolution")
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
    options = parser.parse_args()

    output = flux_pdf(luminosity=1e50, LF="LG", sigma=1,
                      index=2.19, emin=1e4, emax=1e7,
                      density=1e-7, Evolution="HB2006SFR",
                      zmin=0.005, zmax=10., bins=120,
                      LumMin=1e45, LumMax=1e54, nLbins=120,
                      logFMin=-10, logFMax=6, nFluxBins=200)
