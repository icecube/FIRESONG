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
from input_output import output_writer_PDF, get_outputdir, print_config


def flux_pdf(outputdir,
             filename=None,
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

    if filename is None:
        if with_dFdz:
            return logFlux_array, Count_array, zs, Flux_from_fixed_z
        return logFlux_array, Count_array
    
    out = output_writer_PDF(outputdir, filename)
    out.write(logFlux_array, Count_array)
    out.finish()
    if with_dFdz:
        out = output_writer_PDF(outputdir, filename+".dFdz")
        out.write(zs, Flux_from_fixed_z)
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
    parser.add_argument("--eRange", action="store", nargs=2,
                        dest="eRange", type=float, default=[1e4, 1e7],
                        help="Set the minimal and maximal energy boundary of the measured fluxnorm, unit GeV")
    parser.add_argument("--zRange", action="store", nargs=2,
                        dest="zRange", type=float, default=[0.0005, 10.],
                        help="Set the minimal and maximal integration bounderis in redshift")
    parser.add_argument("--zBins", action="store",
                        dest="zBins", type=float, default=120,
                        help="Set number of z bins used in integration")
    parser.add_argument("--lRange", action="store", nargs=2,
                        dest="lRange", type=float, default=[1e45, 1e54],
                        help="Set the minimal and maximal integration bounderis in luminosity, unit erg/yr")
    parser.add_argument("--lBins", action="store",
                        dest="lBins", type=float, default=120,
                        help="Set number of log10(luminosity) bins used in integration")
    parser.add_argument("--fRange", action="store", nargs=2,
                        dest="fRange", type=float, default=[-10, 6],
                        help="Set the minimal and maximal log10(flux) range , unit log10(GeV/cm^2.s.sr)")
    parser.add_argument("--fBins", action="store",
                        dest="fBins", type=float, default=120,
                        help="Set number of log10(flux) bins used in evaluation")
    parser.add_argument("--dFdz", action="store_true", 
                        help="If ggiven a second file is created with dF/dz.")
    options = parser.parse_args()

    output = flux_pdf(outputdir,
                      filename=options.filename,
                      density=options.density,
                      Evolution=options.Evolution,
                      Transient=options.Transient,
                      timescale=options.timescale,
                      fluxnorm=options.fluxnorm,
                      index=options.index,
                      LF=options.LF,
                      sigma=options.sigma,
                      luminosity=options.luminosity,
                      emin=options.eRange[0],
                      emax=options.eRange[1],
                      zmin=options.zRange[0],
                      zmax=options.zRange[1],
                      bins=options.zBins,
                      LumMin=options.lRange[0],
                      LumMax=options.lRange[1],
                      nLbins=options.lBins,
                      logFMin=options.fRange[0],
                      logFMax=options.fRange[1],
                      nFluxBins=options.fBins,
                      with_dFdz=options.dFdz)
