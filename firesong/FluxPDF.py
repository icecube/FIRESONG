#!/usr/bin/python
# Authors: Theo Glauch

"""Calculates the PDF of neutrino fluxes for a given set of parameters"""

# General imports
from __future__ import division
import argparse
# Numpy / Scipy
import numpy as np
# Firesong code
from firesong.Evolution import get_evolution, SourcePopulation, cosmology
from firesong.Evolution import TransientSourcePopulation
from firesong.Luminosity import get_LuminosityFunction
from firesong.input_output import output_writer_PDF, get_outputdir, print_config


def flux_pdf(outputdir,
             filename=None,
             density=1e-9,
             Evolution="MD2014SFR",
             Transient=False,
             timescale=1000.,
             zmin=0.0005,
             zmax=10.,
             bins=10000,
             fluxnorm=1.44e-8,
             index=2.28,
             LF="SC",
             sigma=1.0,
             luminosity=0.0,
             emin=1e4,
             emax=1e7,
             LumMin=1e45, LumMax=1e54, nLbins=500,
             logFMin=-10, logFMax=6, nFluxBins=200,
             with_dFdz=False,
             verbose=True):
    """
    Simulate a universe of neutrino sources and calculate the PDF of fluxes

    Args:
        outputdir (str or None): path to write output. If None, return results
            without writing a file
        filename (str or None): name of the output file. If None, return 
            results without writing a file
        density (float, optional, default=1e-9): local density of neutrino
            sources. Units of Mpc^-3 (Mpc^-3 yr^-1 if Transient=True)
        Evolution (str, optional, default='MD2014SFR'): Redshift evolution model
            assumed. Options include 'MD2014SFR', 'HB2006SFR', 'YMKBH2008SFR',
            'CC2015SNR', and 'NoEvolution'
        Transient (bool, optional, default=False): If true, simulate 
            transient neutrino sources instead of steady sources
        timescale (bool, optional, default=1000): Timescale in seconds
            for transient sources
        zmin (float, optional, default=0.0005): Closest redshift to consider
        zmax (float, optional, default=10.): Farthest redshift to consider
        bins (int, optional, default=1000): Number of bins used when creating
            the redshift PDF
        fluxnorm (float, optional, default=0.9e-8): Normalization on the total
            astrophysical diffuse flux, E^2d Phi/dE. Units of GeV s^-1 sr^-1
        index (float, optional, default=2.13): Spectral index of diffuse flux
        LF (string, optional, default="SC"): Luminosity function, choose
            between standard candle (SC), LogNormal (LG)
        sigma (float, optional, default=1.0): Width of lognormal distribution
            if LF="LG"
        luminosity (float, optional, default=0.0): Manually fix the 
            luminosity of sources if not equal to 0. Overrides fluxnorm. 
            Units of erg/yr
        emin (float, optional, default=1e4): Minimum neutrino energy in GeV
        emax (float, optional, default=1e7): Maximum neutrino energy in GeV
        LumMin (float, optional, default=1e45): Minimum luminosity considered
            in UNITS
        LumMax (float, optional, default=1e54): Max luminosity considered,
            in UNITS
        nLbins (int, optional, default=120): Number of luminosity bins
        logFMin (float, optional, default=-10.): minimum log10 of flux to 
            consider, in UNITS
        logFMax (float, optional, default=6.): max log10 of flux to 
            consider, in UNITS
        nFluxBins (int, optional, default=200): Number of flux bins
        with_dFdz (bool, optional, default=False): Include the derivative
            of the flux pdf vs. redshift
        verbose (bool, optional, default=True): print simulation paramaters
            if True else suppress printed output

    Returns:
        tuple of arrays: Return lists of the fluxes and their counts,
            equivalent to the counts in a histogram with fluxes as bin-edges.
            Choice to include derivative of the PDF as well if with_dFdz
    """

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
    if verbose:
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

    fluxOutOfBounds = []
    if with_dFdz:
        Flux_from_fixed_z = np.zeros(bins)

    # bin-by-bin integration
    if LF == "SC":
        pdf_dL = np.where(np.abs(Ls - np.log10(luminosity_function.mean))
            < deltaL/2., 1.0, 0.0)
    else:
        prepend_val = 10.**(Ls[0] - deltaL) # ensures shapes match
        pdf_dL = luminosity_function.pdf(10**Ls) * \
            np.diff(10.**Ls, prepend=prepend_val)

    pdf_dz = population.RedshiftDistribution(zs) * deltaz / int_norm
    dNs = N_sample * pdf_dz[:,np.newaxis] * pdf_dL

    # flux scales linearly with luminosity for a fixed redshift,
    # so just calculate the z dependence and scale by luminosity after
    flux_from_redshift = population.Lumi2Flux(1.,
            index=index,
            emin=emin,
            emax=emax,
            z=zs)
    fluxes_all = flux_from_redshift[:,np.newaxis] * (10.**Ls)
    logFs = np.log10(fluxes_all)
    
    # Find which flux bins the fluxes fall into
    indices = np.digitize(logFs, bins=logFlux_array) - 1
    mask = (indices >= 0) & (indices < len(logFlux_array) - 1)

    # Calculate the sources out of bounds and the sums for each flux
    fluxOutOfBounds = logFs[~mask]
    Count_array = np.array([np.sum(dNs[indices == ind]) for ind in range(nFluxBins)])

    # Project the sum to isolate z dependence
    Flux_from_fixed_z = np.sum(dNs * (10.**logFs), axis=1)

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
                        help='Simulate transient sources')
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
                        dest="lBins", type=float, default=500,
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
