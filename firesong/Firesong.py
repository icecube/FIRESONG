#!/usr/bin/python
# Authors: Chris Tung
#          Ignacio Taboada
#

"""The main module for simulating neutrino sources."""


import argparse

import numpy as np

# Firesong code
from firesong.Evolution import get_evolution, SourcePopulation
from firesong.Evolution import TransientSourcePopulation, cosmology
from firesong.Luminosity import get_LuminosityFunction
from firesong.input_output import output_writer, print_config, get_outputdir
from firesong.sampling import InverseCDF
from firesong.ConvertFlux import nu_to_gamma

def firesong_simulation(outputdir,
                        filename='Firesong.out',
                        density=1e-9,
                        Evolution="MD2014SFR",
                        Transient=False,
                        timescale=1000.,
                        zmax=10.,
                        bins=10000,
                        fluxnorm=1.44e-8,
                        index=2.28,
                        LF="SC",
                        luminosity=0.0,
                        Gammaflux=False,
                        interaction="pgamma",
                        emin=1e4,
                        emax=1e7,
                        seed=None,
                        verbose=True, 
                        **kwargs):
    """
    Simulate a universe of neutrino sources

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
        timescale (float, optional, default=1000): Timescale in seconds
            for transient sources
        zmax (float, optional, default=10.): Farthest redshift to consider
        bins (int, optional, default=1000): Number of bins used when creating
            the redshift PDF
        fluxnorm (float, optional, default=1.44e-8): Normalization on the total
            astrophysical diffuse flux, E^2dPhi/dE. Units of GeV s^-1 cm^-2 sr^-1
        index (float, optional, default=2.28): Spectral index of diffuse flux
        LF (string, optional, default="SC"): Luminosity function, choose 
            between standard candle (SC), LogNormal (LG), PowerLaw(PL),
            and Broken PowerLaw, (BPL)
        luminosity (float, optional, default=0.0): Manually fix the 
            luminosity of sources if not equal to 0. Overrides fluxnorm. 
            Units of erg/yr
        Gammaflux (bool, optional, default="False"): If true, calculate the
            gamma-ray flux based on the interaction type
        interaction (string, optional, default="pgamma"): Calculate the gamma-ray
            flux from the neutrino flux given an interaction type. Options are
            pgamma and pp
        emin (float, optional, default=1e4): Minimum neutrino energy in GeV
        emax (float, optional, default=1e7): Maximum neutrino energy in GeV
        seed (int or None, optional, default=None): random number seed
        verbose (bool, optional, default=True): print simulation paramaters
            if True else suppress printed output
        **kwargs: the required arguments to be passed to evolution model and 
             luminosity function. 
             The currentyl implemented ones are:
             Log-Normal: lg_width
             Power-Law: pl_lmin, pl_lmax, pl_alpha
             Broken Power-Law: bpl_lmin, bpl_lbreak, bpl_lmax, 
                               bpl_alpha1, bpl_alpha2
             Please refer to the documentation of 

    Returns:
        dict: keys contain simulation results, including the input params
            as well as the sources. Only returned if filename is None
    """

    if Transient:
        population = TransientSourcePopulation(cosmology,
                                               get_evolution(Evolution, **kwargs),
                                               timescale=timescale)
    else:
        population = SourcePopulation(cosmology, get_evolution(Evolution, **kwargs))

    N_sample = int(population.Nsources(density, zmax))

    if luminosity == 0.0:
        ## If luminosity not specified calculate luminosity from diffuse flux
        luminosity = population.StandardCandleLuminosity(fluxnorm,
                                                        density,
                                                        zmax,
                                                        index,
                                                        emin=emin,
                                                        emax=emax)

    luminosity_function = get_LuminosityFunction(luminosity, LF, **kwargs)

    delta_gamma = 2-index
    if verbose:
        print_config(LF, Transient, timescale, Evolution, density, N_sample,
                    luminosity, fluxnorm, delta_gamma, zmax, luminosity, Gammaflux,
                    interaction, mode=" - Calculating Neutrino CDFs ")

    ##################################################
    #        Simulation starts here
    ##################################################

    rng = np.random.RandomState(seed)

    redshift_bins = np.arange(0.0005, zmax, zmax/float(bins))
    RedshiftPDF = [population.RedshiftDistribution(redshift_bins[i])
                   for i in range(0, len(redshift_bins))]
    invCDF = InverseCDF(redshift_bins, RedshiftPDF)

    if filename is None:
        results = {}
        results['header'] = {'LF': LF, 'Transient': Transient, 
            'timescale': timescale, 'fluxnorm': fluxnorm, 
            'delta_gamma': delta_gamma,
            'luminosity': luminosity}
        sources = {}
    else:
        out = output_writer(outputdir, filename)
        out.write_header(LF, Transient, timescale, fluxnorm,
                        delta_gamma, luminosity, Gammaflux, interaction)

    # IMPORTANT notice, in the following "flux" means fluence in
    # Transient mode, but flux in steady source mode,
    # until TotalFlux(TotalFluence) is calculated
    
    # sample sources
    zs = invCDF(rng.uniform(low=0.0, high=1.0, size = N_sample))
    lumis = luminosity_function.sample_distribution(nsources=N_sample, rng=rng)
    if np.ndim(lumis) < 1:
        lumis = np.array([lumis]*N_sample)
    fluxes = population.Lumi2Flux(lumis, index, emin, emax, zs)
    # Random declination over the entire sky
    sinDecs = rng.uniform(-1, 1, size=N_sample)
    declins = np.degrees(np.arcsin(sinDecs))
    # Random ra over the sky
    ras = rng.uniform(0.,360., size=N_sample)

    TotalFlux = np.sum(fluxes)

    # For transient sources, the flux measured on Earth will be
    # red-shifted-fluence/{(1+z)*burst duration}
    if Transient:
        fluxes = population.fluence2flux(fluxes, zs)

    if Gammaflux:
        gammaFluxes = nu_to_gamma(fluxes, index, interaction)

    if filename is None:
        sources['dec'] = declins
        sources['ra'] = ras
        sources['flux'] = fluxes
        sources['z'] = zs
        if Gammaflux:
            sources['gammaFlux'] = gammaFluxes
    else:
        if Gammaflux:
            out.write_gamma(declins, ras, zs, fluxes, gammaFluxes)
        else:
            out.write(declins, ras, zs, fluxes)

    # For transient sources, we calculate the total fluence from all sources,
    # then obtain the diffuse flux by doing a time average over a year
    if Transient:
        TotalFlux = TotalFlux / population.yr2sec
    TotalFlux /= 4*np.pi  # give in per sr

    if filename is None:
        results['total_flux'] = TotalFlux
        results['sources'] = sources
    else:
        out.finish(TotalFlux)
    if verbose:
        print("Actual diffuse flux simulated :")
        log = "E^2 dNdE = {TotalFlux} (E/100 TeV)^({delta_gamma}) [GeV/cm^2.s.sr]"
        print(log.format(**locals()))
    if filename is None:
        return results

if __name__ == "__main__":
    outputdir = get_outputdir()

    # Process command line options
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', action='store', dest='filename',
                        default='Firesong.out', help='Output filename')
    parser.add_argument('-d', action='store', dest='density',
                        type=float, default=1e-9,
                        help='Local neutrino source density [1/Mpc^3] (default 1e-9 / Mpc^3)')
    parser.add_argument("--evolution", action="store",
                        dest="Evolution", default='MD2014SFR',
                        help="Source evolution options:  HB2006SFR, YMKBH2008SFR, CC2015SNR , MD2014SFR (default), NoEvolution")
    parser.add_argument("--transient", action='store_true',
                        dest='Transient', default=False,
                        help='Simulate transient sources (default False)')
    parser.add_argument("--timescale", action='store',
                        dest='timescale', type=float,
                        default=1000., help='time scale of transient sources in seconds (default 1000 seconds)')
    parser.add_argument("--zmax", action="store", type=float,
                        dest="zmax", default=10.,
                        help="Highest redshift to be simulated (default 10.)")
    parser.add_argument("--fluxnorm", action="store", dest='fluxnorm',
                        type=float, default=1.44e-8,
                        help="Astrophysical neutrino flux normalization A on E^2 dN/dE = A (E/100 TeV)^(-index) GeV/cm^2.s.sr (default 1.44e-8)")
    parser.add_argument("--index", action="store", dest='index',
                        type=float, default=2.28,
                        help="Astrophysical neutrino spectral index on E^2 dN/dE = A (E/100 TeV)^(-index) GeV/cm^2.s.sr (default 2.28)")
    parser.add_argument("--LF", action="store", dest="LF", default="SC",
                        help="Luminosity function, SC for standard candles (default), LG for lognormal, PL for powerlaw")
    parser.add_argument("--sigma", action="store",
                        dest="sigma", type=float, default=1.0,
                        help="Width of a log normal Luminosity function in dex (default 1.0)")
    parser.add_argument("--L", action="store",
                        dest="luminosity", type=float, default=0.0,
                        help="Set luminosity for each source [erg/year]. Resets fluxnorm option")
    parser.add_argument("--gammaflux", action='store_true',
                        dest="Gammaflux", default=False,
                        help="Calculate the gamma-ray flux")
    parser.add_argument("--interaction", action="store", dest='interaction',
                        type=str, default="pgamma",
                        help="Interaction type to calculate the gamma-ray flux. Options: pgamma, pp")
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
                        Gammaflux=options.Gammaflux,
                        interaction=options.interaction)
