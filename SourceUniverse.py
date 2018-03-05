#!/usr/bin/python

import numpy as np
import argparse
from Evolution import get_evolution, SourcePopulation, cosmology, FluxPDF 
from Luminosity import get_LuminosityFunction


# Physics Settings
def flux_pdf(luminosity=1e50, LF="LG", sigma=1, 
             index=2.19, emin=1e4, emax=1e7,
             density=1e-7, evol="HB2006SFR",
             zmin=0.005, zmax=10., nzbins=120,
             LumMin=1e45, LumMax=1e54, nLbins=120,
             logFMin=-10, logFMax=6, nFluxBins=200):
    """
    Parameter:
        - density in 1/Mpc^3
        - L_nu in erg/yr
        - sigma in dex
        - gamma
        - flux_to_mu

    Integration Parameters
        - logMu_range = [-10,6]      # expected range in log mu,
        - N_Mu_bins = 200            # number of bins for log nu histogram
        - z_limits = [0.04, 10.]     # Redshift limits
        - nzbins = 120               # number of z bins
        - Lum_limits = [1e45,1e54]   # Luminosity limits
        - nLbins = 120               # number of logLuminosity bins
    """

    population = SourcePopulation(cosmology, get_evolution(evol))
    luminosity_function = get_LuminosityFunction(luminosity, LF=LF, sigma=sigma)

    fPDF = FluxPDF(population, luminosity_function, index, emin, emax, density)

    return fPDF.calc_PDF(zmin=zmin, zmax=zmax, nzbins=nzbins,
                         LumMin=LumMin, LumMax=LumMax, nLbins=nLbins,
                         logFMin=logFMin, logFMax=logFMax, nFluxBins=nFluxBins)


if __name__ == "__main__":

    # Check that the Firesong environmental variable is set
    # This is needed for output and to read exposures, effective areas, etc
    try:
        firesongdir = os.environ['FIRESONG']
    except:
        print "Enviromental variable FIRESONG not set"
        quit()
    outputdir = firesongdir + "/Results/"

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
    
    output = flux_pdf(luminosity=1e50, LF="LG", sigma=1, 
                      index=2.19, emin=1e4, emax=1e7,
                      density=1e-7, evol="HB2006SFR",
                      zmin=0.005, zmax=10., nzbins=120,
                      LumMin=1e45, LumMax=1e54, nLbins=120,
                      logFMin=-10, logFMax=6, nFluxBins=200)
