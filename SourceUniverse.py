#!/usr/bin/python

import numpy as np
import argparse
from Evolution import get_evolution, SourcePopulation, cosmology, FluxPDF 
from Luminosity import get_LuminosityFunction


# Physics Settings
def main(luminosity=1e50, LF="LG", sigma=1, 
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
