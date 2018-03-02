#!/usr/bin/python

import numpy as np
import scipy
import argparse
from scipy.stats import lognorm
from Evolution import get_evolution, SourcePopulation, cosmology
from Luminosity import get_LuminosityFunction
from cosmolopy.distance import diff_comoving_volume

# Physics Settings
def calc_pdf(density=1e-7, L_nu=1e50, sigma=1, gamma=2.19,
             logMu_range=[-10, 6], N_Mu_bins=200,
             z_limits=[0.04, 10.], nzbins=120,
             Lum_limits=[1e45, 1e54], nLbins=120,
             flux_to_mu=10763342917.859608, upper_E=1e7, lower_E=1e4):
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
    
    pop = SourcePopulation(cosmology, get_evolution("HB2006SFR"))

    # Define the Redshift and Luminosity Evolution
    LF = get_LuminosityFunction(argparse.Namespace(LF="LG", sigma=sigma), L_nu)

    int_norm = pop.RedshiftIntegral(z_limits[-1])
    N_tot = pop.Nsources(density, z_limits[-1])
    

    # Setup Arrays
    logMu_array = np.linspace(logMu_range[0], logMu_range[1], N_Mu_bins)
    Flux_from_fixed_z = []

    zs = np.linspace(z_limits[0], z_limits[1], nzbins)
    deltaz = (float(z_limits[1])-float(z_limits[0]))/nzbins

    Ls = np.linspace(np.log10(Lum_limits[0]), np.log10(Lum_limits[1]), nLbins)
    deltaL = (np.log10(Lum_limits[1])-np.log10(Lum_limits[0]))/nLbins

    # Integration
    Count_array = np.zeros(N_Mu_bins)
    muError = []
    tot_bins = nLbins * nzbins
    Flux_from_fixed_z.append([])
    # Loop over redshift bins
    for z_count, z in enumerate(zs):
        
        # Conversion Factor for given z
        u_lim = upper_E/(1+z)   # upper IceCube range
        l_lim = lower_E/(1+z)   # lower IceCube range
        if gamma != 2.0:
            exponent = (2-gamma)
            nenner = 1/exponent*(u_lim**exponent-l_lim**exponent)
        else:
            nenner = (np.log(u_lim)-np.log(l_lim))
        bz = (1/(1e5)**(gamma-2))*1/(nenner)
        
        dlz = pop.LuminosityDistance(z)
        tot_flux_from_z = 0.

        # Loop over Luminosity bins
        for l_count, lum in enumerate(Ls):
                    
                # Number of Sources in
                LF_val = pop.evolution(z)*LF.pdf(10**lum)
                dN = (4*np.pi*LF_val*diff_comoving_volume(z, **pop.cosmology)/int_norm) * deltaL*deltaz*N_tot

                # Flux to Source Strength
                logmu = np.log10(flux_to_mu * 10**lum/pop.GeV_per_sec_2_ergs_per_year /
                                 (4*np.pi*(pop.Mpc2cm*dlz)**2)*bz)

                # Add dN to Histogram
                if logmu < logMu_range[1] and logmu > logMu_range[0]:
                    tot_flux_from_z += dN*10**logmu
                    idx = int((logmu-logMu_range[0])*N_Mu_bins /
                              (logMu_range[1]-logMu_range[0]))
                    Count_array[idx] += dN
                else:
                    muError.append(logmu)

        Flux_from_fixed_z.append(tot_flux_from_z)
        
    print "Total number of sources {:.0f} (All-Sky)".format(N_tot)
    return logMu_array, Count_array, zs, Flux_from_fixed_z


# logMu_array, Count_array, zs, Flux_from_fixed_z = calc_pdf(
#        density=args.density,
#        L_nu=args.luminosity,
#        sigma=args.sigma,
#        gamma=args.index,
#        flux_to_mu=1,
#        logMu_range=[-20, -6],
#        N_Mu_bins=2000,
#        z_limits=[0.04, 10.],  # lower bound in Firesong 0.0005
#        nzbins=1200,
#        Lum_limits=[args.luminosity*1e-5, args.luminosity*1e5],
#        nLbins=1200)
