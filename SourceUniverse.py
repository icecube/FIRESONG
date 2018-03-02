#!/usr/bin/python

import numpy as np
import scipy
from scipy.stats import lognorm
try:
    from cosmolopy.distance import diff_comoving_volume, luminosity_distance
    from cosmolopy.distance import comoving_volume, set_omega_k_0
except:
    pass

# StarFormationHistory from Hopkins & Beacom 2006, unit = M_sun/yr/Mpc^3
@np.vectorize
def HopkinsBeacom2006StarFormationRate(z):
    x = np.log10(1.+z)
    if x < 0.30963:
        return np.power(10, 3.28*x-1.82)
    if (x >= 0.30963)and(x < 0.73878):
        return np.power(10, -0.26*x-0.724)
    if x >= 0.73878:
        return np.power(10, -8.0*x+4.99)

# Use local source Density to deduce total number of sources in the universe
def tot_num_src(redshift_evolution, cosmology, zmax, density):
    integrand = lambda z: redshift_evolution(z) * \
        diff_comoving_volume(z, **cosmology)
    norm = 4 * np.pi * scipy.integrate.quad(integrand, 0, zmax)[0]
    area = 4 * np.pi * scipy.integrate.quad(integrand, 0, 0.01)[0]
    vlocal = comoving_volume(0.01, **cosmology)
    Ntotal = density * vlocal / (area/norm)
    return Ntotal

# Physics Settings
def calc_pdf(density=1e-7, L_nu=1e50, sigma=1, gamma=2.19,
             logMu_range=[-10, 6], N_Mu_bins=200,
             z_limits=[0.04, 10.], nzbins=120,
             Lum_limits=[1e45, 1e54], nLbins=120,
             flux_to_mu=10763342917.859608):
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
    # Conversion Factors
    Mpc_to_cm = 3.086e+24
    erg_to_GeV = 624.151
    year2sec = 365*24*3600

    cosmology = {'omega_M_0': 0.308, 'omega_lambda_0': 0.692, 'h': 0.678}
    cosmology = set_omega_k_0(cosmology)  # Flat universe

    # Define the Redshift and Luminosity Evolution
    redshift_evolution = lambda z: HopkinsBeacom2006StarFormationRate(z)
    LF = lambda logL: np.log(10)*10**logL * \
        lognorm.pdf(10**logL, np.log(10)*sigma,
                    scale=L_nu*np.exp(-0.5*(np.log(10)*sigma)**2))

    N_tot = tot_num_src(redshift_evolution, cosmology, z_limits[-1], density)
    integrand = lambda z: redshift_evolution(z) * diff_comoving_volume(z, **cosmology)
    int_norm = 4 * np.pi * scipy.integrate.quad(integrand, 0, z_limits[-1])[0]
    print "Total number of sources {:.0f} (All-Sky)".format(N_tot)

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
    N_sum = 0
    Flux_from_fixed_z.append([])
    # Loop over redshift bins
    for z_count, z in enumerate(zs):
        # Conversion Factor for given z
        upper_E=1e7
        lower_E=1e4
        u_lim = upper_E/(1+z)   # upper IceCube range
        l_lim = lower_E/(1+z)   # lower IceCube range

        if gamma != 2.0:
            exponent = (2-gamma)
            nenner = 1/exponent*(u_lim**exponent-l_lim**exponent)

        else:
            nenner = (np.log(u_lim)-np.log(l_lim))

        bz = (1/(1e5)**(gamma-2))*1/(nenner)
        dlz = luminosity_distance(z, **cosmology)
        tot_flux_from_z = 0.

        # Loop over Luminosity bins
        for l_count, lum in enumerate(Ls):
                run_id = z_count*nLbins+l_count
                if run_id % (tot_bins/10) == 0.:
                    print "{}%".format(100*run_id/tot_bins)
                # Number of Sources in
                LF_val = redshift_evolution(z)*LF(lum)
                dN = (4*np.pi*LF_val*diff_comoving_volume(z, **cosmology)/int_norm) * deltaL*deltaz*N_tot

                N_sum += dN

                # Flux to Source Strength
                logmu = np.log10(flux_to_mu * erg_to_GeV*10**lum/year2sec /
                                 (4*np.pi*(Mpc_to_cm*dlz)**2)*bz)

                # Add dN to Histogram
                if logmu < logMu_range[1] and logmu > logMu_range[0]:
                    tot_flux_from_z += dN*10**logmu
                    idx = int((logmu-logMu_range[0])*N_Mu_bins /
                              (logMu_range[1]-logMu_range[0]))
                    Count_array[idx] += dN
                else:
                    muError.append(logmu)

        Flux_from_fixed_z.append(tot_flux_from_z)

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
