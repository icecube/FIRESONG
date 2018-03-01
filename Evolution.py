#!/usr/bin/python
#
# x = log10(1+z)

import numpy as np
import scipy
import cosmolopy
cosmology = {'omega_M_0': 0.308, 'omega_lambda_0': 0.692, 'h': 0.678}


def get_evolution(evol):
    evolutions = {"NoEvolution": NoEvolution,
                  "HB2006SFR": HopkinsBeacom2006StarFormationRate,
                  "YMKBH2008SFR": YukselEtAl2008StarFormationRate,
                  "CC2015SNR": CandelsClash2015SNRate,
                  }
    if not evol in evolutions.keys():
        raise NotImplementedError("Source evolution " +
                                  evol + " not implemented.")

    return evolutions[evol]()


class Evolution(object):
    def __init__(self):
        pass

    def parametrization(self, x):
        raise NotImplementedError("Abstract")

    def __call__(self, z):
        return self.parametrization(np.log10(1.+z))


class NoEvolution(Evolution):
    def parametrization(self, x):
        return 1.


class HopkinsBeacom2006StarFormationRate(Evolution):
    """ StarFormationHistory (SFR), from Hopkins and Beacom 2006,
    unit = M_sun/yr/Mpc^3 """

    def parametrization(self, x):
        x = np.atleast_1d(x)
        result = np.zeros_like(x)
        m0 = x < 0.30963
        m1 = np.logical_and(x >= 0.30963, x < 0.73878)
        m2 = x >= 0.73878
        result[m0] = np.power(10, 3.28*x[m0]-1.82)
        result[m1] = np.power(10, -0.26*x[m1]-0.724)
        result[m2] = np.power(10, -8.0*x[m2]+4.99)
        if len(result) == 1:
            return np.asscalar(result)
        return result


class YukselEtAl2008StarFormationRate(Evolution):
    """ Star Formation Rate in units of M_sun/yr/Mpc^3
    arXiv:0804.4008  Eq.5
    """

    def __call__(self, z):
        return self.parametrization(1.+z)

    def parametrization(self, x):
        a = 3.4
        b = -0.3
        c = -3.5
        # z1 = 1
        # z2 =4
        # precomputed B = (1+z1)**(1-a/b)
        B = 5160.63662037
        # precomputed C = (1+z1)**((b-a)/c) * (1 + z2)**(1-b/c)
        C = 9.06337604231
        eta = -10
        r0 = 0.02
        return r0 * (x**(a*eta) + (x/B)**(b*eta) +
                     (x/C)**(c*eta))**(1./eta)


class CandelsClash2015SNRate(Evolution):
    def parametrization(self, x):
        a = 0.015
        b = 1.5
        c = 5.0
        d = 6.1
        density = a*(10.**x)**c / ((10.**x / b)**d+1.)
        return density


class SourcePopulation(object):
    def __init__(self, cosmology, evolution):
        self._zlocal = 0.01
        self.Mpc2cm = 3.086e24                    # Mpc / cm
        self.GeV_per_sec_2_ergs_per_year = 50526  # (GeV/sec) / (ergs/yr)
        self.evolution = evolution

        # Flat universe
        self.cosmology = cosmolopy.distance.set_omega_k_0(cosmology)
        self.dL1 = self.LuminosityDistance(1.)

    def RedshiftDistribution(self, z):
        """ can remove 4*pi becaue we just use this in a normalized way """
        return 4 * np.pi * self.evolution(z) * \
            cosmolopy.distance.diff_comoving_volume(z, **self.cosmology)

    def RedshiftIntegral(self, zmax):
        """ $$ \int_0^{z_\mathrm{max}} \frac{\mathrm{d}N}{\mathrm{d}z}
        \,\mathrm{d}V_c(z) \,\mathrm{d}z $$ """

        integrand = lambda z: self.RedshiftDistribution(z)
        return scipy.integrate.quad(integrand, 0, zmax)[0]

    def LuminosityDistance(self, z):
        # Wrapper function - so that cosmolopy is only imported here.
        return cosmolopy.distance.luminosity_distance(z, **self.cosmology)

    def Nsources(self, density, zmax):
        """ Total number of sources within $z_\mathrm{max}$:

        $$ N_\mathrm{tot} = \rho\cdot V_c(z=0.01)
        \frac{\int_0^{z_\mathrm{max}} \frac{\mathrm{d}N}{\mathrm{d}z}
        V_c(z) \,\mathrm{d}z}{\int_0^{0.01}
        \frac{\mathrm{d}N}{\mathrm{d}z} V_c(z) \,\mathrm{d}z} $$
        """
        vlocal = cosmolopy.distance.comoving_volume(self._zlocal,
                                                    **self.cosmology)
        Ntotal = density * vlocal / \
            (self.RedshiftIntegral(self._zlocal) /
             self.RedshiftIntegral(zmax))
        return Ntotal

    def Flux2Lumi(self, fluxnorm, index, emin, emax, E0=1e5):
        """
        $$ L_\nu = \frac{ \Phi_{z=1}^{PS} }{E_0^2}
        \int_{E_\mathrm{min}}^{E_\mathrm{max}} E
        \left(\frac{E}{E_0}\right)^{-\gamma}\,
        \mathrm{d}E\,4\pi d_L^2(z=1) $$

        Note fluxnorm is E0^2*fluxnorm
        fluxnorm units are []
        """
        integrand = lambda E: E*(E/E0)**(-abs(index))
        flux_integral = scipy.integrate.quad(integrand, emin, emax)[0]
        luminosity = fluxnorm / E0**2. * flux_integral *  \
            self.GeV_per_sec_2_ergs_per_year * \
            4. * np.pi * (self.dL1*self.Mpc2cm)**2.
        # CHECK why (1+z)/2
        return luminosity * 2.**(index-2)

    def Lumi2Flux(self, luminosity, index, emin, emax, E0=1.e5):
        """
        $$ L_\nu = \frac{ \Phi_{z=1}^{PS} }{E_0^2}
        \int_{E_\mathrm{min}}^{E_\mathrm{max}} E
        \left(\frac{E}{E_0}\right)^{-\gamma}\,
        \mathrm{d}E\,4\pi d_L^2(z=1) $$

        Lumi given in ergs/yr
        Note fluxnorm is E0^2*fluxnorm
        fluxnorm units are []
        """
        integrand = lambda E: E*(E/E0)**(-abs(index))
        flux_integral = scipy.integrate.quad(integrand, emin, emax)[0]
        fluxnorm = luminosity / 4. / np.pi / \
            (self.dL1*self.Mpc2cm)**2. / \
            self.GeV_per_sec_2_ergs_per_year / flux_integral * E0**2.
        # CHECK why (1+z)/2
        return fluxnorm * 2.**(2.-index)

    def StandardCandleSources(self, fluxnorm, density, zmax, index):
        """ $$ \Phi_{z=1}^{PS} = \frac{4 \pi \Phi_\mathrm{diffuse}}
        {N_\mathrm{tot}\,d_L^2(z=1)\, \int_0^{10}
        \frac{ (1+z)^{-\gamma+2} }{d_L(z)^2}
        \frac{\frac{\mathrm{d}N}{\mathrm{d}z} V_c(z)}
        { \int_0^{z_\mathrm{max}} \frac{\mathrm{d}N}{\mathrm{d}z'}
        V_c(z') \,\mathrm{d}z'} \,\mathrm{d}z} $$
        """
        norm = self.RedshiftIntegral(zmax)
        Ntotal = self.Nsources(density, zmax)
        all_sky_flux = 4 * np.pi * fluxnorm

        # Here the integral on redshift is done from 0 to 10.
        # This insures proper normalization even if zmax is not 10.
        # CHECK why (1+z)/2.
        Fluxnorm = all_sky_flux / Ntotal / self.dL1**2. / \
            scipy.integrate.quad(lambda z: ((1.+z)/2.)**(-abs(index)+2) /
                                 self.LuminosityDistance(z)**2. *
                                 self.RedshiftDistribution(z) / norm,
                                 0, 10.)[0]

        return Fluxnorm

    def fluxFromRelative(self, flux_z1, z, index):
        dL = self.LuminosityDistance(z)
        # CHECK why (1+z)/2.
        flux = flux_z1 * (self.dL1**2)/(dL**2) * ((1.+z)/2.)**(-index+2)
        return flux


class TransientSourcePopulation(SourcePopulation):
    def __init__(self, cosmology, evolution, timescale):
        super(TransientSourcePopulation, self).__init__(cosmology, evolution)
        self.timescale = timescale
        self.yr2sec = 86400*365

    def RedshiftDistribution(self, z):
        return super(TransientSourcePopulation, self).RedshiftDistribution(z) / (1.+z)

    def StandardCandleSources(self, fluxnorm, density, zmax, index):
        # For transient source, Fluxnorm will be the fluence of a
        # standard candle at z=1, with unit GeV/cm^2 given that the
        # burst rate density is measured in per year.

        norm = self.RedshiftIntegral(zmax)
        Ntotal = self.Nsources(density, zmax)
        all_sky_flux = 4 * np.pi * fluxnorm * self.yr2sec

        # As above, the integral is done from redshift 0 to 10.
        # CHECK why (1+z)/2.
        fluence = all_sky_flux / Ntotal / self.dL1**2. / \
            scipy.integrate.quad(lambda z: ((1.+z)/2.)**(-abs(index)+3) /
                                 (self.LuminosityDistance(z)**2.) *
                                 self.RedshiftDistribution(z) / norm,
                                 0, 10.)[0]

        return fluence

    def Flux2Lumi(self, fluxnorm, index, emin, emax, E0=1e5):
        luminosity = super(TransientSourcePopulation, self).Flux2Lumi(fluxnorm,
                                                                      index,
                                                                      emin,
                                                                      emax,
                                                                      E0=E0)
        return luminosity / self.timescale

    def Lumi2Flux(self, luminosity, index, emin, emax, E0=1e5):
        flux = super(TransientSourcePopulation, self).Lumi2Flux(luminosity,
                                                                index,
                                                                emin,
                                                                emax,
                                                                E0=E0)
        return flux * self.timescale

    def fluxFromRelative(self, flux_z1, z, index):
        flux = super(TransientSourcePopulation, self).fluxFromRelative(flux_z1,
                                                                       z,
                                                                       index)
        # CHECK why (1+z)/2.
        return flux * (1.+z) / 2.

    def fluence2flux(self, fluence, z):
        # For transient sources, the flux measured on Earth will be
        # red-shifted-fluence/{(1+z)*burst duration}
        flux = fluence / ((1.+z)*self.timescale)
        return flux


class Simulation(object):
    def __init__(self, population, luminosity_function, index, zmax,
                 seed=None, zmin=0.0005, bins=10000):
        self.population = population
        self.luminosity_function = luminosity_function
        self.index = abs(index)
        self.zmax = zmax
        self.zmin = zmin
        self.bins = bins
        self.rng = np.random.RandomState(seed)
        self.setup()

    def setup(self):
        redshift_bins = np.arange(self.zmin, self.zmax,
                                  self.zmax/float(self.bins))

        # RedshiftCDF is used for inverse transform sampling
        RedshiftPDF = [self.population.RedshiftDistribution(redshift_bins[i])
                       for i in range(0, len(redshift_bins))]
        RedshiftCDF = np.cumsum(RedshiftPDF)
        RedshiftCDF = RedshiftCDF / RedshiftCDF[-1]

        self.redshift_bins = redshift_bins
        self.redshift_cdf = RedshiftCDF

    def sample_redshift(self, N=None):
        # Generate a histogram to store redshifts.
        # Starts at z = 0.0005 and increases in steps of 0.001
        rand_cdf = self.rng.uniform(0, 1, N)
        bin_index = np.searchsorted(self.redshift_cdf, rand_cdf)
        z = self.redshift_bins[bin_index]
        return z

    def sample_flux(self, N=None):
        z = self.sample_redshift(N)

        flux_z1 = self.luminosity_function.sample_distribution(N, rng=self.rng)
        flux = self.population.fluxFromRelative(flux_z1, z, self.index)

        # Random declination over the entire sky
        sinDec = self.rng.uniform(-1, 1, N)
        declin = np.degrees(np.arcsin(sinDec))

        return flux, z, declin
