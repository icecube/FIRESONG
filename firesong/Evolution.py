#!/usr/bin/python

"""
Calculates features of various cosmic evolution models

Cosmological parameters are hardcoded to Planck (2018) results: 
\(\Omega_{M} = 0.315\), \(\Omega_{L} = 1 - \Omega_{M}\), \(h = 0.674 \)

Planck Collaboration A&A 641, A6 (2020)
arXiv:1807.06209

"""

import numpy as np
import scipy
import cosmolopy
# These are Planck 2015 values
#cosmology = {'omega_M_0': 0.308, 'omega_lambda_0': 0.692, 'h': 0.678}
cosmology = {'omega_M_0': 0.315, 'omega_lambda_0': 0.685, 'h': 0.674}

def get_evolution(evol):
    """
    Get specific evolution model

    Args:
        evol (str): Name of evolution model, options are "NoEvolution",
            "HB2006SFR", "YMKBH2008SFR", "CC2015SNR", "MD2014SFR". See specific
            classes for more details of each model

    Returns:
        Evolution: relevant Evolution object
    """
    evolutions = {"NoEvolution": NoEvolution,
                  "HB2006SFR": HopkinsBeacom2006StarFormationRate,
                  "YMKBH2008SFR": YukselEtAl2008StarFormationRate,
                  "CC2015SNR": CandelsClash2015SNRate,
                  "MD2014SFR": MadauDickinson2014CSFH
                  }
    if not evol in list(evolutions.keys()):
        raise NotImplementedError("Source evolution " +
                                  evol + " not implemented.")

    return evolutions[evol]()

class Evolution(object):
    """
    Abstract class to handle all evolution models
    """
    def __init__(self):
        pass

    def parametrization(self, x):
        raise NotImplementedError("Abstract")

    def __call__(self, z):
        return self.parametrization(np.log10(1.+z))


class NoEvolution(Evolution):
    """
    Evolution model that is flat over cosmic history
    """
    def parametrization(self, x):
        return 1.

    def __str__(self):
        return "No Evolution"


class HopkinsBeacom2006StarFormationRate(Evolution):
    """ 
    StarFormationHistory (SFR), from Hopkins and Beacom 2006,
    unit = M_sun/yr/Mpc^3 

    Model is a piecewise linear fit with the following segments in
    log10(1+z) - log10(rho) space:

    intercepts, slopes, domain: 
    
    -1.82, 3.28, z <= 1.04 

    -0.724, -0.26, 1.04 <= z <= 4.48 

    4.99, -8.0, 4.48 <= z

    Reference: doi:10.1086/506610
    """

    def parametrization(self, x):
        """
        Star formation rate at a given redshift

        Args:
            x (array or float): 1 + z values

        Returns:
            Array or float: Star formation rate
        """
        x = np.atleast_1d(x)
        result = np.zeros_like(x)
        m0 = x < 0.30963
        m1 = np.logical_and(x >= 0.30963, x < 0.73878)
        m2 = x >= 0.73878
        result[m0] = np.power(10, 3.28*x[m0]-1.82)
        result[m1] = np.power(10, -0.26*x[m1]-0.724)
        result[m2] = np.power(10, -8.0*x[m2]+4.99)
        if len(result) == 1:
            return result.item()
        return result

    def __str__(self):
        return "Hopkins and Beacom (2006)"

class YukselEtAl2008StarFormationRate(Evolution):
    r""" 
    Star Formation Rate in units of \(\frac{M_{sun}}{yr Mpc^3}\)

    Model is a continuous broken power law,
    $$ \dot{\rho}_{*}(z)=\dot{\rho}_{0}\left[(1+z)^{a \eta}
    +\left(\frac{1+z}{B}\right)^{b \eta}
    +\left(\frac{1+z}{C}\right)^{c \eta}\right]^{1 / \eta} $$

    with a = 3.4, b=-0.3, c=-3.5, B=5160.63662037,
    C=9.06337604231, \(\dot{\rho}\)=0.02, eta=10

    The given function results in breaks around z=1,4

    Reference: arXiv:0804.4008  Eq.5
    """

    def __call__(self, z):
        return self.parametrization(1.+z)

    def parametrization(self, x):
        """
        Star formation rate at a given redshift

        Args:
            x (array or float): 1 + z values

        Returns:
            Array or float: Star formation rate
        """
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

    def __str__(self):
            return "Yuksel et al. (2008)"

class CandelsClash2015SNRate(Evolution):
    r"""
    This is the implied SFR from Goods/Candels/Clash (2015)
    derive from CC SNe rate and assuming one rate is proportional to the other.
    They use the same functional form as Madau and Dickinson (2014)
    unit = M_sun/yr/Mpc^3 

    Model takes the functional form of 
    $$ \psi(z)=\frac{A(1+z)^{C}}{((1+z) / B)^{D}+1} $$
    with best-fit values A = 0.015, B = 1.5, C = 5.0, D = 6.1

    Reference: arXiv:1509.06574
    """
    
    def __call__(self, z):
        return self.parametrization(1.+z)
    
    def parametrization(self, x):
        """
        Star formation rate at a given redshift

        Args:
            x (array or float): 1 + z values

        Returns:
            Array or float: Star formation rate
        """
        a = 0.015
        b = 1.5
        c = 5.0
        d = 6.1
        density = a*(x**c) / (1. + ( x / b)**d )
        return density
    
    def __str__(self):
        return "Strolger et al. (2015)"


class MadauDickinson2014CSFH(Evolution):
    r""" 
    StarFormationHistory (SFR), from Madau and Dickinson (2014),
    unit = M_sun/yr/Mpc^3 

    Model takes the same functional form as Candels/Clash,  
    $$ \psi(z)=\frac{A(1+z)^{C}}{((1+z) / B)^{D}+1} $$
    but with best-fit parameters A = 0.015, B = 2.7, C = 2.9, D = 5.6

    Reference: arXiv:1403.0007
    """

    def __call__(self, z):
        return self.parametrization(1.+z)
        
    def parametrization(self, x):
        """
        Star formation rate at a given redshift

        Args:
            x (array or float): 1 + z values

        Returns:
            Array or float: Star formation rate
        """
        a = 0.015
        b = 2.7
        c = 2.9
        d = 5.6
        density = a*(x**b) / (1. + (x/c)**d ) 
        return density

    def __str__(self):
        return "Madau and Dickinson (2014)"


class SourcePopulation(object):
    """
    Given an evolution to follow, create a population
    of neutrino sources

    Args:
        cosmology (dict): kwargs to pass to cosmolopy, defaults are
            'omega_M_0': 0.308, 'omega_lambda_0': 0.692, 'h': 0.678
        evolution (Evolution instance): Evolution model for neutrino
            source population

    Attributes:
        _zlocal (float): Describes limit of nearby sources
        Mpc2cm (float): Conversion factor
        GeV_per_sec_2_ergs_per_year (float): Conversion factor
        evolution (Evolution): Evolution model for neutrino source population
        cosmology (cosmolopy instance)
    """

    def __init__(self, cosmology, evolution):
        """
        """
        self._zlocal = 0.01
        self.Mpc2cm = 3.086e24                     # Mpc / cm
        self.GeV_per_sec_2_ergs_per_year = 50526.  # (GeV/sec) / (ergs/yr)
        self.evolution = evolution

        # Flat universe
        self.cosmology = cosmolopy.distance.set_omega_k_0(cosmology)

    def RedshiftDistribution(self, z):
        r""" 
        Provides the unnormalized PDF of number of sources vs. redshift
        by multiplying the \(\frac{dN}{dz} = \frac{d\rho}{dz} \times \frac{dV}{dz}\)
        Note: can remove 4*pi becaue we just use this in a normalized way 

        Args:
            z (array or float): Redshift values

        Returns
            Array of float: Unnormalized PDF of number vs. redshift
        """
        return 4 * np.pi * self.evolution(z) * \
            cosmolopy.distance.diff_comoving_volume(z, **self.cosmology)

    def RedshiftIntegral(self, zmax):
        r""" 
        Integrates the redshift distribution to find the total
        number of sources (before accounting for density) out to zmax

        $$ \int_0^{z_\mathrm{max}} \frac{\mathrm{d}N}{\mathrm{d}z}
        \,\mathrm{d}V_c(z) \,\mathrm{d}z $$ 

        Args:
            zmax (float): upper bound of integral

        Returns:
            float: Number of sources from z=0 to z=zmax
        """

        integrand = lambda z: self.RedshiftDistribution(z)
        return scipy.integrate.quad(integrand, 0, zmax)[0]

    def LuminosityDistance(self, z):
        """
        Convert redshift to luminosity distance.
        If passing many redshifts, a 1d spline is used as cosmolopy can be slow

        Args:
            z (array or float): redshift(s)

        Returns:
            array or float: Luminosity distance(s) in Mpc
        """
        # Wrapper function - so that cosmolopy is only imported here.
        if np.ndim(z) > 0:
            if len(z) > 1000:
                zz = np.linspace(0., 10., 500)
                spl = scipy.interpolate.UnivariateSpline(zz, 
                        cosmolopy.distance.luminosity_distance(zz, 
                            **self.cosmology))
                return spl(z)
        return cosmolopy.distance.luminosity_distance(z, **self.cosmology)

    def Nsources(self, density, zmax):
        r""" Total number of sources within \(z_{\mathrm{max}}\):

        $$ N_\mathrm{tot} = \rho\cdot V_c(z=0.01)
        \frac{\int_0^{z_\mathrm{max}} \frac{\mathrm{d}N}{\mathrm{d}z}
        V_c(z) \,\mathrm{d}z}{\int_0^{0.01}
        \frac{\mathrm{d}N}{\mathrm{d}z} V_c(z) \,\mathrm{d}z} $$

        Args:
            density (float): local density of neutrino sources in Mpc^-3
            zmax (float): Maximal redshift to consider

        Returns:
            float: total number of sources within z_max
        """
        vlocal = cosmolopy.distance.comoving_volume(self._zlocal,
                                                    **self.cosmology)
        Ntotal = density * vlocal / \
            (self.RedshiftIntegral(self._zlocal) /
             self.RedshiftIntegral(zmax))
        return Ntotal

    def Flux2Lumi(self, fluxnorm, index, emin, emax, z=1, E0=1e5):
        r"""
        Converts a flux to a luminosity

        $$ L_\nu = \frac{ \Phi_{z=1}^{PS} }{E_0^2}
        \int_{E_\mathrm{min}}^{E_\mathrm{max}} E
        \left(\frac{E}{E_0}\right)^{-\gamma}\,
        \mathrm{d}E\ \times 4\pi d_L^2(z=1) $$

        Note fluxnorm is E0^2*fluxnorm
        fluxnorm units are [UNITS]

        Args:
            fluxnorm (array or float): Flux of a source in UNITS
            index (float): Spectral index of the flux
            emin (float): Minimum neutrino energy in GeV
            emax (float): Maximum neutrino energy in GeV
            z (array or float, optional, default=1): Redshifts
            E0 (float, optional, default=1): pivot energy in GeV

        Returns:
            float: luminosity in ergs/yr
        """
        flux_integral = self.EnergyIntegral(index, emin, emax, z, E0)
        luminosity = fluxnorm / E0**2. * flux_integral *  \
            self.GeV_per_sec_2_ergs_per_year * \
            4. * np.pi * (self.LuminosityDistance(z)*self.Mpc2cm)**2.
        return luminosity

    def Lumi2Flux(self, luminosity, index, emin, emax, z=1, E0=1e5):
        r"""
        Converts a luminosity to a flux

        $$ L_\nu = \frac{ \Phi_{z=1}^{PS} }{E_0^2}
        \int_{E_\mathrm{min}}^{E_\mathrm{max}} E
        \left(\frac{E}{E_0}\right)^{-\gamma}\,
        \mathrm{d}E\ \times 4\pi d_L^2(z=1) $$

        Note fluxnorm is E0^2*fluxnorm
        fluxnorm units are [UNITS]

        Args:
            luminosity (array or float): luminosity of sources in ergs/yr
            index (float): Spectral index of the flux
            emin (float): Minimum neutrino energy in GeV
            emax (float): Maximum neutrino energy in GeV
            z (array or float, optional, default=1): Redshifts
            E0 (float, optional, default=1): pivot energy in GeV
        
        Returns:
            fluxnorm (array or float): flux of source(s) in UNITS
        """
        flux_integral = self.EnergyIntegral(index, emin, emax, z, E0)
        fluxnorm = luminosity / 4. / np.pi / \
            (self.LuminosityDistance(z)*self.Mpc2cm)**2. / \
            self.GeV_per_sec_2_ergs_per_year / flux_integral * E0**2.
        return fluxnorm

    def EnergyIntegral(self, index, emin, emax, z=1, E0=1e5):
        r""" 
        Calculates energy content in a neutrino flux

        $$\int_{emin/(1+z)}^{emax/(1+z)} E*(E/E0)^{-index} dE$$ 
        
        Args:
            index (float): Spectral index of the flux
            emin (float): Minimum neutrino energy in GeV
            emax (float): Maximum neutrino energy in GeV
            z (array or float, optional, default=1): Redshifts
            E0 (float, optional, default=1): pivot energy in GeV

        Returns:
            float: Energy flux between emin and emax
        """
        if index != 2.0:
            denom = (1+z)**(index-2)
            integral = denom * (emax**(2-index)-emin**(2-index)) / (2-index)
        else:
            integral = np.ones_like(z) * np.log(emax/emin)
        return E0**index * integral

    def StandardCandleSources(self, fluxnorm, density, zmax, index, z0=1.):
        r""" 
        Given a total diffuse neutrino flux, calculate the individual 
        flux contribution from a single source
        
        $$ \Phi_{z=1}^{PS} = \frac{4 \pi \Phi_\mathrm{diffuse}}
        {N_\mathrm{tot}\,d_L^2(z=1)\, \int_0^{10}
        \frac{ (1+z)^{-\gamma+2} }{d_L(z)^2}
        \frac{\frac{\mathrm{d}N}{\mathrm{d}z} V_c(z)}
        { \int_0^{z_\mathrm{max}} \frac{\mathrm{d}N}{\mathrm{d}z'}
        V_c(z') \,\mathrm{d}z'} \,\mathrm{d}z} $$

        Args:
            fluxnorm (float): diffuse astrophysical neutrino flux in UNITS
            density (float): local density of neutrino sources in Mpc^-3
            zmax (float): Maximum redshift considered
            index (float): Spectral index of the flux
            z0 (float, optional, default=1.): Redshift of the source in 
                question

        Returns:
            float: fluxnorm of a source at redshift z0
        """
        norm = self.RedshiftIntegral(zmax)
        Ntotal = self.Nsources(density, zmax)
        all_sky_flux = 4 * np.pi * fluxnorm

        # Here the integral on redshift is done from 0 to 10.
        # This insures proper normalization even if zmax is not 10.
        Fluxnorm = all_sky_flux / Ntotal / self.LuminosityDistance(z0)**2. / \
            scipy.integrate.quad(lambda z: ((1.+z)/(1.+z0))**(2-index) /
                                 self.LuminosityDistance(z)**2. *
                                 self.RedshiftDistribution(z) / norm,
                                 0, 10.)[0]

        return Fluxnorm

    def StandardCandleLuminosity(self, fluxnorm, density, zmax, index,
                                 emin, emax, E0=1e5):
        """ 
        Calculates the standard candle luminosity that characterizes a 
        population of sources which have a fixed total flux

        Args:
            fluxnorm (float): diffuse astrophysical neutrino flux in UNITS
            density (float): local density of neutrino sources in Mpc^-3
            zmax (float): Maximum redshift considered
            index (float): Spectral index of the flux
            emin (float): Minimum neutrino energy in GeV
            emax (float): Maximum neutrino energy in GeV
            E0 (float, optional, default=1): pivot energy in GeV

        Returns:
            float: characteristic luminosity of the population
        """

        flux = self.StandardCandleSources(fluxnorm, density, zmax, index, z0=1)
        luminosity = self.Flux2Lumi(flux, index, emin, emax, z=1, E0=E0)
        return luminosity


class TransientSourcePopulation(SourcePopulation):
    """
    Given an evolution to follow, create a population
    of neutrino sources that only emit for a finite period of time

    See also: :class:`SourcePopulation`

    Args:
        cosmology (dict): kwargs to pass to cosmolopy, defaults are
            'omega_M_0': 0.308, 'omega_lambda_0': 0.692, 'h': 0.678
        evolution (Evolution instance): Evolution model for neutrino
            source population
        timescale (float): Duration (in seconds) of emission

    Attributes:
        timescale (float): Duration (in seconds) of emission
        yr2sec (float): Conversion factor
    """

    def __init__(self, cosmology, evolution, timescale):
        """
        """
        super(TransientSourcePopulation, self).__init__(cosmology, evolution)
        self.timescale = timescale
        self.yr2sec = 86400*365

    def RedshiftDistribution(self, z):
        r"""
        Provides the unnormalized PDF of number of sources vs. redshift
        by multiplying the \(\frac{dN}{dz} = \frac{d\rho}{dz} \times \frac{dV}{dz}\). Corrects for 
        time-dilation with extra factor of 1/1+z

        Args:
            z (array or float): Redshift values

        Returns
            Array of float: Unnormalized PDF of number vs. redshift
        """
        return super(TransientSourcePopulation, self).RedshiftDistribution(z) / (1.+z)

    def StandardCandleSources(self, fluxnorm, density, zmax, index, z0=1.):
        r""" 
        Given a total diffuse neutrino flux, calculate the individual 
        fluence contribution from a single standard candle source,
        given that the burst rate density is measured in per year
        
        $$ \Phi_{z=1}^{PS} = \frac{4 \pi \Phi_\mathrm{diffuse}}
        {N_\mathrm{tot}\,d_L^2(z=1)\, \int_0^{10}
        \frac{ (1+z)^{-\gamma+2} }{d_L(z)^2}
        \frac{\frac{\mathrm{d}N}{\mathrm{d}z} V_c(z)}
        { \int_0^{z_\mathrm{max}} \frac{\mathrm{d}N}{\mathrm{d}z'}
        V_c(z') \,\mathrm{d}z'} \,\mathrm{d}z} $$

        Args:
            fluxnorm (float): diffuse astrophysical neutrino flux in UNITS
            density (float): local density of neutrino sources in Mpc^-3
            zmax (float): Maximum redshift considered
            index (float): Spectral index of the flux
            z0 (float, optional, default=1.): Redshift of the source in 
                question

        Returns:
            float: fluence of a source at redshift z0 in GeV/cm^2
        """
        norm = self.RedshiftIntegral(zmax)
        Ntotal = self.Nsources(density, zmax)
        all_sky_flux = 4 * np.pi * fluxnorm * self.yr2sec

        # As above, the integral is done from redshift 0 to 10.
        fluence = all_sky_flux / Ntotal / self.LuminosityDistance(z0)**2. / \
            scipy.integrate.quad(lambda z: ((1.+z)/(1.+z0))**(3-index) /
                                 (self.LuminosityDistance(z)**2.) *
                                 self.RedshiftDistribution(z) / norm,
                                 0, 10.)[0]

        return fluence

    def Flux2Lumi(self, fluxnorm, index, emin, emax, z=1, E0=1e5):
        r"""
        Converts a fluence to a luminosity. Transient sources require
        fluence to be divided by timescale so that luminosity has
        proper units

        $$ L_\nu = \frac{ \Phi_{z=1}^{PS} }{E_0^2}
        \int_{E_\mathrm{min}}^{E_\mathrm{max}} E
        \left(\frac{E}{E_0}\right)^{-\gamma}\,
        \mathrm{d}E\,4\pi d_L^2(z=1) $$

        Note fluxnorm is E0^2*fluxnorm
        fluxnorm units are [UNITS]

        Args:
            fluxnorm (array or float): Flux of a source in UNITS
            index (float): Spectral index of the flux
            emin (float): Minimum neutrino energy in GeV
            emax (float): Maximum neutrino energy in GeV
            z (array or float, optional, default=1): Redshifts
            E0 (float, optional, default=1): pivot energy in GeV

        Returns:
            float: luminosity in UNITS
        """
        luminosity = super(TransientSourcePopulation, self).Flux2Lumi(fluxnorm,
                                                                      index,
                                                                      emin,
                                                                      emax,
                                                                      z=z,
                                                                      E0=E0)
        return luminosity / self.timescale

    def Lumi2Flux(self, luminosity, index, emin, emax, z=1, E0=1e5):
        r"""
        Converts a luminosity to a fluence

        $$ L_\nu = \frac{ \Phi_{z=1}^{PS} }{E_0^2}
        \int_{E_\mathrm{min}}^{E_\mathrm{max}} E
        \left(\frac{E}{E_0}\right)^{-\gamma}\,
        \mathrm{d}E\ \times 4\pi d_L^2(z=1) $$

        Note fluxnorm is E0^2*fluxnorm
        fluence units are [UNITS]

        Args:
            luminosity (array or float): luminosity of sources in ergs/yr
            index (float): Spectral index of the flux
            emin (float): Minimum neutrino energy in GeV
            emax (float): Maximum neutrino energy in GeV
            z (array or float, optional, default=1): Redshifts
            E0 (float, optional, default=1): pivot energy in GeV
        
        Returns:
            array or float: fluence of source(s) in UNITS
        """
        flux = super(TransientSourcePopulation, self).Lumi2Flux(luminosity,
                                                                index,
                                                                emin,
                                                                emax,
                                                                z=z,
                                                                E0=E0)
        return flux * self.timescale

    def fluence2flux(self, fluence, z):
        """
        Calculates flux measured on Earth, which is red-shifted fluence
        divided by (1+z)*timescale

        Args:
            fluence (array or float): fluence of source(s) in UNITS
            z (array or float): redshift of source(s)

        Returns:
            array or float: fluxes of the sources in UNITS
        """
        # For transient sources, the flux measured on Earth will be
        # red-shifted-fluence/{(1+z)*burst duration}
        flux = fluence / ((1.+z)*self.timescale)
        return flux

#############
#LEGEND AREA#
#############

def get_LEvolution(le_model, lmin, lmax):
    """
    Get specific LuminosityEvolution model (a luminosity distribution
    that is a function of z)

    Args:
        le_model (str): Name of luminosity-evolution model, only supported
            optioin is "HA2014BL"
        lmin (float): log10 of Minimum luminosity considered in erg/s
        lmax (float): log10 of Maximum luminosity considered in erg/s

    Returns:
        LuminosityEvolution: relevant luminosity-evolution object
    """
    evolutions = {"HA2014BL": HardingAbazajian(lmin, lmax)
                  }
    if not le_model in list(evolutions.keys()):
        raise NotImplementedError("Luminosity Evolution " +
                                  le_model + " not implemented.")
    return evolutions[le_model]

class LuminosityEvolution(object):
    """
    Abstract class for the a Luminosity Distribution that depends on z

    Args:
        lmin (float): log10 of Minimum luminosity considered in erg/s
        lmax (float): log10 of Maximum luminosity considered in erg/s
        cosmology (dict, optional, default=cosmology): kwargs to pass 
            to cosmolopy, defaults are 'omega_M_0': 0.308, 
            'omega_lambda_0': 0.692, 'h': 0.678

    Attributes:
        lmin (float): log10 of Minimum luminosity considered in erg/s
        lmax (float): log10 of Maximum luminosity considered in erg/s
        _zlocal (float): Describes limit of nearby sources
        Mpc2cm (float): Conversion factor
        GeV_per_sec_2_ergs_per_year (float): Conversion factor
        cosmology (cosmolopy instance)
    """

    def __init__(self, lmin, lmax, cosmology=cosmology):
        """
        Constructor
        """
        self.cosmology = cosmolopy.distance.set_omega_k_0(cosmology)
        self.lmin = lmin
        self.lmax = lmax
        self._zlocal = 0.01
        self.Mpc2cm = 3.086e24                     # Mpc / cm
        self.GeV_per_sec_2_ergs_per_sec = 1.60218e-3  # (GeV/sec) / (ergs/s)

    def LF(self, L, z):
        """
        Luminosity functions should be implemented by inherited classes
        """
        raise NotImplementedError("Please Specify Model")

    def LuminosityDistance(self, z):
        """
        Convert redshift to luminosity distance.
        If passing many redshifts, a 1d spline is used as cosmolopy can be slow

        Args:
            z (array or float): redshift(s)

        Returns:
            array or float: Luminosity distance(s) in Mpc
        """
        # Wrapper function - so that cosmolopy is only imported here.
        if np.ndim(z) > 0:
            if len(z) > 1000:
                zz = np.linspace(0., 10., 500)
                spl = scipy.interpolate.UnivariateSpline(zz,
                        cosmolopy.distance.luminosity_distance(zz,
                            **self.cosmology))
                return spl(z)
        return cosmolopy.distance.luminosity_distance(z, **self.cosmology)

    def RedshiftDistribution(self, z):
        r"""
        Provides the unnormalized PDF of number of sources vs. redshift
        by multiplying the \(\frac{dN}{dz} = \frac{d\rho}{dz} \times \frac{dV}{dz}\), 
        accounting for the luminosity dependence on z

        $$ P(z) = \int_{Lmin}^{Lmax} LF(L,z) \,dL \,dV_c(z) \,4\pi $$ 

        Args:
            z (array or float): Redshift values

        Returns
            Array of float: Unnormalized PDF of number vs. redshift
        """
        integral = scipy.integrate.quad(lambda L: self.LF(L, z), self.lmin, self.lmax)[0]
        return integral * cosmolopy.distance.diff_comoving_volume(z, **self.cosmology) * \
            4*np.pi

    def L_CDF(self, redshift_bins, luminosity_bins):
        """
        Creates a 2-dimensional cumulative distribution function
        of the number of sources as a function of redshift and luminosity

        Args:
            redshift_bins (array): redshift bin-edges for evaluating the 
                CDF
            luminosity_bins (array): luminosity bin-edges for evaluating the 
                CDF 

        Attributes:
            redshift_bins (array): redshift bin-edges for evaluating the 
                CDF
            luminosity_bins (array): luminosity bin-edges for evaluating the 
                CDF
            Lcdf (2d array): 2D CDF of number of sources vs. redshift and 
                luminosity
        """
        # 2D phase space scan of L and z
        l, z = np.meshgrid(luminosity_bins, redshift_bins)
        L_PDF = self.LF(l, z)
        L_CDF = np.cumsum(L_PDF, axis=1)
        norm = L_CDF[:,-1].reshape((len(redshift_bins),1))
        L_CDF = (1/norm) * L_CDF

        self.redshift_bins = redshift_bins
        self.luminosity_bins = luminosity_bins
        self.Lcdf = L_CDF

    def Luminosity_Sampling(self, z):
        """
        Samples luminosities of sources given their redshifts

        Args:
            z (array or float): redshift(s) of source(s)

        Returns:
            array or float: Sampled luminosities
        """
        lumi = []
        z = np.atleast_1d(z)
        test = np.random.rand(z.shape[0])
        index_1 = np.searchsorted(self.redshift_bins, z)
        for test, index in zip(test, index_1):
            index_2 = np.searchsorted(self.Lcdf[index], test)
            lumi.append(self.luminosity_bins[index_2])
        return np.array(lumi)

    def Nsources(self, zmax):
        r"""
        Integrates full 2-dimensional source count distribution over 
            redshift and luminosity

        $$ N_{tot} = \int_0^{z_{max}} P(z) dz$$ 

        Args:
            zmax (float): Maximum redshift to consider

        Returns:
            float: Total number of sources for the luminosity-evolution model
        """
        return scipy.integrate.quad(lambda z: self.RedshiftDistribution(z), 0, zmax)[0]

    def Lumi2Flux(self, luminosity, index, emin, emax, z=1, E0=1e5):
        r"""
        Converts a luminosity to a fluence

        $$ L_\nu = \frac{ \Phi_{z=1}^{PS} }{E_0^2}
        \int_{E_\mathrm{min}}^{E_\mathrm{max}} E
        \left(\frac{E}{E_0}\right)^{-\gamma}\,
        \mathrm{d}E\ \times 4\pi d_L^2(z=1) $$

        Note fluxnorm is E0^2*fluxnorm
        fluence units are [UNITS]

        Args:
            luminosity (array or float): luminosity of sources in ergs/yr
            index (float): Spectral index of the flux
            emin (float): Minimum neutrino energy in GeV
            emax (float): Maximum neutrino energy in GeV
            z (array or float, optional, default=1): Redshifts
            E0 (float, optional, default=1): pivot energy in GeV
        
        Returns:
            array or float: fluence of source(s) in UNITS
        """
        flux_integral = self.EnergyIntegral(index, emin, emax, z, E0)
        fluxnorm = luminosity / 4. / np.pi / \
            (self.LuminosityDistance(z)*self.Mpc2cm)**2. / \
            self.GeV_per_sec_2_ergs_per_sec / flux_integral * E0**2.
        return fluxnorm

    def EnergyIntegral(self, index, emin, emax, z=1, E0=1e5):
        r""" 
        Calculates energy content in a neutrino flux

        $$\int_{emin/(1+z)}^{emax/(1+z)} E*(E/E0)^{-index} dE$$ 
        
        Args:
            index (float): Spectral index of the flux
            emin (float): Minimum neutrino energy in GeV
            emax (float): Maximum neutrino energy in GeV
            z (array or float, optional, default=1): Redshifts
            E0 (float, optional, default=1): pivot energy in GeV

        Returns:
            float: Energy flux between emin and emax
        """
        if index != 2.0:
            denom = (1+z)**(index-2)
            integral = denom * (emax**(2-index)-emin**(2-index)) / (2-index)
        else:
            integral = np.ones_like(z) * np.log(emax/emin)
        return E0**index * integral

class HardingAbazajian(LuminosityEvolution):
    """ 
    Luminosity dependent density evolution for gamma-ray blazars based
        on X-ray AGN luminosity function

    See also: :class:`LuminosityEvolution`

    Reference: arXiv:1206.4734
               arXiv:1012.1247
               arXiv:0308140
    """
    def __str__(self):
        return "Harding and Abazajian (2012)"
        
    def LF(self, L, z):
        """
        Luminosity function based on X-ray AGN

        Args:
            L (float): log10 of luminosity in erg/s
            z (float): redshift

        Returns:
            float: local PDF value of source count vs. luminosity
                and redshift
        """
        A = 5.04e-6
        gamma1 = 0.43
        L0 = 10**43.94
        gamma2 = 2.23
        zc0 = 1.9
        p10 = 4.23
        p20 = -1.5
        alpha = 0.335
        La = 44.6
        beta1 = 0.
        beta2 = 0.
        
        L = np.atleast_1d(L)
        z = np.atleast_1d(z)
        zc = np.zeros_like(L)
        LF_F = np.zeros_like(L)
        # luminosity distribution at z=0
        LF_L = A*((10**L/L0)**gamma1 + (10**L/L0)**gamma2)**-1
        # density indices 1 and 2 --> constant in this model
        p1 = p10 + beta1 * (L-44.0)
        p2 = p20 + beta2 * (L-44.0)
        # zc, where peak evolution happens
        zc[L>=La] = zc0
        zc[L<La] = zc0*10**((L[L<La]-La)*alpha)
        # density evolution
        LF_F[z<zc] = (1+z[z<zc])**p1[z<zc]
        LF_F[z>=zc] = (1+zc[z>=zc])**p1[z>=zc]*((1+z[z>=zc])/(1+zc[z>=zc]))**p2[z>=zc]
        # total evolution
        return LF_L*LF_F

    def Nsources(self, zmax):
        """
        Calculates total number of sources in the universe out to zmax

        Args:
            zmax (float): Maximum redshift to consider

        Returns:
            float: Total number of sources
        """
        kappa = 9.54e-6         #model specific
        nsource = super(HardingAbazajian, self).Nsources(zmax)
        return nsource*kappa

    def Luminosity_Sampling(self, z):
        """
        Samples luminosities of sources given their redshifts, with 
            appropriately applied unit conversion

        Args:
            z (array or float): redshift(s) of source(s)

        Returns:
            array or float: Sampled luminosities
        """
        L_x_to_rad = 4.21           #model specific
        L = super(HardingAbazajian, self).Luminosity_Sampling(z)
        return 10**(L+L_x_to_rad)
