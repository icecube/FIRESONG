#!/usr/bin/python

import numpy as np
import scipy
import cosmolopy
cosmology = {'omega_M_0': 0.308, 'omega_lambda_0': 0.692, 'h': 0.678}


def get_evolution(evol):
    evolutions = {"NoEvolution": NoEvolution,
                  "HB2006SFR": HopkinsBeacom2006StarFormationRate,
                  "YMKBH2008SFR": YukselEtAl2008StarFormationRate,
                  "CC2015SNR": CandelsClash2015SNRate,
                  "MD2014SFR": MadauDickinson2014CSFH
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
    """This is the implied SRF from Goods/Candels/Clash (2015)
    derive from CC SNe rate and assuming one rate is proportional to the other.
    They use the same functional form as Madau and Dickinson (2014)
    unit = M_sun/yr/Mpc^3 """
    
    def __call__(self, z):
        return self.parametrization(1.+z)
    
    def parametrization(self, x):
        a = 0.015
        b = 1.5
        c = 5.0
        d = 6.1
        density = a*(x**c) / (1. + ( x / b)**d )
        return density


class MadauDickinson2014CSFH(Evolution):
    """ StarFormationHistory (SFR), from Madau and Dickinson (2014),
    unit = M_sun/yr/Mpc^3 """

    def __call__(self, z):
        return self.parametrization(1.+z)
        
    def parametrization(self, x):
        a = 0.015
        b = 2.7
        c = 2.9
        d = 5.6
        density = a*(x**b) / (1. + (x/c)**d ) 
        return density


class SourcePopulation(object):
    def __init__(self, cosmology, evolution):
        self._zlocal = 0.01
        self.Mpc2cm = 3.086e24                     # Mpc / cm
        self.GeV_per_sec_2_ergs_per_year = 50526.  # (GeV/sec) / (ergs/yr)
        self.evolution = evolution

        # Flat universe
        self.cosmology = cosmolopy.distance.set_omega_k_0(cosmology)

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

    def Flux2Lumi(self, fluxnorm, index, emin, emax, z=1, E0=1e5):
        """
        $$ L_\nu = \frac{ \Phi_{z=1}^{PS} }{E_0^2}
        \int_{E_\mathrm{min}}^{E_\mathrm{max}} E
        \left(\frac{E}{E_0}\right)^{-\gamma}\,
        \mathrm{d}E\,4\pi d_L^2(z=1) $$

        Note fluxnorm is E0^2*fluxnorm
        fluxnorm units are []
        """
        flux_integral = self.EnergyIntegral(index, emin, emax, z, E0)
        luminosity = fluxnorm / E0**2. * flux_integral *  \
            self.GeV_per_sec_2_ergs_per_year * \
            4. * np.pi * (self.LuminosityDistance(z)*self.Mpc2cm)**2.
        return luminosity

    def Lumi2Flux(self, luminosity, index, emin, emax, z=1, E0=1e5):
        """
        $$ L_\nu = \frac{ \Phi_{z=1}^{PS} }{E_0^2}
        \int_{E_\mathrm{min}}^{E_\mathrm{max}} E
        \left(\frac{E}{E_0}\right)^{-\gamma}\,
        \mathrm{d}E\,4\pi d_L^2(z=1) $$

        Lumi given in ergs/yr
        Note fluxnorm is E0^2*fluxnorm
        fluxnorm units are []
        """
        flux_integral = self.EnergyIntegral(index, emin, emax, z, E0)
        fluxnorm = luminosity / 4. / np.pi / \
            (self.LuminosityDistance(z)*self.Mpc2cm)**2. / \
            self.GeV_per_sec_2_ergs_per_year / flux_integral * E0**2.
        return fluxnorm

    def EnergyIntegral(self, index, emin, emax, z=1, E0=1e5):
        """ integal_{emin/(1+z)}^{emax/(1+z)} E*(E/E0)^(-index) dE """
        l_lim = emin/(1.+z)
        u_lim = emax/(1.+z)
        if index != 2.0:
            integral = (u_lim**(2-index)-l_lim**(2-index)) / (2-index)
        else:
            integral = np.log(u_lim) - np.log(l_lim)
        return E0**index * integral

    def StandardCandleSources(self, fluxnorm, density, zmax, index, z0=1.):
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
        Fluxnorm = all_sky_flux / Ntotal / self.LuminosityDistance(z0)**2. / \
            scipy.integrate.quad(lambda z: ((1.+z)/(1.+z0))**(2-index) /
                                 self.LuminosityDistance(z)**2. *
                                 self.RedshiftDistribution(z) / norm,
                                 0, 10.)[0]

        return Fluxnorm

    def StandardCandleLuminosity(self, fluxnorm, density, zmax, index,
                                 emin, emax, E0=1e5):
        """ """

        flux = self.StandardCandleSources(fluxnorm, density, zmax, index, z0=1)
        luminosity = self.Flux2Lumi(flux, index, emin, emax, z=1, E0=E0)
        return luminosity


class TransientSourcePopulation(SourcePopulation):
    def __init__(self, cosmology, evolution, timescale):
        super(TransientSourcePopulation, self).__init__(cosmology, evolution)
        self.timescale = timescale
        self.yr2sec = 86400*365

    def RedshiftDistribution(self, z):
        return super(TransientSourcePopulation, self).RedshiftDistribution(z) / (1.+z)

    def StandardCandleSources(self, fluxnorm, density, zmax, index, z0=1.):
        # For transient source, Fluxnorm will be the fluence of a
        # standard candle at z=1, with unit GeV/cm^2 given that the
        # burst rate density is measured in per year.

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
        luminosity = super(TransientSourcePopulation, self).Flux2Lumi(fluxnorm,
                                                                      index,
                                                                      emin,
                                                                      emax,
                                                                      z=z,
                                                                      E0=E0)
        return luminosity / self.timescale

    def Lumi2Flux(self, luminosity, index, emin, emax, z=1, E0=1e5):
        flux = super(TransientSourcePopulation, self).Lumi2Flux(luminosity,
                                                                index,
                                                                emin,
                                                                emax,
                                                                z=z,
                                                                E0=E0)
        return flux * self.timescale

    def fluence2flux(self, fluence, z):
        # For transient sources, the flux measured on Earth will be
        # red-shifted-fluence/{(1+z)*burst duration}
        flux = fluence / ((1.+z)*self.timescale)
        return flux

#############
#LEGEND AREA#
#############

def get_LEvolution(le_model, lmin, lmax):
    evolutions = {"HA2014BL": HardingAbazajian(lmin, lmax)
                  }
    if not le_model in evolutions.keys():
        raise NotImplementedError("Luminosity Evolution " +
                                  le_model + " not implemented.")
    return evolutions[le_model]

class LuminosityEvolution(object):
    """
    For the Luminosity Distribution that depends on z
    """

    def __init__(self, lmin, lmax, cosmology=cosmology):
        self.cosmology = cosmolopy.distance.set_omega_k_0(cosmology)
        self.lmin = lmin
        self.lmax = lmax
        self._zlocal = 0.01
        self.Mpc2cm = 3.086e24                     # Mpc / cm
        self.GeV_per_sec_2_ergs_per_sec = 1.60218e-3  # (GeV/sec) / (ergs/s)

    def LF(self, L, z):
        raise NotImplementedError("Please Specify Model")

    def LuminosityDistance(self, z):
        # Wrapper function - so that cosmolopy is only imported here.
        return cosmolopy.distance.luminosity_distance(z, **self.cosmology)

    def RedshiftDistribution(self, z):
        """
        $$ P(z) = \int_{Lmin}^{Lmax} LF(L,z) \,dL \,dV_c(z) \,4\pi $$ 
        """
        integral = scipy.integrate.quad(lambda L: self.LF(L, z), self.lmin, self.lmax)[0]
        return integral * cosmolopy.distance.diff_comoving_volume(z, **self.cosmology) * \
            4*np.pi

    def L_CDF(self, redshift_bins, luminosity_bins):
        # 2D phase space scan of L and z
        l, z = np.meshgrid(luminosity_bins, redshift_bins)
        L_PDF = self.LF(l, z)
        L_CDF = np.cumsum(L_PDF, axis=1)
        norm = L_CDF[:,-1].reshape((10000,1))
        L_CDF = (1/norm) * L_CDF

        self.redshift_bins = redshift_bins
        self.luminosity_bins = luminosity_bins
        self.Lcdf = L_CDF

    def Luminosity_Sampling(self, z):
        z = np.atleast_1d(z)
        test = np.random.rand(z.shape[0])
        index_1 = np.searchsorted(self.redshift_bins, z)
        for test, index in zip(test, index_1):
            index_2 = np.searchsorted(self.Lcdf[index], test)
        return self.luminosity_bins[index_2]

    def Nsources(self, zmax):
        """
        $$ N_{tot} = \int_0^z_{max} P(z) dz$$ 
        """
        return scipy.integrate.quad(lambda z: self.RedshiftDistribution(z), 0, zmax)[0]

    def Lumi2Flux(self, luminosity, index, emin, emax, z=1, E0=1e5):
        """
        $$ L_\nu = \frac{ \Phi_{z=1}^{PS} }{E_0^2}
        \int_{E_\mathrm{min}}^{E_\mathrm{max}} E
        \left(\frac{E}{E_0}\right)^{-\gamma}\,
        \mathrm{d}E\,4\pi d_L^2(z=1) $$
        Lumi given in ergs/yr
        Note fluxnorm is E0^2*fluxnorm
        fluxnorm units are []
        """
        flux_integral = self.EnergyIntegral(index, emin, emax, z, E0)
        fluxnorm = luminosity / 4. / np.pi / \
            (self.LuminosityDistance(z)*self.Mpc2cm)**2. / \
            self.GeV_per_sec_2_ergs_per_sec / flux_integral * E0**2.
        return fluxnorm

    def EnergyIntegral(self, index, emin, emax, z=1, E0=1e5):
        """ integal_{emin/(1+z)}^{emax/(1+z)} E*(E/E0)^(-index) dE """
        l_lim = emin/(1.+z)
        u_lim = emax/(1.+z)
        if index != 2.0:
            integral = (u_lim**(2-index)-l_lim**(2-index)) / (2-index)
        else:
            integral = np.log(u_lim) - np.log(l_lim)
        return E0**index * integral

class HardingAbazajian(LuminosityEvolution):
    """
    Luminosity Evolution model for Blazars
    """
    def LF(self, L, z):
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
        LF_L = A*((10**L/L0)**gamma1 + (10**L/L0)**gamma2)**-1
        p1 = p10 + beta1 * (L-44.0)
        p2 = p20 + beta2 * (L-44.0)
        zc[L>=La] = zc0
        zc[L<La] = zc0*10**((L[L<La]-La)*alpha)
        LF_F[z-zc<=0] = (1+z[z-zc<=0])**p1[z-zc<=0]
        LF_F[z-zc>0] = (1+zc[z-zc>0])**p1[z-zc>0]*((1+z[z-zc>0])/(1+zc[z-zc>0]))**p2[z-zc>0]
        return LF_L*LF_F

    def Nsources(self, zmax):
        kappa = 9.54e-6
        nsource = super(HardingAbazajian, self).Nsources(zmax)
        return nsource*kappa

    def Luminosity_Sampling(self, z):
        L_x_to_rad = 4.21
        L = super(HardingAbazajian, self).Luminosity_Sampling(z)
        return 10**(L+L_x_to_rad)
