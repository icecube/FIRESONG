#!/usr/bin/python

'''
Distance calculation substitute for cosmolopy
All formula used are from 
https://arxiv.org/pdf/astro-ph/9905116.pdf
'''
import numpy as np
from scipy.integrate import quad

class cosmo_distance(object):
    def __init__(self, **cosmology):
        '''
        To initiate, cosmological parameters must be supplied
        in the form of a dictionary with the following keys
        {'omega_M_0', 'omega_lambda_0', 'h'}

        Input:
        **kwargs = {'omega_M_0', 'omega_lambda_0', 'h'}
        '''
        self.c = 299792458 # speed of light in ms^-1
        # initialize class
        for k in ['omega_M_0', 'omega_lambda_0', 'h']:
            if k not in cosmology.keys():
                raise Exception('Cosmological parameter {} must be supplied'.format(k))
        self.load_param(**cosmology)

    def load_param(self, **cosmology):
        '''
        set up the cosmological parameters
        density of curvature is defined as 1-om0-ode
        unless otherwise specified.
        '''
        self.om0 = cosmology.get('omega_M_0')
        self.ode = cosmology.get('omega_lambda_0')
        self.ok0 = cosmology.get('omega_k_0', 1-self.om0-self.ode)
        self.h = cosmology.get('h')

        # Hubble distance D_h (unit=Mpc)
        self.D_h = self.c/1000./100./self.h

    def E(self, z):
        '''
        time derivative of the logarithm of the scale factor
        E(z) = \frac{\dot{a}(t)}{a(t)}

        Input:
        z     : float or array of float, redshift
        kwargs:  'omega_M_0'     
                 'omega_lambda_0'
                 'omega_k_0'
                 density parameters
        output:
        E      : float or array of float
        '''
        EE = np.sqrt(self.om0*(1.+z)**3.+self.ok0*(1.+z)**2.+self.ode)
        return EE

    def comoving_distance(self, z, z_0=0):
        '''
        line-of-sight comoving distance
        D_c = D_h \int^z_0 \frac{dz'}{E(z')}
        unit = Mpc
        Input:
        z     : float or array of float, redshift
        z0    : float, default is 0, where integration begins
        kwargs:  'omega_M_0'     
                 'omega_lambda_0'
                 'omega_k_0'
                 density parameters
        output:
        D_c      : float or array of float

        '''
        z = np.atleast_1d(z)
        # epsabs can be tweaked to achieve higher precision
        D_c = np.array([self.D_h*quad(lambda x:1./self.E(x), 
                                 z_0, lim, epsabs=1.e-5)[0] for lim in z])
        # return a float if only one redshift
        if np.size(D_c) > 1:
            return D_c
        else:
            return D_c[0]

    def transverse_comoving_distance(self, z):
        '''
        The comoving distance between two events at the same redshift or distance 
        but separated on the sky by some angle is D_m * d\theta
        D_m is the transverse comoving distance
        if omega_R > 0:
            D_m = D_h /\sqrt(Omega_R) \sinh(\sqrt(Omega_R)*D_c/D_h)
        if omega_R = 0:
            D_m = D_c
        if omega_R < 0:
            D_m = D_h /\sqrt(\abs(Omega_R)) \sin(\sqrt(\abs(Omega_R))*D_c/D_h)

        unit = Mpc

        Input:
        z     : float or array of float, redshift
        kwargs:  'omega_M_0'     
                 'omega_lambda_0'
                 'omega_k_0'
                 density parameters
        output:
        D_m      : float or array of float
        '''
        if self.ok0==0:
            return self.comoving_distance(z)
        elif self.ok0>0:
            return self.D_h/np.sqrt(self.ok0)*np.sinh(np.sqrt(self.ok0)\
                   *self.comoving_distance(z)/self.D_h)
        else:
            return self.D_h/np.sqrt(np.abs(self.ok0))*np.sin(np.sqrt(np.abs(self.ok0))\
                   *self.comoving_distance(z)/self.D_h)


    def angular_diameter_distance(self, z):
        '''
        ratio of an object’s physical transverse size to its angular size
        D_A = \frac{D_M}{1+z}
        Input:
        z     : float or array of float, redshift
        kwargs:  'omega_M_0'     
                 'omega_lambda_0'
                 'omega_k_0'
                 density parameters
        output:
        D_A      : float or array of float
        '''
        return self.transverse_comoving_distance(z)/(1.+z)

    def luminosity_distance(self, z):
        '''
        for a flat universe
        D_L = (1+z)*D_M
        unit = Mpc

        Input:
        z     : float or array of float, redshift
        kwargs:  'omega_M_0'     
                 'omega_lambda_0'
                 'omega_k_0'
                 density parameters
        output:
        D_L      : float or array of float
        '''
        return (1.+z)*self.transverse_comoving_distance(z)

    def comoving_volume(self, z):
        '''
        comoving volume up to redshift z
        for flat universe
        V_c = \frac{4\pi}{3}D_M^3 = \frac{4\pi}{3}D_c^3
        other definitions please refer to the reference
        unit = Mpc^3

        Input:
        z     : float or array of float, redshift
        kwargs:  'omega_M_0'     
                 'omega_lambda_0'
                 'omega_k_0'
                 density parameters
        output:
        V_c      : float or array of float
        '''
        if self.ok0 == 0:
            return 4*np.pi/3 * self.comoving_distance(z)**3.
        elif ok0>0:
            D_m = self.transverse_comoving_distance(z)
            D_ratio = D_m/self.D_h
            prefactor = 4.*np.pi*self.D_h**3./2./self.ok0
            return prefactor*(D_ratio*np.sqrt(1+self.ok0*D_ratio**2.)\
            -np.arcsinh(np.sqrt(np.abs(self.ok0))*D_ratio)/np.sqrt(np.abs(self.ok0)))
        else:
            D_m = transverse_comoving_distance(z)
            D_ratio = D_m/D_h
            return prefactor*(D_ratio*np.sqrt(1+self.ok0*D_ratio**2.)\
             -np.arcsin(np.sqrt(np.abs(self.ok0))*D_ratio)/np.sqrt(np.abs(self.ok0)))

    def diff_comoving_volume(self, z):
        '''
        differential comoving volume at redshift z
        dV_c = D_h * \frac{(1+z)^2 * D_A^2}{E(z)} d\Omega dz
        unit = Mpc^3 sr^{-1} 
        '''
        dV_c = self.D_h*(1.+z)**2.*self.angular_diameter_distance(z)**2.\
               /self.E(z)
        return dV_c




