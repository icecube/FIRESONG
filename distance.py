#!/usr/bin/python

'''
Distance calculation substitute for cosmolopy
All formula used are from 
https://cds.cern.ch/record/387177/files/9905116.pdf

Assumed flat universe i.e. omega_R = 0.
'''
import numpy as np
from scipy.integrate import quad

def E(z, **kwargs):
    '''
    time derivative of the logarithm of the scale factor
    E(z) = \frac{\dot{a}(t)}{a(t)}
    '''
    om0 = kwargs.get('omega_M_0', 0.315)
    ode = kwargs.get('omega_lambda_0', 0.685)
    EE = np.sqrt(om0*(1.+z)**3.+ ode)
    return EE

def comoving_distance(z,z_0=0, **kwargs):
    '''
    line-of-sight comoving distance
    D_c = D_h \int^z_0 \frac{dz'}{E(z')}
    unit = Mpc
    '''
    h = kwargs.get('h', 0.674)
    # Hubble distance D_h (unit=Mpc)
    D_h = 299792458./1000./100./h

    z = np.atleast_1d(z)
    D_c = np.array([D_h*quad(lambda x:1./E(x, **kwargs), z_0, lim, epsabs=1.e-5)[0] for lim in z])
    if np.size(D_c) > 1:
        return D_c
    else:
        return D_c[0]

def angular_diameter_distance(z, **kwargs):
    '''
    ratio of an object’s physical transverse size to its angular size
    D_A = \frac{D_M}{1+z} = \frac{D_c}{1+z}
    '''
    return comoving_distance(z, **kwargs)/(1.+z)

def luminosity_distance(z, **kwargs):
    '''
    for a flat universe
    D_L = (1+z)*D_M = (1+z)*D_c
    unit = Mpc
    '''
    return (1.+z)*comoving_distance(z, **kwargs)

def comoving_volume(z, **kwargs):
    '''
    comoving volume up to redshift z
    for flat universe
    V_c = \frac{4\pi}{3}D_M^3 = \frac{4\pi}{3}D_c^3
    unit = Mpc^3
    '''
    return 4*np.pi/3 * comoving_distance(z, **kwargs)**3.

def diff_comoving_volume(z, **kwargs):
    '''
    differential comoving volume at redshift z
    dV_c = D_h * \frac{(1+z)^2 * D_A^2}{E(z)} d\Omega dz
    unit = Mpc^3 sr^{-1} 
    '''
    h = kwargs.get('h', 0.674)
    # Hubble distance D_h (unit=Mpc)
    D_h = 299792458./1000./100./h

    dV_c = D_h*(1.+z)**2.*angular_diameter_distance(z, **kwargs)**2./E(z, **kwargs)
    return dV_c




