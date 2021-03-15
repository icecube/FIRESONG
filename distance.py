#!/usr/bin/python

'''
Distance calculation substitute for cosmolopy
All formula used are from 
https://arxiv.org/pdf/astro-ph/9905116.pdf

Default in all functions are using Planck15
'''
import numpy as np
from scipy.integrate import quad

c = 299792458 # speed of light in ms^-1

def E(z, **kwargs):
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
    om0 = kwargs.get('omega_M_0', 0.315)
    ode = kwargs.get('omega_lambda_0', 0.685)
    ok0 = kwargs.get('omega_k_0', 1-om0-ode)
    EE = np.sqrt(om0*(1.+z)**3.+ok0*(1.+z)**2.+ode)
    return EE

def comoving_distance(z,z_0=0, **kwargs):
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
    h = kwargs.get('h', 0.674)
    # Hubble distance D_h (unit=Mpc)
    D_h = c/1000./100./h

    z = np.atleast_1d(z)
    D_c = np.array([D_h*quad(lambda x:1./E(x, **kwargs), 
                             z_0, lim, epsabs=1.e-5)[0] for lim in z])
    # return a float if only one redshift
    if np.size(D_c) > 1:
        return D_c
    else:
        return D_c[0]

def transverse_comoving_distance(z, **kwargs):
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
    h = kwargs.get('h', 0.674)
    om0 = kwargs.get('omega_M_0', 0.315)
    ode = kwargs.get('omega_lambda_0', 0.685)
    ok0 = kwargs.get('omega_k_0', 1-om0-ode)
    # Hubble distance D_h (unit=Mpc)
    D_h = c/1000./100./h

    if ok0==0:
        return comoving_distance(z, **kwargs)
    elif ok0>0:
        return D_h/np.sqrt(ok0)*np.sinh(np.sqrt(ok0)*comoving_distance(z, **kwargs)/D_h)
    else:
        return D_h/np.sqrt(np.abs(ok0))*np.sin(np.sqrt(np.abs(ok0))\
               *comoving_distance(z, **kwargs)/D_h)


def angular_diameter_distance(z, **kwargs):
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
    return transverse_comoving_distance(z, **kwargs)/(1.+z)

def luminosity_distance(z, **kwargs):
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
    return (1.+z)*transverse_comoving_distance(z, **kwargs)

def comoving_volume(z, **kwargs):
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
    h = kwargs.get('h', 0.674)
    om0 = kwargs.get('omega_M_0', 0.315)
    ode = kwargs.get('omega_lambda_0', 0.685)
    ok0 = kwargs.get('omega_k_0', 1-om0-ode)
    # Hubble distance D_h (unit=Mpc)
    D_h = c/1000./100./h

    if ok0 == 0:
        return 4*np.pi/3 * comoving_distance(z, **kwargs)**3.
    elif ok0>0:
        D_m = transverse_comoving_distance(z, **kwargs)
        D_ratio = D_m/D_h
        return 4.*np.pi*D_h**3./2./ok0*(D_ratio*np.sqrt(1+ok0*D_ratio**2.)\
        -np.arcsinh(np.sqrt(np.abs(ok0))*D_ratio)/np.sqrt(np.abs(ok0)))
    else:
        D_m = transverse_comoving_distance(z, **kwargs)
        D_ratio = D_m/D_h
        return 4.*np.pi*D_h**3./2./ok0*(D_ratio*np.sqrt(1+ok0*D_ratio**2.)\
         -np.arcsin(np.sqrt(np.abs(ok0))*D_ratio)/np.sqrt(np.abs(ok0)))

def diff_comoving_volume(z, **kwargs):
    '''
    differential comoving volume at redshift z
    dV_c = D_h * \frac{(1+z)^2 * D_A^2}{E(z)} d\Omega dz
    unit = Mpc^3 sr^{-1} 
    '''
    h = kwargs.get('h', 0.674)
    # Hubble distance D_h (unit=Mpc)
    D_h = c/1000./100./h

    dV_c = D_h*(1.+z)**2.*angular_diameter_distance(z, **kwargs)**2./E(z, **kwargs)
    return dV_c




