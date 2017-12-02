## This is module for using blazar's luminosity evolution by Harding and Abazajian
## Any other model can follow this as a template
## variables that will be read in LEGEND.py are LF(L,z) and kappa
## LF is the distribution of source as a function of L and z
## kappa is the fraction of source we are interested in the LF

## !!Important!! ##
## L is X-ray luminosity in this model

A = 5.04e-6
gamma1 = 0.43
L0 = 10**43.94
gamma2 = 2.23
zc0 = 1.9
p10 = 4.23
p20 = -1.5
alpha = 0.335
La = 10**44.6
beta1 = 0.
beta2 = 0.
kappa = 9.54e-6
L_x_to_rad = 4.21

def LF_ueda_L(L):
	return A*((10**L/L0)**gamma1 + (10**L/L0)**gamma2)**-1

def zc(L):
	if 10**L >= La:
		return zc0
	if 10**L < La:
		return zc0*(10**L/La)**alpha

def p1(L):
	return p10 + beta1 * (L-44.0)

def p2(L):
	return p20 + beta2 * (L-44.0)

def LF_ueda_f(L,z):
	if z <= zc(L):
		return (1+z)**p1(L)
	if z > zc(L):
		return (1+zc(L))**p1(L)*((1+z)/(1+zc(L)))**p2(L)

def LF(L,z):
	return LF_ueda_f(L,z)*LF_ueda_L(L)
