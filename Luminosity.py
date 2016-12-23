#!/usr/bin/python
#

import numpy as np


# Is this a real luminosty ?
# It'd be nice to have something with the correct units here
def LuminosityFunction(options, nsource, candleflux):
    if options.LF == "SC":
        return candleflux
    if options.LF == "LG":
        return LognormalFlux(candleflux, nsource, options.sigma)
    if options.LF == "PL":
        return PowerLawFlux(candleflux, -2, nsource, options.sigma)

    
def PowerLawFlux(meanflux, index, nsource, width):
## index = -2, f_mean ~= 2ln(10)*F_min
## index < -2, f_mean ~= (index+1)/(index+2)*F_min
        if index == -2:
                width = 10**width
                F_min = meanflux / (width*np.log(10))
                F_max = F_min*10**width
                
                F_bin  = np.linspace(F_min, F_max, 10000)
                pdf = F_bin**index
                cdf = np.cumsum(pdf)
                cdf = cdf / cdf[-1]
                test = np.random.rand(nsource)
                F_id = np.searchsorted(cdf, test)
                flux = np.array([F_bin[i] for i in F_id])
        if index < -2:
                width = 10**width
                F_min = (index+2.)/(index+1.)*meanflux
                F_max = F_min*10**width

                F_bin  = np.linspace(F_min, F_max, 10000)
                pdf = F_bin**index
                cdf = np.cumsum(pdf)
                cdf = cdf / cdf[-1]
                test = np.random.rand(nsource)
                F_id = np.searchsorted(cdf, test)
                flux = np.array([F_bin[i] for i in F_id])
        return flux

def LognormalFlux(meanflux, nsource, width):
        logmean = np.log(meanflux)
        sigma = np.log(10**width)
        lognormaldistribution = np.exp(np.random.normal(logmean-0.5*(sigma**2), sigma, nsource))
        return lognormaldistribution
