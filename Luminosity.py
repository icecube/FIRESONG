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
                F_min = meanflux / (width*np.log(10))
                F_max = F_min*10**width
                
                F_bin  = np.linspace(F_min, F_max, 10000)
                pdf = F_bin**index
                cdf = np.cumsum(pdf)
                cdf = cdf / cdf[-1]
                test = np.random.rand(nsource)
                F_id = np.searchsorted(cdf, test)
                flux = np.array([F_bin[i] for i in F_id])
                return flux

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

## Not really useful at this moment
def LuminosityPDF(options, meanflux):
    if options.LF == "SC":
        return 1
    if options.LF == "LG":
        sigma = np.log(10**options.sigma)
        mu = np.log(meanflux)-sigma**2./2.
        f1_list = np.exp(np.linspace(mu-3*sigma, mu+3*sigma, 10000))
        pdf = np.exp(-(np.log(f1_list)-mu)**2./(2.*sigma**2.))/sigma/np.sqrt(2*np.pi)/f1_list
    if options.LF == "PL":
        #for power-law index = -2 only !!!!!!!
        F_min = meanflux / (options.sigma*np.log(10))
        F_max = F_min*10**options.sigma
        f1_list = 10**np.linspace(np.log10(F_min), np.log10(F_max), 10000)
        pdf = f1_list**-2

    nPDF_f = f1_list*pdf
    nCDF_f = np.cumsum(nPDF_f)
    nCDF_f = nCDF_f / nCDF_f[-1]
    return f1_list, nCDF_f

    
