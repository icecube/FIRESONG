import numpy as np

def PowerLawFlux(meanflux, index, nsource, width):
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

	return flux

def LognormalFlux(meanflux, nsource, width):
	logmean = np.log(meanflux)
	sigma = np.log(10**width)
	lognormaldistribution = np.exp(np.random.normal(logmean-0.5*(sigma**2), sigma, nsource))

	return lognormaldistribution
