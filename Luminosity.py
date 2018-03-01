#!/usr/bin/python
#

import numpy as np
from scipy.stats import lognorm


class LuminosityFunction():
    def __init__(self, candleflux):
        self.meanflux = candleflux              # mean flux

    def sample_distribution(self, nsources, rng=None):
        raise NotImplementedError("Abstract Class")

    def pdf(self):
        raise NotImplementedError("Abstract Class")

    def cdf(self):
        raise NotImplementedError("Abstract Class")

    def cdf_eval(self):
        raise NotImplementedError("Abstract Class")


class SC_LuminosityFunction(LuminosityFunction):

    def sample_distribution(self, nsources=None, rng=None):
        return self.meanflux

    def pdf(self, lumi):
        raise NotImplementedError("""This is a delta function.
            PDF not implemented.""")

    def cdf(self, lumi):
        raise NotImplementedError("""This is a delta function.
            CDF not implemented.""")

    def cdf_eval(self):
        return 1


class LG_LuminosityFunction(LuminosityFunction):

    def __init__(self, candleflux, width):
        self.meanflux = candleflux
        self.width = width                     # width is given in log10
        self.logmean = np.log(self.meanflux)   # log mean flux
        self.sigma = np.log(10**self.width)    # sigma is given in ln
        self.mu = self.logmean-self.sigma**2./2.  # log median flux

    def sample_distribution(self, nsource=None, rng=None):
        """ Samples from the Luminosity Function nsource times

        Parameters:
            number of sources
        """
        if rng is None:
            rng = np.random.RandomState()
        return rng.lognormal(self.mu, self.sigma, nsource)

    def pdf(self, lumi):
        """ Gives the value of the PDF at lumi.

        Parameters:
            lumi: float or array-like, point where PDF is evaluated.

        Notes:
            PDF given by:
                     1                 /     (ln(x) - mu)^2   \
            -------------------- * exp | -  ----------------  |
             x sigma sqrt(2 pi)        \       2 sigma^2      /
        """
        return lognorm.pdf(lumi, s=self.sigma, scale=np.exp(self.mu))

    def cdf(self, lumi):
        """ Gives the value of the CDF at lumi.

        Parameters:
            lumi: float or array-like, point where CDF is evaluated.

        Notes:
            CDF given by:
             1     1       /  (ln(x) - mu)^2   \
            --- + --- erf |  ----------------   |
             2     2       \   sqrt(2) sigma   /
        """
        return lognorm.cdf(lumi, s=self.sigma, scale=np.exp(self.mu))

    def cdf_eval(self):
        """ Evaluates CDF at within a `central` range.
        Returns luminosity list and CDF values
        corresponding to these luminosity values.
        """

        f1_list = np.exp(np.linspace(self.mu-3*self.sigma,
                                     self.mu+3*self.sigma, 10000))
        return f1_list, self.cdf(f1_list)


class PL_LuminosityFunction(LuminosityFunction):
    """
        PDF:
                       (1-alpha)
            ------------------------------------ * x^(-alpha)
            (x_max^(1-alpha) - x_min^(1-alpha))


        CDF:
              x^(1-alpha) - x_min^(1-alpha)
            ---------------------------------
            x_max^(1-alpha) - x_min^(1-alpha)


        inv.CDF:
            P^-1(x) = (x_min^beta + (x_max^beta -x_min^beta)*x)^(1/beta)


        Formulars from: http://up-rs-esp.github.io/bpl/
        and: Clauset, A., Shalizi, C. R., Newman, M. E. J.
        "Power-law Distributions in Empirical Data".
        SFI Working Paper: 2007-12-049 (2007)  arxiv:0706.1062
    """

    def __init__(self, candleflux, index, width):
        """
            index = -2, f_mean ~= width*ln(10)*F_min
            index < -2, f_mean ~= (index+1)/(index+2)*F_min
        """
        self.meanflux = candleflux
        self.index = index
        self.width = width

        if self.index == -2:
            self.Fmin = self.meanflux / (self.width*np.log(10))
        else:
            self.Fmin = (self.index+2.)/(self.index+1.)*self.meanflux

        self.Fmax = self.Fmin*10**self.width

    def sample_distribution(self, nsource=None, rng=None):
        """
        inv.CDF:
        a = 1-index
        P^{-1}(x)=(x_min^a+(x_max^a-x_min^a)*x)^(1/a)
        """
        if rng is None:
            rng = np.random.RandomState()
        x = rng.uniform(0, 1, nsource)
        beta = (1.+self.index)
        return (self.Fmin**beta + (self.Fmax**beta -
                                   self.Fmin**beta)*x)**(1./beta)

    def pdf(self, lumi):
        """
        PDF:
                       (1-alpha)
            ------------------------------------ * x^(-alpha)
            (x_max^(1-alpha) - x_min^(1-alpha))
        """

        norm = (1+self.index)/(self.Fmax**(1+self.index) -
                               self.Fmin**(1+self.index))
        pdf = norm*lumi**self.index
        pdf[np.logical_or(lumi < self.Fmin, lumi > self.Fmax)] = 0
        return pdf

    def cdf(self, lumi):
        """
        CDF:
             x^(1-alpha) - x_min^(1-alpha)
            ---------------------------------
            x_max^(1-alpha) - x_min^(1-alpha)
        """
        return (lumi**(1+self.index) - self.Fmin**(1+self.index)) / \
               ((self.Fmax**(1+self.index) - self.Fmin**(1+self.index)))

    def cdf_eval(self):
        f1_list = np.logspace(np.log10(self.Fmin), np.log10(self.Fmax), 10000)
        return f1_list, self.cdf(f1_list)


# Is this a real luminosty ?
# It'd be nice to have something with the correct units here
def get_LuminosityFunction(options, candleflux):
    if options.LF == "SC":
        return SC_LuminosityFunction(candleflux)
    if options.LF == "LG":
        return LG_LuminosityFunction(candleflux,
                                     options.sigma)
    if options.LF == "PL":
        return PL_LuminosityFunction(candleflux,
                                     options.index,
                                     options.sigma)
    raise NotImplementedError("The luminosity function " +
                              options.LF + " is not implemented.")
