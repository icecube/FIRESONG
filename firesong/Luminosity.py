#!/usr/bin/python

"""Classes to help with various luminosity functions"""

# imports
import numpy as np
from scipy.stats import lognorm


class LuminosityFunction(object):
    """Luminosity Function class.
    """

    def __init__(self, mean_luminosity):
        """Luminosity function.

        Args:
            - mean luminosity (Not not median and not log(mean))
        """
        self.mean = mean_luminosity

    def sample_distribution(self, nsources, rng=None):
        raise NotImplementedError("Abstract Class")

    def pdf(self, lumi):
        raise NotImplementedError("Abstract Class")

    def cdf(self, lumi):
        raise NotImplementedError("Abstract Class")


class SC_LuminosityFunction(LuminosityFunction):
    """ Standard Candle Luminosity functions.
    All sources have same luminosity.
    """

    def sample_distribution(self, nsources=None, rng=None):
        """ Samples from the Luminosity Function nsource times

        Args:
            number of sources
        """
        return self.mean

    def pdf(self, lumi):
        """ Gives the value of the PDF at lumi.

        Args:
            lumi: float or array-like, point where PDF is evaluated.

        Notes:
            PDF for SC not defined. This just works if an array is given.
            PDF is zero every where except where bin edges include the SC flux.
        """
        lumi = np.atleast_1d(lumi)
        if len(lumi) < 2:
            raise NotImplementedError("""This is a delta function.
                PDF not defined.
                We only can return something meaning full for a nice array.
                """)
        pdf = np.zeros_like(lumi)
        pdf[:-1][np.logical_and(self.mean > lumi[:-1],
                 self.mean < lumi[1:])] = 1
        return pdf

    def cdf(self, lumi):
        """ Gives the value of the CDF at lumi.

        Args:
            lumi: float or array-like, point where CDF is evaluated.

        Notes:
            0    lumi < mean
            0.5  lumi = mean
            1    lumi > mean
        """
        lumi = np.atleast_1d(lumi)
        cdf = np.zeros_like(lumi)
        cdf[lumi == self.mean] = .5
        cdf[lumi > self.mean] = 1.
        if len(lumi) == 1:
            return cdf[0]
        return cdf


class LG_LuminosityFunction(LuminosityFunction):

    def __init__(self, mean_luminosity, width):
        """ Log Normal Luminosity function.

        Args:
            - mean luminosity (Not median and not log(mean))
            - width in dex (width in log10)

        Note:
            Relation between mean and median luminosity:
            log(median) = log(mean) - sigma^2/2
            Relation between width and simga:
            sigma = log(10^width)
        """
        self.mean = mean_luminosity
        self.width = width                     # width is given in log10
        self.logmean = np.log(self.mean)       # log mean luminosity
        self.sigma = np.log(10**self.width)    # sigma is given in ln
        self.mu = self.logmean-self.sigma**2./2.  # log median luminosity

    def sample_distribution(self, nsources=None, rng=None):
        """ Samples from the Luminosity Function nsource times

        Args:
            number of sources
        """
        if rng is None:
            rng = np.random.RandomState()
        return rng.lognormal(self.mu, self.sigma, nsources)

    def pdf(self, lumi):
        r""" Gives the value of the PDF at lumi.

        Args:
            lumi: float or array-like, point where PDF is evaluated.

        Notes:
            PDF given by:

            $$ \frac{1}{x\sigma \sqrt{2\pi}} \times 
            \exp{-\frac{(\ln(x)-\mu)^2}{2\sigma^1}}$$
        """
        return lognorm.pdf(lumi, s=self.sigma, scale=np.exp(self.mu))

    def cdf(self, lumi):
        r""" Gives the value of the CDF at lumi.

        Args:
            lumi: float or array-like, point where CDF is evaluated.

        Notes:
            CDF given by:

            $$ \frac{1}{2} + \frac{1}{2} \times 
            \mathrm{erf}\left( \frac{(\ln(x)-\mu)^2}{\sqrt{2}\sigma}\right)$$
        """
        return lognorm.cdf(lumi, s=self.sigma, scale=np.exp(self.mu))


class PL_LuminosityFunction(LuminosityFunction):
    r""" Power-law distribution defined between x_min and x_max.

        PDF:
        $$\frac{1-\alpha}{x_{max}^{1-\alpha}-x_{min}^{1-\alpha}}\times
        x^{-\alpha}$$

        CDF:
        $$\frac{x^{1-\alpha}-x_{min}^{1-\alpha}}
        {x_{max}^{1-\alpha}-x_{min}^{1-\alpha}}$$

        inv.CDF:
        $$ P^{-1}(x) = (x_{min}^{\beta} + (x_{max}^{\beta}
         -x_{min}^{\beta})*x)^{1/\beta}$$

        Formulars from: http://up-rs-esp.github.io/bpl/
        and: Clauset, A., Shalizi, C. R., Newman, M. E. J.
        "Power-law Distributions in Empirical Data".
        SFI Working Paper: 2007-12-049 (2007)  arxiv:0706.1062
    """

    def __init__(self, mean_luminosity, index, width):
        """ Power-law distribution defined between x_min and x_max.
        x_min and x_max are calculated from mean and width of distribution.

        Note:
            f_mean ~= width*ln(10)*F_min,          index = -2
            f_mean ~= (index+1)/(index+2)*F_min    index < -2
        """
        self.mean = mean_luminosity
        self.index = index
        self.width = width

        if self.index == -2:
            self.Fmin = self.mean / (self.width*np.log(10))
        else:
            self.Fmin = (self.index+2.)/(self.index+1.)*self.mean

        self.Fmax = self.Fmin*10**self.width

    def sample_distribution(self, nsources=None, rng=None):
        """Samples from the Luminosity Function nsource times

        Args:
            number of sources
        """
        if rng is None:
            rng = np.random.RandomState()
        x = rng.uniform(0, 1, nsources)
        beta = (1.+self.index)
        return (self.Fmin**beta + (self.Fmax**beta -
                                   self.Fmin**beta)*x)**(1./beta)

    def pdf(self, lumi):
        r"""Gives the value of the PDF at lumi.

        Args:
            lumi: float or array-like, point where PDF is evaluated.

        Notes:
            PDF given by:
            $$\frac{1-\alpha}{x_{max}^{1-\alpha}-x_{min}^{1-\alpha}}
            \times x^{-\alpha}$$
        """

        norm = (1+self.index)/(self.Fmax**(1+self.index) -
                               self.Fmin**(1+self.index))
        pdf = norm*lumi**self.index
        pdf[np.logical_or(lumi < self.Fmin, lumi > self.Fmax)] = 0
        return pdf

    def cdf(self, lumi):
        r"""Gives the value of the CDF at lumi.

        Args:
            lumi: float or array-like, point where CDF is evaluated.

        Note:
            CDF given by:
            $$\frac{x^{1-\alpha}-x_{min}^{1-\alpha}}
            {x_{max}^{1-\alpha}-x_{min}^{1-\alpha}}$$
        """
        return (lumi**(1+self.index) - self.Fmin**(1+self.index)) / \
               ((self.Fmax**(1+self.index) - self.Fmin**(1+self.index)))


def get_LuminosityFunction(mean_luminosity, LF, **kwargs):
    """Returns a Luminosity Function based on an abbreviation.
    Known abbreviations are:
        - SC - Standard Candle
        - LG - Log-Normal
        - PL - Power-Law

    Args:
        - options: Namespace with LF and needed parameters to
                   construct luminosity function.
        - mean_luminosity: mean luminosity.

    Raises 'NotImplementedError' for unknown abbreviation.
    """

    if LF == "SC":
        return SC_LuminosityFunction(mean_luminosity)
    if LF == "LG":
        return LG_LuminosityFunction(mean_luminosity,
                                     kwargs["sigma"])
    if LF == "PL":
        return PL_LuminosityFunction(mean_luminosity,
                                     kwargs["index"],
                                     kwargs["sigma"])
    raise NotImplementedError("The luminosity function " +
                              LF + " is not implemented.")
