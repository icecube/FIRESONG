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
        """
        pass

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
    def __init__(self, mean_luminosity):
        self.mean = mean_luminosity

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

class BoundedPowerLaw:
    r"""Defines a power law with arbitrary scale in the form:

    $$ f(x) = A x^{-\alpha} $$

    Args:
        - A : float, scale factor
        - alpha : float, power law index
        - x0 : float, lower bound of the power law
        - x1 : float, upper bound of the power law
    """

    def __init__(self, A, alpha, x0, x1):
        self.A = A
        self.alpha = alpha
        self.beta = 1 - alpha
        self.x0 = x0
        self.x1 = x1
        self.integral = self._F(self.x1) - self._F(self.x0)

    def _f(self, x):
        """Unbounded power law"""
        return self.A * x ** (-self.alpha)

    def f(self, x):
        """Bounded power law"""
        return np.piecewise(
            x, [x < self.x0, x >= self.x1], [0, 0, lambda x: self._f(x)]
        )

    def _F(self, x):
        """Unbounded indefinite integral"""
        return self.A * x**self.beta / self.beta

    def F(self, x):
        """Definite integral [x0, x]"""
        return np.piecewise(
            x,
            [x < self.x0, x > self.x1],
            [0, self.integral, lambda x: self._F(x) - self._F(self.x0)],
        )

    def inv_F(self, F):
        """Inverse of function of F"""
        return (F * (self.beta / self.A) + self.x0**self.beta) ** (1 / self.beta)


class PL_LuminosityFunction(LuminosityFunction):
    r"""Power-law distribution defined between Lmin and Lmax.

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

    def __init__(self, Lmin, Lmax, alpha):
        """Power-law distribution defined between x_min and x_max.
        x_min and x_max are calculated from mean and width of distribution.

        Note:
            f_mean ~= width*ln(10)*F_min,          index = -2
            f_mean ~= (index+1)/(index+2)*F_min    index < -2
        """
        beta = 1 - alpha
        print(Lmax, Lmin, beta)
        D = Lmax**beta - Lmin**beta
        A = beta / D
        self.PL = BoundedPowerLaw(A=beta / D, alpha=alpha, x0=Lmin, x1=Lmax)

    def sample_distribution(self, nsources, rng=None):
        """Samples from the Luminosity Function nsource times

        Args:
            number of sources
        """
        if rng is None:
            rng = np.random.RandomState()
        F = rng.uniform(0, 1, nsources)
        return self.inv_cdf(F)

    def pdf(self, lumi):
        r"""Gives the value of the PDF at lumi.

        Args:
            lumi: float or array-like, point where PDF is evaluated.

        Notes:
            PDF given by:
            $$\frac{1-\alpha}{x_{max}^{1-\alpha}-x_{min}^{1-\alpha}}
            \times x^{-\alpha}$$
        """
        return self.PL.f(lumi)

    def cdf(self, lumi):
        r"""Gives the value of the CDF at lumi.

        Args:
            lumi: float or array-like, point where CDF is evaluated.

        Note:
            CDF given by:
            $$\frac{x^{1-\alpha}-x_{min}^{1-\alpha}}
            {x_{max}^{1-\alpha}-x_{min}^{1-\alpha}}$$
        """
        return self.PL.F(lumi)

    def inv_cdf(self, F):
        """
        Gives the inverse of the CDF at F
        """
        return self.PL.inv_F(F)

class BPL_LuminosityFunction(LuminosityFunction):
    """
    Broken power law consisting defined between Lmin and Lmax, consisting of two power laws joining at Lbreak.
    """

    def __init__(self, alpha1, alpha2, Lmin, Lbreak, Lmax):
        # auxiliary exponents
        beta1 = 1 - alpha1
        beta2 = 1 - alpha2
        delta = alpha2 - alpha1

        # normalisation factors
        D_ab1 = Lbreak**beta1 - Lmin**beta1
        D_ab2 = Lmax**beta2 - Lbreak**beta2

        # scale factors for power laws
        # calculated from conditions A1 + A2 = 1 and pdf1(Lbreak) = pdf2(Lbreak)
        A1 = 1 / (D_ab1 / beta1 + Lbreak**delta * D_ab2 / beta2)
        A2 = A1 * Lbreak**delta

        self.PL1 = BoundedPowerLaw(A=A1, alpha=alpha1, x0=Lmin, x1=Lbreak)
        self.PL2 = BoundedPowerLaw(A=A2, alpha=alpha2, x0=Lbreak, x1=Lmax)

    def sample_distribution(self, nsources, rng=None):
        if rng is None:
            rng = np.random.RandomState()
        F = rng.uniform(0, 1, nsources)
        return self.inv_cdf(F)

    def pdf(self, lumi):
        return self.PL1.f(lumi) + self.PL2.f(lumi)

    def cdf(self, lumi):
        return self.PL1.F(lumi) + self.PL2.F(lumi)

    def inv_cdf(self, F):
        return np.piecewise(
            F,
            [F <= self.PL1.integral, np.logical_and(F > self.PL1.integral, F < 1)],
            [
                lambda F: self.PL1.inv_F(F),
                lambda F: self.PL2.inv_F(F - self.PL1.integral),
            ],
        )

def get_LuminosityFunction(mean_luminosity, LF, **kwargs):
    """Returns a Luminosity Function based on a configuration string.
    Supported LFs are:
        - SC - Standard Candle
        - LG - Log-Normal
        - PL - Power-Law
        - BPL - Broken-Power-Law

    Args:
        - LF_conf: string, configuration string for the Luminosity Function. 
          The string contains the LF type.

        - kwargs: The specific arguments for each luminosity function
          Log-Normal: lg_width
          Power-Law: pl_lmin, pl_lmax, pl_alpha
          Broken Power-Law: bpl_lmin, bpl_lbreak, bpl_lmax, bpl_alpha1, bpl_alpha2

    Raises 'NotImplementedError' for unknown abbreviation.
    """

    # TODO: this would better go into a try-except block
    LF_type = LF
    if LF_type == "SC":
        return SC_LuminosityFunction(mean_luminosity=mean_luminosity, )
    if LF_type == "LG":
        width = kwargs.get('lg_width', 1.0)
        return LG_LuminosityFunction(mean_luminosity=mean_luminosity, 
                                     width=width)
    if LF_type == "PL":
        lmin = kwargs.get('pl_lmin')
        lmax = kwargs.get('pl_lmax')
        alpha = kwargs.get('pl_alpha')
        return PL_LuminosityFunction(Lmin=lmin, Lmax=lmax, alpha=alpha)
    if LF_type == "BPL":
        lmin = kwargs.get('bpl_lmin')
        lmax = kwargs.get('bpl_lmax')
        lbreak = kwargs.get('bpl_lbreak')
        alpha1 = kwargs.get('bpl_alpha1')
        alpha2 = kwargs.get('bpl_alpha2')
        return BPL_LuminosityFunction(Lmin=lmin, Lbreak=lbreak, Lmax=lmax, 
                                      alpha1=alpha1, alpha2=alpha2)

    raise NotImplementedError("The luminosity function " + LF + " is not implemented.")
