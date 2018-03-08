import numpy as np
from scipy.interpolate import UnivariateSpline


class InverseCDF(object):
    """ Class that constructs the inverse CDF out of a given pdf(x) array.
    The object is callable and returns the inverse of the CDF.
    The object has a method 'sample' that gives random numbers according
    to the pdf.
    """

    def __init__(self, x, pdf, seed=None):
        """ Constructs a InverseCDF object

        Parameters:
            - x    array like, x points at which probability density
                   function is given
            - pdf  array like, probability density function values for x points
            - (seed) seed of random number generator to sample distribution
        """

        self.rng = np.random.RandomState(seed)
        cdf = np.cumsum(pdf)/np.sum(pdf)
        mask = np.diff(cdf) > 0
        self.invCDF = UnivariateSpline(cdf[1:][mask], x[1:][mask], k=1, s=0)

    def __call__(self, x):
        """ Returns the inverse CDF at point(s) x """
        return self.invCDF(x)

    def sample(self, N=None):
        """ Returns random number variable(s) distributed according to
        PDF used at construction. At construction time the seed can be set.

        Parameters:
            - (N) number of random numbers, default=None

        If N is not given (or None) a single value is returned.
        """

        try:
            temp = self.invCDF(self.rng.uniform(0.0, 1.0, size=N), ext=0)
        except:
            temp = self.invCDF(self.rng.uniform(0.0, 1.0, size=N))
        return temp
