"""Implementation of inverse CDF sampling"""

import numpy as np
from scipy.interpolate import interp1d


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

        # Pad the arrays to ensure the cdf starts/ends at the right points
        cdf = np.concatenate([[0,], cdf, [1,]])
        x = np.concatenate([[x.min()-np.finfo(x.dtype).eps,],
                            x,
                            [x.max()+np.finfo(x.dtype).eps,]])

        # Also prepend a 0 here to ensure the shape stays the same
        mask = np.diff(cdf, prepend=0) >= 0
        self.invCDF = interp1d(cdf[mask], x[mask], kind='linear',
                               bounds_error=False, fill_value=(x.min(), x.max()))

    def __call__(self, x):
        x = np.atleast_1d(x)
        """ Returns the inverse CDF at point(s) x """
        if any(x < 0) or any(x > 1):
            raise ValueError("x out of bounds. x has to be in range 0..1")
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
