#!/usr/bin/python

"""
This script contains unittests for sampling.py

Execute either as:
    python test_sampling.py
or:
    python -m unittest test_sampling
"""

import unittest
import sampling
import numpy as np

from scipy.stats import truncnorm

class TestInverseCDF(unittest.TestCase):
    """ Tests get_evolution function, Evolution class and
        all implemented Evolutions.
    """

    @classmethod
    def setUpClass(cls):
        """ once before all tests """
        pass

    @classmethod
    def tearDownClass(cls):
        """ once after all tests """
        pass

    def setUp(self):
        "before each test"
        mean, sigma = 50, 5
        x = np.linspace(0, 100., 10001).astype(np.float32)
        
        bounds = [(x.min() - mean) / sigma,
                  (x.max() - mean) / sigma]
        self.truncnorm = truncnorm(bounds[0], bounds[1], loc=50, scale=5)

        pdf = self.truncnorm.pdf(x)
        self.invCDF = sampling.InverseCDF(x, pdf, seed=1)


    def tearDown(self):
        "after each test"
        pass

    ### tests start here ###

    def test_InverseCDF_evaluation(self):
        self.assertAlmostEqual(self.invCDF(1.0)[0], self.truncnorm.ppf(1.0), 2)
        self.assertAlmostEqual(self.invCDF(0.0)[0], self.truncnorm.ppf(0.0), 2)
        self.assertAlmostEqual(self.invCDF(0.5)[0], self.truncnorm.ppf(0.5), 2)

    def test_InverseCDF_evaluation_out_of_bound(self):
        with self.assertRaises(ValueError):
            self.invCDF(1.1)
        with self.assertRaises(ValueError):
            self.invCDF(-0.1)

    def test_random_default(self):
        self.assertEqual(self.invCDF.sample(), 52.914028811506654)

    def test_random_given_N_1(self):
        self.assertEqual(self.invCDF.sample(N=1), 52.914028811506654)

    def test_random_given_N_5(self):
        x = self.invCDF.sample(N=5)
        self.assertEqual(len(x), 5)
        self.assertEqual(x[0], 43.362544905234685)
        self.assertNotEqual(x[0], x[1])

    def test_random_wrong_parameter(self):
        with self.assertRaises(TypeError):
            x = self.invCDF.sample("Test")
        with self.assertRaises(TypeError):
            x = self.invCDF.sample(1.5)
        #self.assertEqual(len(x), 1)

if __name__ == "__main__":
    unittest.main()
