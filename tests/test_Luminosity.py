#!/usr/bin/python

"""
This script contains unittests for Luminosity.py

Execute either as:
    python test_Luminosity.py
or:
    python -m unittest test_Luminosity
"""

import unittest
import firesong.Luminosity as Luminosity
import numpy as np


class TestLuminosity(unittest.TestCase):
    """ Tests Luminosity classes, with focus on standard
        candles and lognormal distributions
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
        pass

    def tearDown(self):
        "after each test"
        pass

    ### tests start here ###

    def test_get_LuminosityFunction(self):
        lf = Luminosity.get_LuminosityFunction(1e52, "SC")
        self.assertEqual(lf.mean, 1e52)
        with self.assertRaises(NotImplementedError):
            lf = Luminosity.get_LuminosityFunction(1e52, "Test")

    def test_LuminosityFunctionBaseClass_NotImplemented(self):
        lf = Luminosity.LuminosityFunction(1e52)
        with self.assertRaises(NotImplementedError):
            lf.pdf(10)
        with self.assertRaises(NotImplementedError):
            lf.sample_distribution(10)
        with self.assertRaises(NotImplementedError):
            lf.cdf(10)

    def test_SC_LuminosityFunction_CDF(self):
        lf = Luminosity.get_LuminosityFunction(1e50, "SC")
        lumis = np.logspace(49, 51, 3)
        self.assertEqual(lf.cdf(1e49), 0.0)
        self.assertEqual(lf.cdf(1e51), 1.0)
        self.assertEqual(lf.cdf(1e50), 0.5)
        cdf = lf.cdf(lumis)
        for l, c in zip(lumis, cdf):
            self.assertEqual(c, lf.cdf(l))

    def test_SC_LuminosityFunction_PDF(self):
        lf = Luminosity.get_LuminosityFunction(1e50, "SC")
        lumis = np.logspace(49, 51, 30)
        pdf = lf.pdf(lumis)
        self.assertEqual(np.sum(pdf == 0.0), len(pdf)-1)
        self.assertEqual(np.sum(pdf == 1.0), 1)

    def test_SC_LuminosityFunction_sample_distribution(self):
        lf = Luminosity.get_LuminosityFunction(1e50, "SC")
        lumi = lf.sample_distribution(nsources=10)
        self.assertEqual(lumi, 1e50)

    def test_LG_LuminosityFunction_sample_distribution(self):
        rng = np.random.RandomState(seed=1)
        lf = Luminosity.LG_LuminosityFunction(1e50, 1.)
        self.assertEqual(lf.sample_distribution(nsources=10, rng=rng)[0],
                         2.9720274559094399e+50)

    def test_LG_LuminosityFunction_PDF(self):
        lf = Luminosity.LG_LuminosityFunction(1e50, 1.)
        self.assertAlmostEqual(lf.pdf(1e49), 1.712868408774468e-50, places=10)

    def test_LG_LuminosityFunction_CDF(self):
        lf = Luminosity.LG_LuminosityFunction(1e50, 1.)
        self.assertEqual(lf.cdf(1e49), 0.56012752568024438)
        self.assertEqual(lf.cdf(1e60), 1.)

    def test_PL_LuminosityFunction(self):
        pass

if __name__ == "__main__":
    unittest.main()
