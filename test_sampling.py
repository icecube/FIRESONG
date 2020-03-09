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
        x = np.linspace(0, 100., 101)
        pdf = np.exp(-0.5*((x-50.)/5.)**2)
        self.invCDF = sampling.InverseCDF(x, pdf, seed=1)

    def tearDown(self):
        "after each test"
        pass

    ### tests start here ###

    def test_InverseCDF_evaluation(self):
        self.assertTrue(self.invCDF(1.0) <= 100)
        self.assertTrue(self.invCDF(0.0) >= 0)
        self.assertAlmostEqual(self.invCDF(0.5), 50)

    def test_InverseCDF_evaluation_out_of_bound(self):
        with self.assertRaises(ValueError):
            self.invCDF(1.1)
        with self.assertRaises(ValueError):
            self.invCDF(-0.1)

    def test_random_default(self):
        self.assertEqual(self.invCDF.sample(), 48.949116836587841)

    def test_random_given_N_1(self):
        self.assertEqual(self.invCDF.sample(N=1), 48.949116836587841)

    def test_random_given_N_5(self):
        x = self.invCDF.sample(N=5)
        self.assertEqual(len(x), 5)
        self.assertEqual(x[0], 48.949116836587841)
        self.assertNotEqual(x[0], x[1])

    def test_random_wrong_parameter(self):
        with self.assertRaises(TypeError):
            x = self.invCDF.sample("Test")
        with self.assertRaises(TypeError):
            x = self.invCDF.sample(1.5)
        #self.assertEqual(len(x), 1)

if __name__ == "__main__":
    unittest.main()
