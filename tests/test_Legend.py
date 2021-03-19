#!/usr/bin/python

"""
This script contains unittests for Legend.py

Execute either as:
    python test_Legend.py
or:
    python -m unittest test_Legend
"""

import os
import unittest
import firesong.Legend as Legend
import firesong.Evolution as Evolution
import numpy as np

class TestLegendSimulation(unittest.TestCase):
    """ Tests Legend.legend_simulation.
    Also checks to make sure Legend and Evolution
    are working properly together
    """

    @classmethod
    def setUpClass(cls):
        """ once before all tests """
        cls.evol_name = 'HA2014BL'
        cls.legend_uni = Legend.legend_simulation(None,
        filename=None,
        L_Evolution=cls.evol_name,
        verbose=False,
        lmin=40,
        lmax=48,
        seed=1234)

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
    def test_legend_number(self):
        self.assertEqual(len(self.legend_uni['z']),
            3209)

    def test_median_redshift(self):
        self.assertAlmostEqual(np.median(self.legend_uni['z']),
            2.396673768245595,
            places=1)

    def test_total_flux(self):
        self.assertTrue(np.abs(np.log10(np.sum(self.legend_uni['flux'])/7.117604312572643e-06)) < 0.1)

if __name__ == "__main__":
    unittest.main()
