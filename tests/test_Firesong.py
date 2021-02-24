#!/usr/bin/python

"""
This script contains unittests for Firesong.py

Execute either as:
    python test_Firesong.py
or:
    python -m unittest test_Firesong
"""

import unittest
import Firesong, Evolution
import numpy as np

class TestFiresongSimulation(unittest.TestCase):
    """ Tests Firesong.firesong_simulation.
    Also checks to make sure Firesong and Evolution
    are working properly together
    """

    @classmethod
    def setUpClass(cls):
        """ once before all tests """
        cls.evol_name = 'HB2006SFR'
        cls.firesong_uni = Firesong.firesong_simulation(None,
        filename=None,
        Evolution=cls.evol_name,
        verbose=False,
        seed=1234)
        cls.evol = Evolution.get_evolution(cls.evol_name)

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

    def test_firesong_luminosity(self):
        self.assertAlmostEqual(np.log10(self.firesong_uni['header']['luminosity']),
            np.log10(3.6471202781753894e+52),
            places=3)

    def test_firesong_number(self):
        self.assertEqual(len(self.firesong_uni['sources']['z']),
            18460)

    def test_median_redshift(self):
        self.assertAlmostEqual(np.median(self.firesong_uni['sources']['z']),
            2.81456604,
            places=1)

    def test_total_flux(self):
        # Compares both the quoted sum to the direct
        # sum of sources, compares these both to input flux
        self.assertAlmostEqual(np.log10(self.firesong_uni['total_flux']/9e-9),
            0.0,
            places=2)
        self.assertAlmostEqual(np.log10(np.sum(self.firesong_uni['sources']['flux'])),
            np.log10(self.firesong_uni['total_flux']*4.*np.pi),
            places=3)
        

if __name__ == "__main__":
    unittest.main()
