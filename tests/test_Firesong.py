#!/usr/bin/python

"""
This script contains unittests for Firesong.py

Execute either as:
    python test_Firesong.py
or:
    python -m unittest test_Firesong
"""

import os
import unittest
import firesong.Firesong as Firesong
import firesong.Evolution as Evolution
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
            np.log10(6.099042878096284e+52),
            places=3)

    def test_firesong_number(self):
        self.assertEqual(len(self.firesong_uni['sources']['z']),
            18370)

    def test_median_redshift(self):
        self.assertAlmostEqual(np.median(self.firesong_uni['sources']['z']),
            2.81456604,
            places=1)

    def test_total_flux(self):
        # Compares both the quoted sum to the direct
        # sum of sources, compares these both to input flux
        self.assertAlmostEqual(np.log10(self.firesong_uni['total_flux']/1.44e-8),
            0.0,
            places=1)
        self.assertAlmostEqual(np.log10(np.sum(self.firesong_uni['sources']['flux'])),
            np.log10(self.firesong_uni['total_flux']*4.*np.pi),
            places=3)

class TestTransientFiresongSimulation(unittest.TestCase):
    """ Tests Firesong.firesong_simulation.
    Also checks to make sure Firesong and Evolution
    are working properly together, also ensures that 
    the file writing is working properly
    """

    @classmethod
    def setUpClass(cls):
        """ once before all tests """
        cls.evol_name = 'MD2014SFR'
        firesong_uni = Firesong.firesong_simulation('./',
            filename='test_results_tmp.out',
            Evolution=cls.evol_name,
            verbose=True,
            seed=1234,
            Transient=True)
        cls.res = np.loadtxt('./test_results_tmp.out')
        with open('./test_results_tmp.out', 'r') as f:
            header_info = f.readlines()
        header_info = header_info[:11] + header_info[-1:]
        cls.header_info = header_info
        cls.evol = Evolution.get_evolution(cls.evol_name)

    @classmethod
    def tearDownClass(cls):
        """ once after all tests """
        os.remove('./test_results_tmp.out')

    def setUp(self):
        "before each test"
        pass

    def tearDown(self):
        "after each test"
        pass

    def test_summed_fluxes(self):
        # Check that fluxes add up properly
        self.assertAlmostEqual(
            np.sum(self.res.T[3] * (1.+self.res.T[2])) \
                / (86400. * 365.) / (4.*np.pi) * 1000.,
            float(self.header_info[-1].split()[-1]),
            5
        )

    def test_sc_luminosity(self):
        # Check that luminosity is correct
        self.assertEqual(
            float(self.header_info[7].split()[-2]),
            4.9456e+57
        )

    def test_total_simulated_flux(self):
        # Check that simulated flux is roughly the input flux
        sim_flux = float(self.header_info[-1].split()[-1])
        input_flux = float(self.header_info[4].split()[4])
        self.assertAlmostEqual(np.log10(sim_flux),
            np.log10(input_flux),
            0)


if __name__ == "__main__":
    unittest.main()
