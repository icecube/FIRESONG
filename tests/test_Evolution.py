#!/usr/bin/python

"""
This script contains unittests for Evolution.py

Execute either as:
    python test_Evolution.py
or:
    python -m unittest test_Evolution
"""

import unittest
import firesong.Evolution as Evolution


class TestEvolution(unittest.TestCase):
    """ Tests get_evolution function, Evolution class and
        all implemented Evolutions.
    """

    @classmethod
    def setUpClass(cls):
        """ once before all tests """
        cls.evol1 = Evolution.get_evolution("NoEvolution")
        cls.evol2 = Evolution.get_evolution("HB2006SFR")
        cls.evol3 = Evolution.get_evolution("MD2014SFR")
        cls.evol4 = Evolution.get_evolution("YMKBH2008SFR")
        cls.evol5 = Evolution.get_evolution("CC2015SNR")

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

    def test_get_evolution(self):
        with self.assertRaises(NotImplementedError):
            Evolution.get_evolution("Test")

    def test_EvolutionBaseClass(self):
        evol = Evolution.Evolution()
        with self.assertRaises(NotImplementedError):
            evol(1.)

    def test_NoEvolution(self):
        self.assertEqual(self.evol1(1.), self.evol1(2.))
        self.assertEqual(str(self.evol1), 
            "No Evolution")

    def test_HopkinsBeacom2006StarFormationRate(self):
        self.assertAlmostEqual(self.evol2(1.),
            0.14702066600558594,
            3)
        self.assertEqual(str(self.evol2), 
            "Hopkins and Beacom (2006)")

    def test_YukselEtAl2008StarFormationRate(self):
        self.assertAlmostEqual(self.evol4(1.),
            0.19698310613518824,
            3)
        self.assertEqual(str(self.evol4), 
            "Yuksel et al. (2008)")

    def test_CandelsClash2015SNRate(self):
        self.assertAlmostEqual(self.evol5(1.),
            0.0707688850689447,
            3)
        self.assertEqual(str(self.evol5), 
            "Strolger et al. (2015)")

    def test_MadauDickinson2014CSFH(self):
        self.assertAlmostEqual(self.evol3(1.),
            0.08665290222624604,
            3)
        self.assertEqual(str(self.evol3), 
            "Madau and Dickinson (2014)")


class TestSourcePopulation(unittest.TestCase):
    """ Tests SourcePopulation class """

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
        evol = Evolution.get_evolution("NoEvolution")
        self.pop = Evolution.SourcePopulation(Evolution.cosmology,
                                              evol)

    def tearDown(self):
        "after each test"
        pass

    ### tests start here ###

    def test_Nsources(self):
        self.assertEqual(self.pop.Nsources(1e-9, zmax=10.),
                         3745.1796382799967)

    def test_Flux_Lumi(self):
        L = self.pop.Lumi2Flux(1e50, 2.0, 1e3, 1e7, z=1)
        self.assertAlmostEqual(1e50 / self.pop.Flux2Lumi(L,
                                                         2.0,
                                                         1e3,
                                                         1e7,
                                                         z=1),
                               1)

    def test_Diffuse2SC(self):
        self.assertEqual(self.pop.StandardCandleSources(1e-8,
                                                        1e-9,
                                                        zmax=10.,
                                                        index=2.0),
                         9.244753729666578e-11)


class TestTransientSourcePopulation(unittest.TestCase):
    """ Tests TransientSourcePopulation class """

    @classmethod
    def setUpClass(cls):
        """ once before all tests """
        evol = Evolution.get_evolution("NoEvolution")
        cls.transient_pop = Evolution.TransientSourcePopulation(
            Evolution.cosmology,
            evol,
            1000.)
        cls.pop = Evolution.SourcePopulation(
            Evolution.cosmology,
            evol)

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

    def test_redshift_dist(self):
        self.assertEqual(
            self.transient_pop.RedshiftDistribution(2.) * (1.+2.),
            self.pop.RedshiftDistribution(2.)
        )


if __name__ == "__main__":
    unittest.main()
