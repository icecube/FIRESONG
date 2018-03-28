#!/usr/bin/python

"""
This script contains unittests for Evolution.py

Execute either as:
    python test_Evolution.py
or:
    python -m unittest test_Evolution
"""

import unittest
import Evolution


class TestEvolution(unittest.TestCase):
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
        self.evol1 = Evolution.get_evolution("NoEvolution")
        self.evol2 = Evolution.get_evolution("HB2006SFR")

    def tearDown(self):
        "after each test"
        pass

    ### tests start here ###

    def test_get_evolution(self):
        with self.assertRaises(NotImplementedError):
            Evolution.get_evolution("Test")

    def test_EvolutionBaseClass(self):
        pass

    def test_NoEvolution(self):
        self.assertEqual(self.evol1(1.), self.evol1(2.))

    def test_HopkinsBeacom2006StarFormationRate(self):
        pass

    def test_YukselEtAl2008StarFormationRate(self):
        pass

    def test_CandelsClash2015SNRate(self):
        pass

    def test_MadauDickinson2014CSFH(self):
        pass


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
                         3769.7731306195333)

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
                         9.2911752992097727e-11)


class TestTransientSourcePopulation(unittest.TestCase):
    """ Tests TransientSourcePopulation class """

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

    def test_0(self):
        pass


if __name__ == "__main__":
    unittest.main()
