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
            evol = Evolution.get_evolution("Test")

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
        pass

    def tearDown(self):
        "after each test"
        pass

    ### tests start here ###

    def test_0(self):
        pass


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
