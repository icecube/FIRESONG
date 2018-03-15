#!/usr/bin/python

"""
This script contains unittests for Luminosity.py

Execute either as:
    python test_Luminosity.py
or:
    python -m unittest test_Luminosity
"""

import unittest
import sampling


class TestLuminosity(unittest.TestCase):
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
        pass

    def tearDown(self):
        "after each test"
        pass

    ### tests start here ###

    def test_get_LuminosityFunction(self):
        pass

    def test_LuminosityFunctionBaseClass(self):
        pass

    def test_SC_LuminosityFunction(self):
        pass

    def test_LG_LuminosityFunction(self):
        pass

    def test_PL_LuminosityFunction(self):
        pass

if __name__ == "__main__":
    unittest.main()
