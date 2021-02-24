#!/usr/bin/python

"""
This script contains unittests for SingleNeutrinoCDF.py

Execute either as:
    python test_SingleNeutrinoCDF.py
or:
    python -m unittest test_SingleNeutrinoCDF
"""

import unittest
import SingleNeutrinoCDF


class TestSingleNeutrinoCDF(unittest.TestCase):
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

    def test_calc_SingleNeutrinoCDF(self):
        pass

if __name__ == "__main__":
    unittest.main()
