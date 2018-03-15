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
        pass

    def tearDown(self):
        "after each test"
        pass

    ### tests start here ###

    def test_InverseCDF(self):
        pass

if __name__ == "__main__":
    unittest.main()
