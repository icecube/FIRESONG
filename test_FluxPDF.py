#!/usr/bin/python

"""
This script contains unittests for FluxPDF.py

Execute either as:
    python test_FluxPDF.py
or:
    python -m unittest test_FluxPDF
"""

import unittest
import FluxPDF


class TestFiresongSimulation(unittest.TestCase):
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

    def test_flux_pdf(self):
        pass

if __name__ == "__main__":
    unittest.main()
