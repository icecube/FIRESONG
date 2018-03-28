#!/usr/bin/python

"""
This script contains unittests for input_output.py

Execute either as:
    python test_input_output.py
or:
    python -m unittest test_input_output
"""

import unittest
import input_output


class TestOutputWriter(unittest.TestCase):
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

    def test_output_writer(self):
        pass

if __name__ == "__main__":
    unittest.main()
