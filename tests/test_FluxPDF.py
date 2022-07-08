#!/usr/bin/python

"""
This script contains unittests for FluxPDF.py

Execute either as:
    python test_FluxPDF.py
or:
    python -m unittest test_FluxPDF
"""

import os
import unittest
import numpy as np
import firesong.FluxPDF as FluxPDF


class TestFluxPDFSimulation(unittest.TestCase):
    """ Tests flux pdf, which integrates the redshift
        and luminosity distributions and returns a 
        PDF of number vs. flux
    """

    @classmethod
    def setUpClass(cls):
        """ once before all tests """
        cls.sc_flux_pdf = FluxPDF.flux_pdf(None,
            filename=None,
            LF='SC',
            Evolution="MD2014SFR",
            verbose=False,
            logFMin=-20)
        cls.lg_flux_pdf = FluxPDF.flux_pdf(None,
            filename=None,
            LF='LG',
            Evolution="MD2014SFR",
            verbose=False,
            logFMin=-20,
            lg_width=0.01)

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

    def test_sc_flux_pdf(self):
        num_sources = np.sum(self.sc_flux_pdf[1])
        expectation = 13295.051834264095
        diff = (num_sources - expectation) / expectation
        self.assertTrue(diff < 0.1)

    def test_lg_flux_pdf(self):
        num_sources = np.sum(self.lg_flux_pdf[1])
        expectation = 13295.051834264095
        diff = (num_sources - expectation) / expectation
        self.assertTrue(np.abs(diff) < 0.1)

    def test_lg_sc_consistency(self):
        # Narrow LG should act similar to SC
        msk1 = (~np.isnan(np.array(self.sc_flux_pdf[1]))) & \
            (np.array(self.sc_flux_pdf[1]) != 0.)
        msk2 = (~np.isnan(np.array(self.lg_flux_pdf[1]))) & \
            (np.array(self.lg_flux_pdf[1]) != 0.)
        msk = msk1 * msk2
        ratio_arr = np.array(self.sc_flux_pdf[1])[msk] \
            / np.array(self.lg_flux_pdf[1])[msk]
        median_ratio = np.median(ratio_arr)
        diff = np.abs(median_ratio - 1.)
        self.assertTrue(diff < 0.05)


class TestTransientFluxPDFSimulation(unittest.TestCase):
    """ Tests flux pdf, which integrates the redshift
        and luminosity distributions and returns a 
        PDF of number vs. flux
    """

    @classmethod
    def setUpClass(cls):
        """ once before all tests """
        cls.sc_flux_pdf = FluxPDF.flux_pdf('./',
            filename='test_flux_pdf.out',
            LF='SC',
            Evolution="MD2014SFR",
            verbose=False,
            logFMin=-20,
            Transient=True,
            timescale=1000.,
            LumMax=1e59)
        cls.results = np.loadtxt('./test_flux_pdf.out').T

    @classmethod
    def tearDownClass(cls):
        """ once after all tests """
        os.remove('./test_flux_pdf.out')

    def setUp(self):
        "before each test"
        pass

    def tearDown(self):
        "after each test"
        pass

    ### tests start here ###

    def test_total_num(self):
        num_sources = np.sum(self.results[1])
        expectation = 3936.511006153128
        diff = (num_sources - expectation) / expectation
        self.assertTrue(diff < 0.01)


if __name__ == "__main__":
    unittest.main()
