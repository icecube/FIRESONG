#!/usr/bin/python

"""Classed to handle the writing of simulation results to files"""

import re
import os
import gzip
import numpy as np


class output_writer(object):
    """
    Class for handling the output of FIRESONG simulations and writing
    to files

    Args:
        outputdir (str): Path for saving results
        filename (str): Desired output filename
        zNEAR (float, optional, default=0): Only consider nearby sources
            out to redshift=zNEAR if zNEAR is not 0

    Attributes:
        output (file): File being written
        z_near (float): Only consider nearby sources
            out to redshift=zNEAR if zNEAR is not 0
    """
    def __init__(self, outputdir, filename, zNEAR=0):
        """
        """
        self.output = self.open_file(outputdir, filename)
        self.z_near = zNEAR
        if (self.z_near > 0):
            self.near_output = self.open_file(outputdir, "Near_" + filename)

    def open_file(self, outputdir, filename):
        """
        Open file and prepare for writing

        Args:
            outputdir (str): Path for saving results
            filename (str): Desired output filename
        
        Returns:
            file: File being written to
        """
        if re.search('.gz$', filename):
            output = gzip.open(os.path.join(outputdir, str(filename)), 'wb')
        else:
            output = open(os.path.join(outputdir, str(filename)), "w")
        return output

    def write_header(self, LF, Transient, timescale, fluxnorm,
                     delta_gamma, luminosity, Gammaflux, interaction):
        """
        Record relevant simulation parameters in the top of file

        Args:
            LF (str): Luminosity function
            Transient (bool): Transient or steady sources
            timescale (float): Timescale for transient sources
            fluxnorm (float): Normalization on total astrophysical 
                diffuse flux
            delta_gamma (float): Spectral index - 2.
            luminosity (float): Luminosity of sources, 0 if just saturating
                diffuse flux
            Gammaflux (bool): include gamma-ray flux
            interaction (string): type of interaction used to calculate gamma-ray flux
        """
        temp = "# FIRESONG Output description\n"
        if LF == "SC":
            temp += "# Standard candle sources\n"
        if LF == "LG":
            temp += "# Lognormal distributed sources\n"
        if LF == "PL":
            temp += "# PowerLaw distributed sources\n"
        if Transient:
            temp += "# Transient Sources, time scale = {timescale}s \n"
        temp += "# Desired neutrino diffuse flux:\n"
        temp += "#      E^2 dN_diffuse/dE = {fluxnorm} (E/100 TeV)^({delta_gamma}) [GeV/cm^2.s.sr]\n"
        temp += "# Neutrino point source fluxes listed below are of \'A\' where the flux is:\n"
        temp += "#      E^2 dN_PS/dE = A * (E/100 TeV)^({delta_gamma}) [GeV/cm^2.s]\n"
        temp += "# Standard Candle Luminosity: {luminosity:.4e} erg/yr \n"
        temp += "# Note that using 4 years, IceCube sensitivity in the northern hemisphere\n"
        temp += "# is approximately 10^-9 in the units used for A\n"
        if Gammaflux:
            temp += "# Corresponding gamma-ray fluxes listed below are of \'B\' where the flux is:\n"
            temp += "#      dN_PS/dE = B * (E_gamma/1 TeV)^({delta_gamma}-2) [1/TeV.cm^2.s]\n"
            temp += "# Gamma-ray fluxes calculated assuming {interaction} interactions\n"
            temp += "# Dec(deg) RA(deg) Redshift A(GeV/cm^2.s) B(1/TeV.cm^2.s)\n"
        else:
            temp += "# Dec(deg) RA(deg) Redshift A\n"
        self.output.write(temp.format(**locals()))

    def write(self, declin, ra, redshift, flux):
        """
        Write the sources to the output file

        Args:
            declin (array or float): source declination(s)
            ra (array or float): source ra(s)
            redshift (array or float): source redshift(s)
            flux (array or float): source flux(es)
        """
        for d, r, z, f in zip(np.atleast_1d(declin),
                           np.atleast_1d(ra),
                           np.atleast_1d(redshift),
                           np.atleast_1d(flux)):
            self.output.write('{:.4f} {:.4f} {:.4f} {:.6e}\n'.format(d, r, z, f))
            # CHECK why different order
            if z < self.z_near:
                self.near_output.write('{:.4e} {:.4f} {:.4f} {:.4f}\n'.format(f, d, r, z))

    def write_gamma(self, declin, ra, redshift, flux, gammaflux):
        """
        Write the sources to the output file, including the gamma-ray flux

        Args:
            declin (array or float): source declination(s)
            ra (array or float): source ra(s)
            redshift (array or float): source redshift(s)
            flux (array or float): source flux(es)
            gammaflux (array or float): equivalent gamma-ray flux(es)
        """
        for d, r, z, f, g in zip(np.atleast_1d(declin),
                              np.atleast_1d(ra),
                              np.atleast_1d(redshift),
                              np.atleast_1d(flux),
                              np.atleast_1d(gammaflux)):
            self.output.write('{:.4f} {:.4f} {:.4f} {:.6e} {:.6e}\n'.format(d, r, z, f, g))
            if z < self.z_near:
                self.near_output.write('{:.4e} {:.4f} {:.4f} {:.4f} {.4f}\n'.format(f, d, r, z, g))

    def finish(self, tot_flux=None):
        """
        Calculate total flux per steradian and close file
        """
        if tot_flux is not None:
            self.output.write("# E^2 dNdE = {tot_flux}\n".format(**locals()))
        self.output.close()
        if (self.z_near > 0):
            self.near_output.close()


class output_writer_CDF(output_writer):
    """
    Write the results of calculating the neutrino CDF instead of 
    the general firesong_simulation
    """
    def write(self, z, flux, nuCDF):
        self.output.write('{:.4f} {:.6e} {:.6e}\n'.format(float(z), flux, nuCDF))


class output_writer_PDF(output_writer):
    """
    Write the results of calculating the neutrino flux PDF instead of 
    the general firesong_simulation
    """
    def write(self, x, y):
        np.savetxt(self.output, np.array([x, y]).T)


def get_outputdir():
    """ 
    Check that the Firesong environmental variable is set
    This is needed for output and to read exposures, effective areas, etc

    Returns:
        str: path to result folder
    """
    try:
        firesongdir = os.environ['FIRESONG']
    except:
        print("Enviromental variable FIRESONG not set")
        print("set to ./")
        firesongdir = "./"


    return firesongdir


def print_config(LF, Transient, timescale, Evolution, density,
                 N_sample, luminosity_default, fluxnorm, delta_gamma,
                 zmax, luminosity, Gammaflux, interaction, mode="", **kwargs):
    """
    Prints the configuration to the screen.
    """
    str = "##############################################################################\n"
    str += "##### FIRESONG initializing {mode}#####\n"
    if LF == "SC":
        str += "Standard candle sources+\n"
    if LF == "LG":
        str += "Lognormal distributed sources\n"
    if LF == "PL":
        str += "PowerLaw distributed sources\n"
    if Transient:
        str += "Transient Sources, time scale = {timescale}s\n"
    str += "Source evolution assumed: {Evolution}\n"
    str += "Local density of neutrino sources: {density}/Mpc^3\n"
    str += "Total number of neutrinos sources in the Universe: {N_sample}\n"
    if luminosity_default == 0.0:
        str += "Desired neutrino diffuse flux: E^2 dN/dE = {fluxnorm} (E/100 TeV)^({delta_gamma}) GeV/cm^2.s.sr\n"
    if Gammaflux:
        str += "Calculate gamma-ray flux, assuming {interaction} interactions\n"
    str += "Redshift range: 0 - {zmax}\n"
    str += "Standard Candle Luminosity: {luminosity:.4e} erg/yr\n"
    str += "##### FIRESONG initialization done #####\n"
    print((str.format(**locals())))

def print_config_LEGEND(L_Evolution, lmin, lmax, N_sample):
    """
    Prints the configuration to the screen for the luminosity evolution model.
    """
    str = "##############################################################################\n"
    str += "##### LEGEND initializing        #####\n"
    str += "Luminosity Evolution model = {L_Evolution}\n"
    str += "Source Luminosity Limit = {lmin} - {lmax} erg/s \n"
    str += "Total number of sources in the universe: {N_sample} \n"
    str += "##### LEGEND initialization done #####"
    print((str.format(**locals())))
