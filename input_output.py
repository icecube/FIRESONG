#!/usr/bin/python
import re
import os
import gzip
import numpy as np


class output_writer(object):
    def __init__(self, outputdir, filename, zNEAR=0):
        self.output = self.open_file(outputdir, filename)
        self.z_near = zNEAR
        if (self.z_near > 0):
            self.near_output = self.open_file(outputdir, "Near_" + filename)

    def open_file(self, outputdir, filename):
        if re.search('.gz$', filename):
            output = gzip.open(outputdir+str(filename), 'wb')
        else:
            output = open(outputdir+str(filename), "w")
        return output

    def write_header(self, LF, Transient, timescale, fluxnorm,
                     delta_gamma, luminosity):
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
        temp += "# Dec(deg) Redshift A\n"
        self.output.write(temp.format(**locals()))

    def write(self, declin, redshift, flux):
        for d, z, f in zip(np.atleast_1d(declin),
                           np.atleast_1d(redshift),
                           np.atleast_1d(flux)):
            self.output.write('{:.4f} {:.4f} {:.6e}\n'.format(d, z, f))
            # CHECK why different order
            if z < self.z_near:
                self.near_output.write('{:.4e} {:.4f} {:.4f}\n'.format(f, d, z))

    def finish(self, tot_flux=None):
        """give tot_flux per sr """
        if tot_flux is not None:
            self.output.write("# E^2 dNdE = {tot_flux}\n".format(**locals()))
        self.output.close()
        if (self.z_near > 0):
            self.near_output.close()


class output_writer_CDF(output_writer):
    def write(self, z, flux, nuCDF):
        self.output.write('{:.4f} {:.6e} {:.6e}\n'.format(float(z), flux, nuCDF))

class output_writer_PDF(output_writer):
    def write(self, x, y):
        np.savetxt(self.output, np.array([x, y]).T)

def get_outputdir():
    """ Check that the Firesong environmental variable is set
    This is needed for output and to read exposures, effective areas, etc

    Returns path to result folder
    """
    try:
        firesongdir = os.environ['FIRESONG']
    except:
        print "Enviromental variable FIRESONG not set"
        print "set to ./"
        firesongdir = "./"
    return os.path.join(firesongdir, "Results/")


def print_config(LF, Transient, timescale, Evolution, density,
                 N_sample, luminosity_default, fluxnorm, delta_gamma,
                 zmax, luminosity, mode="", **kwargs):
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
    str += "Redshift range: 0 - {zmax}\n"
    str += "Standard Candle Luminosity: {luminosity:.4e} erg/yr\n"
    str += "##### FIRESONG initialization done #####\n"
    print(str.format(**locals()))
