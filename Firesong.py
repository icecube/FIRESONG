#!/usr/bin/python
# Authors: Chris Tung
#          Igncio Taboada
#

# General imports
from __future__ import division
import os
import gzip
import re
import argparse
# Numpy / Scipy
import numpy as np
# Firesong code
from Evolution import get_evolution, SourcePopulation
from Evolution import TransientSourcePopulation, cosmology, Simulation
from Luminosity import get_LuminosityFunction


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
            self.output.write('{:.4f} {:.4f} {:.4e}\n'.format(d, z, f))
            # CHECK why different order
            if z < self.z_near:
                self.near_output.write('{:.4e} {:.4f} {:.4f}\n'.format(f, d, z))

    def finish(self, tot_flux):
        """give tot_flux per sr """
        self.output.write("# E^2 dNdE = {tot_flux}\n".format(**locals()))
        self.output.close()
        if (self.z_near > 0):
            self.near_output.close()


def print_str(LF, Transient, timescale, Evolution, density,
              N_sample, luminosity_default, fluxnorm, delta_gamma,
              zmax, candleflux, luminosity):
    str = "##############################################################################\n"
    str += "##### FIRESONG initializing #####\n"
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
    str += "CandleFlux at z=1: {candleflux:.4e} GeV/cm^2.s\n"
    str += "Standard Candle Luminosity: {luminosity:.4e} erg/yr\n"
    str += "##### FIRESONG initialization done #####\n"
    print(str.format(**locals()))


def firesong_simulation(options, outputdir):

    if options.Transient is True:
        population = TransientSourcePopulation(cosmology,
                                               get_evolution(options.Evolution),
                                               timescale=options.timescale)
    else:
        population = SourcePopulation(cosmology,
                                      get_evolution(options.Evolution))

    N_sample = int(population.Nsources(options.density, options.zmax))

    z0 = 1.
    if options.luminosity == 0.0:
        ## If luminosity not specified calculate candleflux from diffuse flux
        candleflux = population.StandardCandleSources(options.fluxnorm,
                                                      options.density,
                                                      options.zmax,
                                                      options.index,
                                                      z0=z0)
        luminosity = population.Flux2Lumi(candleflux,
                                          options.index,
                                          emin=1.e4,
                                          emax=1.e7,
                                          z=z0,
                                          E0=1e5)
    else:
        ## If luminosity of the sources is specified, calculate candleflux
        candleflux = population.Lumi2Flux(options.luminosity,
                                          options.index,
                                          emin=1.e4,
                                          emax=1.e7,
                                          z=z0,
                                          E0=1.e5)
        luminosity = options.luminosity

    delta_gamma = 2-options.index
    print_str(options.LF, options.Transient, options.timescale,
              options.Evolution, options.density, N_sample,
              options.luminosity, options.fluxnorm, delta_gamma,
              options.zmax, candleflux, luminosity)

    ##################################################
    #        Simulation starts here
    ##################################################

    simulation = Simulation(population,
                            get_LuminosityFunction(options, candleflux),
                            index=options.index,
                            zmax=options.zmax,
                            z0=z0)

    out = output_writer(outputdir, options.filename)
    out.write_header(options.LF, options.Transient, options.timescale,
                     options.fluxnorm, delta_gamma, options.luminosity)

    # CHECK loop can be more efficient
    TotalFlux = 0
    for i in range(0, N_sample):
        # IMPORTANT notice, in the following "flux" means fluence in
        # Transient mode, but flux in steady source mode,
        # until TotalFlux(TotalFluence) is calculated
        flux, z, declin = simulation.sample_flux()
        TotalFlux = TotalFlux + flux

        # For transient sources, the flux measured on Earth will be
        # red-shifted-fluence/{(1+z)*burst duration}
        if options.Transient:
            flux = population.fluence2flux(flux, z)

        out.write(declin, z, flux)
        if i % 100000 == 0 and i > 0:
            print "Generated ", i, " neutrino sources"

    # For transient source, we calculate the total fluence from all sources,
    # then obtain the diffuse flux by doing a time average over a year
    if options.Transient:
        TotalFlux = TotalFlux / population.yr2sec
    TotalFlux /= 4*np.pi  # give in per sr

    out.finish(TotalFlux)
    print "Actual diffuse flux simulated :  E^2 dNdE = {TotalFlux} (E/100 TeV)^({delta_gamma}) [GeV/cm^2.s.sr]".format(**locals())


def firesong_simulation_by_hand(outputdir,
                                filename='Firesong.out',
                                density=1e-9,
                                Evolution="HB2006SFR",
                                Transient=False,
                                timescale=1000.,
                                zmax=10.,
                                fluxnorm=0.9e-8,
                                index=2.13,
                                LF="SC",
                                sigma=1.0,
                                zNEAR=-1,
                                luminosity=0.0):

    options = argparse.Namespace(filename=filename,
                                 density=density,
                                 Evolution=Evolution,
                                 Transient=Transient,
                                 timescale=timescale,
                                 zmax=zmax,
                                 fluxnorm=fluxnorm,
                                 index=index,
                                 LF=LF,
                                 sigma=sigma,
                                 zNEAR=zNEAR,
                                 luminosity=luminosity)
    firesong_simulation(options, outputdir)

if __name__ == "__main__":

    #
    # Check that the Firesong environmental variable is set
    # This is needed for output and to read exposures, effective areas, etc
    try:
        firesongdir = os.environ['FIRESONG']
    except:
        print "Enviromental variable FIRESONG not set"
        quit()
    outputdir = firesongdir + "/Results/"

    #
    # Process command line options
    #
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', action='store',
                        dest='filename',
                        default='Firesong.out',
                        help='Output filename')
    parser.add_argument('-d', action='store', dest='density',
                        type=float, default=1e-9,
                        help='Local neutrino source density [1/Mpc^3]')
    parser.add_argument("--evolution", action="store",
                        dest="Evolution", default='HB2006SFR',
                        help="Source evolution options:  HB2006SFR (default),  NoEvolution")
    parser.add_argument("--transient", action='store_true',
                        dest='Transient', default=False,
                        help='Simulate transient sources, NOT TESTED YET!')
    parser.add_argument("--timescale", action='store',
                        dest='timescale', type=float,
                        default=1000., help='time scale of transient sources, default is 1000sec.')
    parser.add_argument("--zmax", action="store", type=float,
                        dest="zmax", default=10.,
                        help="Highest redshift to be simulated")
    parser.add_argument("--fluxnorm", action="store", dest='fluxnorm',
                        type=float, default=0.9e-8,
                        help="Astrophysical neutrino flux normalization A on E^2 dN/dE = A (E/100 TeV)^(-index+2) GeV/cm^2.s.sr")
    parser.add_argument("--index", action="store", dest='index',
                        type=float, default=2.13,
                        help="Astrophysical neutrino spectral index on E^2 dN/dE = A (E/100 TeV)^(-index+2) GeV/cm^2.s.sr")
    parser.add_argument("--LF", action="store", dest="LF", default="SC",
                        help="Luminosity function, SC for standard candles, LG for lognormal, PL for powerlaw")
    parser.add_argument("--sigma", action="store",
                        dest="sigma", type=float, default=1.0,
                        help="Width of a log normal Luminosity function in dex, default: 1.0")
    parser.add_argument("--zNEAR", action="store", dest="zNEAR",
                        type=float, default=-1,
                        help="Write down a sepaarate file for sources closer than specified redshift. If nothing is specfied, no file is written.")
    parser.add_argument("--L", action="store",
                        dest="luminosity", type=float, default=0.0,
                        help="Set luminosity for each source, will reset fluxnorm option, unit erg/yr")

    options = parser.parse_args()

    firesong_simulation(options, outputdir)
