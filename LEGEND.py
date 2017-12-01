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
import scipy.integrate
# Firesong code
from Evolution import LuminosityEvolution, LuminosityDistance


def legend_simulation(options, output):
    if re.search('.gz$', options.filename):
        output = gzip.open(outputdir+str(options.filename), 'wb')
    else:
        output = open(outputdir+str(options.filename),"w")

    if (options.zNEAR>0):
        if re.search('.gz$', options.filename):
            near_output = gzip.open(outputdir + "Near_" + options.filename,"w")
        else:
            near_output = open(outputdir + "Near_" + options.filename,"w")

    print ("##############################################################################")
    print ("##### FIRESONG-LEGEND initializing #####")
    print ("##### Luminosity Evolution Model: " +str(options.le_model))
    print ("##### Redshift range: 0 - " + str(options.zmax)) 
    print ("##### FIRESONG initialization done #####")

    ##################################################
    #        Simulation starts here
    ##################################################

    output.write("# Desired neutrino diffuse flux:\n")
    output.write("#      E^2 dN_{diffuse}/dE = " + str(options.fluxnorm) + " (E/100 TeV)^(" + str(-(options.index-2.)) + ") [GeV/cm^2.s.sr]\n") 
    output.write("# Neutrino point source fluxes listed below are of \'A\' where the flux is:\n")
    output.write("#      E^2 dN_{PS}/dE = A * (E/100 TeV)^(" + str(-(options.index-2.)) + ") [GeV/cm^2.s]\n") 
    output.write("# Dec(deg) Redshift  A\n")

    ## Luminosity Evolution returns the redshift bin, luminosity bin, number of source at each z bin, and luminosity CDF
    redshift_bins, luminosity_bins, nz, L_cdf = LuminosityEvolution(options)
    Nsource = sum(nz)
    print ("##### {} Blazars in the sky".format(Nsource))

    TotalFlux = 0
    ## Generate declin, L and z
    for index in range(len(nz)):
        if nz[index] > 0:
            #generate nz[index] declination
            sinDec = 2*np.random.rand(nz[index])-1
            declin = 180*np.arcsin(sinDec)/np.pi
            #generate nz[index] luminosity(according to CDF) and z
            test = np.random.rand(nz[index])
            bin_index_l = np.searchsorted(L_cdf[index], test)
            L = np.array(luminosity_bins[bin_index_l])
            z = np.array([redshift_bins[index]]*nz[index])
            #photon flux from radiation luminosity, pivoted at 1GeV, doppler shifted
            pflux = (1+z)**(2-options.index)*(10**L * 624.151)*(2-options.index)/((100)**(2-options.index)-(0.1)**(2-options.index))/4./np.pi/(LuminosityDistance(z)*3.086e24)**2.
            #neutrino flux from radiation photon flux, and change the pivot energy to 100GeV as IceCube convention, assuming 1 to 1 flux ratio
            flux = pflux*(1e5)**(2-options.index)

            for i in range(nz[index]):
                output.write('{:.4f} {:.4f} {:.4e}\n'.format(declin[i], z[i], flux[i]))
                if (z[i]<options.zNEAR):
                    near_output.write('{:.4e} {:.4f} {:.4f}\n'.format(declin[i], flux[i], z[i]))
                TotalFlux += flux[i]

    print "Actual diffuse flux simulated :  E^2 dNdE = " + str(TotalFlux/(4*np.pi)) + " (E/100 TeV)^(" + str(-(options.index-2.)) + ") [GeV/cm^2.s.sr]"
    print "IC Diffuse flux conversion Factor =", options.fluxnorm/(TotalFlux/(4*np.pi))
    output.write("# Diffuse flux E^2 dNdE = " + str(TotalFlux/(4*np.pi)) + "\n")
    output.write('# IC Diffuse flux conversion Factor = '+str(options.fluxnorm/(TotalFlux/(4*np.pi))))
    output.close()
    if (options.zNEAR>0):
        near_output.close()


if __name__ == '__main__':
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
    parser.add_argument('-o', action='store', dest='filename',default= 'Legend.out',
                        help='Output filename')
    parser.add_argument("--model", action="store",
                        dest="le_model", default='blazar',
                        help="evolution model, default is blazar")
    parser.add_argument("--Lmin", action="store", dest="lmin", type=float, default=40,
                        help="Minimum log(luminosity) of source, default is 40")
    parser.add_argument("--Lmax", action="store", dest="lmax", type=float, default=48,
                        help="Maximum log(luminosity) of source, default is 48")
    #parser.add_argument("--transient", action='store_true',
    #                    dest='Transient', default=False,
    #                    help='Simulate transient sources, NOT TESTED YET!')
    #parser.add_argument("--timescale", action='store', dest='timescale', type=float,
    #                    default=1000., help='time scale of transient sources, default is 1000sec.')
    parser.add_argument("--zmax", action="store", type=float,
                        dest="zmax", default=10.,
                        help="Highest redshift to be simulated")
    parser.add_argument("--fluxnorm", action="store", dest='fluxnorm', type=float, default=1.01e-8,
                        help="Astrophysical neutrino flux normalization A on E^2 dN/dE = A (E/100 TeV)^(-index+2) GeV/cm^2.s.sr")
    parser.add_argument("--index", action="store", dest='index', type=float, default=2.19,
                        help="Astrophysical neutrino spectral index on E^2 dN/dE = A (E/100 TeV)^(-index+2) GeV/cm^2.s.sr")
    parser.add_argument("--zNEAR", action="store",dest="zNEAR", type=float,
                        default=-1, help="Write down a separate file for sources closer than specified redshift. If nothing is specfied, no file is written.")
    options = parser.parse_args()

    legend_simulation(options, outputdir)

