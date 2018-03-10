# Firesong
FIRst Extragalactic Simulation Of Neutrinos and Gamma-rays

# Citation
There's no dedicated publication to FIRESONG, though we are
considering writing one. Presently the best reference for FIRESONG is:

https://arxiv.org/abs/1801.09545

Constrains on the extragalactic origin of IceCube's neutrinos using
HAWC
Ignacio Taboada, Chun Fai Tung, Joshua Wood for the HAWC collaboration
Proc. of the 35th International Cosmic Ray Conference (ICRC2017),
Bexco, Busan, Korea.

# Instructions
Set up the enviromental error FIRESONG so
export FIRESONG=/location/of/firesong/

This is needed to read data files (e.g. exposure data from IceCube)
and to write output data.

Several scripts are provided:
* Firesong.py - Generates an instance of all neutrino sources in
  the Universe according to the parameters provided (e.g. local
  neutrino source density). The transient functionality has not been
  tested. 
* NeutrinoAlert.py - Generates a list of neutrinos that have the
  proper statisitical properties of neutrino alerts produced by
  IceCube. Currently this has ONLY been tested for Standard Candle
  Sources.
* SingleNeutrinoCDF.py - Generates the cummulative distribution
  function of a neutrino alert vs. redshift. This has ONLY been tested
  for Standard Candle sources, so it is equivalent to CDF of neutrino
  alert vs. flux.
* FluxPDF.py - Generated the flux probability density distribution of a 
  source.
* CTA/CTA-Firesong.py - Reads in neutrinos generated with
  NeutrinoAlert.py and calculates a significance for 0.5, 5 and 50
  hour observation. Currenly only CTA-North is implemented. Though
  adding CTA-South is straight forward.
* CTA/CTA-Sensitivity.py - Command line version of CTA/CTA-Firesong.py
  but for a single source and a single observation time. 

# Authors as of May 3, 2017
Chris Tung
Ignacio Taboada

We acknowledge help and ideas by Markus Ahlers, Georga Japaridze,
Konstancja Satalecka and Rene Reimann

# Tagged Versions
v1.2 - June 7, 2017

Different luminosity functions can be used for NeutrinoAlert.py

v1.1 - May 9, 2017

Added the option [--L] to specify luminosity for source. Please input the luminosity in unit of erg/yr.

v1.0 - May 3, 2017
Public Release



v0.2 - beta - January 9, 2017
Major functionality is in place. 

First version ready for public realease.

There are two modes of operation:

Firesong.py : It creates a random instance of all the neutrino sources in the Universe. Steady sources have been the most tested. Transient source functionality is present, but not verified. All luminosity functions and evolution options have been tested
NeutrinoAlert.py : The desired number of IceCube detected neutrinos can be simulated. Steady sources have been the most tested. Transient source functionality is present, but not verified. Currently it only works with standard candle sources.
The CTA/ folder provides an example use case of the output of NeutrinoAlert.py



v0.1 - alpha - December 16, 2016
Major functionality is in place.
Problems to be solved:
* There is private IceCube information that needs to be removed before a
beta release. In particular to be able to share with VERITAS/Magic.
* The processing loop is vectorized into numpy. This is ~30% faster than
a "for" loop in python. However it means that data I/O is done at the
end of the run. Data I/O is a major limitation for the code at high
densities and with many simultanous runs in the cluster. This needs
to be fixed before beta.

# Desired changes
* Make CTA/CTA-Sensitivity.py into a module that can be imported into CTA/CTA-Firesong.py
