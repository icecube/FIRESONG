# Tagged Versions
v1.5 - March 28, 2018

Rewritten, added new model, new mode

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

