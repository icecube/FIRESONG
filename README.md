# Firesong
FIRst Extragalactic Simulation Of Neutrinos and Gamma-rays

# Instructions
Set up the enviromental error FIRESONG so
export FIRESONG=/location/of/firesong/

This is needed to read data files (e.g. exposure data from IceCube)
and to write output data.

# Authors as of December 5, 2016
Chris Tung
Ignacio Taboada

We acknowledge help and ideas by Markus Ahlers and Georga Japaridze

# Tagged Versions

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
