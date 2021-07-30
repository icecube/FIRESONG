# Firesong
FIRst Extragalactic Simulation Of Neutrinos and Gamma-rays.

[See the Docs](https://icecube.github.io/FIRESONG/)

[See the FIRESONG paper](https://joss.theoj.org/papers/10.21105/joss.03194)

Documentation for the astrophysics and cosmology of neutrino sources
can be found on:
- [René Rieman's PhD thesis](http://publications.rwth-aachen.de/record/773297),
  section 10.
- [Theo Glauch's Master thesis](https://www.institut3b.physik.rwth-aachen.de/global/show_document.asp?id=aaaaaaaaaavmddj),
  chapter 6.
- [Chun Fai Tung's PhD thesis](https://smartech.gatech.edu/handle/1853/64745),
  chapter 4.

# Set up
This package is developed for Python 3, and can be installed via pip:
```
pip install firesong
```
Or by downloading the repository and running:
```
python setup.py install
```
If you want to execute scripts of FIRESONG directly via shell command, you can specify where you would like the output of simulations to go, by default. In bash: `export FIRESONG=/location/of/FIRESONG/`.

# Basic usage
Several scripts are provided:
* `firesong/Firesong.py` - Generates an instance of all neutrino sources in
  the Universe according to the parameters provided (e.g. local
  neutrino source density). Both the flux and the redshift of neutrino
  sources are calculated.
* `firesong/FluxPDF.py` - Generates the flux probability density distribution of a 
source. It complements Firesong.py because is it much faster for large
source densities. Only the flux of neutrino sources is calculated.
* `firesong/Legend.py` - Generates an instance of gamma-ray sources in the universe
  according to a luminosity dependent density evolution (LDDE). Both the 
  redshift and the gamma-ray flux (without attenuation) are calculated.

Examples:

* A muon neutrino diffuse flux saturation example:

If installed via pip, in in the python console,
```
from firesong.Firesong import firesong_simulation
firesong_simulation('./', density=1e-6, Evolution='CC2015SNR', zmax=4.0,
                    fluxnorm=1.44e-8, index=2.28, LF='SC')
```
or with the repository downloaded

```
python Firesong.py -d 1e-6 --evolution CC2015SNR --zmax 4.0
--fluxnorm 1.44e-8 --index 2.28 --LF SC
```

wlll simulate neutrino sources with a local density of 10^-6 Mpc^-3
with a source density evolution that follows the Clash and Candels 2015
Supernova Rate (CC2015SNR). The simulation will be done up to a
redshift of 4.0. The neutrino luminosity, because it is not specified
as an option, will be calculated internally to saturate a muon neutrino diffuse
flux with a normalization, at 100 TeV, of E^2d\phi/dE = 1.44 x 10^-8
GeV.cm^-2.s^-1.sr^-1 and with a spectral index of -2.28. Neutrino
luminosity is distributed as a delta function, i.e., standard candle
(SC). The result will be output as a text file `firesong.out` in the current 
directoy, or in the directory of environment variable  `FIRESONG` if it is set.

* An exploration of the luminosity vs. local density plane (aka
Kowalski plot) example:
in the console
```
firesong_simulation('./', density=1e-6, Evolution='MD2016SFR', zmax=8.0,
                    index=2.28, LF='SC', luminosity=1e51)
```
or

```
python Firesong.py -d 1e-6 --evolution MD2016SFR --zmax 8.0
--index 2.28 --LF SC -L 1e51
```

will simulate neutrino sources with a local density of 10^-6 Mpc^-3
with a source density evolution that follows the Madau and Dickinson 2016
Star Formation Rate History (MD2014SFR). The simulation will be done up to a
redshift of 8.0.  The neutrino power law spectral index is -2.28. Neutrino
luminosity is distributed as a delta function, i.e., standard candle (SC)
and is set to be 10^51 erg/year. Note that muon neutrino diffuse flux
normalization is ignored when a luminosity is specified, but the
spectral index should still be provided. This mode of operation allows
the exploration of the luminosity vs. local density plane (aka
Kowalski plot).

More examples are included in the `notebooks` directory.
`Jupyter notebook` and `matplotlib` are required to run the examples.

# Tests
All unittest could be run by

```
python -m unittest discover tests/
```

If you would like to suppress the printed output, you may add a `-b` flag to this command. If you want to run a test for a certain file separately use either

```
python -m unittest tests/test_<...>
```

or 

```
python tests/test_<...>.py
```
# How to request support

Questions about support for FIRESONG can be sent to one of the
development leads for FIRESONG. See AUTHORS.md

# How to contribute

Community contributions to Firesong are accepted and welcome. Issues
can be reported by any user, even if not a member of IceCube. Pull
requests can be requested by any user, even if not a member of
IceCube.

