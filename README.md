# Firesong
FIRst Extragalactic Simulation Of Neutrinos and Gamma-rays.

[See the Docs](https://icecube.github.io/FIRESONG/)

Documentation for the astrophysics and cosmology of neutrino sources
can be found on:
- [René Rieman's PhD thesis](http://publications.rwth-aachen.de/record/773297),
  section 10.
- [Theo Glauch's Master thesis](https://www.institut3b.physik.rwth-aachen.de/global/show_document.asp?id=aaaaaaaaaavmddj),
  chapter 6.
- Chun Fai Tung PhD thesis (upcoming)

# Set up
This package is developed for Python3, and can be installed via pip:
```
python -m pip install --editable FIRESONG
```

You can also specify where you would like the output of simulations to go, by default. In bash: `export FIRESONG=/location/of/FIRESONG/`.

# Basic usage
Several scripts are provided:
* `firesong/Firesong.py` - Generates an instance of all neutrino sources in
  the Universe according to the parameters provided (e.g. local
  neutrino source density). Both the flux and the redshift of neutrino
  sources are calculated.
* `firsong/FluxPDF.py` - Generated the flux probability density distribution of a 
source. It complements FIRESONG.py because is it much faster for large
source densities. Only the flux of neutrino sources is calculated.
* `firesong/Legend.py` - Generates an instance of gamma-ray sources in the universe
  according to a luminosity denpendent density evolution (LDDE). Both the 
  redshift and the gamma-ray flux (without attenuation) are calculated.

Examples:

* A muon neutrino diffuse flux saturation example:

```
python Firesong.py -d 1e-6 --evolution CC2015SNR --zmax 4.0
--fluxnorm 1.44e-8 --index 2.28 --LF SC
```

wlll simulate neutrino sources with a local density of 10^-6 Mpc^-3
with a source denstiy evolution that folows the Clash and Candels 2015
Supernova Rate (CC2015SNR). The simulation will be done up to a
redshift of 4.0. The neutrino luminosity, because it is not specified
as an option, will be calculated internally to saturate a muon neutrino diffuse
flux with a normalization, at 100 TeV, of E^2d\phi/dE = 1.44 x 10^-8
GeV.cm^-2.s^-1.sr^-1 and with a spectral index of -2.28. Neutrino
luminosoty is distributed as a delta function, i.e., standard candle
(SC).

* An exploration of the luminosity vs. local density plane (aka
Kowalski plot) example:

```
python Firesong.py -d 1e-6 --evolution MD2016SFR --zmax 8.0
--index 2.28 --LF SC -L 1e51
```

will simulate neutrino sources with a local density of 10^-6 Mpc^-3
with a source denstiy evolution that folows the Madau and Dickinson 2016
Star Formation Rate History (MD2014SFR). The simulation will be done up to a
redshift of 8.0.  The neutrino power law spectral index is -2.28. Neutrino
luminosoty is distributed as a delta function, i.e., standard candle (SC)
and is set to be 10^51 erg/year. Note that muon neutrino diffuse flux
normalization is ignored when a luminosity is specified, but the
spectral index should still be provided. This mode of operation allows
the exploration of the luminosity vs. local density plane (aka
Kowalski plot).

More examples are included in the `notebooks` directory.

# Tests
All unittest could be run by

```
python -m unittest discover tests/
```

If you would like to suppress the printed output, you may add a `-b` flag to this command. If you want to run a test for a certain file seperatly use either

```
python -m unittest tests/test_<...>
```

or 

```
python tests/test_<...>.py
```
