# Firesong
FIRst Extragalactic Simulation Of Neutrinos and Gamma-rays.

Documentation for the astrophysics and cosmology of neutrino sources
can be found on:
- [René Rieman's PhD thesis](http://publications.rwth-aachen.de/record/773297),
  section 10.
- [Theo Glauch's Master thesis](https://www.institut3b.physik.rwth-aachen.de/global/show_document.asp?id=aaaaaaaaaavmddj),
  chapter 6.
- Chun Fai Tung PhD thesis (upcoming)

# Set up
FIRESONG depends on cosmolopy, numpy and scipy. Set up the enviromental variable FIRESONG. In bash: `export FIRESONG=/location/of/firesong/`.

# Basic usage
Several scripts are provided:
* Firesong.py - Generates an instance of all neutrino sources in
  the Universe according to the parameters provided (e.g. local
  neutrino source density). Both the flux and the redshift of neutrino
  sources are calculated.
* FluxPDF.py - Generated the flux probability density distribution of a 
source. It complements FIRESONG.py because is it much faster for large
source densities. Only the flux of neutrino sources is calculated.
* Legend.py - 

# Tests
All unittest could be run by

```
python -m unittest discover tests/
```

If you want to run a test for a certain file seperatly user either

```
python -m unittest tests/test_<...>
```

or 

```
python tests/test_<...>.py
```
