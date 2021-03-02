# Firesong
FIRst Extragalactic Simulation Of Neutrinos and Gamma-rays

# Set up
Set up the enviromental error FIRESONG so (say bash)
export FIRESONG=/location/of/firesong/

This is needed to read data files (e.g. exposure data from IceCube)
and to write output data.

# Basic usage
Several scripts are provided:
* Firesong.py - Generates an instance of all neutrino sources in
  the Universe according to the parameters provided (e.g. local
  neutrino source density). The transient functionality has not been
  tested. 
* FluxPDF.py - Generated the flux probability density distribution of a 
  source.

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
