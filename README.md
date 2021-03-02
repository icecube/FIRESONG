# Firesong
FIRst Extragalactic Simulation Of Neutrinos and Gamma-rays

# Set up
Set up the enviromental variable FIRESONG. In bash: `export FIRESONG=/location/of/firesong/`.

# Basic usage
Several scripts are provided:
* Firesong.py - Generates an instance of all neutrino sources in
  the Universe according to the parameters provided (e.g. local
  neutrino source density). The transient functionality has not been
  tested. 
* FluxPDF.py - Generated the flux probability density distribution of a 
source.
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
