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
Set up the enviromental error FIRESONG so (say bash)
export FIRESONG=/location/of/firesong/

This is needed to read data files (e.g. exposure data from IceCube)
and to write output data.

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
