---
title:  'FIRESONG: A python package to simulate populations of extragalactic neutrino sources'
tags:
   - Python
   - Neutrinos
   - Neutrino Sources
   - Cosmic Rays
   - Multi Messenger Astrophysics
   - Cosmology
authors:
  - name: Chun Fai Tung
    affiliation: 1
  - name: Theo Glauch
    affiliation: 2
  - name: Michael Larson
    orcid: 0000-0002-6996-1155
    affiliation: 3
  - name: Alex Pizzuto
    orcid: 0000-0002-8466-8168
    affiliation: 4
  - name: Rene Reimann
    orcid: 0000-0002-1983-8271
    affiliation: 5
  - name: Ignacio Taboada
    orcid: 0000-0003-3509-3457
    affiliation: 1
affiliations:
  - name: School of Physics. Georgia Institute of Technology. Atlanta, GA 30332, USA
    index: 1
  - name: TUM
    index: 2
  - name: Dept. of Physics, University of Maryland, College Park, MD 20742, USA
    index: 3
  - name: Dept. of Physics and Wisconsin IceCube Particle Astrophysics Center, University of Wisconsin–Madison, Madison, WI 53706, USA
    index: 4
  - name: Johannes Gutenberg University Mainz, Institute of Physics - QUANTUM, 55128 Mainz, Germany
    index: 5
date: 26 February 2021
bibliography: paper.bib
---

# Summary

Neutrinos provide a new perspective on the universe. Due to their weak
interaction with matter, neutrinos carry information from places where
electromagnetic radiation, e.g., gamma rays, cannot
escape. Though astrophysical neutrinos have been detected, the class
or classes of objects that produce them have not been unequivocally
identified. ``FIRESONG`` simulations populate the universe with
neutrino sources. These simulations can be used, among other things,
to study if a given class of astronomical sources is viable
to explain measured astrophysical neutrinos. 

# Background

The IceCube neutrino observatory has discoved an all-sky neutrino flux
in the 10 TeV to 10 PeV energy range
[@IceCube:2019a;@IceCube:2019b;@IceCube:2020a]. IceCube
finds that a power law in energy is a good description of the flux,
with an spectral index ranging from -2.28 to -2.89, depending on the
observation channel used. This flux is apparently isotropic,
consistent with an extragalactic origin for these neutrinos. The flux
is also consistent with equal flux for each of the three neutrino
flavors [@IceCube:2019c], as expected for standard neutrino oscillations over astrophysical
baselines. The origin of this flux is of great scientific interest as
it is expected that neutrino sources are also sources of
ultra-high-energy cosmic rays, which also have an unknown origin. 
IceCube has identified the blazar, a sub-type of Active Galactic
Nuclei (AGN), TXS 0506+056 as a candidate neutrino source
[@IceCube:2018a;@IceCube:2018b]. However, there's also evidence, by
IceCube that gamma-ray bright blazars contribute to no more than 
approximately 27% of the diffuse flux [@IceCube:2017a]. More recently
IceCube has found a neutrino point source hot-spot, just below the 3
sigma threshold normally assigned to evidence, correlated with the
Seyfert II galaxy, another subtype of AGN, NGC 1068 [@IceCube:2020b]. Over
the past 30 years, AGNs and Gamma Ray Bursts (GRBs) were among the
most prominent proposed extragalactic neutrino sources. IceCube has
ruled out GRBs as constributing more than 1% of the diffuse flux
[@IceCube:2015].

The properties of various proposed extragalactic neutrino
source and/or reservoir classes, such as starburst galaxies,
blazars, low luminosity GRBs, Flat Spectrum Radio Quasars, BL Lacs and
galaxy clusters can be summarized in terms of the local density (or
density rate for transient sources) as a function of luminosity (or
per-burst equivalent isotropic energy for transient sources)
[@Kowalski:2014;@MuraseWaxman:2016]. The correct description of each of these classes of 
objects depends on, e.g., the redshift evolution of the density of
sources; but more generally on the luminosity function of the
objects. The existence of a diffuse extragalactic neutrino flux can be
described as a inverse relationship between density (density-rate) and
luminosity (isotropic equivalent energy) [@Kowalski:2014]. This relationship also
depends on the evolution assumed.

The identification of the main sources of the diffuse neutrino flux remains an
open research topic.

# Statement of Need

``FIRESONG`` is a python package to be used by researchers interested in
simulating populations of neutrino sources in the universe and to put
these simulations in the context of IceCube's observation of a diffuse
neutrino flux. The calculations needed to conduct these simulations
are well established 
but also cumbersome and error prone. Indeed several authors have similar
(usually private) code. ``FIRESONG`` provides an open source 
mantained framework for these simulations.  ``FIRESONG`` depends on
the ``cosmolopy`` package [@cosmolopy] 
for cosmological calculations. ``FIRESONG`` also depends on ``numpy``
and ``scipy``. ``FIRESONG`` has already been used on scientific
publication by several observatories of neutrinos or gamma rays:
IceCube [@IceCube:2019d], HAWC and IceCube [@HAWCIceCube:2021], 
HAWC [@HAWC:2018] and CTA [@FiresongCTA:2019]. Though originally concieved
as a stand alone project, maintanance of ``FIRESONG`` is currently
provided by IceCube collaboration members.

# Usage

``FIRESONG`` can be invoked from the command line as ``Firesong.py`` and
configured via command line options outputing a file with a simulated list of
neutrino sources. Alternatively ``FIRESONG`` can also be imported and
produce a python dictionary of the simulated neutrino
sources. ``FIRESONG`` can be used to simulate steady or transient
sources. If no luminosity (isotropic equivalent energy) is provided,
``FIRESONG`` calculates it, as a function of local density (density
rate) and other parameters, so that the IceCube neutrino diffuse flux is fully
saturated. Lack of knowledge of the properties of neutrino sources
motivate simple choices for implemented luminosity distributions: a
delta function (standard candle), a lognormal distribution or a power
law distribution. Various models of Star Formation History are
implemented as well as no evolution.

Legend.py can also be invoked from the command line as ```Lengend.py``` 
and configured in a similar way as ```Firesong.py```. It can also be 
executed in the Python console by importing the function 
```legend_simulation``` from ```Legend```. If invoked as a function,
the output will be a dictionary if the filename option is set to ```None```.
The output dictionary contains the declinations, distances, and fluxes 
of the simulated sources. Simulation of transient sources is not supported 
by ```Legend```. Since ```Legend``` is motivated by Luminosity Dependent 
Density Evolution (LDDE), the distribution of luminosities is decided by 
the evolution model.  


Luminosity functions provide the source density as a function of
source luminosity and cosmological redshift. For observational
purposes, however, we usually care about the source count
distribution, i.e., a function giving the total number of sources with
a specific flux at Earth. The ``FluxPDF.py`` of ``FIRESONG``
calculates a smooth source count distribution by marginalising over
any luminosity function and summing up all the contributions after
accounting for their distance. Once generated, this 1D distribution
can be further used to generate specific realisations of the
luminosity function. This is extremely fast, but doesn’t provide any
information on the sources original redshifts. In that sense it is
complementary to the sampling of ``Firesong.py`` and specifically
useful for cases where the density of sources is extremely high, when
``Firesong.py`` is CPU intensive.

# Acknowledgements

We acknowldge comments, support and ideas by:  Markus Ahlers, George
Japaridze, Konstancja Satalecka and the IceCube collaboration. 
CFT, ML, AP and IT acknowledge support by NSF grant PHY-1913607.

# References
