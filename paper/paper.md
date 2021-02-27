---
title: 'FIRESONG: A python package to simulate populations of extragalactic neutrino sources'
tags:
  - Python
  - Neutrinos
  - Multi-messenger Astrophysics
  - Cosmology
authors:
  - name: Chun Fai Tung
     orcid: XXXXXXXX
     affiliation: 1
  - name: Theo Glauch
     orcid: XXXXXXXX
     affiliation: 2
  - name: Michael Larson
     orcid: XXXXXXXX
     affiliation: 3
  - name: Alex Pizzuto
     orcid: XXXXXXXX
     affiliation: 4
  - name: René Reiman
     orcid: XXXXXXXX
     affiliation: 5
  - name: Ignacio Taboada
     orcid: 0000-0003-3509-3457
     affiliation:1
affiliations:
  - name: School of Physics. Georgia Institute of Technology. Atlanta, GA 30332, USA
     index: 1
  - name: TUM
     index: 2
  - name: UMd
     index: 3
  - name: UW Madison
     index: 4
  - name: I don't know
date: 26 February 2021
bibliography: paper.bib

---

# Background

The IceCube neutrino observatory has discoved an all-sky neutrino flux
in the 10 TeV to 10 PeV energy range
[@PhysRevLett.125.121104,@2019ICRC...36.1004S;@2019ICRC...36.1017S]. IceCube
finds that a power law in energy is a good description of the flux,
with an spectral index ranging from -2.28 to -2.89, depending on the
observation channel used. This flux is apparently isotropic,
consistent with an extragalactic origin for these neutrinos. The flux
is also consistent with equal flux for each of the three neutrino
flavors [@somebody], as consistent with 
expectations for standard neutrino oscillations over astrophysical
baselines. The origin of this flux is of great scientific interest as
it is expected that neutrino sources are also sources of
ultra-high-energy cosmic rays, which also have an unknown origin. 
IceCube has identified the blazar, a sub-type of Active Galactic
Nuclei (AGN), TXS 0506+056 as a candidate neutrino source
[@somebody]. However, there's also evidence, also by IceCube that gamma-ray bright
blazars contribute to no more than 20% of the diffuse flux
[@somebody]. More recently IceCube has found a neutrino point source hot-spot, just below
the 3 sigma threshold normally assigned to evidence, correlated with
the Seyfert II galaxy, another subtype of AGN, NGC
1068 [@somebody]. Over the past 30 years, AGNs and Gamma Ray Bursts
(GRBs) were among the most prominent proposed extragalactic neutrino
sources. IceCube has ruled out GRBs as constributing more than 1% of
the diffuse flux [@somebody].

The properties of various proposed extragalactic neutrino
source (or also reservoir) classes, such as starburst galaxies,
blazars, low luminosity GRBs, Flat Spectrum Radio Quasars, BL Lacs and
galaxy clusters can be summarized in terms of the local density (or
density rate for transient sources) as a function of luminosity (or
per-burst equivalent isotropic energy for transient sources)
[@somebody]. The correct description of each of these classes of 
objects depends on, e.g., the redshift evolution of the density of
sources; but more generally on the luminosity function of the
objects. The exitance of a diffuse extragalactic neutrino flux can be
described as a inverse relationship between density (density-rate) and
luminosity (isotropic energy) [@somebody]. This relationship also
depends on the evolution assumed. Note by Ignacio: We can provide an
example Kowalski plot - but this is not critically needed.

The identification of the main sources of the diffuse flux remains an
open research topic.

# Summary

Summary FIRESONG

Code is tailored for muon neutrinos, but can trivially be adapted by
a user to include all flavors. 

# Usage

Give a brief description on how to use Firesong.

# Statement of Need

Why do we need this

# Acknowledgements

We acknowledge: Markus Ahlers (he gave Ignacio the idea, though
probably he has forgotten), Konstancja Satalecka, George Japaridze and
the IceCube collaboration.
Note, we should consult these people!
IT acknowledges grant Blah.

