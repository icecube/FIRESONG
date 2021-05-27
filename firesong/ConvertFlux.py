import numpy as np

def nu_to_gamma(nuFlux, index, interaction):

    """
    Convert the neutrino flux to gamma-ray flux

    Args:
        nuFlux (float or array): neutrino flux in GeV cm^-2 s^-1
        index (float): spectral index of a power law spectrum
        interaction (string): the interaction type. Options: pgamma, pp

    Returns:
        gammaFlux (float or array): gamma-ray flux in TeV^-1 cm^-2 s^-1
    """

    if interaction == "pgamma":
        gammaFlux = nuFlux * 10**(-7) * 200**index
    elif interaction == "pp":
        gammaFlux = nuFlux * 0.5 * 10**(-7) * 200**index
    #else:
        # only have pgamma and pp interactions types for now

    return gammaFlux

def gamma_to_nu(gammaFlux, index, interaction):

    """
    Convert gamma-ray flux to neutrino flux

    Args:
        gammaFlux (float or array): gamma-ray flux in TeV^-1 cm^-2 s^-1
        index (float): spectral index of a power law spectrum
        interaction (string): the interaction type. Options: pgamma, pp

    Returns:
        nuFlux (float or array): neutrino flux in GeV cm^-2 s^-1
    """

    if interaction == "pgamma":
        nuFlux = gammaFlux * 10**7 * 200**(-index)
    elif interaction == "pp":
        nuFlux = gammaFlux * 2 * 10**7 * 200**(-index)
#    else:
        # only have pgamma and pp interaction types for now

    return nuFlux
        
