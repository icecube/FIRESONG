#!/usr/bin/python
#
# Trivial sample interface to HAWC. 
# Assume that all sources at z<0.1 and within
# the field of view are interesting


def write(z,declin,flux,file):
    # The following sources have z<0.1
    # and are in the declination range relevant to HAWC
    if declin>-26. and declin < 64. and z <0.1:
        file.write('{:.4f} {:.4f} {:.4e}\n'.format(declin, z, flux))
