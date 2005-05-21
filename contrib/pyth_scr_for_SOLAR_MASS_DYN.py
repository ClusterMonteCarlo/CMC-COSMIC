#!/usr/bin/python2.2
# a python script to calculate SOLAR_MASS_DYN for cmc, assuming Salpeter IMF.
#
# SOLAR_MASS_DYN is the mass of Sun in code units,
# it will be equal to 1/<m> where m is given in solar mass.
# e.g., for 0.2-120 Salpeter IMF, <m>~0.69
# therefore SOLAR_MASS_DYN~1.45

from math import *

mmin = 0.2
mmax = 120.0

num = (-2.857142857/pow(mmax, 7.0/20)) - (-2.857142857/pow(mmin, 7.0/20))
denom = (-0.7407407407/pow(mmax, 27.0/20)) - (-0.7407407407/pow(mmin, 27.0/20))

mave = num/denom

SMD = 1.0/mave

print 'average mass =', mave, 'Msolar'
print 'SOLAR_MASS_DYN =', SMD
