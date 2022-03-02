import numpy as np
import matplotlib.pyplot as plt
import juliet
import os

# This file is to analyse the individual TESS sectors.
# The approach is first to analyse out of transit/eclipse data (masked both transit and eclipses)
# and then use the posteriors (on GP, instrumental parameters) from this analysis as
# priors in analyse the whole dataset. The analysis of of out-of-transit/eclipse has been performed
# elsewhere, and in this file, we will take posteriors from that analysis as priors.

# Data to be used:
tim, fl, fle = np.loadtxt(os.getcwd() + '/Data/KELT-20_sector_TESS14.dat', usecols=(0,1,2), unpack=True)


# First planetary parameters:
P, P_err = 	3.4741085, 0.0000019                                 # From Lund et al. 2017
tc, tc_err = 2457503.120049, 0.000190                            # From Lund et al. 2017
t14 = 3.5755/24                                                  # From Lund et al. 2017

# Masking eclipses:
