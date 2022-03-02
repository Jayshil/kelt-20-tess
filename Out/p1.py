import numpy as np
import matplotlib.pyplot as plt
import juliet
import os

# This file is to analyse the individual TESS sectors.
# In this file we first analyse the out of transit/eclipse data (masked both transit and eclipses)
# And the posteriors from this analysis can be used as priors in the 
# analysis of whole dataset.

# Data to be used:
tim, fl, fle = np.loadtxt(os.getcwd() + '/Data/KELT-20_sector_TESS14.dat', usecols=(0,1,2), unpack=True)


# First planetary parameters:
P, P_err = 	3.4741085, 0.0000019                                 # From Lund et al. 2017
tc, tc_err = 2457503.120049, 0.000190                            # From Lund et al. 2017
t14 = 3.5755/24                                                  # From Lund et al. 2017
cycle = round((tim[0] - tc)/P)
T_0 = tc + (cycle*P)
T_ec = T_0 + (P/2)

# Masking transits:
phs1 = juliet.utils.get_phases(tim, P, T_0)
mask_tra = np.where(np.abs(phs1*P) >= 1.5*t14)[0]

# Masking eclipses:
phs2 = juliet.utils.get_phases(tim, P, T_0 + (P/2))
mask_ecl = np.where(np.abs(phs1*P) >= t14)[0]

# Total mask
mask = np.concatenate((mask_tra, mask_ecl))

tim2, fl2, fle2 = tim[mask], fl[mask], fle[mask]