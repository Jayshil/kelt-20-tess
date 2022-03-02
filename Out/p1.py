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
instrument = 'TESS14'


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
tim2, fl2, fle2 = tim[mask], fl[mask], fle[mask]                # This is our dataset!!

# Convert the dataset into dictionaries so that juliet can understand it
tim_oot, fl_oot, fle_oot = {}, {}, {}
tim_oot[instrument], fl_oot[instrument], fle_oot[instrument] = tim2, fl2, fle2

# Priros:
## Instrumental priors
par_ins = ['mdilution_' + instrument, 'mflux_' + instrument, 'sigma_w_' + instrument]
dist_ins = ['fixed', 'normal', 'loguniform']
hyper_ins = [1., [0., 0.1], [0.1, 10000.]]
## GP priros: ExM kernel
#"""
par_gp = ['GP_sigma_' + instrument, 'GP_timescale_' + instrument, 'GP_rho_' + instrument]
dist_gp = ['loguniform', 'loguniform', 'loguniform']
hyper_gp = [[1e-5, 10000.], [1e-3,1e2], [1e-3,1e2]]
"""
# QP kernel
par_gp = ['GP_B_' + instrument, 'GP_C_' + instrument, 'GP_L_' + instrument, 'GP_Prot_' + instrument]
dist_gp = ['loguniform', 'loguniform', 'loguniform','loguniform']
hyper_gp = [[1e-5,1e3], [1e-5,1e4], [1e-3, 1e3], [1.,1e2]]
"""

priors_oot = juliet.utils.generate_priors(par_ins+par_gp, dist_ins+dist_gp, hyper_ins+hyper_gp)

data = juliet.load(priors=priors_oot, t_lc=tim_oot, y_lc=fl_oot, yerr_lc=fle_oot, GP_regressors_lc=tim_oot, out_folder=os.getcwd() + '/Out/Analysis/' + instrument)
res = data.fit(sampler='dynesty', n_live_points=500)