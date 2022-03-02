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
instrument = 'TESS14'

# First planetary parameters:
P, P_err = 	3.4741085, 0.0000019                                 # From Lund et al. 2017
tc, tc_err = 2457503.120049, 0.000190                            # From Lund et al. 2017
t14 = 3.5755/24                                                  # From Lund et al. 2017
cycle = round((tim[0] - tc)/P)
T_0 = tc + (cycle*P)
T_ec = T_0 + (P/2)

# Masking eclipses:
phs2 = juliet.utils.get_phases(tim, P, T_0 + (P/2))
mask_ecl = np.where(np.abs(phs2*P) >= t14)[0]
tim2, fl2, fle2 = tim[mask_ecl], fl[mask_ecl], fle[mask_ecl]

# Priros:
data_oot = juliet.load(input_folder=os.getcwd() + '/Out/Analysis/' + instrument)
res_oot = data_oot.fit()
post_oot = res_oot.posteriors['posterior_samples']
## Instrumental: from Out-of-transit/eclipse analysis
par_ins = ['mdilution_' + instrument, 'mflux_' + instrument, 'sigma_w_' + instrument]
dist_ins = ['fixed', 'normal', 'normal']
hyper_ins = [1.0, [np.median(post_oot['mflux_' + instrument]), np.std(post_oot['mflux_' + instrument])],\
     [np.median(post_oot['sigma_w_' + instrument]), np.std(post_oot['sigma_w_' + instrument])]]
## GP Priors: from out-of-transit/eclipse analysis
par_gp, dist_gp, hyper_gp = [], [], []
for i in post_oot.keys():
    if i[0:2] == 'GP':
        par_gp.append(i)
        dist_gp.append('normal')
        hyper_gp.append([np.median(post_oot(i)), np.std(post_oot(i))])

## Planetary priors
par_P = ['P_p1', 't0_p1', 'r1_p1', 'r2_p1', 'ecc_p1', 'omega_p1', 'a_p1', 'q1_' + instrument, 'q2_' + instrument]
dist_P = ['normal', 'normal', 'uniform', 'uniform', 'fixed', 'fixed', 'loguniform', 'uniform', 'uniform']
hyper_P = [[P, P_err], [T_0, 0.1], [0., 1.], [0., 1.], 0., 90., [1., 100.], [0.,1.], [0.,1.]]

priors = juliet.utils.generate_priors(par_P+par_ins+par_gp, dist_P+dist_ins+dist_gp, hyper_P+hyper_ins+hyper_gp)

data = juliet.load(priors=priors, t_lc=tim2, y_lc=fl2, yerr_lc=fle2, GP_regressors_lc=tim2, out_folder=os.getcwd() + '/Transits/Analysis/' + instrument)