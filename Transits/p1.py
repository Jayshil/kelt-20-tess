import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gd
import juliet
import os

# This file is to analyse the individual TESS sectors.
# The approach is first to analyse out of transit/eclipse data (masked both transit and eclipses)
# and then use the posteriors (on GP, instrumental parameters) from this analysis as
# priors in analyse the whole dataset. The analysis of of out-of-transit/eclipse has been performed
# elsewhere, and in this file, we will take posteriors from that analysis as priors.

# Data to be used:
tim, fl, fle = np.loadtxt(os.getcwd() + '/Data/KELT-20_sector_TESS41.dat', usecols=(0,1,2), unpack=True)
instrument = 'TESS41'

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

# Data so that juliet can understand it!
tim_full, fl_full, fle_full = {}, {}, {}
tim_full[instrument], fl_full[instrument], fle_full[instrument] = tim2, fl2, fle2

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
        hyper_gp.append([np.median(post_oot[i]), np.std(post_oot[i])])

## Planetary priors
par_P = ['P_p1', 't0_p1', 'r1_p1', 'r2_p1', 'ecc_p1', 'omega_p1', 'a_p1', 'q1_' + instrument, 'q2_' + instrument]
dist_P = ['normal', 'normal', 'uniform', 'uniform', 'fixed', 'fixed', 'loguniform', 'uniform', 'uniform']
hyper_P = [[P, P_err], [T_0, 0.1], [0., 1.], [0., 1.], 0., 90., [1., 100.], [0.,1.], [0.,1.]]

priors = juliet.utils.generate_priors(par_P+par_ins+par_gp, dist_P+dist_ins+dist_gp, hyper_P+hyper_ins+hyper_gp)

data = juliet.load(priors=priors, t_lc=tim_full, y_lc=fl_full, yerr_lc=fle_full, GP_regressors_lc=tim_full, out_folder=os.getcwd() + '/Transits/Analysis/' + instrument)
res = data.fit(sampler='dynesty', n_live_points=500)

model = res.lc.evaluate(instrument)

# Let's make sure that it works:
fig = plt.figure(figsize=(16,9))
gs = gd.GridSpec(2,1, height_ratios=[2,1])

# Top panel
ax1 = plt.subplot(gs[0])
ax1.errorbar(tim_full[instrument], fl_full[instrument], yerr=fle_full[instrument], fmt='.', alpha=0.3)
ax1.plot(tim_full[instrument], model, c='k', zorder=100)
ax1.set_ylabel('Relative Flux')
ax1.set_xlim(np.min(tim_full[instrument]), np.max(tim_full[instrument]))
ax1.xaxis.set_major_formatter(plt.NullFormatter())

# Bottom panel
ax2 = plt.subplot(gs[1])
ax2.errorbar(tim_full[instrument], (fl_full[instrument]-model)*1e6, yerr=fle_full[instrument]*1e6, fmt='.', alpha=0.3)
ax2.axhline(y=0.0, c='black', ls='--')
ax2.set_ylabel('Residuals (ppm)')
ax2.set_xlabel('Time (BJD)')
ax2.set_xlim(np.min(tim_full[instrument]), np.max(tim_full[instrument]))

plt.savefig(os.getcwd()+'/Transits/Analysis/' + instrument + '/full_model.png')