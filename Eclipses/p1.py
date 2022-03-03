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
P, P_err = 3.4741003123, (0.0000002756+0.0000002819)/2           # From Transit analysis
tc, tc_err = 2459441.6682182574, (0.0000237897+0.0000237324)/2   # From Transit analysis
t14 = 0.14741540226270058                                        # From Transit analysis
rprs = 0.1155448549                                              # From Transit analysis
aR, aR_err = 7.4609895799, (0.0209709596+0.0215888273)/2         # From Transit analysis
bb, bb_err = 0.5155183998, (0.0049245760+0.0051356091)/2         # From Transit analysis
ac = 0.00031362                                                  # From Transit analysis

cycle = round((tim[0] - tc)/P)
T_0 = tc + (cycle*P)
T_ec = T_0 + (P/2)

# Masking transits:
phs2 = juliet.utils.get_phases(tim, P, T_0)
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
for j in post_oot.keys():
    if j[0:2] == 'GP':
        par_gp.append(j)
        dist_gp.append('normal')
        hyper_gp.append([np.median(post_oot[j]), np.std(post_oot[j])])
## Planetary priors
par_P = ['P_p1', 't0_p1', 'p_p1', 'b_p1', 'ecc_p1', 'omega_p1', 'a_p1', 'fp_p1', 'ac_p1']
dist_P = ['fixed', 'normal', 'fixed', 'normal', 'fixed', 'fixed', 'normal', 'uniform', 'fixed']
hyper_P = [P, [T_0, tc_err], rprs, [bb, bb_err], 0., 90., [aR, aR_err], [1.e-6, 500.e-6], ac]

priors = juliet.utils.generate_priors(par_P+par_ins+par_gp, dist_P+dist_ins+dist_gp, hyper_P+hyper_ins+hyper_gp)

data = juliet.load(priors=priors, t_lc=tim_full, y_lc=fl_full, yerr_lc=fle_full, GP_regressors_lc=tim_full, out_folder=os.getcwd() + '/Eclipses/Analysis/' + instrument)
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

plt.savefig(os.getcwd()+'/Eclipses/Analysis/' + instrument + '/full_model.png')