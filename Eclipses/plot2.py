import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gd
import pickle
import juliet
from glob import glob
import os

instrument = 'TESS41'

p1 = os.getcwd() + '/Eclipses/Analysis/' + instrument

data = juliet.load(input_folder=p1)
res1 = data.fit(sampler='dynesty')

fl1 = glob(os.getcwd() + '/Eclipses/Analysis/' + instrument + '/*pkl')[0]
post = pickle.load(open(fl1,'rb'), encoding='latin1')
post1 = post['posterior_samples']

tim, fl, fle = data.times_lc[instrument], data.data_lc[instrument], data.errors_lc[instrument]

full_model = res1.lc.evaluate(instrument)
gp_model = res1.lc.model[instrument]['GP']
transit_model = res1.lc.model[instrument]['deterministic']

fac = 1/np.max(transit_model)

tim1, fl1, fle1 = tim, (fl-gp_model)*fac, fle*fac
phase = juliet.utils.get_phases(tim1, 3.4741003123, post1['t0_p1']+(3.4741003123/2))
ii = np.where(phase>=0.9)[0]
phase[ii] = phase[ii]-1.0
idx = np.argsort(phase)

# Making a plot
fig = plt.figure(figsize=(16,9))
gs = gd.GridSpec(2,1, height_ratios=[2,1])

ax1 = plt.subplot(gs[0])
ax1.errorbar(phase, fl1, yerr=fle1, fmt='.', alpha=0.3)
ax1.plot(phase[idx], transit_model[idx]*fac, c='k', zorder=100)
ax1.set_ylabel('Normalised Flux')
#ax1.set_xlim(np.min(tim[instrument]), np.max(tim[instrument]))
ax1.set_xlim([-0.05, 0.05])
ax1.set_ylim([0.999, 1.001])
ax1.xaxis.set_major_formatter(plt.NullFormatter())

# Bottom panel
ax2 = plt.subplot(gs[1])
ax2.errorbar(phase, ((fl1/fac)-transit_model)*1e6*fac, yerr=fle1*1e6, fmt='.', alpha=0.3, zorder=1)
ax2.axhline(y=0.0, c='black', ls='--', zorder=5)
ax2.set_ylabel('Residuals (ppm)')
ax2.set_xlabel('Phase')
ax2.set_xlim([-0.05, 0.05])
ax2.set_ylim([-1000,1000])
#ax2.set_xlim(np.min(tim[instrument]), np.max(tim[instrument]))


#plt.savefig('phase_folded.png')
plt.show()