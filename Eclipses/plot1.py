import numpy as np
import matplotlib.pyplot as plt
import pickle
import os
from glob import glob
import juliet

# This file is to plot the distribution of eclipse depth across the sectors.

# Loding the distributions first
## Sector 14
fl1 = glob(os.getcwd() + '/Eclipses/Analysis/TESS14/*pkl')[0]
post14 = pickle.load(open(fl1, 'rb'), encoding='latin1')
fp14 = (post14['posterior_samples']['fp_p1'])*1e6
qua14 = juliet.utils.get_quantiles(fp14)
L14, L14_u, L14_d = qua14[0], qua14[1]-qua14[0], qua14[0]-qua14[2]
L14, L14_u, L14_d = np.around(L14,4), np.around(L14_u,4), np.around(L14_d, 4)
## Sector 40
fl2 = glob(os.getcwd() + '/Eclipses/Analysis/TESS40/*pkl')[0]
post40 = pickle.load(open(fl2, 'rb'), encoding='latin1')
fp40 = (post40['posterior_samples']['fp_p1'])*1e6
qua40 = juliet.utils.get_quantiles(fp40)
L40, L40_u, L40_d = qua40[0], qua40[1]-qua40[0], qua40[0]-qua40[2]
L40, L40_u, L40_d = np.around(L40,4), np.around(L40_u,4), np.around(L40_d, 4)
## Sector 41
fl3 = glob(os.getcwd() + '/Eclipses/Analysis/TESS41/*pkl')[0]
post41 = pickle.load(open(fl3, 'rb'), encoding='latin1')
fp41 = (post41['posterior_samples']['fp_p1'])*1e6
qua41 = juliet.utils.get_quantiles(fp41)
L41, L41_u, L41_d = qua41[0], qua41[1]-qua41[0], qua41[0]-qua41[2]
L41, L41_u, L41_d = np.around(L41,4), np.around(L41_u,4), np.around(L41_d,4)

plt.figure(figsize=(16,9))
# Sector 14
plt.hist(fp14, bins=50, density=True, histtype='step', color='orangered', label='Sector 14', lw=1.5, zorder=10)
plt.axvline(x=L14, ls='--', color='orangered', alpha=0.7, zorder=5)
plt.axvspan(qua14[1], qua14[2], color='orangered', alpha=0.15, zorder=1)

# Sector 40
plt.hist(fp40, bins=50, density=True, histtype='step', color='cornflowerblue', label='Sector 40', lw=1.5, zorder=10)
plt.axvline(x=L40, ls='--', color='cornflowerblue', alpha=0.7, zorder=5)
plt.axvspan(qua40[1], qua40[2], color='cornflowerblue', alpha=0.15, zorder=1)

# Sector 41
plt.hist(fp41, bins=50, density=True, histtype='step', color='darkgreen', label='Sector 41', lw=1.5, zorder=10)
plt.axvline(x=L41, ls='--', color='darkgreen', alpha=0.7, zorder=5)
plt.axvspan(qua41[1], qua41[2], color='darkgreen', alpha=0.15, zorder=1)

# Text
plt.text(225, 0.0145,\
     r'Eclipse depth in Sector 14: ' + str(L14) + '$^{+' + str(L14_u) + '}_{-' + str(L14_d) + '}$ ppm',\
     color='orangered', fontweight='bold')
plt.text(225, 0.0138,\
     r'Eclipse depth in Sector 40: ' + str(L40) + '$^{+' + str(L40_u) + '}_{-' + str(L40_d) + '}$ ppm',\
     color='cornflowerblue', fontweight='bold')
plt.text(225, 0.0131,\
     r'Eclipse depth in Sector 41: ' + str(L41) + '$^{+' + str(L41_u) + '}_{-' + str(L41_d) + '}$ ppm',\
     color='darkgreen', fontweight='bold')

plt.xlabel(r'Eclipse depth ($L$), in ppm')
plt.ylabel('Density')
plt.title('Posterior distributions of Eclipse depth of KELT-20b across TESS Sectors 14, 40 and 41')
plt.legend()
#plt.show()
plt.savefig(os.getcwd() + '/Eclipses/fp.png')