from cProfile import label
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
t014 = np.median(post14['posterior_samples']['t0_p1'])
fp14 = (post14['posterior_samples']['fp_p1'])*1e6
qua14 = juliet.utils.get_quantiles(fp14)
L14, L14_u, L14_d = qua14[0], qua14[1]-qua14[0], qua14[0]-qua14[2]
L14, L14_u, L14_d = np.around(L14,4), np.around(L14_u,4), np.around(L14_d, 4)
## Sector 40
fl2 = glob(os.getcwd() + '/Eclipses/Analysis/TESS40/*pkl')[0]
post40 = pickle.load(open(fl2, 'rb'), encoding='latin1')
t040 = np.median(post40['posterior_samples']['t0_p1'])
fp40 = (post40['posterior_samples']['fp_p1'])*1e6
qua40 = juliet.utils.get_quantiles(fp40)
L40, L40_u, L40_d = qua40[0], qua40[1]-qua40[0], qua40[0]-qua40[2]
L40, L40_u, L40_d = np.around(L40,4), np.around(L40_u,4), np.around(L40_d, 4)
## Sector 41
fl3 = glob(os.getcwd() + '/Eclipses/Analysis/TESS41/*pkl')[0]
post41 = pickle.load(open(fl3, 'rb'), encoding='latin1')
t041 = np.median(post41['posterior_samples']['t0_p1'])
fp41 = (post41['posterior_samples']['fp_p1'])*1e6
qua41 = juliet.utils.get_quantiles(fp41)
L41, L41_u, L41_d = qua41[0], qua41[1]-qua41[0], qua41[0]-qua41[2]
L41, L41_u, L41_d = np.around(L41,4), np.around(L41_u,4), np.around(L41_d,4)

plt.figure(figsize=(16/1.5,9/1.5))
plt.errorbar(np.array([t014]), np.array([L14]), yerr=np.array([L14_d, L14_u]).reshape((2,1)), fmt='o', c='orangered', mfc='white', label='Sector 14')
plt.errorbar(np.array([t040]), np.array([L40]), yerr=np.array([L40_d, L40_u]).reshape((2,1)), fmt='o', c='cornflowerblue', mfc='white', label='Sector 40')
plt.errorbar(np.array([t041]), np.array([L41]), yerr=np.array([L41_d, L41_u]).reshape((2,1)), fmt='o', c='darkgreen', mfc='white', label='Sector 41')

plt.xlabel('Transit time (in BJD)')
plt.ylabel('Eclipse depth (in ppm)')

plt.legend(loc='upper left')
plt.title('TESS eclipse depths of KELT-20b across Sector 14, 40, and 41')
plt.grid()

#plt.show()
plt.savefig(os.getcwd() + '/Eclipses/fp1.png')