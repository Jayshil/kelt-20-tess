import numpy as np
import juliet
import os
import pickle
from glob import glob

# This file is to do a multisector fit
# We do only for in-transit points

sec = ['TESS14', 'TESS40', 'TESS41']

# Importing data
tim, fl, fle = {}, {}, {}
tim14, fl14, fle14 = np.loadtxt(os.getcwd() + '/Data/KELT-20_sector_TESS14.dat', usecols=(0,1,2), unpack=True)
tim40, fl40, fle40 = np.loadtxt(os.getcwd() + '/Data/KELT-20_sector_TESS40.dat', usecols=(0,1,2), unpack=True)
tim41, fl41, fle41 = np.loadtxt(os.getcwd() + '/Data/KELT-20_sector_TESS41.dat', usecols=(0,1,2), unpack=True)

# Planetary parameters
P, P_err = 	3.4741085, 0.0000019                                 # From Lund et al. 2017
tc, tc_err = 2457503.120049, 0.000190                            # From Lund et al. 2017
t14 = 3.5755/24                                                  # From Lund et al. 2017

# Only taking in-transit points
## Sector 14
n14 = round((tim14[0]-tc)/P)
tc14 = tc + (n14+1)*P
phs14 = juliet.utils.get_phases(tim14, P, tc14)
mask14 = np.where(np.abs(phs14*P) <= t14)[0] 
## Sector 40
n40 = round((tim40[0]-tc)/P)
tc40 = tc + (n40+1)*P
phs40 = juliet.utils.get_phases(tim40, P, tc40)
mask40 = np.where(np.abs(phs40*P) <= t14)[0]
## Sector 41
n41 = round((tim41[0]-tc)/P)
tc41 = tc + (n41+1)*P
phs41 = juliet.utils.get_phases(tim41, P, tc41)
mask41 = np.where(np.abs(phs41*P)<= t14)[0]

tim['TESS14'], fl['TESS14'], fle['TESS14'] = tim14[mask14], fl14[mask14], fle14[mask14]
tim['TESS40'], fl['TESS40'], fle['TESS40'] = tim40[mask40], fl40[mask40], fle40[mask40]
tim['TESS41'], fl['TESS41'], fle['TESS41'] = tim41[mask41], fl41[mask41], fle41[mask41]

# For planetary priros
cycle = round((tim['TESS41'][-1]-tc)/P)
T_0 = tc + (cycle-1)*P

# Priors: instrumental and GP
par_ins, dist_ins, hyper_ins = [], [], []
par_gp, dist_gp, hyper_gp = [], [], []
## LDCs
q1, q2 = 'q1_', 'q2_'

for i in range(len(sec)):
    par_ins.append('mdilution_' + sec[i])
    dist_ins.append('fixed')
    hyper_ins.append(1.0)
    q1 = q1 + sec[i] + '_'
    q2 = q2 + sec[i] + '_'
q1, q2 = q1[:-1], q2[:-1]

for i in range(len(sec)):
    pkl = glob(os.getcwd() + '/Transits/Analysis/' + sec[i] + '/*pkl')[0]
    post = pickle.load(open(pkl, 'rb'), encoding='latin1')
    post1 = post['posterior_samples']
    for j in post1.keys():
        if j[0:2] == 'GP':
            par_gp.append(j)
            dist_gp.append('normal')
            hyper_gp.append([np.median(post1[j]), np.std(post1[j])])
        elif j[0:2] == 'mf':
            par_ins.append(j)
            dist_ins.append('normal')
            hyper_ins.append([np.median(post1[j]), np.std(post1[j])])
        elif j[0:2] == 'si':
            par_ins.append(j)
            dist_ins.append('normal')
            hyper_ins.append([np.median(post1[j]), np.std(post1[j])])

par_P = ['P_p1', 't0_p1', 'r1_p1', 'r2_p1', 'ecc_p1', 'omega_p1', 'a_p1', q1, q2]
dist_P = ['normal', 'normal', 'uniform', 'uniform', 'fixed', 'fixed', 'loguniform', 'uniform', 'uniform']
hyper_P = [[P, P_err], [T_0, 0.1], [0., 1.], [0., 1.], 0., 90., [1., 100.], [0.,1.], [0.,1.]]

priors = juliet.utils.generate_priors(par_P+par_ins+par_gp, dist_P+dist_ins+dist_gp, hyper_P+hyper_ins+hyper_gp)

data = juliet.load(priors=priors, t_lc=tim, y_lc=fl, yerr_lc=fle, out_folder=os.getcwd() + '/Transits/Analysis/MultiSector')
res = data.fit(sampler='dynesty', n_live_points=500)