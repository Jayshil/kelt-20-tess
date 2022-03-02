import numpy as np
import matplotlib.pyplot as plt
import juliet
import os

tim, fl, fle = juliet.utils.get_all_TESS_data('KELT-20', radius="0.01 deg")

lst = list(tim.keys())

for i in range(len(lst)):
    tim1, fl1, fle1 = tim[lst[i]], fl[lst[i]], fle[lst[i]]
    f1 = open(os.getcwd() + '/Data/KELT-20_sector_' + lst[i] + '.dat', 'w')
    for j in range(len(tim1)):
        f1.write(str(tim1[j]) + '\t' + str(fl1[j]) + '\t' + str(fle1[j]) + '\n')
    f1.close()
    plt.errorbar(tim1, fl1, yerr=fle1, fmt='.')
    plt.show()