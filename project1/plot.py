#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


caseid = 3
fname = 'p%d.dat' % caseid

x, r, u, p = np.loadtxt(fname, dtype = np.float64, skiprows=1, unpack=True)


fname = 'p%d.pdf' % caseid
pp = PdfPages(fname)

fig, axarr = plt.subplots(3, sharex=True)
axarr[0].set_title('Problem %d' % caseid)
axarr[0].plot(x, r, 'r')
axarr[0].set_ylabel(r'\rho')

axarr[1].plot(x, u, 'g')
axarr[1].set_ylabel(r'u')

axarr[2].plot(x, p, 'b')
axarr[2].set_ylabel(r'p')
axarr[2].set_xlabel(r'x')

if caseid == 1:
    axarr[0].set_ylim([0, 9])
    axarr[1].set_ylim([-1, 9])
    axarr[2].set_ylim([-50, 550])
    for i in range(3):
        axarr[i].axvline(x= 0.3329,color='k',ls='-', label='Contact Surface')
        axarr[i].axvline(x=-0.3995,color='k',ls='--',label='Left Rarefaction Head')
        axarr[i].axvline(x= 0.0455,color='k',ls='-.',label='Left Rarefaction Tail')
        axarr[i].axvline(x= 0.4505,color='k',ls=':', label='Right Shockwave')
    axarr[0].legend(loc=3, fontsize = 'xx-small')
if caseid == 2:
    axarr[0].set_ylim([-0.1, 1.1])
    axarr[1].set_ylim([-2.25, 2.25])
    axarr[2].set_ylim([-0.05, 0.450])
    for i in range(3):
        axarr[i].axvline(x= 0.0000,color='k',ls='-', label='Contact Surface')
        axarr[i].axvline(x=-0.4105,color='k',ls='--',label='Left Rarefaction Head')
        axarr[i].axvline(x=-0.0505,color='k',ls='-.',label='Left Rarefaction Tail')
        axarr[i].axvline(x= 0.0525,color='k',ls=':', label='Right Rarefaction Tail')
        line = axarr[i].axvline(x= 0.4125,color='k', label='Right Rarefaction Head')
        line.set_dashes([8, 4, 2, 4, 2, 4])
    axarr[0].legend(loc=3, fontsize = 'xx-small')
if caseid == 3:
    axarr[0].set_ylim([0, 7])
    axarr[1].set_ylim([-22, 5])
    axarr[2].set_ylim([-100, 1100])
    for i in range(3):
        axarr[i].axvline(x= 0.3000,color='k',ls='-', label='Contact Surface')
        axarr[i].axvline(x=-0.3830,color='k',ls='--',label='Left Rarefaction Head')
        axarr[i].axvline(x=-0.1025,color='k',ls='-.',label='Left Rarefaction Tail')
        axarr[i].axvline(x= 0.3475,color='k',ls=':', label='Right Shockwave')
    axarr[0].legend(loc=2, fontsize = 'xx-small')


plt.tight_layout()
pp.savefig(bbx_inches = 'tight')
pp.close()
