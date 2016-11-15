#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


fname = ['Godunov','HLLC','MacCormack']

inname = 'exact.dat'
xe, re, ue, pe = np.loadtxt(inname, dtype = np.float64, skiprows=1, unpack=True)

for caseid in range(3):
    inname = fname[caseid] + '.dat'
    outname = fname[caseid] + '.pdf'

    x, r, u, p = np.loadtxt(inname, dtype = np.float64, skiprows=1, unpack=True)
    pp = PdfPages(outname)

    fig, axarr = plt.subplots(3, sharex=True)
    axarr[0].set_title(fname[caseid])

    axarr[0].plot(xe, re, '-k')
    axarr[0].plot(x, r, 'or',mec='r',mfc='none',ms=5)
    axarr[0].set_ylabel(r'\rho')

    axarr[1].plot(xe, ue, '-k')
    axarr[1].plot(x, u, 'og',mec='g',mfc='none',ms=5)
    axarr[1].set_ylabel(r'u')

    axarr[2].plot(xe, pe, '-k')
    axarr[2].plot(x, p, 'ob',mec='b',mfc='none',ms=5)
    axarr[2].set_ylabel(r'p')
    axarr[2].set_xlabel(r'x')

    axarr[0].set_ylim([0.75, 1.6])
    axarr[1].set_ylim([-0.02, 0.40])
    axarr[2].set_ylim([0.9, 2.1])

    for i in range(3):
        axarr[i].axvline(x= 8.655,color='k',ls='-', label='Contact Surface')
        axarr[i].axvline(x= -41.75,color='k',ls='--',label='Left Rarefaction Head')
        axarr[i].axvline(x= -31.35,color='k',ls='-.',label='Left Rarefaction Tail')
        axarr[i].axvline(x= 35.35,color='k',ls=':', label='Right Shockwave')
    axarr[0].legend(loc=2, fontsize = 'xx-small')

    plt.tight_layout()
    pp.savefig(bbx_inches = 'tight')
    pp.close()
