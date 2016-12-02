#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


fname = ['ZERO','AVG','MINMOD']
fname = ['MINMOD','SUPERBEE', 'VANLEERNS']
col   = ['b','r','g','c']
symb  = ['o','s','v','^']

inname = 'exact.dat'
xe, re, ue, pe = np.loadtxt(inname, dtype = np.float64, skiprows=1, unpack=True)
inname = 'exact100.dat'
xe100, re100, ue100, pe100 = np.loadtxt(inname, dtype = np.float64, skiprows=1, unpack=True)

pp = PdfPages('task12.pdf')
fig, axarr = plt.subplots(3, sharex=True,figsize=(8, 8))
for caseid in range(3):
    inname = fname[caseid] + '.dat'
    outname = fname[caseid] + '.pdf'

    x, r, u, p = np.loadtxt(inname, dtype = np.float64, skiprows=1, unpack=True)

#   axarr[0].set_title(fname[caseid])

    axarr[0].plot(xe, re, '-k')
    axarr[0].plot(x, r, symb[caseid],mec=col[caseid],mfc='none',ms=3, label=fname[caseid])
    axarr[0].set_ylabel(r'\rho')

    axarr[1].plot(xe, ue, '-k')
    axarr[1].plot(x, u, symb[caseid],mec=col[caseid],mfc='none',ms=3)
    axarr[1].set_ylabel(r'u')

    axarr[2].plot(xe, pe, '-k')
    axarr[2].plot(x, p, symb[caseid],mec=col[caseid],mfc='none',ms=3)
    axarr[2].set_ylabel(r'p')
    axarr[2].set_xlabel(r'x')

    axarr[0].set_ylim([0.75, 1.6])
    axarr[1].set_ylim([-0.02, 0.40])
    axarr[2].set_ylim([0.9, 2.1])
    axarr[0].set_xlim([-50, 50])
    axarr[1].set_xlim([-50, 50])
    axarr[2].set_xlim([-50, 50])

    print 'error'
    print np.linalg.norm( (re100-r) ,ord=2)/ (104**(0.5))
    print np.linalg.norm( (ue100-u) ,ord=2)/ (104**(0.5))
    print np.linalg.norm( (pe100-p) ,ord=2)/ (104**(0.5))

#for i in range(3):
#    axarr[i].axvline(x= 8.655,color='k',ls='-', label='Contact Surface')
#    axarr[i].axvline(x= -41.75,color='k',ls='--',label='Left Rarefaction Head')
#    axarr[i].axvline(x= -31.35,color='k',ls='-.',label='Left Rarefaction Tail')
#    axarr[i].axvline(x= 35.35,color='k',ls=':', label='Right Shockwave')
axarr[0].legend(loc=2, fontsize = 'xx-small')

#axarr[0].plot(xe100, re100, '-^c')
#axarr[0].set_ylabel(r'\rho')
#
#axarr[1].plot(xe100, ue100, '-^c')
#axarr[1].set_ylabel(r'u')
#
#axarr[2].plot(xe100, pe100, '-^c')
#axarr[2].set_ylabel(r'p')
#axarr[2].set_xlabel(r'x')

plt.tight_layout()
pp.savefig(bbx_inches = 'tight')
pp.close()

######################

nt = 16
cmap = plt.get_cmap('jet')
colors = cmap(np.linspace(0, 1.0, nt))

fname = ['sol']
pp = PdfPages('task21.pdf')
fig, axarr = plt.subplots(3, sharex=True,figsize=(8, 8))
fname = ['']*nt
print fname
for i in range(nt):
    fname[i] = 'time%d'%(i+1)
    
print fname
for caseid in range(nt):
    inname = fname[caseid] + '.dat'

    x, r, u, p, tf = np.loadtxt(inname, dtype = np.float64, skiprows=1, unpack=True)

#   axarr[0].set_title(fname[caseid])

    axarr[0].plot(x, r, '-',color=colors[caseid],mfc='none',ms=3, label='t = %f'%tf[0])
    axarr[0].set_ylabel(r'\rho')

    axarr[1].plot(x, u, '-',color=colors[caseid],mfc='none',ms=3)
    axarr[1].set_ylabel(r'u')

    axarr[2].plot(x, p, '-',color=colors[caseid],mfc='none',ms=3)
    axarr[2].set_ylabel(r'p')
    axarr[2].set_xlabel(r'x')

#   axarr[0].set_ylim([0, 30])
#   axarr[1].set_ylim([-10, 10])
#   axarr[2].set_ylim([0, 400])
    axarr[0].set_xlim([0, 1.6])
    axarr[1].set_xlim([0, 1.6])
    axarr[2].set_xlim([0, 1.6])

axarr[0].legend(loc=2, ncol=2, fontsize = 'xx-small')

plt.tight_layout()
pp.savefig(bbx_inches = 'tight')
pp.close()

######################
nt = 6

pp = PdfPages('task22.pdf')
fig = plt.figure(figsize=(8, 4))
fname = ['']*nt
print fname
for i in range(nt):
    fname[i] = 'size%d'%(i+1)
    
print fname
nx = 75
for caseid in range(nt):
    inname = fname[caseid] + '.dat'

    t, p, x = np.loadtxt(inname, dtype = np.float64, skiprows=1, unpack=True)

#   axarr[0].set_title(fname[caseid])
    nx = nx*2

    plt.plot(t, p, '-',mfc='none',ms=3, label='size = %d'%nx)
    plt.ylabel(r'p')
    plt.xlabel(r't')


plt.legend(loc=2, ncol=2, fontsize = 'xx-small')

plt.tight_layout()
pp.savefig(bbx_inches = 'tight')
pp.close()

pp = PdfPages('task22_peak.pdf')
fig = plt.figure(figsize=(8, 4))
fname = ['']*nt
print fname
for i in range(nt):
    fname[i] = 'size%d'%(i+1)
    
gridsize = np.empty(nt,dtype=np.float64)
maxp = np.empty(nt,dtype=np.float64)

nx = 75
for caseid in range(nt):
    inname = fname[caseid] + '.dat'

    t, p, x = np.loadtxt(inname, dtype = np.float64, skiprows=1, unpack=True)
    nx = nx*2
    dx = 1.5/nx
    gridsize[caseid] = nx
    maxp[caseid] = max(p)
print gridsize
plt.plot(gridsize, maxp, '-o',mfc='none',ms=3, label='size = %d'%caseid)
plt.xlabel(r'nx')
plt.ylabel(r'max($p(x=1.3,t)$)')

plt.tight_layout()
pp.savefig(bbx_inches = 'tight')
pp.close()
