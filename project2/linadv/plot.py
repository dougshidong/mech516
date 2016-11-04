#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

caseid = 1
exercise = 3

fname = 'exercise1%d.pdf' %exercise
pp = PdfPages(fname)

symbols = ['o','s','^','>','<']
xe = np.linspace(-10,240,25000)
u1 = [0.5 if xi <= 200 or xi>=220 else 0.5 + 0.075*(xi-200) for xi in xe]
u2 = [0.5 if xi <= 200 or xi>=220 else 1.0 - 0.5 * np.cos(np.pi*(xi-200)/10.0) for xi in xe]

if exercise == 1:
    cfl = [0.8, 1.02, 2.0]
    for i in range(3):
        plt.figure(figsize=(8,3))
        for t in range(5):
            fname = 'p1cfl%dt%d.dat' % (i, t)
            x, u = np.loadtxt(fname, dtype = np.float64, unpack=True)
            plt.plot(x,u,symbols[t],label='Step = %d' % (t+1),mec='none')
        plt.legend(fontsize='x-small')
        plt.xlim([-5,25])
        plt.ylabel(r'$u$')
        plt.xlabel(r'$x$')
        plt.tight_layout()
        pp.savefig(bbx_inches = 'tight')
    pp.close()
if exercise == 2:
    nxl = [1000, 2000, 4000]
    for prob in range(2):
        for scheme in range(2):
            plt.figure(figsize=(8,6))
            for inx, nx in enumerate(nxl):
                fname = 'p2case%df%dnx%d.dat' % (prob+1, scheme+1,nx)
                x, u = np.loadtxt(fname, dtype = np.float64, unpack=True)
                plt.plot(x,u,symbols[inx],label='\Delta x = %6.4f' % (250.0/nx), ms=3,mec='none')
            if prob == 0:
                plt.plot(xe,u1,ls = '-',label='Exact',color = 'k')
                plt.ylim([0.4,2.5])
            if prob == 1:
                plt.plot(xe,u2,ls = '-',label='Exact',color = 'k')
                plt.ylim([0.4,1.6])
            plt.legend(loc=2)
            plt.xlim([195,225])
            plt.ylabel(r'$u$')
            plt.xlabel(r'$x$')
            plt.tight_layout()
            pp.savefig(bbx_inches = 'tight')
    pp.close()
if exercise == 3:
    nxl = np.array([1000, 2000, 4000, 8000])
    scheme_name = ['CIR', 'LW']
    for prob in range(2):
        plt.figure(figsize=(6,4.5))
        for scheme in range(2):
            error = np.empty(len(nxl),dtype=np.float64)
            for inx, nx in enumerate(nxl):
                fname = 'p2case%df%dnx%d.dat' % (prob+1, scheme+1,nx)
                x, u = np.loadtxt(fname, dtype = np.float64, unpack=True)
                if prob == 0:
                    ue = [0.5 if xi <= 200 or xi>=220 else 0.5 + 0.075*(xi-200) for xi in x]
                if prob == 1:
                    ue = [0.5 if xi <= 200 or xi>=220 else 1.0 - 0.5 * np.cos(np.pi*(xi-200)/10.0) for xi in x]
                error[inx] = np.linalg.norm(ue-u,ord=2) / nx
            print nxl
            print error
            slope, intercept = np.polyfit(np.log(nxl/250.0), np.log(error), 1)
            print(slope)
            plt.loglog(nxl/250.0,error,symbols[scheme],label=scheme_name[scheme],ls='-')
            print 'test'
        plt.legend(loc=1)
        plt.xlabel(r'$1/\Delta x$')
        plt.ylabel(r'$\|u_{numerical}-u_{exact}\|_2/nx$')
        plt.tight_layout()
        pp.savefig(bbx_inches = 'tight')
    pp.close()
