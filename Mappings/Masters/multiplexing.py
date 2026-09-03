import matplotlib.pyplot as plt
import numpy as np
from classes import Mapping, linear, quadratic
from scipy.optimize import curve_fit
from copy import copy
import pickle
import functions as ft
from sklearn.neighbors import KernelDensity
import pandas as pd
import os

import matplotlibcolors_v2
plt.style.use('Mappings/Masters/matplotlibrc_v2')


oct = 1
Q = 50e3
chi = 4
ticks = np.array([1e-5, 1e-4, 1e-3, 1e-2])
# yld = .97

nr = 100
Ns = [1000, 2000, 4000]
sigmas = np.logspace(-5, -2, nr)
p0 = .95
cs = 'bop'
fig, ax = plt.subplots(constrained_layout=True, figsize=(18.5/2/2.54, 8/2.54))
for i, N in enumerate(Ns):
    Delta = 2**(oct/(N-1)) - 1
    print(Delta*Q)
    ylds = []
    for sigma in sigmas:
        ylds.append(ft.p0(Q, chi, sigma, Delta))
    idx = np.argmin(np.abs(np.array(ylds)-p0))
    ax.semilogx(sigmas, ylds, label='N=%d' % (N), c=cs[i], lw=2)
    ax.annotate('$\lambda_D=%.f dF$, $N=%d$' % (Delta*Q, N), xy=(4e-4, ylds[-1]-0.01), color=cs[i], ha='left', va='top')
    ax.axvline(sigmas[idx], c=cs[i], ls='--', label=f'$\sigma=%.1f\\times 10^{-4}$' % (sigmas[idx]*1e4))
    n = -3
    while sigmas[idx] < 10**n:
        n -= 1
    else:
        ax.annotate('$\sigma=%.1f \\times 10^{%d}$' % (sigmas[idx]*10**(-n), n), 
                    xy=(.9*sigmas[idx], i*.1), 
                    color=cs[i], 
                    ha='right',
                    va='bottom',
                    rotation=90,
                    bbox=dict(facecolor='none', alpha=0.5, edgecolor='none'))

    ax.set_ylim(0,1)
    ax.set_xlim(1e-5, 1e-2)
# ax.legend(loc='lower right', ncols=3, handlelength=1, columnspacing=0.5)
# ax.axvline(1e-3, c='k', ls='--')

ax.axhline(p0, c='k', ls='--', lw=2, zorder=0)
ax.annotate('$P_0=%d\%%$' % (p0*100), xy=(2e-5, p0-.03), color='k', ha='center', va='top')
# ax.set_xticks(ticks)
# ax.set_xticklabels(ticks)
# ax.set_title('$Q=%d,\chi=%d$' % (Q,chi))
ax.set_xlabel('$\sigma$ (-)')
ax.set_ylabel('$P_0$ (-)')
# Create a secondary x-axis
ax.grid(False, which='minor')
def toQ(x):
    return x*Q

def fromQ(x):
    return x/Q

ax2 = ax.secondary_xaxis('top', functions=(toQ, fromQ))
ax2.set_xticks(ticks*Q)
ax2.set_xticklabels([0.5, 5, 50, 500])
ax2.set_xlabel('$\sigma\\times Q_l$ ($\mathrm{MHz}$)')\

plt.savefig('Mappings/Masters/figures/required_scatter.pdf')
plt.savefig('Mappings/Masters/figures/required_scatter.svg', transparent=True)
plt.show()
