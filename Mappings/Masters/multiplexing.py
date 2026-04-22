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
req_yld = .95
nr = 200
chis = np.logspace(-5, -3, nr)
Ns = np.linspace(2, 10000, nr)
Deltas = 2**(oct/(Ns-1)) - 1
sigmas = np.logspace(-5, -3, nr)
redo=0
file = r'Mappings/Masters/req_Ns.npy'
if os.path.exists(file) and not redo:
    req_Ns = np.load(file)
else:
    req_Ns = np.zeros((len(sigmas), len(chis)))
    for i, sigma in enumerate(sigmas):
        for j, chi in enumerate(chis):
            ylds = []
            for Delta in Deltas:
                ylds.append(ft.p1(chi, sigma, Delta))
            ylds = np.array(ylds)
            id = np.argmin(np.abs(ylds - req_yld))
            req_Ns[j, i] = Ns[id]
    np.save(file, req_Ns)
    print('Calculated req_Ns and saved to file.')


Q = 50e3
fig, axes = plt.subplot_mosaic('a', constrained_layout=True, figsize=(18.5/2/2.54, 8/2.54), sharex=False)
# ax.pcolormesh(Ns, sigmas, req_chi)
ax = axes['a']
ax.set_yscale('log')
ax.set_xscale('log')
c = ax.pcolormesh(sigmas, chis, req_Ns, shading='auto', cmap='viridis', norm='linear')
# cbar = fig.colorbar(c, ax=ax, pad=0)
# cbar.set_label('$\mathrm{MUX}$')
ax2 = ax.secondary_yaxis('right')
ax2.set_yscale('log')
labels = [1, 2, 5, 10, 20, 50]
ax2.set_yticks(np.array(labels)/Q)
ax2.set_yticklabels(labels)

ax3 = ax.secondary_xaxis('top')
ax3.set_xscale('log')
ax3.set_xticks(np.array(labels)/Q)
ax3.set_xticklabels(labels)

contour_levels = [500, 1000, 2000, 4000, 6000, 8000]
CS = ax.contour(sigmas, chis, req_Ns, levels=contour_levels, colors='w')
ax.clabel(CS, inline=True, fontsize=8)
ax.set_aspect('equal')
ax.set_xlim(chis[0], chis[-1])
ax.set_ylim(sigmas[0], sigmas[-1])
ax.set_xlabel('$\sigma$')  
ax.set_ylabel('$\lambda_M^\mathrm{min}$')  
ax2.set_ylabel('$\lambda_\mathrm{min}\\times Q$ $(\mathrm{lw})$')
ax3.set_xlabel('$\sigma\\times Q$ $(\mathrm{lw})$')
ax.annotate('$\mathrm{MUX-factor}$', xy=(0.5, 0.9), xycoords='axes fraction', color='w', fontsize='small', ha='center')
# ax.axvline(2/Q, color='o', linestyle='--')
# ax.plot(chis, sigmas, c='o', ls='--')
plt.savefig('Mappings/Masters/figures/multiplexing_map.pdf')


fig, ax = plt.subplots(figsize=(18.5/2/2.54, 8/2.54), constrained_layout=True)
chi = 4e-5
for chi in [4e-5, 8e-5, 16e-5]:
    id = np.argmin(np.abs(chis - chi))
    lw = chi*Q
    ax.semilogx(sigmas, req_Ns[id, :], label='$\lambda_\mathrm{min}=%.2G\\times10^{-4}(=%d\ \mathrm{lw}/Q)$' % (chi*1e4, lw))
ax.legend()
ax.set_xlabel('$\sigma$')
ax.set_ylabel('$\mathrm{MUX-factor}$')
ax.set_xlim(sigmas[0], sigmas[-1])
ax.set_ylim(0, 10000)
plt.savefig('Mappings/Masters/figures/multiplexing.pdf')
plt.show()