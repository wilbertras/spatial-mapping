import matplotlib.pyplot as plt
import numpy as np
from classes import Mapping, linear, quadratic
from scipy.optimize import curve_fit
from copy import copy
import pickle
import functions as ft
from sklearn.neighbors import KernelDensity
import pandas as pd
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.legend_handler import HandlerTuple
import matplotlibcolors_v2
plt.style.use('Mappings/Masters/matplotlibrc_v2')

Q = 50e3
chi = 2
dir = r'Mappings/Masters/chips/'
offset = True
deg = 2
LT402chip6 = Mapping(dir + 'LT402chip6_master.pkl', Q=Q, min_lw_spacing=chi, deg=deg, offset=offset, mask_edges=True)
LT402chip6_trim = Mapping(dir + 'LT402chip6_master.pkl', type='trim', Q=Q, min_lw_spacing=chi, deg=deg, offset=offset, mask_edges=True)

chips = [LT402chip6]
labels = ['A', 'B']
cs = 'bb'
cs2 = 'oo'
cmap = cm.get_cmap('Oranges', 256)
cmap = LinearSegmentedColormap.from_list("custom_orange", cmap(np.linspace(.1, 1, 256)))
flim = [3.5, 8.5]
dflim = [-20, 80]
dfticks = [-20, 0, 20, 40, 60, 80]
# clims = [[-20, 5], [-2, 2], [-25, 25]]
clim = [-20, 30]
ylims = [[0, 150]]
bins = np.linspace(-20, 80, 40)
cticks = [-20, -10, 0, 10, 20, 30]
pos = 'abcdef'
fig, axes = plt.subplot_mosaic('lll;abc;def', figsize=(18.5/2.54, 18.5/1.75/2.54), gridspec_kw={"height_ratios": [0.025, 1, 1]}, constrained_layout=True)
for i, a in enumerate(pos):
    if i==2 or i==5:
        x = -0.2
    else:
        x = -0.3
    axes[a].annotate('(%s)' % a, xy=(x, 1), xycoords='axes fraction', ha='left', va='top', fontsize='small')
for i, chip in enumerate(chips):
    if offset:
        fit = np.polyval(chip.coeff, chip.fd)
    else:
        if deg==1:
            fit = linear(chip.fd, *chip.coeff)
        elif deg==2:
            fit = quadratic(chip.fd, *chip.coeff)
    ax = axes[pos[i*3]]
    df_f_fit = (fit - chip.fd)/chip.fd
    ax.scatter(chip.fd, chip.df_f*1e3, color=cs[i], label='$\epsilon$', alpha=.75, s=2)
    ax.plot(chip.fd, df_f_fit*1e3, linestyle='--', color='k', label='$\mathrm{fit}$')
    # ax.scatter(chip.fd, chip.df_f_fit*1e3, c=chip.df_f_fit*1e3, vmin=clim[0], vmax=clim[1], label='$\epsilon-\mathrm{fit}$', alpha=.75, s=2)
    ax.scatter(chip.fd, chip.df_f_fit*1e3, fc=cs2[i], label='$\epsilon-\mathrm{fit}$', alpha=.75, s=2)
    # ax.plot(chip.fd, np.zeros(chip.fd.shape), linestyle='-', color='k', label='data')
    ax.set_xlim(flim)
    ax.set_ylim(dflim)
    ax.set_yticks(dfticks)
    ax.set_ylabel('$\epsilon\\times 10^{3}$')
    ax.set_xlabel('$F_D$ $(\mathrm{GHz})$')

    ax.annotate('$\epsilon_M$', xy=(4.2, 60), ha='center', va='center', color='b', fontsize='small')
    ax.annotate('$\mathrm{fit}$', xy=(8, 50), ha='center', va='center', color='k', fontsize='small')
    ax.annotate('$\epsilon$', xy=(8, 20), ha='center', va='center', color=cs2[i], fontsize='small')
    ax.annotate('$\epsilon^\mathrm{trim}$', xy=(8, -10), ha='center', va='center', color='p', fontsize='small')

    # ax.legend(loc='upper right', ncols=3, handlelength=1, frameon=False, fontsize='small')
    ax = axes[pos[2+i*3]]
    empty_array = chip.df_f_fit[chip.map]*1e3
    y_indices, x_indices = np.where(np.isnan(empty_array))
    if i:
        ax.scatter(x_indices, y_indices, color='k', marker='x', s=10, label='NaN', alpha=1, linewidth=.5)
    ax.imshow(empty_array, origin='lower', cmap=cmap, vmin=clim[0], vmax=clim[1])
    # empty_array = chip.df_f[chip.map]*1e3
    y_indices, x_indices = np.where(np.isnan(empty_array))
    ax.scatter(x_indices, y_indices, color='k', marker='x', s=10, label='NaN', alpha=1, linewidth=.5)
    # Set NaN values to a specific color (black)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    for spine in ax.spines.values():
        spine.set_edgecolor('white')

    ax = axes[pos[1+i*3]]
    ax.hist(chip.df_f*1e3, facecolor=cs[i], alpha=.75, label='data', orientation='horizontal', bins=bins)
    ax.hist(chip.df_f_fit*1e3, facecolor='o', alpha=.75, label='data', orientation='horizontal', bins=bins)
    ax.set_xlabel('$\mathrm{Counts}$')
    ax.annotate('$\epsilon_M$', xy=(100, 50), ha='center', va='center', color='b', fontsize='small')
    ax.annotate('$\sigma_M=%.1f\\times 10^{-2}$' % (np.nanstd(chip.df_f)*1e2), xy=(100, 40), ha='center', va='center', color=cs[i], fontsize='small')
    ax.annotate('$\epsilon$', xy=(100, 25), ha='center', va='center', color='o', fontsize='small')
    ax.annotate('$\sigma=%.1f\\times 10^{-2}$' % (np.nanstd(chip.df_f_fit)*1e2), xy=(100, 15), ha='center', va='center', color='o', fontsize='small')
    ax.annotate('$\epsilon^\mathrm{trim}$→', xy=(120, -10), ha='center', va='center', color='p', fontsize='small')
    ax.set_ylim(dflim)
    ax.set_yticks(dfticks)
    ax.set_xlim(ylims[i])
    ax.set_ylabel('$\epsilon\\times 10^{3}$')
cbar = fig.colorbar(axes['c'].imshow(empty_array, origin='lower', cmap=cmap, clim=clim), ax=axes['c'], orientation='vertical', location='right', pad=0.1)
cbar.set_label('$\epsilon\\times 10^{3}$')
cbar.set_ticks(cticks)


axes['l'].axis("off")
handle1 = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='b', markersize=8)
handle2 = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='o', markersize=8)
handle3 = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='p', markersize=8)
axes['l'].legend(
    [(handle1, handle2), handle3],
    ["Before trimming", "After trimming"],
    handler_map={tuple: HandlerTuple(ndivide=2)}, 
    loc='upper center', 
    handletextpad=0.5, frameon=True, ncols=2,
)

chip = LT402chip6_trim
c = 'p'
axes['a'].scatter(chip.fd, chip.df_f*1e3, facecolor=c, alpha=.75, s=2)
axes['a'].set_xlim(flim)
axes['b'].hist(chip.df_f*1e3, facecolor=c, bins=bins, alpha=1, zorder=3, orientation='horizontal')
ylim = 0, 200
clim = -2.5, .5
cticks = [-2, -1, 0]
binsize = 0.075
s = 5
bins=np.arange(ylim[0], ylim[1], binsize)
cmap = LinearSegmentedColormap.from_list('custom_green', ['w', c])
cmap = LinearSegmentedColormap.from_list("custom_orange", cmap(np.linspace(1, 0.2, 256)))
bins = np.arange(clim[0], clim[1], binsize)
ax = axes['d']
ax.scatter(chip.fd, chip.df_f*1e3, alpha=.5, color=c, s=s, zorder=1)
secax = ax.twiny()
secax.scatter(chip.fd, chip.df_f*1e3, color='None')
secax.set_xlim(ax.get_xlim())  # Set limits to match the original x-axis
Lf = chip.Lf[~np.isnan(chip.df_f)]
fingers = np.linspace(10, 10*97, 11, endpoint=True)
for k, finger in enumerate(fingers):
    id = np.argmin(np.abs(Lf - finger))
    fd_fingers = []
    if id:
        fd_finger = chip.fd[id]
        secax.axvline(fd_finger, c='k', ls='--', alpha=.75, lw=.5, zorder=0)
        fd_fingers.append(fd_fingers)
secax.set_xticks(fd_fingers)
secax.set_xticklabels(np.arange(len(fd_fingers)))
ax.set_ylim(clim)
ax.set_xlim(flim)
ax.set_yticks(cticks)
ax.set_ylabel('$\epsilon^\mathrm{trim}\\times 10^{3}$')
ax.set_xlabel('$F_D^\mathrm{trim}$ $(\mathrm{GHz})$')
ax = axes['f']
empty_array = chip.df_f[chip.map]*1e3
y_indices, x_indices = np.where(np.isnan(empty_array))
ax.scatter(x_indices, y_indices, color='k', marker='x', s=10, label='NaN', alpha=1, linewidth=.5)
cbar = fig.colorbar(ax.imshow(empty_array, cmap=cmap, origin='lower', clim=clim), ax=ax, orientation='vertical', location='right', pad=.1)
cbar.set_label('$\epsilon^\mathrm{trim}\\times 10^{3}$')
cbar.set_ticks(cticks)
ax.set_xticks([])
ax.set_yticks([])
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
for spine in ax.spines.values():
    spine.set_edgecolor('white')
ax = axes['e']
ax.hist(chip.df_f*1e3, bins=bins, alpha=.75, density=False, zorder=len(chips)-i, orientation='horizontal', color=c)
ax.annotate('$\epsilon^\mathrm{trim}$', xy=(0.55, 0.85), xycoords='axes fraction', ha='center', va='center', color='p', fontsize='small')
ax.annotate('$\sigma^\mathrm{trim}=%.1f\\times 10^{-4}$' % (np.nanstd(LT402chip6_trim.df_f)*1e4), xy=(0.55, 0.75), xycoords='axes fraction', ha='center', va='center', color='p', fontsize='small')
ax.set_ylim(clim)
ax.set_ylabel('$\epsilon^\mathrm{trim}\\times 10^{3}$')
ax.set_xlim(ylim)
ax.set_xlabel('$\mathrm{Counts}$')

plt.savefig('Mappings/Masters/figures/scatter_before_after.pdf')
plt.savefig('Mappings/Masters/figures/scatter_before_after.svg', transparent=True, bbox_inches='tight')
plt.show()