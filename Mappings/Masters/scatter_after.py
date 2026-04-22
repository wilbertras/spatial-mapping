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
import matplotlibcolors_v2
from matplotlib.colors import LinearSegmentedColormap
plt.style.use('Mappings/Masters/matplotlibrc_v2')

# Q = 50e3
# chi = 2
# dir = 'Mappings/Masters/chips/'
# LT402chip6 = Mapping(dir + 'LT402chip6_master.pkl', Q=Q, min_lw_spacing=chi)
# LT402chip6_trim = Mapping(dir + 'LT402chip6_master.pkl', type='trim', Q=Q, min_lw_spacing=chi)
# f, s21 = np.load(dir + 'LT402chip6_S21_dark.npy')
# f_trim, s21_trim = np.load(dir + 'LT402chip6_S21_dark_trim.npy')

# chips = [LT402chip6, LT402chip6_trim]
# labels = ['C', 'C, trimmed']
# fig, axes = plt.subplot_mosaic('abcc;decc', constrained_layout=True, figsize=(18.5/2.54, 8/2.54))
# ylim = -20, 80
# clim = -2.5, .5
# cticks = [-2, -1, 0]
# flim = 4.5, 5
# binsize = 2.5
# s = 5
# bins=np.arange(ylim[0], ylim[1], binsize)
# cs = 'bp'
# from matplotlib.colors import LinearSegmentedColormap
# cmap = LinearSegmentedColormap.from_list('custom_green', ['w', cs[-1]])
# for i, chip in enumerate(chips):
#     ax = axes['a']
#     ax.scatter(chips[0].fd, chip.df_f*1e3, alpha=.5, label='Device %s' % labels[i], s=s, c=cs[i], edgecolor='None')
#     ax.set_ylim(ylim)
#     # ax.set_xlim(flim)
#     ax.set_ylabel('$\epsilon \\times 10^{3}$')
#     ax = axes['b']
#     ax.hist(chip.df_f*1e3, bins=bins, alpha=.75, label='Device %s' % labels[i], zorder=len(chips)-i, orientation='horizontal', facecolor=cs[i], edgecolor='None')
#     ax.set_ylim(ylim)
#     ax.set_xscale('linear')
#     ax.set_yticklabels([])
#     if i==1:
#         binsize = 0.1
#         bins = np.arange(clim[0], clim[1], binsize)
#         ax = axes['d']
#         ax.scatter(chips[0].fd, chip.df_f*1e3, alpha=.5, label='Device %s' % labels[i], color=cs[i], s=s, zorder=1)
#         secax = ax.twiny()
#         secax.scatter(chips[0].fd, chip.df_f*1e3, color='None')
#         secax.set_xlim(ax.get_xlim())  # Set limits to match the original x-axis
#         Lf = chip.Lf[~np.isnan(chip.df_f_fit)]
#         fingers = np.linspace(10*1, 10*97, 11, endpoint=True)
#         for k, finger in enumerate(fingers):
#             id = np.argmin(np.abs(Lf - finger))
#             fd_fingers = []
#             if id:
#                 fd_finger = chips[0].fd[id]
#                 secax.axvline(fd_finger, c='k', ls='--', alpha=.75, lw=.5, zorder=0)
#                 fd_fingers.append(fd_fingers)
#         secax.set_xticks(fd_fingers)
#         secax.set_xticklabels(np.arange(len(fd_fingers)))
#         ax.set_ylim(clim)
#         ax.set_yticks(cticks)
#         ax.set_ylabel('$\epsilon \\times 10^{3}$')
#         ax.set_xlabel('$f$ $(\mathrm{GHz})$')
#         ax = axes['e']
#         ax.hist(chip.df_f*1e3, bins=bins, alpha=.75, density=False, label='Device %s' % labels[i], zorder=len(chips)-i, orientation='horizontal', color=cs[i])
#         ax.set_xscale('linear')
#         ax.set_ylim(clim)
#         ax.set_yticks(cticks)
#         ax.set_yticklabels([])
#         ax.set_xlabel('$\mathrm{Counts}$')
#         ax = axes['c']
#         df_f = chip.df_f_fit*1e3
#         # ax.imshow(df_f[chip.map], origin='lower', cmap='purples', clim=clim)
#         cbar = fig.colorbar(ax.imshow(df_f[chip.map], cmap=cmap, origin='lower', clim=clim), ax=ax, orientation='vertical', location='right')
#         cbar.set_label('$\epsilon \\times 10^{3}$')
#         cbar.set_ticks(cticks)
#         ax.set_xticks([])
#         ax.set_yticks([])
#         ax.set_yticklabels([])
#         ax.set_xticklabels([])
#         ax.set_xlabel('$x$')
#         ax.set_ylabel('$y$')

# axes['b'].annotate('$\sigma_M=%.1f\\times 10^{-3}$' % (np.nanstd(chips[0].df_f)*1e3), xy=(0.5, 0.9), xycoords='axes fraction', ha='center', va='top', color=cs[0], fontsize='small')
# axes['e'].annotate('$\sigma_M=%.1f\\times 10^{-3}$' % (np.nanstd(chips[1].df_f)*1e3), xy=(0.5, 0.9), xycoords='axes fraction', ha='center', va='top', color=cs[1], fontsize='small')
# plt.show()

Q = 50e3
chi = 2
dir = 'Mappings/Masters/chips/'
LT402chip6 = Mapping(dir + 'LT402chip6_master.pkl', Q=Q, min_lw_spacing=chi)
LT402chip6_trim = Mapping(dir + 'LT402chip6_master.pkl', type='trim', Q=Q, min_lw_spacing=chi)
f, s21 = np.load(dir + 'LT402chip6_S21_dark.npy')
f_trim, s21_trim = np.load(dir + 'LT402chip6_S21_dark_trim.npy')

chips = [LT402chip6, LT402chip6_trim]
labels = ['C', 'C, trimmed']
fig, axes = plt.subplot_mosaic('dce', constrained_layout=True, figsize=(18.5/2.54, 8/2.54))
ylim = -20, 80
clim = -2.5, .5
cticks = [-2, -1, 0]
flim = 4.5, 5
binsize = 2.5
s = 5
bins=np.arange(ylim[0], ylim[1], binsize)
cs = 'bp'
from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list('custom_green', ['w', cs[-1]])
for i, chip in enumerate(chips):
    if i==1:
        binsize = 0.1
        bins = np.arange(clim[0], clim[1], binsize)
        ax = axes['d']
        ax.scatter(chips[0].fd, chip.df_f*1e3, alpha=.5, label='Device %s' % labels[i], color=cs[i], s=s, zorder=1)
        secax = ax.twiny()
        secax.scatter(chips[0].fd, chip.df_f*1e3, color='None')
        secax.set_xlim(ax.get_xlim())  # Set limits to match the original x-axis
        Lf = chip.Lf[~np.isnan(chip.df_f_fit)]
        fingers = np.linspace(10*1, 10*97, 11, endpoint=True)
        for k, finger in enumerate(fingers):
            id = np.argmin(np.abs(Lf - finger))
            fd_fingers = []
            if id:
                fd_finger = chips[0].fd[id]
                secax.axvline(fd_finger, c='k', ls='--', alpha=.75, lw=.5, zorder=0)
                fd_fingers.append(fd_fingers)
        secax.set_xticks(fd_fingers)
        secax.set_xticklabels(np.arange(len(fd_fingers)))
        ax.set_ylim(clim)
        ax.set_yticks(cticks)
        ax.set_ylabel('$\epsilon \\times 10^{3}$')
        ax.set_xlabel('$f$ $(\mathrm{GHz})$')
        ax = axes['e']
        ax.hist(chip.df_f*1e3, bins=bins, alpha=.75, density=False, label='Device %s' % labels[i], zorder=len(chips)-i, orientation='vertical', color=cs[i])
        ax.set_xscale('linear')
        ax.set_xlim(clim)
        ax.set_yticks(cticks)
        ax.set_yticklabels([])
        ax.set_xlabel('$\mathrm{Counts}$')
        ax = axes['c']
        df_f = chip.df_f_fit*1e3
        # ax.imshow(df_f[chip.map], origin='lower', cmap='purples', clim=clim)
        cbar = fig.colorbar(ax.imshow(df_f[chip.map], cmap=cmap, origin='lower', clim=clim), ax=ax, orientation='horizontal', location='bottom')
        cbar.set_label('$\epsilon \\times 10^{3}$')
        cbar.set_ticks(cticks)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')

# axes['b'].annotate('$\sigma_M=%.1f\\times 10^{-3}$' % (np.nanstd(chips[0].df_f)*1e3), xy=(0.5, 0.9), xycoords='axes fraction', ha='center', va='top', color=cs[0], fontsize='small')
axes['e'].annotate('$\sigma_M=%.1f\\times 10^{-3}$' % (np.nanstd(chips[1].df_f)*1e3), xy=(0.3, 0.9), xycoords='axes fraction', ha='center', va='top', color=cs[1], fontsize='small')
plt.show()
# plt.savefig('after_trimming.png')