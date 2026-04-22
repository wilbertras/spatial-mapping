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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import matplotlibcolors_v2
plt.style.use('Mappings/Masters/matplotlibrc_v2')

Q=50e3
ylds_at = 4
dir = 'Mappings/Masters/chips/'
LT402chip6 = Mapping(dir + 'LT402chip6_master.pkl')
LT402chip6_trim = Mapping(dir + 'LT402chip6_master.pkl', type='trim')
f, s21 = np.load(dir + 'LT402chip6_S21_dark.npy')
f_trim, s21_trim = np.load(dir + 'LT402chip6_S21_dark_trim.npy')
dflim = [0, 100]
dfticks = np.arange(0, 101, 25)
df_flim = [0, 100/Q]
df_fticks = dfticks/Q
bins = np.linspace(dflim[0], dflim[1], 40)
chips = [LT402chip6, LT402chip6_trim]
labels = ['C', 'C, trimmed']
cs = 'bp'
slim = [-40, 5]
inset_slim = [-35, 0]
flim = [4, 8.2]
inset_flimt = [4.7, 4.9]
fmt = LT402chip6_trim.fm
fdt  = LT402chip6_trim.fd
ldt = np.nanmean(np.diff(fdt[:500])/fdt[:499]*Q)
ldt2 = np.nanmean(np.diff(fdt[600:1000])/fdt[600:999]*Q)
print(ldt, ldt2)
fm = LT402chip6.fm
kid_start = np.nanargmin(np.abs(fmt-inset_flimt[0]))
kid_end = np.nanargmin(np.abs(fmt-inset_flimt[1]))
inset_flim = [fm[kid_start], fm[kid_end]]
pos = 'adbce'
fig, axes = plt.subplot_mosaic('aaabbb;ccddee', figsize=(18.5/2.54,12/2.54), constrained_layout=True)
for i, a in enumerate(pos):
    if a in 'cde':
        x = -0.3
    else:
        x = -0.1
    axes[a].annotate('(%s)' % a, xy=(x, 1), xycoords='axes fraction', ha='left', va='top', fontsize='small')
ax = axes['a']
f_interp = np.linspace(f[0], f[-1], len(f), endpoint=True)
s21_fit = np.interp(f_interp, [f[0], f[-1]], [s21[0], s21[-1]])
ax.plot(f, s21-s21_fit, label='S21', lw=.01, color=cs[0], zorder=0)
ax.set_xlim(flim)
ax.set_ylim(slim)
ax.set_ylabel('$|S_{21}|^2$')
ax.set_xlabel('$f$ $(\mathrm{GHz})$')
ax.set_title('Before trimming')
ax_inset = inset_axes(ax, width="40%", height="50%", loc='lower right')
ax_inset.plot(f, s21-s21_fit, label='S21', lw=.2, color=cs[0])
ax_inset.set_xlim(inset_flim)
ax_inset.set_ylim(slim)
ax_inset.set_xticks([])
ax_inset.set_yticks([])
ax.add_patch(plt.Rectangle((inset_flim[0], inset_slim[0]), inset_flim[1] - inset_flim[0], inset_slim[1] - inset_slim[0], 
                                   ec='k', fc='None', alpha=1, transform=ax.transData, lw=1))  # Changed alpha to 0.5
ax = axes['b']
f_interp = np.linspace(f_trim[0], f_trim[-1], len(f_trim), endpoint=True)
s21_fit = np.interp(f_interp, [f_trim[0], f_trim[-1]], [s21_trim[0], s21_trim[-1]])
ax.plot(f_trim, s21_trim-s21_fit, label='S21', c=cs[1], lw=.01, zorder=0)
ax.set_ylim(slim)
ax.set_xlim(flim)
ax.set_xlabel('$f$ $(\mathrm{GHz})$')
ax.set_ylabel('$|S_{21}|^2$')
ax.set_title('After trimming')
ax_inset = inset_axes(ax, width="40%", height="50%", loc='lower right')
ax_inset.plot(f, s21_trim-s21_fit, label='S21', lw=.2, color=cs[1])
ax_inset.set_xlim(inset_flimt)
ax_inset.set_ylim(slim)
ax_inset.set_xticks([])
ax_inset.set_yticks([])
ax.add_patch(plt.Rectangle((inset_flimt[0], inset_slim[0]), inset_flimt[1] - inset_flimt[0], inset_slim[1] - inset_slim[0], 
                                   ec='k', fc='None', alpha=1, transform=ax.transData, lw=1))  # Changed alpha to 0.5

ax = axes['c']
_ = ax.hist(chips[0].rel_spacings, bins, facecolor=cs[0], alpha=.75)
ax.set_xlim(dflim)
ax.set_xticks(dfticks)
ax.axvline(ylds_at, c='r', ls='-', lw=1, label='Min. spacing', zorder=5)
ax.set_ylabel('$\mathrm{Counts}$')   
ax.set_xlabel('$\lambda_M \\times Q$ $(\mathrm{lw})$')
ax2 = ax.secondary_xaxis('top')
ax2.set_xticks(dfticks)
ax2.set_xticklabels(df_fticks)
ax2.set_xlabel('$\lambda_M$')

ax = axes['e']
lambdas, _, _ = chips[1].spacings(Q, ylds_at)
lo = np.nanargmax(lambdas)
print(lo, lambdas[lo])
_ = ax.hist(lambdas[:lo], bins, facecolor=cs[1], edgecolor='None', alpha=.75, label='lower band')
_ = ax.hist(lambdas[lo:], bins, facecolor='None', edgecolor=cs[1], alpha=.75, label='upper band')
ax.axvline(ylds_at, c='r', ls='-', lw=1, zorder=5)
_ = ax.set_xlim(dflim)
# _ = ax.set_ylim([0, 120])
ax.set_xlabel('$\lambda_M^\mathrm{trim} \\times Q$ $(\mathrm{lw})$')
ax.set_ylabel('$\mathrm{Counts}$')   
ax2 = ax.secondary_xaxis('top')
ax2.set_xticks(dfticks)
ax2.set_xticklabels(df_fticks)
ax2.set_xlabel('$\lambda_M$')
sigma = np.nanstd(chips[1].df_f)
ft.p0_numeric(1, 50e3, 1, sigma, ylds_at, nr_iter=100, fs=chips[1].fd)
x, kernel, fwhm = ft.kde(lambdas[:lo], bw=1)
std = fwhm/Q/2.355
ax.plot(x, kernel*len(lambdas[:lo])*(bins[1]-bins[0]), c='k', ls='--', lw=1, label='KDE')
ax.annotate('$\sigma_\mathrm{KDE}^\mathrm{trim}=%.2G\\times 10^{-5}$' % (std*1e5), xy=(0.5, .95), xycoords='axes fraction', ha='center', va='top', color='k')
# ax.axvline(ldt)
# ax.axvline(ldt2)
ax.set_ylim(0, 125)
ax.legend(loc='center right', handlelength=1)

ax = axes['d']
thresholds = np.arange(1, 11)
ylds_before = []
ylds_after = []
for threshold in thresholds:
    ylds_before.append(chips[0].spacings(Q, threshold)[2])
    ylds_after.append(chips[1].spacings(Q, threshold)[2])
ax.plot(thresholds, np.array(ylds_before), lw=1, label='Before trimming', c=cs[0])
ax.plot(thresholds, np.array(ylds_after), lw=1, label='After trimming', c=cs[1])	
yld_before = ylds_before[ylds_at-1]
yld_after = ylds_after[ylds_at-1]
ax.axvline(ylds_at, c='r', ls='-', lw=1, zorder=0, label='$\lambda_M^\mathrm{min}$')
ax.scatter(ylds_at, yld_before, c=cs[0], s=10)
ax.scatter(ylds_at, yld_after, c=cs[1], s=10)
ax.annotate('%d\%%' % (yld_before*1e2), xy=(ylds_at, yld_before), xytext=(-10, -15), textcoords='offset points', ha='center', va='bottom', color=cs[0])
ax.annotate('%d\%%' % (yld_after*1e2), xy=(ylds_at, yld_after), xytext=(-10, -15), textcoords='offset points', ha='center', va='bottom', color=cs[1])
print('yields: %.2f %% and %.2f %% ' % (yld_before, yld_after))
lim = [0, 10]
ax.set_xlim(lim)
ticks = np.arange(0, 11, 5)
ax.set_xticks(ticks)
ax.set_xticklabels(ticks)
tickss = np.arange(0, 11, 1)
ax.set_xticks(tickss, minor=True)
ax.set_ylim(0.5, 1)
ax2 = ax.secondary_xaxis('top')
ax2.set_xticks(ticks)
ax2.set_xticklabels(ticks/Q)
ax2.set_xticks(tickss, minor=True)
ax2.set_xlabel('$\lambda_M^\mathrm{min}$')
ax.set_xlabel('$\lambda_M^\mathrm{min} \\times Q$ $(\mathrm{lw})$')
ax.set_ylabel('$P_0$')
# ax.legend()


# ax.legend(loc='upper right')
plt.savefig('Mappings/Masters/figures/yield.pdf')
plt.show()