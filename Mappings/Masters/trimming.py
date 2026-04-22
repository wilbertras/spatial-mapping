import matplotlib.pyplot as plt
import numpy as np
from classes import Mapping, linear, quadratic
from scipy.optimize import curve_fit
from copy import copy
import pickle
import functions as ft
from sklearn.neighbors import KernelDensity
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlibcolors_v2
plt.style.use('Mappings/Masters/matplotlibrc_v2')

def poly(x, a, b, c, d, e, f):
    return a*x**5 + b*x**4 + c*x**3+d*x**2+e*x+f


dir = r'Mappings\Masters\chips\LT402chip6_master.pkl'
chip = Mapping(dir)
chiptrim = Mapping(dir, type='trim')

meas = chip.fm
design = chip.fd
trimdesign = chiptrim.fd
fingers = chip.Lf
newfingers = chiptrim.Lf
trims = newfingers - fingers
nanmask = np.isnan(meas)

popt, pcov = curve_fit(poly, meas[~nanmask], fingers[~nanmask])
fit_fingers = poly(meas, *popt)   

df_name = 'LT402_design_simulations'
dir = r'C:\Users\wilbertr\ownCloud2\PhD\PhD\PythonProjects\sonnet\LT402 Design/'
Lcs_at_Q_design = np.load(dir + df_name + '_Lc.npy')
f0s_at_Q_design = np.load(dir + df_name + '_f0.npy')*1e-9
df = pd.read_pickle(dir + '/' + df_name)
Lfs = np.unique(df['L_fingers'])


fig, axes = plt.subplot_mosaic('aabee;ccdee', figsize=(18.5/2.54, 18.5/2/2.54), constrained_layout=True)
pos = 'abcde'
for i, a in enumerate(pos):
    if a in 'bd':
        x = -0.4
    elif a in 'e':
        x = -0.1
    else:
        x = -0.2
    axes[a].annotate('(%s)' % a, xy=(x, 1), xycoords='axes fraction', ha='left', va='top', fontsize='small')
s = 2
ax = axes['a']
# ax.set_title('$F_M \\rightarrow F_D^\mathrm{trim}$')
ax.plot(np.sort(meas), c='o', ls='None', marker='o', alpha=.5, ms=s, zorder=0)
ax.plot(np.sort(trimdesign), ls='None', marker='o',c='p', alpha=.5, ms=s, zorder=0)
ax_inset = inset_axes(ax, width="45%", height="45%", loc='lower right')
s = 3
x_min, x_max = 330, 370  # Example x range
y_min, y_max = 5.15, 5.30    # Example y range
ax.add_patch(plt.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min, 
                                   ec='k', fc='None', alpha=1, transform=ax.transData, lw=1))  # Changed alpha to 0.5
                                # Draw diagonal lines connecting the rectangle corners to the inset axes corners
ax_inset.plot(np.sort(meas), c='o', ls='None', marker='o', alpha=.5, ms=s)
ax_inset.plot(np.sort(trimdesign), ls='None', marker='o', c='p', alpha=.5, ms=s)
ax_inset.set_xlim(x_min, x_max)
ax_inset.set_ylim(y_min, y_max)
ax_inset.set_xticks([])
ax_inset.set_yticks([])
ax.set_ylabel('$F$ $(\mathrm{GHz})$')
ax.set_xlabel('Sorted KID index')

s = 5
a=.5
ax = axes['c']
# ax.set_title('$F_D^\mathrm{trim} \\rightarrow L_\mathrm{IDC}^\mathrm{trim}$')
fs = np.linspace(0, 10*97, 11, endpoint=True)
for i, finger in enumerate(fs):
    if not i:
        ax.axhline(finger, c='gray', ls='--', alpha=1, lw=.5, label='Fingers', zorder=0)
    else:
        ax.axhline(finger, c='gray', ls='--', alpha=1, lw=.5, zorder=0)
ax.scatter(design, fingers, c='b', alpha=a, s=s, label='$F_D$', zorder=1)
ax.plot(f0s_at_Q_design, Lfs, c='k', label='Simulation', ls='-', zorder=2)
x = np.linspace(np.nanmin(f0s_at_Q_design), np.nanmax(f0s_at_Q_design), 100)
ax.plot(x, poly(x, *popt), c='k', ls='--', label='Fit $F_M$', zorder=2)
ax.scatter(meas, fingers, alpha=a, s=s, label='$F_M$', c='o', zorder=1)
ax.scatter(trimdesign, newfingers, c='p', ls='-', alpha=a, s=s, label='$F_D^\mathrm{trim}$', zorder=1)
# ax.legend(bbox_to_anchor=(0.5, 1.05), loc='lower center',
#           ncols=2, mode="expand", borderaxespad=0., ax=axes['a'])
ax.set_ylabel('$L_{IDC}$ $(\mathrm{\mu m})$')
ax.set_xlabel('$F$ $(\mathrm{GHz})$')
ax.set_xlim(3.5, 8.5)
ax.grid(False)
ax.set_ylim(0, fs[-1])
ax.set_yticks(fs)
ax.set_yticklabels(['%d' % int(f) for f in fs])

# Create inset axes
x_min, x_max = 5.15, 5.30  # Example x range
y_min, y_max = 380, 440    # Example y range
ax.add_patch(plt.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min, 
                                   ec='k', fc='None', alpha=1, transform=ax.transData, lw=1))  # Changed alpha to 0.5

s = 10
ax_inset = inset_axes(ax, width="55%", height="55%", loc='upper right')
for i, finger in enumerate(fs):
    ax_inset.axhline(finger, c='gray', ls='--', alpha=1, lw=.5, zorder=0)
ax_inset.plot(f0s_at_Q_design, Lfs, c='k', ls='-', zorder=2)
ax_inset.scatter(design, fingers, c='b', alpha=.5, s=s, zorder=1)
ax_inset.plot(x, poly(x, *popt), c='k', ls='--', zorder=2)
ax_inset.scatter(meas, fingers, alpha=.5, s=s, c='o', zorder=1)
ax_inset.scatter(trimdesign, newfingers, c='p', alpha=.5, s=s, zorder=1)
ax_inset.set_xlim(x_min, x_max)
ax_inset.set_ylim(y_min, y_max)
ax_inset.set_xticks([])
ax_inset.set_yticks([])

ax = axes['b']
# ax.set_title('$\mathrm{Difference}$')
_ = ax.hist(trims, bins='auto', fc='p')
ax.set_xlabel('$L_\mathrm{IDC}^\mathrm{trim}-L_\mathrm{IDC}$ $(\mathrm{\mu m})$')
ax.set_ylabel('$\mathrm{Counts}$')
ax = axes['d']
# ax.set_title('$\mathrm{Difference}$')
ax.hist((trimdesign-meas)*1e3, bins='auto', fc='p')
ax.set_xlabel('$F_D^\mathrm{trim}-F_D$ $(\mathrm{MHz})$')
ax.set_ylabel('$\mathrm{Counts}$')
ax = axes['e']
img = plt.imread(r'Mappings\Masters/figures/trim.png')
ax.imshow(img, aspect='equal')
# ax.imshow(img, aspect='auto')
ax.set_xticks([])
ax.set_yticks([])
fig.legend(loc='upper center', bbox_to_anchor=(.83, .2), ncol=3, frameon=True, columnspacing=0.5, handlelength=1)
plt.savefig(r'Mappings\Masters/figures/trimming.pdf')
plt.show()
