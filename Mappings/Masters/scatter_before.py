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

Q = 50e3
chi = 2
dir = r'Mappings/Masters/chips/'
offset = True
deg = 2
LT361chip7 = Mapping(dir + 'LT361chip7_master.pkl', Q=Q, min_lw_spacing=chi, deg=deg, offset=offset)
LT402chip6 = Mapping(dir + 'LT402chip6_master.pkl', Q=Q, min_lw_spacing=chi, deg=deg, offset=offset, mask_edges=True)
LT402chip6_trim = Mapping(dir + 'LT402chip6_master.pkl', type='trim', Q=Q, min_lw_spacing=chi, deg=deg, offset=offset, mask_edges=True)
chips = [LT402chip6, LT361chip7]
# chips = [LT402chip1, LT361chip7, LT402chip6]
labels = ['A', 'B']
cs = 'bb'
cs2 = 'oo'
cmap = cm.get_cmap('Oranges', 256)
cmap = LinearSegmentedColormap.from_list("custom_orange", cmap(np.linspace(0.05, 1.0, 256)))
  # Constrain the range starting at 10%
# Create a custom colormap from black to orange
# custom_cmap = LinearSegmentedColormap.from_list("custom_orange", ['w', 'o'])
# cmap = custom_cmap
flim = [3.5, 8.5]
dflim = [-20, 80]
dfticks = [-20, 0, 20, 40, 60, 80]
# clims = [[-20, 5], [-2, 2], [-25, 25]]
clim = [-20, 30]
ylims = [[0, 300], [0, 80], [0, 300]]
bins=np.linspace(-20, 80, 80)
cticks = [-20, -10, 0, 10, 20, 30]
# pos = 'abcdefghi'
# fig, axes = plt.subplot_mosaic('abc;def', figsize=(18.5/2.54, 13/2.54), constrained_layout=True)
# for i, chip in enumerate(chips):
#     if offset:
#         fit = np.polyval(chip.coeff, chip.fd)
#     else:
#         if deg==1:
#             fit = linear(chip.fd, *chip.coeff)
#         elif deg==2:
#             fit = quadratic(chip.fd, *chip.coeff)
#     ax = axes[pos[i*3]]
#     df_f_fit = (fit - chip.fd)/chip.fd
#     ax.scatter(chip.fd, chip.df_f*1e3, color=cs[i], label='$\epsilon$', alpha=.75, s=2)
#     ax.plot(chip.fd, df_f_fit*1e3, linestyle='--', color='k', label='$\mathrm{fit}$')
#     # ax.scatter(chip.fd, chip.df_f_fit*1e3, c=chip.df_f_fit*1e3, vmin=clim[0], vmax=clim[1], label='$\epsilon-\epsilon_f$', alpha=.75, s=2)
#     ax.scatter(chip.fd, chip.df_f_fit*1e3, fc=cs2[i], label='$\epsilon-\epsilon_f$', alpha=.75, s=2)
#     # ax.plot(chip.fd, np.zeros(chip.fd.shape), linestyle='-', color='k', label='data')
#     ax.set_xlim(flim)
#     ax.set_ylim(dflim)
#     ax.set_yticks(dfticks)
#     ax.set_ylabel('$\epsilon\\times 10^{3}$')

#     if not i:
#         ax.annotate('$\epsilon_M$', xy=(4.2, 60), ha='center', va='center', color='b', fontsize='small', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
#         ax.annotate('$\epsilon_f$', xy=(8, 45), ha='center', va='center', color='k', fontsize='small', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
#         ax.annotate('$\epsilon_M-\epsilon_f$', xy=(7.75, 20), ha='center', va='center', color=cs2[i], fontsize='small', bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
#         ax.annotate('$\epsilon_M^\mathrm{trim}$', xy=(0.9, 0.1), xycoords='axes fraction', ha='center', va='center', color='p', fontsize='small')

#     # ax.legend(loc='upper right', ncols=3, handlelength=1, frameon=False, fontsize='small')
#     ax = axes[pos[1+i*3]]
#     if i==0:
#         empty_array = np.full((32, 32), np.nan)  # Create a 32x32 array filled with NaNs
#         center_x, center_y = (empty_array.shape[0] // 2, empty_array.shape[1] // 2)
#         half_size = chip.df_f_fit[chip.map].shape[0] // 2
#         # Place df_f_spatial[chips[1].map] in the center
#         empty_array[center_x - half_size:center_x + half_size, center_y - half_size:center_y + half_size] = chip.df_f_fit[chip.map]*1e3
#     else:
#         empty_array = chip.df_f_fit[chip.map]*1e3
#     y_indices, x_indices = np.where(np.isnan(empty_array))
#     if not i:
#         ax.scatter(x_indices, y_indices, color='k', marker='x', s=20, label='NaN', alpha=1, linewidth=.5)
#     ax.imshow(empty_array, origin='lower', cmap=cmap, vmin=clim[0], vmax=clim[1])
#     # Set NaN values to a specific color (black)
#     ax.set_xticks([])
#     ax.set_yticks([])
#     ax.set_xlabel('$x$')
#     ax.set_ylabel('$y$')

#     # cbar = fig.colorbar(ax.imshow(empty_array, origin='lower', cmap='viridis', clim=clims[i]), ax=ax, orientation='horizontal', location='top', pad=.01)
#     ax = axes[pos[2+i*3]]
#     ax.hist(chip.df_f*1e3, facecolor=cs[i], alpha=.75, label='data', orientation='vertical', bins=bins)
#     chip.remove_spatial()
#     ax.hist(chip.df_f_spatial*1e3, facecolor='gray', alpha=.75, orientation='vertical', bins=bins)
#     # ax.hist(chip.df_f_fit*1e3, facecolor='gray', alpha=.5, label='data', orientation='horizontal')
#     # ax.hist(chip.df_f_fit, facecolor='None', edgecolor=cs[i], alpha=.5, label='data', orientation='horizontal')
#     # ax.annotate('$\mu=%.0f\\times 10^{-3}$' % (np.nanmean(chip.df_f)*1e3), xy=(0.5, 0.9), xycoords='axes fraction', ha='center', va='center', color=cs[i], fontsize='small')
#     if not i:
#         ax.annotate('$\epsilon_M^\mathrm{trim}$', xy=(0.6, 0.9), xycoords='axes fraction', ha='center', va='center', color='p', fontsize='small')
#         ax.annotate('$\epsilon_M^*=\epsilon_M-\epsilon_f-\epsilon_{xy}$', xy=(0.6, 0.6), xycoords='axes fraction', ha='center', va='center', color='gray', fontsize='small')
#         ax.annotate('$\epsilon_M$', xy=(0.65, 0.3), xycoords='axes fraction', ha='center', va='center', color='b', fontsize='small')
#         ax.annotate('$\sigma_M^\mathrm{trim}=%.2G\\times 10^{-3}$' % (np.nanstd(LT402chip6_trim.df_f)*1e3), xy=(0.6, 0.8), xycoords='axes fraction', ha='center', va='center', color='p', fontsize='small')
#     ax.annotate('$\sigma_M^*=%.2G\\times 10^{-3}$' % (np.nanstd(chip.df_f_spatial)*1e3), xy=(0.6, 0.5+i*.1), xycoords='axes fraction', ha='center', va='center', color='gray', fontsize='small')
#     ax.annotate('$\sigma_M=%.2G\\times 10^{-3}$' % (np.nanstd(chip.df_f)*1e3), xy=(0.65, 0.2+i*.1), xycoords='axes fraction', ha='center', va='center', color=cs[i], fontsize='small')
#     ax.set_xlim(dflim)
#     ax.set_xticks(dfticks)
#     secax = ax.twinx()
#     ax.set_yticklabels([])
#     secax.set_ylabel('$\mathrm{Counts}$')
#     secax.set_ylim(ylims[i])
#     ax.set_ylim(ylims[i])
# axes['d'].set_xlabel('$F_D$ $(\mathrm{GHz})$')
# axes['f'].set_xlabel('$\epsilon\\times 10^{3}$')
# axes['a'].set_title('Frequency dependence')
# axes['b'].set_title('Spatial dependence')
# axes['c'].set_title('Frequency scatter')
# cbar = fig.colorbar(axes['e'].imshow(empty_array, origin='lower', cmap=cmap, clim=clim), ax=axes['e'], orientation='horizontal', location='bottom', pad=-.1)
# cbar.set_label('$(\epsilon_M-\epsilon_f)\\times 10^{3}$')
# cbar.set_ticks(cticks)

# axes['c'].hist(LT402chip6_trim.df_f*1e3, facecolor='p', bins=bins, alpha=1, zorder=0)
# axes['a'].scatter(LT402chip6_trim.fd, LT402chip6_trim.df_f*1e3, facecolor='p', alpha=.75, s=2)
# # plt.savefig('Mappings/Masters/figures/scatter_before.png')
# plt.show()

