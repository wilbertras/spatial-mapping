import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import pickle
import functions as ft
import matplotlibcolors
plt.style.use('matplotlibrc')

class Mapping:
    """
    Mapping class to handle the mapping of data.
    """

    def __init__(self, file, mask_edges=False, type='', Q=20e3, min_lw_spacing=4, deg=2):
        self.file = file
        with open(self.file, 'rb') as f:
            arr = pickle.load(f)
        self.Q = Q
        self.min_lw_spacing = min_lw_spacing
        self.deg = deg
        self.type = type
        self.mask_edges = mask_edges
        self.initialize(arr)

    def initialize(self, arr):
        design = self.type + 'design'
        measured = self.type + 'measured'

        self.row = arr['design']['row'].astype(int)
        self.col = arr['design']['col'].astype(int)
        self.M = self.row.max() + 1
        self.N = self.col.max() + 1
        self.nr = self.M * self.N
        self.idx = np.arange(self.nr)
        if 'fingers' in arr[design]:
            self.Lf = arr[design]['fingers']
        if 'couplers' in arr[design]:
            self.Lc = arr[design]['couplers']
        self.fd = arr[design]['f0']
        if np.any(self.fd > 1e9):
            self.fd = self.fd / 1e9  # Convert to GHz
        self.fm = arr[measured]['f0']
        if np.any(self.fm > 1e9):
            self.fm = self.fm / 1e9  # Convert to GHz
        self.map = self.make_map(self.M, self.N)
        self.df_f = self.comp_df_f(self.fd, self.fm)
        edge_mask = np.zeros(self.nr, dtype=bool)
        ids_center = self.map[1:-1, :-1].flatten()
        edge_mask[ids_center] = True
        self.edge_mask = edge_mask
        self.fit_func = self.fit(self.fd, self.fm)
        self.fd_fit = self.fit_func(self.fd)
        self.std = np.nanstd(self.df_f)
        self.df_f_fit = self.comp_df_f(self.fd_fit, self.fm)
        self.std_fit = np.nanstd(self.df_f_fit)
        self.fab_yield = np.sum(~np.isnan(self.fm)) / self.nr
        self.rel_spacings, self.suff_spaced, self.usable_yield = self.spacings(self.Q, self.min_lw_spacing)


    def remap(self, ids):
        with open(self.file, 'rb') as f:
            arr = pickle.load(f)
        arr[self.type + 'measured']['f0'][ids] = np.nan
        self.initialize(arr)

    def nan_edges(self):
        with open(self.file, 'rb') as f:
            arr = pickle.load(f)
        arr[self.type + 'measured'][~self.edge_mask] = np.nan
        self.initialize(arr)

    def make_map(self, M, N):
        """
        Create kid_id_board
        """
        map = np.ones((M, N), dtype=int)
        map[self.row, self.col] = self.idx
        return map
    
    def comp_df_f(self, design, meas):
        """
        Compute fractional frequency error
        """
        df_f = (meas - design) / design
        return df_f
    
    def fit(self, design, meas):
        nanmask = ~(np.isnan(meas) | np.isnan(design))
        if self.mask_edges:
           mask = self.edge_mask & nanmask
        else:
            mask = nanmask
        coeff = np.polyfit(design[mask], meas[mask], self.deg)
        return np.poly1d(coeff)

    def spacings(self, Q, nr_lw_spacing):
        argnans = np.isnan(self.fm)
        argsort = np.argsort(self.fm[~argnans])
        sorted_f0s = self.fm[~argnans][argsort]
        diffs = np.diff(sorted_f0s)
        fwhms = sorted_f0s / Q
        rel_diffs = diffs/fwhms[:-1]
        too_close = rel_diffs<nr_lw_spacing
        suff_spaced = np.ones(sorted_f0s.shape)
        suff_spaced[1:] -= too_close
        suff_spaced[:-1] -= too_close
        suff_spaced = (suff_spaced==True)
        return rel_diffs, suff_spaced, np.sum(suff_spaced) / self.nr

    def plot_design_vs_meas(self, ax=None, flim=[None, None]):
        if not ax:
            fig, ax = plt.subplots(figsize=(5,5))
        ax.scatter(self.fd, self.fd, label='Design')
        ax.scatter(self.fd, self.fm, label='Measured')
        x = np.linspace(np.nanmin(self.fd), np.nanmax(self.fd), 100)
        ax.plot(x, self.fit_func(x), label='Fit, deg=%d' % self.deg, c='k', ls='--')
        ax.legend()
        ax.set_xlim(flim)
        ax.set_ylim(flim)
        ax.set_xlabel('Frequency [GHz]')
        ax.set_ylabel('Frequency [GHz]')
        ax.set_yticks(ax.get_xticks())

    def plot_scatter_histogram(self, ax=None, clim=[None, None]):
        if not ax:
            fig, ax = plt.subplots(figsize=(5,5))
        ax.hist(self.df_f_fit, bins='auto', label='$\sigma_f=%.1e$' % self.std_fit, facecolor='o')
        ax.legend()
        ax.set_xlim(clim)
        ax.set_xlabel('Frac. freq. error')
        ax.set_ylabel('Counts')

    def plot_scatter_vs_freq(self, ax=None, clim=[None, None], flim=[None, None]):
        if not ax:
            fig, ax = plt.subplots(figsize=(5,5))
        ax.scatter(self.fd, self.df_f_fit, c=self.df_f_fit, cmap='viridis', vmin=clim[0], vmax=clim[1])
        ax.set_xlabel('Frequency [GHz]')
        ax.set_ylabel('Frac. freq. error')
        ax.set_xlim(flim)
        ax.set_ylim(clim)
    
    def plot_scatter_map(self, ax=None, clim=[None, None]):
        if not ax:
            fig, ax = plt.subplots(figsize=(5,5))
        im = ax.imshow(self.df_f_fit[self.map], origin='lower', vmin=clim[0], vmax=clim[1])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xlim(-.5, self.N-.5)
        ax.set_ylim(-.5, self.M-.5)
        ax.grid(False, which='both')
        ax.set_xlabel('$x$ $[px]$')
        ax.set_ylabel('$y$ $[px]$')
        cbar = plt.colorbar(im, ax=ax)
        _ = cbar.ax.set_ylabel('Frac. freq. error')
    
    def plot_spacings(self, ax=None, lwlim=[None, None]):
        if not ax:
            fig, ax = plt.subplots()
        lo_id = np.argmax(self.rel_spacings)
        _ = ax.hist(self.rel_spacings[:lo_id], bins='auto', alpha=.5, label='lower band')
        _ = ax.hist(self.rel_spacings[lo_id+1:], bins='auto', alpha=.5, label='upper band')
        ax.axvline(self.min_lw_spacing, c='r', ls='--', lw=1, label='min. spacing (%d lw, Q=%dk)' % (self.min_lw_spacing, self.Q/1e3))
        ax.set_xlim(lwlim)
        ax.set_xlabel('Freq. spacing (# lw)')
        ax.set_ylabel('# KIDs')
        ax.set_title('Usable yield: %.2f%%' % (self.usable_yield*1e2))
        ax.legend()

    def plot_overview(self, clim=[None, None], flim=[None, None], lwlim=[None, None]):
            fig, axes = plt.subplot_mosaic('abcd', constrained_layout=True, figsize=(13,3))
            self.plot_design_vs_meas(ax=axes['a'], flim=flim)
            # self.plot_scatter_histogram(ax=axes['b'], clim=clim)
            self.plot_scatter_vs_freq(ax=axes['b'], clim=clim, flim=flim)
            self.plot_scatter_map(ax=axes['c'], clim=clim)
            self.plot_spacings(ax=axes['d'], lwlim=lwlim)


    def remove_spatial(self):
        df_f_map = self.df_f_fit[self.map]
        if self.mask_edges:
            x = np.arange(self.N-1)
            y = np.arange(1, self.M-1)
            Z = df_f_map[1:-1, :-1]
        else:
            x = np.arange(self.N)
            y = np.arange(self.M)
            Z = df_f_map

        X, Y = np.meshgrid(x, y)
        degree = 2
        coeffs, powers = ft.fit_poly2d_nan_safe(X, Y, Z, degree)
        X_fit, Y_fit = np.meshgrid(np.arange(self.M), np.arange(self.N))
        Z_fit = ft.eval_poly2d_on_grid(X_fit, Y_fit, coeffs, powers)
        self.df_f_spatial = np.zeros_like(self.df_f_fit)
        self.df_f_spatial[self.map.flatten()] = (df_f_map - Z_fit).flatten()
        