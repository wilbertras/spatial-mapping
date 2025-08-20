import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import matplotlibcolors
plt.style.use('matplotlibrc')

def func(x, a, b):
    return a*x+b


class Mapping:
    """
    Mapping class to handle the mapping of data.
    """

    def __init__(self, file):
        self.file = file
        arr = np.loadtxt(self.file, skiprows=1, delimiter=',', unpack=True)
        self.initialize(arr)

    def initialize(self, arr):
        self.idx = arr[0].astype(int)
        self.row = arr[1].astype(int)
        self.col = arr[2].astype(int)
        self.M = self.row.max() + 1
        self.N = self.col.max() + 1
        self.nr = self.idx.max() + 1
        self.Lf = arr[3]
        self.Lc = arr[4]
        self.fd = arr[5]
        self.fm = arr[6]
        self.map = self.make_map(self.M, self.N)
        self.df_f = self.comp_df_f(self.fd, self.fm)
        self.fd_fit, self.popt = self.fit(self.fd, self.fm)
        self.std = np.nanstd(self.df_f)
        self.df_f_fit = self.comp_df_f(self.fd_fit, self.fm)
        self.std_fit = np.nanstd(self.df_f_fit)

    def remap(self, ids):
        arr = np.loadtxt(self.file, skiprows=1, delimiter=',', unpack=True)
        arr[6:, ids] = np.nan
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
        nanmask = np.isnan(meas)
        popt, pcov = curve_fit(func, design[~nanmask], meas[~nanmask])
        design_fit = func(design, *popt)
        return design_fit, popt

    
    def plot(self):
            fig, axes = plt.subplot_mosaic('abcd', constrained_layout=True, figsize=(13,3))
            ax = axes['a']
            ax.scatter(self.fd, self.fd, label='Design: $f_0^d$')
            ax.scatter(self.fd, self.fm, label='Measured: $f_0^m$')
            x = np.linspace(np.amin(self.fd), np.amax(self.fd), 100)
            ax.plot(x, func(x, *self.popt), label='Corr. design: $f_0^{d^*}$', c='k', ls='--')
            ax.legend()
            min = np.nanmin((self.fd, self.fm))-.1
            max = np.nanmax((self.fd, self.fm))+.1
            xlim = [np.floor(min * 2) / 2, np.ceil(max * 2) / 2]
            ylim = xlim
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_xlabel('Frequency [GHz]')
            ax.set_ylabel('Frequency [GHz]')
            ax.set_yticks(ax.get_xticks())
            ax = axes['b']
            ax.hist(self.df_f_fit, bins='auto', label='$\sigma_f=%.1e$' % self.std_fit, facecolor='o')
            ax.legend()
            max = np.ceil(3/np.nanmax(self.df_f_fit))
            cmin = np.floor(np.nanmin(self.df_f_fit) * max) / max
            cmax = np.ceil(np.nanmax(self.df_f_fit) * max) / max
            ax.set_xlim([cmin, cmax])
            ax.legend()
            ax.set_xlabel('Frac. freq. error')
            ax.set_ylabel('Counts')
            ax = axes['c']
            ax.scatter(self.fd, self.df_f_fit, c=self.df_f_fit, cmap='viridis', vmin=cmin, vmax=cmax, label='$\delta f/f$')
            ax.set_xlim(xlim)
            ax.set_ylim([cmin, cmax])
            ax.set_xlabel('Frequency [GHz]')
            ax.set_ylabel('Frac. freq. error')
            ax = axes['d']
            im = ax.imshow(self.df_f_fit[self.map], origin='lower', vmin=cmin, vmax=cmax)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_xlim(-.5, self.N-.5)
            ax.set_ylim(-.5, self.M-.5)
            ax.grid(False, which='both')
            ax.set_xlabel('$x$ $[px]$')
            ax.set_ylabel('$y$ $[px]$')
            cbar = fig.colorbar(im, ax=ax)
            _ = cbar.ax.set_ylabel('$(f_0^m - f_0^{d^*})/f_0^{d^*}$')
    

class TrimmedMapping(Mapping):
    """
    TrimmedMapping class for handling trimmed data.
    """

    def __init__(self, file):
        super().__init__(file)
        arr = np.loadtxt(file, skiprows=1, delimiter=',', unpack=True)
        self.trimitialize(arr)
    
    def trimitialize(self, arr):
        self.fd = arr[7]
        self.Lf = arr[8]
        self.fm = arr[9]
        self.df_f = self.comp_df_f(self.fd, self.fm)
        self.std = np.nanstd(self.df_f)
        self.fd_fit, self.popt = self.fit(self.fd, self.fm)
        self.df_f_fit = self.comp_df_f(self.fd_fit, self.fm)
        self.std_fit = np.nanstd(self.df_f_fit)

    def remap(self, ids):
        arr = np.loadtxt(self.file, skiprows=1, delimiter=',', unpack=True)
        arr[7:, ids] = np.nan
        self.trimitialize(arr)