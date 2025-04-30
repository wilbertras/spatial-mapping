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
        arr = np.loadtxt(file, skiprows=1, delimiter=',', unpack=True)
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
        self.fd2 = arr[7]
        self.Lf2 = arr[8]
        self.fm2 = arr[9]
        self.trimmed = np.any(~np.isnan(self.fm2))
        if self.trimmed:
            self.nanmask = np.isnan(self.fm2)
        else:
            self.nanmask = np.isnan(self.fm)
        self.map = self.make_map(self.M, self.N)
        self.df_f = self.comp_df_f(self.fd, self.fm)
        self.fd_fit, self.popt = self.fit(self.fd, self.fm)
        self.std = np.nanstd(self.df_f)
        self.df_f_fit = self.comp_df_f(self.fd_fit, self.fm)
        self.std_fit = np.nanstd(self.df_f_fit)
        if self.trimmed:
            self.df_f2 = self.comp_df_f(self.fd2, self.fm2)
            self.std2 = np.nanstd(self.df_f2)
            self.fd2_fit, self.popt2 = self.fit(self.fd2, self.fm2)
            self.df_f2_fit = self.comp_df_f(self.fd2_fit, self.fm2)
            self.std2_fit = np.nanstd(self.df_f2_fit)

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
        popt, pcov = curve_fit(func, design[~self.nanmask], meas[~self.nanmask])
        design_fit = func(design, *popt)
        return design_fit, popt
    
