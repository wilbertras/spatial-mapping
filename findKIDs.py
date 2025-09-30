import numpy as np
from functions import *
from scipy.signal import find_peaks

if __name__ == "__main__":
    n = 5         # smoothing window length must be even
    deg = 3
    mph = 1    # minimal peak height 
    arr = open_numpy_array()

    freqs, s21 = arr
    smooth_s21 = smooth_s21(s21, n, deg)
    d2s21 = sec_diff(smooth_s21)
    locs, props = find_peaks(-smooth_s21, height=0, prominence=.5)

    fig = plot_pks(freqs, s21, d2s21, locs, mph)
    plot_dfs(locs, freqs)
    plt.show()
    file_path = save_numpy_array(freqs[locs].T)
    if file_path:
        fig.savefig(file_path.split('.')[0] + '.png')