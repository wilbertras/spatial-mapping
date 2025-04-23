import numpy as np
from functions import *
from scipy.signal import find_peaks

if __name__ == "__main__":
    sw = 9         # smoothing window length must be even
    mph = 0.01      # minimal peak height 

    arr = open_numpy_array()
    freqs, s21 = arr
    d2s21 = sec_diff(s21, sw)
    locs, props = find_peaks(d2s21, height=mph, prominence=mph)

    fig = plot_pks(freqs, s21, d2s21, locs, mph, sw)
    plot_dfs(locs, freqs)
    plt.show()
    file_path = save_numpy_array(freqs[sw:][locs].T)
    if file_path:
        fig.savefig(file_path.split('.')[0] + '.png')