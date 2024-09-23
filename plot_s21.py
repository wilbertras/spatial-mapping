import numpy as np
from functions import *
from scipy.signal import find_peaks

if __name__ == "__main__":
    arr = open_numpy_array(8)
    f = arr[:, 0]
    dB = arr[:, 1]
    phase = arr[:, 2]
    plot_s21(f, dB)
    plt.show()