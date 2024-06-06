from functions import *
import time

if __name__ == "__main__":
    fstart = 5          # GHz
    fstop = 5.1           # GHz
    scanbw = 100        # MHz
    nr_points = 6401
    power = -110        # dBm
    ifbw = 1000         # Hz

    st = time.time()
    freqs, s21 = get_s21(fstart, fstop, scanbw, nr_points, power, ifbw)
    et = time.time()
    elapsed_time = et - st
    print('Elapsed time = %d seconds' % elapsed_time)

    file_path = save_numpy_array(np.stack((freqs, s21), axis=-1).T)

    fig = plot_s21(freqs, s21)
    fig.savefig(file_path.split('.')[0] + '.png')