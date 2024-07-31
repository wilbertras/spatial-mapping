from functions import *
import time

if __name__ == "__main__":
    fstart = 4.9         # GHz
    fstop = 7.1           # GHz
    scanbw = 100        # MHz
    nr_points = 3201
    power = -112        # dBm
    ifbw = 500         # Hz

    st = time.time()
    freqs, s21 = get_s21(fstart, fstop, scanbw, nr_points, power, ifbw)
    et = time.time()
    elapsed_time = et - st
    print('Elapsed time = %d seconds' % elapsed_time)

    fig = plot_s21(freqs, s21)
    file_path = save_numpy_array(np.stack((freqs, s21), axis=-1).T)
    fig.savefig(file_path.split('.')[0] + '.png')