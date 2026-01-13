from functions import *
import time

if __name__ == "__main__":
    fstart = 4       # GHz
    fstop = 4.5      # GHz
    scanbw = 500       # MHz
    nr_points = 6401
    power = -110        # dBm
    ifbw = 10000         # Hz

    st = time.time()
    freqs, s21 = get_s21(fstart, fstop, scanbw, nr_points, power, ifbw, calfile=False)
    et = time.time()
    elapsed_time = et - st
    print('Elapsed time = %d seconds' % elapsed_time)

    fig = plot_s21(freqs, s21)
    plt.show()
    file_path = save_numpy_array(np.stack((freqs, s21), axis=-1).T)
    if file_path:
        fig.savefig(file_path.split('.')[0] + '.png')


# if __name__ == "__main__":
#     fstart = 8       # GHz
#     fstop = 8.5      # GHz
#     scanbw = 500       # MHz
#     nr_points = 6401
#     power = -110        # dBm
#     ifbw = 10000         # Hz

#     st = time.time()
#     freqs, s21 = get_s21(fstart, fstop, scanbw, nr_points, power, ifbw, calfile=False)
#     et = time.time()
#     elapsed_time = et - st
#     print('Elapsed time = %d seconds' % elapsed_time)