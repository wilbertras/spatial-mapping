from functions import *
import time

if __name__ == "__main__":
    fstart = 5       # GHz
    fstop = 6    # GHz
    scanbw = 100        # MHz
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


# Q=50k, fwhm = 5e9/50e3=.1 MHz, 5 points/fwhm = 5000 points/100MHz 
# 1 GHz, 50 MHz, 3201 points, 10k IFBW: 73s
# 1 GHz, 100 MHz, 6401 points, 10k IFBW: 53s
# 1 GHz, 100 MHz, 3201 points, 10k IFBW: 36s
# 1 GHz, 100 MHz, 3201 points, 1k IFBW: 46s
# 1 GHz, 200 MHz, 6401 points, 1k IFBW: 43s
# 1 GHz, 500 MHz, 16001 points, 1k IFBW: 42s