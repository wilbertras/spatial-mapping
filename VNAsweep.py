from functions import *
import time


fstart = 5          # GHz
fstop = 8           # GHz
scanbw = 100        # MHz
nr_points = 6401
power = -110        # dBm
ifbw = 1000         # Hz

st = time.time()
freqs, s21 = get_s21(fstart, fstop, scanbw, nr_points, power, ifbw)
et = time.time()
elapsed_time = et - st
print('Elapsed time = %d seconds' % elapsed_time)
plot_s21(freqs, s21)
np.save('VNAsweep.npy', np.stack((freqs, s21), axis=-1).T)
print('saved: VNAsweep.npy')
plt.savefig('VNAsweep.png')
print('saved: VNAsweep.png')