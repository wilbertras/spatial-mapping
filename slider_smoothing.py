import matplotlib.pyplot as plt
import numpy as np
from functions import *
from scipy.signal import find_peaks, savgol_filter

from matplotlib.widgets import Button, Slider


n0 = 5         # smoothing window length must be even
deg0 = 4      # minimal peak height 

arr = open_numpy_array()
freqs, s21 = arr

smooth_s21 = savgol_filter(s21, n0, deg0)


fig, ax = plt.subplots()
fig.subplots_adjust(bottom=0.25)
ax.plot(freqs, s21, color='gray')
l, = ax.plot(freqs, smooth_s21, lw=1)

ax_n = fig.add_axes([0.25, 0.1, 0.65, 0.03])
ax_deg = fig.add_axes([0.25, 0.15, 0.65, 0.03])

# create the sliders
sn = Slider(
    ax_n, "n", 1, 10,
    valinit=n0, valstep=1,
    color="green"
)

sdeg = Slider(
    ax_deg, "deg", 1, 10,
    valinit=deg0, valstep=1,
    initcolor='none'  # Remove the line marking the valinit position.
)


def update(val):
    n = sn.val
    deg = sdeg.val
    l.set_ydata(savgol_filter(s21, n, deg))
    fig.canvas.draw_idle()


sn.on_changed(update)
sdeg.on_changed(update)

ax_reset = fig.add_axes([0.8, 0.025, 0.1, 0.04])
button = Button(ax_reset, 'Reset', hovercolor='0.975')


def reset(event):
    sn.reset()
    sdeg.reset()
button.on_clicked(reset)


plt.show()