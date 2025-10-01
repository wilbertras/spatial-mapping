import matplotlib.pyplot as plt
import numpy as np
from functions import *
from scipy.signal import find_peaks, savgol_filter

from matplotlib.widgets import Button, Slider


n0 = 5         # smoothing window length must be even
deg0 = 4      # minimal peak height 
mph0 = 0
mpp0 = 0

arr = open_numpy_array()
freqs, s21 = arr

smooth_s21 = savgol_filter(s21, n0, deg0)
locs = find_peaks(-smooth_s21, height=mph0, prominence=mpp0)[0]

fig, ax = plt.subplots(figsize=(8,4))
fig.subplots_adjust(bottom=0.4)
ax.plot(freqs, s21, color='gray')
l, = ax.plot(freqs, smooth_s21, lw=1)
s, = ax.plot(freqs[locs], smooth_s21[locs], color='red', ls='None', marker='o', label='%d kids' % len(locs))    
ax.legend()

ax_n = fig.add_axes([0.25, 0.1, 0.65, 0.03])
ax_deg = fig.add_axes([0.25, 0.15, 0.65, 0.03])
ax_mph = fig.add_axes([0.25, 0.2, 0.65, 0.03])
ax_mpp = fig.add_axes([0.25, 0.25, 0.65, 0.03])

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

smph = Slider(
    ax_mph, "mph", -50, 0,
    valinit=mph0, valstep=1,
    initcolor='none'  # Remove the line marking the valinit position.
)

smpp = Slider(
    ax_mpp, "mpp", 0, 30,
    valinit=mpp0, valstep=1,
    initcolor='none'  # Remove the line marking the valinit position.
)

def update(val):
    n = sn.val
    deg = sdeg.val
    mph = smph.val
    mpp = smpp.val
    smooth_s21 = savgol_filter(s21, n, deg)
    l.set_ydata(smooth_s21)
    locs = find_peaks(-smooth_s21, height=mph, prominence=mpp)[0]
    s.set_xdata(freqs[locs])
    s.set_ydata(smooth_s21[locs])
    ax.legend(handles=[s], labels=['%d kids' % len(locs)])
    fig.canvas.draw_idle()


sn.on_changed(update)
sdeg.on_changed(update)
smph.on_changed(update)
smpp.on_changed(update)

ax_reset = fig.add_axes([0.8, 0.025, 0.1, 0.04])
button = Button(ax_reset, 'Reset', hovercolor='0.975')


def reset(event):
    sn.reset()
    sdeg.reset()
    smph.reset()
    smpp.reset()
button.on_clicked(reset)


plt.show()