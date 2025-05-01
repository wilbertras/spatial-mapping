import numpy as np
import matplotlib.pyplot as plt
from functions import *
from scipy.signal import find_peaks, convolve
from scipy.optimize import curve_fit
from glob import glob
from natsort import natsorted
import matplotlibcolors
import time
plt.style.use('matplotlibrc')


def line(x, a, b):
    return a * x + b 

def inv_line(x, a, b):
    return (x - b) / a

if __name__ == "__main__":
    fstart = 4         # GHz
    fstop = 5         # GHz
    scanbw = 200       # MHz
    nr_points = 6401
    power = -110        # dBm
    ifbw = 10000         # Hz
    sw = 3 
    mph = 0.2 
    window = np.ones(sw)/sw
    vna = False

    if vna: 
        try: 
            vna = connect2vi("GPIB0::16::INSTR", timeout=3000000)
            weinschell = connect2vi("GPIB0::10::INSTR", timeout=300000)
            vna = True
        except:
            vna = False

    arr = np.loadtxt('LTD25/LT361w2chip9_master.txt', skiprows=1, delimiter=',')    
    idx = arr[:, 0].astype(int)
    row = arr[:, 1].astype(int)
    col = arr[:, 2].astype(int)
    M = row.max() + 1
    N = col.max() + 1
    nr = idx.max() + 1
    f0_mapped = arr[:, 6]
    nanmask = np.isnan(f0_mapped)
    map = np.ones((M, N), dtype=int)
    map[row, col] = idx

    if vna:
        f_dark, s21_dark = get_s21(fstart, fstop, scanbw, nr_points, power, ifbw, calfile='D:\KIDS\KIDS.csa')
    else:    
        f_dark, s21_dark = np.load('Mappings\LT361w2chip9\S21s_20241211_17h7\S21_dark.npy')
        files = natsorted(glob('Mappings\LT361w2chip9\S21s_20241211_17h7\*_*y0.npy'))
    
    frac_bw = .1
    end = int(frac_bw * len(f_dark))
    f_dark = f_dark[:end]
    s21_dark = s21_dark[:end]
    smooth_s21_dark = convolve(s21_dark, window, mode='same')
    d2s21 = sec_diff(s21_dark, sw)
    locs, props = find_peaks(d2s21, height=mph, prominence=mph)
    dw = len(s21_dark)-len(d2s21) - sw
    locs += dw
    f0_dark = f_dark[locs]
    s21_dark_min = s21_dark[locs]

    f0_mapped_corr = f0_mapped
    tones = np.empty(nr)
    tones[:] = np.nan
    tone_locs = np.empty(nr)
    tone_locs[:] = np.nan
    for id, f0 in enumerate(f0_mapped_corr):
        if ~np.isnan(f0):
            closest =  np.argmin(np.abs(f0 - f0_dark))
            tone = f0_dark[closest]
            tone_loc = locs[closest]
            if id > 0:
                prev_tone = tones[~np.isnan(tones)][-1]
                if tone == prev_tone:
                    curr_diff = np.abs(tone - f0)
                    prev_diff = np.abs(prev_tone - f0)
                    if curr_diff < prev_diff:
                        tones[tones==prev_tone] = np.nan
                        tone_locs[tones==prev_tone] = np.nan
                    else:
                        tone = np.nan
                        tone_loc = np.nan
                else:
                    pass
            else:
                pass
        else:
            tone = np.nan
            tone_loc = np.nan
        tones[id] = tone
        tone_locs[id] = tone_loc
    nantones = np.isnan(tone_locs)
    s21_dark_min = np.empty(nr)
    s21_dark_min[:] = np.nan
    s21_dark_min[~nantones] = smooth_s21_dark[tone_locs[~nantones].astype(int)]


    fig, ax = plt.subplots()
    i = 0
    while plt.fignum_exists(fig.number):
        if vna:
            f, s21 = get_s21(fstart, fstop, scanbw, nr_points, power, ifbw, calfile='D:\KIDS\KIDS.csa')
        else:
            f, s21 = np.load(files[i])
            i += 1
        smooth_s21 = convolve(s21, window, mode='same')
        s21_min = smooth_s21[tone_locs[~nantones].astype(int)]
        
        ds21_min = np.empty(nr)
        ds21_min[:] = np.nan
        ds21_min[~nantones] = s21_min - s21_dark_min[~nantones] 
        ax.clear()
        ax.imshow(ds21_min[map], origin='lower', cmap='viridis', vmin=np.nanmin(ds21_min), vmax=np.nanmax(ds21_min))
        plt.draw()
        if i == len(files):
            i = 0
        else:
            plt.pause(0.05)
    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(f, smooth_s21, lw=.5, label='s21')
    ax.plot(f_dark, smooth_s21_dark, lw=.5, label='dark')
    ax.scatter(f[tone_locs[~nantones].astype(int)], s21_min, lw=.5, label='image')
    ax.scatter(f0_mapped, np.zeros(nr), marker='v', color='None', edgecolor='y', label='mapped')
    ax.scatter(f_dark[locs], smooth_s21_dark[locs], marker='^', color='None', edgecolor='r', label='dark')
    ax.scatter(tones, s21_dark_min, marker='o', color='g', edgecolor='None', label='tones')
    ax.set_ylim([np.amin(smooth_s21_dark), np.amax(smooth_s21_dark)])
    ax.legend()
    plt.show()