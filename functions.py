import pyvisa
import matplotlib.pyplot as plt
import numpy as np
import win32com.client
from scipy.interpolate import interp1d


def get_s21(vna, fstart, fstop, subscanbw, num_points, kidpower, ifbw):
    BFout = np.load('calibration files/BFout.npy')
    BFin = np.load('calibration files/BFin.npy')
    RFout = np.load('calibration files/RFout.npy')
    bfout = interp1d(BFout[:, 0], BFout[:, 1])
    bfin = interp1d(BFin[:, 0], BFin[:, 1])
    rfout = interp1d(RFout[:, 0], RFout[:, 1])

    totscanbw = fstop - fstart
    num_subscans = int(np.ceil(totscanbw / subscanbw))
    realfstart = fstart
    realfstop = fstart + num_subscans * subscanbw
    f0start = realfstart + subscanbw / 2
    freqs = np.linspace(realfstart, realfstop, num_points*num_subscans)
    s21 = []
    for i in range(num_subscans):
        f0 = f0start + i*subscanbw
        if i == 0:
            KID_cryoOUt = bfout(f0)
            GainRFbox = rfout(f0)
            PcryoOUt = kidpower + KID_cryoOUt + GainRFbox
            PVNAin = PcryoOUt - 2
            if PVNAin > 62:
                att = 62
                print('PVNAin > 62')
            elif PVNAin < 2:
                att = 2
                print('PVNAin < 2')
            else:
                att = PVNAin
            weinschell = np.round(att / 2) * 2
        # set_weinschell(weinschell)
        KID_cryoIn = bfin(f0)
        vna_power = kidpower - KID_cryoIn
        subscan = vna_scan(vna, f0, subscanbw, num_points, vna_power, ifbw)
        s21.append(subscan)
    s21 = np.array(s21).flatten()
    return freqs, s21


def connect2vi(VISA, timeout=300000):
    rm = pyvisa.ResourceManager() 
    vi = rm.open_resource(VISA)
    vi.timeout = timeout
    try:
        print('Connected to VI: ' + vi.query('*IDN?'))
    except:
        print('Could not connect to VI')
    return vi


def vna_scan(vna, f0, subscanbw, num_points, vna_power, ifbw):
    # Initialize VNA

    # Set sweep params
    vna.write(f'SENS1:FREQ:CENT {f0}GHz')
    vna.write(f'SENS1:FREQ:SPAN {subscanbw}GHz')
    vna.write(f'SENS1:BWID {ifbw}Hz')
    vna.write(f'SENS1:SWE:POIN {num_points}') 
    vna.write(f'SENS1:PAR:S21:FORM MLOG')
    vna.write(f'SOUR1:POW1:LEV {vna_power}')
    vna.write('FORM:DATA ASCII')
    # Trigger a single sweep
    vna.write(f'INIT:IMM;')
    # Wait for the measurement to complete (you may need to adjust the wait time)
    try:
        response = vna.query_ascii_values(f'CALC1:DATA? FDATA')
        s21 = np.array(response)
        if vna.query(f'*OPC?'):
            vna.close()
    except:
        s21 = np.zeros(num_points)
        print('No data gathered')
        vna.close()
    return s21


def set_weinschell(weinschell, attn):
    attn_chan2 = 2
    weinschell.write(f'CHAN1;ATTN {attn}')
    weinschell.write(f'CHAN2;ATTN {attn_chan2}')
    print(weinschell.query(f'CHAN 1;ATTN?'))
    print(weinschell.query(f'CHAN 2;ATTN?'))
    weinschell.query(f'*OPC?')
            

def plot_s21(freqs, s21):
    fig, ax = plt.subplots()
    _ = ax.plot(freqs, s21)
    ax.set_xlabel('F [GHz]')
    ax.set_ylabel('|$S_{21}$| [dB]')
    ax.set_title('S21 from VNA')
    ax.grid(which='Major')
    plt.show()


# run
freqs, s21 = get_s21(4, 6, 0.1, 3201, -110, 10000)
print(s21.shape)
