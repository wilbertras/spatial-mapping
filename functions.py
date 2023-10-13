import pyvisa
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d


def get_s21(fstart, fstop, subscanbw, num_points, kidpower, ifbw):
    bfout, bfin, rfout = power_calibration()

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
            elif PVNAin < 2:
                att = 2
            else:
                att = PVNAin
            att = np.round(att / 2) * 2
        set_weinschell(att)
        KID_cryoIn = bfin(f0)
        vna_power = kidpower - KID_cryoIn
        subscan = vna_scan(f0, subscanbw, num_points, vna_power, ifbw, i)
        s21.append(subscan)
    s21 = np.array(s21).flatten()
    print('S21 completed')
    return freqs, s21


def power_calibration():
    BFout = np.load('calibration files/BFout.npy')
    BFin = np.load('calibration files/BFin.npy')
    RFout = np.load('calibration files/RFout.npy')
    bfout = interp1d(BFout[:, 0], BFout[:, 1])
    bfin = interp1d(BFin[:, 0], BFin[:, 1])
    rfout = interp1d(RFout[:, 0], RFout[:, 1])
    return bfout, bfin, rfout

def connect2vi(VISA, timeout=300000):
    rm = pyvisa.ResourceManager() 
    vi = rm.open_resource(VISA)
    vi.timeout = timeout
    try:
        vi.query('*IDN?')
    except:
        print('Could not connect to VI')
    return vi

def init_vna(vna):
    vna.write('OUTP ON')
    vna.write('MMEMORY:LOAD "D:\KIDS\KIDs.csa";')
    vna.query(f'*OPC?')
    vna.write('SENS1:SWE:TRIG:POIN OFF;')
    vna.write('TRIG:SCOP curr;')
    vna.write('INIT1:CONT ON;')

def vna_scan(f0, subscanbw, num_points, vna_power, ifbw, id):
    # Reconnect with VNA
    vna = connect2vi("GPIB0::16::INSTR", timeout=3000000)
    # Initialize VNA
    init_vna(vna)
    # Set sweep params
    session = 'Scan%d' % id
    vna.write(f'DISP:WIND1:TRAC1:DEL;')
    vna.write(f'DISP:WIND1:TRAC1:feed "{session}";')
    vna.write(f'DISP:WIND1:TRAC1:Y:')
    vna.write(f'CALC1:PAR:DEF "{session}", S21;')
    vna.write(f'CALC1:PAR:SEL "{session}";')
    vna.write(f'SENS1:FREQ:CENT {f0}GHz;')
    vna.write(f'SENS1:FREQ:SPAN {subscanbw}GHz;')
    vna.write(f'SENS1:BWID {ifbw}Hz;')
    vna.write(f'SENS1:SWE:POIN {num_points};') 
    vna.write(f'SENS1:PAR:S21:FORM MLOG;')
    vna.write(f'SOUR1:POW1:LEV {vna_power};')
    power = float(vna.query(f'SOUR1:POW1:LEV?'))
    print('VNA power = %.1f' % power)
    vna.write('FORM:DATA ASCII;')
    # Trigger a single sweep
    vna.write(f'INIT:IMM;')
    # Wait for the measurement to complete (you may need to adjust the wait time)
    try:
        response = vna.query_ascii_values(f'CALC1:DATA? FDATA;')
        s21 = np.array(response)
        if vna.query(f'*OPC?'):
            vna.close()
            print('Subscan %d complete' % (id))
    except:
        s21 = np.zeros(num_points)
        print('No data gathered')
        vna.close()
    return s21


def set_weinschell(attn):
    weinschell = connect2vi("GPIB0::10::INSTR", timeout=300000)
    attn_chan2 = 2
    weinschell.write(f'CHAN1;ATTN {attn}')
    weinschell.write(f'CHAN2;ATTN {attn_chan2}')
    if weinschell.query(f'*OPC?'):
        weinschell.close()
            

def plot_s21(freqs, s21):
    fig, ax = plt.subplots()
    _ = ax.plot(freqs, s21, lw=0.2)
    ax.set_xlabel('F [GHz]')
    ax.set_ylabel('|$S_{21}$| [dB]')
    ax.set_title('S21 from VNA')
    ax.grid(which='Major')
    plt.show()



# run
freqs, s21 = get_s21(4, 5, 0.5, 3201, -110, 10000)
plot_s21(freqs, s21)
