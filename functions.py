import pyvisa
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import time

def get_s21(fstart, fstop, subscanbw, num_points, kidpower, ifbw):
    bfout, bfin, rfout = power_calibration()

    totscanbw = fstop - fstart
    num_subscans = int(np.ceil(totscanbw / subscanbw))
    realfstart = fstart
    realfstop = fstart + num_subscans * subscanbw
    f0start = realfstart + subscanbw / 2
    freqs = np.linspace(realfstart, realfstop, num_points*num_subscans)
    scans = []
    # Connect with VI's
    vna = connect2vi("GPIB0::16::INSTR", timeout=3000000)
    weinschell = connect2vi("GPIB0::10::INSTR", timeout=300000)
    # Initialize VNA
    init_vna(vna)
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
        set_weinschell(weinschell, att)
        KID_cryoIn = bfin(f0)
        vna_power = kidpower - KID_cryoIn
        subscan = vna_scan(vna, f0, subscanbw, num_points, vna_power, ifbw, i)
        scans.append(subscan)
    s21 = np.array(scans).flatten()
    print('S21 completed')
    vna.close()
    weinschell.close()
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
    vna.write(f'SYST:PRES')
    vna.write(f'CONT:AUX:OUTP2:VOLT 0')
    vna.write(f'CONT:AUX:OUTP1:VOLT 5')
    
    vna.write('OUTP ON')
    # vna.write('MMEMORY:LOAD "D:\KIDS\KIDs.csa";')
    vna.query(f'*OPC?')
    vna.write('SENS1:SWE:TRIG:POIN OFF;')
    vna.write('TRIG:SCOP CURR;')
    vna.write('INIT1:CONT ON;')

def vna_scan(vna, f0, subscanbw, num_points, vna_power, ifbw, id):

    # Set sweep params
    session = 'Scan%d' % id
    vna.write(f'DISP:WIND:TRAC1:DEL;')
    vna.write(f'CALC1:PAR:DEF "{session}", S21;')
    vna.write(f'CALC1:PAR:SEL "{session}";')
    vna.write(f'DISP:WIND:TRAC1:FEED "{session}";')  
    vna.write(f'DISP:WIND:TRAC1:Y:AUTO')
    vna.write(f'SENS1:FREQ:CENT {f0}GHz;')
    vna.write(f'SENS1:FREQ:SPAN {subscanbw}GHz;')
    vna.write(f'SOUR1:POW1:LEV {vna_power};')
    vna.write(f'SENS1:BWID {ifbw}Hz;')
    vna.write(f'SENS1:SWE:POIN {num_points};') 
    vna.write(f'MMEM:STOR:TRAC:FORM:SNP DB;')
    vna.write('FORM:DATA ASCII;')
    vna.write(f'sens1:swe:trig:poin off;')
    vna.write(f'sens1:swe:time:auto on;')          
    vna.write(f'TRIG:SOUR MAN;')
    vna.write(f'INIT:CONT OFF;')
    power = float(vna.query(f'SOUR1:POW1:LEV?'))
    print('VNA power = %.2f' % power)
    
    # Trigger a single sweep
    vna.write(f'TRIG:SCOP CURR;')
    vna.write(f'INIT:IMM;')
    if vna.query('*OPC?'):
        response = vna.query_ascii_values(f'CALC1:DATA? FDATA;')
        s21 = np.array(response)
        if vna.query(f'*OPC?'):
            print('Subscan %d complete' % (id))
        
    return s21


def set_weinschell(weinschell, attn):
    attn_chan2 = 2
    weinschell.write(f'CHAN1;ATTN {attn}')
    weinschell.write(f'CHAN2;ATTN {attn_chan2}')
            

def plot_s21(freqs, s21):
    fig, ax = plt.subplots()
    _ = ax.plot(freqs, s21, lw=0.2)
    ax.set_xlim([np.amin(freqs), np.amax(freqs)])
    ax.set_xlabel('F [GHz]')
    ax.set_ylabel('|$S_{21}$| [dB]')
    ax.set_title('S21 from VNA')
    ax.grid(which='Major')
    plt.show()



# run
# st = time.time()
# freqs, s21 = get_s21(4, 6, 0.1, 3201, -110, 1000)
# et = time.time()
# elapsed_time = et - st
# print('Elapsed time = %d seconds' % elapsed_time)
# plot_s21(freqs, s21)
