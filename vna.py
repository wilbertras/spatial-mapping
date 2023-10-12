import pyvisa
import matplotlib.pyplot as plt
import numpy as np
import win32com.client
from scipy.interpolate import interp1d


def get_s21(fstart, fstop, subscanbw, num_points, kidpower, ifbw, timeout=3000000):
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
        subscan = vna_scan(f0, subscanbw, num_points, vna_power, ifbw)
        s21.append(subscan)
    s21 = np.array(s21)
    return freqs, s21


def vna_scan(f0, subscanbw, num_points, vna_power, ifbw, timeout=300000):
    VISA = "GPIB0::16::INSTR"
    rm = pyvisa.ResourceManager() 
    vna = rm.open_resource(VISA)
    vna.timeout = timeout
    if vna.query('*IDN?'):
        print("Connected to VNA:", vna.query('*IDN?;'))
        vna.write(f'SENS1:FREQ:CENT {f0}GHz')
        vna.write(f'SENS1:FREQ:SPAN {subscanbw}GHz')
        vna.write(f'SENS1:BWID {ifbw}Hz;')
        vna.write(f'SENS1:SWE:POIN {num_points};')  # Set the number of sweep points
        vna.write(f'SENS1:PAR:S21:FORM MLOG;')
        vna.write(f'SOUR1:POW1:LEV {vna_power};')
        vna.write('FORM:DATA ASCII')
        # Trigger a single sweep
        vna.write(f'INIT:IMM;')
        # Wait for the measurement to complete (you may need to adjust the wait time)
        try:
            response = vna.query_ascii_values(f'CALC1:DATA? FDATA')
            s21 = np.array(response)
        except:
            s21 = np.zeros(num_points)
            print('No data gathered')
        vna.close()
    else:
        s21 = np.zeros(num_points)
        vna.close()
    return s21


def set_weinschell(attn):
    VISA = "GPIB0::10::INSTR"
    attn_chan2 = 2
    try:
        rm = pyvisa.ResourceManager() 
        try_visa = 1
    except:
        try_visa = 0
        print('NO Pyvisa')
    if try_visa:
        weinshell = rm.open_resource(VISA)
        print(weinshell.query(f'*IDN?'))
        weinshell.write(f'CHAN1;ATTN {attn}')
        weinshell.write(f'CHAN2;ATTN {attn_chan2}')
        # print(weinshell.query(f'CHAN 1;ATTN?'))
        # print(weinshell.query(f'CHAN 2;ATTN?'))
        weinshell.query(f'*OPC?')
            

def plot_s21(freqs, s21):
    fig, ax = plt.subplots()
    _ = ax.plot(freqs, s21)
    ax.set_xlabel('F [GHz]')
    ax.set_ylabel('|$S_{21}$| [dB]')
    ax.set_title('S21 from VNA')
    ax.grid(which='Major')
    plt.show()


def run_labview():
    # Create a LabVIEW Automation object
    lv = win32com.client.Dispatch("LabVIEW.Application")

    # Load a VI (Virtual Instrument) file
    vi_path = r'D:\JBtests_labview\VNASWEEP2.vi'
    vi = lv.GetVIReference(vi_path)
    # vi.SetControlValue("Cooler", 'Bluefors')
    # vi.SetControlValue("AmpNew", 'BF LNF 4-8 + Miteq')
    # vi.SetControlValue("Fstart GHz", 4)
    # vi.SetControlValue("Fstop GHz", 6)
    vi.SetControlValue("BW_per_Scan_[MHz]", 100)
    # vi.SetControlValue("# pts per scan", 6401)
    # vi.SetControlValue("KID Power", -110)
    # vi.SetControlValue("IF Bandwidth (Hz)", 10000)
    # Close LabVIEW
    # Run the VI
    vi.Call()
    # Set input parameter values (assuming your VI has a numeric control named "InputParameter")

    lv.Quit()


# run_labview()
freqs, s21 = get_s21(4, 6, 0.1, 3201, -110, 10000)
print(s21.shape)
