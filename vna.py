import pyvisa
import matplotlib.pyplot as plt
import numpy as np
import win32com.client


def get_s21(f0, bw, num_points, power, ifbw, timeout=3000000):
    freqs = np.linspace(f0-bw/2, f0+bw/2, num_points)
    VISA = "GPIB0::16::INSTR"
    try:
        rm = pyvisa.ResourceManager() 
        vna = rm.open_resource(VISA)
        vna.timeout = timeout
        if vna.query('*IDN?;'):
            print("Connected to VNA:", vna.query('*IDN?;'))
            # Configure VNA parameters (example)
            kidname = 'kidx, timex'
            par = 'S21'
            # _, status = vna.write(f'OUTP ON;')
            # _, status = vna.write(f'MMEMORY:LOAD "D:\KIDS\KIDs.csa";')  # Define S21 measurement
            # _, status = vna.write(f'*OPC;')
            # _, status = vna.write(f'SENS1:SWE:TRIG:POIN OFF; TRIG:SCOP CURR; INIT1:CONT ON;')
            # _, status = vna.write(f'disp:wind1:trac1:DEL;')
            # _, status = vna.write(f'CALC1:PAR:DEF "{kidname}", "{par}"')
            # _, status = vna.write(f'CALC1:PAR:SEL "{kidname}";')
            # _, status = vna.write(f'disp:wind1:trac1:feed "{kidname}";')
            # _, status = vna.write(f'DISP:WIND1:TRAC1:Y:AUTO')
            _, status = vna.write(f'SENS1:FREQ:CENT {f0}GHz; SENS1:FREQ:SPAN {bw}GHz; SENS1:BWID {ifbw}Hz; SENS1:SWE:POIN {num_points}; SENS1:PAR:S21:FORM MLOG;SOUR1:POW1:LEV {power};')
            # vna.write(f'SENS1:FREQ:SPAN {bw}GHz;')
            # vna.write(f'SENS1:BWID {ifbw}Hz;')
            # vna.write(f'SENS1:SWE:POIN {num_points};')  # Set the number of sweep points
            # vna.write(f'SENS1:PAR:S21:FORM MLOG;')
            # vna.write(f'SOUR1:POW1:LEV {power};')
            vna.write('FORM:DATA ASCII')
            # Trigger a single sweep
            vna.write(f'INIT:IMM;')

            # Wait for the measurement to complete (you may need to adjust the wait time)
            response = vna.write(f'CALC1:DATA? FDATA;')
            s21 = np.array(response)
            vna.close()
        else:
            s21 = np.zeros(num_points)
    except Exception as err: 
        print("Exception")  
        s21 = np.zeros(freqs.shape)
    finally:
        print("Sweep complete")
    return freqs, s21


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
# # freqs, s21 = get_s21(6, 0.1, 101, -59, 10000)
