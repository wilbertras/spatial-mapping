import pyvisa
import matplotlib.pyplot as plt
import numpy as np


def get_s21(f0, bw, num_points, power, ifbw, timeout=3000000):
    freqs = np.linspace(f0-bw/2, f0+bw/2, num_points)
    VISA = "GPIB0::16::INSTR"
    try:
        rm = pyvisa.ResourceManager() 
        try_visa = 1
    except:
        try_visa = 0
        print('No pyvisa')
    if try_visa:
        vna = rm.open_resource(VISA)
        vna.timeout = timeout
        if vna.query('*IDN?'):
            print("Connected to VNA:", vna.query('*IDN?'))
            # Configure VNA parameters (example)
            vna.write(f'CALC1:PAR:S21:DEF')  # Define S21 measurement
            vna.write(f'SENS1:FREQ:CENT {f0}GHz')
            vna.write(f'SENS1:FREQ:SPAN {bw}GHz')
            vna.write(f'SENS1:BWID {ifbw}Hz')
            vna.write(f'SENS1:SWE:POIN {num_points}')  # Set the number of sweep points
            vna.write(f'SENS1:PAR:S21:FORM MLOG')
            vna.write(f'SOUR1:POW1:LEV {power}')

            # Trigger a single sweep
            vna.write(f'INIT:IMM')

            # Wait for the measurement to complete (you may need to adjust the wait time)
            
           
            try:
                # vna.query(f'*OPC?')
                # Perform your query with the adjusted timeout
                response = vna.query_ascii_values(f'CALC1:DATA? FDATA')
                s21 = np.array(response)
                # If the query succeeds, you can process the response here
                print("Query result positive")
            except pyvisa.VisaIOError as e:
                s21 = np.zeros(num_points)
                # input(timeout)
                # vna.timeout = timeout
                if "Timeout" in str(e):
                    print("Query timed out. Consider increasing the timeout.")
                else:
                    print("An error occurred:", e)

            vna.close()
        else:
            s21 = np.zeros(num_points)
    else:
        s21 = np.zeros(num_points)
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
