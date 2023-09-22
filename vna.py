# import pyvisa
import matplotlib.pyplot as plt
import numpy as np


def get_s21(start, stop, res, num_points):
    freqs = np.linspace(start, stop, num_points)
    VISA = "GPIB0::16::INSTR"
    try:
        rm = pyvisa.ResourceManager() 
        try_visa = 1
    except:
        try_visa = 0
    if try_visa:
        vna = rm.open_resource(VISA)
        vna.read_termination = '\n'
        vna.write_termination = '\n'
        if vna.query('*IDN?'):
            print("Connected to VNA:", vna.query('*IDN?'))
            # Configure VNA parameters (example)
            vna.write(f'CALC1:PAR:S21:DEF')  # Define S21 measurement
            vna.write(f'SENS1:FREQ:START {start}GHz')
            vna.write(f'SENS1:FREQ:STOP {stop}GHz')
            vna.write(f'SENS1:BAND:RES {res}kHz')
            vna.write(f'SENS1:SWE:POIN {num_points}')  # Set the number of sweep points
            vna.write(f'SENS1:PAR:S21:FORM MLOG')

            # Trigger a single sweep
            vna.write(f'INIT:IMM')

            # Wait for the measurement to complete (you may need to adjust the wait time)
            vna.query(f'*OPC?')

            # Retrieve S21 data
            # s21_data = np.abs(vna.query_ascii_values(f'CALC1:DATA? SDATA')[:num_points])
            s21 = np.array(vna.query_ascii_values(f'CALC1:DATA? FDATA'))

            # frequency_data = vna.query_ascii_values('CALC1:X?')  # Query the frequency values

            # frequency_data = vna.query_ascii_values('SENS1:FREQ:DATA?')
            vna.close()
        else:
            s21 = np.zeros(len(freqs))
    else:
        s21 = np.zeros(len(freqs))
    return freqs, s21


def plot_s21(freqs, s21):
    fig, ax = plt.subplots()
    _ = ax.plot(freqs, s21)
    ax.set_xlabel('F [GHz]')
    ax.set_ylabel('|$S_{21}$| [dB]')
    ax.set_title('S21 from VNA')
    ax.grid(which='Major')
    plt.show()