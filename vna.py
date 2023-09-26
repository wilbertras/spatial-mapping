import pyvisa
import matplotlib.pyplot as plt
import numpy as np


def get_s21(start, stop, num_points, timeout=3000):
    freqs = np.linspace(start, stop, num_points)
    VISA = "GPIB0::16::INSTR"
    try:
        rm = pyvisa.ResourceManager() 
        try_visa = 1
    except:
        try_visa = 0
        print('NOT connected to VNA')
    if try_visa:
        vna = rm.open_resource(VISA)
        vna.timeout = timeout
        vna.read_termination = '\n'
        vna.write_termination = '\n'
        if vna.query('*IDN?'):
            print("Connected to VNA:", vna.query('*IDN?'))
            # Configure VNA parameters (example)
            vna.write(f'CALC1:PAR:S21:DEF')  # Define S21 measurement
            vna.write(f'SENS1:FREQ:START {start}GHz')
            vna.write(f'SENS1:FREQ:STOP {stop}GHz')
            vna.write(f'SENS1:SWE:POIN {num_points}')  # Set the number of sweep points
            vna.write(f'SENS1:PAR:S21:FORM MLOG')
            # vna.write(f'SENS1:POW:ATT AREC,{att}')

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


def plot_s21(freqs, s21):
    fig, ax = plt.subplots()
    _ = ax.plot(freqs, s21)
    ax.set_xlabel('F [GHz]')
    ax.set_ylabel('|$S_{21}$| [dB]')
    ax.set_title('S21 from VNA')
    ax.grid(which='Major')
    plt.show()