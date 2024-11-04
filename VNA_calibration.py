from functions import *
import time
import pyvisa


if __name__ == "__main__":
    vna = connect2vi("GPIB0::16::INSTR", timeout=3000000)

    init_vna(vna, calibfile='D:\KIDS\KIDs.csa')

    vna_scan(vna, 5, 2, 3201, -110, 1000, 1)

    calibrate_vna(vna)

    if vna.query('*OPC?'):
        response = vna.query_ascii_values('CALC1:DATA? FDATA;')
        s21 = np.array(response)

    fig, ax = plt.subplots()
    ax.plot(s21)
    
    vna.close()
    plt.show()