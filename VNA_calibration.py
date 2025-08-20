from functions import *
import time
import pyvisa


if __name__ == "__main__":
    vna = connect2vi("GPIB0::16::INSTR", timeout=3000000)
    file='D:\KIDS\KIDs.csa'
    init_vna(vna, calibfile=file)

    vna_scan(vna, 3, 9, 3201, -110, 10000, 1)

    calibrate_vna(vna)

    if vna.query('*OPC?'):
        response = vna.query_ascii_values('CALC1:DATA? FDATA;')
        s21 = np.array(response)
    # freqs = np.array(vna.query_ascii_values('CALC1:MEAS1:MIX:XAX?;'))
    init_vna(vna, calibfile='D:\KIDS\KIDs.csa')
    vna_scan(vna, 3, 9, 3201, -110, 10000, 1)
    vna.write(f'MMEMORY:STOR:CSAR "{file}";')  # Opslaan met dezelfde bestandsnaam

    fig, ax = plt.subplots()
    ax.plot(s21)

    vna.close()
    plt.show()