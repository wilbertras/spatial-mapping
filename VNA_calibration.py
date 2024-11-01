from functions import *
import time
import pyvisa


if __name__ == "__main__":
    vna = connect2vi("GPIB0::16::INSTR", timeout=3000000)
    # weinschell = connect2vi("GPIB0::10::INSTR", timeout=300000)
    init_vna(vna, calibfile='D:\KIDS\KIDs.csa')
    calibrate_vna(vna)
    vna.close()
