Calibration procedure as understood from the various labview files
We want to input desired power at KID and compute what power that translated to at VNA
KIDpower + KID_cryoOUT + GainRFbox(freq) - Weinschell = P@VNAin
# KIDpower + KID_cryoOUTp + GainRFbox(freq)-P@VNAin=W

1. KIDpower obtained from BleuforsOutput2022_4-8GHz.txt 
2. Interpolate at F0 to get the correct value
3. KID_cryoOUT is obtained from RFbox_Entropy.txt 
4. Interpolate at F0 to get the correct value
5. Add: PcryoOUt=KIDpower + KID_cryoOUT
6. Substract: P@VNAin = PcryoOUt --2
7. Check if PryoOUt <= 62
8. Check if PryoOUt > 0
9. Clip PcryoOUt to even number, Weinschell=round(PryoOUt / 2)*2, write to VNA
10. Power into VNA@F0=KIDpower + KID_cryoOUT - Weinschell

