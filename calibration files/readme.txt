Calibration procedure as understood from the various labview files
We want to input desired power at KID and compute what power that translated to at VNA
# KIDpower + KID_cryoOUT + GainRFbox(freq) - Weinschell = P@VNAin
# KIDpower + KID_cryoOUTp + GainRFbox(freq)-P@VNAin=W

0. Input KIDpower
1. KID_cryoOUt obtained from BleuforsOutput2022_4-8GHz.txt 
2. Interpolate at F0 to get the correct value
3. GainRFbox is obtained from RFbox_Entropy.txt 
4. Interpolate at F0 to get the correct value
5. Add: PcryoOUt=KIDpower + KID_cryoOUT + GainRFbox
6. Substract: P@VNAin = PcryoOUt --2
7. Check if PcryoOUt <= 62, else: PcryoOUt=62, give error
8. Check if PcryoOUt > 2, else: 
9. Clip PcryoOUt to even number, Weinschell=round(PryoOUt / 2)*2, write to VNA
10. Power into VNA@F0=KIDpower + KID_cryoOUT + GainRFbox - Weinschell

