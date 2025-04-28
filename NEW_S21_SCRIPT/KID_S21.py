import numpy as np
import glob
import pandas as pd
import re
import io
import os
import matplotlib.pyplot as plt

from Khalil import KhalilModel_magspace

def Khalil_func_magspace(f, f0, Ql, Qc_re, dw, a, b):
    return np.abs(1 - (Ql/Qc_re *(1 + 2j*Ql*dw/f0) / (1 + 2j * Ql * (f - f0) / f0)))*np.abs((b*(f-f0)+a))

# Function to add results dynamically

def def_Pint(Q,Qe,Pread):
        Pint = 10*np.log10((1/np.pi)*(Q**2/Qe)) + Pread
        return Pint

def create_result_pd():
    df_results = pd.DataFrame(columns=["KID", "Power", "Temperature",
                                       "f0", 'f0_std',
                                       "Ql", "Ql_std",
                                       "Qc", "Qc_std",
                                       "Qi", "Qi_std",
                                       "Pint"])
    return df_results
    
def add_result(df_results, kid, power, temperature, fit_result, Pint):
    
    new_entry = pd.DataFrame({"KID": [kid], "Power": [power], "Temperature": [temperature],
                              "f0": [fit_result.params['f0'].value],
                              "f0_std": [fit_result.params['f0'].stderr],
                              "Ql": [fit_result.params['Ql'].value],
                              "Ql_std": [fit_result.params['Ql'].stderr],
                              "Qc": [fit_result.params['Qc_re'].value],
                              "Qc_std": [fit_result.params['Qc_re'].stderr],
                              "Qi": [fit_result.params['Qi'].value],
                              "Qi_std": [fit_result.params['Qi'].stderr],
                              "Pint": [Pint],
                              "redchisqr": [fit_result.redchi]})
    
    df_results = pd.concat([df_results, new_entry], ignore_index=True)
    
    return df_results

def find_S21_files(path, kid='', pread=''):

    string = path + 'KID%s_%sdBm*.dat' % (kid,pread)    
    files = sorted(glob.glob(string)) # a sorted lis of filenames in path
    nr_files = len(files)
    
    kid = []
    for file in files:
        kid.append(int(re.findall(r"KID(\d+)_", file)[0]))
        # kid.append(int(re.findall(path + "KID(\d+)_(\d{2,3})dBm_", file)[0]))
        
    kids = np.unique(kid)
    
    return files, kids
    # optional print some info like path name and nr of KIDs
    
def loop_over_S21_files(path, kid=None, pread=None, plot=False):
    if not kid:
        kid = '*'
    else:
        kid = str(kid)
    if not pread:
        pread = '*'
    else:
        pread = str(pread)

    filenames, kids = find_S21_files(path, kid, pread)
    
    # output_dir = "temperature_data"
    # os.makedirs(output_dir, exist_ok=True)
    
    df_results = create_result_pd()

    for file_path in filenames:
        [kid, power] = re.findall("KID(\d+)_(\d+)dBm_", file_path)[0]
        kid_id = int(kid)
        Pread = -float(power)
        
        with open(file_path, 'r') as file:
            file_contents = file.readlines()

        # Extract data organized by temperature
        data_by_temperature = preprocess_file_fast(file_contents)
                
        all_temperatures = list(data_by_temperature.keys())
        plot_temperatures = all_temperatures[::5]
        
        for temperature, df in data_by_temperature.items():
            #print(f"Temperature: {temperature} K")

            # Extract Frequency and S21 data
            frequencies = df['Frequency'].values
            s21_values = df['dB'].values
            rad_values = df['Rad'].values
        
            if (np.min(frequencies) != np.max(frequencies)):
        
                # here we start fitting the data
                result = Fit_S21(frequencies, s21_values)

                S21_fit_line = result.eval()
                if plot:
                    if len(filenames) % plot == 0:
                        fig, ax = plt.subplots()
                        result.plot_fit(ax)
                Pint = def_Pint(result.params['Ql'].value, result.params['Qc_re'].value, Pread)
                
                df_results = add_result(df_results, kid_id, Pread, temperature, result, Pint)
                
                df['Mag'] = 10**((s21_values - np.mean(s21_values[0:100]))/20)
                df['Fit'] = S21_fit_line
                # fig, ax = plt.subplots()
                # ax.plot(frequencies, df['Mag'])
                # ax.plot(frequencies, df['Fit'])
                # file_name = "KID" + str(kid_id) + "_" + str(int(Pread)) + f"dBm_{temperature*1e3:.2f}mK.csv"
                # output_path = os.path.join(output_dir, file_name)

                # Save the DataFrame (Frequency, S21, Rad) to CSV
                # df.to_csv(output_path, index=False)
        
        # print("Saved data for KID " + str(kid_id) +", Pread -" +str(int(-Pread)) + ' dBm', end='\r')
        
    return df_results

            
def Fit_S21(f, S21_dB):
    
    # start with a very basic version that does nothing special....
    
    # The matlab versions applies a smotth here but i skip that for now
    S21_dB = S21_dB - np.mean(S21_dB[0:100]) # normalize the data
    
    S21_mag = 10**(S21_dB/20)
               
    # load the models
    Model_mag = KhalilModel_magspace
    model_mag = Model_mag(f, S21_mag)
    # a first estimate of the fit
    result_pre = model_mag.fit(S21_mag, f=f, params = model_mag.guess)
    # model_mag.plot_fit()
    
    
        
    return result_pre
    #print('fitted the data! yeah!')
    
    
def extract_data_by_temperature(file_contents):
    data_by_temperature = {}
    current_temperature = None
    current_data = []

    for line in file_contents:
        # Check if the line contains temperature information
        temp_match = re.match(r'Temperature in K:(\d+\.\d+)', line)
        if temp_match:
            if current_temperature is not None and current_data:
                # Store the current data block before moving to the next temperature
                data_by_temperature[current_temperature] = current_data
            
            # Start a new temperature section
            current_temperature = float(temp_match.group(1))
            current_data = []
        
        # Check if the line contains tab-separated data (GHz, dB, Rad)
        elif re.match(r'\d+\.\d+\t', line):
            current_data.append(line.strip().split('\t'))
    
    # Save the last block of data
    if current_temperature is not None and current_data:
        data_by_temperature[current_temperature] = current_data

    return data_by_temperature

def preprocess_file_fast(file_contents):
    data_by_temperature = {}
    current_temperature = None
    current_data = []

    # Process each line, reducing Python overhead
    for line in file_contents:
        # If temperature is found
        temp_match = re.match(r'Temperature in K:(\d+\.\d+)', line)
        if temp_match:
            if current_temperature is not None and current_data:
                # Convert current_data to a pandas DataFrame
                data_str = "\n".join(current_data)
                df = pd.read_csv(io.StringIO(data_str), sep='\t', header=None, names=['Frequency', 'dB', 'Rad'])
                data_by_temperature[current_temperature] = df
            
            # Set new temperature and reset current data
            current_temperature = float(temp_match.group(1))
            current_data = []
        
        # Collect tab-separated data (GHz, dB, Rad)
        elif re.match(r'\d+\.\d+\t', line):
            current_data.append(line.strip())
    
    # Save the last block of data
    if current_temperature is not None and current_data:
        data_str = "\n".join(current_data)
        df = pd.read_csv(io.StringIO(data_str), sep='\t', header=None, names=['Frequency', 'dB', 'Rad'])
        data_by_temperature[current_temperature] = df

    return data_by_temperature