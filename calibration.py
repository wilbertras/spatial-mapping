import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d




def read_txt(file_path):
    # Define the delimiter used to separate values in each line
    delimiter = '\t'  # Change this to the appropriate delimiter

    # Initialize a list to store the extracted data
    data = []

    # Open the file and read its contents line by line
    with open(file_path, 'r') as file:
        for line in file:
            # Split the line into fields using the specified delimiter
            fields = line.strip().split(delimiter)
            data.append(fields)
    np_data = np.array(data, dtype=np.float64)
    return np_data


def make_npy():
    BF_calib = read_txt(r'calibration files\BlueforsOutput2020_4-8GHz.txt')
    RF_calib = read_txt(r'calibration files\RFbox_Entropy.txt')
    np.save('calibration files/BF_calib.npy', BF_calib)
    np.save('calibration files/RF_calib.npy', RF_calib)


def plot_calib_files():
    BF_calib = np.load('calibration files/BF_calib.npy')
    RF_calib = np.load('calibration files/RF_calib.npy')
    fig, ax = plt.subplots()
    ax.plot(BF_calib[:, 0], BF_calib[:, 1], label='BF')
    ax.plot(RF_calib[:, 0], RF_calib[:, 1], label='RFcol1')
    ax.plot(RF_calib[:, 0], RF_calib[:, 2], label='RFcol2')
    bf = interp1d(BF_calib[:, 0], BF_calib[:, 1])
    rf = interp1d(RF_calib[:, 0], RF_calib[:, 1])
    ax.set_title('Calibration curves')
    ax.set_xlabel('F [GHz]')
    ax.set_ylabel('P [dBm]')
    ax.legend()
    plt.show()


## Uncomment the two functions below to generate and plot the calbration files
# make_npy()
plot_calib_files()
