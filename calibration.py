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
    BFout = read_txt(r'calibration files\BlueforsOutput2020_4-8GHz.txt')
    BFin = read_txt(r'calibration files\BlueforsInput2020.txt')
    RFout = read_txt(r'calibration files\RFbox_Entropy.txt')

    np.save('calibration files/BFout.npy', BFout)
    np.save('calibration files/BFin.npy', BFin)
    np.save('calibration files/RFout.npy', RFout)


def plot_calib_files():
    BFout = np.load('calibration files/BFout.npy')
    BFin = np.load('calibration files/BFin.npy')
    RFout = np.load('calibration files/RFout.npy')
    fig, ax = plt.subplot_mosaic('abc', figsize=(9, 3), constrained_layout=True, sharex=True)
    ax['a'].plot(BFout[:, 0], BFout[:, 1], label='Bluefors out', c='tab:blue')
    ax['b'].plot(BFin[:, 0], BFin[:, 1], label='Bluefors in', c='tab:red')
    ax['c'].plot(RFout[:, 0], RFout[:, 1], label='Readout in', c='tab:green')
    # ax.plot(RFout[:, 0], RFout[:, 2], label='RFcol2')
    bf = interp1d(BFout[:, 0], BFout[:, 1])
    rf = interp1d(RFout[:, 0], RFout[:, 1])
    ax['a'].set_title('Bluefors out')
    ax['a'].set_xlabel('F [GHz]')
    ax['a'].set_ylabel('P [dBm]')
    ax['a'].set_ylim((-80, 80))
    ax['b'].set_xlabel('F [GHz]')
    ax['b'].set_title('Bluefors in')
    ax['b'].set_ylabel('P [dBm]')
    ax['b'].set_ylim((-80, 80))
    ax['c'].set_xlabel('F [GHz]')
    ax['c'].set_title('Readout gain')
    ax['c'].set_ylabel('P [dBm]')
    ax['c'].set_ylim((-80, 80))
    # ax.legend()
    plt.show()


## Uncomment the two functions below to generate and plot the calbration files
make_npy()
plot_calib_files()
