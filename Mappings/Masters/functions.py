import numpy as np
from itertools import product
from scipy.special import erf
from scipy.signal import welch


def build_design_matrix(x, y, degree):
    """
    Build polynomial design matrix for 2D total-degree polynomial.
    x, y are 1D arrays of valid (non-NaN) data.
    """
    terms = []
    powers = []
    for i, j in product(range(degree+1), repeat=2):
        if i + j <= degree:
            terms.append((x**i) * (y**j))
            powers.append((i, j))

    return np.vstack(terms).T, powers


def fit_poly2d_nan_safe(X, Y, Z, degree):
    """
    Fit polynomial surface z = f(x,y)
    X, Y, Z are 2D arrays of equal shape, Z may contain NaNs.
    """
    # Flatten the arrays
    x = X.ravel()
    y = Y.ravel()
    z = Z.ravel()

    # Mask out NaN values
    mask = ~np.isnan(z)
    x_valid = x[mask]
    y_valid = y[mask]
    z_valid = z[mask]

    # Build design matrix using only valid points
    M, powers = build_design_matrix(x_valid, y_valid, degree)

    # Least squares solution
    coeffs, *_ = np.linalg.lstsq(M, z_valid, rcond=None)

    return coeffs, powers


def eval_poly2d_on_grid(X, Y, coeffs, powers):
    """
    Evaluate polynomial surface back on a 2D grid.
    NaNs in Z do not matter for evaluation.
    """
    Z_fit = np.zeros_like(X, dtype=float)
    for c, (i, j) in zip(coeffs, powers):
        Z_fit += c * (X**i) * (Y**j)
    return Z_fit


def rel_freqs(f0, f1, nr_kids):
    oct = np.log2(f1/f0)
    spacing = 2**(oct/nr_kids)
    powers = np.arange(nr_kids)
    f0s = f0 * (spacing)**powers
    return f0s


def yld(Q, chi, sigma, Delta, fs=None):
    if Delta is None and fs is not None:
        f0 = fs[0]
        f1 = fs[-1]
        oct = np.log2(f1/f0)
        N = len(fs)
        Delta = 1-2**(oct/N)
    else:
        pass
    return p0(Q, chi, sigma, Delta)


def p0(Q, chi, sigma, Delta):
    w = 1/Q
    n = np.arange(1,10000)
    return np.prod(1-(erf((n*Delta + chi*w)/(np.sqrt(2)*sigma)) - erf((n*Delta - chi*w)/(np.sqrt(2)*sigma)))/2)**2


def yld_map(Ns, sigmas, Q, chi, oct):
    map = np.zeros((len(sigmas), len(Ns)))
    for i, N in enumerate(Ns):
        Delta = 1-2**(oct/N)
        for j, sigma in enumerate(sigmas):
            map[i, j] = yld(Q, chi, sigma, Delta)
    return map

def bin2mat(file_path):
    data = np.fromfile(file_path, dtype='>f8', count=-1)
    data = data.reshape((-1, 2))

    I = data[:, 0]
    Q = data[:, 1]

    # From I and Q data to Radius/Magnitude and Phase
    r = np.sqrt(I**2 + Q**2)
    I /= np.mean(r) # Normalize I to 1
    Q /= np.mean(r) # Normalize Q to 1
    R = np.sqrt(I**2 + Q**2)

    P = np.arctan2(Q, I) 
    P = np.pi - P % (2 * np.pi) # Convert phase to be taken from the negative I axis
    return R, P


def logsmooth(fs, sxx):
    # Number of points per decade
    points_per_decade = 10

    # Define logarithmic bins
    log_min = np.log10(fs[1:].min())
    log_max = np.log10(fs.max())
    log_bins = np.logspace(log_min, log_max, int((log_max - log_min) * points_per_decade))

    # Bin the data
    binned_frequencies = []
    binned_powers = []

    for i in range(len(log_bins) - 1):
        # Find indices of frequencies within the current bin
        bin_indices = (fs >= log_bins[i]) & (fs < log_bins[i + 1])
        
        # Calculate the mean frequency and power spectrum value for the bin
        if np.any(bin_indices):
            binned_frequencies.append(np.mean(fs[bin_indices]))
            binned_powers.append(np.mean(sxx[bin_indices]))

    # Convert to arrays
    binned_frequencies = np.array(binned_frequencies)
    binned_powers = np.array(binned_powers)
    return binned_frequencies, binned_powers


def get_noise_psd(file, wl=256, fs=1e6):
    _, phase = bin2mat(file)
    phase -= np.mean(phase)
    fxx, nxx = welch(phase, fs=int(fs), nperseg=wl, window='flattop', return_onesided=True)
    return logsmooth(fxx, nxx)


