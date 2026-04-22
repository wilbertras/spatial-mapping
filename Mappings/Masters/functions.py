import numpy as np
from itertools import product
from scipy.special import erf
from scipy.signal import welch
from sklearn.neighbors import KernelDensity


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
    spacing = 2**(oct/(nr_kids-1))
    powers = np.arange(nr_kids)
    f0s = f0 * (spacing)**powers
    return f0s


def yld(Q, chi, sigma, Delta, fs=None):
    if Delta is None and fs is not None:
        f0 = fs[0]
        f1 = fs[-1]
        oct = np.log2(f1/f0)
        N = len(fs)
        Delta = 2**(oct/(N-1)) - 1
    else:
        pass
    return p0(Q, chi, sigma, Delta)


def p0(Q, chi, sigma, Delta):
    w = 1/Q
    n = np.arange(1,10000)
    return np.prod(1-(erf((n*Delta + chi*w)/(np.sqrt(2)*sigma)) - erf((n*Delta - chi*w)/(np.sqrt(2)*sigma)))/2)**2


def p1(chi, sigma, Delta):
    n = np.arange(1,10000)
    return np.prod(1-(erf((n*Delta + chi)/(np.sqrt(2)*sigma)) - erf((n*Delta - chi)/(np.sqrt(2)*sigma)))/2)**2


def p0_numeric(oct, Q, nr_kids, sigma, chi, nr_iter=10, fs=None):
    if fs is not None:
        f_d = fs
        nr_kids = len(f_d)
    else:
        f0 = 1
        f1 = f0 * 2**oct
        f_d = rel_freqs(f0, f1, nr_kids)
    f_m = np.random.normal(f_d, sigma*f_d, (nr_iter, nr_kids))
    f_m = np.sort(f_m, axis=1)
    df_left = np.absolute(f_m[:, 1:] - f_m[:, :-1]) / f_m[:, :-1]
    df_right = np.absolute(f_m[:, :-1] - f_m[:, 1:]) / f_m[:, 1:]
    too_close = np.zeros((nr_iter, nr_kids))
    too_close_left = df_left < (chi / Q)
    too_close_right = df_right < (chi / Q)
    too_close[:, 1:] += too_close_left
    too_close[:, :-1] += too_close_right 
    yields = np.sum(too_close==0, axis=1) / nr_kids
    return np.mean(yields, axis=0)


def yld_map(Ns, sigmas, Q, chi, oct, type='analytic'):
    map = np.zeros((len(sigmas), len(Ns)))
    for i, N in enumerate(Ns):
        if type=='analytic':
            Delta = 2**(oct/(N-1)) - 1
            for j, sigma in enumerate(sigmas):
                map[i, j] = yld(Q, chi, sigma, Delta)
        elif type=='numeric':
            for j, sigma in enumerate(sigmas):
                for i, N in enumerate(Ns):
                    map[i, j] = p0_numeric(oct, Q, N, sigma, chi)
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
    return phase, logsmooth(fxx, nxx)


# def comp_yield(f0s, Q, threshold):
#     diffs = np.diff(f0s)
#     fwhms = f0s/Q
#     rel_diffs = diffs/fwhms[:-1]
#     lo_id = np.argmax(diffs)
#     too_close = rel_diffs<threshold
#     good = np.ones(f0s.shape)
#     good[1:] -= too_close
#     good[:-1] -= too_close
#     spaced = np.sum(good==True)
#     total = len(f0s)
#     return rel_diffs, lo_id, spaced/total


# def plot_spacings(f0s, Q, threshold, kernel=False, ax=None, title='', c='b', bins='auto'):
#     rel_diffs, lo_id, yld = comp_yield(f0s, Q, threshold)
#     if not ax:
#         fig, ax = plt.subplots()
#     else:
#         lower = rel_diffs[:lo_id]
#         upper = rel_diffs[lo_id+1:]
#         _ = ax.hist(lower, bins=bins, alpha=.5, label='lower band', facecolor=c, edgecolor=c, density=True)
#         _ = ax.hist(upper, bins=bins, alpha=.5, label='upper band', facecolor='w', edgecolor=c, density=True)
#         if kernel:
#             x, density, std = kde(lower, bw=1)
#             ax.plot(x, density, c='k', ls='--')
#             print(std/np.sqrt(2)/Q)
#             x, density, std = kde(upper, bw=1)
#             ax.plot(x, density, c='k', ls='--')
#             print(std/np.sqrt(2)/Q)
#         # ax.axvline(threshold, c='r', ls='--', lw=1, label='minimal spacing')
#         # ax.set_title(title + 'usable yield = %.1f%%' % (yld*1e2), fontsize=10)
#         ax.set_ylabel('$\#$')
#         ax.legend()


def kde(y, ylim=[], bw=.2):
    y = y[~np.isnan(y)].reshape(-1,1)
    if ylim:
        y = y[(y >= ylim[0]) & (y <= ylim[1])].reshape(-1,1)
    kde = KernelDensity(kernel='gaussian', bandwidth=bw)
    kde.fit(y)

    # Evaluate the density
    x = np.linspace(y.min(), y.max(), 1000).reshape(-1,1)
    log_density = kde.score_samples(x)
    density = np.exp(log_density)  # Convert log density to actual density
    # density = density * (y.max() - y.min()) + y.min()  # Convert density back to original scale
    half_height = density.max() / 2
    indices = np.where(density >= half_height)[0]
    width = x[indices[-1]] - x[indices[0]]
    return x, density, width[0]