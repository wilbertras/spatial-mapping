import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm


f0 = 4e9
oct = 1
Q = 50e3
delta = 20  # fwhm

# oct1: 1.6e3
# oct2: 3.2e3
# oct3: 5.2e3
# oct4: 6.5e3

def carlo(sigma, df, plot=False):
    a = np.random.normal(f0, sigma*f0, N)
    b = np.random.normal(f0+df, sigma*(f0+df), N)
    Pr = np.sum(b - a > delta * a/Q) / N
    return Pr


def semi_analytic(sigma, df):
    f = np.linspace(f0-nr_sigma*sigma*f0, f0+nr_sigma*sigma*f0, N)
    Q_func = norm.sf(f+delta*f0/Q, f0+df, sigma*(f0+df))
    pdf = norm.pdf(f, f0, sigma*f0)
    Pr = np.sum(pdf*Q_func*(2*nr_sigma*sigma*f0)/N)
    return Pr

# N = 10000
# df = 2e6
# nr_sigma = 10
# sigmas = np.logspace(-5, 0)
# Pr_carlo = []
# Pr_semi = []
# for sigma in sigmas:
#     Pr_semi.append(semi_analytic(sigma, df))
#     Pr_carlo.append(carlo(sigma, df))
# fig, ax = plt.subplots()
# ax.semilogx(sigmas, Pr_semi, label='semi-analytic') 
# ax.semilogx(sigmas, Pr_carlo, label='carlo') 
# ax.legend()
# plt.show()

def res_yield(sigma, nr_f0s):
    f_d = rel_freqs(nr_f0s)
    nr_fs = len(f_d)
    f_m = np.random.normal(f_d, sigma*f_d, (N, nr_f0s))
    nr = np.sum((f_m[:, 1:] - f_m[:, :-1]) > (delta*f_d[1:]/Q).T, axis=1)
    return np.mean(nr / (nr_fs-1), axis=0), np.std(nr / (nr_fs-1), axis=0)

def rel_freqs(nr):
    f_end = 2**oct * f0
    x = Q*(2**(oct/nr)-1)
    powers = np.arange(nr)
    freqs = f0 * (1+x/Q)**powers
    return freqs


N = 100
sigmas = np.logspace(-5, 0)
nrs = np.logspace(1, 4, 4, endpoint=True)
fig, ax = plt.subplots()
for nr in nrs:
    yields = []
    for sigma in sigmas:
        yields.append(res_yield(sigma, int(nr)))
    yields = np.stack(yields)    
    ax.semilogx(sigmas, yields[:,0], label=f'nr_f0s={nr}')
    ax.fill_between(sigmas, yields[:,0]-yields[:,1], yields[:,0]+yields[:,1], alpha=.5)
ax.legend()
ax.set_xlabel(['# kids'])
ax.set_ylabel(['expected yield'])
ax.set_title('yield vs sigma for various nr_f0s')
plt.show()

# N = 10
# nrs = np.logspace(1, 5, 40, endpoint=True)
# sigmas = np.logspace(-5, -2, 4, endpoint=True)
# yields2 = []
# fig, ax = plt.subplots()
# for sigma in sigmas:
#     yields2 = []
#     for nr in nrs:
#         yields2.append(res_yield(sigma, int(nr)))
#     yields2 = np.stack(yields2)    
#     ax.semilogx(nrs, yields2[:,0], label=f'$\sigma/f$=1E$^%d$' % (np.log10(sigma)))
#     ax.fill_between(nrs, yields2[:,0]-yields2[:,1], yields2[:,0]+yields2[:,1], alpha=.5)
# ax.axvline(oct/np.log2(delta/Q+1), c='k', ls='--', label='$oct/(\mathrm{log}_2(\lambda/Q+1))$')
# ax.set_xlabel('# kids')
# ax.set_ylabel('expected yield')
# ax.set_ylim(0, 1.01)
# ax.set_xlim(nrs[0], nrs[-1])
# ax.legend()
# ax.set_title('yield vs nr kids for various sigma')
# plt.show()
