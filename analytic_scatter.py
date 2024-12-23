import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm


df = 4e6
f0 = 6e9
f1 = f0 + df
Q = 50e3
delta = 20
N = 1000
Pr_carlo = []
Pr_semi = []
nr_sigma = 10
yields = []

def carlo(sigma, plot=False):
    f = np.linspace(f0-nr_sigma*sigma*f0, f0+nr_sigma*sigma*f0, N)
    a = np.random.normal(f0, sigma*f0, N)
    b = np.random.normal(f1, sigma*f1, N)
    Pr = np.sum(b - a > delta) / N
    if plot:
        fig, ax = plt.subplots()
        ax.hist(a, bins='auto')
        ax.hist(b, bins='auto')
        ax.set_title(f'Probability is {Pr}')
        plt.show()
    return Pr


def semi_analytic(sigma):
    f = np.linspace(f0-nr_sigma*sigma*f0, f0+nr_sigma*sigma*f0, N)
    Q = norm.sf(f+delta, f1, sigma*f1)
    p = norm.pdf(f, f0, sigma*f0)
    Pr = np.sum(p*Q*(2*nr_sigma*sigma*f0)/N)
    return Pr

def res_yield(sigma, nr_f0s):
    f0s = np.linspace(4e9+df, 6e9-df, nr_f0s, endpoint=True)
    nr_f0s = len(f0s)
    fs = np.random.normal(f0s, sigma*f0s, (N, nr_f0s))
    nr = np.sum((fs[:, 1:] - fs[:, :-1]) > (delta*f0s[1:]/Q).T, axis=1)
    # nr -= np.sum(fs < 4e9, axis=1)
    # nr -= np.sum(fs > 6e9, axis=1)
    return np.mean(nr / (nr_f0s-1), axis=0)


# sigmas = np.logspace(-5, 0)
# for sigma in sigmas:
#     Pr_semi.append(semi_analytic(sigma))
#     Pr_carlo.append(carlo(sigma))
# fig, ax = plt.subplots()
# ax.semilogx(sigmas, Pr_semi, label='semi-analytic') 
# ax.semilogx(sigmas, Pr_carlo, label='carlo') 
# ax.legend()
# plt.show()


sigmas = np.logspace(-5, 0)
nrs = np.logspace(1, 4, 4, endpoint=True)
fig, ax = plt.subplots()
for nr in nrs:
    yields = []
    for sigma in sigmas:
        yields.append(res_yield(sigma, int(nr)))
    ax.semilogx(sigmas, yields, label=f'nr_f0s={nr}')
ax.legend()
ax.set_title('yield vs sigma for various nr_f0s')
plt.show()


# nrs = np.logspace(1, 4, 30, endpoint=True)
# yields2 = []
# fig, ax = plt.subplots()
# for sigma in np.logspace(-5, -2, 4, endpoint=True):
#     yields2 = []
#     for nr in nrs:
#         yields.append(res_yield2(sigma, int(nr)))
#     ax.semilogx(nrs, yields2, label=f'sigma={sigma}')
# ax.legend()
# ax.set_title('yield vs nr kids for various sigma')
# plt.show()