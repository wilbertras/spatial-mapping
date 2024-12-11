import numpy as np
import matplotlib.pyplot as plt
plt.style.use('matplotlibrc')
import matplotlibcolors

std_df = []
std_df_f = []
f_d = np.linspace(4, 8, 1000)
kids = np.arange(len(f_d))
scatters = np.logspace(-6,0, 1000)
nr_swapped = 0
hline=scatters[-1]
min_df = 1-3
yields = []
for scatter in scatters:
    sf = np.random.normal(0, scatter, f_d.shape)
    a = 1.1
    f_m = (a+sf)*f_d
    
    df_f = (f_d - f_m) / f_d
    std_df_f.append(np.std(df_f))

    sort = np.argsort(f_m)
    swapped = np.sum(kids - sort != 0)
    ids = np.arange(0, 400, 10)
    removed = np.ones(f_m.shape, dtype=bool)
    # removed[ids] = 0
    f_m_sorted = f_m[sort][removed]
    if swapped and not nr_swapped:
        print(scatter, swapped)
        hline = scatter
        nr_swapped += 1
    spacing = f_m_sorted[:-1] - f_m_sorted[1:]
    df = (spacing) / 6 / np.sqrt(2)
    std_df.append(np.std(df))
    too_close = spacing < min_df
    left = np.hstack((too_close, 0))
    right = np.hstack((0, too_close))
    yields.append(np.sum((left+right)==0)/ len(df))
    # fig, ax = plt.subplots()
    # ax.hist(df, bins='auto', facecolor='o')
    # ax.set_xlabel('%.1e' % (scatter))

std_df_f = np.array(std_df_f)
std_df = np.array(std_df)
yields = np.array(yields)
fig, ax = plt.subplot_mosaic('acdb', constrained_layout=True)
ax['a'].scatter(f_d, f_m)
ax['a'].scatter(f_d[removed], f_m_sorted)
ax['a'].scatter(f_d, f_d)
ax['d'].hist(df_f, bins='auto', label='$\sigma=%.1e$' % (std_df_f[-1]))
ax['d'].hist(df, bins='auto', label='$\sigma=%.1e$' % (std_df[-1]))
ax['d'].legend()
ax['c'].loglog(scatters, std_df_f, label='$\sigma_{\\chi_n}$')
ax['c'].loglog(scatters, std_df, label='$\sigma_{\\Upsilon_i}$')
ax['b'].semilogx(scatters, yields*1e2, label='yields')
ax['c'].set_xlabel('In: fractional frequency scatter $\sigma$')
ax['c'].set_ylabel('Out: fractional frequency scatter $\sigma$')
# ax['c'].loglog(scatters, np.sqrt(2)*(5e-3/((f_d[0]+f_d[-1])/2)+1)*scatters, label='theory')
ax['c'].axvline(hline, color='k', linestyle='--', label='KIDs start to swap')
ax['c'].legend()
plt.show()