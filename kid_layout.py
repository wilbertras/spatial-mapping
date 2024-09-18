import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
import matplotlib as mpl

if __name__ == "__main__":
    N = 32
    
    filename = './kid id boards/kid_id_board_%dx%d' % (N, N)

    ids = np.zeros((N, N))

    for i in range(N):
        id = np.arange(0, N**2, N, dtype=int) + i
        ids[i] = np.roll(id, i*(int(np.ceil(N/3))))

    fig, axes = plt.subplot_mosaic('a', figsize=(6, 6), constrained_layout=True)
    ax = axes['a']
    ax.imshow(ids, origin='lower', cmap='viridis')
    for i in range(N):
        for j in range(N):
            idx = ids[i, j]
            if ~np.isnan(idx):
                color = 'w'
                text = ax.text(j, i, '%d' % idx, ha="center", va="center", color=color, fontsize='xx-small') 
    ax.set_xticks(np.arange(0, N, 1))
    ax.set_yticks(np.arange(0, N, 1))
    ax.set_xticklabels(np.arange(1, N+1, 1))
    ax.set_yticklabels(np.arange(1, N+1, 1))
    ax.set_xticks(np.arange(0.5, N, 1), minor=True)
    ax.set_yticks(np.arange(0.5, N, 1), minor=True)
    ax.set_xlabel('$\it x$ $[px]$')
    ax.set_ylabel('$\it y$ $[px]$')
    ax.grid(True, which='minor', color='None', linestyle='-', linewidth=0.5)
    ax.grid(False, which='major')

    # np.save(filename + '.npy', ids)
    plt.savefig(filename + '.png')
    plt.show()