#!/usr/bin/env python

import copy

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


def parse(fname):
    nY, nX = np.loadtxt(fname, max_rows=1, skiprows=3, unpack=True, dtype=int)
    data = np.loadtxt(fname, skiprows=4)
    data = np.reshape(data, ((nX, nY, -1)))
    Yran = np.array(data[0, :, 0])
    Xran = np.array(data[:, 0, 1])
    data = np.swapaxes(data, 0, 1)
    return data, Yran, Xran


with open('eos_plotter.dat') as f:
    title = f.readline().strip()
    xlabel = f.readline().strip()
    ylabel = f.readline().strip()

eosDT, Yran, Xran = parse('eos_plotter.dat')

# set up plot and labels
fig, ax = plt.subplots(figsize=(5, 4))
ax.set_title(title)
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_xlim(Xran.min(), Xran.max())
ax.set_ylim(Yran.min(), Yran.max())

# set up color map
cmap = copy.copy(mpl.cm.get_cmap("viridis"))
cmap.set_over('white')
cmap.set_under('black')
cmap.set_bad('grey')

# set color bar limits
# None will auto-set limits
cbar_min = None
cbar_max = None

pcol = ax.pcolormesh(Xran, Yran, eosDT[..., 2], shading='nearest', cmap=cmap,
                     vmin=cbar_min, vmax=cbar_max)
pcol.set_edgecolor('face')
cax = fig.colorbar(pcol, extend='both')
cax.set_label('')

# save figure
fig.savefig('eos_plotter.png', dpi=150)

plt.show()
