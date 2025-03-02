#!/usr/bin/env python

import copy
import numpy as np
import matplotlib.pyplot as plt


def parse(fname):
    nY, nX = np.loadtxt(fname, max_rows=1, skiprows=3, unpack=True, dtype=int)
    data = np.loadtxt(fname, skiprows=4)
    data = np.reshape(data, ((nX, nY, -1)))
    Yran = np.array(data[0, :, 0])
    Xran = np.array(data[:, 0, 1])
    data = np.swapaxes(data, 0, 1)
    return data, Yran, Xran


with open("kap_plotter.dat") as f:
    title = f.readline().strip()
    xlabel = f.readline().strip()
    ylabel = f.readline().strip()

kapDT, Yran, Xran = parse("kap_plotter.dat")

# set up plot and labels
fig, ax = plt.subplots(figsize=(4, 3))
ax.set_title(title)
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_xlim(Xran.min(), Xran.max())
ax.set_ylim(Yran.min(), Yran.max())

# set up color map
cmap = copy.copy(plt.get_cmap("coolwarm"))
cmap.set_over("tab:red")
cmap.set_under("black")

# set color bar limits
# None will auto-set limits
cbar_min = -14
cbar_max = 2

pcol = ax.pcolormesh(
    Xran,
    Yran,
    kapDT[..., 2],
    shading="nearest",
    cmap=cmap,
    vmin=cbar_min,
    vmax=cbar_max,
)
pcol.set_edgecolor("face")
cax = fig.colorbar(pcol, extend="both")
cax.set_label("")

# save figure
fig.savefig("kap_plotter.png", dpi=300)

plt.show()
