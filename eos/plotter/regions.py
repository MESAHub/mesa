#!/usr/bin/env python

import copy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def parse(fname):
    nY, nX = np.loadtxt(fname, max_rows=1, skiprows=3, unpack=True, dtype=int)
    data = np.loadtxt(fname, skiprows=4)
    data = np.reshape(data, ((nX, nY, -1)))
    Yran = data[0,:,0]
    Xran = data[:,0,1]
    data = np.swapaxes(data, 0, 1)
    return data, Yran, Xran


with open('eos_plotter.dat') as f:
    title = f.readline().strip()
    xlabel = f.readline().strip()
    ylabel = f.readline().strip()

# overwrite with fancier labels
xlabel = r'$\log_{10}(\rho/{\rm g\,cm^{-3}})$'
ylabel = r'$\log_{10}(T/{\rm K})$'
title = r'MESA EOS Regions ($X=0.7$, $Z=0.02$)'

eosDT, Yran, Xran = parse('eos_plotter.dat')

# set up plot and labels
fig, ax = plt.subplots(figsize=(5,4))
ax.set_title(title)
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_xlim(Xran.min(), Xran.max())
ax.set_ylim(Yran.min(), Yran.max())

# set up color map
cmap = mpl.cm.get_cmap("Accent")
bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

pcol = ax.pcolormesh(Xran, Yran, eosDT[...,2], shading='nearest', cmap=cmap, norm=norm)
pcol.set_edgecolor('face')
cax = fig.colorbar(pcol, ticks=[0, 1, 2, 3, 4, 5, 6, 7])
cax.set_label('')
cax.ax.set_yticklabels(['blend', 'HELM', 'OPAL/SCVH', 'FreeEOS', 'PC', 'Skye', 'CMS', 'ideal'])

# save figure
fig.savefig('eos_regions.png', dpi=300)

# plt.show()
