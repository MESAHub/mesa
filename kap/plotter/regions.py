#!/usr/bin/env python

import copy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors


plt.style.use('mesa_kap_regions.mplstyle')

def parse(fname):
    nY, nX = np.loadtxt(fname, max_rows=1, skiprows=3, unpack=True, dtype=int)
    data = np.loadtxt(fname, skiprows=4)
    data = np.reshape(data, ((nX, nY, -1)))
    Yran = data[0,:,0]
    Xran = data[:,0,1]
    data = np.swapaxes(data, 0, 1)
    return data, Yran, Xran


with open('kap_plotter.dat') as f:
    title = f.readline().strip()
    xlabel = f.readline().strip()
    ylabel = f.readline().strip()

# overwrite with fancier labels
xlabel = r'$\log(\rho/{\rm g\,cm^{-3}})$'
ylabel = r'$\log(T/{\rm K})$'
title = r'MESA kap Regions'

kapDT, Yran, Xran = parse('kap_plotter.dat')

apjcolwidth=3.38
# set up plot and labels
#fig, ax = plt.subplots(figsize=(apjcolwidth,apjcolwidth*4./5.)) # for paper figures
fig, ax = plt.subplots(figsize=(4,3)) # for website pngs
ax.set_title(title)
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_xlim(Xran.min(), Xran.max())
ax.set_ylim(Yran.min(), Yran.max())

# set up color map (slightly customized to make Skye blue)
my_colors = np.array(mpl.colormaps["Accent"].colors) # array so that entries are editable
#tmp = my_colors[4].copy()
#my_colors[4] = my_colors[3]
#my_colors[5] = tmp
cmap = colors.ListedColormap(my_colors)
bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
norm = colors.BoundaryNorm(bounds, cmap.N)

pcol = ax.pcolormesh(Xran, Yran, kapDT[...,2], shading='nearest', cmap=cmap, norm=norm, rasterized=True)
pcol.set_edgecolor('face')
cax = fig.colorbar(pcol, ticks=[0, 1, 2, 3, 4,5 ])
cax.set_label('')
cax.ax.tick_params(labelsize='x-small')
cax.ax.minorticks_off()
cax.ax.set_yticklabels(['blend', 'Low T', 'Type 2', 'High T', 'Compton', 'No tables'])

# save figure
fig.savefig('kap_regions.png', dpi=150)

plt.show()
