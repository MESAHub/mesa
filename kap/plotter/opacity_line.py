#!/usr/bin/env python

import copy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import interp1d

#xpts = np.array([0,0.1,0.2,0.35,0.5,0.7,0.8,0.9,0.95,1.0])
xpts = np.array([0.1,0.2,0.35,0.5,0.7,0.8,0.9,0.95])

def parse(fname):
    nY, nX = np.loadtxt(fname, max_rows=1, skiprows=3, unpack=True, dtype=int)
    data = np.loadtxt(fname, skiprows=4)
    data = np.reshape(data, ((nX, nY, -1)))
    Yran = np.array(data[0,:,0])
    Xran = np.array(data[:,0,1])
    data = np.swapaxes(data, 0, 1)
    return data, Yran, Xran, nX, nY

#files = ['opacity_linear_bigrange.dat','opacity_cubic_bigrange.dat']
files = ['opacity_cubic_bigrange.dat']

with open(files[0]) as f:
    title = f.readline().strip()
    xlabel = f.readline().strip()
    ylabel = f.readline().strip()


# set up plot and labels
fig, ax = plt.subplots()
ax.set_xlabel(xlabel)
ax.set_ylabel(title)
    
for infile in files:
    kapDT, Yran, Xran, nX, nY = parse(infile)
    
    logkap_logTmid = kapDT[int(nY/2),:,2]
    logkap_logTmid_m1 = kapDT[int(nY/2)-10,:,2]
    logkap_logTmid_p1 = kapDT[int(nY/2)+10,:,2]

    logTmid_m1 = kapDT[int(nY/2)-10,0,0]
    logTmid = kapDT[int(nY/2),0,0]
    logTmid_p1 = kapDT[int(nY/2)+10,0,0]
    
    # check temperature
    print("logT =",logTmid_m1,logTmid,logTmid_p1)

    ax.plot(Xran,logkap_logTmid_m1,label='logT = {:.2f}'.format(logTmid_m1))
    ax.plot(Xran,logkap_logTmid,label='logT = {:.2f}'.format(logTmid))
    ax.plot(Xran,logkap_logTmid_p1,label='logT = {:.2f}'.format(logTmid_p1))
    #ax.plot(Xran,logkap_logTmid_p1,ls='--')
    ax.set_xlim(Xran.min(), Xran.max())

    finterp = interp1d(Xran,logkap_logTmid)
    ax.scatter(xpts,finterp(xpts))

#ax.set_xlim(0.8,1.0)
#ax.set_ylim(-1,-0.95)


ax.legend()

# save figure
fig.savefig('opacity_line.png', dpi=300)

plt.show()
