#!/usr/bin/env python

import copy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import interp1d

xpts = np.array([0.0,0.1,0.2,0.35,0.5,0.7,0.8,0.9,0.95,1.0])
#xpts = np.array([0.1,0.2,0.35,0.5,0.7,0.8,0.9,0.95])

def parse(fname):
    nY, nX = np.loadtxt(fname, max_rows=1, skiprows=3, unpack=True, dtype=int)
    data = np.loadtxt(fname, skiprows=4)
    data = np.reshape(data, ((nX, nY, -1)))
    Yran = np.array(data[0,:,0])
    Xran = np.array(data[:,0,1])
    data = np.swapaxes(data, 0, 1)
    return data, Yran, Xran, nX, nY

files = ['linear.dat','cubic.dat']

with open(files[0]) as f:
    title = f.readline().strip()
    xlabel = f.readline().strip()
    ylabel = f.readline().strip()


# set up plot and labels
fig, ax = plt.subplots()
ax.set_title(r'$\log_{10} T = 6.5$, $\log_{10} \rho = -5$, $Z=0$')
ax.set_xlabel(xlabel)
ax.set_ylabel(title)
    
for infile in files:
    kapDT, Yran, Xran, nX, nY = parse(infile)
    
    logkap_logTmid = kapDT[int(nY/2),:,2]
    logkap_logTmid_m1 = kapDT[int(nY/2)-1,:,2]
    logkap_logTmid_p1 = kapDT[int(nY/2)+1,:,2]

    ax.plot(Xran,logkap_logTmid,label='interpolated derivative')
    #ax.plot(Xran,logkap_logTmid_p1,ls='--')
    ax.set_xlim(Xran.min(), Xran.max())

    finterp = interp1d(Xran,logkap_logTmid,fill_value='extrapolate')
    ax.scatter(xpts,finterp(xpts))

# test out python interpolation over the xpts too
#ypts = finterp(xpts)
#f = interp1d(xpts,ypts,kind='cubic',fill_value='extrapolate')
#x = np.linspace(0,1,1000)
#ax.plot(x,f(x),c='tab:purple')
    
# now try to evaluate numerically from opacity files
#files = ['opacity_linear.dat','opacity_cubic.dat']
files = ['opacity_cubic.dat']
for infile in files:
    kapDT, Yran, Xran, nX, nY = parse(infile)
    logkap_logTmid = kapDT[int(nY/2),:,2]
    logkap_logTmid_m1 = kapDT[int(nY/2)-1,:,2]
    logkap_logTmid_p1 = kapDT[int(nY/2)+1,:,2]
    
    # natural logs for derivatives
    logTmid = np.log(10**kapDT[int(nY/2),0,0])
    logTmid_m1 = np.log(10**kapDT[int(nY/2)-1,0,0])
    logTmid_p1 = np.log(10**kapDT[int(nY/2)+1,0,0])
    
    print(np.log10(np.exp(logTmid_m1)),np.log10(np.exp(logTmid)),np.log10(np.exp(logTmid_p1)))
    
    nderiv_m1 = (logkap_logTmid - logkap_logTmid_m1)/(logTmid - logTmid_m1)
    nderiv_p1 = (logkap_logTmid - logkap_logTmid_p1)/(logTmid - logTmid_p1)

    #ax.plot(Xran,nderiv_m1,ls=':',c='k')
    #ax.plot(Xran,nderiv_p1,ls='--',c='k')
    ax.plot(Xran,(nderiv_p1+nderiv_m1)/2,ls='--',c='k',label='numeric derivative')

ax.legend()
    
# save figure
fig.savefig('line_kap.png', dpi=300)

