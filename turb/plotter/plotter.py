#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# set up plot and labels
fig, ax = plt.subplots(figsize=(4,3))

data = np.genfromtxt('turb_plotter.dat')
R0 = data[:,1]
# these are dimensionless diffusion coefficients Dth/kappa_C
D0 = data[:,2]
D1 = data[:,3]
D2 = data[:,4]

D1_frg = data[:,5]
D2_frg = data[:,6]

plt.plot(R0,D2,label=r'$H_B = H_{B2}$ from inlist')
plt.plot(R0,D1,label=r'$H_B = H_{B1}$ from inlist')
plt.plot(R0,D0,label=r'$H_B = 0$')

plt.scatter(R0,D2_frg,s=10)
plt.scatter(R0,D1_frg,s=7)

plt.yscale('log')
plt.ylim(bottom=3e-3)

plt.legend(frameon=True)

plt.xlabel(r'$R_0$')
plt.ylabel(r'$D_C/\kappa_C$')

# save figure
fig.savefig('turb_plotter.png', dpi=300)
