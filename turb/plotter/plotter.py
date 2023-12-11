#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# set up plot and labels
fig, ax = plt.subplots(figsize=(4,3))

data = np.genfromtxt('turb_plotter.dat',names=True)
R0 = data['R0']

# these are dimensionless diffusion coefficients Dth/kappa_C
D0 = data['Dth_HG19_HB0']
D1 = data['Dth_HG19_HB1']
D2 = data['Dth_HG19_HB2']

D0_frg = data['Dth_FRG24_HB0']
D1_frg = data['Dth_FRG24_HB1']
D2_frg = data['Dth_FRG24_HB2']

plt.plot(R0,D2,label=r'$H_B = H_{B2}$ from inlist',ls='--',dashes=(1,1),lw=1.2,c='tab:blue')
plt.plot(R0,D1,label=r'$H_B = H_{B1}$ from inlist',ls='--',dashes=(1,1),lw=1.2,c='tab:green')
plt.plot(R0,D0,label=r'$H_B = 0$',ls='--',dashes=(4,2),lw=1.2,c='tab:gray')

plt.scatter(R0,D2_frg,s=10,color='tab:blue')
plt.scatter(R0,D1_frg,s=7,color='tab:green')
plt.scatter(R0,D0_frg,s=7,color='tab:gray')

plt.yscale('log')
plt.ylim(bottom=3e-3)

plt.legend(frameon=True)

plt.xlabel(r'$R_0$')
plt.ylabel(r'$D_C/\kappa_C$')

# save figure
fig.savefig('turb_plotter.png', dpi=300)
