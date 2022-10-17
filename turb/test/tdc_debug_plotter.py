'''
This script can be used to plot the debug output from TDC.

TDC can write data to out.data. When it does so, it writes Q(Y).
The job of the TDC solver is to solve for Y such that Q(Y)==0, so
visualizing this function can be helpful.

Here we plot (-Y, Q(Y)), just because the non-trivial cases all have Y < 0
and it's easier to see using a log scale.
'''
import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('out.data')

plt.plot(-data[:, 0], data[:, 1])
plt.xscale('log')
plt.yscale('symlog', linthreshy=1e30)
plt.show()
