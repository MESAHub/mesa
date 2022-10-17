#!/usr/bin/env python

"""Script to get targets for MESA's `surface_effects` test case.
Running the test case should write the frequency calculation data to
`freqs.dat`.  We then use Python to independently fit the different
surface effects and write the target values to the terminal.

"""

import numpy as np
from scipy.optimize import leastsq

leastsq_kwargs = {'ftol': 1e-12, 'xtol': 1e-12}


def get_Qnl(l, mdl, Ek):
    Qnl = np.zeros_like(Ek)
    Qnl[l == 0] = 1
    Qnl[l > 0] = Ek[l > 0] / np.interp(mdl[l > 0], mdl[l == 0], Ek[l == 0])
    return Qnl


def cubic(obs, err, mdl, Ek):
    X = (np.array(mdl) ** 3 / (np.array(Ek) * np.array(err))).reshape((-1, 1))
    return np.linalg.lstsq(X, (np.array(obs) - np.array(mdl)) / np.array(err),
                           rcond=-1)[0]


def combined(obs, err, mdl, Ek):
    X = np.power(np.array(mdl).reshape((-1, 1)), [[-1, 3]]) / (
                np.array(Ek) * np.array(err)).reshape((-1, 1))
    return np.linalg.lstsq(X, (np.array(obs) - np.array(mdl)) / np.array(err),
                           rcond=-1)[0]


def power_law(l, obs, err, mdl, Ek):
    Qnl = get_Qnl(l, mdl, Ek)
    return leastsq(lambda z: (z[0] * mdl ** z[1] + Qnl * (mdl - obs)) / err,
                   [-10.0, 3.0], **leastsq_kwargs)[0]


def f_sonoi(nu, numax, a, b):
    return a * numax * (1.0 - 1.0 / (1.0 + (nu / numax) ** b))


def sonoi(numax, l, obs, err, mdl, Ek):
    Qnl = get_Qnl(l, mdl, Ek)
    return leastsq(
        lambda z: (f_sonoi(mdl, numax, z[0], z[1]) + Qnl * (mdl - obs)) / err,
        [-0.01, 8.0], **leastsq_kwargs)[0]


def print_target(name, value):
    s = '         target_{:s} = {:23.16e}'.format(name, value)
    s = s.replace('e+', 'd+').replace('e-', 'd-')
    print(s)


l, obs, err, mdl, Ek = np.loadtxt('freqs.dat', skiprows=1).T

a3 = cubic(obs, err, mdl, Ek)[0]
b1, b3 = combined(obs, err, mdl, Ek)
p0, p1 = power_law(l, obs, err, mdl, Ek)
s0, s1 = sonoi(3090.0, l, obs, err, mdl, Ek)
p0 = p0 * 3090. ** p1

print_target('a3', a3)
print_target('b1', b1)
print_target('b3', b3)
print_target('p0', p0)
print_target('p1', p1)
print_target('s0', s0)
print_target('s1', s1)

# finally, do Kjeldsen correction, which is more complicated

obs0 = obs[l == 0]
mdl0 = mdl[l == 0]

numax = 3090.0
Dnu = np.median(np.diff(obs[l == 0]))

# # carefully copy how MESA does it
# numax_sun = 3100.0
# Dnu_sun = 135.0
# nmax = numax/Dnu*(Dnu_sun/numax_sun)*22.6-1.6
# norders = int((obs0[-1]-obs0[0])/Dnu + 0.5) + 1
# n = np.ones_like(obs0)
# n[0] = np.floor(nmax - (norders-1)//2)
# n[1:] = n[0] + np.floor((obs0[1:]-obs0[0])/Dnu + 0.5)
# print(n)

# or do it your own way
n = np.floor(obs0 / Dnu) - 1.0

Dnu_obs = np.sum((obs0 - obs0.mean()) * (n - n.mean())) \
          / np.sum((n - n.mean()) ** 2)  # KBCD (8)
Dnu_mdl = np.sum((mdl0 - mdl0.mean()) * (n - n.mean())) \
          / np.sum((n - n.mean()) ** 2)  # KBCD (9)

b = 4.9
r = (b - 1) / (b * mdl0.mean() / obs0.mean() - Dnu_mdl / Dnu)  # KBCD (6)
a = (obs0.mean() - r * mdl0.mean()) * len(obs0) / np.sum((obs0 / numax) ** b)

print_target('a_div_r', a / r)
