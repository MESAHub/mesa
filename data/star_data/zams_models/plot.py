#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize=(4, 3))
Lsun = 3.839e33  # erg/s

masses, zones = np.genfromtxt(
    "zams_z2m2_y28.data", skip_header=14, max_rows=35, unpack=True
)

Teffs = []
Ls = []
h1_center = []
for i, mass in enumerate(masses):
    nz = zones[i]
    lines_to_skip = 50 + 5 * (i + 1) + int(np.sum(zones[:i]))
    lnT, L, h1 = np.genfromtxt(
        "zams_z2m2_y28.data",
        skip_header=lines_to_skip,
        max_rows=nz,
        usecols=(2, 4, 6),
        unpack=True,
    )

    Teffs.append(np.exp(lnT[0]))
    Ls.append(L[0] / Lsun)
    h1_center.append(h1[-1])

plt.scatter(masses, h1_center)
plt.xscale("log")

plt.xlabel(r"Mass [$M_\odot$]")
plt.ylabel(r"Center $X$")

plt.ylim(0.6965, 0.6975)

plt.tight_layout()

plt.savefig("CenterXfrac.png", dpi=200)

plt.clf()
plt.scatter(Teffs, Ls)
plt.scatter([5778], [1], marker="*")
plt.xscale("log")
plt.yscale("log")
plt.gca().invert_xaxis()

plt.xlabel(r"$T_{\rm eff}$ [K]")
plt.ylabel(r"$L$ [$L_\odot$]")

plt.tight_layout()

plt.savefig("HR.png", dpi=200)
