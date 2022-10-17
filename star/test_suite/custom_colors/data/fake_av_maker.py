import numpy as np

# This makes some fake data for demonstrating MESA's ability to interpolate
# over arbitrary data in the colors module


tmin = 3.0
tmax = 6.0
nt = 21

gmin = -5.0
gmax = 5.0
ng = 21

amin = 0.0
amax = 10.0
na = 21

with open("fake_av_v.txt", 'w') as f:
    print("#T logg Av av_v", file=f)
    for i in np.linspace(tmin, tmax, nt):
        for j in np.linspace(gmin, gmax, ng):
            for k in np.linspace(amin, amax, na):
                rr = i / (tmax - tmin) + j / (gmax - gmin) + k / (amax - amin)
                print("{:10.4f}".format(10 ** i), "{:10.4f}".format(j),
                      "{:10.4f}".format(k), "{:10.4f}".format(rr), file=f)
