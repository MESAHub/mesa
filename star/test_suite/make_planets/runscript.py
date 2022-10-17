#!/usr/bin/env python

import math
import os
import shutil

import mysubprograms as my

msun = 1.9892e33
rsun = 6.9598e10
mjup = 1.8986e30
rjup = 6.9911e9
mearth = 5.97e27
sigma = 5.67e-5
au = 1.496e13

# fiducial parameters
mp = 1.0  # planet mass in mjup
rp = 2.0  # inital planet radius in rjup
mcore = 10.0  # planet core mass in mearth
rhocore = 10.0  # core density in g/cc
mp_wo_core = mp - mcore * mearth / mjup  # mass of planet before core put in
z = 0.02  # metallicity of both planet and star
y = 0.24  # helium fraction of both planet and (initial) star
maxage = 1.e10  # ending age

# parameters having to do with irradiation
ms = 1.0  # star mass in msun
rs = 1.0  # star radius in rsun
Teff_star = 5800  # stellar Teff
orb_sep = 0.05  # orbital sepration in AU (to set day-side flux)
irrad_col = 300.0  # column depth for depositing stellar radiation as heat
flux_dayside = sigma * Teff_star ** 4 * (
            rs * rsun / orb_sep / au) ** 2  # flux hitting planet's dayside
Teq = (flux_dayside / 4.0 / sigma) ** 0.25

# flags to skip steps
do_create_planet = True
do_put_in_core = True
do_evolve_planet = True

f = open('logfile', 'w')

lgmmin = -1.0
lgmmax = math.log10(80.0)
imax = 100
for i in range(1, imax + 1, 1):
    
    mp = math.pow(10.0, lgmmin + (lgmmax - lgmmin) * i / imax)
    mp_wo_core = mp - mcore * mearth / mjup
    
    my.print_parameters(mp, rp, mcore, rhocore, mp_wo_core, irrad_col,
                        flux_dayside, Teq, y, z, maxage)
    
    createmodel = "planet_create_" + str(mp_wo_core)[0:6] + "_MJ_" + str(
        rp) + "_RJ.mod"
    coremodel = "planet_core_" + str(mp)[0:6] + "_MJ_" + str(mcore)[
                                                         0:6] + "_ME_" + str(
        rp) + "_RJ.mod"
    evolvemodel = "planet_evolve_" + str(mp)[0:6] + "_MJ_" + str(mcore)[
                                                             0:6] + "_ME_" +\
                  str(
        rp) + "_RJ.mod"
    
    if do_create_planet:
        inlist1 = "inlist_create_" + str(mp)[0:6] + "_MJ_" + str(mcore)[
                                                             0:6] + "_ME_" +\
                  str(
            rp) + "_RJ"
        run_time = my.create_planet(mp_wo_core, rp, y, z, inlist1, createmodel)
    
    success = True
    if not os.path.exists(createmodel):
        success = False
    k = open('LOGS/history.data', 'r')
    for line in k.readlines():
        pass
    last_temp = line
    k.close()
    last = last_temp.split()
    print
    "final model number in create=", last[0]
    print
    "last[0]==1000", last[0] == "1000"
    if last[0] == "1000":
        success = False
    outstring = '%6.3f\t%6.3f\t%6.3f\t%s\n' % (mp, rp, run_time, success)
    f.write(outstring)
    if not success:
        continue
    
    if do_put_in_core:
        if mcore > 0.0:
            inlist2 = "inlist_core_" + str(mp)[0:6] + "_MJ_" + str(mcore)[
                                                               0:6] + "_ME_"\
                      + str(
                rp) + "_RJ"
            run_time = my.put_core_in_planet(mcore, rhocore, inlist2,
                                             createmodel, coremodel)
        else:
            shutil.copyfile(createmodel, coremodel)
    
    if not os.path.exists(coremodel):
        continue
    
    if do_evolve_planet:
        inlist3 = "inlist_evolve_" + str(mp)[0:6] + "_MJ_" + str(mcore)[
                                                             0:6] + "_ME_" + str(
            rp) + "_RJ"
        run_time = my.evolve_planet(irrad_col, flux_dayside, maxage, inlist3,
                                    coremodel, evolvemodel)
    
    if not os.path.exists(evolvemodel):
        continue

f.close()
