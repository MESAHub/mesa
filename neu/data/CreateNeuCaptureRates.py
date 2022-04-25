# ***********************************************************************
#
#   Copyright (C) 2022  The MESA Team
#
#   this file is part of mesa.
#
#   mesa is free software; you can redistribute it and/or modify
#   it under the terms of the gnu general library public license as published
#   by the free software foundation; either version 2 of the license, or
#   (at your option) any later version.
#
#   mesa is distributed in the hope that it will be useful, 
#   but without any warranty; without even the implied warranty of
#   merchantability or fitness for a particular purpose.  see the
#   gnu library general public license for more details.
#
#   you should have received a copy of the gnu library general public license
#   along with this software; if not, write to the free software
#   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
#
# ***********************************************************************
 
# This computes the neutrino capture rates using SKYNET
# This was last run with SKYNET version: 49ff6c91d58d7bcc97143b649329ed186b0958ae
# With the following patches applied:
# 0001-Make-it-easier-to-call-the-neutrino-capture-rates-fr.patch

from SkyNet import *
import numpy as np
import itertools
import os
import subprocess


# Setup Skynet
nuclib = NuclideLibrary.CreateFromWebnucleoXML(SkyNetRoot
    + "/data/webnucleo_nuc_v2.0.xml")

reaclib = REACLIB(SkyNetRoot + "/data/reaclib")

opts = NetworkOptions()
opts.ConvergenceCriterion = NetworkConvergenceCriterion.Mass
opts.MassDeviationThreshold = 1.0E-10
opts.IsSelfHeating = True
opts.CalculateStrongInverseRates = True
opts.DisableStdoutOutput = True

out = NetworkOutput.CreateNew("precompute_reaction_libs_with_neutrinos",
    nuclib, False)

nl = NeutrinoReactionLibrary(
    SkyNetRoot + "/data/neutrino_reactions.dat",
    "Neutrino interactions", nuclib, opts, False, True)
print("Loaded neutrino library")


def compute_skynet(T9,mu):
    tstate = ThermodynamicState()
    tstate.SetT9(T9)
    tstate.SetEtaElectron(mu)

    # Creates Fermi-dirac function for the neutrino potential with mu=0
    nd = NeutrinoDistributionFermiDirac_Create(
        T9,
        0.0
    )

    tstate.SetNeutrinoDistribution(nd)

    # unused parameters but need setting to null
    partFunc = []
    expArgCap = 0.0
    pZs = std_vector_int([0])
    pScreen = std_vector_double([0.0])

    nl.CalculateRates(
        tstate,
        partFunc,
        expArgCap,
        pZs,
        pScreen,
    )

    return *nl.Rates(),*nl.InverseRates(),*nl.HeatingRates(), *nl.HeatingInverseRates()


minlogT = 6.0
maxlogT = 11.0
numT = 101 

minMu = -10.0
maxMu = 5.0
numMu = 101

logTs = np.linspace(minlogT,maxlogT,numT)
mus = np.linspace(minMu,maxMu,numMu)

filename = './neutrino_captures.txt'

with open(filename,'w') as f:
    print(f"# {numT*numMu} {numT} {numMu}",file=f)
    print("# logT mu prot_e_cap neut_pos_cap neut_neu_cap prot_aneu_cap Q_prot_e_cap Q_neut_pos_cap Q_neut_neu_cap Q_prot_aneu_cap",file=f)
    for idx,(logT,mu) in enumerate(itertools.product(logTs,mus)):
        t = (10**logT)/10**9
        x = compute_skynet(t,mu)
        print(f"{logT:.16e} {mu:.16e}",*['{:.16e}'.format(j) for j in x],file=f)

subprocess.run(['xz','-zk',filename])
