import matplotlib.pyplot as plt
import numpy as np
import os
import mesa_reader as mr

# pip install git+https://github.com/wmwolf/py_mesa_reader.git

# if the directory plt_out/ does not exitst, make it
if not os.path.exists("plt_out"):
    os.makedirs("plt_out")

# enable Latex for figure labels
plt.rcParams["text.usetex"] = True

# read in the history data from MESA
history = mr.MesaData("LOGS/history.data")

# load the last profile saved (largest model number)
load_dir = mr.MesaLogDir("./LOGS")
profile = load_dir.profile_data()

# Plot HR diagram
HR_logT_min = 3.55
HR_logT_max = 3.85
HR_logL_min = 0.1
HR_logL_max = 1.0
plt.figure()
plt.plot(history.log_Teff, history.log_L)
plt.plot(history.log_Teff[-1], history.log_L[-1], "ro")
plt.xlabel(r"$\log_{10}(T_{\rm eff})$")
plt.ylabel(r"$\log_{10}(L/L_{\odot})$")
plt.title("HR Diagram")
plt.xlim(HR_logT_min, HR_logT_max)
plt.ylim(HR_logL_min, HR_logL_max)
plt.gca().set_box_aspect(1)
plt.gca().invert_xaxis()
plt.savefig("plt_out/HR_diagram.svg", bbox_inches="tight", pad_inches=0)
# plt.show()

# Plot Profile of metal mass fraction
plt.figure()
plt.plot(profile.mass, profile.z_mass_fraction_metals)
plt.xlabel(r"$\log_{10}(M(r)/M_{\odot})$")
plt.ylabel(r"metal mass fraction")
plt.title("Metal mass fraction profile")
plt.xlim(0, profile.mass[0])
plt.ylim(0.038, 0.042)
plt.savefig("plt_out/z_mass_fraction_metals.svg", bbox_inches="tight", pad_inches=0)
# plt.show()
