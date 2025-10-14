import mesa_reader as mr
import matplotlib.pyplot as plt

# A simple plotting script for MESA data with mesa_reader
# mesa_reader offers a lot of functionality!
# Check out the documentation at:
# https://billwolf.space/py_mesa_reader/

# Specify the directory containing MESA logs
logs = mr.MesaLogDir("LOGS")

# Read the history and profile data
history = logs.history
profile = logs.profile_data(model_number=-1)  # -1 for the last model

# Make an HR diagram from the history data
plt.figure()
plt.plot(history.log_Teff, history.log_L)
plt.xlabel(r"$\log_{10}(T_{\mathrm{eff}})\ \mathrm{[K]}$")
plt.ylabel(r"$\log_{10}(L)\ [L_\odot]$")
plt.gca().invert_xaxis()

# Plot the density profile as a function of the mass coordinate
plt.figure()
plt.plot(profile.mass, profile.Rho)
plt.xlabel(r"$M(r)\ [M_\odot]$")
plt.ylabel(r"$\rho\ [\mathrm{g/cm}^3]$")
plt.yscale("log")

# Show plots
plt.show()
