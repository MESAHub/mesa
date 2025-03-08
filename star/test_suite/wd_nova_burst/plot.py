import matplotlib.pyplot as plt
import numpy as np
import os
import mesa_reader as mr

# if the directory plt_out/ does not exits, make it
if not os.path.exists("plt_out"):
    os.makedirs("plt_out")

# "MESA" styles for plotting
# note that this requires having a LaTeX installation locally.
# If you don't have that, you can comment out the "text.usetex" line.
plt.rcParams["figure.figsize"] = (3.38, 2.535)
plt.rcParams["lines.markersize"] = 4
plt.rcParams["lines.linewidth"] = 1.5
plt.rcParams["text.usetex"] = True
plt.rcParams["font.size"] = 10
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "Computer Modern Roman"
plt.rcParams["axes.titlesize"] = "medium"
plt.rcParams["axes.labelsize"] = "medium"
plt.rcParams["legend.fontsize"] = 8
plt.rcParams["legend.frameon"] = False
plt.rcParams["figure.dpi"] = 300

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.top"] = True
plt.rcParams["ytick.right"] = True
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True

plt.rcParams["savefig.bbox"] = "tight"
plt.rcParams["savefig.pad_inches"] = 0.1
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["savefig.format"] = "svg"

plt.rcParams["axes.formatter.use_mathtext"] = True

# read in the history data from MESA
logs = mr.MesaLogDir("LOGS")
h = logs.history
p = logs.profile_data()

# Create a figure
fig = plt.figure(figsize=(10, 6))

# Create a GridSpec with 6 rows and 3 columns
gs_outer = fig.add_gridspec(1, 2, width_ratios=[1, 2])
gs_left = gs_outer[0].subgridspec(2, 1, hspace=0.3)
gs_right = gs_outer[1].subgridspec(3, 2, hspace=0, wspace=0.75)

# First column: Two plots of roughly equal height
ax1 = fig.add_subplot(gs_left[0])
ax2 = fig.add_subplot(gs_left[1])

# Second column: Three stacked axes sharing the same x-axis
ax3 = fig.add_subplot(gs_right[0, 0])
ax4 = fig.add_subplot(gs_right[1, 0])
ax5 = fig.add_subplot(gs_right[2, 0])

# Third column: Three stacked axes sharing a different x-axis
ax6 = fig.add_subplot(gs_right[0, 1])
ax7 = fig.add_subplot(gs_right[1, 1])
ax8 = fig.add_subplot(gs_right[2, 1])

# HR diagram
ax1.loglog(h.Teff, h.L)
ax1.set_xlabel(r"Effective Temperature [K]")
ax1.set_ylabel(r"Luminosity [$L_{\odot}$]")
ax1.invert_xaxis()
# add marker to last point
ax1.plot(h.Teff[-1], h.L[-1], "o", color="C3")

# T-Rho Profile
# load data from MESA_DIR and save to arrays in T-Rho space
plot_info_dir = os.path.join(os.environ["MESA_DIR"], "data/star_data/plot_info")

# degeneracy data
psi4_file = os.path.join(plot_info_dir, "psi4.data")
psi4_data = np.genfromtxt(psi4_file)
psi4_xs = 10 ** np.array(psi4_data.T[0])
psi4_ys = 10 ** np.array(psi4_data.T[1])

# hydrogen burn line
hydrogen_burn_file = os.path.join(plot_info_dir, "hydrogen_burn.data")
hydrogen_burn_data = np.genfromtxt(hydrogen_burn_file)
hydrogen_burn_xs = 10 ** np.array(hydrogen_burn_data.T[0])
hydrogen_burn_ys = 10 ** np.array(hydrogen_burn_data.T[1])

# helium burn line
helium_burn_file = os.path.join(plot_info_dir, "helium_burn.data")
helium_burn_data = np.genfromtxt(helium_burn_file)
helium_burn_xs = 10 ** np.array(helium_burn_data.T[0])
helium_burn_ys = 10 ** np.array(helium_burn_data.T[1])

# plot raw T-Rho data; sets limits for us
ax2.loglog(p.Rho, p.T)

# plot psi4 and hydrogen burn lines
# preserve limits... shouldn't have to do it like this, but
# autoscale commands always seem to fail
x_left, x_right = ax2.get_xlim()
y_bottom, y_top = ax2.get_ylim()
ax2.plot(psi4_xs, psi4_ys, color="lightgrey", ls=":", zorder=-5)
ax2.plot(hydrogen_burn_xs, hydrogen_burn_ys, ls="--", color="lightgrey", zorder=-5)
ax2.plot(helium_burn_xs, helium_burn_ys, ls="--", color="lightgrey", zorder=-5)
ax2.set_xlim(x_left, x_right)
ax2.set_ylim(y_bottom, y_top)
ax2.set_xlabel(r"Density [g/cm$^3$]")
ax2.set_ylabel(r"Temperature [K]")

# plot strong burning zones
high_burning = np.where(p.eps_nuc > 1e7, True, False)
mid_burning = np.where((p.eps_nuc > 1e3) & (p.eps_nuc < 1e7), True, False)
low_burning = np.where((p.eps_nuc > 1) & (p.eps_nuc < 1e3), True, False)
# temporarily set cap style to round
save_capstyle = plt.rcParams["lines.solid_capstyle"]
plt.rcParams["lines.solid_capstyle"] = "round"
ax2.plot(
    p.Rho[high_burning],
    p.T[high_burning],
    marker="o",
    ls="",
    ms=6,
    color="C3",
    zorder=-1,
)
ax2.plot(
    [],
    [],
    lw=6,
    color="C3",
    label=r"$\epsilon_{\mathrm{nuc}} > 10^7\,\mathrm{erg/g/s}$",
)
ax2.plot(
    p.Rho[mid_burning],
    p.T[mid_burning],
    marker="o",
    ls="",
    ms=4.5,
    color="C1",
    zorder=-2,
)
ax2.plot(
    [],
    [],
    color="C1",
    lw=4.5,
    label=r"$\epsilon_{\mathrm{nuc}} > 10^3\,\mathrm{erg/g/s}$",
)
ax2.plot(
    p.Rho[low_burning],
    p.T[low_burning],
    marker="o",
    ls="",
    ms=3,
    color="goldenrod",
    zorder=-3,
)
ax2.plot(
    [],
    [],
    color="goldenrod",
    lw=3,
    label=r"$\epsilon_{\mathrm{nuc}} > 1\,\mathrm{erg/g/s}$",
)
ax2.legend(loc="lower right")
# restore capstyle
plt.rcParams["lines.solid_capstyle"] = save_capstyle

# time series for last 5 years
window = 5

# top panel: photospheric and H-burning luminosities
mask = h.star_age > (max(h.star_age) - window)
ax3.set_title("Last 5 Years")
ax3.semilogy(h.star_age[mask], h.L[mask], label=r"$L$")
ax3.semilogy(h.star_age[mask], h.LH[mask], ls="--", label=r"$L_{\mathrm{H}}$")
ax3.set_ylabel(r"Luminosity [$L_{\odot}$]")
ax3.legend(loc="upper left")
ax3.set_xticklabels([])

# temperature and radius
ax4.semilogy(h.star_age[mask], h.Teff[mask], label=r"$T_{\mathrm{eff}}$")
ax4.semilogy([], [], label=r"$R$", ls="--", color="C1")
ax4.legend(loc="upper left")
ax4b = ax4.twinx()
ax4b.semilogy(h.star_age[mask], h.R[mask], ls="--", color="C1")
ax4.set_ylabel(r"Eff. Temp. [K]", color="C0")
ax4b.set_ylabel(r"Radius [$R_{\odot}$]", color="C1")
ax4.set_xticklabels([])

# mass loss and envelope mass
ax5.semilogy(
    h.star_age[mask],
    h.star_mass[mask] - h.he_core_mass[mask],
    label=r"$\Delta M_{\mathrm{H}}$",
)
ax5.semilogy([], [], ls="--", label=r"$|\dot{M}|$")
ax5.legend(loc="center left")
ax5.set_ylabel(r"Envelope Mass [$M_{\odot}$]", color="C0")
ax5b = ax5.twinx()
ax5b.semilogy(h.star_age[mask], h.abs_mdot[mask], ls="--", color="C1")
ax5b.set_ylabel(r"$\left\vert\dot{M}\right\vert$ [$M_{\odot}$/yr]", color="C1")

ax5.set_xlabel(r"Star Age [yr]")

# profile at the end of the simulation
# first up, common isotope mass fractions
xm = p.star_mass - p.mass
for iso in ["h1", "he4", "c12", "n14", "o16"]:
    element = "".join([i for i in iso if not i.isdigit()]).capitalize()
    mass_number = "".join([i for i in iso if i.isdigit()])
    tex_iso = r"${}^{" + mass_number + r"}$" + element
    ax6.loglog(xm, p.data(iso), label=tex_iso)

ax6.set_title("Final Profile")
ax6.set_ylim(4e-5, 1.5)
ax6.legend(loc="lower left")
ax6.set_ylabel(r"Mass Fraction")
ax6.set_xticklabels([])

# nuclear energy generation rates
ax7.loglog(xm, p.eps_nuc, label=r"$\epsilon_{\mathrm{nuc}}$")
ax7.loglog(xm, p.pp, ls="--", label=r"$\epsilon_{\mathrm{pp}}$")
ax7.loglog(xm, p.cno, ls=":", label=r"$\epsilon_{\mathrm{CNO}}$")
ax7.legend(loc="upper right")
ax7.set_ylabel(r"Specific Power [erg/g/s]")
ax7.set_xticklabels([])
ax7.set_ylim(1e-6, 1e10)

# basic thermodynamic quantities
ax8.loglog(xm, p.Rho, label="$\\rho$")
ax8.loglog([], [], ls="--", label="$T$")
ax8.set_ylim(2e-3, 2e6)
ax8.set_ylabel(r"Density [g/cm$^3$]", color="C0")
ax8b = ax8.twinx()
ax8b.loglog(xm, p.T, ls="--", color="C1")
ax8b.set_ylim(2e6, 1.25e8)
ax8b.set_ylabel(r"Temperature [K]", color="C1")
ax8.legend(loc="upper right")

for ax in [ax6, ax7, ax8]:
    ax.set_xlim(8e-4, 1e-10)
ax8.set_xlabel(r"Exterior Mass Coordinate [M$_{\odot}$]")


fig.savefig("plt_out/wd_nova_burst_grid.svg")
