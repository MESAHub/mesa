#!/usr/bin/env python3
"""Plot nonlinear growth and outer-boundary diagnostics for this test case."""

from __future__ import annotations

import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


CASE = Path(__file__).resolve().parents[1]
LOGS = CASE / "LOGS"
OUT = Path(__file__).resolve().parent

PHOTO_COLUMNS = [
    "photosphere_black_body_T",
    "photosphere_cell_T",
    "photosphere_cell_log_T",
    "photosphere_cell_density",
    "photosphere_cell_log_density",
    "photosphere_cell_opacity",
    "photosphere_cell_log_opacity",
    "photosphere_L",
    "photosphere_log_L",
    "photosphere_r",
    "photosphere_log_r",
    "photosphere_m",
    "photosphere_v_km_s",
    "photosphere_cell_k",
    "photosphere_column_density",
    "photosphere_csound",
    "photosphere_log_column_density",
    "photosphere_opacity",
    "photosphere_v_div_cs",
    "photosphere_xm",
    "photosphere_cell_free_e",
    "photosphere_cell_log_free_e",
    "photosphere_logg",
    "photosphere_T",
]

OLD_EXTRAS = [
    "num_periods",
    "period",
    "growth",
    "max_v_div_cs",
    "max_v_div_vesc",
    "delta_R",
    "delta_Teff",
    "delta_logTeff",
    "delta_logL",
    "delta_Mag",
    "KE_growth",
    "KE_growth_avg",
]

NEW_EXTRAS = [
    "num_periods",
    "period",
    "growth",
    "max_v_div_cs",
    "max_v_div_vesc",
    "delta_R",
    "delta_rphot",
    "delta_V",
    "delta_Vphot",
    "delta_Teff",
    "delta_logTeff",
    "delta_logL",
    "delta_Mag",
    "KE_growth",
    "KE_growth_avg",
]


def file_line(path: Path, number: int) -> list[str]:
    with path.open() as stream:
        for i, line in enumerate(stream, start=1):
            if i == number:
                return line.split()
    raise RuntimeError(f"missing line {number} in {path}")


def current_history_names(path: Path) -> list[str]:
    header = file_line(path, 6)
    last = subprocess.run(
        ["tail", "-n", "1", str(path)],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.split()
    ndata = len(last)
    if len(header) == ndata:
        return header

    # Native history columns are read without recompiling.  This case therefore
    # can append the newly enabled photosphere columns to an existing history
    # file whose original header does not contain them.
    if "num_periods" not in header:
        raise RuntimeError("cannot identify the extra-history column block")
    standard = header[: header.index("num_periods")]
    insert_at = standard.index("log_Teff") + 1
    standard = standard[:insert_at] + PHOTO_COLUMNS + standard[insert_at:]
    nextras = ndata - len(standard)
    if nextras == len(OLD_EXTRAS):
        return standard + OLD_EXTRAS
    if nextras == len(NEW_EXTRAS):
        return standard + NEW_EXTRAS
    raise RuntimeError(
        f"history has {ndata} values, but its stale header has {len(header)} names"
    )


def read_history_tail(path: Path, nlines: int = 30000) -> dict[str, np.ndarray]:
    names = current_history_names(path)
    wanted = [
        "model_number",
        "day",
        "log_Teff",
        "log_L",
        "radius",
        "v_surf_km_s",
        "v_div_csound_surf",
        "photosphere_log_L",
        "photosphere_r",
        "photosphere_v_km_s",
        "photosphere_v_div_cs",
        "num_periods",
        "period",
        "growth",
        "max_v_div_cs",
        "delta_R",
        "delta_logL",
        "KE_growth",
        "KE_growth_avg",
    ]
    wanted = [name for name in wanted if name in names]
    indices = [names.index(name) for name in wanted]
    values = {name: [] for name in wanted}

    process = subprocess.Popen(
        ["tail", "-n", str(nlines), str(path)],
        stdout=subprocess.PIPE,
        text=True,
    )
    assert process.stdout is not None
    for line in process.stdout:
        fields = line.split()
        if len(fields) != len(names):
            continue
        try:
            selected = [float(fields[i]) for i in indices]
        except ValueError:
            continue
        for name, value in zip(wanted, selected):
            values[name].append(value)
    if process.wait() != 0:
        raise RuntimeError("tail failed while reading history.data")

    data = {name: np.asarray(column) for name, column in values.items()}
    if not len(data["model_number"]):
        raise RuntimeError("no history rows were read")

    # Keep the newest branch if history.data contains an appended restart.
    rollback = np.flatnonzero(np.diff(data["model_number"]) < 0)
    if rollback.size:
        start = rollback[-1] + 1
        data = {name: column[start:] for name, column in data.items()}
    return data


def cycle_end_indices(history: dict[str, np.ndarray]) -> np.ndarray:
    cycles = history["num_periods"].astype(int)
    # The first row in a tail or restarted branch need not be a cycle boundary.
    changed = np.r_[False, cycles[1:] != cycles[:-1]]
    indices = np.flatnonzero(changed & (cycles > 1))
    return indices


def plot_cycle_history(history: dict[str, np.ndarray]) -> tuple[int, int, int]:
    ends = cycle_end_indices(history)
    if ends.size < 3:
        raise RuntimeError("not enough complete cycles in the history tail")
    ends = ends[-45:]
    cycle = history["num_periods"][ends]
    growth = history["growth"][ends]
    ke_growth = history["KE_growth"][ends]
    radius_max = history["radius"][ends]
    delta_r = history["delta_R"][ends]
    radius_min = radius_max - delta_r
    radius_mean = 0.5 * (radius_max + radius_min)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    ax = axes[0, 0]
    ax.plot(cycle, growth, "o-", label=r"$\Delta R$ growth")
    ax.plot(cycle, ke_growth, "o-", label="KE growth")
    ax.plot(cycle, (1 + growth) ** 2 - 1, "--", label=r"$(1+g_R)^2-1$")
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_ylabel("fractional growth per cycle")
    ax.legend(frameon=False)

    ax = axes[0, 1]
    ax.plot(cycle, radius_max, label=r"$R_{\max}$")
    ax.plot(cycle, radius_mean, label=r"$(R_{\max}+R_{\min})/2$")
    ax.plot(cycle, radius_min, label=r"$R_{\min}$")
    ax.set_ylabel(r"surface radius ($R_\odot$)")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    ax.plot(cycle, history["max_v_div_cs"][ends], "o-", label="max surface Mach")
    ax.plot(cycle, delta_r / radius_mean, "o-", label=r"$\Delta R/\bar R$")
    ax.axhline(1, color="0.4", lw=0.8, ls="--")
    ax.set_xlabel("completed cycle")
    ax.set_ylabel("nonlinearity measure")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    ax.plot(cycle, history["delta_logL"][ends], "o-", label=r"$\Delta\log L$")
    ax2 = ax.twinx()
    ax2.plot(cycle, history["period"][ends], color="tab:red", label="period")
    ax.set_xlabel("completed cycle")
    ax.set_ylabel(r"$\Delta\log L$")
    ax2.set_ylabel("period (day)", color="tab:red")

    fig.suptitle("Nonlinear cycle diagnostics")
    fig.savefig(OUT / "history_cycle_growth.png", dpi=180)
    plt.close(fig)

    # A completed cycle runs between two successive positive-to-negative
    # surface-velocity crossings.  The last change index is its right edge.
    all_ends = cycle_end_indices(history)
    left = int(all_ends[-2])
    right = int(all_ends[-1])
    return left, right, int(history["num_periods"][right])


def plot_recent_history(history: dict[str, np.ndarray], right: int) -> None:
    ends = cycle_end_indices(history)
    start = int(ends[-4]) if ends.size >= 4 else 0
    stop = min(right + 1, len(history["day"]))
    day = history["day"][start:stop]
    day = day - day[0]

    fig, axes = plt.subplots(4, 1, figsize=(11, 9), sharex=True, constrained_layout=True)
    axes[0].plot(day, history["radius"][start:stop], label="surface")
    if "photosphere_r" in history:
        axes[0].plot(day, history["photosphere_r"][start:stop], label="photosphere")
    axes[0].set_ylabel(r"radius ($R_\odot$)")
    axes[0].legend(frameon=False)

    axes[1].plot(day, history["v_surf_km_s"][start:stop], label="surface")
    if "photosphere_v_km_s" in history:
        axes[1].plot(day, history["photosphere_v_km_s"][start:stop], label="photosphere")
    axes[1].axhline(0, color="0.4", lw=0.8)
    axes[1].set_ylabel("velocity (km/s)")

    axes[2].plot(day, history["log_L"][start:stop], label="surface")
    if "photosphere_log_L" in history:
        axes[2].plot(day, history["photosphere_log_L"][start:stop], label="photosphere")
    axes[2].set_ylabel(r"$\log L/L_\odot$")

    axes[3].plot(day, history["log_Teff"][start:stop])
    axes[3].set_ylabel(r"$\log T_{\rm eff}$")
    axes[3].set_xlabel("time since start of displayed interval (day)")
    fig.suptitle("Last three completed cycles")
    fig.savefig(OUT / "history_recent_cycles.png", dpi=180)
    plt.close(fig)


def read_profile(path: Path) -> tuple[dict[str, float], dict[str, np.ndarray]]:
    global_names = file_line(path, 2)
    global_values = []
    for value in file_line(path, 3):
        try:
            global_values.append(float(value))
        except ValueError:
            global_values.append(np.nan)
    names = file_line(path, 6)
    values = np.loadtxt(path, skiprows=6)
    if values.ndim == 1:
        values = values[np.newaxis, :]
    global_data = dict(zip(global_names, global_values))
    data = {name: values[:, i] for i, name in enumerate(names)}
    return global_data, data


def profile_index() -> tuple[np.ndarray, np.ndarray]:
    values = np.loadtxt(LOGS / "profiles.index", skiprows=1)
    if values.ndim == 1:
        values = values[np.newaxis, :]
    return values[:, 0].astype(int), values[:, 2].astype(int)


def nearest_profile(model: int, models: np.ndarray, numbers: np.ndarray) -> Path:
    number = numbers[np.argmin(np.abs(models - model))]
    return LOGS / f"profile{number}.data"


def plot_phase_profiles(
    history: dict[str, np.ndarray], left: int, right: int
) -> tuple[float, float, float, float, float, float, float, float, float]:
    models, numbers = profile_index()
    model_lo = int(history["model_number"][left])
    model_hi = int(history["model_number"][right])
    use = (models >= model_lo) & (models <= model_hi)
    cycle_models = models[use]
    cycle_numbers = numbers[use]
    if cycle_models.size < 4:
        raise RuntimeError("not enough profiles in the last completed cycle")

    sl = slice(left, right + 1)
    phase_models = history["model_number"][sl].astype(int)
    radius = history["radius"][sl]
    velocity = history["v_surf_km_s"][sl]
    targets = [
        ("maximum radius", phase_models[np.argmax(radius)]),
        ("maximum contraction", phase_models[np.argmin(velocity)]),
        ("minimum radius", phase_models[np.argmin(radius)]),
        ("maximum expansion", phase_models[np.argmax(velocity)]),
    ]

    selected = []
    for label, model in targets:
        path = nearest_profile(int(model), cycle_models, cycle_numbers)
        selected.append((label, path, *read_profile(path)))

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.08, 0.9, len(selected)))
    for color, (label, path, global_data, data) in zip(colors, selected):
        logt = data["logT"]
        outer = (logt >= 3.45) & (logt <= 5.3)
        pressure = 10 ** data["logP"]
        axes[0, 0].plot(logt[outer], data["v_div_csound"][outer], color=color, label=label)
        axes[0, 1].plot(logt[outer], data["lum_conv_div_L"][outer], color=color)
        axes[1, 0].plot(logt[outer], data["log_opacity"][outer], color=color)
        uq = np.maximum(np.abs(data["Uq"]) / pressure, 1e-20)
        pvsc = np.maximum(np.abs(data["Pvsc"]) / pressure, 1e-20)
        axes[1, 1].plot(logt[outer], uq[outer], color=color)
        axes[1, 1].plot(logt[outer], pvsc[outer], color=color, ls="--")

    axes[0, 0].axhline(0, color="0.4", lw=0.8)
    axes[0, 0].set_ylabel(r"$v/c_s$")
    axes[0, 0].legend(frameon=False, fontsize=9)
    axes[0, 1].set_ylabel(r"$L_{\rm conv}/L$")
    axes[1, 0].set_ylabel(r"$\log\kappa$")
    axes[1, 0].set_xlabel(r"$\log T$")
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_ylabel(r"$|U_q|/P$ (solid), $P_{\rm vsc}/P$ (dashed)")
    axes[1, 1].set_xlabel(r"$\log T$")
    fig.suptitle("Outer-envelope structure through the last completed cycle")
    fig.savefig(OUT / "profile_last_cycle.png", dpi=180)
    plt.close(fig)

    # Estimate an upper-scale atmospheric work using both a hydrostatic
    # P_rad + g*tau/kappa estimate and the much larger surface-cell pressure.
    # The former is only a proxy because the exact atmosphere integration is
    # not written to the profile.
    g_cgs = 6.67430e-8
    msun = 1.988409870698051e33
    rsun = 6.957e10
    lsun = 3.828e33
    crad = 7.5657e-15
    volumes = []
    p_proxy = []
    p_cell = []
    kinetic = []
    luminosity = []
    hydrogen = []
    times = []
    max_uq_ratio = 0.0
    max_pvsc_ratio = 0.0
    max_pvsc_logt = np.nan
    for model, number in zip(cycle_models, cycle_numbers):
        global_data, data = read_profile(LOGS / f"profile{number}.data")
        r = data["radius"][0] * rsun
        mass = data["mass"][0] * msun
        temperature = 10 ** data["logT"][0]
        opacity = 10 ** data["log_opacity"][0]
        tau = data["tau"][0]
        gravity = g_cgs * mass / r**2
        prad = crad * temperature**4 / 3
        volumes.append(4 * np.pi * r**3 / 3)
        p_proxy.append(prad + gravity * tau / opacity)
        p_cell.append(10 ** data["logP"][0])
        dm = data["dq"] * global_data["star_mass"] * msun
        kinetic.append(np.sum(0.5 * dm * data["velocity"] ** 2))
        luminosity.append(data["luminosity"][0] * lsun)
        hydrogen.append(data["x_mass_fraction_H"][0])
        times.append(global_data["time_seconds"])
        pressure = 10 ** data["logP"]
        max_uq_ratio = max(max_uq_ratio, np.max(np.abs(data["Uq"]) / pressure))
        pvsc_ratio = np.abs(data["Pvsc"]) / pressure
        if np.max(pvsc_ratio) > max_pvsc_ratio:
            max_pvsc_ratio = np.max(pvsc_ratio)
            max_pvsc_logt = data["logT"][np.argmax(pvsc_ratio)]

    volumes = np.asarray(volumes)
    p_proxy = np.asarray(p_proxy)
    p_cell = np.asarray(p_cell)
    kinetic = np.asarray(kinetic)
    luminosity = np.asarray(luminosity)
    times = np.asarray(times)
    dvolume = np.diff(volumes)
    work_proxy = -np.sum(0.5 * (p_proxy[1:] + p_proxy[:-1]) * dvolume)
    work_cell = -np.sum(0.5 * (p_cell[1:] + p_cell[:-1]) * dvolume)
    work_abs = np.sum(0.5 * (p_cell[1:] + p_cell[:-1]) * np.abs(dvolume))
    radiant_energy = np.trapz(luminosity, times)
    hydrogen_span = np.ptp(hydrogen)
    return (
        work_proxy,
        work_cell,
        work_abs,
        np.max(kinetic),
        radiant_energy,
        hydrogen_span,
        max_uq_ratio,
        max_pvsc_ratio,
        max_pvsc_logt,
    )


def plot_lna_work() -> tuple[float, float, float, float]:
    path = CASE / "LNA" / "star_LNA_work_1.data"
    work = np.loadtxt(path, comments="#")
    logt = work[:, 1]
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True, constrained_layout=True)
    labels = [
        (3, "pressure"),
        (5, "eddy viscosity"),
        (6, "radiative luminosity"),
        (7, "convective luminosity"),
        (9, "total"),
    ]
    for column, label in labels:
        axes[0].plot(logt, work[:, column], label=label)
    axes[0].set_yscale("symlog", linthresh=1e-4)
    axes[0].set_ylabel("local normalized work")
    axes[0].legend(frameon=False, ncol=2)

    cumulative = [
        (10, "pressure"),
        (12, "eddy viscosity"),
        (13, "radiative luminosity"),
        (14, "convective luminosity"),
        (16, "total"),
    ]
    for column, label in cumulative:
        axes[1].plot(logt, work[:, column], label=label)
    axes[1].axhline(0, color="0.4", lw=0.8)
    axes[1].set_xlabel(r"$\log T$")
    axes[1].set_ylabel("cumulative normalized work")
    axes[1].legend(frameon=False, ncol=2)
    fig.suptitle("Linear fundamental-mode work audit")
    fig.savefig(OUT / "lna_fundamental_work.png", dpi=180)
    plt.close(fig)

    period_data = np.loadtxt(CASE / "LNA" / "star_LNA_period_growth.data", comments="#")
    mode = period_data[0]
    surface_cell_work = work[0, 9]
    return mode[6], mode[13], mode[11], surface_cell_work


def main() -> None:
    history = read_history_tail(LOGS / "history.data")
    left, right, cycle = plot_cycle_history(history)
    plot_recent_history(history, right)
    (
        work_proxy,
        work_cell,
        work_abs,
        ke_max,
        e_rad,
        x_span,
        max_uq_ratio,
        max_pvsc_ratio,
        max_pvsc_logt,
    ) = plot_phase_profiles(history, left, right)
    lna_period, lna_amp_growth, lna_ke_growth, surface_work = plot_lna_work()

    growth = history["growth"][right]
    ke_growth = history["KE_growth"][right]
    print(f"last completed cycle: {cycle}")
    print(f"nonlinear radius-amplitude growth/cycle: {growth:.6f}")
    print(f"nonlinear KE growth/cycle: {ke_growth:.6f}")
    print(f"square-amplitude prediction: {(1 + growth)**2 - 1:.6f}")
    print(f"LNA fundamental period: {lna_period:.3f} day")
    print(f"LNA amplitude growth/cycle: {lna_amp_growth:.6f}")
    print(f"LNA KE growth/cycle: {lna_ke_growth:.6f}")
    print(f"LNA outermost-cell normalized work: {surface_work:.3e}")
    print(f"atmosphere-work proxy over cycle: {work_proxy:.3e} erg")
    print(f"surface-cell-pressure work over cycle: {work_cell:.3e} erg")
    print(f"absolute surface-cell P dV scale: {work_abs:.3e} erg")
    print(f"maximum sampled pulsation KE: {ke_max:.3e} erg")
    print(f"radiated energy over sampled cycle: {e_rad:.3e} erg")
    print(f"maximum phase change in profile hydrogen fraction: {x_span:.3e}")
    print(f"maximum |Uq|/P in sampled cycle: {max_uq_ratio:.3e}")
    print(
        f"maximum Pvsc/P in sampled cycle: {max_pvsc_ratio:.3e} "
        f"at logT={max_pvsc_logt:.3f}"
    )


if __name__ == "__main__":
    main()
