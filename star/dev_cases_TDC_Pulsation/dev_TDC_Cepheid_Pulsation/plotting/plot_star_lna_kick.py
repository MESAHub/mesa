#!/usr/bin/env python3

from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from mesa_reader import MesaData


PLOT_DIR = Path(__file__).resolve().parent
CASE_DIR = PLOT_DIR.parent
LOGS_DIR = CASE_DIR / "LOGS"
LNA_DIR = CASE_DIR / "LNA"


def read_integer_control(name):
    text = (CASE_DIR / "inlist_pulses").read_text()
    match = re.search(rf"^\s*{name}\s*=\s*([-+]?\d+)", text, re.MULTILINE)
    if match is None:
        raise ValueError(f"could not read {name} from inlist_pulses")
    return int(match.group(1))


def read_metadata(path):
    metadata = {}
    with path.open() as stream:
        for line in stream:
            if not line.startswith("#"):
                break
            fields = line[1:].split()
            if len(fields) == 2:
                metadata[fields[0]] = fields[1]
    return metadata


def load_eigenfunction(mode):
    path = LNA_DIR / f"star_LNA_eigenfunction_{mode}.data"
    return np.loadtxt(path), read_metadata(path)


def displacement_shapes(data):
    radius = data[:, 3]
    dlnr = data[:, 6] + 1j * data[:, 7]
    relative_dlnr = dlnr / dlnr[0]

    phase = np.angle(relative_dlnr)
    sign = np.ones_like(phase)
    sign[np.abs(phase) > 0.5 * np.pi] = -1.0
    current = sign * np.abs(radius * dlnr) / np.abs(radius[0] * dlnr[0])

    # This is the real displacement slice used by the former GYRE kick.
    real_slice = np.real(radius * dlnr / (radius[0] * dlnr[0]))
    return current, real_slice, phase


def count_nodes(shape, threshold=0.0, mask=None):
    crossings = shape[:-1] * shape[1:] < 0.0
    crossings &= np.maximum(np.abs(shape[:-1]), np.abs(shape[1:])) > threshold
    if mask is not None:
        crossings &= mask[:-1]
    return int(np.count_nonzero(crossings))


def kinematic_quality(data, metadata):
    sigma = float(metadata["sigma_real"]) + 1j * float(metadata["sigma_imag"])
    radius = data[:, 3]
    dlnr = data[:, 6] + 1j * data[:, 7]
    velocity = data[:, 8] + 1j * data[:, 9]

    velocity_term = velocity / radius
    radius_term = sigma * dlnr
    residual = velocity_term - radius_term
    scale = np.abs(velocity_term) + np.abs(radius_term) + np.finfo(float).tiny
    balance_residual = np.abs(residual) / scale
    amplitude_ratio = np.abs(velocity_term) / np.maximum(
        np.abs(radius_term), np.finfo(float).tiny
    )
    return balance_residual, amplitude_ratio


def profile_number_for_model(model_number):
    index = np.loadtxt(LOGS_DIR / "profiles.index", dtype=int, skiprows=1)
    row = index[np.argmin(np.abs(index[:, 0] - model_number))]
    return int(row[2])


def profile_for_model(model_number):
    number = profile_number_for_model(model_number)
    return MesaData(str(LOGS_DIR / f"profile{number}.data"))


def plot_history(history, kick_model):
    kick_index = np.flatnonzero(history.model_number == kick_model)[0]
    kick_day = history.day[kick_index]
    time = history.day - kick_day
    use = time >= -50.0

    figure, axes = plt.subplots(3, 1, figsize=(9.0, 8.2), sharex=True)
    axes[0].semilogy(time[use], history.time_step_sec[use], color="#1f4e79")
    axes[0].set_ylabel("Timestep (s)")

    axes[1].plot(time[use], history.v_surf_km_s[use], color="#b33a3a", label=r"$v_{\rm surf}$")
    axes[1].plot(
        time[use],
        history.v_div_csound_surf[use],
        color="#2e7d32",
        label=r"$v_{\rm surf}/c_s$",
    )
    axes[1].axhline(0.0, color="0.75", linewidth=0.8)
    axes[1].set_ylabel("Surface velocity")
    axes[1].legend(frameon=False, ncol=2)

    retry_increment = np.diff(history.num_retries, prepend=history.num_retries[0])
    axes[2].step(time[use], history.num_retries[use], where="post", color="#6a3d9a", label="cumulative retries")
    axes[2].scatter(
        time[use][retry_increment[use] > 0],
        history.num_retries[use][retry_increment[use] > 0],
        s=14,
        color="#d17c00",
        label="accepted after retry",
        zorder=3,
    )
    axes[2].set_ylabel("Retries")
    axes[2].set_xlabel("Days from Star LNA kick")
    axes[2].legend(frameon=False)

    for axis in axes:
        axis.axvline(0.0, color="0.2", linestyle="--", linewidth=0.9)
        axis.grid(alpha=0.2)

    figure.suptitle("Nonlinear response to the Star LNA velocity kick")
    figure.tight_layout()
    figure.savefig(PLOT_DIR / "kick_history.png", dpi=180)
    plt.close(figure)


def plot_kick_profile(kick_model, selected_mode):
    before = profile_for_model(kick_model)
    after = profile_for_model(kick_model + 10)
    later = profile_for_model(kick_model + 200)
    data, metadata = load_eigenfunction(selected_mode)
    current, real_slice, phase = displacement_shapes(data)
    requested_velocity = float(metadata["star_LNA_kick_vsurf_km_per_sec"])

    figure, axes = plt.subplots(2, 2, figsize=(11.0, 7.8))

    axes[0, 0].plot(before.logT, requested_velocity * current, color="#b33a3a", label="current kick")
    axes[0, 0].plot(before.logT, requested_velocity * real_slice, color="#1f4e79", linestyle="--", label="real displacement slice")
    axes[0, 0].axhline(0.0, color="0.75", linewidth=0.8)
    axes[0, 0].set_ylabel("Velocity shape (km/s)")
    axes[0, 0].set_title("Selected eigenfunction")
    axes[0, 0].legend(frameon=False)

    outer = before.logT <= 4.2
    axes[0, 1].plot(before.logT[outer], requested_velocity * current[outer], color="#b33a3a", label="imposed at model 200")
    axes[0, 1].plot(after.logT[outer], after.vel_km_per_s[outer], color="#d17c00", label="model 210")
    axes[0, 1].plot(later.logT[outer], later.vel_km_per_s[outer], color="#2e7d32", label="model 400")
    axes[0, 1].axhline(0.0, color="0.75", linewidth=0.8)
    axes[0, 1].set_ylabel("Velocity (km/s)")
    axes[0, 1].set_title("Near-surface response")
    axes[0, 1].legend(frameon=False)

    axes[1, 0].plot(before.logT, phase / np.pi, color="#6a3d9a")
    axes[1, 0].axhline(0.5, color="0.65", linestyle=":", linewidth=0.8)
    axes[1, 0].axhline(-0.5, color="0.65", linestyle=":", linewidth=0.8)
    axes[1, 0].set_ylabel(r"Phase of $\delta\ln R$ ($\pi$ rad)")
    axes[1, 0].set_xlabel(r"$\log_{10}(T/{\rm K})$")

    axes[1, 1].plot(before.logT, before.vel_km_per_s, color="0.45", label="model 200")
    axes[1, 1].plot(after.logT, after.vel_km_per_s, color="#d17c00", label="model 210")
    axes[1, 1].plot(later.logT, later.vel_km_per_s, color="#2e7d32", label="model 400")
    axes[1, 1].axhline(0.0, color="0.75", linewidth=0.8)
    axes[1, 1].set_ylabel("Velocity (km/s)")
    axes[1, 1].set_xlabel(r"$\log_{10}(T/{\rm K})$")
    axes[1, 1].set_title("Profile evolution")
    axes[1, 1].legend(frameon=False)

    for axis in axes.flat:
        axis.grid(alpha=0.2)

    figure.suptitle(f"Star LNA kick from selected mode {selected_mode}")
    figure.tight_layout()
    figure.savefig(PLOT_DIR / "kick_profile.png", dpi=180)
    plt.close(figure)


def collect_mode_quality(background):
    period_growth = np.loadtxt(LNA_DIR / "star_LNA_period_growth.data")
    records = []
    for row in period_growth:
        mode = int(row[0])
        data, metadata = load_eigenfunction(mode)
        current, real_slice, _ = displacement_shapes(data)
        balance_residual, amplitude_ratio = kinematic_quality(data, metadata)
        records.append(
            {
                "mode": mode,
                "period": row[6],
                "growth": row[9],
                "selected": bool(int(row[2])),
                "nodes": count_nodes(real_slice),
                "nodes_001": count_nodes(real_slice, threshold=1.0e-2),
                "outer_nodes": count_nodes(
                    real_slice, threshold=1.0e-3, mask=background.logT < 5.0
                ),
                "median_balance_residual": float(np.median(balance_residual)),
                "median_amplitude_ratio": float(np.median(amplitude_ratio)),
                "current_min": float(np.min(current)),
                "current_max": float(np.max(current)),
            }
        )
    return records


def plot_mode_quality(records):
    mode = np.array([record["mode"] for record in records])
    nodes = np.array([record["nodes_001"] for record in records])
    ratio = np.array([record["median_amplitude_ratio"] for record in records])
    growth = np.array([record["growth"] for record in records])
    selected = np.array([record["selected"] for record in records])

    figure, axes = plt.subplots(3, 1, figsize=(9.0, 8.4), sharex=True)
    axes[0].bar(mode, nodes, color=np.where(selected, "#b33a3a", "#1f4e79"))
    axes[0].set_ylabel("Nodes above 1%\nof surface amplitude")

    axes[1].semilogy(mode, ratio, "o-", color="#6a3d9a")
    axes[1].axhline(1.0, color="0.4", linestyle="--", linewidth=0.8)
    axes[1].set_ylabel(r"Median $|\delta v/r|/|\sigma\delta\ln R|$")

    axes[2].plot(mode, growth, "o-", color="#2e7d32")
    axes[2].axhline(0.0, color="0.65", linewidth=0.8)
    axes[2].scatter(mode[selected], growth[selected], s=55, facecolors="none", edgecolors="#b33a3a", linewidths=1.5, label="kick mode")
    axes[2].set_ylabel("log KE per cycle")
    axes[2].set_xlabel("One-based selected mode")
    axes[2].legend(frameon=False)

    for axis in axes:
        axis.set_xticks(mode)
        axis.grid(alpha=0.2)

    figure.suptitle("Selected-mode quality checks")
    figure.tight_layout()
    figure.savefig(PLOT_DIR / "mode_quality.png", dpi=180)
    plt.close(figure)


def write_diagnostics(history, kick_model, selected_mode, records):
    kick_index = np.flatnonzero(history.model_number == kick_model)[0]
    next_index = np.flatnonzero(history.model_number == kick_model + 1)[0]
    background = profile_for_model(kick_model)
    data, metadata = load_eigenfunction(selected_mode)
    current, real_slice, _ = displacement_shapes(data)
    balance_residual, amplitude_ratio = kinematic_quality(data, metadata)
    crossings = np.flatnonzero(current[:-1] * current[1:] < 0.0)
    first_crossing = crossings[0]
    sigma_imag = float(metadata["sigma_imag"])
    period_days = 2.0 * np.pi / sigma_imag / 86400.0
    requested_velocity = float(metadata["star_LNA_kick_vsurf_km_per_sec"])

    lines = [
        "Star LNA kick diagnostics",
        "=========================",
        "",
        f"kick model: {kick_model}",
        f"selected one-based mode: {selected_mode}",
        f"selected zero-based control index: {int(metadata['control_mode_index'])}",
        f"selected period: {period_days:.8g} days",
        f"requested surface velocity: {requested_velocity:.8g} km/s",
        f"surface Mach number after kick: {history.v_div_csound_surf[next_index]:.8g}",
        f"accepted timestep before kick: {history.time_step_sec[kick_index]:.8g} s",
        f"first accepted timestep after kick: {history.time_step_sec[next_index]:.8g} s",
        f"timestep reduction factor: {history.time_step_sec[kick_index]/history.time_step_sec[next_index]:.8g}",
        f"retries before first accepted post-kick step: {int(history.num_retries[next_index]-history.num_retries[kick_index])}",
        f"elapsed simulated time after kick: {history.day[-1]-history.day[kick_index]:.8g} days",
        "",
        "Selected displacement shape",
        "---------------------------",
        f"all real-slice sign changes: {count_nodes(real_slice)}",
        f"sign changes above 1 percent of surface amplitude: {count_nodes(real_slice, threshold=1.0e-2)}",
        f"sign changes above 0.1 percent for logT < 5: {count_nodes(real_slice, threshold=1.0e-3, mask=background.logT < 5.0)}",
        f"first current-kick sign change: zones {first_crossing+1} to {first_crossing+2}",
        f"first sign-change logT interval: {background.logT[first_crossing]:.8g} to {background.logT[first_crossing+1]:.8g}",
        f"first sign-change velocities: "
        f"{requested_velocity*current[first_crossing]:.8g} to "
        f"{requested_velocity*current[first_crossing+1]:.8g} km/s",
        f"maximum difference from the former GYRE real displacement slice: {np.max(np.abs(current-real_slice)):.8g}",
        "",
        "Kinematic eigenvector check",
        "---------------------------",
        "The radius row requires delta_v/r = sigma*delta_lnR on the static background.",
        f"median normalized balance residual: {np.median(balance_residual):.8g}",
        f"maximum normalized balance residual: {np.max(balance_residual):.8g}",
        f"median amplitude ratio |delta_v/r|/|sigma*delta_lnR|: {np.median(amplitude_ratio):.8g}",
        "",
        "Selected mode table",
        "-------------------",
        "mode period_days logKE nodes nodes_gt_0.01 outer_nodes median_balance median_amplitude_ratio",
    ]
    for record in records:
        lines.append(
            f"{record['mode']:4d} {record['period']:12.5e} {record['growth']:12.5e} "
            f"{record['nodes']:5d} {record['nodes_001']:5d} {record['outer_nodes']:5d} "
            f"{record['median_balance_residual']:12.5e} {record['median_amplitude_ratio']:12.5e}"
        )
    (PLOT_DIR / "kick_diagnostics.txt").write_text("\n".join(lines) + "\n")


def main():
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
        }
    )

    kick_model = read_integer_control("star_LNA_model_number")
    selected_mode = read_integer_control("star_LNA_kick_mode_1")
    history = MesaData(str(LOGS_DIR / "history.data"))
    background = profile_for_model(kick_model)
    records = collect_mode_quality(background)

    plot_history(history, kick_model)
    plot_kick_profile(kick_model, selected_mode)
    plot_mode_quality(records)
    write_diagnostics(history, kick_model, selected_mode, records)


if __name__ == "__main__":
    main()
