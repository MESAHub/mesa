#!/usr/bin/env python3
"""Audit per-step energy errors and mesh changes in a profile sequence."""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mesa_reader as mr
import numpy as np


CASE = Path(__file__).resolve().parents[1]
LOGS = CASE / "LOGS"
OUT = Path(__file__).resolve().parent
G_CGS = 6.67430e-8


def cell_masses(profile: mr.MesaData) -> np.ndarray:
    """Return cell masses using MESA's envelope-mass normalization of dq."""
    xmstar = float(profile.star_mass) * float(profile.msun) - float(profile.M_center)
    return np.asarray(profile.dq) * xmstar


def history_names(path: Path) -> list[str]:
    with path.open() as stream:
        for number, line in enumerate(stream, start=1):
            if number == 6:
                names = line.split()
                break
        else:
            raise RuntimeError(f"missing history header in {path}")
    last = subprocess.run(
        ["tail", "-n", "1", str(path)], check=True, capture_output=True, text=True
    ).stdout.split()
    if len(names) != len(last):
        raise RuntimeError(
            f"history header has {len(names)} columns, but data have {len(last)}"
        )
    return names


def read_history_tail(path: Path, models: set[int], nlines: int = 5000):
    names = history_names(path)
    wanted = [
        "model_number",
        "num_zones",
        "time_step_sec",
        "log_rel_error_in_energy_conservation",
        "log_rel_run_E_err",
    ]
    missing = [name for name in wanted if name not in names]
    if missing:
        raise RuntimeError(f"missing history columns: {', '.join(missing)}")
    indices = [names.index(name) for name in wanted]
    rows = {}
    process = subprocess.Popen(
        ["tail", "-n", str(nlines), str(path)], stdout=subprocess.PIPE, text=True
    )
    assert process.stdout is not None
    for line in process.stdout:
        fields = line.split()
        if len(fields) != len(names):
            continue
        try:
            model = int(float(fields[indices[0]]))
            if model in models:
                rows[model] = {
                    name: float(fields[index])
                    for name, index in zip(wanted, indices)
                }
        except ValueError:
            continue
    if process.wait() != 0:
        raise RuntimeError("tail failed while reading history.data")
    return rows


def infer_abs_total_energy(profile: mr.MesaData) -> float:
    error = np.asarray(profile.ergs_error)
    log_rel = np.asarray(profile.log_rel_E_err)
    good = (
        np.isfinite(error)
        & np.isfinite(log_rel)
        & (error != 0.0)
        & (log_rel > -300.0)
    )
    estimates = np.abs(error[good]) / np.power(10.0, log_rel[good])
    return float(np.median(estimates)) if estimates.size else np.nan


def integrated_energies(profile: mr.MesaData) -> dict[str, float]:
    """Reproduce the nonrotating u-flag quadrature in star_utils.f90."""
    dm = cell_masses(profile)
    radius = np.asarray(profile.radius) * float(profile.rsun)
    mass = np.asarray(profile.mass) * float(profile.msun)
    grav_face = -G_CGS * mass / radius
    grav_inner = np.empty_like(grav_face)
    grav_inner[:-1] = grav_face[1:]
    grav_inner[-1] = (
        -G_CGS * float(profile.M_center) / float(profile.R_center)
        if float(profile.R_center) > 0.0
        else 0.0
    )
    specific_pe = 0.5 * (grav_face + grav_inner)
    specific_ke = 0.5 * np.square(np.asarray(profile.velocity))

    mlt_vc = np.asarray(profile.mlt_vc)
    specific_turb = 0.75 * np.square(mlt_vc)
    specific_turb[:-1] += 0.75 * np.square(mlt_vc[1:])

    internal = float(np.sum(dm * np.asarray(profile.energy)))
    potential = float(np.sum(dm * specific_pe))
    kinetic = float(np.sum(dm * specific_ke))
    turbulent = float(np.sum(dm * specific_turb))
    return {
        "internal": internal,
        "potential": potential,
        "kinetic": kinetic,
        "turbulent": turbulent,
        "total": internal + potential + kinetic + turbulent,
    }


def mechanical_balance_cells(
    previous: mr.MesaData, profile: mr.MesaData, dt_sec: float
) -> np.ndarray:
    """Approximate the zone contribution outside the internal-energy solve."""
    dm = cell_masses(profile)

    def specific_pe(data: mr.MesaData) -> np.ndarray:
        radius = np.asarray(data.radius) * float(data.rsun)
        mass = np.asarray(data.mass) * float(data.msun)
        grav_face = -G_CGS * mass / radius
        grav_inner = np.empty_like(grav_face)
        grav_inner[:-1] = grav_face[1:]
        grav_inner[-1] = -G_CGS * float(data.M_center) / float(data.R_center)
        return 0.5 * (grav_face + grav_inner)

    delta_pe = dm * (specific_pe(profile) - specific_pe(previous))
    delta_ke = 0.5 * dm * (
        np.square(np.asarray(profile.velocity))
        - np.square(np.asarray(previous.velocity))
    )
    work = dt_sec * dm * np.asarray(profile.dwork_dm)
    eq_face = np.asarray(profile.Eq)
    eq_cell = 0.5 * eq_face
    eq_cell[:-1] += 0.5 * eq_face[1:]
    eq = dt_sec * dm * eq_cell
    return delta_pe + delta_ke - work + eq


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("first", nargs="?", type=int, default=301)
    parser.add_argument("last", nargs="?", type=int, default=314)
    args = parser.parse_args()

    profiles = []
    for number in range(args.first, args.last + 1):
        profiles.append((number, mr.MesaData(str(LOGS / f"profile{number}.data"))))

    models = {int(profile.model_number) for _, profile in profiles}
    history = read_history_tail(LOGS / "history.data", models)

    rows = []
    previous_dq = None
    previous_profile = None
    previous_energies = None
    for number, profile in profiles:
        model = int(profile.model_number)
        error = np.asarray(profile.ergs_error)
        dq = np.asarray(profile.dq)
        abs_total_energy = infer_abs_total_energy(profile)
        relative_error = error / abs_total_energy
        imax = int(np.nanargmax(np.abs(relative_error)))
        if previous_dq is None or previous_dq.size != dq.size:
            dq_l1_change = np.nan
            dq_max_change = np.nan
        else:
            delta_dq = dq - previous_dq
            dq_l1_change = float(np.sum(np.abs(delta_dq)))
            dq_max_change = float(np.max(np.abs(delta_dq)))
        previous_dq = dq.copy()
        hrow = history.get(model, {})
        energies = integrated_energies(profile)
        history_log_rel_error = hrow.get(
            "log_rel_error_in_energy_conservation", np.nan
        )
        global_error_abs = (
            abs_total_energy * 10.0**history_log_rel_error
            if np.isfinite(history_log_rel_error)
            else np.nan
        )
        dt_sec = hrow.get(
            "time_step_sec", float(profile.time_step) * 365.25 * 86400.0
        )
        if previous_profile is None:
            delta_internal = np.nan
            delta_potential = np.nan
            delta_kinetic = np.nan
            delta_turbulent = np.nan
            delta_total = np.nan
            radiation = np.nan
            surface_work = np.nan
            center_work = np.nan
            reconstructed_error = np.nan
            eddy_transfer_mismatch = np.nan
            integrated_work = np.nan
            integrated_eq = np.nan
            cell_energy_closure = np.nan
            mechanical_closure = np.nan
        else:
            delta_internal = energies["internal"] - previous_energies["internal"]
            delta_potential = energies["potential"] - previous_energies["potential"]
            delta_kinetic = energies["kinetic"] - previous_energies["kinetic"]
            delta_turbulent = energies["turbulent"] - previous_energies["turbulent"]
            delta_total = energies["total"] - previous_energies["total"]
            luminosity_surface = float(profile.luminosity[0]) * float(profile.lsun)
            radiation = dt_sec * (luminosity_surface - float(profile.L_center))
            surface_speed = (
                (float(profile.radius[0]) - float(previous_profile.radius[0]))
                * float(profile.rsun)
                / dt_sec
            )
            surface_pressure = 10.0 ** float(profile.logP[0])
            surface_radius = float(profile.radius[0]) * float(profile.rsun)
            surface_work = (
                dt_sec
                * 4.0
                * np.pi
                * surface_radius**2
                * surface_pressure
                * surface_speed
            )
            center_speed = (
                (float(profile.R_center) - float(previous_profile.R_center)) / dt_sec
            )
            center_pressure = 10.0 ** float(profile.logP[-1])
            center_work = (
                dt_sec
                * 4.0
                * np.pi
                * float(profile.R_center) ** 2
                * center_pressure
                * center_speed
            )
            reconstructed_error = delta_total + radiation + surface_work - center_work
            dm = cell_masses(profile)
            eq_face = np.asarray(profile.Eq)
            eq_cell = 0.5 * eq_face
            eq_cell[:-1] += 0.5 * eq_face[1:]
            integrated_work = dt_sec * float(
                np.sum(dm * np.asarray(profile.dwork_dm))
            )
            integrated_eq = dt_sec * float(np.sum(dm * eq_cell))
            cell_energy_closure = (
                delta_internal
                + delta_turbulent
                + radiation
                + integrated_work
                - integrated_eq
            )
            mechanical_closure = (
                delta_potential
                + delta_kinetic
                - integrated_work
                + integrated_eq
                + surface_work
                - center_work
            )
            velocity_mid = 0.5 * (
                np.asarray(profile.velocity) + np.asarray(previous_profile.velocity)
            )
            eddy_transfer_mismatch = dt_sec * float(
                np.sum(dm * (eq_cell + velocity_mid * np.asarray(profile.Uq)))
            )
        rows.append(
            {
                "profile_number": number,
                "model_number": model,
                "num_zones": int(profile.num_zones),
                "dt_sec": dt_sec,
                "history_log_rel_E_err": history_log_rel_error,
                "history_log_rel_run_E_err": hrow.get("log_rel_run_E_err", np.nan),
                "abs_total_energy_for_rel_error_erg": abs_total_energy,
                "global_error_abs_erg": global_error_abs,
                "cell_residual_signed_rel": float(np.sum(relative_error)),
                "cell_residual_abs_rel": float(np.sum(np.abs(relative_error))),
                "max_cell_residual_rel": float(relative_error[imax]),
                "max_error_zone": int(profile.zone[imax]),
                "max_error_mass_Msun": float(profile.mass[imax]),
                "max_error_q": float(profile.q[imax]),
                "max_error_logR": float(profile.logR[imax]),
                "max_error_logT": float(profile.logT[imax]),
                "max_error_logRho": float(profile.logRho[imax]),
                "max_error_velocity_cm_s": float(profile.velocity[imax]),
                "max_error_dwork_dm": float(profile.dwork_dm[imax]),
                "max_error_Eq": float(profile.Eq[imax]),
                "max_error_Uq": float(profile.Uq[imax]),
                "dq_l1_change": dq_l1_change,
                "dq_max_change": dq_max_change,
                "total_energy_erg": energies["total"],
                "delta_internal_erg": delta_internal,
                "delta_potential_erg": delta_potential,
                "delta_kinetic_erg": delta_kinetic,
                "delta_turbulent_erg": delta_turbulent,
                "delta_turbulent_div_global_error_abs": (
                    delta_turbulent / global_error_abs
                    if global_error_abs > 0.0
                    else np.nan
                ),
                "delta_total_erg": delta_total,
                "radiation_loss_erg": radiation,
                "surface_work_erg": surface_work,
                "center_work_erg": center_work,
                "reconstructed_error_rel": reconstructed_error / abs_total_energy,
                "eddy_transfer_mismatch_erg": eddy_transfer_mismatch,
                "eddy_transfer_mismatch_rel": eddy_transfer_mismatch / abs_total_energy,
                "integrated_work_erg": integrated_work,
                "integrated_Eq_erg": integrated_eq,
                "cell_energy_closure_rel": cell_energy_closure / abs_total_energy,
                "mechanical_closure_rel": mechanical_closure / abs_total_energy,
            }
        )
        previous_profile = profile
        previous_energies = energies

    csv_path = OUT / f"energy_error_profiles_{args.first}_{args.last}.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    model = np.asarray([row["model_number"] for row in rows])
    history_log = np.asarray([row["history_log_rel_E_err"] for row in rows])
    residual_signed = np.asarray([row["cell_residual_signed_rel"] for row in rows])
    residual_abs = np.asarray([row["cell_residual_abs_rel"] for row in rows])
    dq_l1 = np.asarray([row["dq_l1_change"] for row in rows])
    reconstructed = np.asarray([row["reconstructed_error_rel"] for row in rows])
    global_error_abs = np.asarray([row["global_error_abs_erg"] for row in rows])
    delta_turbulent = np.asarray([row["delta_turbulent_erg"] for row in rows])

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    ax = axes[0, 0]
    ax.plot(model, history_log, "o-", label="global step error")
    ax.plot(
        model,
        np.log10(np.maximum(np.abs(residual_signed), 1e-30)),
        "s--",
        label="signed sum of cell residuals",
    )
    ax.plot(
        model,
        np.log10(np.maximum(residual_abs, 1e-30)),
        "^:",
        label="sum of absolute cell residuals",
    )
    ax.plot(
        model,
        np.log10(np.maximum(np.abs(reconstructed), 1e-30)),
        "d-.",
        label="profile-reconstructed balance",
    )
    ax.set_ylabel(r"$\log_{10}|\Delta E/E|$")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)

    ax = axes[0, 1]
    ax.semilogy(model[1:], np.maximum(dq_l1[1:], 1e-30), "o-")
    ax.set_ylabel(r"$\sum_k|dq_k^{n}-dq_k^{n-1}|$")
    ax.grid(alpha=0.25)

    selected_indices = [1, len(profiles) // 2, len(profiles) - 1]
    ax = axes[1, 0]
    for index in selected_indices:
        number, profile = profiles[index]
        _, previous = profiles[index - 1]
        balance = mechanical_balance_cells(previous, profile, rows[index]["dt_sec"])
        ax.plot(
            profile.mass,
            balance / infer_abs_total_energy(profile),
            label=f"profile {number}",
        )
    ax.set_xlabel(r"enclosed mass ($M_\odot$)")
    ax.set_ylabel(r"approx. mechanical balance per zone / $|E|$")
    ax.set_yscale("symlog", linthresh=1e-8)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)

    ax = axes[1, 1]
    for index in selected_indices:
        number, profile = profiles[index]
        _, previous = profiles[index - 1]
        balance = mechanical_balance_cells(previous, profile, rows[index]["dt_sec"])
        cumulative = np.cumsum(balance) / infer_abs_total_energy(profile)
        ax.plot(profile.mass, cumulative, label=f"profile {number}")
    ax.set_xlabel(r"enclosed mass ($M_\odot$)")
    ax.set_ylabel(r"cumulative approximate mechanical balance / $|E|$")
    ax.grid(alpha=0.25)

    figure_path = OUT / f"energy_error_profiles_{args.first}_{args.last}.png"
    fig.savefig(figure_path, dpi=180)
    plt.close(fig)

    worst_global = int(np.nanargmax(history_log))
    worst_residual = int(np.nanargmax(residual_abs))
    max_mesh_change = float(np.nanmax(dq_l1[1:]))
    valid_turbulent = np.isfinite(global_error_abs) & np.isfinite(delta_turbulent)
    turbulent_error_correlation = np.corrcoef(
        global_error_abs[valid_turbulent], delta_turbulent[valid_turbulent]
    )[0, 1]
    turbulent_error_ratio = np.nanmedian(
        delta_turbulent[valid_turbulent] / global_error_abs[valid_turbulent]
    )
    summary_path = OUT / f"energy_error_profiles_{args.first}_{args.last}.txt"
    with summary_path.open("w") as stream:
        stream.write(f"Profiles {args.first}--{args.last}\n")
        stream.write(
            f"Worst global step: model {model[worst_global]}, "
            f"log10(|dE/E|)={history_log[worst_global]:.6e}\n"
        )
        stream.write(
            f"Largest absolute cell-residual sum: model {model[worst_residual]}, "
            f"sum|e_k|/|E|={residual_abs[worst_residual]:.6e}\n"
        )
        stream.write(f"Maximum dq L1 change between profiles: {max_mesh_change:.6e}\n")
        stream.write(
            "Median profile-reconstructed |dE/E|: "
            f"{np.nanmedian(np.abs(reconstructed)):.6e}\n"
        )
        stream.write(
            "Median eddy-transfer mismatch |dE/E|: "
            f"{np.nanmedian(np.abs([row['eddy_transfer_mismatch_rel'] for row in rows])):.6e}\n"
        )
        stream.write(
            "Correlation of |global error| with accepted-state delta Eturb: "
            f"{turbulent_error_correlation:.6e}\n"
        )
        stream.write(
            "Median accepted-state delta Eturb / |global error|: "
            f"{turbulent_error_ratio:.6e}\n"
        )
        stream.write(f"CSV: {csv_path.name}\n")
        stream.write(f"Figure: {figure_path.name}\n")

    print(summary_path.read_text(), end="")


if __name__ == "__main__":
    main()
