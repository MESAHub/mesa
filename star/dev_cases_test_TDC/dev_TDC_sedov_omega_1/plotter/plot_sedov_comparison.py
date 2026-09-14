#!/usr/bin/env python3

"""Compare MESA Sedov profiles with the Kamm-Timmes solution."""

from __future__ import annotations

import argparse
import contextlib
import io
from pathlib import Path
import warnings

import matplotlib.pyplot as plt
import mesa_reader as mr
import numpy as np

from globalvars import comvars as gv
from sedov_1d import sed_1d


HERE = Path(__file__).resolve().parent
CASE_DIR = HERE.parent

GAMMA = 1.4
GEOMETRY = 3.0
OMEGA = 1.0
BLAST_ENERGY = 1.464276
RHO_0 = 1.0
P_AMBIENT = 1.0e-6
NZ_BASELINE = 2000


def profile_header(path: Path) -> dict[str, float | str]:
    with path.open(encoding="utf-8") as handle:
        handle.readline()
        names = handle.readline().split()
        values = handle.readline().split()
    header = {}
    for name, value in zip(names, values):
        try:
            header[name] = float(value)
        except ValueError:
            header[name] = value.strip('"')
    return header


def nearest_profile(
    log_dir: Path, target_time: float
) -> tuple[Path, dict[str, float | str]]:
    rows = np.loadtxt(log_dir / "profiles.index", skiprows=1, dtype=int, ndmin=2)
    candidates = []
    for model_number, _, profile_number in rows:
        path = log_dir / f"profile{profile_number}.data"
        header = profile_header(path)
        time_seconds = float(header["time_seconds"])
        candidates.append((abs(time_seconds - target_time), path, header))
    _, path, header = min(candidates, key=lambda item: item[0])
    return path, header


def exact_solution(time_seconds: float, r_max: float) -> dict[str, np.ndarray | float]:
    gv.its = 20
    gv.xgeom = GEOMETRY
    gv.omega = OMEGA
    gv.gamma = GAMMA

    radius = np.linspace(1.0e-4, r_max, 1800)
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        density, velocity, pressure, _, _, _, _, radius = sed_1d(
            time_seconds,
            radius.size,
            radius,
            BLAST_ENERGY,
            RHO_0,
            0.0,
            0.0,
            0.0,
            1.0,
            gv,
        )

    shock_radius = gv.r2
    shock_speed = 2.0 * shock_radius / ((GEOMETRY + 2.0 - OMEGA) * time_seconds)
    rho_upstream = RHO_0 * shock_radius ** (-OMEGA)
    rho_postshock = (GAMMA + 1.0) * rho_upstream / (GAMMA - 1.0)
    velocity_postshock = 2.0 * shock_speed / (GAMMA + 1.0)
    pressure_postshock = 2.0 * rho_upstream * shock_speed**2 / (GAMMA + 1.0)
    return {
        "radius": radius,
        "density": density,
        "velocity": velocity,
        "pressure": pressure,
        "shock_radius": shock_radius,
        "rho_postshock": rho_postshock,
        "velocity_postshock": velocity_postshock,
        "pressure_postshock": pressure_postshock,
    }


def mesa_profile(path: Path) -> dict[str, np.ndarray]:
    profile = mr.MesaData(str(path))
    order = np.argsort(profile.radius_cm)
    return {
        "radius": profile.radius_cm[order],
        "density": profile.density[order],
        "velocity": profile.velocity[order],
        "pressure": profile.pressure[order],
    }


def precursor_diagnostics(
    mesa: dict[str, np.ndarray], exact: dict[str, np.ndarray | float]
) -> dict[str, float]:
    radius = mesa["radius"]
    shock_radius = float(exact["shock_radius"])
    k_shock = int(np.argmin(abs(radius - shock_radius)))
    k_lo = max(0, k_shock - 5)
    k_hi = min(radius.size - 1, k_shock + 5)
    local_dr = np.median(np.diff(radius[k_lo : k_hi + 1]))
    ahead = (radius > shock_radius + 5.0 * local_dr) & (radius < shock_radius + 0.1)
    if not np.any(ahead):
        raise RuntimeError("profile does not extend ahead of the shock")

    ambient_density = RHO_0 * radius[ahead] ** (-OMEGA)
    rho_error = mesa["density"][ahead] / ambient_density - 1.0
    velocity_error = mesa["velocity"][ahead] / float(exact["velocity_postshock"])
    pressure_error = (mesa["pressure"][ahead] - P_AMBIENT) / float(
        exact["pressure_postshock"]
    )

    rho_k = int(np.argmax(abs(rho_error)))
    velocity_k = int(np.argmax(abs(velocity_error)))
    pressure_k = int(np.argmax(abs(pressure_error)))
    ahead_radius = radius[ahead]

    behind = (radius < shock_radius - 5.0 * local_dr) & (radius > shock_radius - 0.2)
    exact_density = np.interp(radius[behind], exact["radius"], exact["density"])
    density_postshock_error = (mesa["density"][behind] - exact_density) / float(
        exact["rho_postshock"]
    )
    density_postshock_k = int(np.argmax(abs(density_postshock_error)))
    behind_radius = radius[behind]
    density_postshock_full_k = np.flatnonzero(behind)[density_postshock_k]
    exact_velocity = np.interp(
        radius[density_postshock_full_k], exact["radius"], exact["velocity"]
    )
    exact_pressure = np.interp(
        radius[density_postshock_full_k], exact["radius"], exact["pressure"]
    )
    target_dr = radius[-1] / NZ_BASELINE
    return {
        "local_dr": local_dr,
        "rho_error": rho_error[rho_k],
        "rho_error_radius": ahead_radius[rho_k],
        "velocity_error": velocity_error[velocity_k],
        "velocity_error_radius": ahead_radius[velocity_k],
        "pressure_error": pressure_error[pressure_k],
        "pressure_error_radius": ahead_radius[pressure_k],
        "density_postshock_error": density_postshock_error[density_postshock_k],
        "density_postshock_error_radius": behind_radius[density_postshock_k],
        "velocity_at_density_error": (
            mesa["velocity"][density_postshock_full_k] - exact_velocity
        )
        / float(exact["velocity_postshock"]),
        "pressure_at_density_error": (
            mesa["pressure"][density_postshock_full_k] - exact_pressure
        )
        / float(exact["pressure_postshock"]),
        "inner_dr_at_density_error": (
            radius[density_postshock_full_k] - radius[density_postshock_full_k - 1]
        )
        / target_dr,
        "outer_dr_at_density_error": (
            radius[density_postshock_full_k + 1] - radius[density_postshock_full_k]
        )
        / target_dr,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--logs",
        nargs="+",
        default=["LOGS_high_res", "LOGS_sedov"],
        help="LOG directories relative to the Sedov case",
    )
    parser.add_argument(
        "--labels",
        nargs="+",
        help="plot labels corresponding to the LOG directories",
    )
    parser.add_argument(
        "--time", type=float, default=0.36, help="target time in seconds"
    )
    parser.add_argument("--output", type=Path, default=HERE / "sedov_comparison.pdf")
    parser.add_argument("--max-long", type=float, default=1.25)
    parser.add_argument("--max-short", type=float, default=2.5)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    plt.style.use(HERE / "mesa.mplstyle")
    if args.labels is not None and len(args.labels) != len(args.logs):
        raise ValueError("--labels must have one entry for each --logs entry")
    labels = args.logs if args.labels is None else args.labels

    selected = []
    for name, label in zip(args.logs, labels):
        log_dir = CASE_DIR / name
        path, header = nearest_profile(log_dir, args.time)
        mesa = mesa_profile(path)
        exact = exact_solution(header["time_seconds"], max(mesa["radius"]))
        selected.append((label, path, header, mesa, exact))

    fields = (
        ("density", r"$\rho\;[\mathrm{g\,cm^{-3}}]$"),
        ("velocity", r"$u\;[\mathrm{cm\,s^{-1}}]$"),
        ("pressure", r"$P\;[\mathrm{erg\,cm^{-3}}]$"),
    )
    fig, axes = plt.subplots(3, 2, figsize=(15, 12), sharex=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    scales = {
        "density": "rho_postshock",
        "velocity": "velocity_postshock",
        "pressure": "pressure_postshock",
    }
    residual_labels = {
        "density": r"$\Delta\rho/\rho_2$",
        "velocity": r"$\Delta u/u_2$",
        "pressure": r"$\Delta P/P_2$",
    }

    for row, (field, ylabel) in enumerate(fields):
        profile_axis = axes[row, 0]
        residual_axis = axes[row, 1]
        for run_number, (name, _, header, mesa, exact) in enumerate(selected):
            color = colors[run_number % len(colors)]
            profile_axis.plot(
                mesa["radius"],
                mesa[field],
                color=color,
                label=f"{name}: model {int(header['model_number'])}",
            )
            profile_axis.axvline(float(exact["shock_radius"]), color=color, alpha=0.25)

            exact_at_mesa = np.interp(mesa["radius"], exact["radius"], exact[field])
            residual = (mesa[field] - exact_at_mesa) / float(exact[scales[field]])
            k_shock = int(np.argmin(abs(mesa["radius"] - float(exact["shock_radius"]))))
            k_lo = max(0, k_shock - 5)
            k_hi = min(mesa["radius"].size - 1, k_shock + 5)
            local_dr = np.median(np.diff(mesa["radius"][k_lo : k_hi + 1]))
            away_from_shock = (
                abs(mesa["radius"] - float(exact["shock_radius"])) > 5.0 * local_dr
            )
            residual_axis.plot(
                mesa["radius"][away_from_shock],
                residual[away_from_shock],
                color=color,
                label=name,
            )
        reference = selected[0][4]
        profile_axis.plot(
            reference["radius"],
            reference[field],
            color="black",
            linestyle="--",
            label="Kamm-Timmes",
        )
        profile_axis.set_ylabel(ylabel)
        profile_axis.legend(fontsize=12)
        residual_axis.axhline(0.0, color="black", linewidth=1.0)
        residual_axis.set_ylabel(residual_labels[field])
        residual_axis.legend(fontsize=12)

    axes[-1, 0].set_xlabel(r"$r\;[\mathrm{cm}]$")
    axes[-1, 1].set_xlabel(r"$r\;[\mathrm{cm}]$")
    axes[-1, 0].set_xlim(0.4, 0.75)
    axes[-1, 1].set_xlim(0.4, 0.75)
    fig.savefig(args.output)
    fig.savefig(args.output.with_suffix(".png"), dpi=200)

    mesh_fig, mesh_axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
    for run_number, (name, _, _, mesa, exact) in enumerate(selected):
        color = colors[run_number % len(colors)]
        exact_density = np.interp(mesa["radius"], exact["radius"], exact["density"])
        density_residual = (mesa["density"] - exact_density) / float(
            exact["rho_postshock"]
        )
        mesh_axes[0].plot(mesa["radius"], density_residual, color=color, label=name)

        dr = np.diff(mesa["radius"])
        rmid = 0.5 * (mesa["radius"][1:] + mesa["radius"][:-1])
        target_dr = mesa["radius"][-1] / NZ_BASELINE
        mesh_axes[1].plot(rmid, dr / target_dr, color=color, label=name)
        for axis in mesh_axes:
            axis.axvline(float(exact["shock_radius"]), color=color, alpha=0.25)

    mesh_axes[0].axhline(0.0, color="black", linewidth=1.0)
    mesh_axes[0].set_ylabel(r"$\Delta\rho/\rho_2$")
    mesh_axes[0].legend(fontsize=12)
    mesh_axes[1].axhline(args.max_long, color="black", linestyle=":", linewidth=1.0)
    mesh_axes[1].axhline(
        1.0 / args.max_short, color="black", linestyle=":", linewidth=1.0
    )
    mesh_axes[1].set_ylabel(r"$\Delta r/\Delta r_{\rm target}$")
    mesh_axes[1].set_xlabel(r"$r\;[\mathrm{cm}]$")
    mesh_axes[1].legend(fontsize=12)
    mesh_axes[1].set_xlim(0.45, 0.65)
    mesh_output = args.output.with_name(args.output.stem + "_mesh")
    mesh_fig.savefig(mesh_output.with_suffix(".pdf"))
    mesh_fig.savefig(mesh_output.with_suffix(".png"), dpi=200)

    for name, path, header, mesa, exact in selected:
        diagnostic = precursor_diagnostics(mesa, exact)
        print(
            f"{name}: {path.name}, model={int(header['model_number'])}, "
            f"t={header['time_seconds']:.8f} s, "
            f"r_shock={float(exact['shock_radius']):.8f} cm"
        )
        print(
            "  upstream maxima beyond five local cells: "
            f"dr={diagnostic['local_dr']:.3e} cm, "
            f"dRho/Rho={diagnostic['rho_error']:.3e} "
            f"at {diagnostic['rho_error_radius']:.6f} cm, "
            f"u/u2={diagnostic['velocity_error']:.3e} "
            f"at {diagnostic['velocity_error_radius']:.6f} cm, "
            f"dP/P2={diagnostic['pressure_error']:.3e} "
            f"at {diagnostic['pressure_error_radius']:.6f} cm"
        )
        print(
            "  postshock density maximum away from the front: "
            f"dRho/Rho2={diagnostic['density_postshock_error']:.3e} "
            f"at {diagnostic['density_postshock_error_radius']:.6f} cm, "
            f"du/u2={diagnostic['velocity_at_density_error']:.3e}, "
            f"dP/P2={diagnostic['pressure_at_density_error']:.3e}, "
            "adjacent dr/target="
            f"({diagnostic['inner_dr_at_density_error']:.3f}, "
            f"{diagnostic['outer_dr_at_density_error']:.3f})"
        )


if __name__ == "__main__":
    main()
