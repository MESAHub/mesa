#!/usr/bin/env python3
"""Plot the cell energy-balance diagnostic from the Noh test."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_MODELS = (1, 10, 50, 100, 250, 500, 750, 1000)
GAMMA = 1.66667
NOH_DENSITY = 64.0
NOH_PRESSURE = 64.0 / 3.0
NOH_ENERGY = 0.5
NOH_ENTROPY_PROXY = NOH_PRESSURE / NOH_DENSITY**GAMMA
DATA_COLUMNS = (
    "model",
    "k",
    "age_s",
    "dt_s",
    "r_cm",
    "dm_g",
    "rho",
    "u",
    "P",
    "energy",
    "dU_erg",
    "dK_erg",
    "dW_erg",
    "balance_erg",
    "equation_error_erg",
    "balance_minus_error_erg",
    "relative_balance",
)


def parse_args() -> argparse.Namespace:
    case_dir = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description="Plot noh_energy_balance.txt without loading the full file into memory."
    )
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        default=case_dir / "noh_energy_balance.txt",
        help="energy-balance file (default: ../noh_energy_balance.txt)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=case_dir / "plotter" / "plots",
        help="directory for PNG and CSV output",
    )
    parser.add_argument(
        "--models",
        nargs="+",
        type=int,
        default=list(DEFAULT_MODELS),
        help="models for detailed radial plots",
    )
    parser.add_argument(
        "--radius-max",
        type=float,
        default=None,
        help="maximum radius for detailed profiles; defaults to 2.5 shock radii",
    )
    return parser.parse_args()


def parse_model_header(line: str) -> dict[str, float | int] | None:
    fields = line.split()
    if len(fields) != 15 or fields[:2] != ["#", "model"]:
        return None
    return {
        "model": int(fields[2]),
        "nz": int(fields[4]),
        "age_s": float(fields[6]),
        "total_dU_erg": float(fields[8]),
        "total_dK_erg": float(fields[10]),
        "total_dW_erg": float(fields[12]),
        "total_balance_erg": float(fields[14]),
    }


def new_aggregate(header: dict[str, float | int]) -> dict[str, float | int]:
    return {
        **header,
        "rows_read": 0,
        "sum_equation_error_erg": 0.0,
        "sum_balance_minus_error_erg": 0.0,
        "max_abs_balance_erg": 0.0,
        "max_balance_k": -1,
        "max_balance_r_cm": np.nan,
        "max_abs_balance_minus_error_erg": 0.0,
        "max_balance_minus_error_k": -1,
        "max_balance_minus_error_r_cm": np.nan,
        "max_relative_balance": 0.0,
        "max_relative_balance_minus_error": 0.0,
        "max_abs_dU_near_shock_erg": 0.0,
        "max_dU_near_shock_r_cm": np.nan,
        "center_rho": np.nan,
        "center_P": np.nan,
        "center_energy": np.nan,
        "center_u": np.nan,
        "center_r_cm": np.inf,
        "min_postshock_rho": np.inf,
        "min_postshock_r_cm": np.nan,
    }


def update_aggregate(record: dict[str, float | int], values: list[float]) -> None:
    k = int(values[1])
    radius = values[4]
    rho = values[6]
    balance = values[13]
    equation_error = values[14]
    balance_minus_error = values[15]
    dU = values[10]
    local_energy_scale = abs(dU) + abs(values[11]) + abs(values[12])

    record["rows_read"] += 1
    record["sum_equation_error_erg"] += equation_error
    record["sum_balance_minus_error_erg"] += balance_minus_error

    if abs(balance) > record["max_abs_balance_erg"]:
        record["max_abs_balance_erg"] = abs(balance)
        record["max_balance_k"] = k
        record["max_balance_r_cm"] = radius
    if abs(balance_minus_error) > record["max_abs_balance_minus_error_erg"]:
        record["max_abs_balance_minus_error_erg"] = abs(balance_minus_error)
        record["max_balance_minus_error_k"] = k
        record["max_balance_minus_error_r_cm"] = radius
    if local_energy_scale > 0.0:
        record["max_relative_balance"] = max(
            record["max_relative_balance"], abs(balance) / local_energy_scale
        )
        record["max_relative_balance_minus_error"] = max(
            record["max_relative_balance_minus_error"],
            abs(balance_minus_error) / local_energy_scale,
        )
    if radius < record["center_r_cm"]:
        record["center_r_cm"] = radius
        record["center_rho"] = rho
        record["center_P"] = values[8]
        record["center_energy"] = values[9]
        record["center_u"] = values[7]

    shock_radius = float(record["age_s"]) / 3.0
    if radius <= 1.5 * shock_radius:
        if abs(dU) > record["max_abs_dU_near_shock_erg"]:
            record["max_abs_dU_near_shock_erg"] = abs(dU)
            record["max_dU_near_shock_r_cm"] = radius
    if shock_radius > 0.0 and radius <= 0.8 * shock_radius:
        if rho < record["min_postshock_rho"]:
            record["min_postshock_rho"] = rho
            record["min_postshock_r_cm"] = radius


def read_balance_file(
    path: Path, selected_models: set[int]
) -> tuple[list[dict[str, float | int]], dict[int, np.ndarray]]:
    summaries: dict[int, dict[str, float | int]] = {}
    profiles: dict[int, list[list[float]]] = {model: [] for model in selected_models}
    current_model: int | None = None

    with path.open("r", encoding="ascii", errors="strict") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("# model "):
                header = parse_model_header(line)
                if header is not None:
                    current_model = int(header["model"])
                    summaries[current_model] = new_aggregate(header)
                    if current_model in profiles:
                        profiles[current_model] = []
                continue
            if not line or line.startswith("#"):
                continue

            fields = line.split()
            if len(fields) != len(DATA_COLUMNS):
                # A concurrently running model can leave one incomplete final line.
                continue
            try:
                values = [float(field) for field in fields]
            except ValueError as error:
                raise ValueError(f"invalid data on line {line_number}") from error

            model = int(values[0])
            if current_model != model or model not in summaries:
                raise ValueError(
                    f"data for model {model} has no matching header on line {line_number}"
                )
            update_aggregate(summaries[model], values)
            if model in profiles:
                profiles[model].append(values)

    if not summaries:
        raise ValueError(f"no model data found in {path}")

    summary_rows = [
        summaries[model]
        for model in sorted(summaries)
        if summaries[model]["rows_read"] == summaries[model]["nz"]
    ]
    if not summary_rows:
        raise ValueError(f"no complete model data found in {path}")
    complete_models = {int(row["model"]) for row in summary_rows}
    profile_arrays = {
        model: np.asarray(rows, dtype=float)
        for model, rows in profiles.items()
        if rows and model in complete_models
    }
    for row in summary_rows:
        if not np.isfinite(row["min_postshock_rho"]):
            row["min_postshock_rho"] = np.nan
    return summary_rows, profile_arrays


def write_summary_csv(path: Path, rows: list[dict[str, float | int]]) -> None:
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def symlog_threshold(*arrays: np.ndarray) -> float:
    values = np.concatenate([np.abs(array[np.isfinite(array)]) for array in arrays])
    values = values[values > 0.0]
    if values.size == 0:
        return 1.0
    return max(float(np.nanmax(values)) * 1.0e-10, np.finfo(float).tiny)


def save_balance_history(output: Path, rows: list[dict[str, float | int]]) -> None:
    age = np.asarray([row["age_s"] for row in rows], dtype=float)
    dU = np.asarray([row["total_dU_erg"] for row in rows], dtype=float)
    dK = np.asarray([row["total_dK_erg"] for row in rows], dtype=float)
    dW = np.asarray([row["total_dW_erg"] for row in rows], dtype=float)
    balance = np.asarray([row["total_balance_erg"] for row in rows], dtype=float)
    equation_error = np.asarray(
        [row["sum_equation_error_erg"] for row in rows], dtype=float
    )
    closure_difference = np.asarray(
        [row["sum_balance_minus_error_erg"] for row in rows], dtype=float
    )
    global_energy_scale = np.abs(dU) + np.abs(dK) + np.abs(dW)
    relative_balance = np.zeros_like(balance)
    relative_equation_error = np.zeros_like(equation_error)
    relative_closure_difference = np.zeros_like(closure_difference)
    np.divide(
        np.abs(balance),
        global_energy_scale,
        out=relative_balance,
        where=global_energy_scale > 0.0,
    )
    np.divide(
        np.abs(equation_error),
        global_energy_scale,
        out=relative_equation_error,
        where=global_energy_scale > 0.0,
    )
    np.divide(
        np.abs(closure_difference),
        global_energy_scale,
        out=relative_closure_difference,
        where=global_energy_scale > 0.0,
    )
    max_relative_balance = np.asarray(
        [row["max_relative_balance"] for row in rows], dtype=float
    )
    max_relative_difference = np.asarray(
        [row["max_relative_balance_minus_error"] for row in rows], dtype=float
    )

    figure, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
    floor = np.finfo(float).tiny
    axes[0].semilogy(age, np.maximum(np.abs(dU), floor), label=r"$|\sum\Delta U|$")
    axes[0].semilogy(age, np.maximum(np.abs(dK), floor), label=r"$|\sum\Delta K|$")
    axes[0].semilogy(age, np.maximum(np.abs(dW), floor), label=r"$|\sum\Delta W|$")
    axes[0].set_ylabel("Term magnitude (erg per step)")
    axes[0].legend(ncol=3)

    axes[1].semilogy(
        age, np.maximum(relative_balance, floor), label=r"$|\sum B_k|/E_{\rm step}$"
    )
    axes[1].semilogy(
        age,
        np.maximum(relative_equation_error, floor),
        label=r"$|\sum$ equation error$|/E_{\rm step}$",
    )
    axes[1].semilogy(
        age,
        np.maximum(relative_closure_difference, floor),
        label="closure difference",
    )
    axes[1].set_ylabel("Relative global closure")
    axes[1].legend(ncol=3)

    axes[2].semilogy(
        age,
        np.maximum(max_relative_balance, floor),
        label=r"max $|B_k|/E_{{\rm step},k}$",
    )
    axes[2].semilogy(
        age,
        np.maximum(max_relative_difference, floor),
        label=r"max $|B_k-\mathrm{error}_k|/E_{{\rm step},k}$",
    )
    axes[2].set_ylabel("Relative local closure")
    axes[2].set_xlabel("Time (s)")
    axes[2].legend()

    figure.suptitle("Noh discrete energy balance")
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def save_wall_heating_history(output: Path, rows: list[dict[str, float | int]]) -> None:
    age = np.asarray([row["age_s"] for row in rows], dtype=float)
    center_rho = np.asarray([row["center_rho"] for row in rows], dtype=float)
    min_rho = np.asarray([row["min_postshock_rho"] for row in rows], dtype=float)
    center_pressure = np.asarray([row["center_P"] for row in rows], dtype=float)
    center_energy = np.asarray([row["center_energy"] for row in rows], dtype=float)
    min_rho_radius = np.asarray(
        [row["min_postshock_r_cm"] for row in rows], dtype=float
    )
    max_dU_radius = np.asarray(
        [row["max_dU_near_shock_r_cm"] for row in rows], dtype=float
    )
    shock_radius = age / 3.0

    figure, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
    axes[0].plot(age, center_rho / NOH_DENSITY, label="innermost cell")
    axes[0].plot(age, min_rho / NOH_DENSITY, label=r"minimum for $r<0.8r_s$")
    axes[0].axhline(1.0, color="black", linestyle="--", label="analytic")
    axes[0].set_yscale("log")
    axes[0].set_ylabel(r"$\rho/\rho_{\rm Noh}$")
    axes[0].legend()

    axes[1].plot(age, center_pressure / NOH_PRESSURE, label="central pressure")
    axes[1].plot(age, center_energy / NOH_ENERGY, label="central specific energy")
    axes[1].axhline(1.0, color="black", linestyle="--", label="analytic")
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Central value / analytic value")
    axes[1].legend(ncol=2)

    normalized_min_rho_radius = np.full_like(age, np.nan)
    normalized_max_dU_radius = np.full_like(age, np.nan)
    np.divide(
        min_rho_radius,
        shock_radius,
        out=normalized_min_rho_radius,
        where=shock_radius > 0.0,
    )
    np.divide(
        max_dU_radius,
        shock_radius,
        out=normalized_max_dU_radius,
        where=shock_radius > 0.0,
    )
    axes[2].plot(age, normalized_min_rho_radius, label=r"radius of minimum $\rho$")
    axes[2].plot(
        age,
        normalized_max_dU_radius,
        label=r"radius of max $|\Delta U|$ near shock",
    )
    axes[2].axhline(1.0, color="black", linestyle="--", label="shock")
    axes[2].set_ylabel(r"$r/r_s$")
    axes[2].set_xlabel("Time (s)")
    axes[2].legend()

    figure.suptitle("Noh wall-heating diagnostics")
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def save_wall_heating_profiles(output: Path, profiles: dict[int, np.ndarray]) -> None:
    usable = {
        model: data
        for model, data in profiles.items()
        if np.any(data[:, 4] <= 0.8 * data[0, 2] / 3.0)
    }
    if not usable:
        return

    figure, axes = plt.subplots(4, 1, figsize=(10, 13), sharex=True)
    colors = plt.get_cmap("viridis")(np.linspace(0.0, 1.0, len(usable)))
    for color, model in zip(colors, sorted(usable)):
        data = usable[model]
        shock_radius = data[0, 2] / 3.0
        view = data[data[:, 4] <= 1.2 * shock_radius]
        view = view[np.argsort(view[:, 4])]
        radius = view[:, 4] / shock_radius
        entropy_proxy = view[:, 8] / view[:, 6] ** GAMMA
        label = f"{model}: t={data[0, 2]:.3e} s"

        axes[0].plot(radius, view[:, 6] / NOH_DENSITY, color=color, label=label)
        axes[1].plot(radius, view[:, 8] / NOH_PRESSURE, color=color)
        axes[2].plot(radius, view[:, 9] / NOH_ENERGY, color=color)
        axes[3].plot(radius, entropy_proxy / NOH_ENTROPY_PROXY, color=color)

    labels = (
        r"$\rho/\rho_{\rm Noh}$",
        r"$P/P_{\rm Noh}$",
        r"$e/e_{\rm Noh}$",
        r"$(P/\rho^\gamma)/(P/\rho^\gamma)_{\rm Noh}$",
    )
    for axis, label in zip(axes, labels):
        axis.axhline(1.0, color="black", linestyle="--")
        axis.axvline(1.0, color="gray", linestyle=":")
        axis.set_yscale("log")
        axis.set_ylabel(label)
        axis.set_xlim(0.0, 1.2)
    axes[0].legend(ncol=2)
    axes[-1].set_xlabel(r"$r/r_s$")

    figure.suptitle("Noh wall-heating structure")
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def save_profile_plot(
    output: Path, model: int, data: np.ndarray, radius_max: float | None
) -> None:
    age = data[0, 2]
    shock_radius = age / 3.0
    xmax = radius_max
    if xmax is None:
        xmax = min(0.2, max(0.02, 2.5 * shock_radius))
    mask = data[:, 4] <= xmax
    view = data[mask]
    if view.size == 0:
        return

    radius = view[:, 4]
    dm = view[:, 5]
    specific = np.zeros((view.shape[0], 6))
    np.divide(view[:, 10:16], dm[:, None], out=specific, where=dm[:, None] > 0.0)

    figure, axes = plt.subplots(4, 1, figsize=(10, 13), sharex=True)
    axes[0].plot(radius, view[:, 6], label="density")
    axes[0].axhline(64.0, color="black", linestyle="--", label="analytic 64")
    axes[0].set_ylabel(r"$\rho$ (g cm$^{-3}$)")
    axes[0].legend()

    axes[1].plot(radius, specific[:, 0], label=r"$\Delta U/dm$")
    axes[1].plot(radius, specific[:, 1], label=r"$\Delta K/dm$")
    axes[1].plot(radius, specific[:, 2], label=r"$\Delta W/dm$")
    axes[1].set_yscale(
        "symlog",
        linthresh=symlog_threshold(specific[:, 0], specific[:, 1], specific[:, 2]),
    )
    axes[1].set_ylabel("Specific energy per step (erg/g)")
    axes[1].legend(ncol=3)

    axes[2].plot(radius, specific[:, 3], label=r"$B_k/dm$")
    axes[2].plot(radius, specific[:, 4], label="equation error/dm")
    axes[2].plot(radius, specific[:, 5], label="difference/dm")
    axes[2].set_yscale(
        "symlog",
        linthresh=symlog_threshold(specific[:, 3], specific[:, 4], specific[:, 5]),
    )
    axes[2].set_ylabel("Specific closure error (erg/g)")
    axes[2].legend(ncol=3)

    axes[3].plot(radius, view[:, 16], label="relative balance")
    axes[3].plot(radius, view[:, 7], label="velocity")
    axes[3].set_yscale("symlog", linthresh=1.0e-10)
    axes[3].set_ylabel("Relative balance or velocity")
    axes[3].set_xlabel("Radius (cm)")
    axes[3].legend()

    for axis in axes:
        axis.axvline(shock_radius, color="gray", linestyle=":")
        axis.set_xlim(0.0, xmax)

    figure.suptitle(f"Noh energy balance: model {model}, t={age:.6e} s")
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def main() -> None:
    args = parse_args()
    if not args.input.is_file():
        raise FileNotFoundError(args.input)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows, profiles = read_balance_file(args.input, set(args.models))
    write_summary_csv(args.output_dir / "energy_balance_summary.csv", rows)
    save_balance_history(args.output_dir / "energy_balance_history.png", rows)
    save_wall_heating_history(args.output_dir / "wall_heating_history.png", rows)
    save_wall_heating_profiles(args.output_dir / "wall_heating_profiles.png", profiles)
    for model in sorted(profiles):
        save_profile_plot(
            args.output_dir / f"energy_balance_model_{model:06d}.png",
            model,
            profiles[model],
            args.radius_max,
        )

    available = ", ".join(str(model) for model in sorted(profiles))
    print(f"Read {len(rows)} accepted models from {args.input}")
    print(f"Wrote plots and summary to {args.output_dir}")
    print(f"Detailed profiles: {available or 'none of the requested models'}")


if __name__ == "__main__":
    main()
