#!/usr/bin/env python
"""Plot the selected spectrum and one radial LNA eigenfunction."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


CASE = Path(__file__).resolve().parents[1]


def read_table(path):
    names = None
    with path.open() as stream:
        for line in stream:
            if line.startswith(("# mode control_mode_index ", "# k q ", "# k logT ")):
                names = line[1:].split()
                break
        if names is None:
            raise ValueError(f"Column header missing: {path}")
        rows = [line for line in stream if line.strip() and not line.startswith("#")]
    if not rows:
        raise ValueError(f"No selected data in {path}; inspect the LNA run output.")
    data = np.loadtxt(rows, ndmin=2)
    if data.shape[1] != len(names):
        raise ValueError(f"Column count does not match the header: {path}")
    return {name: data[:, i] for i, name in enumerate(names)}


def plot_spectrum(modes, output, selected_index):
    fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True, layout="constrained")
    period = modes["period_days"]
    growth = modes["logKE_per_cycle"]
    colors = np.where(growth > 0, "tab:red", "tab:blue")
    axes[0].scatter(period, growth, c=colors)
    axes[0].axhline(0, color="0.4", linewidth=0.8)
    axes[0].set_ylabel("Log kinetic energy change per cycle")
    axes[0].set_title("Selected radial LNA modes")
    for x, y, index in zip(period, growth, modes["control_mode_index"]):
        if len(period) > 20 and int(index) != selected_index:
            continue
        axes[0].annotate(
            str(int(index)),
            (x, y),
            xytext=(3, 4),
            textcoords="offset points",
            fontsize=8,
        )
    residual = np.maximum(modes["max_eigenvector_residual"], np.finfo(float).tiny)
    axes[1].scatter(period, residual, c=colors)
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Eigenvector residual")
    axes[1].set_xlabel("Period (days)")
    for ax in axes:
        ax.set_xscale("log")
        ax.grid(alpha=0.2)
    fig.savefig(output / "spectrum.png", dpi=180)
    plt.close(fig)


def plot_mode(directory, modes, index, output):
    selected = np.flatnonzero(modes["control_mode_index"].astype(int) == index)
    if len(selected) != 1:
        raise ValueError(f"Mode {index} is not in the selected table.")
    row = selected[0]
    file_mode = int(modes["mode"][row])
    data = read_table(directory / f"star_LNA_eigenfunction_{file_mode}.data")
    radius = data["r"] / np.max(data["r"])
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
    labels = [r"$\delta\ln r$", r"$\delta\ln T$", r"$\delta L/\max|\delta L|$"]
    for ax, variable, label in zip(axes.flat, ("lnR", "lnT", "L"), labels):
        value = data[f"re_{variable}"] + 1j * data[f"im_{variable}"]
        scale = np.max(np.abs(value)) if variable == "L" else 1.0
        if scale == 0:
            scale = 1.0
        ax.plot(radius, value.real / scale, label="Real")
        ax.plot(radius, value.imag / scale, label="Imaginary")
        ax.set_ylabel(label)
        ax.legend(fontsize=8)
    work_file = directory / f"star_LNA_work_{file_mode}.data"
    ax = axes[1, 1]
    if work_file.exists():
        work = read_table(work_file)
        for name, label in (
            ("pressure", "Pressure"),
            ("turb_pressure", "Turbulent pressure"),
            ("eddy_visc", "Eddy viscosity"),
            ("rad_lum", "Radiative"),
            ("conv_lum", "Convective"),
            ("turb_lum", "Turbulent transport"),
        ):
            ax.plot(work["r_div_R"], work[f"c_{name}_work"], label=label)
        ax.legend(fontsize=7)
        ax.set_ylabel("Cumulative diagnostic work / kinetic energy")
    else:
        ax.text(
            0.5, 0.5, "Work output not available", ha="center", transform=ax.transAxes
        )
    for ax in axes.flat:
        ax.set_xlabel(r"$r/R$")
        ax.grid(alpha=0.2)
    fig.suptitle(
        f"Mode {index}: P = {modes['period_days'][row]:.6g} days, "
        f"logKE/cycle = {modes['logKE_per_cycle'][row]:.5g}"
    )
    fig.savefig(output / f"mode_{index}.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lna", type=Path, default=CASE / "LNA")
    parser.add_argument("--output", type=Path, default=CASE / "plots")
    parser.add_argument(
        "--mode", type=int, default=0, help="Zero-based terminal-table mode index"
    )
    args = parser.parse_args()
    try:
        modes = read_table(args.lna / "star_LNA_period_growth.data")
        args.output.mkdir(parents=True, exist_ok=True)
        plot_spectrum(modes, args.output, args.mode)
        plot_mode(args.lna, modes, args.mode, args.output)
    except (OSError, ValueError) as error:
        parser.exit(1, f"{error}\n")
    print(f"Plots written to {args.output}")


if __name__ == "__main__":
    main()
