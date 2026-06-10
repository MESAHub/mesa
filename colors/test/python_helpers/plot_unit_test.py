#!/usr/bin/env python3
"""
plot_unit_test_pgstar_style.py
==============================
PGSTAR-styled diagnostic plots for the MESA colors module unit test.

This keeps the enhanced unit-test diagnostics but restyles the figure so it
looks much closer to a PGSTAR panel layout:
  - black background
  - white axes, ticks, and labels
  - cyan / yellow / orange accent colours
  - inward ticks and boxed panels
  - compact in-panel diagnostic text

It reads colors/test/test_output and produces a four-panel figure:
  1. Colour indices vs [M/H]
  2. Colour indices vs log g
  3. Colour indices vs Teff
  4. Solar SED sanity panel with Johnson transmissions, a scaled blackbody,
     Wien-law peak marker, and sampled SED peak marker.

Run from colors/test/:
    python3 python_helpers/plot_unit_test_pgstar_style.py

Requires: matplotlib, numpy
"""

from __future__ import annotations

import os
import sys
from typing import Dict, Iterable, List, Tuple

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

# -----------------------------------------------------------------------------
# paths
# -----------------------------------------------------------------------------

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_TEST_DIR = os.path.dirname(SCRIPT_DIR)


def resolve_test_output() -> tuple[str, str]:
    """Return (test_output_path, output_pdf_path) using a few sensible fallbacks."""
    candidates = []

    env_path = os.environ.get("TEST_OUTPUT_PATH", "").strip()
    if env_path:
        candidates.append(env_path)

    candidates.extend(
        [
            os.path.join(DEFAULT_TEST_DIR, "test_output"),
            os.path.join(SCRIPT_DIR, "test_output"),
            os.path.join(os.getcwd(), "test_output"),
        ]
    )

    for cand in candidates:
        if cand and os.path.exists(cand):
            out_dir = os.path.dirname(cand)
            return cand, os.path.join(out_dir, "unit_test_diagnostic_pgstar.pdf")

    # default, even if not found yet, so downstream error messages are clear
    fallback = os.path.join(DEFAULT_TEST_DIR, "test_output")
    return fallback, os.path.join(DEFAULT_TEST_DIR, "unit_test_diagnostic_pgstar.pdf")


TEST_OUTPUT, OUT_FILE = resolve_test_output()

MESA_DIR = os.environ.get("MESA_DIR", "")
FILTER_DIR = (
    os.path.join(MESA_DIR, "data", "colors_data", "filters", "Generic", "Johnson")
    if MESA_DIR
    else ""
)

# -----------------------------------------------------------------------------
# PGSTAR-like style constants
# -----------------------------------------------------------------------------

PG_BG = "#000000"
PG_FG = "#f2f2f2"
PG_GRID = "#00d7d7"
PG_CYAN = "#00c8d7"
PG_SKY = "#73d2ff"
PG_BLUE = "#4c6fff"
PG_GREEN = "#31d843"
PG_YELLOW = "#ffd400"
PG_ORANGE = "#ff9d00"
PG_RED = "#ff4040"
PG_MAGENTA = "#ff62d5"
PG_GREY = "#8a8a8a"

COLOR_INDEX_STYLES = {
    "U-B": {"color": PG_YELLOW, "marker": "o"},
    "B-V": {"color": PG_CYAN, "marker": "s"},
    "V-I": {"color": PG_ORANGE, "marker": "^"},
    "V-J": {"color": PG_MAGENTA, "marker": "D"},
}

FILTER_COLORS = {
    "U": PG_MAGENTA,
    "B": PG_BLUE,
    "V": PG_GREEN,
    "R": PG_ORANGE,
    "I": PG_RED,
    "J": PG_CYAN,
    "M": PG_GREY,
}

GROUP1_TEMPS = {
    "solar": 5778.0,
    "hot_ms": 15000.0,
    "cool_giant": 4000.0,
}

# -----------------------------------------------------------------------------
# parsing
# -----------------------------------------------------------------------------


def parse_test_output(path: str) -> Dict[str, object]:
    """Parse the colors test output file into structured blocks."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"test output not found: {path}")

    data: Dict[str, object] = {
        "group1": {},
        "group2a": {},
        "group2b": {},
        "group2c": {},
        "sed_wav": [],
        "sed_flux": [],
    }

    section = None
    current_key = None
    current_block: Dict[str, float] = {}

    def flush(dest: Dict, key, block: Dict[str, float]) -> None:
        if key is not None and block:
            dest[key] = dict(block)

    with open(path, encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue

            if line.startswith("# group1"):
                section = "group1"
                current_key = None
                current_block = {}
                continue
            if line.startswith("# group2a"):
                section = "group2a"
                current_key = None
                current_block = {}
                continue
            if line.startswith("# group2b"):
                section = "group2b"
                current_key = None
                current_block = {}
                continue
            if line.startswith("# group2c"):
                section = "group2c"
                current_key = None
                current_block = {}
                continue
            if line.startswith("# SED sample"):
                if section in ("group1", "group2a", "group2b", "group2c"):
                    flush(data[section], current_key, current_block)
                section = "sed"
                continue

            if line.startswith("# case:") and section == "group1":
                flush(data["group1"], current_key, current_block)
                current_key = line.split(":", 1)[-1].strip()
                current_block = {}
                continue
            if line.startswith("# MH=") and section == "group2a":
                flush(data["group2a"], current_key, current_block)
                current_key = float(line.split("=", 1)[-1].strip())
                current_block = {}
                continue
            if line.startswith("# logg=") and section == "group2b":
                flush(data["group2b"], current_key, current_block)
                current_key = float(line.split("=", 1)[-1].strip())
                current_block = {}
                continue
            if line.startswith("# Teff=") and section == "group2c":
                flush(data["group2c"], current_key, current_block)
                current_key = float(line.split("=", 1)[-1].strip())
                current_block = {}
                continue

            if line.startswith("#"):
                continue

            if section == "sed":
                parts = line.split()
                if len(parts) == 2:
                    try:
                        data["sed_wav"].append(float(parts[0]))
                        data["sed_flux"].append(float(parts[1]))
                    except ValueError:
                        pass
                continue

            if (
                section in ("group1", "group2a", "group2b", "group2c")
                and current_key is not None
            ):
                parts = line.split()
                if len(parts) == 2:
                    current_block[parts[0]] = float(parts[1])

    if section in ("group1", "group2a", "group2b", "group2c"):
        flush(data[section], current_key, current_block)

    data["sed_wav"] = np.array(data["sed_wav"], dtype=float)
    data["sed_flux"] = np.array(data["sed_flux"], dtype=float)
    return data


# -----------------------------------------------------------------------------
# filters
# -----------------------------------------------------------------------------


def load_filters(filter_dir: str) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """Load Johnson filter curves if available. Return an empty dict otherwise."""
    filters: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    if not filter_dir or not os.path.isdir(filter_dir):
        return filters

    for fname in os.listdir(filter_dir):
        if not fname.endswith(".dat"):
            continue
        name = fname[:-4]
        wav: List[float] = []
        trans: List[float] = []
        with open(os.path.join(filter_dir, fname), encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    wav.append(float(parts[0]))
                    trans.append(float(parts[1]))
        if wav:
            filters[name] = (np.array(wav, dtype=float), np.array(trans, dtype=float))

    if not filters:
        return filters

    eff_wav = {
        name: np.average(w, weights=np.maximum(t, 0.0))
        for name, (w, t) in filters.items()
    }
    return dict(sorted(filters.items(), key=lambda item: eff_wav[item[0]]))


# -----------------------------------------------------------------------------
# physics helpers
# -----------------------------------------------------------------------------


def mag_color(block: Dict[str, float], blue: str, red: str) -> float:
    if blue not in block or red not in block:
        return np.nan
    return block[blue] - block[red]


def planck_lambda_aa(wav_aa: np.ndarray, temperature_k: float) -> np.ndarray:
    """Planck function B_lambda in arbitrary cgs-consistent units."""
    h = 6.62607015e-27  # erg s
    c = 2.99792458e10  # cm s^-1
    k = 1.380649e-16  # erg K^-1

    lam_cm = wav_aa * 1.0e-8
    expo = (h * c) / (lam_cm * k * temperature_k)
    expo = np.clip(expo, 1.0e-12, 700.0)
    return (2.0 * h * c**2) / (lam_cm**5) / np.expm1(expo)


def wien_peak_aa(temperature_k: float) -> float:
    return 2.897771955e7 / temperature_k


# -----------------------------------------------------------------------------
# styling helpers
# -----------------------------------------------------------------------------


def set_pgstar_rcparams() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": PG_BG,
            "savefig.facecolor": PG_BG,
            "axes.facecolor": PG_BG,
            "axes.edgecolor": PG_FG,
            "axes.labelcolor": PG_FG,
            "axes.titlecolor": PG_FG,
            "xtick.color": PG_FG,
            "ytick.color": PG_FG,
            "text.color": PG_FG,
            "legend.facecolor": PG_BG,
            "legend.edgecolor": PG_FG,
            "font.size": 11,
            "font.family": "DejaVu Sans",
            "mathtext.default": "regular",
        }
    )


def style_axis(ax, title: str | None = None) -> None:
    ax.set_facecolor(PG_BG)
    for spine in ax.spines.values():
        spine.set_color(PG_FG)
        spine.set_linewidth(1.0)
    ax.tick_params(
        axis="both",
        which="major",
        colors=PG_FG,
        direction="in",
        top=True,
        right=True,
        length=6,
        width=1.0,
        labelsize=9,
    )
    ax.tick_params(
        axis="both",
        which="minor",
        colors=PG_FG,
        direction="in",
        top=True,
        right=True,
        length=3,
        width=0.8,
    )
    ax.minorticks_on()
    ax.grid(
        True, which="major", color=PG_GRID, alpha=0.18, linestyle="--", linewidth=0.6
    )
    if title is not None:
        ax.set_title(title, fontsize=10, color=PG_FG, pad=7)


def style_secondary_axis(ax) -> None:
    for spine in ax.spines.values():
        spine.set_color(PG_FG)
        spine.set_linewidth(1.0)
    ax.tick_params(
        axis="both",
        which="major",
        colors=PG_FG,
        direction="in",
        top=True,
        right=True,
        length=6,
        width=1.0,
        labelsize=9,
    )
    ax.tick_params(
        axis="both", which="minor", colors=PG_FG, direction="in", top=True, right=True
    )
    ax.set_facecolor("none")


def add_panel_tag(ax, tag: str) -> None:
    ax.text(
        0.02,
        0.98,
        tag,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        color=PG_YELLOW,
    )


# -----------------------------------------------------------------------------
# plotting helpers
# -----------------------------------------------------------------------------


def plot_color_sweep(
    ax,
    sweep_data: Dict[float, Dict[str, float]],
    x_vals: Iterable[float],
    x_label: str,
    panel_tag: str,
    bol_on_right: bool = True,
) -> None:
    x_vals = list(x_vals)
    color_defs = [("U", "B"), ("B", "V"), ("V", "I"), ("V", "J")]

    style_axis(ax)
    add_panel_tag(ax, panel_tag)

    for blue, red in color_defs:
        label = f"{blue}-{red}"
        y = [mag_color(sweep_data[x], blue, red) for x in x_vals]
        style = COLOR_INDEX_STYLES[label]
        ax.plot(
            x_vals,
            y,
            color=style["color"],
            marker=style["marker"],
            linestyle="-",
            linewidth=1.6,
            markersize=4.5,
            label=label,
        )

    ax.set_xlabel(x_label, fontsize=10)
    ax.set_ylabel("colour index (mag)", fontsize=10)

    if bol_on_right:
        ax2 = ax.twinx()
        style_secondary_axis(ax2)
        ref = sweep_data[x_vals[0]].get("Mag_bol", np.nan)
        y_bol = [sweep_data[x].get("Mag_bol", np.nan) - ref for x in x_vals]
        ax2.plot(
            x_vals,
            y_bol,
            color=PG_FG,
            linestyle="--",
            linewidth=1.1,
            marker="x",
            markersize=4,
            label=r"$\Delta$Mag$_{bol}$",
        )
        ax2.set_ylabel(r"$\Delta$Mag$_{bol}$ (mag)", fontsize=9, color=PG_FG)


def make_group1_summary(group1: Dict[str, Dict[str, float]]) -> str:
    lines = ["anchor cases"]
    for case in ("cool_giant", "solar", "hot_ms"):
        if case not in group1:
            continue
        block = group1[case]
        temp = GROUP1_TEMPS.get(case)
        lines.append(
            f"{case:>10s}  T={temp:5.0f}  B-V={mag_color(block, 'B', 'V'):+5.2f}"
        )
    return "\n".join(lines)


# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------


def main() -> None:
    set_pgstar_rcparams()

    data = parse_test_output(TEST_OUTPUT)
    filters = load_filters(FILTER_DIR)

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.0))
    fig.patch.set_facecolor(PG_BG)

    # ------------------------------------------------------------------
    # panel 1
    # ------------------------------------------------------------------
    ax = axes[0, 0]
    g2a = data["group2a"]
    x_mh = sorted(g2a.keys())
    plot_color_sweep(ax, g2a, x_mh, "[M/H]", panel_tag="colour sweep :: metallicity")
    ax.set_title(r"fixed $T_{eff}=5778$ K, $\log g=4.44$", fontsize=10, color=PG_FG)
    ax.set_xticks(x_mh)

    # ------------------------------------------------------------------
    # panel 2
    # ------------------------------------------------------------------
    ax = axes[0, 1]
    g2b = data["group2b"]
    x_logg = sorted(g2b.keys())
    plot_color_sweep(ax, g2b, x_logg, r"$\log g$", panel_tag="colour sweep :: gravity")
    ax.set_title(r"fixed $T_{eff}=5778$ K, [M/H] = 0.0", fontsize=10, color=PG_FG)
    ax.set_xticks(x_logg)

    # shared legend for top row / left-bottom diagnostics
    sweep_handles = [
        Line2D(
            [0],
            [0],
            color=style["color"],
            marker=style["marker"],
            linestyle="-",
            linewidth=1.6,
            markersize=5,
            label=label,
        )
        for label, style in COLOR_INDEX_STYLES.items()
    ]
    sweep_handles.append(
        Line2D(
            [0],
            [0],
            color=PG_FG,
            marker="x",
            linestyle="--",
            linewidth=1.1,
            markersize=5,
            label=r"$\Delta$Mag$_{bol}$",
        )
    )
    legend = ax.legend(
        handles=sweep_handles,
        fontsize=8,
        loc="lower left",
        framealpha=1.0,
        ncol=1,
        title="tracks",
        title_fontsize=8,
    )
    legend.get_frame().set_facecolor(PG_BG)
    legend.get_frame().set_edgecolor(PG_FG)
    for text in legend.get_texts():
        text.set_color(PG_FG)
    legend.get_title().set_color(PG_YELLOW)

    # ------------------------------------------------------------------
    # panel 3
    # ------------------------------------------------------------------
    ax = axes[1, 0]
    g2c = data["group2c"]
    x_teff = sorted(g2c.keys())
    plot_color_sweep(
        ax, g2c, x_teff, r"$T_{eff}$ (K)", panel_tag="colour sweep :: temperature"
    )
    ax.set_title(r"fixed $\log g=4.0$, [M/H] = 0.0", fontsize=10, color=PG_FG)
    ax.set_xticks(x_teff)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

    # small in-panel text summary like PGSTAR annotations
    summary_text = make_group1_summary(data["group1"])
    ax.text(
        0.03,
        0.05,
        summary_text,
        transform=ax.transAxes,
        fontsize=7.8,
        color=PG_YELLOW,
        ha="left",
        va="bottom",
        family="DejaVu Sans Mono",
        bbox={
            "boxstyle": "square,pad=0.25",
            "facecolor": PG_BG,
            "edgecolor": PG_FG,
            "alpha": 0.95,
        },
    )

    # ------------------------------------------------------------------
    # panel 4
    # ------------------------------------------------------------------
    ax = axes[1, 1]
    style_axis(ax, title="solar SED sanity check")
    add_panel_tag(ax, "sampled atmosphere  /  blackbody  /  filters")
    ax2 = ax.twinx()
    style_secondary_axis(ax2)

    wav = data["sed_wav"]
    flux = data["sed_flux"]
    if wav.size == 0 or flux.size == 0:
        raise ValueError("No SED sample found in test_output.")

    mask = flux > 0.0
    wav_pos = wav[mask]
    flux_pos = flux[mask]

    sed_line = ax.semilogy(
        wav_pos,
        flux_pos,
        color=PG_CYAN,
        linewidth=2.0,
        label="solar SED",
        zorder=4,
    )[0]

    t_sun = 5778.0
    bb = planck_lambda_aa(wav_pos, t_sun)
    bb *= flux_pos.max() / bb.max()
    bb_line = ax.semilogy(
        wav_pos,
        bb,
        color=PG_YELLOW,
        linestyle="--",
        linewidth=1.5,
        label="scaled blackbody",
        zorder=3,
    )[0]

    lambda_wien = wien_peak_aa(t_sun)
    sampled_peak_idx = int(np.argmax(flux_pos))
    sampled_peak_wav = wav_pos[sampled_peak_idx]
    sampled_peak_flux = flux_pos[sampled_peak_idx]

    wien_line = ax.axvline(
        lambda_wien,
        color=PG_ORANGE,
        linestyle=":",
        linewidth=1.2,
        label="Wien peak",
        zorder=2,
    )
    peak_point = ax.scatter(
        [sampled_peak_wav],
        [sampled_peak_flux],
        color=PG_RED,
        marker="o",
        s=24,
        label="sampled peak",
        zorder=5,
    )

    filter_handles: List[Line2D] = []
    if filters:
        for fname, (fwav, ftrans) in filters.items():
            peak = np.max(ftrans)
            if peak <= 0.0:
                continue
            norm_t = ftrans / peak
            color = FILTER_COLORS.get(fname, PG_GREY)
            ax2.fill_between(fwav, norm_t, alpha=0.13, color=color)
            handle = ax2.plot(
                fwav,
                norm_t,
                color=color,
                linewidth=1.0,
                alpha=0.95,
                label=fname,
            )[0]
            filter_handles.append(handle)
    else:
        ax.text(
            0.03,
            0.92,
            "Johnson filter files not found; plotting SED-only diagnostic",
            transform=ax.transAxes,
            fontsize=7.5,
            color=PG_ORANGE,
            ha="left",
            va="top",
            family="DejaVu Sans Mono",
        )

    ax.set_xlabel("wavelength (A)", fontsize=10)
    ax.set_ylabel(r"flux (erg s$^{-1}$ cm$^{-2}$ A$^{-1}$)", fontsize=10)
    ax2.set_ylabel("normalised transmission", fontsize=9, color=PG_FG)
    ax2.set_ylim(0.0, 1.35)
    ax.set_xlim(wav_pos.min() * 0.9, min(wav_pos.max(), 1.1e4))

    # white/amber annotation block instead of matplotlib arrows all over the place
    info = (
        f"T_solar = {t_sun:7.1f} K\n"
        f"lambda_Wien = {lambda_wien:7.0f} A\n"
        f"lambda_peak(sampled) = {sampled_peak_wav:7.0f} A\n"
        f"delta = {sampled_peak_wav - lambda_wien:+7.0f} A"
    )
    ax.text(
        0.03,
        0.05,
        info,
        transform=ax.transAxes,
        fontsize=7.8,
        color=PG_YELLOW,
        ha="left",
        va="bottom",
        family="DejaVu Sans Mono",
        bbox={
            "boxstyle": "square,pad=0.25",
            "facecolor": PG_BG,
            "edgecolor": PG_FG,
            "alpha": 0.95,
        },
    )

    sed_handles = [sed_line, bb_line, wien_line, peak_point] + filter_handles
    sed_labels = [handle.get_label() for handle in sed_handles]
    sed_legend = ax.legend(
        sed_handles,
        sed_labels,
        fontsize=7.6,
        loc="lower right",
        framealpha=1.0,
        title="legend",
        title_fontsize=8,
    )
    sed_legend.get_frame().set_facecolor(PG_BG)
    sed_legend.get_frame().set_edgecolor(PG_FG)
    for text in sed_legend.get_texts():
        text.set_color(PG_FG)
    sed_legend.get_title().set_color(PG_YELLOW)

    fig.tight_layout(rect=[0.01, 0.02, 0.99, 0.95])
    fig.savefig(OUT_FILE, dpi=180, bbox_inches="tight")
    print(f"saved: {OUT_FILE}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        sys.exit(f"ERROR: {exc}")
