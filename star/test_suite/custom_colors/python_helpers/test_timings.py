#!/usr/bin/env python3
"""
timing_benchmark.py
===================
Benchmarks the wall-clock overhead of every user-facing MESA colors
configuration against a no-colors baseline.

Each test runs the same stellar evolution model (7 M☉, pre-MS → ZAMS)
for a fixed number of timesteps.  Wall times are measured, normalised
to the baseline, and summarised in a formatted terminal table and a
set of publication-quality plots saved to ``test_results/timing/``.

Configurations tested
---------------------
  baseline          – use_colors = .false.  (pure MESA, no photometry)
  colors_minimal    – use_colors = .true.   (Vega, no CSV, no Newton)
  colors_AB         – mag_system = 'AB'
  colors_ST         – mag_system = 'ST'
  colors_csv        – + make_csv = .true.
  colors_csv_permod – + make_csv + sed_per_model = .true.
  colors_newton     – + colors_per_newton_step = .true.
  colors_full       – make_csv + sed_per_model + newton_step (everything on)

Run from the custom_colors test suite directory:
    python3 python_helpers/timing_benchmark.py

Requirements
------------
  - ./mk must have been run first
  - matplotlib, numpy must be importable
  - MESA_DIR environment variable must be set
"""

import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# =============================================================================
# SETUP – run from the work directory regardless of where the script lives
# =============================================================================

_SCRIPT_DIR = Path(__file__).resolve().parent
_WORK_DIR   = _SCRIPT_DIR.parent
os.chdir(_WORK_DIR)

MESA_DIR = os.environ.get("MESA_DIR", "")
if not MESA_DIR:
    sys.exit("ERROR: MESA_DIR is not set.")

# =============================================================================
# BENCHMARK PARAMETERS
# =============================================================================

# Number of converged timesteps each run must complete before timing stops.
# 60 steps gives ~30–90 s per config on a modern workstation, which is long
# enough to average out startup noise while keeping the full suite tractable.
N_STEPS = 1000

# Wall-clock timeout per run (seconds).  Must be >> expected run time.
RUN_TIMEOUT = 600

# Where plots and the CSV summary are written.
OUT_DIR = Path("test_results/timing")

# Inlist file that carries the &colors namelist (and &star_job / &controls).
INLIST_COLORS = "inlist_colors"
INLIST_MAIN   = "inlist"          # top-level inlist read by ./rn

# =============================================================================
# TEST CONFIGURATIONS
# =============================================================================
# Each entry is (label, description, dict-of-colors-parameters).
# Parameters not listed here keep the values already in inlist_colors.
# Boolean values must be given as the Fortran strings '.true.' / '.false.'.

CONFIGS = [
    (
        "baseline",
        "No colors (pure MESA)",
        {
            "use_colors":             ".false.",
        },
    ),
    (
        "colors_minimal",
        "Colors on – Vega, no CSV, no Newton",
        {
            "use_colors":             ".true.",
            "mag_system":             "'Vega'",
            "make_csv":               ".false.",
            "sed_per_model":          ".false.",
            "colors_per_newton_step": ".false.",
        },
    ),
    (
        "colors_AB",
        "Colors on – AB magnitude system",
        {
            "use_colors":             ".true.",
            "mag_system":             "'AB'",
            "make_csv":               ".false.",
            "sed_per_model":          ".false.",
            "colors_per_newton_step": ".false.",
        },
    ),
    (
        "colors_ST",
        "Colors on – ST magnitude system",
        {
            "use_colors":             ".true.",
            "mag_system":             "'ST'",
            "make_csv":               ".false.",
            "sed_per_model":          ".false.",
            "colors_per_newton_step": ".false.",
        },
    ),
    (
        "colors_csv",
        "Colors on – Vega + make_csv",
        {
            "use_colors":             ".true.",
            "mag_system":             "'Vega'",
            "make_csv":               ".true.",
            "sed_per_model":          ".false.",
            "colors_per_newton_step": ".false.",
        },
    ),
    (
        "colors_csv_permod",
        "Colors on – Vega + make_csv + sed_per_model",
        {
            "use_colors":             ".true.",
            "mag_system":             "'Vega'",
            "make_csv":               ".true.",
            "sed_per_model":          ".true.",
            "colors_per_newton_step": ".false.",
        },
    ),
    (
        "colors_newton",
        "Colors on – Vega + colors_per_newton_step",
        {
            "use_colors":             ".true.",
            "mag_system":             "'Vega'",
            "make_csv":               ".false.",
            "sed_per_model":          ".false.",
            "colors_per_newton_step": ".true.",
        },
    ),
    (
        "colors_full",
        "Colors on – all options enabled",
        {
            "use_colors":             ".true.",
            "mag_system":             "'Vega'",
            "make_csv":               ".true.",
            "sed_per_model":          ".true.",
            "colors_per_newton_step": ".true.",
        },
    ),
]

# Colour palette for plots (one per config, reused across subplots)
PALETTE = [
    "#4c72b0",  # baseline  – blue
    "#dd8452",  # minimal   – orange
    "#55a868",  # AB        – green
    "#c44e52",  # ST        – red
    "#8172b2",  # csv       – purple
    "#937860",  # csv+per   – brown
    "#da8bc3",  # newton    – pink
    "#8c8c8c",  # full      – grey
]

# =============================================================================
# INLIST PATCHING
# =============================================================================

def _backup(path: str) -> str:
    bak = path + ".timing_bak"
    shutil.copy2(path, bak)
    return bak


def _restore(path: str, bak: str):
    shutil.copy2(bak, path)
    os.remove(bak)


def _set_param(text: str, param: str, value: str) -> str:
    """Replace ``param = value`` inside the &colors namelist block only.

    Extracts the &colors block, substitutes within it, then splices it
    back.  Commented-out occurrences of the parameter elsewhere in the
    file are never touched.  If the parameter is absent it is inserted
    before the closing slash.
    """
    block_re = re.compile(
        r"(&colors\b.*?^/)[ \t]*(?:!.*)?$",
        re.IGNORECASE | re.DOTALL | re.MULTILINE,
    )
    m = block_re.search(text)
    if not m:
        return text
    block = m.group(0)

    # Only match un-commented assignment lines.
    pattern = re.compile(
        r"^(?!\s*!)(\s*" + re.escape(param) + r"\s*=\s*)\S.*$",
        re.IGNORECASE | re.MULTILINE,
    )
    if pattern.search(block):
        new_block = pattern.sub(r"\g<1>" + value, block)
    else:
        new_block = re.sub(
            r"(\n/)",
            r"\n      " + param + " = " + value + r"\1",
            block,
            count=1,
        )
    return text[:m.start()] + new_block + text[m.end():]

def patch_inlist_colors(params: dict):
    """Rewrite inlist_colors with the given &colors parameter values.

    Modifies &colors parameters, disables pgstar, and enforces a hard
    stop at N_STEPS by setting max_model_number and disabling any other
    stopping conditions that could fire first.
    """
    with open(INLIST_COLORS, "r") as fh:
        text = fh.read()

    for param, value in params.items():
        text = _set_param(text, param, value)

    # Disable interactive plotting.
    text = re.sub(r"(?i)(pgstar_flag\s*=\s*)\S+", r"\g<1>.false.", text)

    # Hard stop: max_model_number = N_STEPS.
    # Insert after &controls if not already present, otherwise replace.
    if re.search(r"max_model_number\s*=", text, re.IGNORECASE):
        text = re.sub(
            r"(?im)(max_model_number\s*=\s*)\S+",
            rf"\g<1>{N_STEPS}",
            text,
        )
    else:
        text = re.sub(
            r"(?m)^(&controls)",
            rf"\1\n      max_model_number = {N_STEPS}",
            text,
            count=1,
        )

    # Disable xa_central_lower_limit so He exhaustion can't stop the run.
    text = re.sub(
        r"(?im)(xa_central_lower_limit\(1\)\s*=\s*)\S+",
        r"\g<1>-1d0",
        text,
    )

    # Push max_age far out so it can't fire within 100 steps.
    text = re.sub(
        r"(?im)(max_age\s*=\s*)\S+",
        r"\g<1>1d99",
        text,
    )

    with open(INLIST_COLORS, "w") as fh:
        fh.write(text)


# =============================================================================
# RUN
# =============================================================================

def _count_steps(output: str) -> int:
    """Count the number of converged timestep lines in MESA terminal output.

    MESA step lines look like (9+ leading spaces):
        "         10   5.477806   4013.312 ..."
    Pre-MS relaxation lines use the same format but appear before the
    "finished doing relax num steps" marker.  We strip everything before
    that marker so only real evolution steps are counted.
    """
    # Discard pre-MS relaxation output so its step numbers aren't counted.
    marker = "finished doing relax num steps"
    idx = output.find(marker)
    if idx != -1:
        output = output[idx:]

    # Step lines: 1–20 leading spaces, integer step number, then whitespace
    # and a floating-point value (lg_Tmax).  \s+ covers the full range of
    # MESA's column-width indentation.
    matches = re.findall(r"^\s+(\d+)\s+\d+\.\d", output, re.MULTILINE)
    if not matches:
        return 0
    return max(int(m) for m in matches)


def run_config(label: str, params: dict) -> dict:
    """Patch inlists, run MESA, return timing metadata."""
    bak_colors = _backup(INLIST_COLORS)

    result = {
        "label":    label,
        "wall_s":   None,
        "steps":    0,
        "success":  False,
        "timeout":  False,
        "s_per_step": None,
    }

    try:
        patch_inlist_colors(params)

        t0 = time.perf_counter()
        proc = subprocess.run(
            ["./rn"],
            capture_output=True,
            text=True,
            timeout=RUN_TIMEOUT,
        )
        wall = time.perf_counter() - t0

        output = proc.stdout + proc.stderr
        steps  = _count_steps(output)

        failed = any(m in output for m in [
            "ERROR STOP", "STOP 1", "Backtrace",
            "failed to find", "ERROR: failed",
            "Error: Could not open",
        ])

        if steps == 0 or failed:
            print("  [DEBUG] MESA output (last 40 lines):")
            for line in output.splitlines()[-40:]:
                print("    " + line)

        result["wall_s"]   = wall
        result["steps"]    = steps
        result["success"]  = (not failed) and steps > 0
        result["s_per_step"] = (wall / steps) if steps > 0 else None

    except subprocess.TimeoutExpired:
        result["timeout"] = True
        print(f"  [TIMEOUT] {label} exceeded {RUN_TIMEOUT} s")

    finally:
        _restore(INLIST_COLORS, bak_colors)

    return result


# =============================================================================
# TERMINAL REPORTING
# =============================================================================

_PASS = "\033[92m PASS \033[0m"
_FAIL = "\033[91m FAIL \033[0m"
_TIME = "\033[93mTIMEOUT\033[0m"

_SEP  = "─" * 92

def print_header():
    print()
    print("╔" + "═" * 90 + "╗")
    print("║" + "  MESA COLORS  –  TIMING BENCHMARK".center(90) + "║")
    print("╚" + "═" * 90 + "╝")
    print(f"  Work dir   : {_WORK_DIR}")
    print(f"  MESA_DIR   : {MESA_DIR}")
    print(f"  Steps/run  : {N_STEPS}")
    print(f"  Timeout    : {RUN_TIMEOUT} s")
    print(f"  Configs    : {len(CONFIGS)}")
    print()


def print_running(idx: int, label: str, desc: str):
    print(_SEP)
    print(f"  [{idx}/{len(CONFIGS)}]  {label}")
    print(f"           {desc}")


def print_run_result(r: dict, baseline_sps: float | None):
    if r["timeout"]:
        status = _TIME
    elif r["success"]:
        status = _PASS
    else:
        status = _FAIL

    wall   = r["wall_s"]
    steps  = r["steps"]
    sps    = r["s_per_step"]

    wall_str  = f"{wall:7.1f} s" if wall  is not None else "    N/A  "
    steps_str = f"{steps:4d} steps"
    sps_str   = f"{sps:.3f} s/step" if sps is not None else "     N/A    "

    overhead = ""
    if sps is not None and baseline_sps is not None and baseline_sps > 0:
        pct = 100.0 * (sps - baseline_sps) / baseline_sps
        sign = "+" if pct >= 0 else ""
        overhead = f"  overhead {sign}{pct:.1f}%"

    print(f"  {status}  {wall_str}  {steps_str}  {sps_str}{overhead}")


def print_summary_table(results: list[dict]):
    baseline_sps = None
    for r in results:
        if r["label"] == "baseline" and r["s_per_step"] is not None:
            baseline_sps = r["s_per_step"]
            break

    print()
    print("╔" + "═" * 90 + "╗")
    print("║" + "  SUMMARY TABLE".center(90) + "║")
    print("╠" + "═" * 90 + "╣")
    hdr = f"  {'Config':<24} {'Status':<8} {'Wall (s)':>9} {'Steps':>6} {'s/step':>9} {'Overhead':>10}"
    print("║" + hdr + " " * (90 - len(hdr)) + "║")
    print("╠" + "═" * 90 + "╣")

    for r in results:
        if r["timeout"]:
            st = "TIMEOUT"
        elif r["success"]:
            st = "PASS"
        else:
            st = "FAIL"

        wall_s  = f"{r['wall_s']:.1f}"   if r["wall_s"]     is not None else "N/A"
        sps_s   = f"{r['s_per_step']:.3f}" if r["s_per_step"] is not None else "N/A"

        if r["s_per_step"] is not None and baseline_sps is not None and baseline_sps > 0:
            pct    = 100.0 * (r["s_per_step"] - baseline_sps) / baseline_sps
            oh_str = f"{pct:+.1f}%"
        else:
            oh_str = "–"

        row = (
            f"  {r['label']:<24} {st:<8} {wall_s:>9} "
            f"{r['steps']:>6} {sps_s:>9} {oh_str:>10}"
        )
        print("║" + row + " " * max(0, 90 - len(row)) + "║")

    print("╚" + "═" * 90 + "╝")
    print()


# =============================================================================
# PLOTS
# =============================================================================

def _valid(results: list[dict]) -> list[dict]:
    return [r for r in results if r["success"] and r["s_per_step"] is not None]


def plot_wall_time(results: list[dict], out: Path):
    valid = _valid(results)
    if not valid:
        print("  [plot_wall_time]  no successful runs – skipping")
        return

    labels  = [r["label"]  for r in valid]
    walls   = [r["wall_s"] for r in valid]
    colors  = [PALETTE[results.index(r)] for r in valid]

    fig, ax = plt.subplots(figsize=(11, 5))
    bars = ax.bar(labels, walls, color=colors, edgecolor="white", linewidth=0.8, zorder=3)

    # Annotate bar tops
    for bar, w in zip(bars, walls):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            w + max(walls) * 0.01,
            f"{w:.1f} s",
            ha="center", va="bottom", fontsize=8.5, fontweight="bold",
        )

    ax.set_ylabel("Wall-clock time (s)", fontsize=12)
    ax.set_title(
        f"Total wall-clock time per configuration  ({N_STEPS} steps each)",
        fontsize=13, fontweight="bold",
    )
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
    ax.set_ylim(0, max(walls) * 1.18)
    ax.yaxis.grid(True, alpha=0.35, zorder=0)
    ax.set_axisbelow(True)
    fig.tight_layout()

    fpath = out / "timing_wall_time.png"
    fig.savefig(fpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  [plot]  {fpath}")


def plot_overhead(results: list[dict], out: Path):
    valid = _valid(results)
    if not valid:
        return

    baseline_sps = next(
        (r["s_per_step"] for r in valid if r["label"] == "baseline"), None
    )
    if baseline_sps is None:
        print("  [plot_overhead]  baseline missing – skipping overhead plot")
        return

    non_baseline = [r for r in valid if r["label"] != "baseline"]
    if not non_baseline:
        return

    labels   = [r["label"] for r in non_baseline]
    overhead = [
        100.0 * (r["s_per_step"] - baseline_sps) / baseline_sps
        for r in non_baseline
    ]
    colors = [PALETTE[results.index(r)] for r in non_baseline]

    fig, ax = plt.subplots(figsize=(10, 5))
    bars = ax.bar(labels, overhead, color=colors, edgecolor="white", linewidth=0.8, zorder=3)

    for bar, pct in zip(bars, overhead):
        va    = "bottom" if pct >= 0 else "top"
        yoffs = max(overhead) * 0.02 if pct >= 0 else min(overhead) * 0.02
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            pct + yoffs,
            f"{pct:+.1f}%",
            ha="center", va=va, fontsize=8.5, fontweight="bold",
        )

    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Overhead vs baseline (%)", fontsize=12)
    ax.set_title(
        "Per-step overhead relative to no-colors baseline",
        fontsize=13, fontweight="bold",
    )
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
    y_pad = max(abs(v) for v in overhead) * 0.18
    ax.set_ylim(min(min(overhead) - y_pad, -2), max(overhead) + y_pad)
    ax.yaxis.grid(True, alpha=0.35, zorder=0)
    ax.set_axisbelow(True)
    fig.tight_layout()

    fpath = out / "timing_overhead_pct.png"
    fig.savefig(fpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  [plot]  {fpath}")


def plot_s_per_step(results: list[dict], out: Path):
    """Horizontal bar chart of seconds-per-step, sorted ascending."""
    valid = _valid(results)
    if not valid:
        return

    valid_sorted = sorted(valid, key=lambda r: r["s_per_step"])
    labels = [r["label"]      for r in valid_sorted]
    sps    = [r["s_per_step"] for r in valid_sorted]
    clrs   = [PALETTE[results.index(r)] for r in valid_sorted]

    fig, ax = plt.subplots(figsize=(9, 0.7 * len(labels) + 2))
    bars = ax.barh(labels, sps, color=clrs, edgecolor="white", linewidth=0.8, zorder=3)

    for bar, v in zip(bars, sps):
        ax.text(
            v + max(sps) * 0.01,
            bar.get_y() + bar.get_height() / 2,
            f"{v:.3f} s",
            va="center", fontsize=8.5, fontweight="bold",
        )

    ax.set_xlabel("Seconds per converged timestep", fontsize=12)
    ax.set_title("Per-step cost by configuration", fontsize=13, fontweight="bold")
    ax.set_xlim(0, max(sps) * 1.20)
    ax.xaxis.grid(True, alpha=0.35, zorder=0)
    ax.set_axisbelow(True)
    fig.tight_layout()

    fpath = out / "timing_s_per_step.png"
    fig.savefig(fpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  [plot]  {fpath}")


def plot_feature_breakdown(results: list[dict], out: Path):
    """Grouped bar chart decomposing the effect of each feature flag."""
    valid  = _valid(results)
    bline  = next((r for r in valid if r["label"] == "baseline"),  None)
    minimal = next((r for r in valid if r["label"] == "colors_minimal"), None)

    if bline is None or minimal is None:
        print("  [plot_feature_breakdown]  missing baseline or minimal – skipping")
        return

    # Feature contributions relative to minimal-colors (all flags off)
    def delta(label):
        r = next((x for x in valid if x["label"] == label), None)
        if r is None:
            return None
        return (r["s_per_step"] - minimal["s_per_step"])

    features = {
        "colors\n(vs baseline)": minimal["s_per_step"] - bline["s_per_step"],
        "make_csv":              delta("colors_csv"),
        "sed_per_model\n(on top of csv)": (
            (next((x["s_per_step"] for x in valid if x["label"] == "colors_csv_permod"), None) or 0)
            - (next((x["s_per_step"] for x in valid if x["label"] == "colors_csv"),     None) or 0)
        ),
        "AB system":     delta("colors_AB"),
        "ST system":     delta("colors_ST"),
        "newton_step":   delta("colors_newton"),
    }

    labels = [k for k, v in features.items() if v is not None]
    values = [v for v in features.values() if v is not None]
    clrs   = ["#c44e52" if v > 0 else "#55a868" for v in values]

    fig, ax = plt.subplots(figsize=(10, 5))
    bars = ax.bar(labels, values, color=clrs, edgecolor="white", linewidth=0.8, zorder=3)

    for bar, v in zip(bars, values):
        va    = "bottom" if v >= 0 else "top"
        yoffs = max(abs(x) for x in values) * 0.02
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            v + (yoffs if v >= 0 else -yoffs),
            f"{v:+.4f} s",
            ha="center", va=va, fontsize=8, fontweight="bold",
        )

    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Δ seconds per step vs colors_minimal", fontsize=11)
    ax.set_title(
        "Incremental cost of individual features\n(positive = slower, negative = faster than colors_minimal)",
        fontsize=12, fontweight="bold",
    )
    positive_patch = mpatches.Patch(color="#c44e52", label="Added cost")
    negative_patch = mpatches.Patch(color="#55a868", label="Reduced cost / faster")
    ax.legend(handles=[positive_patch, negative_patch], fontsize=9)
    y_pad = max(abs(v) for v in values) * 0.20
    ax.set_ylim(min(values) - y_pad, max(values) + y_pad)
    ax.yaxis.grid(True, alpha=0.35, zorder=0)
    ax.set_axisbelow(True)
    fig.tight_layout()

    fpath = out / "timing_feature_breakdown.png"
    fig.savefig(fpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  [plot]  {fpath}")


def save_csv(results: list[dict], out: Path):
    fpath = out / "timing_results.csv"
    with open(fpath, "w") as fh:
        fh.write("label,description,wall_s,steps,s_per_step,success,timeout\n")
        for r, (_, desc, _) in zip(results, CONFIGS):
            wall = f"{r['wall_s']:.3f}"   if r["wall_s"]     is not None else ""
            sps  = f"{r['s_per_step']:.5f}" if r["s_per_step"] is not None else ""
            fh.write(
                f"{r['label']},{desc},{wall},{r['steps']},{sps},"
                f"{r['success']},{r['timeout']}\n"
            )
    print(f"  [csv]   {fpath}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print_header()

    results = []
    baseline_sps = None

    for idx, (label, desc, params) in enumerate(CONFIGS, start=1):
        print_running(idx, label, desc)
        r = run_config(label, params)
        results.append(r)
        print_run_result(r, baseline_sps)
        if label == "baseline" and r["s_per_step"] is not None:
            baseline_sps = r["s_per_step"]

    print(_SEP)
    print_summary_table(results)

    print("  Generating plots …")
    plot_wall_time(results, OUT_DIR)
    plot_overhead(results, OUT_DIR)
    plot_s_per_step(results, OUT_DIR)
    plot_feature_breakdown(results, OUT_DIR)
    save_csv(results, OUT_DIR)

    n_pass = sum(1 for r in results if r["success"])
    n_fail = len(results) - n_pass
    print()
    print(f"  Done.  {n_pass} passed, {n_fail} failed.")
    print(f"  Outputs written to: {OUT_DIR.resolve()}")
    print()


if __name__ == "__main__":
    main()
