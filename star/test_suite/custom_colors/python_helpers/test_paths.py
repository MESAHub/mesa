#!/usr/bin/env python3
"""
path_test.py
============
Tests all path resolution cases for the MESA colors module's instrument,
stellar_atm, vega_sed, and colors_results_directory parameters.

For each test the script:
  1. Copies (not symlinks) the data to the appropriate location
  2. Pre-creates the output directory (belt-and-suspenders alongside the
     Fortran mkdir-p fix in open_iteration_file)
  3. Rewrites inlist_colors with the test paths
  4. Runs 'make run' long enough to pass colors initialisation
  5. Saves three proof files to the test output dir:
       proof_plot.png   — CMD + HR diagram + filter light curves
       proof_newton.png — Teff and magnitude vs Newton iteration number
       proof_sed.png    — sample filter SED convolution curves
  6. Does NOT clean up — everything is left on disk for inspection

Run from the custom_colors test suite directory:
    python3 path_test.py

Requirements:
  - Run './mk' first to build
  - MESA_DIR environment variable must be set
  - mesa_reader, matplotlib, pandas, numpy must be importable
"""

import os
import re
import shutil
import subprocess
import textwrap
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    import pandas as pd

    HAVE_PD = True
except ImportError:
    HAVE_PD = False
    print("WARNING: pandas not found — SED plots will be skipped.")

try:
    import mesa_reader as mr

    HAVE_MR = True
except ImportError:
    HAVE_MR = False
    print("WARNING: mesa_reader not found — history-based plots will be skipped.")

# =============================================================================
# CONFIGURATION
# =============================================================================
os.chdir("../")

_mesa = os.environ.get("MESA_DIR", "")

SOURCE_INSTRUMENT = os.path.join(_mesa, "data/colors_data/filters/Generic/Johnson")
SOURCE_STELLAR_ATM = os.path.join(
    _mesa, "data/colors_data/stellar_models/Kurucz2003all"
)
SOURCE_VEGA_SED = os.path.join(_mesa, "data/colors_data/stellar_models/vega_flam.csv")

INLIST_COLORS = "inlist_colors"
LOGS_DIR = "LOGS"
RUN_TIMEOUT = 45

SUCCESS_MARKERS = ["step    1", "step      1", "      1   "]
FAILURE_MARKERS = [
    "Error: Could not open file",
    "colors_utils.f90",
    "ERROR STOP",
    "STOP 1",
    "Backtrace",
    "failed to find",
    "ERROR: failed",
]

# =============================================================================
# HELPERS
# =============================================================================

INLIST_BACKUP = INLIST_COLORS + ".bak"
CWD = os.getcwd()

MESA_DIR = os.environ.get("MESA_DIR", "")
if not MESA_DIR:
    raise EnvironmentError("MESA_DIR is not set.")


def backup_inlist():
    shutil.copy2(INLIST_COLORS, INLIST_BACKUP)


def restore_inlist():
    shutil.copy2(INLIST_BACKUP, INLIST_COLORS)


def resolve_outdir(output_dir: str) -> Path:
    """Return an absolute Path for output_dir (handles relative and ../)."""
    if os.path.isabs(output_dir):
        return Path(output_dir)
    return Path(os.path.normpath(os.path.join(CWD, output_dir)))


def patch_inlist(instrument, stellar_atm, vega_sed, colors_results_directory):
    with open(INLIST_COLORS, "r") as f:
        text = f.read()

    def replace_param(text, key, value):
        pattern = rf"(^\s*{re.escape(key)}\s*=\s*)['\"].*?['\"]"
        replacement = rf"\g<1>'{value}'"
        new_text, n = re.subn(pattern, replacement, text, flags=re.MULTILINE)
        if n == 0:
            raise ValueError(f"Could not find parameter '{key}' in {INLIST_COLORS}")
        return new_text

    text = replace_param(text, "instrument", instrument)
    text = replace_param(text, "stellar_atm", stellar_atm)
    text = replace_param(text, "vega_sed", vega_sed)
    text = replace_param(text, "colors_results_directory", colors_results_directory)

    with open(INLIST_COLORS, "w") as f:
        f.write(text)


def copy_data(src, dst):
    dst_path = Path(dst)
    dst_path.parent.mkdir(parents=True, exist_ok=True)
    if dst_path.exists() or dst_path.is_symlink():
        shutil.rmtree(dst) if dst_path.is_dir() else dst_path.unlink()
    shutil.copytree(src, dst) if Path(src).is_dir() else shutil.copy2(src, dst)


def run_star(timeout):
    try:
        proc = subprocess.Popen(
            ["make", "run"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        try:
            output, _ = proc.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            output, _ = proc.communicate()
    except FileNotFoundError:
        return False, "ERROR: 'make' not found."

    found_failure = any(m in output for m in FAILURE_MARKERS)
    found_success = any(m in output for m in SUCCESS_MARKERS)
    return found_success and not found_failure, output


def print_result(name, success, output):
    status = "✓ PASS" if success else "✗ FAIL"
    print(f"\n{'=' * 70}")
    print(f"  {status}  |  {name}")
    print(f"{'=' * 70}")
    lines = output.strip().splitlines()
    relevant = [
        l
        for l in lines
        if any(
            kw.lower() in l.lower()
            for kw in [
                "color",
                "filter",
                "error",
                "step",
                "instrument",
                "could not open",
                "backtrace",
                "stop",
                "results",
            ]
        )
    ]
    if relevant:
        print("  Relevant output:")
        for l in relevant[:25]:
            print(f"    {l}")
    else:
        print("  (no recognisable output — showing tail)")
        for l in lines[-15:]:
            print(f"    {l}")


# =============================================================================
# PROOF PLOTS
# =============================================================================

_PHASE_MAP = {
    -1: ("Relax", "#C0C0C0"),
    1: ("Starting", "#E6E6FA"),
    2: ("Pre-MS", "#FF69B4"),
    3: ("ZAMS", "#00FF00"),
    4: ("IAMS", "#0000FF"),
    5: ("TAMS", "#FF8C00"),
    6: ("He-Burn", "#8A2BE2"),
    7: ("ZACHeB", "#9932CC"),
    8: ("TACHeB", "#BA55D3"),
    9: ("TP-AGB", "#8B0000"),
    10: ("C-Burn", "#FF4500"),
    14: ("WDCS", "#708090"),
}


def _read_filter_columns(history_file):
    with open(history_file) as f:
        for line in f:
            if "model_number" in line:
                cols = line.split()
                try:
                    idx = cols.index("Interp_rad")
                    return cols[idx + 1 :]
                except ValueError:
                    pass
    return []


def _phase_colors(md_raw, skip):
    if hasattr(md_raw, "phase_of_evolution"):
        codes = md_raw.phase_of_evolution[skip:]
    else:
        codes = np.zeros(len(md_raw.model_number) - skip, dtype=int)
    names = [_PHASE_MAP.get(int(c), ("Unknown", "#808080"))[0] for c in codes]
    colors = [_PHASE_MAP.get(int(c), ("Unknown", "#808080"))[1] for c in codes]
    return names, colors


def _best_color_pair(filter_columns, md_raw, skip):
    fc = {f.lower(): f for f in filter_columns}

    def get(name):
        real = fc.get(name.lower())
        if real is None:
            return None
        try:
            v = getattr(md_raw, real)
        except AttributeError:
            try:
                v = md_raw.data(real)
            except Exception:
                return None
        return v[skip:] if isinstance(v, np.ndarray) and len(v) > skip else v

    for blue, red, mag_name in [
        ("gbp", "grp", "g"),
        ("b", "r", "v"),
        ("b", "v", "v"),
        ("v", "r", "v"),
        ("g", "r", "g"),
    ]:
        b, r, m = get(blue), get(red), get(mag_name)
        if b is not None and r is not None and m is not None:
            return b - r, m, f"{blue.upper()} − {red.upper()}", mag_name.upper()

    if len(filter_columns) >= 2:
        d0, d1 = get(filter_columns[0]), get(filter_columns[-1])
        if d0 is not None and d1 is not None:
            return (
                d0 - d1,
                d0,
                f"{filter_columns[0]} − {filter_columns[-1]}",
                filter_columns[0],
            )

    return None, None, "Color", "Mag"


def plot_history_result(test_name, out_path, success):
    """
    proof_plot.png — 3 panels: CMD, HR diagram, filter light curves.
    Also copies history.data into the output directory.
    """
    if not HAVE_MR:
        return

    history_file = os.path.join(LOGS_DIR, "history.data")
    if not os.path.exists(history_file):
        print(f"  [proof_plot] {history_file} not found — skipping.")
        return

    shutil.copy2(history_file, out_path / "history.data")

    try:
        md_raw = mr.MesaData(history_file)
    except Exception as e:
        print(f"  [proof_plot] mesa_reader failed: {e}")
        return

    skip = 5
    filter_columns = _read_filter_columns(history_file)
    phase_names, p_cols = _phase_colors(md_raw, skip)
    color_idx, mag, xlabel, ylabel = _best_color_pair(filter_columns, md_raw, skip)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle(f"{test_name}  —  {'PASS ✓' if success else 'FAIL ✗'}", fontsize=13)

    ax = axes[0]
    if color_idx is not None:
        ax.scatter(color_idx, mag, c=p_cols, s=15, alpha=0.7, edgecolors="none")
        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.invert_yaxis()
        ax.xaxis.set_ticks_position("top")
        ax.xaxis.set_label_position("top")
    else:
        ax.text(
            0.5, 0.5, "No colour pair", ha="center", va="center", transform=ax.transAxes
        )
    ax.set_title("Colour–Magnitude", fontsize=11)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    try:
        ax.scatter(
            md_raw.Teff[skip:],
            md_raw.log_L[skip:],
            c=p_cols,
            s=15,
            alpha=0.7,
            edgecolors="none",
        )
        ax.set_xlabel("Teff (K)", fontsize=11)
        ax.set_ylabel("log L/L☉", fontsize=11)
        ax.invert_xaxis()
        ax.xaxis.set_ticks_position("top")
        ax.xaxis.set_label_position("top")
    except Exception as e:
        ax.text(
            0.5,
            0.5,
            f"HR error:\n{e}",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=8,
        )
    ax.set_title("HR Diagram", fontsize=11)
    ax.grid(True, alpha=0.3)

    ax = axes[2]
    try:
        age = md_raw.star_age[skip:]
        plotted = 0
        for filt in filter_columns:
            try:
                data = getattr(md_raw, filt)[skip:]
            except AttributeError:
                try:
                    data = md_raw.data(filt)[skip:]
                except Exception:
                    continue
            ax.plot(age, data, lw=1.2, label=filt, alpha=0.85)
            plotted += 1
        if plotted:
            ax.invert_yaxis()
            ax.legend(fontsize=8, ncol=2)
        else:
            ax.text(
                0.5,
                0.5,
                "No filter data",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
        ax.set_xlabel("Age (yr)", fontsize=11)
        ax.set_ylabel("Magnitude", fontsize=11)
    except Exception as e:
        ax.text(
            0.5,
            0.5,
            f"LC error:\n{e}",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=8,
        )
    ax.set_title("Light Curves", fontsize=11)
    ax.grid(True, alpha=0.3)

    seen = {}
    for name, color in zip(phase_names, p_cols):
        seen.setdefault(name, color)
    legend_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=c,
            markersize=8,
            label=n,
            markeredgecolor="none",
        )
        for n, c in seen.items()
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=len(seen),
        fontsize=9,
        title="Evolutionary phase",
        title_fontsize=9,
        bbox_to_anchor=(0.5, -0.04),
    )

    plt.tight_layout(rect=[0, 0.06, 1, 1])
    out_png = out_path / "proof_plot.png"
    plt.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  [proof_plot]   saved → {out_png}")


def plot_newton_iter_result(test_name, out_path, success):
    """
    proof_newton.png — 2 panels from iteration_colors.data:
      Left:  Newton iteration vs Teff, coloured by model_number
      Right: Newton iteration vs first available filter magnitude, coloured by model_number
    """
    iter_file = out_path / "iteration_colors.data"
    if not iter_file.exists():
        print(f"  [proof_newton] {iter_file} not found — skipping.")
        return

    with open(iter_file) as f:
        header_line = f.readline().strip().lstrip("#").strip()
    col_names = header_line.split()

    try:
        data = np.loadtxt(iter_file, comments="#")
    except Exception as e:
        print(f"  [proof_newton] Failed to load {iter_file}: {e}")
        return
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[0] < 2:
        print(f"  [proof_newton] Too few rows ({data.shape[0]}) — skipping.")
        return

    def col(name):
        try:
            return data[:, col_names.index(name)]
        except ValueError:
            return None

    model_col = col("model")
    iter_col = col("iter")
    teff_col = col("Teff")

    if iter_col is None or teff_col is None or model_col is None:
        print(
            "  [proof_newton] Expected columns (model/iter/Teff) not found — skipping."
        )
        return

    # First non-sentinel filter magnitude column (after Flux_bol)
    mag_col_data, mag_col_label = None, None
    try:
        flux_idx = col_names.index("Flux_bol")
        for candidate in col_names[flux_idx + 1 :]:
            c = col(candidate)
            if c is not None and not np.all(c == -99.0):
                mag_col_data = c
                mag_col_label = candidate
                break
    except ValueError:
        pass

    n_panels = 2 if mag_col_data is not None else 1
    fig, axes = plt.subplots(1, n_panels, figsize=(6 * n_panels, 5))
    if n_panels == 1:
        axes = [axes]
    fig.suptitle(
        f"Newton Iterations  |  {test_name}  —  {'PASS ✓' if success else 'FAIL ✗'}",
        fontsize=12,
    )

    sc = axes[0].scatter(
        iter_col,
        teff_col,
        c=model_col,
        cmap="viridis",
        s=15,
        alpha=0.7,
        edgecolors="none",
    )
    axes[0].set_xlabel("Newton iteration", fontsize=11)
    axes[0].set_ylabel("Teff (K)", fontsize=11)
    axes[0].set_title("Teff per Newton iteration", fontsize=11)
    axes[0].grid(True, alpha=0.3)
    plt.colorbar(sc, ax=axes[0], label="model number")

    if mag_col_data is not None:
        sc2 = axes[1].scatter(
            iter_col,
            mag_col_data,
            c=model_col,
            cmap="plasma",
            s=15,
            alpha=0.7,
            edgecolors="none",
        )
        axes[1].set_xlabel("Newton iteration", fontsize=11)
        axes[1].set_ylabel(f"{mag_col_label} (mag)", fontsize=11)
        axes[1].set_title(f"{mag_col_label} per Newton iteration", fontsize=11)
        axes[1].invert_yaxis()
        axes[1].grid(True, alpha=0.3)
        plt.colorbar(sc2, ax=axes[1], label="model number")

    plt.tight_layout()
    out_png = out_path / "proof_newton.png"
    plt.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  [proof_newton] saved → {out_png}")


def plot_sed_result(test_name, out_path, success):
    """
    proof_sed.png — up to 5 filter SED convolution files (*_SED.csv, excluding VEGA_*)
    overlaid on the full stellar SED.
    """
    if not HAVE_PD:
        return

    sed_files = sorted(
        [
            f
            for f in out_path.iterdir()
            if f.name.endswith("_SED.csv") and not f.name.startswith("VEGA")
        ]
    )

    if not sed_files:
        print(f"  [proof_sed]    No *_SED.csv files in {out_path} — skipping.")
        return

    fig, ax = plt.subplots(figsize=(10, 5))
    fig.suptitle(
        f"SED Sample  |  {test_name}  —  {'PASS ✓' if success else 'FAIL ✗'}",
        fontsize=12,
    )

    sed_plotted = False
    for f in sed_files[:5]:
        try:
            df = (
                pd.read_csv(f, delimiter=",", header=0)
                .rename(columns=str.strip)
                .dropna()
            )
        except Exception as e:
            print(f"  [proof_sed]    Could not read {f.name}: {e}")
            continue

        wavelengths = df.get("wavelengths", pd.Series()).to_numpy()
        flux = df.get("fluxes", pd.Series()).to_numpy()
        convolved_flux = df.get("convolved_flux", pd.Series()).to_numpy()

        if not sed_plotted and len(wavelengths) > 0 and len(flux) > 0:
            ax.plot(
                wavelengths, flux, color="black", lw=1.5, label="Full SED", zorder=5
            )
            sed_plotted = True

        filter_label = f.stem.replace("_SED", "")
        if len(wavelengths) > 0 and len(convolved_flux) > 0:
            ax.plot(
                wavelengths,
                convolved_flux,
                lw=1.2,
                label=f"{filter_label} (convolved)",
                alpha=0.85,
            )

    ax.set_xlabel("Wavelength (Å)", fontsize=11)
    ax.set_ylabel("Flux", fontsize=11)
    ax.set_title("Filter convolution check", fontsize=11)
    ax.legend(loc="best", fontsize=9)
    ax.ticklabel_format(style="plain", useOffset=False)
    ax.set_xlim([0, 60000])
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    out_png = out_path / "proof_sed.png"
    plt.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  [proof_sed]    saved → {out_png}")


def run_all_proof_plots(test_name, output_dir, success):
    out_path = resolve_outdir(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    plot_history_result(test_name, out_path, success)
    plot_newton_iter_result(test_name, out_path, success)
    plot_sed_result(test_name, out_path, success)


# =============================================================================
# DEFAULT INPUT PATHS
# =============================================================================

_DEFAULT_INSTRUMENT = "data/colors_data/filters/Generic/Johnson"
_DEFAULT_STELLAR_ATM = "data/colors_data/stellar_models/Kurucz2003all/"
_DEFAULT_VEGA_SED = "data/colors_data/stellar_models/vega_flam.csv"


def _run_test(name, instrument, stellar_atm, vega_sed, out):
    """Pre-create output dir, patch inlist, run, plot."""
    resolve_outdir(out).mkdir(parents=True, exist_ok=True)
    patch_inlist(instrument, stellar_atm, vega_sed, out)
    success, output = run_star(RUN_TIMEOUT)
    print_result(name, success, output)
    run_all_proof_plots(name, out, success)
    return success


# =============================================================================
# INPUT PATH TESTS
# =============================================================================


def test01_input_mesa_dir_relative():
    return _run_test(
        "INPUT: MESA_DIR-relative",
        "data/colors_data/filters/Generic/Johnson",
        "data/colors_data/stellar_models/Kurucz2003all/",
        "data/colors_data/stellar_models/vega_flam.csv",
        "test_results/test01_output",
    )


def test02_input_cwd_dotslash():
    stage = "./test_staged/test02_input_cwd_dotslash"
    copy_data(SOURCE_INSTRUMENT, f"{stage}/filters/Generic/Johnson")
    copy_data(SOURCE_STELLAR_ATM, f"{stage}/stellar_models/Kurucz2003all")
    copy_data(SOURCE_VEGA_SED, f"{stage}/stellar_models/vega_flam.csv")
    return _run_test(
        "INPUT: CWD-relative (./)",
        f"{stage}/filters/Generic/Johnson",
        f"{stage}/stellar_models/Kurucz2003all",
        f"{stage}/stellar_models/vega_flam.csv",
        "test_results/test02_output",
    )


def test03_input_cwd_dotdotslash():
    stage = "../test_staged/test03_input_cwd_dotdotslash"
    copy_data(SOURCE_INSTRUMENT, f"{stage}/filters/Generic/Johnson")
    copy_data(SOURCE_STELLAR_ATM, f"{stage}/stellar_models/Kurucz2003all")
    copy_data(SOURCE_VEGA_SED, f"{stage}/stellar_models/vega_flam.csv")
    return _run_test(
        "INPUT: CWD-relative (../)",
        f"{stage}/filters/Generic/Johnson",
        f"{stage}/stellar_models/Kurucz2003all",
        f"{stage}/stellar_models/vega_flam.csv",
        "test_results/test03_output",
    )


def test04_input_absolute():
    stage = os.path.join(CWD, "test_staged/test04_input_absolute")
    copy_data(SOURCE_INSTRUMENT, os.path.join(stage, "filters/Generic/Johnson"))
    copy_data(SOURCE_STELLAR_ATM, os.path.join(stage, "stellar_models/Kurucz2003all"))
    copy_data(SOURCE_VEGA_SED, os.path.join(stage, "stellar_models/vega_flam.csv"))
    return _run_test(
        "INPUT: Absolute path",
        os.path.join(stage, "filters/Generic/Johnson"),
        os.path.join(stage, "stellar_models/Kurucz2003all"),
        os.path.join(stage, "stellar_models/vega_flam.csv"),
        "test_results/test04_output",
    )


def test05_input_slash_mesa_fallback():
    return _run_test(
        "INPUT: /-prefixed MESA_DIR fallback",
        "/data/colors_data/filters/Generic/Johnson",
        "/data/colors_data/stellar_models/Kurucz2003all/",
        "/data/colors_data/stellar_models/vega_flam.csv",
        "test_results/test05_output",
    )


# =============================================================================
# OUTPUT PATH TESTS
# =============================================================================


def test06_output_plain_relative():
    return _run_test(
        "OUTPUT: Plain relative",
        _DEFAULT_INSTRUMENT,
        _DEFAULT_STELLAR_ATM,
        _DEFAULT_VEGA_SED,
        "test_results/test06_output",
    )


def test07_output_cwd_dotslash():
    return _run_test(
        "OUTPUT: CWD-relative (./)",
        _DEFAULT_INSTRUMENT,
        _DEFAULT_STELLAR_ATM,
        _DEFAULT_VEGA_SED,
        "./test_results/test07_output",
    )


def test08_output_dotdotslash():
    return _run_test(
        "OUTPUT: CWD-relative (../)",
        _DEFAULT_INSTRUMENT,
        _DEFAULT_STELLAR_ATM,
        _DEFAULT_VEGA_SED,
        "../test_results_parent/test08_output",
    )


def test09_output_absolute():
    return _run_test(
        "OUTPUT: Absolute path",
        _DEFAULT_INSTRUMENT,
        _DEFAULT_STELLAR_ATM,
        _DEFAULT_VEGA_SED,
        os.path.join(CWD, "test_results/test09_output"),
    )


def test10_output_nested_subdir():
    return _run_test(
        "OUTPUT: Nested subdirectory",
        _DEFAULT_INSTRUMENT,
        _DEFAULT_STELLAR_ATM,
        _DEFAULT_VEGA_SED,
        "test_results/test10_output/nested/SED",
    )


# =============================================================================
# MAIN
# =============================================================================

TESTS = [
    ("INPUT: MESA_DIR-relative", test01_input_mesa_dir_relative),
    ("INPUT: CWD-relative (./)", test02_input_cwd_dotslash),
    ("INPUT: CWD-relative (../)", test03_input_cwd_dotdotslash),
    ("INPUT: Absolute path", test04_input_absolute),
    ("INPUT: /-prefixed MESA fallback", test05_input_slash_mesa_fallback),
    ("OUTPUT: Plain relative", test06_output_plain_relative),
    ("OUTPUT: CWD-relative (./)", test07_output_cwd_dotslash),
    ("OUTPUT: CWD-relative (../)", test08_output_dotdotslash),
    ("OUTPUT: Absolute path", test09_output_absolute),
    ("OUTPUT: Nested subdirectory", test10_output_nested_subdir),
]

if __name__ == "__main__":
    print(
        textwrap.dedent(f"""
        MESA Colors — Path Resolution Test Suite
        =========================================
        Work directory : {CWD}
        MESA_DIR       : {MESA_DIR}
        Run timeout    : {RUN_TIMEOUT}s per test

        NOTE: No cleanup is performed. All staged input data and all output
        directories are left on disk for inspection after the run.

        Staged input data : ./test_staged/   and  ../test_staged/
        Output data       : ./test_results/  and  ../test_results_parent/
        Proof plots per test (in each output dir):
          proof_plot.png   — CMD + HR diagram + filter light curves
          proof_newton.png — Teff and filter mag vs Newton iteration
          proof_sed.png    — filter SED convolution sample
    """)
    )

    backup_inlist()

    results = {}
    try:
        for name, fn in TESTS:
            print(f"\nRunning: {name} ...")
            try:
                results[name] = fn()
            except Exception as e:
                print(f"  EXCEPTION: {e}")
                results[name] = False
            finally:
                restore_inlist()
                time.sleep(1)
    finally:
        restore_inlist()
        if os.path.exists(INLIST_BACKUP):
            os.remove(INLIST_BACKUP)

    print(f"\n\n{'=' * 70}")
    print("  SUMMARY")
    print(f"{'=' * 70}")
    for name, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {status}  {name}")
    total = len(results)
    passed = sum(results.values())
    print(f"\n  {passed}/{total} tests passed")
    print(f"{'=' * 70}")
    print(
        textwrap.dedent(f"""
        Directories to inspect:
          Input staging  : {CWD}/test_staged/
                           {os.path.dirname(CWD)}/test_staged/
          Output + plots : {CWD}/test_results/
                           {os.path.dirname(CWD)}/test_results_parent/
    """)
    )

    exit(0 if passed == total else 1)
