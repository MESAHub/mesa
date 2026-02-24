#!/usr/bin/env python3
"""
test_path_resolution.py
=======================
Tests all path resolution cases for the MESA colors module's instrument,
stellar_atm, and vega_sed parameters. For each test case the script:
  1. Stages the data at the appropriate location (symlink or copy)
  2. Rewrites inlist_colors with the test paths
  3. Runs 'make run' long enough to pass colors initialisation
  4. Inspects output for success/failure
  5. Cleans up before the next test

Run from the custom_colors test suite directory:
    python3 test_path_resolution.py

Requirements:
  - Run './mk' first to build
  - MESA_DIR environment variable must be set
  - The canonical source paths below must exist on your system
"""

import os
import re
import shutil
import subprocess
import textwrap
import time
from pathlib import Path

# =============================================================================
# CONFIGURATION — edit these to match your system
# =============================================================================

# Canonical source locations (used as the data to stage for each test).
# All tests use the built-in Johnson+Kurucz data so the filter set is
# consistent with the pgstar inlist (which expects V, B, R, etc.).
_mesa = os.environ.get("MESA_DIR", "")
SOURCE_INSTRUMENT  = os.path.join(_mesa, "data/colors_data/filters/Generic/Johnson")
SOURCE_STELLAR_ATM = os.path.join(_mesa, "data/colors_data/stellar_models/Kurucz2003all")
SOURCE_VEGA_SED    = os.path.join(_mesa, "data/colors_data/stellar_models/vega_flam.csv")

# Path to the inlist_colors file (relative to this script's location)
INLIST_COLORS = "inlist_colors"

# How long to let 'make run' execute before killing it (seconds).
# Colors initialisation happens in the first few seconds.
RUN_TIMEOUT = 45

# Strings that indicate colors initialised and the run reached step 1
SUCCESS_MARKERS = [
    "step    1",
    "step      1",
    "      1   ",   # various MESA step-line formats
]

# Strings that indicate a colors path failure
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

MESA_DIR = os.environ.get("MESA_DIR", "")
if not MESA_DIR:
    raise EnvironmentError("MESA_DIR is not set. Source your MESA SDK before running this script.")


def backup_inlist():
    shutil.copy2(INLIST_COLORS, INLIST_BACKUP)


def restore_inlist():
    shutil.copy2(INLIST_BACKUP, INLIST_COLORS)


def patch_inlist(instrument: str, stellar_atm: str, vega_sed: str):
    """Replace the three path parameters in inlist_colors."""
    with open(INLIST_COLORS, "r") as f:
        text = f.read()

    def replace_param(text, key, value):
        # Matches:   key = 'anything'   or   key = "anything"
        pattern = rf"(^\s*{re.escape(key)}\s*=\s*)['\"].*?['\"]"
        replacement = rf"\g<1>'{value}'"
        new_text, n = re.subn(pattern, replacement, text, flags=re.MULTILINE)
        if n == 0:
            raise ValueError(f"Could not find parameter '{key}' in {INLIST_COLORS}")
        return new_text

    text = replace_param(text, "instrument",  instrument)
    text = replace_param(text, "stellar_atm", stellar_atm)
    text = replace_param(text, "vega_sed",    vega_sed)

    with open(INLIST_COLORS, "w") as f:
        f.write(text)


def make_symlink(src: str, dst: str):
    """Create a symlink at dst pointing to src, creating parent dirs as needed."""
    dst_path = Path(dst)
    dst_path.parent.mkdir(parents=True, exist_ok=True)
    if dst_path.exists() or dst_path.is_symlink():
        dst_path.unlink()
    os.symlink(os.path.abspath(src), dst)


def remove_path(path: str):
    p = Path(path)
    if p.is_symlink() or p.is_file():
        p.unlink()
    elif p.is_dir():
        shutil.rmtree(path)


def run_star(timeout: int) -> tuple[bool, str]:
    """
    Run 'make run', kill after timeout seconds, return (success, output).
    Success means a SUCCESS_MARKER was found and no FAILURE_MARKER was found.
    """
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

    success = found_success and not found_failure
    return success, output


def print_result(name: str, success: bool, output: str):
    status = "✓ PASS" if success else "✗ FAIL"
    print(f"\n{'='*70}")
    print(f"  {status}  |  {name}")
    print(f"{'='*70}")
    lines = output.strip().splitlines()
    relevant = [l for l in lines if any(
        kw.lower() in l.lower() for kw in
        ["color", "filter", "error", "step", "instrument", "could not open", "backtrace", "stop"]
    )]
    if relevant:
        print("  Relevant output:")
        for l in relevant[:25]:
            print(f"    {l}")
    else:
        print("  (no recognisable output — showing tail)")
        for l in lines[-15:]:
            print(f"    {l}")


# =============================================================================
# TEST CASES
# =============================================================================

def test_mesa_dir_relative():
    """
    Case 1: No leading / or ./ — resolved relative to $MESA_DIR.
    Uses the Johnson filters that ship with MESA so no staging needed.
    """
    patch_inlist(
        instrument  = "data/colors_data/filters/Generic/Johnson",
        stellar_atm = "data/colors_data/stellar_models/Kurucz2003all/",
        vega_sed    = "data/colors_data/stellar_models/vega_flam.csv",
    )
    success, output = run_star(RUN_TIMEOUT)
    print_result("MESA_DIR-relative (default built-in paths)", success, output)
    return success


def test_cwd_relative_dotslash():
    """
    Case 2: Paths beginning with ./ — resolved relative to the work directory.
    Symlinks source data into a local ./test_data/ subdirectory.
    """
    make_symlink(SOURCE_INSTRUMENT,  "./test_data/filters/Generic/Johnson")
    make_symlink(SOURCE_STELLAR_ATM, "./test_data/stellar_models/Kurucz2003all")
    make_symlink(SOURCE_VEGA_SED,    "./test_data/stellar_models/vega_flam.csv")

    patch_inlist(
        instrument  = "./test_data/filters/Generic/Johnson",
        stellar_atm = "./test_data/stellar_models/Kurucz2003all",
        vega_sed    = "./test_data/stellar_models/vega_flam.csv",
    )
    success, output = run_star(RUN_TIMEOUT)
    print_result("CWD-relative with ./", success, output)

    remove_path("./test_data")
    return success


def test_cwd_relative_dotdotslash():
    """
    Case 3: Paths beginning with ../ — resolved relative to the parent of the
    work directory. Symlinks data one level up.
    """
    make_symlink(SOURCE_INSTRUMENT,  "../colors_test_data/filters/Generic/Johnson")
    make_symlink(SOURCE_STELLAR_ATM, "../colors_test_data/stellar_models/Kurucz2003all")
    make_symlink(SOURCE_VEGA_SED,    "../colors_test_data/stellar_models/vega_flam.csv")

    patch_inlist(
        instrument  = "../colors_test_data/filters/Generic/Johnson",
        stellar_atm = "../colors_test_data/stellar_models/Kurucz2003all",
        vega_sed    = "../colors_test_data/stellar_models/vega_flam.csv",
    )
    success, output = run_star(RUN_TIMEOUT)
    print_result("CWD-relative with ../", success, output)

    remove_path("../colors_test_data")
    return success


def test_absolute_path():
    """
    Case 4: True absolute paths — used as-is when they exist on disk.
    No staging needed; canonical source paths are used directly.
    """
    patch_inlist(
        instrument  = os.path.abspath(SOURCE_INSTRUMENT),
        stellar_atm = os.path.abspath(SOURCE_STELLAR_ATM),
        vega_sed    = os.path.abspath(SOURCE_VEGA_SED),
    )
    success, output = run_star(RUN_TIMEOUT)
    print_result("Absolute path", success, output)
    return success


def test_slash_prefixed_mesa_fallback():
    """
    Case 5: A path starting with / that does NOT exist as a standalone
    absolute path — should fall back to $MESA_DIR-prepended resolution
    via the inquire() check in resolve_path.
    """
    patch_inlist(
        instrument  = "/data/colors_data/filters/Generic/Johnson",
        stellar_atm = "/data/colors_data/stellar_models/Kurucz2003all/",
        vega_sed    = "/data/colors_data/stellar_models/vega_flam.csv",
    )
    success, output = run_star(RUN_TIMEOUT)
    print_result("/-prefixed MESA_DIR-relative (fallback via inquire)", success, output)
    return success


# =============================================================================
# MAIN
# =============================================================================

TESTS = [
    ("MESA_DIR-relative",            test_mesa_dir_relative),
    ("CWD-relative (./)",            test_cwd_relative_dotslash),
    ("CWD-relative (../)",           test_cwd_relative_dotdotslash),
    ("Absolute path",                test_absolute_path),
    ("/-prefixed MESA_DIR fallback", test_slash_prefixed_mesa_fallback),
]

if __name__ == "__main__":
    print(textwrap.dedent(f"""
        MESA Colors — Path Resolution Test Suite
        =========================================
        Work directory : {os.getcwd()}
        MESA_DIR       : {MESA_DIR}
        Run timeout    : {RUN_TIMEOUT}s per test
    """))

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
                time.sleep(1)  # brief pause between runs
    finally:
        restore_inlist()
        # Remove backup
        if os.path.exists(INLIST_BACKUP):
            os.remove(INLIST_BACKUP)

    # Summary
    print(f"\n\n{'='*70}")
    print("  SUMMARY")
    print(f"{'='*70}")
    for name, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {status}  {name}")
    total = len(results)
    passed = sum(results.values())
    print(f"\n  {passed}/{total} tests passed")
    print(f"{'='*70}\n")

    exit(0 if passed == total else 1)