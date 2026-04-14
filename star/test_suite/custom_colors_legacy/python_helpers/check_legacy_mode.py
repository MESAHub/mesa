#!/usr/bin/env python3
"""Validate that a MESA colors test run used legacy precomputed photometry.

Expected layout when run from python_helpers/:
    ../inlist_colors
    ../LOGS/history.data
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except FileNotFoundError:
        raise SystemExit(f"ERROR: missing file: {path}")


def parse_photo_precompute(inlist_text: str) -> bool | None:
    match = re.search(r"^\s*photo_precompute\s*=\s*\.(true|false)\.\s*$",
                      inlist_text, re.IGNORECASE | re.MULTILINE)
    if match is None:
        return None
    return match.group(1).lower() == "true"


def parse_history_table(history_text: str) -> tuple[list[str], list[list[float]]]:
    lines = [line.rstrip() for line in history_text.splitlines() if line.strip()]

    names_idx = None
    for i, line in enumerate(lines):
        if "model_number" in line.split():
            names_idx = i
            break

    if names_idx is None:
        raise SystemExit("ERROR: could not locate history column names line containing 'model_number'")

    names = lines[names_idx].split()
    rows: list[list[float]] = []

    for line in lines[names_idx + 1 :]:
        parts = line.split()
        if len(parts) != len(names):
            continue
        try:
            row = [float(x.replace("D", "E").replace("d", "e")) for x in parts]
        except ValueError:
            continue
        rows.append(row)

    if not rows:
        raise SystemExit("ERROR: no numeric history rows found after the column-name line")

    return names, rows


def col_values(names: list[str], rows: list[list[float]], name: str) -> list[float]:
    try:
        idx = names.index(name)
    except ValueError:
        raise SystemExit(f"ERROR: required history column '{name}' not found")
    return [row[idx] for row in rows]


def is_all_close(values: list[float], target: float, tol: float = 1e-12) -> bool:
    return all(math.isfinite(v) and abs(v - target) <= tol for v in values)


def any_finite_non_sentinel(values: list[float], sentinel: float = -1.0) -> bool:
    for v in values:
        if math.isfinite(v) and abs(v - sentinel) > 1e-12:
            return True
    return False


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    inlist_path = root / "inlist_colors"
    history_path = root / "LOGS" / "history.data"

    inlist_text = read_text(inlist_path)
    photo_precompute = parse_photo_precompute(inlist_text)
    if photo_precompute is None:
        print("FAIL: could not find 'photo_precompute' in inlist_colors")
        return 1
    if not photo_precompute:
        print("FAIL: inlist_colors does not enable legacy mode (photo_precompute is not .true.)")
        return 1

    history_text = read_text(history_path)
    names, rows = parse_history_table(history_text)

    mag_bol = col_values(names, rows, "Mag_bol")
    flux_bol = col_values(names, rows, "Flux_bol")
    interp_rad = col_values(names, rows, "Interp_rad")

    if not any_finite_non_sentinel(mag_bol):
        print("FAIL: Mag_bol is missing or never populated with finite values")
        return 1

    if not is_all_close(flux_bol, -1.0):
        print("FAIL: Flux_bol is not the expected legacy sentinel (-1) for all rows")
        return 1

    if not is_all_close(interp_rad, -1.0):
        print("FAIL: Interp_rad is not the expected legacy sentinel (-1) for all rows")
        return 1

    try:
        interp_idx = names.index("Interp_rad")
    except ValueError:
        print("FAIL: Interp_rad column missing")
        return 1

    phot_cols = names[interp_idx + 1 :]
    if not phot_cols:
        print("FAIL: no photometric columns found after Interp_rad")
        return 1

    populated_phot_cols = []
    for name in phot_cols:
        values = col_values(names, rows, name)
        if any_finite_non_sentinel(values):
            populated_phot_cols.append(name)

    if not populated_phot_cols:
        print("FAIL: no populated photometric columns were found")
        return 1

    print("PASS: legacy colors signature detected")
    print("  - photo_precompute = .true.")
    print("  - Mag_bol populated")
    print("  - Flux_bol = -1 for all rows")
    print("  - Interp_rad = -1 for all rows")
    print(f"  - populated photometric columns found: {', '.join(populated_phot_cols[:8])}")
    if len(populated_phot_cols) > 8:
        print(f"    ... and {len(populated_phot_cols) - 8} more")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
