#!/usr/bin/env python3
"""
stitch_history.py — appends colors columns onto the original MESA history file.

Reads history_file and output_file from the existing &post_proc namelist in
the inlist one directory above this script. The colors output_file is expected
to have been produced by colors_post_proc. The original history file is renamed
with an _old suffix and replaced in-place with the stitched result.

Can be run standalone:
    python stitch_history.py

Or imported and called from another script:
    from stitch_history import stitch
    stitch(root_dir)
"""

from pathlib import Path
import re
import sys

HISTORY_HEADER_LINES = 6   # 5 MESA header lines + 1 column-name line
COLORS_SKIP_COLS     = 4   # star_age, Teff, log_g, feh (already in history)

SCRIPT_DIR = Path(__file__).resolve().parent


def stitch(root_dir: Path = None):
    """
    Stitch the colors post-proc output onto the MESA history file.

    Parameters
    ----------
    root_dir : Path, optional
        Project root (directory containing the inlist). Defaults to one
        directory above this script, matching standalone behaviour.

    Raises
    ------
    SystemExit on any fatal error so callers can catch it if needed.
    """
    if root_dir is None:
        root_dir = SCRIPT_DIR.parent

    root_dir = Path(root_dir).resolve()
    inlist_path = root_dir / "inlist"

    def resolve_from_root(path_str: str) -> Path:
        p = Path(path_str)
        return p if p.is_absolute() else root_dir / p

    try:
        text = inlist_path.read_text()
    except FileNotFoundError:
        sys.exit(f"error: inlist not found: {inlist_path}")

    m_hist = re.search(r"history_file\s*=\s*'([^']+)'", text)
    m_out  = re.search(r"output_file\s*=\s*'([^']+)'", text)

    history_file = resolve_from_root(m_hist.group(1) if m_hist else "LOGS/history.data")
    colors_file  = resolve_from_root(m_out.group(1)  if m_out  else "post_proc_output.data")

    # --- read original history file ---
    try:
        hist_lines = history_file.read_text().splitlines(keepends=True)
    except FileNotFoundError:
        sys.exit(f"error: history file not found: {history_file}")

    hist_header = hist_lines[:HISTORY_HEADER_LINES]
    hist_data   = [l for l in hist_lines[HISTORY_HEADER_LINES:] if l.strip()]
    hist_cols   = hist_header[HISTORY_HEADER_LINES - 1].split()

    # --- read colors output file ---
    try:
        col_lines = colors_file.read_text().splitlines(keepends=True)
    except FileNotFoundError:
        sys.exit(f"error: colors output file not found: {colors_file}")

    if not col_lines:
        sys.exit(f"error: colors output file is empty: {colors_file}")

    col_header = col_lines[0]
    col_data   = [l for l in col_lines[1:] if l.strip()]
    col_names  = col_header.split()[COLORS_SKIP_COLS:]

    # --- check row counts match ---
    if len(hist_data) != len(col_data):
        sys.exit(
            f"error: row count mismatch — "
            f"history has {len(hist_data)} models, "
            f"colors output has {len(col_data)}"
        )

    # --- build stitched header ---
    stitched_col_header = (
        hist_header[HISTORY_HEADER_LINES - 1].rstrip()
        + "   " + "   ".join(col_names) + "\n"
    )
    stitched_header = hist_header[:HISTORY_HEADER_LINES - 1] + [stitched_col_header]

    # --- stitch rows ---
    stitched_rows = []
    for h_line, c_line in zip(hist_data, col_data):
        color_vals = c_line.split()[COLORS_SKIP_COLS:]
        stitched_rows.append(h_line.rstrip() + "   " + "   ".join(color_vals) + "\n")

    # --- rename original and write stitched file in its place ---
    old_history = history_file.with_name(history_file.name + "_old")
    history_file.rename(old_history)

    with history_file.open("w") as f:
        f.writelines(stitched_header)
        f.writelines(stitched_rows)

    print(f"root    : {root_dir}")
    print(f"inlist  : {inlist_path}")
    print(f"history : {history_file} ({len(hist_data)} models, {len(hist_cols)} cols)")
    print(f"colors  : {colors_file} (+{len(col_names)} color cols)")
    print(f"renamed : {history_file.name} -> {old_history.name}")
    print(f"wrote   : {history_file}")


if __name__ == "__main__":
    stitch()
