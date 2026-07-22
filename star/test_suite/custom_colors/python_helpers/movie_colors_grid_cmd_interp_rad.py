#!/usr/bin/env python3
"""
Two-panel movie for the MESA colors module.

Top panel:
    3D stellar-atmosphere grid in (Teff, logg, metallicity), using the real
    SED nodes from lookup_table.csv, plus the stellar evolution track.

Bottom panel:
    CMD coloured by Interp_rad.

This script is intentionally conservative:
    - it uses the actual lookup_table.csv nodes, not a fake Cartesian product
    - it uses history.data for Teff/logg/Interp_rad/CMD columns
    - if history.data has no metallicity column, it reads initial_z from the
      inlist and uses that as a constant metallicity coordinate
    - you can override the metallicity coordinate explicitly with
      --metallicity-coordinate
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import re
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import mesa_reader as mr
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from plot_history import MesaView, read_header_columns, setup_hr_diagram_params


FLOAT_RE = re.compile(
    r"""
    (?P<value>
        [+-]?
        (?:
            (?:\d+\.\d*) |
            (?:\.\d+) |
            (?:\d+)
        )
        (?:
            [dDeE][+-]?\d+
        )?
    )
    """,
    re.VERBOSE,
)


def fortran_float(value: str) -> float:
    """Parse a Fortran-style float such as 0.02d0."""
    return float(value.strip().replace("D", "e").replace("d", "e"))


def parse_inlist_value(inlist_file: str | Path, name: str) -> float | None:
    """Return the first numeric assignment for `name` in a MESA inlist."""
    path = Path(inlist_file)
    if not path.exists():
        return None

    pattern = re.compile(rf"^\s*{re.escape(name)}\s*=", re.IGNORECASE)

    for raw_line in path.read_text(errors="replace").splitlines():
        line = raw_line.split("!", 1)[0]
        if not pattern.search(line):
            continue

        rhs = line.split("=", 1)[1]
        match = FLOAT_RE.search(rhs)
        if match:
            return fortran_float(match.group("value"))

    return None


def read_lookup_table(lookup_file: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read real atmosphere grid nodes from lookup_table.csv."""
    lookup_file = Path(lookup_file)
    if not lookup_file.exists():
        raise FileNotFoundError(f"lookup_table.csv not found: {lookup_file}")

    teff: list[float] = []
    logg: list[float] = []
    meta: list[float] = []

    with lookup_file.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"teff", "logg", "metallicity"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"{lookup_file} is missing required columns: {sorted(missing)}"
            )

        for row in reader:
            try:
                t = float(row["teff"])
                g = float(row["logg"])
                m = float(row["metallicity"])
            except (TypeError, ValueError):
                continue

            if not (np.isfinite(t) and np.isfinite(g) and np.isfinite(m)):
                continue

            # Ignore dummy/sentinel values if they ever appear in real axes.
            if t in (999.0, -999.0) or g in (999.0, -999.0) or m in (999.0, -999.0):
                continue

            teff.append(t)
            logg.append(g)
            meta.append(m)

    if not teff:
        raise ValueError(f"No usable grid nodes found in {lookup_file}")

    return np.asarray(teff), np.asarray(logg), np.asarray(meta)


def mesa_colnames(md) -> set[str]:
    names: set[str] = set()

    for attr in ("bulk_names", "header_names"):
        vals = getattr(md, attr, None)
        if vals is not None:
            names.update(str(v) for v in vals)

    # MesaView wrappers may not expose all mesa_reader metadata.
    for key in getattr(md, "__dict__", {}):
        names.add(str(key))

    return names


def has_mesa_column(md, name: str) -> bool:
    if name in mesa_colnames(md):
        return True
    try:
        _ = getattr(md, name)
        return True
    except Exception:
        pass
    try:
        _ = md.data(name)
        return True
    except Exception:
        return False


def get_mesa_column(md, name: str) -> np.ndarray:
    try:
        return np.asarray(getattr(md, name), dtype=float)
    except Exception:
        pass

    try:
        return np.asarray(md.data(name), dtype=float)
    except Exception as exc:
        raise KeyError(f"Could not read history column {name!r}") from exc


def first_existing_column(md, candidates: Iterable[str]) -> tuple[str, np.ndarray] | None:
    for name in candidates:
        if has_mesa_column(md, name):
            return name, get_mesa_column(md, name)
    return None


def get_teff(md) -> np.ndarray:
    if has_mesa_column(md, "Teff"):
        return get_mesa_column(md, "Teff")
    if has_mesa_column(md, "log_Teff"):
        return 10.0 ** get_mesa_column(md, "log_Teff")
    raise KeyError("history.data needs either Teff or log_Teff")


def get_metallicity_track(
    md,
    n: int,
    *,
    inlist_file: str | Path | None,
    metallicity_coordinate: float | None,
    convert_initial_z_to_feh: bool,
    z_ref: float,
) -> tuple[np.ndarray, str]:
    """Return the metallicity coordinate to plot in the top panel."""
    if metallicity_coordinate is not None:
        return (
            np.full(n, metallicity_coordinate, dtype=float),
            f"constant command-line metallicity = {metallicity_coordinate:g}",
        )

    direct = first_existing_column(
        md,
        (
            "metallicity",
            "Metallicity",
            "feh",
            "FeH",
            "FEH",
            "mh",
            "MH",
            "log_metallicity",
        ),
    )
    if direct is not None:
        name, values = direct
        return values, f"history column {name}"

    zcol = first_existing_column(md, ("initial_z", "Initial_z", "surface_z", "Surface_z", "Z"))
    if zcol is not None:
        name, zvals = zcol
        if convert_initial_z_to_feh:
            values = np.log10(np.maximum(zvals, 1.0e-99) / z_ref)
            return values, f"log10({name}/{z_ref:g})"
        return zvals, f"history column {name}"

    initial_z = None
    if inlist_file is not None:
        initial_z = parse_inlist_value(inlist_file, "initial_z")

    if initial_z is None:
        raise KeyError(
            "No metallicity column found in history.data and no initial_z found in "
            "the supplied inlist. Use --metallicity-coordinate VALUE."
        )

    if convert_initial_z_to_feh:
        value = math.log10(initial_z / z_ref)
        return (
            np.full(n, value, dtype=float),
            f"constant log10(initial_z/{z_ref:g}) = {value:g}",
        )

    return (
        np.full(n, initial_z, dtype=float),
        f"constant initial_z from inlist = {initial_z:g}",
    )


def normalized_nearest_indices(
    path_teff: np.ndarray,
    path_logg: np.ndarray,
    path_meta: np.ndarray,
    grid_teff: np.ndarray,
    grid_logg: np.ndarray,
    grid_meta: np.ndarray,
    chunk_size: int = 256,
) -> np.ndarray:
    """
    Nearest grid node using the same normalized Euclidean idea as
    compute_interp_radius in bolometric.f90.
    """
    eps = 1.0e-12

    t_min, t_max = np.nanmin(grid_teff), np.nanmax(grid_teff)
    g_min, g_max = np.nanmin(grid_logg), np.nanmax(grid_logg)
    m_min, m_max = np.nanmin(grid_meta), np.nanmax(grid_meta)

    t_range = max(t_max - t_min, eps)
    g_range = max(g_max - g_min, eps)
    m_range = max(m_max - m_min, eps)

    grid_norm = np.column_stack(
        (
            (grid_teff - t_min) / t_range,
            (grid_logg - g_min) / g_range,
            (grid_meta - m_min) / m_range,
        )
    )

    path_norm = np.column_stack(
        (
            (path_teff - t_min) / t_range,
            (path_logg - g_min) / g_range,
            (path_meta - m_min) / m_range,
        )
    )

    nearest = np.empty(len(path_norm), dtype=int)

    for lo in range(0, len(path_norm), chunk_size):
        hi = min(lo + chunk_size, len(path_norm))
        diff = path_norm[lo:hi, None, :] - grid_norm[None, :, :]
        dist2 = np.sum(diff * diff, axis=2)
        nearest[lo:hi] = np.argmin(dist2, axis=1)

    return nearest


def set_3d_point(artist, x: float, y: float, z: float) -> None:
    artist._offsets3d = ([x], [y], [z])


def make_colors_grid_movie(
    history_file: str | Path = "../LOGS/history.data",
    lookup_file: str | Path = (
        "$MESA_DIR/data/colors_data/stellar_models/"
        "Kurucz2003all__alpha_04/lookup_table.csv"
    ),
    inlist_file: str | Path | None = "../inlist_colors",
    outfile: str | Path = "colors_grid_cmd_interp_rad.mp4",
    fps: int = 30,
    total_seconds: int = 20,
    mesa_view_stride: int = 5,
    metallicity_coordinate: float | None = None,
    convert_initial_z_to_feh: bool = False,
    z_ref: float = 0.02,
    log_teff_axis: bool = False,
):
    history_file = Path(history_file)
    lookup_file = Path(str(lookup_file).replace("$MESA_DIR", str(Path.home() / "MESA/mesa")))
    if "$MESA_DIR" in str(lookup_file):
        import os

        mesa_dir = os.environ.get("MESA_DIR")
        if mesa_dir:
            lookup_file = Path(str(lookup_file).replace("$MESA_DIR", mesa_dir))

    if inlist_file is not None:
        inlist_file = Path(inlist_file)

    md_raw = mr.MesaData(str(history_file))
    md = MesaView(md_raw, mesa_view_stride) if mesa_view_stride > 1 else md_raw

    all_cols, filter_columns = read_header_columns(str(history_file))
    hr_color, hr_mag, hr_xlabel, hr_ylabel, color_index = setup_hr_diagram_params(
        md, filter_columns
    )

    teff = get_teff(md)
    logg = get_mesa_column(md, "log_g")
    interp_rad = get_mesa_column(md, "Interp_rad")

    n = min(len(teff), len(logg), len(interp_rad), len(hr_color), len(hr_mag))
    teff = teff[:n]
    logg = logg[:n]
    interp_rad = interp_rad[:n]
    hr_color = np.asarray(hr_color, dtype=float)[:n]
    hr_mag = np.asarray(hr_mag, dtype=float)[:n]

    meta, meta_source = get_metallicity_track(
        md,
        n,
        inlist_file=inlist_file,
        metallicity_coordinate=metallicity_coordinate,
        convert_initial_z_to_feh=convert_initial_z_to_feh,
        z_ref=z_ref,
    )
    meta = np.asarray(meta, dtype=float)[:n]

    grid_teff, grid_logg, grid_meta = read_lookup_table(lookup_file)

    finite = (
        np.isfinite(teff)
        & np.isfinite(logg)
        & np.isfinite(meta)
        & np.isfinite(interp_rad)
        & np.isfinite(hr_color)
        & np.isfinite(hr_mag)
    )
    teff = teff[finite]
    logg = logg[finite]
    meta = meta[finite]
    interp_rad = interp_rad[finite]
    hr_color = hr_color[finite]
    hr_mag = hr_mag[finite]

    if len(teff) == 0:
        raise ValueError("No finite history points available for plotting.")

    nearest = normalized_nearest_indices(teff, logg, meta, grid_teff, grid_logg, grid_meta)

    if log_teff_axis:
        x_grid = np.log10(grid_teff)
        x_path = np.log10(teff)
        x_label = r"$\log_{10}(T_{\rm eff}/{\rm K})$"
    else:
        x_grid = grid_teff
        x_path = teff
        x_label = r"$T_{\rm eff}$ [K]"

    nframes = max(2, int(fps * total_seconds))
    frame_to_i = np.linspace(0, len(teff) - 1, nframes).astype(int)

    fig = plt.figure(figsize=(10, 10))
    ax_grid = fig.add_subplot(2, 1, 1, projection="3d")
    ax_cmd = fig.add_subplot(2, 1, 2)

    # Top: real atmosphere grid nodes and stellar path.
    ax_grid.scatter(
        x_grid,
        grid_logg,
        grid_meta,
        marker="s",
        s=5,
        alpha=0.16,
        linewidths=0,
        label="real atmosphere SED nodes",
    )

    full_path, = ax_grid.plot(
        x_path,
        logg,
        meta,
        lw=1.0,
        alpha=0.35,
        label="stellar path",
    )

    prog_path, = ax_grid.plot([], [], [], lw=2.2, label="path so far")
    current_grid = ax_grid.scatter([], [], [], marker="o", s=65, depthshade=True)
    nearest_grid = ax_grid.scatter([], [], [], marker="s", s=75, depthshade=True)
    radius_line, = ax_grid.plot([], [], [], lw=1.8, ls="--", alpha=0.9)

    ax_grid.set_xlabel(x_label)
    ax_grid.set_ylabel(r"$\log g$")
    ax_grid.set_zlabel("metallicity coordinate")

    if not log_teff_axis:
        ax_grid.set_xlim(np.nanmax(x_grid), np.nanmin(x_grid))
    else:
        ax_grid.set_xlim(np.nanmax(x_grid), np.nanmin(x_grid))

    # Bottom: CMD coloured by interpolation radius.
    cmd_scatter = ax_cmd.scatter(
        hr_color,
        hr_mag,
        c=interp_rad,
        s=9,
        alpha=0.85,
    )
    cmd_path, = ax_cmd.plot(hr_color, hr_mag, lw=0.8, alpha=0.25)
    cmd_current = ax_cmd.scatter(
        [],
        [],
        c=[],
        s=85,
        edgecolors="black",
        linewidths=0.8,
        zorder=5,
    )

    cbar = fig.colorbar(cmd_scatter, ax=ax_cmd)
    cbar.set_label("Interp_rad: normalized distance to nearest real SED node")

    ax_cmd.set_xlabel(hr_xlabel)
    ax_cmd.set_ylabel(hr_ylabel)
    ax_cmd.invert_yaxis()

    info = fig.text(
        0.02,
        0.015,
        "",
        ha="left",
        va="bottom",
        fontsize=9,
    )


    fig.subplots_adjust(left=0.08, right=0.92, bottom=0.08, top=0.94, hspace=0.22)

    def update(frame: int):
        i = int(frame_to_i[frame])
        j = int(nearest[i])

        # Top path/progressive point.
        prog_path.set_data_3d(x_path[: i + 1], logg[: i + 1], meta[: i + 1])

        set_3d_point(current_grid, x_path[i], logg[i], meta[i])
        set_3d_point(nearest_grid, x_grid[j], grid_logg[j], grid_meta[j])

        radius_line.set_data_3d(
            [x_path[i], x_grid[j]],
            [logg[i], grid_logg[j]],
            [meta[i], grid_meta[j]],
        )

        # A slow view drift makes the 3D geometry easier to read.
        ax_grid.view_init(elev=45, azim=35 + 80.0 * frame / max(nframes - 1, 1))

        # Bottom current point.
        cmd_current.set_offsets([[hr_color[i], hr_mag[i]]])
        cmd_current.set_array(np.asarray([interp_rad[i]]))
        cmd_current.set_clim(cmd_scatter.get_clim())

        info.set_text(
            f"history: {history_file}    "
            f"grid: {lookup_file}    "
            f"metallicity source: {meta_source}    "
            f"Teff={teff[i]:.0f} K, logg={logg[i]:.3f}, "
            f"metallicity={meta[i]:.4g}, Interp_rad={interp_rad[i]:.4g}"
        )

        return prog_path, current_grid, nearest_grid, radius_line, cmd_current, info

    anim = FuncAnimation(
        fig,
        update,
        frames=nframes,
        interval=1000 / fps,
        blit=False,
    )

    writer = FFMpegWriter(fps=fps, bitrate=3500)
    anim.save(str(outfile), writer=writer)
    plt.close(fig)

    print(f"Wrote {outfile}")
    print(f"Metallicity source: {meta_source}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Make a two-panel colors-module interpolation movie."
    )
    parser.add_argument("--history", default="../LOGS/history.data")
    parser.add_argument(
        "--lookup",
        default="$MESA_DIR/data/colors_data/stellar_models/Kurucz2003all__alpha_04/lookup_table.csv",
        help="Path to atmosphere lookup_table.csv",
    )
    parser.add_argument(
        "--inlist",
        default="../inlist_colors",
        help="Inlist to parse for initial_z if history.data has no metallicity column",
    )
    parser.add_argument("--outfile", default="colors_grid_cmd_interp_rad.mp4")
    parser.add_argument("--fps", type=int, default=30)
    parser.add_argument("--seconds", type=int, default=20)
    parser.add_argument(
        "--mesa-view-stride",
        type=int,
        default=5,
        help="Downsample history through MesaView, matching the existing movie script.",
    )
    parser.add_argument(
        "--metallicity-coordinate",
        type=float,
        default=None,
        help="Explicit constant metallicity coordinate for the stellar path.",
    )
    parser.add_argument(
        "--convert-initial-z-to-feh",
        action="store_true",
        help="If using initial_z, plot log10(initial_z / z_ref) instead of initial_z.",
    )
    parser.add_argument(
        "--z-ref",
        type=float,
        default=0.02,
        help="Reference Z used with --convert-initial-z-to-feh.",
    )
    parser.add_argument(
        "--log-teff-axis",
        action="store_true",
        help="Plot log10(Teff) in the 3D panel instead of Teff.",
    )
    return parser


def main():
    args = build_parser().parse_args()

    # Preserve the old behaviour if the user runs from inside python_helpers:
    # find ../LOGS/history.data by default.
    history = args.history
    if history == "../LOGS/history.data":
        matches = glob.glob("../LOGS/history.data")
        if matches:
            history = matches[0]

    make_colors_grid_movie(
        history_file=history,
        lookup_file=args.lookup,
        inlist_file=args.inlist,
        outfile=args.outfile,
        fps=args.fps,
        total_seconds=args.seconds,
        mesa_view_stride=args.mesa_view_stride,
        metallicity_coordinate=args.metallicity_coordinate,
        convert_initial_z_to_feh=args.convert_initial_z_to_feh,
        z_ref=args.z_ref,
        log_teff_axis=args.log_teff_axis,
    )


if __name__ == "__main__":
    main()
