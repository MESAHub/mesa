#!/usr/bin/env python

import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.backends.backend_pdf import PdfPages


DATA = Path("data/eos_composition_partials.csv")
CONTOUR_DATA = Path("data/eos_composition_partials_contours.csv")
DFRIDR_DATA = Path("data/eos_composition_partials_dfridr_contours.csv")
FIGURES = Path("figures")
SUMMARY = Path("data/eos_composition_partials_summary.txt")
COMBINED = FIGURES / "all_composition_partial_contours.pdf"
NO_COVERAGE_COLOR = "#ffffff"

plt.style.use("../mesa_eos_regions.mplstyle")
mpl.rcParams.update({
    "figure.max_open_warning": 0,
    "axes.titlepad": 4.0,
    "legend.handlelength": 1.6,
})

KEY_ROWS = [
    "chiX",
    "mu",
    "lnPgas",
    "lnE",
    "lnS",
    "lnfree_e",
    "eta",
    "grad_ad",
    "chiRho",
    "chiT",
    "Cp",
    "Cv",
    "gamma1",
    "gamma3",
]

KEY_DIRECTIONS = ["he4", "c12", "o16"]
DFRIDR_ROWS = ["Pgas", "mu", "lnE"]
DFRIDR_DIRECTIONS = ["he4", "c12", "o16"]

FRACTION_COLS = [
    "frac_OPAL_SCVH",
    "frac_HELM",
    "frac_Skye",
    "frac_PC",
    "frac_FreeEOS",
    "frac_CMS",
    "frac_ideal",
]

COMPONENTS = [
    ("HELM", "frac_HELM"),
    ("OPAL/SCVH", "frac_OPAL_SCVH"),
    ("FreeEOS", "frac_FreeEOS"),
    ("PC", "frac_PC"),
    ("Skye", "frac_Skye"),
    ("CMS", "frac_CMS"),
    ("ideal", "frac_ideal"),
]

ROW_LABELS = {
    "chiX": r"\chi_{X_i}",
    "Pgas": r"$P_{\rm gas}$",
    "mu": r"$\mu$",
    "lnPgas": r"$\ln P_{\rm gas}$",
    "lnE": r"$\ln E$",
    "lnS": r"$\ln S$",
    "lnfree_e": r"$\ln n_{e,{\rm free}}$",
    "eta": r"$\eta$",
    "grad_ad": r"$\nabla_{\rm ad}$",
    "chiRho": r"$\chi_\rho$",
    "chiT": r"$\chi_T$",
    "Cp": r"$C_P$",
    "Cv": r"$C_V$",
    "gamma1": r"$\Gamma_1$",
    "gamma3": r"$\Gamma_3$",
}

DIR_LABELS = {
    "h1": r"{\rm h1}",
    "he4": r"{\rm he4}",
    "c12": r"{\rm c12}",
    "o16": r"{\rm o16}",
    "mg24": r"{\rm mg24}",
    "fe56": r"{\rm fe56}",
}

SWEEP_LABELS = {
    "cool_T_sweep": r"cool $T$ sweep",
    "dense_T_sweep": r"dense $T$ sweep",
    "hot_T_sweep": r"hot $T$ sweep",
    "rho_sweep": r"$\rho$ sweep",
}

LINE_COLORS = {
    "h1": "#1f77b4",
    "c12": "#d62728",
    "o16": "#2ca02c",
    "mg24": "#9467bd",
    "fe56": "#ff7f0e",
}


def as_float(value):
    return float(value.replace("D", "E"))


def safe_name(value):
    return "".join(ch if ch.isalnum() else "_" for ch in value).strip("_")


def load_rows(path, has_fd):
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        rows = []
        numeric_keys = [
            "logT",
            "logRho",
            *FRACTION_COLS,
        ]
        for key in ("analytic", "raw_partial", "sink_projected_partial"):
            if key in reader.fieldnames:
                numeric_keys.append(key)
        if has_fd:
            numeric_keys.extend([
                "finite_difference",
                "finite_difference_error",
                "abs_error",
                "rel_error",
            ])
        for row in reader:
            for key in numeric_keys:
                row[key] = as_float(row[key])
            if "point" in row:
                row["point"] = int(row["point"])
            row["row"] = int(row["row"])
            rows.append(row)
    return rows


def clean_figures():
    FIGURES.mkdir(parents=True, exist_ok=True)
    for pattern in ("*.pdf", "*.png", "*.svg"):
        for old in FIGURES.glob(pattern):
            old.unlink()


def write_summary(rows):
    stats = defaultdict(lambda: {"max_rel": 0.0, "max_abs": 0.0, "max_fd_err": 0.0, "count": 0})
    for row in rows:
        key = (row["sweep"], row["row_name"])
        stats[key]["max_rel"] = max(stats[key]["max_rel"], row["rel_error"])
        stats[key]["max_abs"] = max(stats[key]["max_abs"], row["abs_error"])
        stats[key]["max_fd_err"] = max(stats[key]["max_fd_err"], row["finite_difference_error"])
        stats[key]["count"] += 1

    with SUMMARY.open("w") as f:
        f.write("EOS composition partial validation summary\n")
        f.write("finite_difference uses constrained Ridders extrapolation with he4 as sink\n\n")
        f.write(f"{'sweep':18s} {'row':20s} {'max_rel':>14s} {'max_abs':>14s} {'max_fd_err':>14s} {'n':>8s}\n")
        for (sweep, row_name), item in sorted(stats.items()):
            f.write(
                f"{sweep:18s} {row_name:20s} "
                f"{item['max_rel']:14.6e} {item['max_abs']:14.6e} "
                f"{item['max_fd_err']:14.6e} {item['count']:8d}\n"
            )
    return stats


def save_figure(fig, name, pdf):
    pdf.savefig(fig)
    fig.savefig(FIGURES / f"{name}.pdf")
    fig.savefig(FIGURES / f"{name}.png", dpi=220)
    plt.close(fig)


def row_math(row_name):
    return ROW_LABELS[row_name].strip("$")


def raw_partial_label(row_name):
    if row_name == "chiX":
        return r"$\chi_{X_i}=(\partial\ln P/\partial X_i)_{\rho,T}$"
    q = row_math(row_name)
    return rf"$\left(\partial {q}/\partial X_i\right)_{{\rho,T}}$"


def contour_title(row_name):
    if row_name == "chiX":
        return rf"Brunt pressure composition partial: {raw_partial_label(row_name)}"
    return rf"raw composition partial: {raw_partial_label(row_name)}"


def direction_partial_label(direction):
    return rf"$\partial/\partial X_{{{DIR_LABELS[direction]}}}$"


def sink_direction_label(direction, sink):
    return rf"$D_{{X_{{{DIR_LABELS[direction]}}}}}^{{X_{{{DIR_LABELS[sink]}}}}}$"


def xkey_for_sweep(sweep):
    return "logRho" if sweep == "rho_sweep" else "logT"


def symlog_linthresh(values):
    vals = np.asarray([abs(v) for v in values if np.isfinite(v) and v != 0.0])
    if vals.size == 0:
        return 1e-12
    return max(1e-14, np.nanpercentile(vals, 5)*0.2)


def plot_eos_fractions(rows, pdf):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["sweep"]].append(row)

    fig, axes = plt.subplots(2, 2, figsize=(6.7, 4.8), sharey=True, constrained_layout=True)
    axes = axes.ravel()
    for ax, sweep in zip(axes, sorted(grouped)):
        points = {}
        for row in grouped[sweep]:
            points.setdefault(row["point"], row)
        ordered = [points[key] for key in sorted(points)]
        xkey = xkey_for_sweep(sweep)
        x = np.array([row[xkey] for row in ordered])
        for col in FRACTION_COLS:
            y = np.array([row[col] for row in ordered])
            ax.plot(x, y, lw=1.0, label=col.replace("frac_", ""))
        ax.set_title(SWEEP_LABELS.get(sweep, sweep))
        ax.set_xlabel(r"$\log \rho$" if xkey == "logRho" else r"$\log T$")
        ax.set_ylim(-0.05, 1.05)
    axes[0].set_ylabel("EOS fraction")
    axes[2].set_ylabel("EOS fraction")
    axes[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=6)
    save_figure(fig, "eos_fractions", pdf)


def plot_line_comparisons(rows, pdf):
    grouped = defaultdict(list)
    for row in rows:
        if row["row_name"] in KEY_ROWS and row["direction"] in KEY_DIRECTIONS:
            grouped[(row["sweep"], row["row_name"], row["direction"])].append(row)

    for sweep in sorted({key[0] for key in grouped}):
        xkey = xkey_for_sweep(sweep)
        xlabel = r"$\log(\rho/{\rm g\,cm^{-3}})$" if xkey == "logRho" else r"$\log(T/{\rm K})$"
        for row_name in KEY_ROWS:
            keys = [(sweep, row_name, direction) for direction in KEY_DIRECTIONS]
            if not any(key in grouped for key in keys):
                continue

            values = []
            for key in keys:
                for row in grouped.get(key, []):
                    values.extend([row["analytic"], row["finite_difference"]])
            linthresh = symlog_linthresh(values)

            fig, (ax, axerr) = plt.subplots(
                2, 1, sharex=True, figsize=(4.8, 3.6),
                gridspec_kw={"height_ratios": [2.2, 1.0]},
                constrained_layout=True,
            )
            for key in keys:
                ordered = sorted(grouped.get(key, []), key=lambda item: item["point"])
                if not ordered:
                    continue
                direction = key[2]
                color = LINE_COLORS[direction]
                x = np.array([item[xkey] for item in ordered])
                analytic = np.array([item["analytic"] for item in ordered])
                fd = np.array([item["finite_difference"] for item in ordered])
                rel = np.array([max(item["rel_error"], 1e-18) for item in ordered])
                ax.plot(x, analytic, color=color, lw=1.0, label=DIR_LABELS[direction])
                ax.plot(x, fd, color=color, lw=0.9, ls=":")
                axerr.semilogy(x, rel, color=color, lw=0.9)

            ax.set_title(rf"{SWEEP_LABELS.get(sweep, sweep)}: {ROW_LABELS[row_name]}")
            ax.set_ylabel(r"$\partial_{X_i-X_{\rm he4}} q$")
            ax.set_yscale("symlog", linthresh=linthresh)
            ax.legend(ncol=3, loc="best", fontsize=6)
            axerr.set_xlabel(xlabel)
            axerr.set_ylabel("rel. err.")
            axerr.set_ylim(1e-16, 2.0)
            save_figure(fig, f"line_{safe_name(sweep)}_{safe_name(row_name)}", pdf)


def grid_from_rows(rows, value_key):
    xs = np.array(sorted({row["logRho"] for row in rows}))
    ys = np.array(sorted({row["logT"] for row in rows}))
    x_index = {value: i for i, value in enumerate(xs)}
    y_index = {value: i for i, value in enumerate(ys)}
    z = np.full((len(ys), len(xs)), np.nan)
    for row in rows:
        z[y_index[row["logT"]], x_index[row["logRho"]]] = row[value_key]
    return xs, ys, z


def masked_grid(z):
    return np.ma.masked_invalid(z)


def partial_cmap(name):
    cmap = mpl.colormaps[name].copy()
    cmap.set_bad(NO_COVERAGE_COLOR, alpha=1.0)
    return cmap


def plot_eos_region_contour(contour_rows, pdf):
    points = {}
    for row in contour_rows:
        key = (row["logT"], row["logRho"])
        if key not in points:
            points[key] = row

    rows = list(points.values())
    xs = np.array(sorted({row["logRho"] for row in rows}))
    ys = np.array(sorted({row["logT"] for row in rows}))
    x_index = {value: i for i, value in enumerate(xs)}
    y_index = {value: i for i, value in enumerate(ys)}
    z = np.full((len(ys), len(xs)), np.nan)
    for row in rows:
        fractions = np.array([row[col] for _, col in COMPONENTS])
        imax = int(np.nanargmax(fractions))
        code = imax + 1
        if fractions[imax] < 0.999:
            code = 0
        z[y_index[row["logT"]], x_index[row["logRho"]]] = code

    my_colors = np.array(mpl.colormaps["Accent"].colors)
    tmp = my_colors[4].copy()
    my_colors[4] = my_colors[3]
    my_colors[5] = tmp
    cmap = colors.ListedColormap(my_colors)
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 5.5, 7.5]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    cmap.set_under("white")
    cmap.set_bad(NO_COVERAGE_COLOR, alpha=1.0)

    fig, ax = plt.subplots(figsize=(4.8, 3.75))
    mesh = ax.pcolormesh(xs, ys, masked_grid(z), shading="nearest", cmap=cmap, norm=norm, rasterized=True)
    ax.set_title(r"MESA EOS Regions ($X=0.7$, $Z=0.02$)")
    ax.set_xlabel(r"$\log(\rho/{\rm g\,cm^{-3}})$")
    ax.set_ylabel(r"$\log(T/{\rm K})$")
    ax.set_xlim(xs.min(), xs.max())
    ax.set_ylim(ys.min(), ys.max())
    ax.text(4.0, 2.3, "no coverage", fontsize="small")
    cbar = fig.colorbar(mesh, ax=ax, ticks=[0, 1, 2, 3, 4.5, 6.5], pad=0.02)
    cbar.set_label("")
    cbar.ax.set_yticklabels(["blend", "HELM", "OPAL/SCVH", "FreeEOS", "Skye", "ideal"])
    cbar.ax.tick_params(labelsize="x-small")
    cbar.ax.minorticks_off()
    fig.tight_layout(pad=0.35)
    save_figure(fig, "contour_eos_regions", pdf)


def contour_norm(all_values):
    vals = np.asarray([v for v in all_values if np.isfinite(v)])
    if vals.size == 0:
        return colors.Normalize(vmin=-1.0, vmax=1.0), partial_cmap("coolwarm")
    vmax = np.nanpercentile(np.abs(vals), 98.0)
    vmax = max(vmax, np.nanmax(np.abs(vals))*0.2, 1e-20)
    linthresh = max(vmax*1e-4, 1e-14)
    return colors.SymLogNorm(linthresh=linthresh, vmin=-vmax, vmax=vmax, base=10), partial_cmap("coolwarm")


def compact_symlog_ticks(norm):
    if not isinstance(norm, colors.SymLogNorm):
        return None

    vmax = max(abs(norm.vmin), abs(norm.vmax))
    if vmax <= 0:
        return [0.0]

    emin = math.ceil(math.log10(norm.linthresh))
    emax = math.floor(math.log10(vmax))
    if emin > emax:
        return [-vmax, 0.0, vmax]

    positives = [10.0**e for e in range(emin, emax + 1)]
    if len(positives) > 2:
        positives = [positives[0], positives[-1]]
    ticks = [-value for value in reversed(positives)] + [0.0] + positives
    return [tick for tick in ticks if norm.vmin <= tick <= norm.vmax]


def plot_partial_contours(contour_rows, pdf):
    grouped = defaultdict(list)
    for row in contour_rows:
        if row["row_name"] in KEY_ROWS and row["direction"] in KEY_DIRECTIONS:
            grouped[(row["row_name"], row["direction"])].append(row)

    for row_name in KEY_ROWS:
        all_values = []
        for direction in KEY_DIRECTIONS:
            all_values.extend(row["raw_partial"] for row in grouped.get((row_name, direction), []))
        norm, cmap = contour_norm(all_values)

        fig = plt.figure(figsize=(6.35, 2.55))
        gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.055], wspace=0.18)
        axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
        cax = fig.add_subplot(gs[0, 3])
        mesh = None
        for ax, direction in zip(axes, KEY_DIRECTIONS):
            records = grouped.get((row_name, direction), [])
            xs, ys, z = grid_from_rows(records, "raw_partial")
            mesh = ax.pcolormesh(xs, ys, masked_grid(z), shading="auto", cmap=cmap, norm=norm, rasterized=True)
            ax.set_title(direction_partial_label(direction), fontsize=8)
            ax.set_xlabel(r"$\log(\rho/{\rm g\,cm^{-3}})$")
            if ax is axes[0]:
                ax.set_ylabel(r"$\log(T/{\rm K})$")
            else:
                ax.tick_params(labelleft=False)
            ax.tick_params(labelsize=7, pad=1)
        fig.suptitle(contour_title(row_name), y=0.985, fontsize=9.5)
        if mesh is not None:
            cbar = fig.colorbar(mesh, cax=cax)
            ticks = compact_symlog_ticks(norm)
            if ticks is not None:
                cbar.set_ticks(ticks)
            cbar.set_label(raw_partial_label(row_name), fontsize=8, labelpad=2)
            cbar.ax.tick_params(labelsize=7, pad=1)
        fig.subplots_adjust(left=0.07, right=0.945, bottom=0.18, top=0.81, wspace=0.22)
        save_figure(fig, f"contour_partials_{safe_name(row_name)}", pdf)


def plot_dfridr_contours(dfridr_rows, pdf):
    grouped = defaultdict(list)
    for row in dfridr_rows:
        if row["row_name"] in DFRIDR_ROWS and row["direction"] in DFRIDR_DIRECTIONS:
            grouped[(row["row_name"], row["direction"])].append(row)

    cmap = partial_cmap("viridis")
    norm = colors.Normalize(vmin=-12.0, vmax=0.0)
    for row_name in DFRIDR_ROWS:
        fig = plt.figure(figsize=(6.35, 2.55))
        gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.055], wspace=0.18)
        axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
        cax = fig.add_subplot(gs[0, 3])
        mesh = None
        for ax, direction in zip(axes, DFRIDR_DIRECTIONS):
            records = grouped.get((row_name, direction), [])
            xs, ys, rel = grid_from_rows(records, "rel_error")
            with np.errstate(divide="ignore", invalid="ignore"):
                z = np.log10(np.maximum(rel, 1e-18))
            mesh = ax.pcolormesh(xs, ys, masked_grid(z), shading="auto", cmap=cmap, norm=norm, rasterized=True)
            sink = records[0]["sink"] if records else "h1"
            ax.set_title(sink_direction_label(direction, sink), fontsize=8)
            ax.set_xlabel(r"$\log(\rho/{\rm g\,cm^{-3}})$")
            if ax is axes[0]:
                ax.set_ylabel(r"$\log(T/{\rm K})$")
            else:
                ax.tick_params(labelleft=False)
            ax.tick_params(labelsize=7, pad=1)
        fig.suptitle(rf"DFRIDR composition partial check: {ROW_LABELS[row_name]}", y=0.985, fontsize=9.5)
        if mesh is not None:
            cbar = fig.colorbar(mesh, cax=cax, extend="max")
            cbar.set_label(r"$\log_{10}|A-F|/\max(|A|,|F|)$", fontsize=8, labelpad=2)
            cbar.ax.tick_params(labelsize=7, pad=1)
        fig.subplots_adjust(left=0.07, right=0.945, bottom=0.18, top=0.81, wspace=0.22)
        save_figure(fig, f"dfridr_relerr_{safe_name(row_name)}", pdf)


def plot_summary_heatmap(stats, pdf):
    sweeps = sorted({key[0] for key in stats})
    row_names = sorted({key[1] for key in stats})
    matrix = np.full((len(row_names), len(sweeps)), np.nan)
    for i, row_name in enumerate(row_names):
        for j, sweep in enumerate(sweeps):
            rel = stats.get((sweep, row_name), {}).get("max_rel", np.nan)
            matrix[i, j] = np.log10(max(rel, 1e-18)) if np.isfinite(rel) else np.nan

    fig, ax = plt.subplots(figsize=(5.5, 5.8), constrained_layout=True)
    image = ax.imshow(matrix, aspect="auto", vmin=-16, vmax=0, cmap="viridis")
    ax.set_xticks(range(len(sweeps)))
    ax.set_xticklabels([SWEEP_LABELS.get(sweep, sweep) for sweep in sweeps], rotation=30, ha="right")
    ax.set_yticks(range(len(row_names)))
    ax.set_yticklabels([ROW_LABELS.get(name, name) for name in row_names])
    ax.set_title(r"max $\log_{10}$ relative error")
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label(r"$\log_{10}\max |A-F|/\max(|A|,|F|)$")
    save_figure(fig, "summary_max_rel_error", pdf)


def main():
    clean_figures()
    rows = load_rows(DATA, has_fd=True)
    contour_rows = load_rows(CONTOUR_DATA, has_fd=False)
    dfridr_rows = load_rows(DFRIDR_DATA, has_fd=True)
    write_summary(rows)

    with PdfPages(COMBINED) as pdf:
        plot_eos_region_contour(contour_rows, pdf)
        plot_partial_contours(contour_rows, pdf)
        plot_dfridr_contours(dfridr_rows, pdf)

    print(f"wrote {SUMMARY}")
    print(f"wrote contour PDF and PNG plots to {FIGURES}")
    print(f"wrote combined contour PDF report to {COMBINED}")


if __name__ == "__main__":
    main()
