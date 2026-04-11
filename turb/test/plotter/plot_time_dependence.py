#!/usr/bin/env python

import sys

sys.dont_write_bytecode = True

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


HERE = Path(__file__).resolve().parent
CSV_PATH = HERE / "time_dependence.csv"
OUTPUT_DIR = HERE / "plots"
STYLE_PATH = HERE / "mesa.mplstyle"
OUTPUT_PATH = OUTPUT_DIR / "time_dependence_conv_vel_panels.pdf"

MODE_ORDER = [
    "mlt",
    "tdc",
    "tdc_with_mlt_corr",
    "tdc_with_arnett_closure",
    "tdc_with_arnett_closure_tdc_ss",
    "tdc_with_arnett_closure_tdc_ss_mlt_corr",
    "tdc_with_acceleration_limit",
]

MODE_LABELS = {
    "mlt": "MLT",
    "tdc": "plain TDC",
    "tdc_with_mlt_corr": "TDC + MLT corr",
    "tdc_with_arnett_closure": "TDC + Arnett (MLT steady state)",
    "tdc_with_arnett_closure_tdc_ss": "TDC + Arnett (TDC steady state)",
    "tdc_with_arnett_closure_tdc_ss_mlt_corr": "TDC + Arnett (TDC steady state, MLT corr)",
    "tdc_with_acceleration_limit": "TDC + acceleration limit",
}

MODE_COLORS = {
    "mlt": plt.cm.tab10(0),
    "tdc": plt.cm.tab10(1),
    "tdc_with_mlt_corr": plt.cm.tab10(2),
    "tdc_with_arnett_closure": plt.cm.tab10(3),
    "tdc_with_arnett_closure_tdc_ss": plt.cm.tab10(4),
    "tdc_with_arnett_closure_tdc_ss_mlt_corr": plt.cm.tab10(5),
    "tdc_with_acceleration_limit": "black",
}

MODE_LINEWIDTHS = {
    "mlt": 5.0,
    "tdc": 5.0,
    "tdc_with_mlt_corr": 4.0,
    "tdc_with_arnett_closure": 2.0,
    "tdc_with_arnett_closure_tdc_ss": 2.0,
    "tdc_with_arnett_closure_tdc_ss_mlt_corr": 2.0,
    "tdc_with_acceleration_limit": 3.0,
}

MODE_LINESTYLES = {
    "mlt": "-",
    "tdc": "-",
    "tdc_with_mlt_corr": "--",
    "tdc_with_arnett_closure": "-",
    "tdc_with_arnett_closure_tdc_ss": "-",
    "tdc_with_arnett_closure_tdc_ss_mlt_corr": "-",
    "tdc_with_acceleration_limit": ":",
}

MODE_ZORDERS = {
    "mlt": 1,
    "tdc": 1,
    "tdc_with_mlt_corr": 2,
    "tdc_with_acceleration_limit": 3,
    "tdc_with_arnett_closure": 4,
    "tdc_with_arnett_closure_tdc_ss": 4,
    "tdc_with_arnett_closure_tdc_ss_mlt_corr": 4,
}

PANEL_GROUPS = [
    "decay_evolution",
    "growth_evolution",
    "envelope_decay_evolution",
    "cburn_growth_evolution",
]

LINE_ALPHA = 0.8


def read_rows():
    groups = {}
    with CSV_PATH.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            group = row["group"]
            mode = row["mode"]
            row["step"] = int(row["step"])
            row["time"] = float(row["time"])
            row["conv_vel"] = float(row["conv_vel"])
            groups.setdefault(group, {}).setdefault(mode, []).append(row)

    for tracks in groups.values():
        for rows in tracks.values():
            rows.sort(key=lambda row: row["step"])

    return groups


if __name__ == "__main__":
    if not CSV_PATH.exists():
        raise SystemExit(f"missing CSV: {CSV_PATH}")

    OUTPUT_DIR.mkdir(exist_ok=True)

    groups = read_rows()

    missing_groups = [group_name for group_name in PANEL_GROUPS if group_name not in groups]
    if missing_groups:
        raise SystemExit(f"missing groups in CSV: {', '.join(missing_groups)}")

    if STYLE_PATH.exists():
        plt.style.use(str(STYLE_PATH))

    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["axes.facecolor"] = "white"
    plt.rcParams["savefig.format"] = "pdf"

    fig, axes = plt.subplots(2, 2, figsize=(24.0, 16.0))
    axes = axes.flatten()

    for ax, group_name in zip(axes, PANEL_GROUPS):
        tracks = groups[group_name]

        conv_vel_values = []
        for mode in MODE_ORDER:
            rows = tracks.get(mode, [])
            if not rows:
                continue

            x = [row["time"] for row in rows]
            y = [row["conv_vel"] for row in rows]
            conv_vel_values.extend(y)

            ax.plot(
                x,
                y,
                label=MODE_LABELS[mode],
                color=MODE_COLORS[mode],
                linewidth=MODE_LINEWIDTHS[mode],
                linestyle=MODE_LINESTYLES[mode],
                alpha=LINE_ALPHA,
                zorder=MODE_ZORDERS[mode],
            )

        ymin = min(conv_vel_values)
        ymax = max(conv_vel_values)
        if ymax == ymin:
            pad = max(1.0, abs(ymax)) * 0.05
        else:
            pad = 0.05 * (ymax - ymin)

        ax.set_title(
            group_name
            .replace("cburn", "carbon burning")
            .replace("_evolution", "")
            .replace("_", " ")
            .title()
        )
        ax.set_xlabel("time [s]")
        ax.set_ylabel("conv_vel [cm/s]")
        ax.set_ylim(ymin - pad, ymax + pad)
        ax.grid(True)

        x_formatter = ScalarFormatter(useMathText=True)
        x_formatter.set_scientific(True)
        x_formatter.set_powerlimits((-3, 3))
        x_formatter.set_useOffset(False)
        ax.xaxis.set_major_formatter(x_formatter)

        y_formatter = ScalarFormatter(useMathText=True)
        y_formatter.set_scientific(True)
        y_formatter.set_powerlimits((-3, 3))
        y_formatter.set_useOffset(False)
        ax.yaxis.set_major_formatter(y_formatter)

    axes[0].legend(loc="upper right")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.97))
    fig.savefig(OUTPUT_PATH)
    plt.close(fig)

    print(f"wrote {OUTPUT_PATH}")
