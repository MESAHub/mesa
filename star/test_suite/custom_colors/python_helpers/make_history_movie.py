#!/usr/bin/env python3
# Make a movie that visually matches the live HISTORY_check viewer.
# One frame per row; each frame shows the full chain up to that row.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from HISTORY_check import HistoryChecker  # uses static_HISTORY_check under the hood


do_tqdm = True
try:
    from tqdm import tqdm
except ImportError:
    do_tqdm = False


# Reuse the live viewer so formatting/logic stays identical

OUT = "history.mp4"
FPS = 24
DPI = 150

# Build once from the user's checker and pull arrays
checker = HistoryChecker(history_file="../LOGS/history.data", refresh_interval=0.1)
checker.update_data()  # read md + compute columns/colors
checker.setup_plot()  # set labels, axes styling, grids

age = np.asarray(checker.Star_Age, float)
Teff = np.asarray(checker.Teff, float)
Log_L = np.asarray(checker.Log_L, float)
hr_x = np.asarray(checker.hr_color, float)
hr_y = np.asarray(checker.hr_mag, float)
phase_colors = (
    np.asarray(checker.phase_colors, object) if len(checker.phase_colors) else None
)
filters = list(checker.filter_columns)
fig, axes = checker.fig, checker.axes

if age.size == 0:
    raise RuntimeError("history.data has zero rows")


# Helper: padded limits
def pad(lo, hi, frac=0.05):
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        return lo - 1, hi + 1
    span = hi - lo
    return lo - frac * span, hi + frac * span


# Fix axis limits once so the camera doesn't jump
age_lo, age_hi = pad(np.nanmin(age), np.nanmax(age))
Teff_lo, Teff_hi = pad(np.nanmin(Teff), np.nanmax(Teff))
LogL_lo, LogL_hi = pad(np.nanmin(Log_L), np.nanmax(Log_L))
hrx_lo, hrx_hi = pad(np.nanmin(hr_x), np.nanmax(hr_x))
hry_lo, hry_hi = pad(np.nanmin(hr_y), np.nanmax(hr_y))
ci_lo, ci_hi = pad(np.nanmin(checker.color_index), np.nanmax(checker.color_index))

# Preload all magnitude series + their assigned colors from the live viewer
mag_series = []
for i, f in enumerate(filters):
    y = np.asarray(getattr(checker.md, f), float)
    if y.size != age.size:
        continue
    color = checker.filter_colors[i] if i < len(checker.filter_colors) else "black"
    mag_series.append((f, y, color))

if mag_series:
    m_all = np.vstack([s for _, s, _ in mag_series])
    mag_lo, mag_hi = pad(np.nanmin(m_all), np.nanmax(m_all))
else:
    mag_lo, mag_hi = 0.0, 1.0

writer = FFMpegWriter(fps=FPS, codec="libx264", extra_args=["-pix_fmt", "yuv420p"])

with writer.saving(fig, OUT, dpi=DPI):
    N = age.size

    if do_tqdm:
        iterator = tqdm(range(N))
    else:
        iterator = range(N)

    for i in iterator:
        sl = slice(0, i + 1)

        # Clear and re-apply the exact same formatting
        for ax in axes.flatten():
            ax.cla()
        checker.setup_plot()

        # Colors for this frame (phase legend if available; else colormap colors are already in phase_colors)
        C = (
            phase_colors[sl]
            if phase_colors is not None and len(phase_colors) == N
            else None
        )

        # === Top-left: CMD/HR (color vs mag) ===
        axes[0, 0].plot(
            hr_x[sl], hr_y[sl], marker=",", linestyle="-", color="k", alpha=0.2
        )
        if C is not None:
            axes[0, 0].scatter(
                hr_x[sl], hr_y[sl], c=C, s=15, alpha=0.9, edgecolors="none"
            )
        else:
            axes[0, 0].scatter(hr_x[sl], hr_y[sl], s=15, alpha=0.9, edgecolors="none")
        axes[0, 0].set_xlim(hrx_lo, hrx_hi)
        axes[0, 0].set_ylim(hry_hi, hry_lo)  # invert mags

        # === Top-right: Teff vs Log L ===
        axes[0, 1].plot(
            Teff[sl], Log_L[sl], marker=",", linestyle="-", color="k", alpha=0.2
        )
        if C is not None:
            axes[0, 1].scatter(
                Teff[sl], Log_L[sl], c=C, s=15, alpha=0.9, edgecolors="none"
            )
        else:
            axes[0, 1].scatter(Teff[sl], Log_L[sl], s=15, alpha=0.9, edgecolors="none")
        axes[0, 1].set_xlim(Teff_hi, Teff_lo)  # Teff inverted x
        axes[0, 1].set_ylim(LogL_lo, LogL_hi)

        # === Bottom-left: Age vs color index ===
        axes[1, 0].plot(
            age[sl],
            checker.color_index[sl],
            marker=",",
            linestyle="-",
            color="k",
            alpha=0.2,
        )
        if C is not None:
            axes[1, 0].scatter(
                age[sl],
                checker.color_index[sl],
                c=C,
                s=15,
                alpha=1.0,
                edgecolors="none",
            )
        else:
            axes[1, 0].scatter(
                age[sl], checker.color_index[sl], s=15, alpha=1.0, edgecolors="none"
            )
        axes[1, 0].set_xlim(age_lo, age_hi)
        axes[1, 0].set_ylim(ci_lo, ci_hi)

        # === Bottom-right: Age vs all filter mags (full chain) ===
        for label, series, color in mag_series:
            axes[1, 1].plot(
                age[sl],
                series[sl],
                marker="o",
                linestyle="-",
                markersize=3,
                alpha=0.8,
                label=label,
                color=color,
            )
        axes[1, 1].set_xlim(age_lo, age_hi)
        axes[1, 1].set_ylim(mag_hi, mag_lo)  # invert mags
        if mag_series:
            axes[1, 1].legend()

        # Phase legend exactly like the live viewer (only when real phases are present)
        if getattr(checker, "has_phase", False):
            legend_elements = checker.create_phase_legend()
            if len(legend_elements) > 0:
                # keep the filter legend on bottom-right and add a centered phase legend
                axes[1, 1].legend()
                n_phases = len(legend_elements)
                ncol = max(1, (n_phases + 1) // 1)
                fig.legend(
                    handles=legend_elements,
                    loc="upper center",
                    bbox_to_anchor=(0.5, 0.53),
                    ncol=ncol,
                    title_fontsize=10,
                    fontsize=10,
                    frameon=True,
                    fancybox=True,
                    shadow=True,
                )

        plt.subplots_adjust(top=0.92)
        writer.grab_frame()

print(f"Wrote {age.size} frames -> {OUT}")
