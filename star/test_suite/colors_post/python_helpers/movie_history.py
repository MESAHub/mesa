#!/usr/bin/env python3
# Make a movie that visually matches the live HISTORY_check viewer.
# One frame per row; each frame shows the full chain up to that row.

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter
from plot_history_live import HistoryChecker  # uses static_HISTORY_check under the hood

do_tqdm = True
try:
    from tqdm import tqdm
except ImportError:
    do_tqdm = False


OUT = "history.mp4"
FPS = 24
DPI = 150

# Build once from the user's checker and pull arrays
checker = HistoryChecker(history_file="../LOGS/history.data", refresh_interval=0.1)
checker.update_data()  # read md + compute columns/colors
checker.setup_plot()   # set labels, axes styling, grids — called ONCE

age = np.asarray(checker.Star_Age, float)
Teff = np.asarray(checker.Teff, float)
Log_L = np.asarray(checker.Log_L, float)
hr_x = np.asarray(checker.hr_color, float)
hr_y = np.asarray(checker.hr_mag, float)
color_index = np.asarray(checker.color_index, float)
phase_colors = (
    np.asarray(checker.phase_colors, object) if len(checker.phase_colors) else None
)
has_phases = phase_colors is not None and len(phase_colors) == age.size
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
ci_lo, ci_hi = pad(np.nanmin(color_index), np.nanmax(color_index))

# Preload all magnitude series + their assigned colors from the live viewer.
# Apply the same physical mask as update_plot so duplicate/non-physical columns
# (e.g. from unused atmosphere grids) are silently skipped.
mag_series = []
seen_filters = set()
for i, f in enumerate(filters):
    if f in seen_filters:
        continue
    y = np.asarray(getattr(checker.md, f), float)
    if y.size != age.size:
        continue
    mask = np.isfinite(y) & np.isfinite(age) & (y < 90.0) & (y > -50.0)
    if not np.any(mask):
        continue
    seen_filters.add(f)
    color = checker.filter_colors[i] if i < len(checker.filter_colors) else "black"
    mag_series.append((f, y, color))

if mag_series:
    m_all = np.vstack([s for _, s, _ in mag_series])
    mag_lo, mag_hi = pad(np.nanmin(m_all), np.nanmax(m_all))
else:
    mag_lo, mag_hi = 0.0, 1.0

# ── Create all artists ONCE with empty data ───────────────────────────────────
# Top-left: CMD/HR
hr_line,   = axes[0, 0].plot([], [], marker=",", linestyle="-", color="k", alpha=0.2)
hr_sc      = axes[0, 0].scatter([], [], s=15, alpha=0.9, edgecolors="none")
axes[0, 0].set_xlim(hrx_lo, hrx_hi)
axes[0, 0].set_ylim(hry_hi, hry_lo)  # invert mags

# Top-right: Teff vs Log L
teff_line, = axes[0, 1].plot([], [], marker=",", linestyle="-", color="k", alpha=0.2)
teff_sc    = axes[0, 1].scatter([], [], s=15, alpha=0.9, edgecolors="none")
axes[0, 1].set_xlim(Teff_hi, Teff_lo)  # Teff inverted x
axes[0, 1].set_ylim(LogL_lo, LogL_hi)

# Bottom-left: Age vs color index
ci_line,   = axes[1, 0].plot([], [], marker=",", linestyle="-", color="k", alpha=0.2)
ci_sc      = axes[1, 0].scatter([], [], s=15, alpha=1.0, edgecolors="none")
axes[1, 0].set_xlim(age_lo, age_hi)
axes[1, 0].set_ylim(ci_lo, ci_hi)

# Bottom-right: Age vs filter mags — one Line2D per filter, created once
mag_lines = []
for label, series, color in mag_series:
    ln, = axes[1, 1].plot([], [], marker="o", linestyle="-", markersize=3,
                           alpha=0.8, label=label, color=color)
    mag_lines.append((ln, series))
axes[1, 1].set_xlim(age_lo, age_hi)
axes[1, 1].set_ylim(mag_hi, mag_lo)  # invert mags
if mag_series:
    axes[1, 1].legend()

# Phase legend is the same every frame — build it once outside the loop
if getattr(checker, "has_phase", False):
    legend_elements = checker.create_phase_legend()
    if len(legend_elements) > 0:
        n_phases = len(legend_elements)
        fig.legend(
            handles=legend_elements,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.53),
            ncol=max(1, (n_phases + 1) // 1),
            title_fontsize=10,
            fontsize=10,
            frameon=True,
            fancybox=True,
            shadow=True,
        )

plt.subplots_adjust(top=0.92)

# ── Render frames ─────────────────────────────────────────────────────────────
writer = FFMpegWriter(fps=FPS, codec="libx264", extra_args=["-pix_fmt", "yuv420p"])

with writer.saving(fig, OUT, dpi=DPI):
    N = age.size
    iterator = tqdm(range(N)) if do_tqdm else range(N)

    for i in iterator:
        sl = slice(0, i + 1)

        # Top-left: CMD/HR
        hr_line.set_data(hr_x[sl], hr_y[sl])
        hr_sc.set_offsets(np.c_[hr_x[sl], hr_y[sl]])
        if has_phases:
            hr_sc.set_facecolors(phase_colors[sl])

        # Top-right: Teff vs Log L
        teff_line.set_data(Teff[sl], Log_L[sl])
        teff_sc.set_offsets(np.c_[Teff[sl], Log_L[sl]])
        if has_phases:
            teff_sc.set_facecolors(phase_colors[sl])

        # Bottom-left: Age vs color index
        ci_line.set_data(age[sl], color_index[sl])
        ci_sc.set_offsets(np.c_[age[sl], color_index[sl]])
        if has_phases:
            ci_sc.set_facecolors(phase_colors[sl])

        # Bottom-right: Age vs filter mags
        for ln, series in mag_lines:
            ln.set_data(age[sl], series[sl])

        writer.grab_frame()

print(f"Wrote {age.size} frames -> {OUT}")
