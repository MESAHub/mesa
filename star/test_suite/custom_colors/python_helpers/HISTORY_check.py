#!/usr/bin/env python3.8
####################################################
#
# Author: M Joyce, Modified by N. Miller,
#         Further modified to automatically use all filters
#         by reading the header of the history file.
#
####################################################
import numpy as np
import glob
import matplotlib.pyplot as plt
import mesa_reader as mr

# Locate the history.data file
f = glob.glob("../LOGS/history.data")[0]

# Read the MESA data (this part is unchanged and works)
md = mr.MesaData(f)
Teff = md.Teff
Log_L = md.log_L
Log_g = md.log_g
Log_R = md.log_R
Star_Age = md.star_age
Mag_bol = md.Mag_bol
Flux_bol = np.log10(md.Flux_bol)  # Convert to log scale

# -------------------------------------------------------------------
# Determine all history file column names by reading the header.
#
# The history file (by your description) has two header blocks.
# The second header block (which lists the history data column names)
# has a line that includes "model_number".
# We read that line and split it into column names.
# -------------------------------------------------------------------
header_line = None
with open(f, "r") as fp:
    for line in fp:
        if "model_number" in line:
            header_line = line.strip()
            break

if header_line is None:
    raise ValueError(
        "Could not find a header line containing 'model_number' in the file."
    )

# Split the header line on whitespace.
all_cols = header_line.split()
# For example, with your sample header, all_cols might be:
# ['model_number', 'num_zones', 'star_age', 'log_dt', 'star_mass', 'Teff',
#  'log_Teff', 'log_L', 'log_R', 'log_g', 'num_retries', 'num_iters',
#  'Mag_bol', 'Flux_bol', 'Gbp_bright', 'Gbp', 'Gbp_faint', 'G', 'Grp', 'Grvs']
# print("All history columns:", all_cols)

# -------------------------------------------------------------------
# Assume that all columns after "Flux_bol" are filters.
# (That is, Flux_bol is the last "fixed" column.)
# -------------------------------------------------------------------
try:
    flux_index = all_cols.index("Flux_bol")
except ValueError:
    raise ValueError("The header does not contain a 'Flux_bol' column.")

filter_names = all_cols[flux_index + 1 :]
print("Detected filter columns:", filter_names)

# -------------------------------------------------------------------
# For the HR diagram, we want a color index.
#
# If Gaia filters ("Gbp", "Grp", and "G") are among the filters,
# we use them (as in your working version).
# Otherwise, if at least two filter columns exist, we use the first two.
# -------------------------------------------------------------------
if "Gbp" in filter_names and "Grp" in filter_names and "G" in filter_names:
    hr_color = md.Gbp - md.Grp
    hr_mag = md.G
    hr_xlabel = "Gbp - Grp"
    hr_ylabel = "G"
    color_index = hr_color
else:
    if len(filter_names) >= 2:
        # Use the first two filters.
        f1 = filter_names[0]
        f2 = filter_names[1]
        # We can retrieve the data using getattr (mesa_reader creates attributes for each column)
        try:
            col1 = getattr(md, f1)
            col2 = getattr(md, f2)
        except AttributeError:
            # Alternatively, use md.data(key) if the attribute isn't present.
            col1 = md.data(f1)
            col2 = md.data(f2)
        hr_color = col1 - col2
        hr_mag = col1
        hr_xlabel = f"{f1} - {f2}"
        hr_ylabel = f1
        color_index = hr_color
    else:
        raise ValueError("Not enough filter columns found to construct a color index.")

# -------------------------------------------------------------------
# Create the plots.
#
# The top-left panel shows the HR diagram (color vs. magnitude),
# the top-right panel is Teff vs. Log_L,
# the bottom-left panel is Age vs. color index,
# and the bottom-right panel shows Age vs. magnitude for every filter.
# -------------------------------------------------------------------
fig, axes = plt.subplots(
    2, 2, figsize=(14, 18), gridspec_kw={"hspace": 0.2, "wspace": 0.2}
)

# Top-left plot: HR Diagram (Color vs. Magnitude)
axes[0, 0].plot(hr_color, hr_mag, "go")
axes[0, 0].set_xlabel(hr_xlabel)
axes[0, 0].set_ylabel(hr_ylabel)
axes[0, 0].invert_yaxis()

# Top-right plot: Teff vs. Log_L
axes[0, 1].plot(Teff, Log_L, "go")
axes[0, 1].set_xlabel("Teff (K)")
axes[0, 1].set_ylabel("Log_L")
axes[0, 1].invert_xaxis()
axes[0, 1].yaxis.set_label_position("right")
axes[0, 1].yaxis.tick_right()

# Bottom-left plot: Age vs. Color Index
axes[1, 0].plot(Star_Age, color_index, "kx")
axes[1, 0].set_xlabel("Age")
axes[1, 0].set_ylabel(f"Color ({hr_xlabel})")

# Bottom-right plot: Age vs. All Filter Magnitudes
for filt in filter_names:
    # Retrieve each filter's data.
    # Since mesa_reader creates properties for each column,
    # we can try using getattr; if that fails, we use md.data(filt)
    try:
        col_data = getattr(md, filt)
    except AttributeError:
        col_data = md.data(filt)
    axes[1, 1].plot(Star_Age, col_data, marker="o", linestyle="-", label=filt)

axes[1, 1].set_xlabel("Age")
axes[1, 1].set_ylabel("Magnitude")
axes[1, 1].invert_yaxis()
axes[1, 1].yaxis.set_label_position("right")
axes[1, 1].yaxis.tick_right()
axes[1, 1].legend()

plt.show()
