#!/usr/bin/env python3.8
####################################################
#
# Author: M Joyce, Modified by N. Miller,
#         Further modified to automatically use all filters
#         by reading the header of the history file.
#         Enhanced with MESA's built-in evolutionary phase identification.
#
####################################################
import glob
import matplotlib.pyplot as plt
import mesa_reader as mr
import numpy as np


def get_mesa_phase_info(phase_code):
    """
    Map MESA's phase_of_evolution integer codes to phase names and colors.
    Based on MESA's exact internal phase definitions from star_data_def.inc
    """
    # MESA phase codes (exact from source code)
    phase_map = {
        -1: ("Relax", "#C0C0C0"),  # Silver - Relaxation phase
        1: ("Starting", "#E6E6FA"),  # Lavender - Starting phase
        2: ("Pre-MS", "#FF69B4"),  # Hot pink - Pre-main sequence
        3: ("ZAMS", "#00FF00"),  # Bright green - Zero-age main sequence
        4: ("IAMS", "#0000FF"),  # Blue - Intermediate-age main sequence
        5: ("TAMS", "#FF8C00"),  # Dark orange - Terminal-age main sequence
        6: ("He-Burn", "#8A2BE2"),  # Blue violet - Helium burning (general)
        7: ("ZACHeB", "#9932CC"),  # Dark orchid - Zero-age core helium burning
        8: ("TACHeB", "#BA55D3"),  # Medium orchid - Terminal-age core helium burning
        9: ("TP-AGB", "#8B0000"),  # Dark red - Thermally pulsing AGB
        10: ("C-Burn", "#FF4500"),  # Orange red - Carbon burning
        11: ("Ne-Burn", "#FF6347"),  # Tomato - Neon burning
        12: ("O-Burn", "#FF8C00"),  # Dark orange - Oxygen burning
        13: ("Si-Burn", "#FFA500"),  # Orange - Silicon burning
        14: ("WDCS", "#708090"),  # Slate gray - White dwarf cooling sequence
    }

    return phase_map.get(phase_code, ("Unknown", "#808080"))


def get_phase_info_from_mesa(md):
    """Get evolutionary phase information using MESA's phase_of_evolution."""

    # Check if phase_of_evolution exists in the data
    if hasattr(md, "phase_of_evolution"):
        phase_codes = md.phase_of_evolution
    else:
        print("Warning: phase_of_evolution not found in history file.")
        print("Make sure to add 'phase_of_evolution' to your history_columns.list")
        # Fallback to unknown phase
        n_models = len(md.model_number)
        phase_codes = np.full(n_models, -1)

    phases = []
    phase_colors = []

    for code in phase_codes:
        phase_name, color = get_mesa_phase_info(int(code))
        phases.append(phase_name)
        phase_colors.append(color)

    return phases, phase_colors


def read_header_columns(history_file):
    """Read column headers from history file."""
    header_line = None
    with open(history_file, "r") as fp:
        for line in fp:
            if "model_number" in line:
                header_line = line.strip()
                break

    if header_line is None:
        print("Warning: Could not find header line with 'model_number'")
        return [], []

    # Split the header line on whitespace
    all_cols = header_line.split()

    # Find the index of Flux_bol
    try:
        flux_index = all_cols.index("Flux_bol")
        filter_columns = all_cols[flux_index + 1 :]
    except ValueError:
        print("Warning: Could not find 'Flux_bol' column in header")
        filter_columns = []

    return all_cols, filter_columns


def setup_hr_diagram_params(md, filter_columns):
    """Set up parameters for HR diagram based on available filters."""
    if "Gbp" in filter_columns and "Grp" in filter_columns and "G" in filter_columns:
        hr_color = md.Gbp - md.Grp
        hr_mag = md.G
        hr_xlabel = "Gbp - Grp"
        hr_ylabel = "G"
        color_index = hr_color
    else:
        if len(filter_columns) >= 2:
            # Use the first two filters
            f1 = filter_columns[0]
            f2 = filter_columns[1]

            # Retrieve the data using getattr or data method
            try:
                col1 = getattr(md, f1)
                col2 = getattr(md, f2)
            except AttributeError:
                col1 = md.data(f1)
                col2 = md.data(f2)

            hr_color = col1 - col2
            hr_mag = col1
            hr_xlabel = f"{f1} - {f2}"
            hr_ylabel = f1
            color_index = hr_color
        else:
            # Default values if not enough filters
            print("Warning: Not enough filter columns to construct color index")
            hr_color = np.zeros_like(md.Teff)
            hr_mag = np.zeros_like(md.Teff)
            hr_xlabel = "Color Index"
            hr_ylabel = "Magnitude"
            color_index = hr_color

    return hr_color, hr_mag, hr_xlabel, hr_ylabel, color_index


def create_phase_plots(history_file="../LOGS/history.data"):
    """
    Create plots with MESA's evolutionary phase color coding.
    """
    # Read the MESA data
    md = mr.MesaData(history_file)

    # Basic stellar parameters
    Teff = md.Teff
    Log_L = md.log_L
    # Log_g = md.log_g
    # Log_R = md.log_R
    Star_Age = md.star_age
    # Mag_bol = md.Mag_bol
    # Flux_bol = np.log10(md.Flux_bol)

    # Read headers and get filter info
    all_cols, filter_columns = read_header_columns(history_file)

    # Set up HR diagram parameters
    hr_color, hr_mag, hr_xlabel, hr_ylabel, color_index = setup_hr_diagram_params(
        md, filter_columns
    )

    # Get evolutionary phase information using MESA's built-in phases
    phases, phase_colors = get_phase_info_from_mesa(md)

    # Print phase statistics
    unique_phases, counts = np.unique(phases, return_counts=True)
    print("\nEvolutionary phases found:")
    for phase, count in zip(unique_phases, counts):
        print(f"  {phase}: {count} models")

    # Create the plots
    fig, axes = plt.subplots(
        2, 2, figsize=(14, 18), gridspec_kw={"hspace": 0.3, "wspace": 0.2}
    )

    # Top-left plot: HR Diagram (Color vs. Magnitude) with phase colors
    axes[0, 0].scatter(
        hr_color, hr_mag, c=phase_colors, s=20, alpha=0.7, edgecolors="none"
    )

    axes[0, 0].set_xlabel(hr_xlabel, fontsize=14)
    axes[0, 0].set_ylabel(hr_ylabel, fontsize=14)
    axes[0, 0].invert_yaxis()
    axes[0, 0].xaxis.set_ticks_position("top")
    axes[0, 0].xaxis.set_label_position("top")
    axes[0, 0].grid(True, alpha=0.3)

    # Top-right plot: Teff vs. Log_L with phase colors
    axes[0, 1].scatter(Teff, Log_L, c=phase_colors, s=20, alpha=0.7, edgecolors="none")

    axes[0, 1].set_xlabel("Teff (K)", fontsize=14)
    axes[0, 1].set_ylabel("Log L/L☉", fontsize=14)
    axes[0, 1].invert_xaxis()
    axes[0, 1].yaxis.set_label_position("right")
    axes[0, 1].yaxis.tick_right()
    axes[0, 1].xaxis.set_ticks_position("top")
    axes[0, 1].xaxis.set_label_position("top")
    axes[0, 1].grid(True, alpha=0.3)

    # Bottom-left plot: Age vs. Color Index with phase colors
    axes[1, 0].scatter(
        Star_Age, color_index, c=phase_colors, s=20, alpha=0.7, edgecolors="none"
    )

    axes[1, 0].set_xlabel("Age (years)", fontsize=14)
    axes[1, 0].set_ylabel(f"Color ({hr_xlabel})", fontsize=14)
    axes[1, 0].grid(True, alpha=0.3)

    # Bottom-right plot: Age vs. All Filter Magnitudes
    for filt in filter_columns:
        try:
            col_data = getattr(md, filt)
        except AttributeError:
            try:
                col_data = md.data(filt)
            except Exception:
                print(f"Warning: Could not retrieve data for filter {filt}")
                continue

        axes[1, 1].plot(
            Star_Age,
            col_data,
            marker="o",
            linestyle="-",
            label=filt,
            markersize=3,
            alpha=0.8,
        )

    axes[1, 1].set_xlabel("Age (years)", fontsize=14)
    axes[1, 1].set_ylabel("Magnitude", fontsize=14)
    axes[1, 1].invert_yaxis()
    axes[1, 1].yaxis.set_label_position("right")
    axes[1, 1].yaxis.tick_right()
    axes[1, 1].grid(True, alpha=0.3)

    # Create a custom legend for evolutionary phases (only show phases that appear)
    unique_phases = []
    unique_colors = []
    for phase, color in zip(phases, phase_colors):
        if phase not in unique_phases:
            unique_phases.append(phase)
            unique_colors.append(color)

    # Add legend to the first subplot
    legend_elements = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=color,
            markersize=8,
            label=phase,
            markeredgecolor="none",
        )
        for phase, color in zip(unique_phases, unique_colors)
    ]

    # Position legend outside the plot
    axes[0, 0].legend(
        handles=legend_elements,
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        title="Evolutionary Phase",
    )

    plt.tight_layout()
    return fig, axes, phases, phase_colors


def main():
    # Locate the history.data file
    try:
        history_file = glob.glob("../LOGS/history.data")[0]
    except IndexError:
        history_file = "../LOGS/history.data"
        print(f"Warning: No history.data file found, will check for {history_file}")

    print("Using MESA's built-in phase_of_evolution for phase identification")
    print("Make sure 'phase_of_evolution' is in your history_columns.list")
    print("\nMESA Phase Definitions:")
    print("  -1: Relax, 1: Starting, 2: Pre-MS, 3: ZAMS, 4: IAMS, 5: TAMS")
    print("  6: He-Burn, 7: ZACHeB, 8: TACHeB, 9: TP-AGB")
    print("  10: C-Burn, 11: Ne-Burn, 12: O-Burn, 13: Si-Burn, 14: WDCS")

    fig, axes, phases, phase_colors = create_phase_plots(history_file)
    plt.show()


if __name__ == "__main__":
    main()
