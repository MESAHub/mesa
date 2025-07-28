#!/usr/bin/env python3
import glob
import os
import matplotlib.pyplot as plt
import mesa_reader as mr
import numpy as np
from matplotlib.animation import FuncAnimation

# Import functions from static version for consistency
from static_HISTORY_check import (
    # get_mesa_phase_info,
    read_header_columns,
    setup_hr_diagram_params,
)


def get_improved_mesa_phase_info(phase_code):
    """
    Map MESA's phase_of_evolution integer codes to phase names and more intuitive colors.
    Colors follow stellar evolution logic: blue (hot/young) -> red (cool/evolved) -> white (remnants)
    """
    # Improved phase colors that make physical sense
    phase_map = {
        -1: ("Relax", "#808080"),  # Gray - Relaxation phase
        1: ("Starting", "#E6E6FA"),  # Lavender - Starting phase
        2: ("Pre-MS", "#9370DB"),  # Medium purple - Pre-main sequence (young)
        3: ("ZAMS", "#4169E1"),  # Royal blue - Zero-age MS (hot, young MS)
        4: ("IAMS", "#00CED1"),  # Dark turquoise - Intermediate-age MS
        5: ("TAMS", "#FF6347"),  # Tomato - Terminal-age MS (cooler, evolved MS)
        6: ("He-Burn", "#FF4500"),  # Orange red - Helium burning (post-MS)
        7: ("ZACHeB", "#FF8C00"),  # Dark orange - Zero-age core helium burning
        8: ("TACHeB", "#FFA500"),  # Orange - Terminal-age core helium burning
        9: ("TP-AGB", "#DC143C"),  # Crimson - Thermally pulsing AGB (very evolved)
        10: ("C-Burn", "#B22222"),  # Fire brick - Carbon burning (massive stars)
        11: ("Ne-Burn", "#8B0000"),  # Dark red - Neon burning
        12: ("O-Burn", "#800000"),  # Maroon - Oxygen burning
        13: ("Si-Burn", "#654321"),  # Dark brown - Silicon burning (pre-collapse)
        14: ("WDCS", "#F5F5F5"),  # White smoke - White dwarf cooling sequence
    }

    return phase_map.get(phase_code, ("Unknown", "#696969"))  # Dim gray for unknown


def get_phase_info_from_mesa(md):
    """Get evolutionary phase information using MESA's phase_of_evolution with improved colors."""

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
        phase_name, color = get_improved_mesa_phase_info(int(code))
        phases.append(phase_name)
        phase_colors.append(color)

    return phases, phase_colors


def get_filter_colors(filter_names):
    """Generate distinct colors for different filters."""
    # Use a colormap that provides good contrast
    cmap = plt.cm.Set1  # Good for distinct colors
    colors = []
    for i, filt in enumerate(filter_names):
        color = cmap(i / max(len(filter_names), 1))
        colors.append(color)
    return colors


class HistoryChecker:
    def __init__(self, history_file="../LOGS/history.data", refresh_interval=0.01):
        """
        Initialize the History checker with auto-refresh capability and MESA phase color coding.

        Args:
            history_file: Path to the MESA history.data file
            refresh_interval: Time in seconds between refresh attempts
        """
        self.history_file = history_file
        self.refresh_interval = refresh_interval
        self.last_modified = None

        # Create the figure and axes
        self.fig, self.axes = plt.subplots(
            2, 2, figsize=(14, 18), gridspec_kw={"hspace": 0.01, "wspace": 0.01}
        )

        self.update_flag = 0
        # Initial setup
        self.filter_columns = []
        self.phases = []
        self.phase_colors = []
        self.filter_colors = []
        self.update_data()
        self.setup_plot()

    def setup_plot(self):
        """Set up the plot with formatting and labels (titles removed, labels enlarged, top x-axis added)."""
        # Top-left plot: HR Diagram (Color vs. Magnitude)
        self.axes[0, 0].set_xlabel(self.hr_xlabel, fontsize=16)
        self.axes[0, 0].set_ylabel(self.hr_ylabel, fontsize=16)
        self.axes[0, 0].invert_yaxis()
        self.axes[0, 0].xaxis.set_ticks_position("top")
        self.axes[0, 0].xaxis.set_label_position("top")
        self.axes[0, 0].grid(True, alpha=0.3)

        # Top-right plot: Teff vs. Log_L
        self.axes[0, 1].set_xlabel("Teff (K)", fontsize=16)
        self.axes[0, 1].set_ylabel("Log L/L☉", fontsize=16)
        self.axes[0, 1].invert_xaxis()
        self.axes[0, 1].yaxis.set_label_position("right")
        self.axes[0, 1].yaxis.tick_right()
        self.axes[0, 1].xaxis.set_ticks_position("top")
        self.axes[0, 1].xaxis.set_label_position("top")
        self.axes[0, 1].grid(True, alpha=0.3)

        # Bottom-left plot: Age vs. Color Index
        self.axes[1, 0].set_xlabel("Age (years)", fontsize=16)
        self.axes[1, 0].set_ylabel(f"Color ({self.hr_xlabel})", fontsize=16)
        self.axes[1, 0].grid(True, alpha=0.3)

        # Bottom-right plot: Age vs. All Filter Magnitudes
        self.axes[1, 1].set_xlabel("Age (years)", fontsize=16)
        self.axes[1, 1].set_ylabel("Magnitude", fontsize=16)
        self.axes[1, 1].invert_yaxis()
        self.axes[1, 1].yaxis.set_label_position("right")
        self.axes[1, 1].yaxis.tick_right()
        self.axes[1, 1].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.subplots_adjust(top=0.92)  # Adjust if needed for legend spacing

    def check_for_changes(self):
        """Check if history file has been modified since last check."""
        if not os.path.exists(self.history_file):
            print(f"Warning: History file {self.history_file} not found")
            return False

        current_mtime = os.path.getmtime(self.history_file)

        if self.last_modified is None or current_mtime > self.last_modified:
            self.last_modified = current_mtime
            return True

        return False

    def update_data(self):
        """Read data from history file and extract relevant columns."""
        if not os.path.exists(self.history_file):
            print(f"Warning: History file {self.history_file} not found")
            return

        try:
            # Read the MESA data
            self.md = mr.MesaData(self.history_file)

            # Basic stellar parameters
            self.Teff = self.md.Teff
            self.Log_L = self.md.log_L
            self.Log_g = self.md.log_g
            self.Log_R = self.md.log_R
            self.Star_Age = self.md.star_age
            self.Mag_bol = self.md.Mag_bol
            self.Flux_bol = np.log10(self.md.Flux_bol)

            # Read header columns using imported function
            self.all_cols, self.filter_columns = read_header_columns(self.history_file)

            # Set up HR diagram parameters using imported function
            (
                self.hr_color,
                self.hr_mag,
                self.hr_xlabel,
                self.hr_ylabel,
                self.color_index,
            ) = setup_hr_diagram_params(self.md, self.filter_columns)

            # Get evolutionary phase information using MESA's built-in phases with improved colors
            self.phases, self.phase_colors = get_phase_info_from_mesa(self.md)

            # Get filter colors for the magnitude plot
            self.filter_colors = get_filter_colors(self.filter_columns)

        except Exception as e:
            print(f"Error reading history data: {e}")

    def create_phase_legend(self):
        """Create legend for evolutionary phases."""
        unique_phases = []
        unique_colors = []
        for phase, color in zip(self.phases, self.phase_colors):
            if phase not in unique_phases:
                unique_phases.append(phase)
                unique_colors.append(color)

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

        return legend_elements

    def update_plot(self, frame):
        """Update the plot with new data if the file has changed."""
        if not self.check_for_changes():
            return

        if self.update_flag < 5:
            print(f"Detected changes in {self.history_file}, updating plot...")
            self.update_flag = self.update_flag + 1
        elif self.update_flag == 5:
            print(
                f"... MESA is clearly running so no more updates about: {self.history_file}"
            )
            self.update_flag = self.update_flag + 1

        # Update data
        self.update_data()

        # Clear all axes
        for ax in self.axes.flatten():
            ax.clear()

        # Reset plot formatting
        self.setup_plot()

        # Top-left plot: HR Diagram with phase colors
        if len(self.phases) > 0:
            self.axes[0, 0].plot(
                self.hr_color,
                self.hr_mag,
                marker=",",
                linestyle="-",
                color="k",
                alpha=0.2,
            )

            self.axes[0, 0].scatter(
                self.hr_color,
                self.hr_mag,
                c=self.phase_colors,
                s=15,
                alpha=0.9,
                linestyle="-",
                edgecolors="none",
            )
        else:
            self.axes[0, 0].plot(self.hr_color, self.hr_mag, "go")

        # Top-right plot: Teff vs. Log_L with phase colors
        if len(self.phases) > 0:
            self.axes[0, 1].plot(
                self.Teff,
                self.Log_L,
                marker=",",
                linestyle="-",
                color="k",
                alpha=0.2,
            )

            self.axes[0, 1].scatter(
                self.Teff,
                self.Log_L,
                c=self.phase_colors,
                s=15,
                alpha=0.9,
                linestyle="-",
                edgecolors="none",
            )
        else:
            self.axes[0, 1].plot(self.Teff, self.Log_L, "go")

        # Bottom-left plot: Age vs. Color Index with phase colors
        if len(self.phases) > 0:
            self.axes[1, 0].plot(
                self.Star_Age,
                self.color_index,
                marker=",",
                linestyle="-",
                color="k",
                alpha=0.2,
            )

            self.axes[1, 0].scatter(
                self.Star_Age,
                self.color_index,
                c=self.phase_colors,
                s=15,
                alpha=1,
                linestyle="-",
                edgecolors="none",
            )

        else:
            self.axes[1, 0].plot(self.Star_Age, self.color_index, "kx")

        # Bottom-right plot: Age vs. All Filter Magnitudes with filter-specific colors
        for i, filt in enumerate(self.filter_columns):
            # Retrieve filter magnitude data
            try:
                col_data = getattr(self.md, filt)
            except AttributeError:
                try:
                    col_data = self.md.data(filt)
                except Exception:
                    print(f"Warning: Could not retrieve data for filter {filt}")
                    continue

            # Use filter-specific color
            color = self.filter_colors[i] if i < len(self.filter_colors) else "black"

            self.axes[1, 1].plot(
                self.Star_Age,
                col_data,
                marker="o",
                linestyle="-",
                label=filt,
                color=color,
                markersize=3,
                alpha=0.8,
            )

        # Add legend to filter plot
        self.axes[1, 1].legend()

        # Add evolutionary phase legend at the top of the figure
        if len(self.phases) > 0:
            legend_elements = self.create_phase_legend()
            # Calculate number of columns for max 2 rows
            n_phases = len(legend_elements)
            ncol = max(1, (n_phases + 1) // 1)  # Ceiling division to get max 2 rows

            self.fig.legend(
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

        # Adjust layout to make room for top legend
        plt.tight_layout()
        plt.subplots_adjust(top=0.92)

        # Update the figure
        self.fig.canvas.draw_idle()

    def run(self):
        """Start the auto-refreshing display."""
        print(f"Monitoring history file: {self.history_file}")
        print(f"Refresh interval: {self.refresh_interval} seconds")
        print("Using MESA's built-in phase_of_evolution for color coding")
        print("Make sure 'phase_of_evolution' is in your history_columns.list")
        print("\nMESA Phase Definitions:")
        print("  -1: Relax, 1: Starting, 2: Pre-MS, 3: ZAMS, 4: IAMS, 5: TAMS")
        print("  6: He-Burn, 7: ZACHeB, 8: TACHeB, 9: TP-AGB")
        print("  10: C-Burn, 11: Ne-Burn, 12: O-Burn, 13: Si-Burn, 14: WDCS")

        # Initial plot
        self.update_plot(0)

        # Set up animation
        self.animation = FuncAnimation(
            self.fig,
            self.update_plot,
            interval=self.refresh_interval * 1000,  # Convert to milliseconds
            cache_frame_data=False,
        )

        # Show the plot (this will block until window is closed)
        plt.show()


def main():
    # Locate the history.data file
    try:
        history_file = glob.glob("../LOGS/history.data")[0]
    except IndexError:
        history_file = "../LOGS/history.data"  # Default path if not found
        print(f"Warning: No history.data file found, will check for {history_file}")

    checker = HistoryChecker(history_file=history_file, refresh_interval=0.1)
    checker.run()


if __name__ == "__main__":
    main()
