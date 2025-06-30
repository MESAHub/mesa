#!/usr/bin/env python3
import glob
import os
import time
import matplotlib.pyplot as plt
import mesa_reader as mr
import numpy as np
from matplotlib.animation import FuncAnimation


class HistoryChecker:
    def __init__(self, history_file="../LOGS/history.data", refresh_interval=1, skip_rows=5):
        """
        Initialize the History checker with auto-refresh capability.
        
        Args:
            history_file: Path to the MESA history.data file
            refresh_interval: Time in seconds between refresh attempts
        """
        self.history_file = history_file
        self.refresh_interval = refresh_interval
        self.last_modified = None
        self.skip_rows = skip_rows

        # Create the figure and axes
        self.fig, self.axes = plt.subplots(
            2, 2, figsize=(14, 18), gridspec_kw={"hspace": 0.2, "wspace": 0.2}
        )
    
        self.update_flag = 0
        # Initial setup
        self.filter_columns = []
        self.update_data()
        self.setup_plot()
    
    def setup_plot(self):
        """Set up the plot with formatting and labels."""
        # Top-left plot: HR Diagram (Color vs. Magnitude)
        self.axes[0, 0].set_xlabel(self.hr_xlabel)
        self.axes[0, 0].set_ylabel(self.hr_ylabel)
        self.axes[0, 0].invert_yaxis()
        
        # Top-right plot: Teff vs. Log_L
        self.axes[0, 1].set_xlabel("Teff (K)")
        self.axes[0, 1].set_ylabel("Log_L")
        self.axes[0, 1].invert_xaxis()
        self.axes[0, 1].yaxis.set_label_position("right")
        self.axes[0, 1].yaxis.tick_right()
        
        # Bottom-left plot: Age vs. Color Index
        self.axes[1, 0].set_xlabel("Age")
        self.axes[1, 0].set_ylabel(f"Color ({self.hr_xlabel})")
        
        # Bottom-right plot: Age vs. All Filter Magnitudes
        self.axes[1, 1].set_xlabel("Age")
        self.axes[1, 1].set_ylabel("Magnitude")
        self.axes[1, 1].invert_yaxis()
        self.axes[1, 1].yaxis.set_label_position("right")
        self.axes[1, 1].yaxis.tick_right()
        
        plt.tight_layout()
    
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
        if not os.path.exists(self.history_file):
            print(f"Warning: History file {self.history_file} not found")
            return

        try:
            self.md = mr.MesaData(self.history_file)
            self.slice = slice(self.skip_rows, None)

            s = self.slice
            self.Teff = self.md.Teff[s]
            self.Log_L = self.md.log_L[s]
            self.Log_g = self.md.log_g[s]
            self.Log_R = self.md.log_R[s]
            self.Star_Age = self.md.star_age[s]
            self.Mag_bol = self.md.Mag_bol[s]
            self.Flux_bol = np.log10(self.md.Flux_bol[s])

            self.read_header_columns()
            self.setup_hr_diagram_params()

        except Exception as e:
            print(f"Error reading history data: {e}")

    
    def read_header_columns(self):
        """Read column headers from history file."""
        header_line = None
        with open(self.history_file, "r") as fp:
            for line in fp:
                if "model_number" in line:
                    header_line = line.strip()
                    break
        
        if header_line is None:
            print("Warning: Could not find header line with 'model_number'")
            self.all_cols = []
            self.filter_columns = []
            return
        
        # Split the header line on whitespace
        self.all_cols = header_line.split()
        
        # Find the index of Flux_bol
        try:
            flux_index = self.all_cols.index("Flux_bol")
            self.filter_columns = self.all_cols[flux_index + 1:]
            #print("Detected filter columns:", self.filter_columns)
        except ValueError:
            print("Warning: Could not find 'Flux_bol' column in header")
            self.filter_columns = []
    

    def setup_hr_diagram_params(self):
        s = self.slice

        if "Gbp" in self.filter_columns and "Grp" in self.filter_columns and "G" in self.filter_columns:
            self.hr_color = self.md.Gbp[s] - self.md.Grp[s]
            self.hr_mag = self.md.G[s]
            self.hr_xlabel = "Gbp - Grp"
            self.hr_ylabel = "G"
            self.color_index = self.hr_color  # DO NOT slice again
        else:
            if len(self.filter_columns) >= 2:
                f1 = self.filter_columns[0]
                f2 = self.filter_columns[1]

                try:
                    col1 = getattr(self.md, f1)[s]
                    col2 = getattr(self.md, f2)[s]
                except AttributeError:
                    col1 = self.md.data(f1)[s]
                    col2 = self.md.data(f2)[s]

                self.hr_color = col1 - col2
                self.hr_mag = col1
                self.hr_xlabel = f"{f1} - {f2}"
                self.hr_ylabel = f1
                self.color_index = self.hr_color
            else:
                print("Warning: Not enough filter columns to construct color index")
                n = len(self.Teff)
                self.hr_color = np.zeros(n)
                self.hr_mag = np.zeros(n)
                self.hr_xlabel = "Color Index"
                self.hr_ylabel = "Magnitude"
                self.color_index = self.hr_color

    
    def update_plot(self, frame):
        """Update the plot with new data if the file has changed."""
        if not self.check_for_changes():
            return
        
        if self.update_flag < 5:
            print(f"Detected changes in {self.history_file}, updating plot...")
            self.update_flag = self.update_flag + 1
        elif self.update_flag == 5:
            print(f"... MESA is clearly running so no more updates about: {self.history_file}")
            self.update_flag = self.update_flag + 1
        
        # Update data
        self.update_data()
        
        # Clear all axes
        for ax in self.axes.flatten():
            ax.clear()
        
        # Reset plot formatting
        self.setup_plot()
        
        # Top-left plot: HR Diagram (Color vs. Magnitude)
        self.axes[0, 0].plot(self.hr_color, self.hr_mag, "go")
        
        # Top-right plot: Teff vs. Log_L
        self.axes[0, 1].plot(self.Teff, self.Log_L, "go")
        
        # Bottom-left plot: Age vs. Color Index
        self.axes[1, 0].plot(self.Star_Age, self.color_index, "kx")
        
        for filt in self.filter_columns:
            try:
                col_data = getattr(self.md, filt)[self.slice]
            except AttributeError:
                try:
                    col_data = self.md.data(filt)[self.slice]
                except:
                    print(f"Warning: Could not retrieve data for filter {filt}")
                    continue

            self.axes[1, 1].plot(self.Star_Age, col_data, marker="o", linestyle="-", label=filt)

        
        # Add legend to filter plot
        self.axes[1, 1].legend()
        
        # Update the figure
        self.fig.canvas.draw_idle()
    
    def run(self):
        """Start the auto-refreshing display."""
        print(f"Monitoring history file: {self.history_file}")
        print(f"Refresh interval: {self.refresh_interval} seconds")
        
        # Initial plot
        self.update_plot(0)
        
        # Set up animation
        self.animation = FuncAnimation(
            self.fig, 
            self.update_plot,
            interval=self.refresh_interval * 1000,  # Convert to milliseconds
            cache_frame_data=False
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
    
    checker = HistoryChecker(history_file=history_file, refresh_interval=1, skip_rows = 200)
    checker.run()


if __name__ == "__main__":
    main()