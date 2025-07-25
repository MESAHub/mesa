import os
import csv
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation

# Use matplotlib's mathtext renderer instead of LaTeX to avoid dependency issues
plt.rcParams['text.usetex'] = False

class SEDChecker:
    def __init__(
        self,
        directory="../SED/",
        xlim=None,
        ylim=None,
        refresh_interval=0.1,
        save_video=False,
        video_filename="sed_animation.mp4",
        video_fps=10,
        video_duration=30,
    ):
        self.directory = directory
        self.xlim = xlim
        self.ylim = ylim
        self.refresh_interval = refresh_interval
        self.file_timestamps = {}
        self.fig, self.ax = plt.subplots(figsize=(12, 6))
        self.full_sed_plotted = False

        self.save_video = save_video
        self.video_filename = video_filename
        self.video_fps = video_fps
        self.video_duration = video_duration
        self.video_frame_count = int(video_fps * video_duration)
        self.frame_counter = 0

        self.setup_plot()

    def setup_plot(self):
        """Set up the plot with formatting and labels."""
        # Use matplotlib's mathtext for math rendering
        self.ax.set_xlabel(r"Wavelength ($\AA$)", fontsize=14)
        self.ax.set_ylabel(r"Flux (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=14)
        # Set logarithmic scale for x-axis to better visualize wide wavelength range
        self.ax.set_xscale('log')
        # Only apply ticklabel formatting to y-axis (log scale incompatible with ticklabel_format)
        self.ax.yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
        if self.xlim:
            self.ax.set_xlim(self.xlim)
        if self.ylim:
            self.ax.set_ylim(self.ylim)
        plt.tight_layout()

    def find_output_files(self):
        return [
            f
            for f in os.listdir(self.directory)
            if "SED.csv" in f and os.path.isfile(os.path.join(self.directory, f))
        ]

    def normalize(self, data):
        if data is None or len(data) == 0:
            return np.array([])
        try:
            if isinstance(data[0], str):
                data = np.array([float(x) for x in data])
            data = data.astype(float)
        except (ValueError, TypeError):
            return np.zeros_like(data, dtype=float)
        try:
            data_min = np.min(data)
            data_max = np.max(data)
            if data_max != data_min:
                return (data - data_min) / (data_max - data_min)
            else:
                return np.zeros_like(data)
        except Exception:
            return np.zeros_like(data)

    def check_for_changes(self):
        try:
            output_files = self.find_output_files()
            changed = False
            if set(output_files) != set(self.file_timestamps.keys()):
                changed = True
            for file_path in output_files:
                full_path = os.path.join(self.directory, file_path)
                current_mtime = os.path.getmtime(full_path)
                if (
                    file_path not in self.file_timestamps
                    or self.file_timestamps[file_path] != current_mtime
                ):
                    self.file_timestamps[file_path] = current_mtime
                    changed = True
            return changed
        except Exception:
            return False

    def calculate_data_range(self, data_list):
        all_min = float("inf")
        all_max = float("-inf")
        for data in data_list:
            if data is None or len(data) == 0:
                continue
            try:
                if isinstance(data[0], str):
                    data = np.array([float(x) for x in data])
                data = np.array(data, dtype=float)
                all_min = min(all_min, np.min(data))
                all_max = max(all_max, np.max(data))
            except Exception:
                continue
        return (None, None) if all_min == float("inf") else (all_min, all_max)

    def ensure_numeric(self, data):
        if data is None or len(data) == 0:
            return np.array([])
        try:
            if isinstance(data[0], str):
                data = np.array([float(x) for x in data])
            return np.array(data, dtype=float)
        except (ValueError, TypeError):
            return np.zeros(len(data))

    def read_csv_file(self, file_path):
        """Read CSV file without pandas"""
        try:
            with open(file_path, 'r') as f:
                reader = csv.DictReader(f)
                
                # Strip whitespace from column names
                fieldnames = [field.strip() for field in reader.fieldnames]
                
                data = {}
                for field in fieldnames:
                    data[field] = []
                
                for row in reader:
                    # Skip rows with any missing values
                    if any(not value.strip() for value in row.values()):
                        continue
                        
                    for field in fieldnames:
                        try:
                            data[field].append(float(row[field.strip()]))
                        except (ValueError, KeyError):
                            data[field].append(np.nan)
                
                # Remove any columns that are all NaN
                for field in list(data.keys()):
                    if all(np.isnan(x) for x in data[field]):
                        del data[field]
                
                return data
        except Exception:
            return {}

    def update_plot(self, frame):
        if self.save_video:
            self.frame_counter += 1
            if self.frame_counter >= self.video_frame_count:
                plt.close(self.fig)
                return

        if not self.check_for_changes() and not self.save_video:
            return

        current_xlim = self.ax.get_xlim() if not self.xlim else self.xlim
        current_ylim = self.ax.get_ylim() if not self.ylim else self.ylim
        self.ax.clear()
        self.full_sed_plotted = False

        output_files = self.find_output_files()
        if not output_files:
            self.setup_plot()
            self.ax.set_xlim(current_xlim)
            self.ax.set_ylim(current_ylim)
            return

        all_flux_data = []
        all_convolved_flux_data = []
        sed_wavelengths = np.array([])  # Track SED wavelengths for smart x-axis clipping
        sed_flux = np.array([])  # Track SED flux for smart x-axis clipping

        # Legend and color mapping for filters
        legend_labels = {
            "Gbp": r"Gaia G$_{\mathrm{BP}}$",
            "Grp": r"Gaia G$_{\mathrm{RP}}$", 
            "G": "Gaia G",
            "Grvs": r"Gaia G$_{\mathrm{RVS}}$",
            "VEGA": "VEGA",
        }
        
        # Color scheme: main filters clear, others from Monet's water lilies palette
        filter_colors = {
            "Gbp": "#0000FF",       # Bright blue for Blue photometer
            "Grp": "#FF0000",       # Bright red for Red photometer  
            "G": "#228B22",         # Forest green (less bright) for broad band
            "Grvs": "#9370DB",      # Monet-inspired soft purple
            "VEGA": "#DDA0DD",      # Monet-inspired soft plum/lavender
        }
        
        displayed_labels = set()
        legend_handles = []
        legend_labels_list = []

        for file_path in output_files:
            linestyle = "--" if "VEGA" in file_path else "-"
            try:
                # Read CSV file without pandas
                data_dict = self.read_csv_file(os.path.join(self.directory, file_path))
                
                if not data_dict:
                    continue
                
                wavelengths = self.ensure_numeric(data_dict.get("wavelengths", []))
                flux = self.ensure_numeric(data_dict.get("fluxes", []))
                convolved_flux = self.ensure_numeric(data_dict.get("convolved_flux", []))

                if len(flux) > 0:
                    all_flux_data.append(flux)
                    # Store SED data for intelligent x-axis clipping - only for the actual SED
                    if len(wavelengths) > 0 and len(flux) > 0 and not self.full_sed_plotted:
                        sed_wavelengths = wavelengths
                        sed_flux = flux
                if len(convolved_flux) > 0:
                    all_convolved_flux_data.append(convolved_flux)

                # Plot full SED only once
                if len(wavelengths) > 0 and len(flux) > 0 and not self.full_sed_plotted:
                    line = self.ax.plot(
                        wavelengths,
                        flux,
                        color="black",
                        linewidth=1.5,
                        linestyle=linestyle,
                    )[0]
                    legend_handles.append(line)
                    legend_labels_list.append("Full SED")
                    self.full_sed_plotted = True

                # Plot convolved flux with proper legend handling and physical colors
                if len(wavelengths) > 0 and len(convolved_flux) > 0:
                    short_label = next(
                        (v for k, v in legend_labels.items() if k in file_path), None
                    )
                    # Get the corresponding color for this filter
                    filter_key = next(
                        (k for k in legend_labels.keys() if k in file_path), None
                    )
                    color = filter_colors.get(filter_key, None) if filter_key else None
                    
                    if short_label and short_label not in displayed_labels:
                        line = self.ax.plot(
                            wavelengths,
                            convolved_flux,
                            linewidth=1,
                            linestyle=linestyle,
                            color=color,
                        )[0]
                        legend_handles.append(line)
                        legend_labels_list.append(short_label)
                        displayed_labels.add(short_label)
                    else:
                        # Plot without adding to legend
                        self.ax.plot(
                            wavelengths,
                            convolved_flux,
                            linewidth=1,
                            linestyle=linestyle,
                            color=color,
                        )
            except Exception:
                continue

        # Create legend with explicit handles and labels
        if legend_handles and legend_labels_list:
            self.ax.legend(legend_handles, legend_labels_list, loc="best", fontsize=10)
            
        self.setup_plot()

        # Smart x-axis scaling: use wavelengths where flux > 0.1% of max flux
        if not self.xlim and len(sed_wavelengths) > 0 and len(sed_flux) > 0:
            max_flux = np.max(sed_flux)
            threshold = 0.001 * max_flux
            
            # Create subset of wavelengths where flux > threshold
            flux_mask = np.array(sed_flux) > threshold
            if np.any(flux_mask):
                subset_wavelengths = np.array(sed_wavelengths)[flux_mask]
                min_wave = np.min(subset_wavelengths)  
                max_wave = np.max(subset_wavelengths)
            
                # Just add 10% padding on each side in log space
                log_min = np.log10(max(min_wave, 50))  # Don't go below 50 Å
                log_max = np.log10(min(max_wave, 100000))  # Don't go above 100,000 Å
                log_range = log_max - log_min
                padding = log_range * 0.00001
                
                self.ax.set_xlim([10**(log_min - padding), 10**(log_max + padding)])

        if not self.ylim and all_convolved_flux_data:
            min_val, max_val = self.calculate_data_range(
                all_convolved_flux_data + all_flux_data
            )
            if min_val is not None and max_val is not None:
                padding = (max_val - min_val) * 0.1
                self.ax.set_ylim([min_val, max_val + padding])
            else:
                self.ax.set_ylim(current_ylim)
        elif not self.ylim:
            self.ax.set_ylim(current_ylim)

        # Keep current x-axis limits if manually set, otherwise auto-scaling is handled above
        if self.xlim:
            self.ax.set_xlim(self.xlim)

        self.fig.canvas.draw_idle()

    def run(self):
        print(f"Monitoring directory: {self.directory}")
        print(f"Refresh interval: {self.refresh_interval} seconds")
        if self.save_video:
            print(f"Saving video to {self.video_filename}")
            try:
                writer = FFMpegWriter(
                    fps=self.video_fps,
                    metadata=dict(title="SED Animation", artist="MESA Colors"),
                    bitrate=5000,
                )
                self.animation = FuncAnimation(
                    self.fig,
                    self.update_plot,
                    frames=self.video_frame_count,
                    interval=self.refresh_interval * 1000,
                    cache_frame_data=False,
                )
                self.animation.save(self.video_filename, writer=writer)
                print(f"Video saved to {self.video_filename}")
            except Exception as e:
                print(f"Error saving video: {e}")
        else:
            try:
                self.animation = FuncAnimation(
                    self.fig,
                    self.update_plot,
                    interval=self.refresh_interval * 1000,
                    cache_frame_data=False,
                )
                plt.show()
            except KeyboardInterrupt:
                print("\nMonitoring stopped by user.")
            except Exception as e:
                print(f"Error in animation: {e}")


def main():
    checker = SEDChecker(
        directory="../SED/",  # Change if needed
        xlim=None,  # Auto-scale x-axis to SED data range
        ylim=None,
        refresh_interval=0.1,
        save_video=False,
        video_filename="sed_animation.mp4",
        video_fps=10,
        video_duration=30,
    )
    checker.run()


if __name__ == "__main__":
    main()