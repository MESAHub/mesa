import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FFMpegWriter, FuncAnimation

class SEDChecker:
    def __init__(
        self,
        directory="../LOGS/SED/",
        xlim=None,
        ylim=None,
        refresh_interval=2,
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
        self.ax.set_xlabel("Wavelength (Å)", fontsize=14)
        self.ax.set_ylabel("Flux (erg s⁻¹ cm⁻² Å⁻¹)", fontsize=14)
        # self.ax.set_title("Combined SEDs and Filters")  # Removed title
        self.ax.ticklabel_format(style="plain", useOffset=False)
        if self.xlim:
            self.ax.set_xlim(self.xlim)
        if self.ylim:
            self.ax.set_ylim(self.ylim)
        plt.tight_layout()

    def find_output_files(self):
        return [
            f for f in os.listdir(self.directory)
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
                    file_path not in self.file_timestamps or
                    self.file_timestamps[file_path] != current_mtime
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

        # Legend simplification
        legend_labels = {
            "Gbp": "Gaia G_BP",
            "Grp": "Gaia G_RP",
            "G": "Gaia G",
            "Grvs": "Gaia G_RVS",
            "VEGA": "VEGA"
        }
        displayed_labels = set()

        for file_path in output_files:
            linestyle = "--" if "VEGA" in file_path else "-"
            try:
                df = pd.read_csv(os.path.join(self.directory, file_path), delimiter=",", header=0)
                df = df.rename(columns=str.strip).dropna()
                wavelengths = self.ensure_numeric(df.get("wavelengths"))
                flux = self.ensure_numeric(df.get("fluxes"))
                convolved_flux = self.ensure_numeric(df.get("convolved_flux"))

                if len(flux) > 0:
                    all_flux_data.append(flux)
                if len(convolved_flux) > 0:
                    all_convolved_flux_data.append(convolved_flux)

                if len(wavelengths) > 0 and len(flux) > 0 and not self.full_sed_plotted:
                    self.ax.plot(
                        wavelengths, flux,
                        label="Full SED", color="black",
                        linewidth=1.5, linestyle=linestyle
                    )
                    self.full_sed_plotted = True

                if len(wavelengths) > 0 and len(convolved_flux) > 0:
                    short_label = next((v for k, v in legend_labels.items() if k in file_path), None)
                    label = short_label if short_label and short_label not in displayed_labels else None
                    if label:
                        displayed_labels.add(label)
                    self.ax.plot(
                        wavelengths, convolved_flux,
                        label=label, linewidth=1, linestyle=linestyle
                    )
            except Exception:
                continue

        self.ax.legend(loc="best", fontsize=10)
        self.setup_plot()

        if not self.ylim and all_convolved_flux_data:
            min_val, max_val = self.calculate_data_range(all_convolved_flux_data + all_flux_data)
            if min_val is not None and max_val is not None:
                padding = (max_val - min_val) * 0.1
                self.ax.set_ylim([min_val, max_val + padding])
            else:
                self.ax.set_ylim(current_ylim)
        elif not self.ylim:
            self.ax.set_ylim(current_ylim)

        if not self.xlim:
            self.ax.set_xlim(current_xlim)

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
                    self.fig, self.update_plot,
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
                    self.fig, self.update_plot,
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
        directory="../LOGS/SED/",
        xlim=[2500, 20000],
        ylim=None,
        refresh_interval=5,
        save_video=False,
        video_filename="sed_animation.mp4",
        video_fps=10,
        video_duration=30,
    )
    checker.run()


if __name__ == "__main__":
    main()
