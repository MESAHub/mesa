import os
import csv
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation

# disable latex because it's a pain in the ass
plt.rcParams["text.usetex"] = False


class SEDChecker:
    """
    Monitors and plots spectral energy distribution (SED) data from MESA colors output

    This thing watches CSV files and makes pretty animated plots of stellar SEDs.
    Works with Gaia filter convolutions and all that jazz.
    """

    # ============================================================================
    # INITIALIZATION AND SETUP
    # ============================================================================

    def __init__(
        self,
        directory=None,
        xlim=None,
        ylim=None,
        refresh_interval=0.1,
        save_video=False,
        video_filename="sed_animation.mp4",
        video_fps=10,
        video_duration=30,
        inlist_file="../inlist_colors",
    ):
        """
        Initialize the SED checker thingy

        Most params are optional - it'll try to figure stuff out from the inlist file
        """

        # grab stellar parameters from the inlist file first
        self.stellar_params = self.parse_inlist_file(inlist_file)

        # basic settings - use inlist directory if none specified
        self.directory = (
            directory if directory is not None else self.stellar_params["sed_directory"]
        )
        self.xlim = xlim
        self.ylim = ylim
        self.refresh_interval = refresh_interval

        # keep track of file changes for monitoring
        self.file_timestamps = {}
        self.full_sed_plotted = False

        # set up the matplotlib figure
        self.fig, self.ax = plt.subplots(figsize=(12, 6))

        # video recording stuff
        self.save_video = save_video
        self.video_filename = video_filename
        self.video_fps = video_fps
        self.video_duration = video_duration
        self.video_frame_count = int(video_fps * video_duration)
        self.frame_counter = 0

        # do initial plot setup
        self.setup_plot()

    def parse_inlist_file(self, inlist_path="../inlist_colors"):
        """
        Parse the MESA inlist file to get stellar parameters

        This is kinda fragile and probably breaks sometimes but works most of the time
        """
        params = {
            "sed_directory": "../SED/",
            "distance": None,
            "initial_mass": None,
            "initial_z": None,
        }

        try:
            with open(inlist_path, "r") as f:
                content = f.read()

            lines = content.split("\n")
            for line in lines:
                line = line.strip()

                # skip comments and empty lines
                if "=" not in line or line.startswith("!"):
                    continue

                # remove inline comments
                if "!" in line:
                    line = line.split("!")[0].strip()

                # look for the stuff we care about
                if "colors_results_directory" in line:
                    value = line.split("=")[1].strip().strip("'\"")
                    params["sed_directory"] = f"../{value}/"

                elif "distance" in line and "distance" == line.split("=")[0].strip():
                    value = (
                        line.split("=")[1].strip().replace("d", "e")
                    )  # fix fortran notation
                    try:
                        params["distance"] = float(value)
                    except ValueError:
                        pass

                elif "initial_mass" in line:
                    value = line.split("=")[1].strip().replace("d", "e")
                    try:
                        params["initial_mass"] = float(value)
                    except ValueError:
                        pass

                elif "initial_z" in line:
                    value = line.split("=")[1].strip().replace("d", "e")
                    try:
                        params["initial_z"] = float(value)
                    except ValueError:
                        pass

            return params

        except Exception as e:
            print(f"Warning: Could not parse inlist file {inlist_path}: {e}")
            return params

    def setup_plot(self):
        """
        Set up the plot with proper labels, scales, and styling
        """
        # basic labels and formatting
        self.ax.set_xlabel(r"Wavelength ($\AA$)", fontsize=14)
        self.ax.set_ylabel(r"Flux (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=14)
        self.ax.set_xscale("log")  # log scale makes way more sense for wavelengths
        self.ax.yaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))

        # add those EM spectrum colored regions
        self.add_em_spectrum_regions()

        # put stellar info in a nice box in the corner
        stellar_info = self.format_stellar_info()
        if stellar_info:
            self.ax.text(
                0.99,
                0.98,
                stellar_info,
                transform=self.ax.transAxes,
                fontsize=10,
                verticalalignment="top",
                horizontalalignment="right",
                bbox=dict(
                    boxstyle="round,pad=0.4",
                    facecolor="white",
                    edgecolor="lightblue",
                    alpha=0.3,
                    linewidth=1.5,
                ),
            )

        # apply manual limits if they were set
        if self.xlim:
            self.ax.set_xlim(self.xlim)
        if self.ylim:
            self.ax.set_ylim(self.ylim)

        plt.tight_layout()

    def format_stellar_info(self):
        """
        Format the stellar parameters into a nice info box string
        """
        info_lines = []

        if self.stellar_params["initial_mass"] is not None:
            info_lines.append(f"Mass: {self.stellar_params['initial_mass']:.1f} M☉")

        if self.stellar_params["initial_z"] is not None:
            info_lines.append(f"Z: {self.stellar_params['initial_z']:.4f}")

        if self.stellar_params["distance"] is not None:
            # convert from cm to parsecs because who the hell uses cm for distances
            distance_pc = self.stellar_params["distance"] / 3.0857e18
            if distance_pc < 1000:
                info_lines.append(f"Distance: {distance_pc:.1f} pc")
            else:
                info_lines.append(f"Distance: {distance_pc / 1000:.1f} kpc")

        return "\n".join(info_lines)

    def add_em_spectrum_regions(self):
        """
        Add faint colored background regions for different parts of the EM spectrum
        Makes it easier to see what wavelength range you're looking at
        """
        em_regions = [
            (1, 100, "X-ray", "#8B008B"),
            (100, 2000, "FUV", "#4B0082"),
            (2000, 4000, "UV", "#9400D3"),
            (4000, 5000, "Blue", "#0000FF"),
            (5000, 5700, "Green", "#008000"),
            (5700, 7500, "Red", "#FF0000"),
            (7500, 25000, "Near-IR", "#8B0000"),
            (25000, 400000, "Mid-IR", "#A0522D"),
            (400000, 1e7, "Far-IR", "#D2691E"),
            (1e7, 1e12, "Radio", "#696969"),
        ]

        # only add regions that are actually visible
        xlim = self.ax.get_xlim()

        for min_wave, max_wave, label, color in em_regions:
            # skip if region is outside current view
            if max_wave < xlim[0] or min_wave > xlim[1]:
                continue

            # clip to visible range
            visible_min = max(min_wave, xlim[0])
            visible_max = min(max_wave, xlim[1])

            # add the faint colored background
            self.ax.axvspan(visible_min, visible_max, alpha=0.05, color=color, zorder=0)

            # add label at the top
            center_wave = np.sqrt(
                visible_min * visible_max
            )  # geometric mean for log scale
            self.ax.text(
                center_wave,
                0.95,
                label,
                transform=self.ax.get_xaxis_transform(),
                horizontalalignment="center",
                verticalalignment="top",
                fontsize=8,
                color=color,
                alpha=0.1,
                weight="bold",
            )

    # ============================================================================
    # FILE HANDLING AND MONITORING
    # ============================================================================

    def find_output_files(self):
        """
        Find all the SED CSV files in the directory
        """
        return [
            f
            for f in os.listdir(self.directory)
            if "SED.csv" in f and os.path.isfile(os.path.join(self.directory, f))
        ]

    def check_for_changes(self):
        """
        Check if any of the CSV files have been modified since last check
        Basic file watching functionality
        """
        try:
            output_files = self.find_output_files()
            changed = False

            # check if file list changed
            if set(output_files) != set(self.file_timestamps.keys()):
                changed = True

            # check modification times
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

    def read_csv_file(self, file_path):
        """
        Read CSV file without pandas to keep dependencies minimal

        Handles the usual CSV parsing headaches like whitespace and missing values
        """
        try:
            with open(file_path, "r") as f:
                reader = csv.DictReader(f)

                # clean up column names
                fieldnames = [field.strip() for field in reader.fieldnames]

                data = {}
                for field in fieldnames:
                    data[field] = []

                for row in reader:
                    # skip rows with missing values
                    if any(not value.strip() for value in row.values()):
                        continue

                    for field in fieldnames:
                        try:
                            data[field].append(float(row[field.strip()]))
                        except (ValueError, KeyError):
                            data[field].append(np.nan)

                # get rid of columns that are all NaN
                for field in list(data.keys()):
                    if all(np.isnan(x) for x in data[field]):
                        del data[field]

                return data
        except Exception:
            return {}

    # ============================================================================
    # DATA PROCESSING UTILITIES
    # ============================================================================

    def ensure_numeric(self, data):
        """
        Convert data to numeric array or return zeros if conversion fails
        Because you never know what weird stuff ends up in CSV files
        """
        if data is None or len(data) == 0:
            return np.array([])
        try:
            if isinstance(data[0], str):
                data = np.array([float(x) for x in data])
            return np.array(data, dtype=float)
        except (ValueError, TypeError):
            return np.zeros(len(data))

    def normalize(self, data):
        """
        Normalize data to 0-1 range
        Not sure if this is actually used anywhere but keeping it just in case
        """
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

    def calculate_data_range(self, data_list):
        """
        Find the overall min/max values across multiple datasets
        Useful for setting consistent y-axis limits
        """
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

    # ============================================================================
    # MAIN PLOTTING FUNCTIONALITY
    # ============================================================================

    def update_plot(self, frame):
        """
        Main plotting function that gets called every refresh interval

        This is where all the magic happens - reads the CSV files and makes the plot
        """
        # handle video frame counting
        if self.save_video:
            self.frame_counter += 1
            if self.frame_counter >= self.video_frame_count:
                plt.close(self.fig)
                return

        # only update if files changed (unless we're making a video)
        if not self.check_for_changes() and not self.save_video:
            return

        # save current axis limits before clearing
        current_xlim = self.ax.get_xlim() if not self.xlim else self.xlim
        current_ylim = self.ax.get_ylim() if not self.ylim else self.ylim

        # clear the plot and reset flags
        self.ax.clear()
        self.full_sed_plotted = False

        # find all the CSV files
        output_files = self.find_output_files()
        if not output_files:
            self.setup_plot()
            self.ax.set_xlim(current_xlim)
            self.ax.set_ylim(current_ylim)
            return

        # collect all data for scaling
        all_flux_data = []
        all_convolved_flux_data = []
        sed_wavelengths = []  # for auto x-axis scaling
        sed_flux_data = []

        # set up filter names and colors - gaia filters get physically meaningful colors
        legend_labels = {
            "Gbp": r"Gaia G$_{\mathrm{BP}}$",
            "Grp": r"Gaia G$_{\mathrm{RP}}$",
            "G": "Gaia G",
            "Grvs": r"Gaia G$_{\mathrm{RVS}}$",
            "VEGA": "VEGA",
        }

        filter_colors = {
            "Gbp": "#0000FF",  # blue for blue photometer
            "Grp": "#FF0000",  # red for red photometer
            "G": "#228B22",  # green for broad band
            "Grvs": "#9370DB",  # purple
            "VEGA": "#DDA0DD",  # light purple
        }

        # keep track of what we've added to the legend
        displayed_labels = set()
        legend_handles = []
        legend_labels_list = []

        # process each CSV file
        for file_path in output_files:
            if (
                "VEGA" not in file_path
                and "bright" not in file_path
                and "faint" not in file_path
            ):
                linestyle = "--" if "VEGA" in file_path else "-"

                try:
                    data_dict = self.read_csv_file(
                        os.path.join(self.directory, file_path)
                    )

                    if not data_dict:
                        continue

                    # extract the data columns
                    wavelengths = self.ensure_numeric(data_dict.get("wavelengths", []))
                    flux = self.ensure_numeric(data_dict.get("fluxes", []))
                    convolved_flux = self.ensure_numeric(
                        data_dict.get("convolved_flux", [])
                    )

                    # collect data for range calculations
                    if len(flux) > 0:
                        all_flux_data.append(flux)
                        if len(wavelengths) > 0:
                            sed_wavelengths.extend(wavelengths)
                            sed_flux_data.extend(flux)

                    if len(convolved_flux) > 0:
                        all_convolved_flux_data.append(convolved_flux)

                    # plot the full SED (black line) only once
                    if (
                        len(wavelengths) > 0
                        and len(flux) > 0
                        and not self.full_sed_plotted
                    ):
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

                    # plot the convolved flux with proper colors and labels
                    if len(wavelengths) > 0 and len(convolved_flux) > 0:
                        # figure out what filter this is
                        short_label = next(
                            (v for k, v in legend_labels.items() if k in file_path),
                            None,
                        )
                        filter_key = next(
                            (k for k in legend_labels.keys() if k in file_path), None
                        )
                        color = (
                            filter_colors.get(filter_key, None) if filter_key else None
                        )

                        # add to legend if we haven't seen this filter yet
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
                            # plot without adding to legend (duplicate filter)
                            self.ax.plot(
                                wavelengths,
                                convolved_flux,
                                linewidth=1,
                                linestyle=linestyle,
                                color=color,
                            )

                except Exception:
                    continue  # skip files that can't be read

        # add the legend
        if legend_handles and legend_labels_list:
            self.ax.legend(
                legend_handles, legend_labels_list, loc="lower right", fontsize=10
            )

        # redo the plot setup (labels, EM regions, etc.)
        self.setup_plot()

        # auto-scale x-axis based on where the SED has significant flux
        if not self.xlim and sed_wavelengths and sed_flux_data:
            sed_waves = np.array(sed_wavelengths)
            sed_fluxes = np.array(sed_flux_data)

            # find wavelengths where flux is above 1% of maximum
            max_sed_flux = np.max(sed_fluxes)
            if max_sed_flux > 0:
                mask = sed_fluxes > 0.01 * max_sed_flux
                active_wavelengths = sed_waves[mask]

                if len(active_wavelengths) > 0:
                    min_wave = np.min(active_wavelengths)
                    max_wave = np.max(active_wavelengths)

                    # add some padding so things don't touch the edges
                    wave_range = max_wave - min_wave
                    padding = wave_range * 0.01

                    self.ax.set_xlim(
                        [min_wave - padding * 0.01, max_wave + padding * 10]
                    )

        # set y-axis limits based on all the data
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

        # override with manual limits if they were set
        if self.xlim:
            self.ax.set_xlim(self.xlim)

        # redraw the canvas
        self.fig.canvas.draw_idle()

    # ============================================================================
    # MAIN RUN METHOD
    # ============================================================================

    def run(self):
        """
        Start the monitoring and animation

        Either shows live plot or saves a video depending on settings
        """
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


# ============================================================================
# MAIN FUNCTION
# ============================================================================


def main():
    """
    Main function - set up the checker with default params and run it
    """
    checker = SEDChecker(
        directory="../SED/",
        xlim=None,  # auto-scale x-axis to SED data range
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
