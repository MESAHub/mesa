import os
import time
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.animation import FuncAnimation, FFMpegWriter


class SEDChecker:
    def __init__(self, directory="../LOGS/SED/", xlim=None, ylim=None, refresh_interval=2, 
                 save_video=False, video_filename="sed_animation.mp4", video_fps=10, video_duration=30):
        """
        Initialize the SED checker with auto-refresh capability.
        
        Args:
            directory: Directory containing SED CSV files
            xlim: Optional limits for x-axis [min, max]
            ylim: Optional limits for y-axis [min, max]
            refresh_interval: Time in seconds between refresh attempts
            save_video: Whether to save the animation as a video
            video_filename: Filename for the saved video
            video_fps: Frames per second for the video
            video_duration: Duration of the video in seconds
        """
        self.directory = directory
        self.xlim = xlim
        self.ylim = ylim
        self.refresh_interval = refresh_interval
        self.file_timestamps = {}
        self.fig = None
        self.full_sed_plotted = False
        
        # Video saving settings
        self.save_video = save_video
        self.video_filename = video_filename
        self.video_fps = video_fps
        self.video_duration = video_duration
        self.video_frame_count = int(video_fps * video_duration)
        self.frame_counter = 0
        
        # Create the figure and axis
        self.fig, self.ax = plt.subplots(figsize=(12, 6))
        self.setup_plot()
    
    def setup_plot(self):
        """Set up the plot with formatting and labels."""
        self.ax.set_xlabel("Wavelengths")
        self.ax.set_ylabel("Flux")
        self.ax.set_title("Combined SEDs and Filters")
        self.ax.ticklabel_format(style="plain", useOffset=False)
        
        # Set fixed axis limits if provided
        if self.xlim:
            self.ax.set_xlim(self.xlim)
        if self.ylim:
            self.ax.set_ylim(self.ylim)
            
        plt.tight_layout()
    
    def find_output_files(self):
        """Finds all files in the directory that contain 'SED.csv' in their name."""
        return [
            f
            for f in os.listdir(self.directory)
            if "SED.csv" in f and os.path.isfile(os.path.join(self.directory, f))
        ]
    
    def normalize(self, data):
        """Normalize data to the range [0, 1] with type safety."""
        if data is None or len(data) == 0:
            return np.array([])  # Return empty array if data is None or empty
        
        # Handle string data by converting to float
        try:
            # If data contains strings, convert to float
            if isinstance(data[0], str):
                data = np.array([float(x) for x in data])
            # Ensure data is numeric (could be object dtype)
            data = data.astype(float)
        except (ValueError, TypeError) as e:
            print(f"Warning: Data conversion error in normalize: {e}")
            return np.zeros_like(data, dtype=float)
        
        # Now do the normalization with numeric data
        try:
            data_min = np.min(data)
            data_max = np.max(data)
            if data_max != data_min:
                return (data - data_min) / (data_max - data_min)
            else:
                return np.zeros_like(data)
        except Exception as e:
            print(f"Warning: Normalization error: {e}")
            return np.zeros_like(data)
    
    def check_for_changes(self):
        """Check if any files have been modified since last check."""
        try:
            output_files = self.find_output_files()
            changed = False
            
            # Check for new files
            if set(output_files) != set(self.file_timestamps.keys()):
                changed = True
            
            # Check for modified files
            for file_path in output_files:
                full_path = os.path.join(self.directory, file_path)
                current_mtime = os.path.getmtime(full_path)
                
                if (file_path not in self.file_timestamps or 
                    self.file_timestamps[file_path] != current_mtime):
                    self.file_timestamps[file_path] = current_mtime
                    changed = True
            
            return changed
        except Exception as e:
            print(f"Error checking for changes: {e}")
            return False
    
    def calculate_data_range(self, data_list):
        """Calculate the min and max values across all datasets for axis scaling."""
        all_min = float('inf')
        all_max = float('-inf')
        
        for data in data_list:
            if data is None or len(data) == 0:
                continue
                
            try:
                # Convert to numeric if needed
                if isinstance(data[0], str):
                    data = np.array([float(x) for x in data])
                data = np.array(data, dtype=float)
                
                data_min = np.min(data)
                data_max = np.max(data)
                all_min = min(all_min, data_min)
                all_max = max(all_max, data_max)
            except Exception as e:
                print(f"Warning: Error in calculate_data_range: {e}")
                continue
                
        if all_min == float('inf') or all_max == float('-inf'):
            return None, None
            
        return all_min, all_max
    
    def ensure_numeric(self, data):
        """Ensure data is numeric, converting strings if needed."""
        if data is None or len(data) == 0:
            return np.array([])
            
        try:
            # If data contains strings, convert to float
            if isinstance(data[0], str):
                data = np.array([float(x) for x in data])
            # Ensure data is numeric
            return np.array(data, dtype=float)
        except (ValueError, TypeError) as e:
            print(f"Warning: Data conversion error: {e}")
            return np.zeros(len(data))
    
    def update_plot(self, frame):
        """Update function for animation - processes files and updates plot."""
        if self.save_video:
            self.frame_counter += 1
            # Print progress periodically
            if self.frame_counter % 10 == 0:
                print(f"Creating video: Frame {self.frame_counter}/{self.video_frame_count}")
            
            # Stop the animation when we've reached the desired frame count
            if self.frame_counter >= self.video_frame_count:
                plt.close(self.fig)
                return
        
        if not self.check_for_changes() and not self.save_video:
            return
        
        # Clear the plot but preserve limits
        current_xlim = self.ax.get_xlim() if not self.xlim else self.xlim
        current_ylim = self.ax.get_ylim() if not self.ylim else self.ylim
        
        self.ax.clear()
        self.full_sed_plotted = False
        
        # Process all files
        output_files = self.find_output_files()
        if not output_files:
            print("No output files found.")
            # Restore axis setup
            self.setup_plot()
            if not self.xlim:
                self.ax.set_xlim(current_xlim)
            if not self.ylim:
                self.ax.set_ylim(current_ylim)
            return
        
        # Lists to store all flux data for auto-scaling if needed
        all_flux_data = []
        all_convolved_flux_data = []
        
        for file_path in output_files:
            if "VEGA" in file_path:
                linestyle = "--"
            else:
                linestyle = "-"
            
            try:
                # Read the CSV file
                df = pd.read_csv(os.path.join(self.directory, file_path), delimiter=",", header=0)
                df = df.rename(columns=str.strip).dropna()
                
                # Extract columns safely
                wavelengths = None
                flux = None
                convolved_flux = None
                filter_wavelengths = None
                filter_trans = None
                
                if "wavelengths" in df.columns:
                    wavelengths = self.ensure_numeric(df["wavelengths"].values)
                
                if "fluxes" in df.columns:
                    flux = self.ensure_numeric(df["fluxes"].values)
                    if len(flux) > 0:
                        all_flux_data.append(flux)
                
                if "convolved_flux" in df.columns:
                    convolved_flux = self.ensure_numeric(df["convolved_flux"].values)
                    if len(convolved_flux) > 0:
                        all_convolved_flux_data.append(convolved_flux)
                
                if "filter_wavelengths" in df.columns:
                    filter_wavelengths = self.ensure_numeric(df["filter_wavelengths"].values)
                
                if "filter_trans" in df.columns:
                    filter_trans = self.ensure_numeric(df["filter_trans"].values)
                
                # Plot full SED only once (assume it is the same across files)
                if wavelengths is not None and len(wavelengths) > 0 and flux is not None and len(flux) > 0 and not self.full_sed_plotted:
                    self.ax.plot(
                        wavelengths,
                        flux,
                        label="Full SED (common)",
                        color="black",
                        linewidth=1.5,
                        linestyle=linestyle,
                    )
                    self.full_sed_plotted = True
                
                # Plot convolved SED if available
                if wavelengths is not None and len(wavelengths) > 0 and convolved_flux is not None and len(convolved_flux) > 0:
                    self.ax.plot(
                        wavelengths,
                        convolved_flux,
                        label=f"Convolved SED ({file_path})",
                        linewidth=1,
                        linestyle=linestyle,
                    )
                
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
                continue
        
        # Add legend
        self.ax.legend(loc="best", fontsize=8)
        
        # Apply axis formatting
        self.setup_plot()
        
        # Handle automatic y-axis scaling if ylim is not explicitly set
        if not self.ylim and all_convolved_flux_data:
            # Calculate data ranges
            min_val, max_val = self.calculate_data_range(all_convolved_flux_data + all_flux_data)
            if min_val is not None and max_val is not None:
                # Add 10% padding
                padding = (max_val - min_val) * 0.1
                self.ax.set_ylim([min_val - padding, max_val + padding])
            else:
                # Restore previous ylim if we can't calculate new ones
                self.ax.set_ylim(current_ylim)
        elif not self.ylim:
            # Restore previous ylim if no data available
            self.ax.set_ylim(current_ylim)
        
        # Ensure x limits are maintained
        if not self.xlim:
            self.ax.set_xlim(current_xlim)
        
        # Update the plot
        self.fig.canvas.draw_idle()
    
    def run(self):
        """Start the auto-refreshing display."""
        print(f"Monitoring directory: {self.directory}")
        print(f"Refresh interval: {self.refresh_interval} seconds")
        
        # Set up animation
        if self.save_video:
            print(f"Saving video to {self.video_filename}")
            print(f"Video settings: {self.video_fps} fps, {self.video_duration} seconds")
            
            try:
                # Create writer
                writer = FFMpegWriter(fps=self.video_fps, metadata=dict(title='SED Animation', 
                                                                        artist='MESA Colors'), 
                                    bitrate=5000)
                
                # Set up animation with save
                self.animation = FuncAnimation(
                    self.fig, 
                    self.update_plot,
                    frames=self.video_frame_count,
                    interval=self.refresh_interval * 1000,  # Convert to milliseconds
                    cache_frame_data=False
                )
                
                # Save animation
                self.animation.save(self.video_filename, writer=writer)
                print(f"Video saved to {self.video_filename}")
            except Exception as e:
                print(f"Error saving video: {e}")
        else:
            try:
                # Normal interactive mode
                self.animation = FuncAnimation(
                    self.fig, 
                    self.update_plot,
                    interval=self.refresh_interval * 1000,  # Convert to milliseconds
                    cache_frame_data=False
                )
                
                # Show the plot (this will block until window is closed)
                plt.show()
            except KeyboardInterrupt:
                print("\nMonitoring stopped by user.")
            except Exception as e:
                print(f"Error in animation: {e}")


def main():
    # Create and run the checker with default settings
    # All parameters are optional
    checker = SEDChecker(
        directory="../LOGS/SED/",
        xlim=[0, 20000],
        ylim=None,
        refresh_interval=5,
        save_video=False,
        video_filename="sed_animation.mp4",
        video_fps=10,
        video_duration=30
    )
    checker.run()


if __name__ == "__main__":
    main()