import pandas as pd
import matplotlib.pyplot as plt

# Load the data
file_path = "LOGS/history.data"

# Define column names based on the sample provided (modify as necessary for exact column names)
columns = [
    "model_number", "num_zones", "star_age", "log_dt", "star_mass", "bc_H", "bc_J_bb", "bc_J",
    "abs_mag_H", "abs_mag_J", "abs_mag_V", "abs_mag_B", "abs_mag_U", "Teff", "log_Teff", "log_L",
    "log_R", "log_g", "num_retries", "num_iters", "av_v"
]

# Read data into DataFrame, skipping any potential metadata rows
data = pd.read_csv(file_path, delim_whitespace=True, skiprows=1, names=columns)

# Convert columns to numeric, coerce errors to NaN, and drop rows with NaN values
data = data.apply(pd.to_numeric, errors='coerce').dropna()

# Define which magnitudes to plot
magnitudes = ["abs_mag_H", "abs_mag_J", "abs_mag_V", "abs_mag_B", "abs_mag_U"]

# Define parameters to plot magnitudes against
parameters = ["star_age", "Teff", "log_R", "log_L"]

# Plot magnitudes against each parameter
for param in parameters:
    plt.figure(figsize=(10, 6))
    for mag in magnitudes:
        plt.plot(data[param], data[mag], label=mag)
    
    plt.xlabel(param)
    plt.ylabel("Magnitude")
    plt.title(f"Magnitudes vs {param}")
    plt.legend()
    plt.grid()
    plt.show()

# Define color indices for CMDs
color_indices = {
    "B-V": ("abs_mag_B", "abs_mag_V"),
    "V-J": ("abs_mag_V", "abs_mag_J"),
    "V-H": ("abs_mag_V", "abs_mag_H"),
    "U-B": ("abs_mag_U", "abs_mag_B")
}

# Plot color-magnitude diagrams
for color, (mag1, mag2) in color_indices.items():
    color_index = data[mag1] - data[mag2]  # Calculate color index
    plt.figure(figsize=(10, 6))
    plt.scatter(color_index, data[mag1], label=f"{mag1} vs {color} color index", s=10)
    
    plt.xlabel(f"{color} Color Index ({mag1} - {mag2})")
    plt.ylabel(mag1)
    plt.gca().invert_yaxis()  # Invert y-axis for magnitudes
    plt.title(f"Color-Magnitude Diagram: {mag1} vs {color}")
    plt.legend()
    plt.grid()
    plt.show()

