import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def find_output_files(directory):
    """Finds all files in the directory that contain 'SED.csv' in their name."""
    return [
        f
        for f in os.listdir(directory)
        if "SED.csv" in f and os.path.isfile(os.path.join(directory, f))
    ]


def normalize(data):
    if data is None or len(data) == 0:
        return np.array([])  # Return empty array if data is None or empty
    return (
        (data - np.min(data)) / (np.max(data) - np.min(data))
        if np.max(data) != np.min(data)
        else data
    )


def process_all_files(directory, xlim=None):
    output_files = find_output_files(directory)
    if not output_files:
        print("No output files found.")
        return

    plt.figure(figsize=(12, 6))

    full_sed_plotted = False

    for file_path in output_files:
        if "VEGA" in file_path:
            linestyle = "--"
        else:
            linestyle = "-"

        df = (
            pd.read_csv(os.path.join(directory, file_path), delimiter=",", header=0)
            .rename(columns=str.strip)
            .dropna()
        )

        # Extract columns safely
        wavelengths = df.get("wavelengths", pd.Series()).to_numpy()
        flux = df.get("fluxes", pd.Series()).to_numpy()
        convolved_flux = df.get("convolved_flux", pd.Series()).to_numpy()
        filter_wavelengths = df.get("filter_wavelengths", pd.Series()).to_numpy()
        filter_trans = df.get("filter_trans", pd.Series()).to_numpy()

        print(f"Processing {file_path}")
        print("Wavelengths shape:", wavelengths.shape)
        print("Flux shape:", flux.shape)
        print("Convolved Flux shape:", convolved_flux.shape)

        # Normalize data safely
        nflux = normalize(flux)
        nconvolved_flux = normalize(convolved_flux)
        nfilter_trans = normalize(filter_trans)

        # Plot full SED only once (assume it is the same across files)
        if len(wavelengths) > 0 and len(flux) > 0 and not full_sed_plotted:
            plt.plot(
                wavelengths,
                flux,
                label="Full SED (common)",
                color="black",
                linewidth=1.5,
                linestyle=linestyle,
            )
            full_sed_plotted = True

        # Plot convolved SED if available
        if len(wavelengths) > 0 and len(convolved_flux) > 0:
            plt.plot(
                wavelengths,
                convolved_flux,
                label=f"Convolved SED ({file_path})",
                linewidth=1,
                linestyle=linestyle,
            )

        # Plot normalized filters if available
        if len(filter_wavelengths) > 0 and len(nfilter_trans) > 0:
            pass  # plt.plot(filter_wavelengths, filter_trans, label=f"Filter ({file_path})", linewidth=1, linestyle=":")

    # Formatting
    plt.xlabel("Wavelengths")
    plt.ylabel("Flux")
    plt.title("Combined SEDs and Filters")
    plt.ticklabel_format(style="plain", useOffset=False)

    # Add legend to differentiate curves
    plt.legend(loc="best", fontsize=8)

    # Set x-axis limits if provided
    if xlim:
        plt.xlim(xlim)

    plt.tight_layout()
    plt.show()


def main():
    directory = "../LOGS/SED/"  # Change if needed
    process_all_files(directory, xlim=[0, 60000])


if __name__ == "__main__":
    main()
