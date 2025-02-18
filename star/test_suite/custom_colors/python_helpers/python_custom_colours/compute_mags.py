import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import os
import sys
from scipy.interpolate import griddata
import argparse

# Physical constants (cgs units)
h = 6.62607015e-27  # Planck constant (erg·s)
c = 2.99792458e10  # Speed of light (cm/s)
k_B = 1.380649e-16  # Boltzmann constant (erg/K)
R_sun = 6.957e10  # Solar radius in cm
pc_in_cm = 3.0857e18  # Parsec in cm


def read_mesa_history(filename):
    """
    Reads the MESA history.data file into a pandas DataFrame.

    Parameters:
        filename (str): Path to the history.data file.

    Returns:
        pd.DataFrame: DataFrame containing the history data.
    """
    try:
        with open(filename, "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: History file '{filename}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading history file: {e}")
        sys.exit(1)

    # Find the index where the column headers start
    header_line_idx = None
    for idx, line in enumerate(lines):
        if line.strip().startswith("model_number"):
            header_line_idx = idx
            break
    if header_line_idx is None:
        print("Error: Header line with column names not found in history file.")
        sys.exit(1)

    # Extract headers
    headers = lines[header_line_idx].strip().split()

    # Read data into DataFrame
    try:
        data = pd.read_csv(
            filename,
            skiprows=header_line_idx + 1,
            sep="\s+",  # Use regex to handle any whitespace
            names=headers,
            comment="#",
            engine="python",
        )
    except Exception as e:
        print(f"Error parsing history file: {e}")
        sys.exit(1)

    return data


def read_filter_data(filter_file):
    """
    Reads the filter transmission data from a file.

    Parameters:
        filter_file (str): Path to the filter data file.

    Returns:
        tuple: Tuple containing:
            - filter_wavelength_angstrom (np.ndarray): Wavelength array in Angstroms.
            - filter_transmission (np.ndarray): Transmission values (0 to 1).
    """
    if not os.path.isfile(filter_file):
        print(f"Error: Filter file '{filter_file}' does not exist.")
        sys.exit(1)

    try:
        filter_data = np.loadtxt(filter_file, comments="#", delimiter=",", skiprows=1)
    except Exception as e:
        print(f"Error loading filter file '{filter_file}': {e}")
        sys.exit(1)

    if filter_data.ndim != 2 or filter_data.shape[1] < 2:
        print(f"Error: Filter file '{filter_file}' has an unexpected format.")
        sys.exit(1)

    filter_wavelength_units = filter_data[:, 0]
    filter_transmission = filter_data[:, 1]

    # Assuming filter wavelengths are in Angstroms. If not, adjust accordingly
    # Uncomment the appropriate line based on your filter file's wavelength units.

    # If wavelengths are in microns (µm):
    # filter_wavelength_angstrom = filter_wavelength_units * 1e4  # µm to Å

    # If wavelengths are in nanometers (nm):
    # filter_wavelength_angstrom = filter_wavelength_units * 10    # nm to Å

    # If wavelengths are already in Angstroms (Å):
    filter_wavelength_angstrom = filter_wavelength_units

    return filter_wavelength_angstrom, filter_transmission


def planck(Teff, wavelength):
    """
    Calculate the spectral radiance using Planck's Law.

    Parameters:
    Teff (float): Effective temperature in Kelvin (K)
    wavelength (float or np.ndarray): Wavelength in meters (m)

    Returns:
    float or np.ndarray: Spectral radiance (flux) in W/m²/m/sr
    """
    # Physical constants in SI units
    h = 6.62607015e-34  # Planck's constant (J·s)
    c = 2.99792458e8  # Speed of light in vacuum (m/s)
    k = 1.380649e-23  # Boltzmann's constant (J/K)

    # Ensure wavelength is a NumPy array for vectorized operations
    wavelength = np.array(wavelength, dtype=np.float64)

    # Avoid division by zero or negative wavelengths
    wavelength = np.where(wavelength > 0, wavelength, np.nan)

    # Calculate the exponent: (h * c) / (wavelength * k * Teff)
    exponent = (h * c) / (wavelength * k * Teff)

    # To prevent overflow in the exponential, cap the exponent at a maximum value
    # np.exp(700) is a safe upper limit to prevent overflow (np.exp(709) ~ 8.218407e+307)
    exponent = np.clip(exponent, a_min=None, a_max=700)

    # Calculate the spectral radiance using Planck's Law
    # Handle cases where exponent is very large by setting flux to 0
    with np.errstate(over="ignore", invalid="ignore"):
        exp_term = np.exp(exponent)
        flux = (2.0 * h * c**2) / (wavelength**5) / (exp_term - 1.0)
        flux = np.where(np.isfinite(flux), flux, 0.0)  # Replace inf and nan with 0.0

    # If the input was a scalar, return a scalar
    if np.isscalar(wavelength):
        return float(flux)
    return flux


def obs_flux(flux_lambda, radius, distance):
    """
    Calculate the observed flux at a distance from the blackbody.

    Parameters:
    flux_lambda float: Spectral radiance (flux) in W/m²/m/sr
    radius (float): Radius of the star in meters (m)
    distance (float): Distance to the observer in meters (m)

    Returns:
    float: Observed flux in W/m²/m
    """
    # Calculate the observed flux using the formula:
    # F_obs(λ) = π * B(λ, T) * (R / d)^2
    flux_obs = np.pi * flux_lambda * (radius / distance) ** 2

    return flux_obs


def compute_filter_magnitude(
    T_eff, log_g, log_R, filter_wavelength_angstrom, filter_transmission
):
    """
    Computes the absolute AB magnitude in a specified filter by convolving the stellar spectrum with the filter transmission.

    Parameters:
        T_eff (float): Effective temperature in Kelvin.
        log_g (float): Logarithm of surface gravity.
        log_R (float): Logarithm (base 10) of stellar radius relative to the Sun.
        filter_wavelength_angstrom (np.ndarray): Wavelength array of the filter in Angstroms.
        filter_transmission (np.ndarray): Transmission values of the filter (0 to 1).

    Returns:
        float: Absolute AB magnitude in the filter.
    """
    # Define wavelength range based on filter
    print("teff", T_eff)
    print("logg", log_g)
    print("logr", log_R)
    wavelength_min = np.min(filter_wavelength_angstrom)
    wavelength_max = np.max(filter_wavelength_angstrom)
    wavelength_length = int(len(filter_wavelength_angstrom) * 1.5)

    # Generate wavelength array around filter range
    wavelength_angstrom = np.linspace(
        wavelength_min * 0.9, wavelength_max * 1.1, wavelength_length
    )
    wavelength_m = wavelength_angstrom * 1e-10
    filter_wavelength_m = filter_wavelength_angstrom * 1e-10
    print(
        "Generated Wavelengths (middle 5):",
        wavelength_angstrom[
            int(len(wavelength_angstrom) / 2) : int(len(wavelength_angstrom) / 2) + 5
        ],
    )
    print(
        "Generated Wavelengths (middle 5):",
        wavelength_m[int(len(wavelength_m) / 2) : int(len(wavelength_m) / 2) + 5],
    )
    print(
        "Filter Wavelengths (middle 5):",
        filter_wavelength_angstrom[
            int(len(filter_wavelength_angstrom) / 2) : int(
                len(filter_wavelength_angstrom) / 2
            )
            + 5
        ],
    )

    # Calculate flux_lambda using Planck function
    flux_lambda = planck(wavelength_m, T_eff)  # erg/s/cm²/Å
    print(
        "Flux Lambda (middle 5):",
        flux_lambda[int(len(flux_lambda) / 2) : int(len(flux_lambda) / 2) + 5],
    )

    distance = 3.086e17
    radii = 10**log_R * 6.957e8

    observed_flux_lambda = obs_flux(flux_lambda, radii, distance)
    print(
        "Observed Flux (middle 5):",
        observed_flux_lambda[
            int(len(observed_flux_lambda) / 2) : int(len(observed_flux_lambda) / 2) + 5
        ],
    )

    filter_wavelength_m = filter_wavelength_angstrom * 1e-10

    # Interpolate filter transmission to match wavelength_m grid
    filter_interp_func = interp1d(
        filter_wavelength_m, filter_transmission, bounds_error=False, fill_value=0
    )
    filter_transmission_interp = filter_interp_func(wavelength_m)
    print(
        "Filter Transmission Interpolated (middle 10):",
        filter_transmission_interp[
            len(filter_transmission_interp) // 2
            - 5 : len(filter_transmission_interp) // 2
            + 5
        ],
    )

    # Compute f_nu
    # f_nu = ∫ f_lambda * T(lambda) * (c / lambda^2) d lambda / ∫ T(lambda) * (c / lambda^2) d lambda
    numerator = np.trapz(
        observed_flux_lambda * filter_transmission_interp * c / wavelength_m**2,
        wavelength_m,
    )
    denominator = np.trapz(
        filter_transmission_interp * c / wavelength_m**2, wavelength_m
    )
    f_nu = numerator / denominator  # W/m²/Hz

    print("f_nu (W/m²/Hz):", f_nu)

    # Convert f_nu to cgs units (erg/s/cm²/Hz)
    f_nu_cgs = f_nu * 1e3  # 1 W/m² = 1e3 erg/s/cm²

    print("f_nu (erg/s/cm²/Hz):", f_nu_cgs)

    # Compute AB magnitude
    m_AB = -2.5 * np.log10(f_nu_cgs) - 48.60
    print("AB Magnitude:", m_AB)

    return (
        m_AB,
        wavelength_angstrom,
        flux_lambda,
        filter_transmission,
        observed_flux_lambda,
    )


def read_color_file(color_file):
    """
    Reads the color file into a pandas DataFrame.

    Parameters:
        color_file (str): Path to the color file.

    Returns:
        pd.DataFrame: DataFrame containing the color file data.
    """
    try:
        # Assuming columns are space-separated and there's a header to skip
        data = pd.read_csv(color_file, delim_whitespace=True, comment="#", header=None)
        # Rename columns based on your example
        data.columns = [
            "Teff",
            "logg",
            "M_div_H",
            "U",
            "B",
            "V",
            "R",
            "I",
            "J",
            "H",
            "K",
            "L",
            "Lprime",
            "M",
        ]
    except Exception as e:
        print(f"Error reading color file '{color_file}': {e}")
        sys.exit(1)

    return data


def compute_bolometric_magnitude_from_color_file(T_eff, log_g, color_data, band="V"):
    """
    Computes the bolometric correction using the color file data.

    Parameters:
        T_eff (float): Effective temperature in Kelvin.
        log_g (float): Logarithm of surface gravity.
        color_data (pd.DataFrame): DataFrame containing the color file data.
        band (str): Filter band for the bolometric correction (default is 'V').

    Returns:
        float: Bolometric correction.
    """
    if band not in color_data.columns:
        raise ValueError(f"Band '{band}' not found in color file.")

    # Prepare the grid and points for interpolation
    points = color_data[["Teff", "logg"]].values
    values = color_data[band].values

    # Interpolate using griddata
    bolometric_correction = griddata(points, values, (T_eff, log_g), method="linear")

    if np.isnan(bolometric_correction):
        print(
            f"Warning: Bolometric correction not found for T_eff={T_eff}, log_g={log_g}."
        )
        return np.nan

    return bolometric_correction


def generate_modified_history(history, bol_mag, filter_mag, filter_used, filter_file):
    """
    Adds new columns to the history DataFrame.

    Parameters:
        history (pd.DataFrame): Original history DataFrame.
        bol_mag (pd.Series): Bolometric magnitudes.
        filter_mag (pd.Series): Filter magnitudes.
        filter_used (str): Name or path of the filter used.
        filter_file (str): Path to the filter file.

    Returns:
        pd.DataFrame: Modified history DataFrame with new columns.
    """
    history["bol_magnitude"] = bol_mag
    history["filter_magnitude"] = filter_mag
    history["filter_used"] = filter_used
    return history


def write_modified_history(history, output_file):
    """
    Writes the modified history DataFrame to a file, ensuring no duplicate headers
    and formatting floats without scientific notation.

    Parameters:
        history (pd.DataFrame): Modified history DataFrame.
        output_file (str): Path to the output file.
    """
    # print(history)
    try:
        with open(output_file, "w") as f:
            # Write header only once
            f.write(",".join(history.columns) + "\n")
            # Write data without scientific notation
            history.to_csv(f, sep=",", index=False, header=False, float_format="%.6f")
        print(f"Updated history data saved to '{output_file}'.")
    except Exception as e:
        print(f"Error writing to output file '{output_file}': {e}")


def plot_spectrum(
    wavelength_angstrom, flux_lambda, filter_response, observed_flux_lambda, filter_file
):
    """
    Plots the stellar spectrum, filter transmission, and convolved spectrum.

    Parameters:
        wavelength_angstrom (np.ndarray): Wavelength array in Angstroms.
        flux_lambda (np.ndarray): Stellar flux per Angstrom.
        filter_response (np.ndarray): Filter transmission interpolated onto the stellar wavelength grid.
        observed_flux_lambda (np.ndarray): Convolved flux.
        filter_file (str): Path to the filter file used.
    """
    plt.figure(figsize=(10, 6))
    plt.plot(wavelength_angstrom, flux_lambda, label="Stellar Spectrum")
    plt.plot(
        wavelength_angstrom,
        filter_response * np.max(flux_lambda),
        label="Filter Transmission (scaled)",
    )
    plt.plot(wavelength_angstrom, observed_flux_lambda, label="Convolved Spectrum")
    plt.xlabel("Wavelength (Å)")
    plt.ylabel("Flux (erg/s/cm²/Å)")
    plt.legend()
    plt.title(f"Spectrum Convolution with Filter: {os.path.basename(filter_file)}")
    plt.grid(True)
    plt.show()


def main():
    """
    Main function to orchestrate the computation of bolometric and filter magnitudes,
    and to generate a modified history file with the new data.
    """
    default_filter = "data/filters/PAN_STARRS/PAN_STARRS/g.dat"
    # User specifies the filter file
    filter_file = ""  # input("Enter the path to the filter file (e.g., default_filter): ").strip()
    if filter_file == "":
        filter_file = default_filter

    if not os.path.isfile(filter_file):
        print(
            f"Error: Filter file '{filter_file}' does not exist. Please check the path."
        )
        sys.exit(1)

    # Read filter data
    filter_wavelength_angstrom, filter_transmission = read_filter_data(filter_file)

    # Read MESA history data
    history_file = "LOGS/history.data"
    if not os.path.isfile(history_file):
        print(
            f"Error: History file '{history_file}' does not exist. Please check the path."
        )
        sys.exit(1)

    history = read_mesa_history(history_file)

    # Initialize lists to store computed magnitudes
    bol_magnitudes = []
    filter_magnitudes = []

    # Process each model in the history
    total_models = len(history)
    print(f"Processing {total_models} models...")

    for idx, row in history.iterrows():
        T_eff = row.get("Teff")
        log_g = row.get("log_g")
        log_R = row.get("log_R")
        metallicity = row.get("initial_feh", 0.0)

        # Check for missing values
        if pd.isnull(T_eff) or pd.isnull(log_g):
            print(
                f"Warning: Missing T_eff or log_g for model index {idx}. Assigning NaN."
            )
            bol_magnitudes.append(np.nan)
            filter_magnitudes.append(np.nan)
            continue

        # Compute bolometric magnitude
        bol_mag = compute_bolometric_magnitude(T_eff, log_g, metallicity)
        bol_magnitudes.append(bol_mag)

        # Compute filter magnitude
        (
            filter_mag,
            wavelength_angstrom,
            flux_lambda,
            filter_response,
            observed_flux_lambda,
        ) = compute_filter_magnitude(
            T_eff, log_g, log_R, filter_wavelength_angstrom, filter_transmission
        )

        filter_magnitudes.append(filter_mag)
        # exit()
        # Optional: Print progress every 10%
        # if (idx + 1) % max(1, total_models // 10) == 0:
        # print(f"Processed {idx + 1}/{total_models} models.")

    # Add new columns to history
    history = generate_modified_history(
        history,
        pd.Series(bol_magnitudes),
        pd.Series(filter_magnitudes),
        os.path.basename(filter_file),
        filter_file,
    )

    # Write the modified history to a new file
    output_file = "history_with_extras.data"
    write_modified_history(history, output_file)

    # Plotting the last model as an example
    last_model_idx = total_models - 1
    if not np.isnan(filter_magnitudes[last_model_idx]):
        print(f"Generating plot for the last model (index {last_model_idx})...")
        (
            filter_mag,
            wavelength_angstrom,
            flux_lambda,
            filter_response,
            observed_flux_lambda,
        ) = compute_filter_magnitude(
            history.at[last_model_idx, "Teff"],
            history.at[last_model_idx, "log_g"],
            history.at[last_model_idx, "log_R"],
            filter_wavelength_angstrom,
            filter_transmission,
        )
        plot_spectrum(
            wavelength_angstrom,
            flux_lambda,
            filter_response,
            observed_flux_lambda,
            filter_file,
        )
    else:
        print("No valid filter magnitude computed for the last model. Skipping plot.")


import os
import sys
import pandas as pd
import numpy as np


def synth_main():
    """
    Main function to orchestrate the computation of bolometric and filter magnitudes,
    and to handle a CSV file input.
    """
    default_filter = "data/filters/PAN_STARRS/PAN_STARRS/g.dat"
    # User specifies the filter file
    filter_file = ""  # input("Enter the path to the filter file (e.g., default_filter): ").strip()
    if filter_file == "":
        filter_file = default_filter

    if not os.path.isfile(filter_file):
        print(
            f"Error: Filter file '{filter_file}' does not exist. Please check the path."
        )
        sys.exit(1)

    # Read filter data
    filter_wavelength_angstrom, filter_transmission = read_filter_data(filter_file)

    # Read CSV data
    csv_file = "input_data.csv"
    if not os.path.isfile(csv_file):
        print(f"Error: CSV file '{csv_file}' does not exist. Please check the path.")
        sys.exit(1)

    data = pd.read_csv(csv_file)

    # Initialize lists to store computed magnitudes
    bol_magnitudes = []
    filter_magnitudes = []

    # Process each row in the CSV
    total_models = len(data)
    print(f"Processing {total_models} entries...")

    for idx, row in data.iterrows():
        T_eff = row.get("teff")
        log_g = row.get("logg")
        metallicity = row.get("meta", 0.0)
        flux = row.get("flux")
        wavelength = row.get("wavelength")

        # Check for missing values
        if (
            pd.isnull(T_eff)
            or pd.isnull(log_g)
            or pd.isnull(flux)
            or pd.isnull(wavelength)
        ):
            print(f"Warning: Missing data for row index {idx}. Assigning NaN.")
            bol_magnitudes.append(np.nan)
            filter_magnitudes.append(np.nan)
            continue

        # Compute bolometric magnitude
        bol_mag = compute_bolometric_magnitude(T_eff, log_g, metallicity)
        bol_magnitudes.append(bol_mag)

        # Compute filter magnitude
        (
            filter_mag,
            wavelength_angstrom,
            flux_lambda,
            filter_response,
            observed_flux_lambda,
        ) = compute_filter_magnitude(
            T_eff, log_g, None, filter_wavelength_angstrom, filter_transmission
        )

        filter_magnitudes.append(filter_mag)

    # Add new columns to the data
    data["bolometric_magnitude"] = pd.Series(bol_magnitudes)
    data["filter_magnitude"] = pd.Series(filter_magnitudes)

    # Write the modified data to a new file
    output_file = "output_data_with_extras.csv"
    data.to_csv(output_file, index=False)
    print(f"Modified data saved to {output_file}")

    # Plotting the last entry as an example
    last_idx = total_models - 1
    if not np.isnan(filter_magnitudes[last_idx]):
        print(f"Generating plot for the last entry (index {last_idx})...")
        (
            filter_mag,
            wavelength_angstrom,
            flux_lambda,
            filter_response,
            observed_flux_lambda,
        ) = compute_filter_magnitude(
            data.at[last_idx, "teff"],
            data.at[last_idx, "logg"],
            None,
            filter_wavelength_angstrom,
            filter_transmission,
        )
        plot_spectrum(
            wavelength_angstrom,
            flux_lambda,
            filter_response,
            observed_flux_lambda,
            filter_file,
        )
    else:
        print("No valid filter magnitude computed for the last entry. Skipping plot.")


if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process data for stellar magnitudes.")
    parser.add_argument(
        "-s",
        "--synth",
        action="store_true",
        help="Run the synthetic data processing instead of the main program.",
    )
    args = parser.parse_args()

    # Decide which function to run
    if args.synth:
        synth_main()
    else:
        main()
