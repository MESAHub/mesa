import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import os
import sys
import matplotlib.pyplot as plt

# Physical constants (cgs units)
h = 6.62607015e-27  # Planck constant (erg·s)
c = 2.99792458e10  # Speed of light (cm/s)
k_B = 1.380649e-16  # Boltzmann constant (erg/K)


def read_mesa_history(filename):
    """
    Reads the MESA history file, neatly prints the content,
    and parses it into a pandas DataFrame.
    """
    print(f"Reading MESA history file: {filename}")
    try:
        with open(filename, "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: History file '{filename}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading history file: {e}")
        sys.exit(1)

    # Print the contents of the history file
    print("\n--- History File Contents (Preview) ---")
    for idx, line in enumerate(lines[:20]):  # Print the first 20 lines for brevity
        # Truncate long lines and print with line numbers
        line_content = line.strip()
        truncated_line = (
            (line_content[:120] + "...") if len(line_content) > 120 else line_content
        )
        print(f"{idx + 1:>4}: {truncated_line}")
    print("--- End of History File Preview ---\n")

    # Find the header line
    header_line_idx = None
    for idx, line in enumerate(lines):
        if line.strip().startswith("model_number"):
            header_line_idx = idx
            break
    if header_line_idx is None:
        print("Error: Header line with column names not found in history file.")
        sys.exit(1)

    print(f"Header found at line {header_line_idx + 1}:")
    print(lines[header_line_idx].strip())  # Print header without truncation

    headers = lines[header_line_idx].strip().split()

    # Parse the file into a DataFrame
    try:
        data = pd.read_csv(
            filename,
            skiprows=header_line_idx + 1,
            sep=r"\s+",
            names=headers,
            comment="#",
            engine="python",
            nrows=1,  # Load only N lines
        )
    except Exception as e:
        print(f"Error parsing history file: {e}")
        sys.exit(1)

    print(f"History file successfully parsed with {len(data)} entries.")
    return data


def get_closest_stellar_models(teff, log_g, lookup_table):
    print(f"Finding closest stellar models for Teff={teff}, log_g={log_g}")
    distances = np.sqrt(
        (lookup_table[" teff"] - teff) ** 2 + (lookup_table[" logg"] - log_g) ** 2
    )
    closest_indices = np.argsort(distances)[:2]
    if len(closest_indices) < 2:
        raise ValueError("Not enough models for interpolation.")
    print(f"Closest models found: {closest_indices}")
    return lookup_table.iloc[closest_indices]


def interpolate_sed(teff, log_g, stellar_model_dir, lookup_table):
    print(f"Interpolating SED for Teff={teff}, log_g={log_g}")
    closest_models = get_closest_stellar_models(teff, log_g, lookup_table)

    print(f"Closest models for interpolation:\n {closest_models}")

    # Check for a perfect match
    if (closest_models.iloc[0][" teff"] == teff) and (
        closest_models.iloc[0][" logg"] == log_g
    ):
        print("Perfect match found for the stellar model.")
        sed_path = os.path.join(stellar_model_dir, closest_models.iloc[0]["#file_name"])
        sed = pd.read_csv(
            sed_path, sep="\s+", comment="#", header=None, names=["wavelength", "flux"]
        )
        return sed["wavelength"], sed["flux"]

    sed1_path = os.path.join(stellar_model_dir, closest_models.iloc[0]["#file_name"])
    sed2_path = os.path.join(stellar_model_dir, closest_models.iloc[1]["#file_name"])
    print(f"Loading SEDs: {sed1_path}, {sed2_path}")

    sed1 = pd.read_csv(
        sed1_path, sep="\s+", comment="#", header=None, names=["wavelength", "flux"]
    )
    sed2 = pd.read_csv(
        sed2_path, sep="\s+", comment="#", header=None, names=["wavelength", "flux"]
    )

    print(f"SED1 head:\n{sed1.head()}")
    print(f"SED2 head:\n{sed2.head()}")

    # Check if the two closest models are identical in Teff
    if closest_models.iloc[0][" teff"] == closest_models.iloc[1][" teff"]:
        print("Closest models have identical Teff; skipping interpolation.")
        return sed1["wavelength"], sed1["flux"]  # Return one model directly

    # Interpolation weights
    weight1 = np.abs(closest_models.iloc[1][" teff"] - teff) / (
        np.abs(closest_models.iloc[1][" teff"] - closest_models.iloc[0][" teff"])
    )
    weight2 = 1.0 - weight1
    print(f"Interpolation weights: weight1={weight1}, weight2={weight2}")

    # Interpolate the flux
    flux = weight1 * sed1["flux"] + weight2 * sed2["flux"]
    return sed1["wavelength"], flux


def calculate_bolometric_magnitude(wavelength, flux):
    print("Calculating bolometric magnitude...")
    bolometric_flux = np.trapz(flux, wavelength)
    bolometric_magnitude = -2.5 * np.log10(bolometric_flux) - 48.6
    print(f"Bolometric magnitude calculated: {bolometric_magnitude}")
    return bolometric_magnitude


def convolve_with_filter(wavelength, flux, filter_dir, vega_sed):
    print(f"Convolving SED with filters in: {filter_dir}")
    filter_files = [f for f in os.listdir(filter_dir) if f.endswith(".dat")]
    print(f"Found {len(filter_files)} filter files: {filter_files}")

    # Convert wavelength to cm and flux to erg/s/cm²/Hz
    wavelength_cm = wavelength * 1e-8
    flux_nu = flux * (wavelength_cm**2) / c
    results = []
    summary = []

    vega_interp = interp1d(
        vega_sed["wavelength"], vega_sed["flux"], bounds_error=False, fill_value=0
    )

    for filter_file in filter_files:
        filter_path = os.path.join(filter_dir, filter_file)
        filter_data = pd.read_csv(
            filter_path, sep=",", header=0, names=["wavelength", "transmission"]
        )

        # Interpolate Vega SED onto the filter wavelength grid
        vega_flux_at_filter = vega_interp(filter_data["wavelength"])
        # Calculate λ_eff using the Vega spectrum
        numerator = np.trapz(
            filter_data["wavelength"]
            * filter_data["transmission"]
            * vega_flux_at_filter,
            filter_data["wavelength"],
        )
        denominator = np.trapz(
            filter_data["transmission"] * vega_flux_at_filter, filter_data["wavelength"]
        )
        lambda_eff = numerator / denominator

        peak_wavelength = filter_data["wavelength"][
            filter_data["transmission"].idxmax()
        ]
        filter_interp = interp1d(
            filter_data["wavelength"],
            filter_data["transmission"],
            bounds_error=False,
            fill_value=0,
        )
        filter_response = filter_interp(wavelength)

        # Restrict integration to valid filter response range
        valid = filter_response > 0
        numerator = np.trapz(flux_nu[valid] * filter_response[valid], wavelength[valid])
        denominator = np.trapz(filter_response[valid], wavelength[valid])

        if denominator == 0:
            m_AB = np.nan
            f_lambda = np.nan
            warning = f"Warning: Filter {filter_file} has no valid transmission in the SED range."
        else:
            f_nu = numerator / denominator
            m_AB = -2.5 * np.log10(f_nu) - 48.6

            # Calculate f_lambda at the peak wavelength
            f_lambda = (
                f_nu * c / (peak_wavelength * 1e-8) ** 2
            )  # Convert peak_wavelength to cm
            warning = None

        results.append((m_AB, peak_wavelength))
        summary.append(
            {
                "Filter": filter_file,
                "λ_eff": lambda_eff,
                "Peak Wavelength": peak_wavelength,
                "Magnitude (m_AB)": m_AB,
                "Flux (f_nu)": f_nu if denominator != 0 else None,
                "Flux (f_lambda)": f_lambda if denominator != 0 else None,
                "Warning": warning,
            }
        )

        def format_value(val):
            """Formats values with up to 10 decimal places; switches to scientific if below 1e-10."""
            if val is None or np.isnan(val):
                return "N/A"
            elif abs(val) < 1e-10:
                return f"{val:.3e}"  # Use scientific notation
            else:
                return f"{val:.10f}".rstrip("0").rstrip(".")  # Remove trailing zeros

    print("\nSummary of Results:")
    print("-" * 50)
    for s in summary:
        print(f"Filter: {s['Filter']}")
        print(f"  Effective Wavelength: {format_value(s['λ_eff'])}")
        print(f"  Peak Wavelength: {format_value(s['Peak Wavelength'])}")
        print(f"  Magnitude (m_AB): {format_value(s['Magnitude (m_AB)'])}")
        print(f"  Flux (f_nu): {format_value(s['Flux (f_nu)'])}")
        print(f"  Flux (f_lambda): {format_value(s['Flux (f_lambda)'])}")
        if s["Warning"]:
            print(f"  {s['Warning']}")
        print("-" * 50)

    return results


def plot_magnitudes(all_filter_results, output_file):
    """
    Plots magnitudes against their corresponding peak wavelengths.

    Parameters:
        all_filter_results (list): List of lists containing (magnitude, peak_wavelength) tuples.
        output_file (str): Path to the output file used to save the history DataFrame.
    """
    print("Plotting magnitudes...")
    plt.figure(figsize=(10, 10))

    for model_idx, filter_results in enumerate(all_filter_results):
        magnitudes, wavelengths = zip(*filter_results)
        plt.scatter(wavelengths, magnitudes, label=f"Model {model_idx + 1}", alpha=0.6)

    plt.xlabel("Wavelength (Angstroms)")
    plt.ylabel("AB Magnitude")
    plt.gca().invert_yaxis()  # Magnitudes are inverted (smaller is brighter)
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize="small")
    plt.tight_layout()
    plot_file = output_file.replace(".csv", "_magnitudes_plot.png")
    plt.savefig(plot_file)
    print(f"Plot saved to {plot_file}")
    plt.show()


def main():
    history_file = "LOGS/history.data"
    stellar_model = "data/stellar_models/Kurucz/"
    instrument = "data/filters/JWST/NIRCam/"
    vega_sed_file = "data/stellar_models/vega_sed.txt"
    vega_sed = pd.read_csv(
        vega_sed_file, sep=",", names=["wavelength", "flux"], comment="#"
    )  # https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/castelli-and-kurucz-atlas
    print(f"Starting main process...")
    lookup_table_file = os.path.join(stellar_model, "lookup_table.csv")
    lookup_table = pd.read_csv(lookup_table_file)
    print(f"Lookup table loaded: {lookup_table_file}")

    bolometric_magnitudes = []
    all_filter_results = []

    history = read_mesa_history(history_file)
    print(f"Processing {len(history)} models from history file...")

    for idx, row in history.iterrows():
        print(f"Processing model {idx + 1}/{len(history)}")
        teff = 3500  # row.get('Teff')
        log_g = 3  # row.get('log_g')
        log_R = row.get("log_R")
        metallicity = -0.5  # row.get('initial_feh', 0.0)

        wavelength, flux = interpolate_sed(teff, log_g, stellar_model, lookup_table)
        bol_mag = calculate_bolometric_magnitude(wavelength, flux)
        bolometric_magnitudes.append(bol_mag)

        filter_results = convolve_with_filter(wavelength, flux, instrument, vega_sed)
        all_filter_results.append(
            filter_results
        )  # Store tuples of (mag, peak_wavelength)

    # Add bolometric magnitudes to history
    history["bolometric_magnitude"] = bolometric_magnitudes

    # Add filter magnitudes as separate columns in the DataFrame
    for idx, filter_results in enumerate(all_filter_results):
        for mag, peak_wavelength in filter_results:
            filter_name = f"mag_{int(peak_wavelength)}"
            if filter_name not in history:
                history[filter_name] = np.nan
            history.at[idx, filter_name] = mag

    # Save results
    output_file = "history_with_magnitudes.csv"
    history.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}.")

    # Plot magnitudes against wavelengths
    plot_magnitudes(all_filter_results, output_file)


if __name__ == "__main__":
    main()
