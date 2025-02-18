from scipy.interpolate import griddata
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import os
import sys


# Physical constants (cgs units)
h = 6.62607015e-27   # Planck constant (erg·s)
c = 2.99792458e10    # Speed of light (cm/s)
k_B = 1.380649e-16   # Boltzmann constant (erg/K)
R_sun = 6.957e10      # Solar radius in cm
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
        with open(filename, 'r') as f:
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
        if line.strip().startswith('model_number'):
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
            sep='\s+',  # Use regex to handle any whitespace
            names=headers,
            comment='#',
            engine='python'
        )
    except Exception as e:
        print(f"Error parsing history file: {e}")
        sys.exit(1)

    return data

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
        data = pd.read_csv(color_file, delim_whitespace=True, comment='#', header=None)
        # Rename columns based on your example
        data.columns = [
            'Teff', 'logg', 'M_div_H', 'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'L', 'Lprime', 'M'
        ]
    except Exception as e:
        print(f"Error reading color file '{color_file}': {e}")
        sys.exit(1)
    
    return data

def compute_bolometric_magnitude_from_color_file(T_eff, log_g, color_data, band='V'):
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
    points = color_data[['Teff', 'logg']].values
    values = color_data[band].values
    
    # Interpolate using griddata
    bolometric_correction = griddata(points, values, (T_eff, log_g), method='linear')
    
    if np.isnan(bolometric_correction):
        print(f"Warning: Bolometric correction not found for T_eff={T_eff}, log_g={log_g}.")
        return np.nan
    
    return bolometric_correction

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
        filter_data = np.loadtxt(filter_file, comments='#', delimiter = ',', skiprows = 1)
    except Exception as e:
        print(f"Error loading filter file '{filter_file}': {e}")
        sys.exit(1)
    
    if filter_data.ndim != 2 or filter_data.shape[1] < 2:
        print(f"Error: Filter file '{filter_file}' has an unexpected format.")
        sys.exit(1)
    
    filter_wavelength_units = filter_data[:, 0]
    filter_transmission = filter_data[:, 1]
    
    # Assuming filter wavelengths are in Angstroms. If not, adjust accordingly.
    # Uncomment the appropriate line based on your filter file's wavelength units.
    
    # If wavelengths are in microns (µm):
    # filter_wavelength_angstrom = filter_wavelength_units * 1e4  # µm to Å
    
    # If wavelengths are in nanometers (nm):
    # filter_wavelength_angstrom = filter_wavelength_units * 10    # nm to Å
    
    # If wavelengths are already in Angstroms (Å):
    filter_wavelength_angstrom = filter_wavelength_units

    return filter_wavelength_angstrom, filter_transmission



def planck(wavelength_angstrom, T_eff):
    # Convert wavelength from Å to cm
    wavelength_cm = wavelength_angstrom * 1e-8  # Å to cm

    # Blackbody flux formula in erg/s/cm²/cm
    exponent = (h * c) / (wavelength_cm * k_B * T_eff)
    flux_per_cm = (2.0 * h * c**2) / (wavelength_cm**5) / (np.exp(exponent) - 1.0)

    # Convert flux from per cm to per Å
    flux_per_angstrom = flux_per_cm * 1e-8  # 1 cm = 1e8 Å

    return flux_per_angstrom
import numpy as np
from scipy.interpolate import interp1d

# Physical constants
c = 2.99792458e10        # Speed of light in cm/s
R_sun = 6.957e10         # Solar radius in cm
pc_in_cm = 3.0857e18     # Parsec in cm
h = 6.62607015e-27       # Planck constant in erg*s

def planck(wavelength_angstrom, T_eff):
    """
    Computes the Planck function B_lambda for given wavelengths and temperature.
    
    Parameters:
        wavelength_angstrom (np.ndarray): Wavelengths in Angstroms.
        T_eff (float): Effective temperature in Kelvin.
    
    Returns:
        np.ndarray: Spectral radiance in erg/s/cm²/Å.
    """
    wavelength_cm = wavelength_angstrom * 1e-8  # Convert Å to cm
    exponent = (h * c) / (wavelength_cm * 1.380649e-16 * T_eff)  # h*c/(k*T)
    B_lambda = (2.0 * h * c**2) / (wavelength_cm**5 * (np.exp(exponent) - 1.0))
    return B_lambda

def compute_filter_magnitude(T_eff, log_g, log_R, filter_wavelength_angstrom, filter_transmission):
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
    print('teff', T_eff)
    print('logg', log_g)
    print('logr', log_R)
    wavelength_min = np.min(filter_wavelength_angstrom)
    wavelength_max = np.max(filter_wavelength_angstrom)
    wavelength_length = int(len(filter_wavelength_angstrom) * 1.5)

    # Generate wavelength array around filter range
    wavelength_angstrom = np.linspace(wavelength_min * 0.9, wavelength_max * 1.1, wavelength_length)
    print("Generated Wavelengths (middle 5):", wavelength_angstrom[int(len(wavelength_angstrom)/2):int(len(wavelength_angstrom)/2)+5])
    print("Filter Wavelengths (middle 5):", filter_wavelength_angstrom[int(len(filter_wavelength_angstrom)/2):int(len(filter_wavelength_angstrom)/2)+5])

    # Calculate flux_lambda using Planck function
    flux_lambda = planck(wavelength_angstrom, T_eff)  # erg/s/cm²/Å
    print("Flux Lambda (middle 5):", flux_lambda[int(len(flux_lambda)/2):int(len(flux_lambda)/2)+5])

    # Convert log_R to R in cm
    R = 10**log_R * R_sun  # cm
    print("Stellar Radius (cm):", R)

    # Define distance to star (10 pc in cm)
    d = 10 * pc_in_cm  # cm
    print("Distance to Star (cm):", d)

    # Apply flux scaling: F_lambda = pi * (R/d)^2 * B_lambda
    scaling_factor = np.pi * (R / d)**2
    flux_lambda *= scaling_factor
    print("Scaling Factor (pi*(R/d)^2):", scaling_factor)
    print("Scaled Flux Lambda (middle 5):", flux_lambda[int(len(flux_lambda)/2):int(len(flux_lambda)/2)+5])

    # Clip stellar spectrum to filter wavelength range
    mask = (wavelength_angstrom >= wavelength_min) & (wavelength_angstrom <= wavelength_max)
    wavelength_angstrom_clipped = wavelength_angstrom[mask]
    flux_lambda_clipped = flux_lambda[mask]
    print("Clipped Flux Lambda (middle 5):", flux_lambda_clipped[int(len(flux_lambda_clipped)/2):int(len(flux_lambda_clipped)/2)+5])

    # Interpolate filter response onto the stellar spectrum wavelength grid
    filter_interp = interp1d(
        filter_wavelength_angstrom, filter_transmission,
        bounds_error=False, fill_value=0.0
    )
    filter_response_clipped = filter_interp(wavelength_angstrom_clipped)
    print("Filter Response Clipped (middle 5):", filter_response_clipped[int(len(filter_response_clipped)/2):int(len(filter_response_clipped)/2)+5])

    # Convolve spectrum with filter response
    observed_flux_lambda = flux_lambda_clipped * filter_response_clipped  # erg/s/cm²/Å
    print("Observed Flux Lambda (middle 5):", observed_flux_lambda[int(len(observed_flux_lambda)/2):int(len(observed_flux_lambda)/2)+5])

    # Integrate F_lambda * lambda over wavelength to get F_nu * c
    integrand = observed_flux_lambda * wavelength_angstrom_clipped * 1e-8  # Convert Å to cm
    flux_nu_c = np.trapz(integrand, wavelength_angstrom_clipped * 1e-8)  # erg/s/cm²
    flux_nu = flux_nu_c / c  # erg/s/cm²/Hz
    print("Flux_nu (erg/s/cm²/Hz):", flux_nu)

    if flux_nu <= 0:
        print(f"Warning: Non-positive F_nu encountered (T_eff={T_eff}, log_g={log_g}).")
        return np.nan

    # AB magnitude zero-point flux
    F_zero_nu = 3631e-23  # erg/s/cm²/Hz

    # Calculate the absolute AB magnitude
    M_filter = -2.5 * np.log10(flux_nu / F_zero_nu)
    print("AB Magnitude:", M_filter)

    return M_filter, wavelength_angstrom_clipped, flux_lambda_clipped, filter_response_clipped, observed_flux_lambda




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
    history['bol_magnitude'] = bol_mag
    history['filter_magnitude'] = filter_mag
    history['filter_used'] = filter_used
    return history


def write_modified_history(history, output_file):
    """
    Writes the modified history DataFrame to a file, ensuring no duplicate headers
    and formatting floats without scientific notation.
    
    Parameters:
        history (pd.DataFrame): Modified history DataFrame.
        output_file (str): Path to the output file.
    """
    #print(history)
    try:
        with open(output_file, 'w') as f:
            # Write header only once
            f.write(','.join(history.columns) + '\n')
            # Write data without scientific notation
            history.to_csv(f, sep=',', index=False, header=False, float_format='%.6f')
        print(f"Updated history data saved to '{output_file}'.")
    except Exception as e:
        print(f"Error writing to output file '{output_file}': {e}")



def plot_spectrum(wavelength_angstrom, flux_lambda, filter_response, observed_flux_lambda, filter_file):
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
    plt.plot(wavelength_angstrom, flux_lambda, label='Stellar Spectrum')
    plt.plot(
        wavelength_angstrom,
        filter_response * np.max(flux_lambda),
        label='Filter Transmission (scaled)'
    )
    plt.plot(wavelength_angstrom, observed_flux_lambda, label='Convolved Spectrum')
    plt.xlabel('Wavelength (Å)')
    plt.ylabel('Flux (erg/s/cm²/Å)')
    plt.legend()
    plt.title(f'Spectrum Convolution with Filter: {os.path.basename(filter_file)}')
    plt.grid(True)
    plt.show()

def main():
    """
    Main function to orchestrate the computation of bolometric and filter magnitudes,
    and to generate a modified history file with the new data.
    """
    default_filter = 'data/filters/PAN_STARRS/PAN_STARRS/g.dat'
    # User specifies the filter file
    filter_file = ''#input("Enter the path to the filter file (e.g., default_filter): ").strip()
    if filter_file == '':
        filter_file = default_filter
    
    if not os.path.isfile(filter_file):
        print(f"Error: Filter file '{filter_file}' does not exist. Please check the path.")
        sys.exit(1)
    
    # Read filter data
    filter_wavelength_angstrom, filter_transmission = read_filter_data(filter_file)
    
    # Read MESA history data
    history_file = 'LOGS/history.data'
    if not os.path.isfile(history_file):
        print(f"Error: History file '{history_file}' does not exist. Please check the path.")
        sys.exit(1)
    
    history = read_mesa_history(history_file)

    # Initialize lists to store computed magnitudes
    bol_magnitudes = []
    filter_magnitudes = []
    
    # Process each model in the history
    total_models = len(history)
    print(f"Processing {total_models} models...")
    
    for idx, row in history.iterrows():
        T_eff = row.get('Teff')
        log_g = row.get('log_g')
        log_R = row.get('log_R')        
        metallicity = row.get('initial_feh', 0.0)
        
        # Check for missing values
        if pd.isnull(T_eff) or pd.isnull(log_g):
            print(f"Warning: Missing T_eff or log_g for model index {idx}. Assigning NaN.")
            bol_magnitudes.append(np.nan)
            filter_magnitudes.append(np.nan)
            continue
        
        # Compute bolometric magnitude
        bol_mag = compute_bolometric_magnitude(T_eff, log_g)
        bol_magnitudes.append(bol_mag)
        
        # Compute filter magnitude
        filter_mag, wavelength_angstrom, flux_lambda, filter_response, observed_flux_lambda = compute_filter_magnitude(T_eff, log_g, log_R, filter_wavelength_angstrom, filter_transmission)
        
        filter_magnitudes.append(filter_mag)
        exit()
        # Optional: Print progress every 10%
        #if (idx + 1) % max(1, total_models // 10) == 0:
            #print(f"Processed {idx + 1}/{total_models} models.")
    
    # Add new columns to history
    history = generate_modified_history(
        history,
        pd.Series(bol_magnitudes),
        pd.Series(filter_magnitudes),
        os.path.basename(filter_file),
        filter_file
    )
    
    # Write the modified history to a new file
    output_file = 'history_with_extras.data'
    write_modified_history(history, output_file)
    
    # Plotting the last model as an example
    last_model_idx = total_models - 1
    if not np.isnan(filter_magnitudes[last_model_idx]):
        print(f"Generating plot for the last model (index {last_model_idx})...")
        filter_mag, wavelength_angstrom, flux_lambda, filter_response, observed_flux_lambda = compute_filter_magnitude(
            history.at[last_model_idx, 'Teff'],
            history.at[last_model_idx, 'log_g'],
            history.at[last_model_idx, 'log_R'],
            filter_wavelength_angstrom,
            filter_transmission
        )
        plot_spectrum(wavelength_angstrom, flux_lambda, filter_response, observed_flux_lambda, filter_file)
    else:
        print("No valid filter magnitude computed for the last model. Skipping plot.")

if __name__ == '__main__':
    main()

