import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Physical constants
h = 6.62607015e-34  # Planck constant (J·s)
c = 2.99792458e8    # Speed of light (m/s)
k_B = 1.380649e-23  # Boltzmann constant (J/K)

def read_mesa_history(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Find the index where the column headers start
    for idx, line in enumerate(lines):
        if line.strip().startswith('model_number'):
            header_line_idx = idx
            break
    else:
        raise ValueError('Header line with column names not found.')

    # Extract the header line
    headers = lines[header_line_idx].strip().split()

    # Data starts after the header line
    data_start_idx = header_line_idx + 1

    # Read the data into a DataFrame
    data = pd.read_csv(
        filename,
        skiprows=data_start_idx,
        delim_whitespace=True,
        names=headers,
        comment='#',
        engine='python'
    )

    return data

def planck(wavelength_m, T):
    exponent = (h * c) / (wavelength_m * k_B * T)
    exponent = np.minimum(exponent, 700)  # Avoid overflow in exponential
    numerator = 2 * h * c**2
    denominator = (wavelength_m**5) * (np.exp(exponent) - 1)
    return numerator / denominator

# Read MESA history data
history = read_mesa_history('LOGS/history.data')

# Display available columns
print("Available columns:", history.columns.tolist())

# Extract the effective temperature
possible_Teff_columns = ['Teff', 'photosphere_Teff', 'log_Teff']

for col in possible_Teff_columns:
    if col in history.columns:
        if col.startswith('log_'):
            T_eff = 10 ** history[col].iloc[-1]
        else:
            T_eff = history[col].iloc[-1]
        print(f"Effective Temperature ({col}): {T_eff} K")
        break
else:
    raise KeyError('Effective temperature column not found in history data.')


# Load filter transmission data
filter_data = np.loadtxt('data/filters/Liverpool_IOI.J.dat')

filter_wavelength_nm = filter_data[:, 0]
filter_transmission = filter_data[:, 1]

# Define the wavelength range
wavelength_nm = np.linspace(min(filter_wavelength_nm), max(filter_wavelength_nm), len(filter_wavelength_nm))  # in nanometers
wavelength_m = wavelength_nm * 1e-9  # Convert to meters

# Calculate blackbody spectral radiance
B_lambda = planck(wavelength_m, T_eff)

# Interpolate filter response
filter_interp = interp1d(filter_wavelength_nm, filter_transmission, bounds_error=False, fill_value=0)
filter_response = filter_interp(wavelength_nm)

# Convolve spectrum with filter
flux_lambda = B_lambda * filter_response

# Integrate over wavelength to get flux
flux = np.trapz(flux_lambda, wavelength_m)

# Zero-point flux for Vega in W/m^2/m
F_zero = 3.58e-8 / 1e-6  # Convert μm to m

# Calculate the absolute magnitude in the Vega system
M_V = -2.5 * np.log10(flux / F_zero)
print(f"Absolute Magnitude (Vega system): {M_V}")

# Plotting the spectrum and filter response (optional)
plt.figure(figsize=(10, 6))
plt.plot(wavelength_nm, B_lambda, label='Blackbody Spectrum')
plt.plot(filter_wavelength_nm, filter_transmission * max(B_lambda), label='Filter Transmission')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Spectral Radiance')
plt.legend()
plt.show()


