import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Physical constants
h = 6.62607015e-27   # Planck constant (erg·s)
c = 2.99792458e10    # Speed of light (cm/s)
k_B = 1.380649e-16   # Boltzmann constant (erg/K)

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

def load_stellar_spectrum(T_eff, log_g, metallicity, wavelength_min, wavelength_max):
    """
    Placeholder function to load stellar atmosphere spectrum.
    Currently uses a blackbody approximation.
    """
    # Adjust the wavelength range to cover the filter's range (~5.6 microns)
    wavelength_angstrom = np.linspace(wavelength_min, wavelength_max, 10000)  # 3 to 8 microns in Å
    wavelength_cm = wavelength_angstrom * 1e-8  # Convert Å to cm

    # Calculate blackbody spectral radiance (erg/s/cm^2/Å)
    exponent = (h * c) / (wavelength_cm * k_B * T_eff)
    exponent = np.minimum(exponent, 700)  # Avoid overflow
    flux_lambda = (2 * h * c**2) / (wavelength_cm**5) / (np.exp(exponent) - 1)
    flux_lambda = flux_lambda * 1e-8  # Convert from per cm to per Å

    return wavelength_angstrom, flux_lambda

# Read MESA history data
history = read_mesa_history('LOGS/history.data')

# Display available columns
print("Available columns:", history.columns.tolist())

# Extract stellar parameters
T_eff = history['Teff'].iloc[-1]
log_g = history['log_g'].iloc[-1]
metallicity = history.get('initial_feh', 0.0)

print(f"Effective Temperature: {T_eff} K")
print(f"log(g): {log_g}")
print(f"Metallicity [Fe/H]: {metallicity}")


# Load filter transmission data
filterloc = 'filters/JWST/MIRI/F560W.dat'
filter_wavelength_micron, filter_transmission = np.loadtxt(
    filterloc, skiprows=1, delimiter=',', unpack=True
)
print(filter_wavelength_micron)

# **Correct the wavelength units**
# Convert filter wavelengths from microns to Angstroms
filter_wavelength_angstrom = filter_wavelength_micron * 1e4  # microns to Å

# Load stellar atmosphere spectrum
wavelength_min = np.min(filter_wavelength_angstrom)
wavelength_max = np.max(filter_wavelength_angstrom)
wavelength_angstrom, flux_lambda = load_stellar_spectrum(T_eff, log_g, metallicity, wavelength_min, wavelength_max)
# Interpolate filter response onto the stellar spectrum wavelength grid
filter_interp = interp1d(
    filter_wavelength_angstrom, filter_transmission, bounds_error=False, fill_value=0
)
filter_response = filter_interp(wavelength_angstrom)

# Convolve spectrum with filter response
observed_flux_lambda = flux_lambda * filter_response

# Integrate over wavelength to get flux (erg/s/cm^2)
flux = np.trapz(observed_flux_lambda, wavelength_angstrom)

# **Calculate zero-point flux for the filter**
# Using the AB magnitude system zero-point flux
# Zero-point flux density for AB system: F_nu = 3631 Jy
F_zero_nu = 3631e-23  # Convert Jy to erg/s/cm^2/Hz

# Convert F_nu to F_lambda at the effective wavelength
lambda_eff = 5.6e4  # 5.6 microns in Å
F_zero_lambda = (F_zero_nu * c) / (lambda_eff**2)  # erg/s/cm^2/Å

# Calculate the absolute magnitude in the AB system
M_filter = -2.5 * np.log10(flux / F_zero_lambda)
print(f"Absolute Magnitude ({filterloc}): {M_filter}")

# Plotting the spectrum, filter response, and convolved spectrum
plt.figure(figsize=(10, 6))
plt.plot(wavelength_angstrom, flux_lambda, label='Stellar Spectrum')
plt.plot(
    filter_wavelength_angstrom,
    filter_transmission * np.max(flux_lambda),
    label='Filter Transmission'
)
plt.plot(wavelength_angstrom, observed_flux_lambda, label='Convolved Spectrum')
plt.xlabel('Wavelength (Å)')
plt.ylabel('Flux (erg/s/cm²/Å)')
plt.legend()
plt.show()


