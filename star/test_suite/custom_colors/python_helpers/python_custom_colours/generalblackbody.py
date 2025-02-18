import numpy as np
import os
from scipy import integrate, interpolate
from synphot import SourceSpectrum

# Constants (as before)
c = 299792458.0
teffsun = 5777.0
rsun = 6.95700 * 10**8
h = 6.6260693 * 10 ** (-34)
k = 1.380658 * 10 ** (-23)
pc = 3.0857 * 10**16
stefan = (2.0 * (np.pi**5) * (k**4)) / (15 * (c**2) * (h**3))
lsun = 4.0 * np.pi * stefan * rsun**2 * teffsun**4
magZero = 0.03
mbolsun = 4.77

# Vega calibration
vega_spectrum = SourceSpectrum.from_vega()
vega_wave = vega_spectrum.waveset.value * 10**-10
vega_flux = vega_spectrum(vega_spectrum.waveset).value * 10**7
vega_interp = interpolate.interp1d(
    vega_wave, vega_flux + 1e-30, bounds_error=False, fill_value=0.0
)


# Blackbody function for wavelength
def bb_wave(wave, teff):
    return ((2 * np.pi * h * c**2) / (wave**5)) * (
        1.0 / (np.exp((h * c) / (wave * k * teff)) - 1.0)
    )


# Constant calculation function
def const(t):
    return (
        mbolsun
        - 2.5 * np.log10((4.0 * np.pi * (10.0 * pc) ** 2 * stefan * t**4) / lsun)
        - magZero
    )


# Modified function for bolometric correction
def do_one(f_in, vega_interp, teff, log_g):
    print(f"Processing {f_in}")
    j = np.genfromtxt(f_in, names=["wave_nm", "flux"], delimiter=",", skip_header=1)
    j_m = j["wave_nm"] * 10**-10
    j_f = j["flux"]
    j_f_interp = interpolate.interp1d(
        j_m, j_f + 1e-30, bounds_error=False, fill_value=0.0
    )

    # Calculate normalization
    j_norm, _ = integrate.quad(j_f_interp, j_m[0], j_m[-1], epsabs=1e-5, epsrel=1e-5)

    # Modified bolometric correction function
    def bc(t, g):
        try:
            top, _ = integrate.quad(
                lambda wave: wave * bb_wave(wave, t) * (j_f_interp(wave) / j_norm),
                j_m[0],
                j_m[-1],
                epsabs=1e-5,
                epsrel=1e-5,
            )
            bot, _ = integrate.quad(
                lambda wave: wave * vega_interp(wave) * (j_f_interp(wave) / j_norm),
                j_m[0],
                j_m[-1],
                epsabs=1e-5,
                epsrel=1e-5,
            )
            if bot == 0:  # Avoid division by zero
                print(
                    f"Warning: Zero denominator encountered in {f_in} for T={t}, log_g={g}."
                )
                return np.nan
            # Example adjustment: scale const(t) with log_g
            return (
                const(t) - 0.1 * g + 2.5 * np.log10(top / bot)
            )  # Adjust as needed based on BC definition
        except Exception as e:
            print(f"Error in bc calculation for {f_in}, T={t}, log_g={g}: {e}")
            return np.nan

    return np.array([[bc(temp, g) for g in log_g] for temp in teff])


# Main processing function with log_g range
def process_filter_families(base_dir, teff, log_g):
    for root, dirs, files in os.walk(base_dir):
        family = os.path.relpath(root, base_dir).split(os.sep)[0]
        if not files:
            continue

        filter_data = []
        for f in files:
            filter_path = os.path.join(root, f)
            filter_data.append((f, do_one(filter_path, vega_interp, teff, log_g)))

        # Write output
        output_file = f"blackbody_{family}.dat"
        with open(output_file, "w") as f_out:
            header = "#Teff logg mdivh " + " ".join(
                [f"{filter_name.split('.')[0]}_bb" for filter_name, _ in filter_data]
            )
            print(header, file=f_out)
            for i, temp in enumerate(teff):
                for j, g in enumerate(log_g):
                    data_row = f"{temp} {g} 0.0 " + " ".join(
                        [str(data[i][j]) for _, data in filter_data]
                    )
                    print(data_row, file=f_out)


# Define temperature and log_g range
teff = np.arange(100.0, 50001.0, 100.0)
log_g = np.arange(0.0, 5.1, 0.1)  # Example log_g range
base_dir = "data/filters/"

# Run processing
process_filter_families(base_dir, teff, log_g)
