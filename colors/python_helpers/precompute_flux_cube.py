#!/usr/bin/env python3
import argparse
import os
import struct

import numpy as np
from tqdm import tqdm


def load_sed(filepath, index):
    """Load a spectral energy distribution file."""
    wavelengths = []
    fluxes = []

    with open(filepath, "r") as f:
        # Skip header lines (comments)
        for line in f:
            if line.strip() and not line.startswith("#"):
                break

        # Read wavelength and flux values
        try:
            values = line.strip().split()
            wavelengths.append(float(values[0]))
            fluxes.append(float(values[1]))
        except (ValueError, IndexError):
            pass  # Skip any malformed lines

        for line in f:
            try:
                values = line.strip().split()
                wavelengths.append(float(values[0]))
                fluxes.append(float(values[1]))
            except (ValueError, IndexError):
                pass  # Skip any malformed lines

    return np.array(wavelengths), np.array(fluxes)


def load_lookup_table(lookup_file):
    """Load the lookup table with model parameters."""
    file_names = []
    teff_values = []
    logg_values = []
    meta_values = []

    with open(lookup_file, "r") as f:
        # Skip header line
        header = f.readline().strip().split(",")

        # Find column indices
        file_col = 0  # Assume first column is filename
        teff_col = None
        logg_col = None
        meta_col = None

        for i, col in enumerate(header):
            col = col.strip().lower()
            if col == "teff":
                teff_col = i
            elif col == "logg":
                logg_col = i
            elif col in ["meta", "feh"]:
                meta_col = i

        # Read data rows
        for line in f:
            if line.strip():
                values = line.strip().split(",")
                if len(values) <= max(
                    file_col, teff_col or 0, logg_col or 0, meta_col or 0
                ):
                    continue  # Skip lines that don't have enough values

                file_names.append(values[file_col].strip())

                try:
                    teff = float(values[teff_col]) if teff_col is not None else 0.0
                    logg = float(values[logg_col]) if logg_col is not None else 0.0
                    meta = float(values[meta_col]) if meta_col is not None else 0.0

                    teff_values.append(teff)
                    logg_values.append(logg)
                    meta_values.append(meta)
                except (ValueError, IndexError):
                    # Skip rows with invalid values
                    file_names.pop()  # Remove the added filename

    print(
        f"Column indices found - teff: {teff_col}, logg: {logg_col}, meta: {meta_col}"
    )
    return (
        file_names,
        np.array(teff_values),
        np.array(logg_values),
        np.array(meta_values),
    )


def get_unique_sorted(values, tolerance=1e-8):
    """Get sorted unique values with tolerance."""
    sorted_vals = np.sort(values)
    unique_indices = [0]

    for i in range(1, len(sorted_vals)):
        if abs(sorted_vals[i] - sorted_vals[unique_indices[-1]]) > tolerance:
            unique_indices.append(i)

    return sorted_vals[unique_indices]


def save_binary_file(
    output_file, wavelengths, teff_grid, logg_grid, meta_grid, flux_cube
):
    with open(output_file, "wb") as f:
        # Write dimensions
        f.write(
            struct.pack(
                "4i", len(teff_grid), len(logg_grid), len(meta_grid), len(wavelengths)
            )
        )

        # Write grid arrays
        teff_grid.astype(np.float64).tofile(f)
        logg_grid.astype(np.float64).tofile(f)
        meta_grid.astype(np.float64).tofile(f)
        wavelengths.astype(np.float64).tofile(f)

        # FIXED: Transpose to match Fortran's column-major order expectations
        # This swaps the dimension order to (wavelength, meta, logg, teff)
        transposed_cube = np.transpose(flux_cube, (3, 2, 1, 0))
        transposed_cube.astype(np.float64).tofile(f)

        # Write flux cube - make sure it's in Fortran column-major order
        flux_cube.astype(np.float64).tofile(f)

        # Also save a text version of the grids for debugging
        debug_file = output_file[:-4] + ".txt"
        with open(debug_file, "w") as df:
            df.write("# Teff grid\n")
            for val in teff_grid:
                df.write(f"{val}\n")
            df.write("\n# Logg grid\n")
            for val in logg_grid:
                df.write(f"{val}\n")
            df.write("\n# Metallicity grid\n")
            for val in meta_grid:
                df.write(f"{val}\n")


def precompute_flux_cube(stellar_model_dir, output_file):
    """Precompute the 3D flux cube for all wavelengths."""
    # Load the lookup table
    lookup_file = os.path.join(stellar_model_dir, "lookup_table.csv")
    file_names, teff_values, logg_values, meta_values = load_lookup_table(lookup_file)

    print(f"Loaded {len(file_names)} model files from lookup table")

    # Get unique parameter values
    teff_grid = get_unique_sorted(teff_values)
    logg_grid = get_unique_sorted(logg_values)
    meta_grid = get_unique_sorted(meta_values)

    print(f"Unique Teff values: {len(teff_grid)}")
    for i, v in enumerate(teff_grid):
        if i < 5 or i > len(teff_grid) - 5:
            print(f"  Teff[{i}] = {v}")

    print(f"Unique logg values: {len(logg_grid)}")
    for i, v in enumerate(logg_grid):
        if i < 5 or i > len(logg_grid) - 5:
            print(f"  logg[{i}] = {v}")

    print(f"Unique metallicity values: {len(meta_grid)}")
    for i, v in enumerate(meta_grid):
        if i < 5 or i > len(meta_grid) - 5:
            print(f"  meta[{i}] = {v}")

    # Load first SED to get wavelength grid
    first_file = os.path.join(stellar_model_dir, file_names[0])
    wavelengths, _ = load_sed(first_file, 0)
    n_lambda = len(wavelengths)

    print(f"Wavelength grid has {n_lambda} points")

    # Create the flux cube
    flux_cube = np.zeros((len(teff_grid), len(logg_grid), len(meta_grid), n_lambda))

    # Create a map to track where we've filled in the cube
    filled_map = np.zeros((len(teff_grid), len(logg_grid), len(meta_grid)), dtype=bool)

    # Process all model files
    model_index = 0
    missing_points = 0

    for file_name, teff, logg, meta in tqdm(
        zip(file_names, teff_values, logg_values, meta_values),
        total=len(file_names),
        desc="Processing models",
    ):
        # Find the indices in the grids
        i_teff = np.searchsorted(teff_grid, teff)
        if i_teff == len(teff_grid) or teff_grid[i_teff] != teff:
            if i_teff > 0 and i_teff < len(teff_grid):
                # Check which one is closer
                if abs(teff_grid[i_teff - 1] - teff) < abs(teff_grid[i_teff] - teff):
                    i_teff = i_teff - 1
            elif i_teff == len(teff_grid):
                i_teff = i_teff - 1

        i_logg = np.searchsorted(logg_grid, logg)
        if i_logg == len(logg_grid) or logg_grid[i_logg] != logg:
            if i_logg > 0 and i_logg < len(logg_grid):
                # Check which one is closer
                if abs(logg_grid[i_logg - 1] - logg) < abs(logg_grid[i_logg] - logg):
                    i_logg = i_logg - 1
            elif i_logg == len(logg_grid):
                i_logg = i_logg - 1

        i_meta = np.searchsorted(meta_grid, meta)
        if i_meta == len(meta_grid) or meta_grid[i_meta] != meta:
            if i_meta > 0 and i_meta < len(meta_grid):
                # Check which one is closer
                if abs(meta_grid[i_meta - 1] - meta) < abs(meta_grid[i_meta] - meta):
                    i_meta = i_meta - 1
            elif i_meta == len(meta_grid):
                i_meta = i_meta - 1

        # Load the SED
        model_path = os.path.join(stellar_model_dir, file_name)
        try:
            model_wavelengths, model_fluxes = load_sed(model_path, model_index)

            # Check if wavelength grids match
            if len(model_wavelengths) == n_lambda and np.allclose(
                model_wavelengths, wavelengths
            ):
                # Same grid, directly store values
                flux_cube[i_teff, i_logg, i_meta, :] = model_fluxes
                filled_map[i_teff, i_logg, i_meta] = True
            else:
                # Interpolate to the common wavelength grid
                flux_cube[i_teff, i_logg, i_meta, :] = np.interp(
                    wavelengths,
                    model_wavelengths,
                    model_fluxes,
                    left=model_fluxes[0],
                    right=model_fluxes[-1],  # Use edge values
                )
                filled_map[i_teff, i_logg, i_meta] = True
        except Exception as e:
            missing_points += 1
            print(f"Error processing {model_path}: {str(e)}")

        model_index += 1

    # Check for empty grid points
    empty_points = np.sum(~filled_map)
    if empty_points > 0:
        user_input = input("fill gaps? [True/False]: ").strip().lower()
        fill_empty = user_input in ("true", "t", "yes", "y", "1")

        if fill_empty:
            print(f"Warning: {empty_points} grid points are empty and need filling")

            # Simple fill algorithm: use nearest neighbor
            for i_teff in tqdm(range(len(teff_grid))):
                for i_logg in range(len(logg_grid)):
                    for i_meta in range(len(meta_grid)):
                        if not filled_map[i_teff, i_logg, i_meta]:
                            # Find nearest filled neighbor
                            best_dist = float("inf")
                            best_i, best_j, best_k = -1, -1, -1

                            for ii in range(len(teff_grid)):
                                for jj in range(len(logg_grid)):
                                    for kk in range(len(meta_grid)):
                                        if filled_map[ii, jj, kk]:
                                            # Simple Euclidean distance, could be improved
                                            d_teff = (
                                                teff_grid[i_teff] - teff_grid[ii]
                                            ) / 1000  # Scale for comparable magnitude
                                            d_logg = logg_grid[i_logg] - logg_grid[jj]
                                            d_meta = meta_grid[i_meta] - meta_grid[kk]
                                            dist = d_teff**2 + d_logg**2 + d_meta**2

                                            if dist < best_dist:
                                                best_dist = dist
                                                best_i, best_j, best_k = ii, jj, kk

                            if best_i >= 0:
                                # Copy data from nearest neighbor
                                flux_cube[i_teff, i_logg, i_meta, :] = flux_cube[
                                    best_i, best_j, best_k, :
                                ]
                                filled_map[i_teff, i_logg, i_meta] = True

    # Save the precomputed data
    save_binary_file(
        output_file, wavelengths, teff_grid, logg_grid, meta_grid, flux_cube
    )

    print(f"Precomputed data saved to {output_file}")
    print(f"Missing points during processing: {missing_points}")

    # Print some stats for verification
    print("\nData summary:")
    print(f"Wavelength range: {wavelengths[0]} to {wavelengths[-1]}")
    print(f"Teff range: {teff_grid[0]} to {teff_grid[-1]}")
    print(f"logg range: {logg_grid[0]} to {logg_grid[-1]}")
    print(f"Metallicity range: {meta_grid[0]} to {meta_grid[-1]}")
    print(f"Flux cube shape: {flux_cube.shape}")
    print(f"Flux cube size: {flux_cube.nbytes / (1024 * 1024):.2f} MB")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Precompute flux cube for stellar atmosphere models"
    )
    parser.add_argument(
        "--model_dir",
        type=str,
        required=True,
        help="Directory containing stellar model files",
    )
    parser.add_argument(
        "--output", type=str, default="flux_cube.bin", help="Output binary file"
    )

    args = parser.parse_args()
    precompute_flux_cube(args.model_dir, args.output)
