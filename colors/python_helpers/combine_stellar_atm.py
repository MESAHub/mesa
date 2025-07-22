#!/usr/bin/env python3
import argparse
import os
import shutil
import struct

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import simps
from scipy.interpolate import interp1d
from scipy.spatial import KDTree
from tqdm import tqdm

SIGMA = 5.670374419e-5  # erg s-1 cm-2 K-4


def renorm_to_sigmaT4(wl, flux, Teff):
    """wl in Å, flux in erg/cm²/s/Å  →  renormalised Fλ."""
    Fbol_model = simps(flux, wl)  # Å cancels
    Fbol_target = SIGMA * Teff**4
    return flux * (Fbol_target / Fbol_model)


def find_stellar_models(base_dir="../data/stellar_models/"):
    """Find all stellar model directories containing lookup tables."""
    model_dirs = []

    for item in os.listdir(base_dir):
        item_path = os.path.join(base_dir, item)
        if os.path.isdir(item_path):
            lookup_file = os.path.join(item_path, "lookup_table.csv")
            if os.path.exists(lookup_file):
                model_dirs.append((item, item_path))

    return model_dirs


def select_models_interactive(model_dirs):
    """Present model options to user and get selection."""
    print("\nAvailable stellar atmosphere models:")
    print("-" * 50)
    for idx, (name, path) in enumerate(model_dirs, start=1):
        # Count models in directory
        lookup_file = os.path.join(path, "lookup_table.csv")
        try:
            df = pd.read_csv(lookup_file, comment="#")
            n_models = len(df)
        except (FileNotFoundError, pd.errors.EmptyDataError, ValueError):
            n_models = "?"
        print(f"{idx}. {name} ({n_models} models)")

    print("\nEnter the numbers of models to combine (comma-separated):")
    print("Example: 1,3,5 or 'all' for all models")

    user_input = input("> ").strip()

    if user_input.lower() == "all":
        return model_dirs

    try:
        indices = [int(x.strip()) - 1 for x in user_input.split(",")]
        selected = [model_dirs[i] for i in indices if 0 <= i < len(model_dirs)]
        return selected
    except (ValueError, IndexError):
        print("Invalid input. Using all models.")
        return model_dirs


def load_model_data(model_path):
    """Load lookup table and extract parameter information."""
    lookup_file = os.path.join(model_path, "lookup_table.csv")

    # Read the CSV, handling the # comment character
    with open(lookup_file, "r") as f:
        # Read header line
        header = f.readline().strip()
        if header.startswith("#"):
            header = header[1:].strip()
        columns = [col.strip() for col in header.split(",")]

    # Read the data
    df = pd.read_csv(lookup_file, comment="#", names=columns, skiprows=1)

    # Extract parameters
    file_col = columns[0]  # Usually 'file_name' or 'filename'

    # Find parameter columns (case-insensitive)
    teff_col = None
    logg_col = None
    meta_col = None

    for col in columns:
        col_lower = col.lower()
        if "teff" in col_lower:
            teff_col = col
        elif "logg" in col_lower or "log(g)" in col_lower:
            logg_col = col
        elif "meta" in col_lower or "feh" in col_lower or "m/h" in col_lower:
            meta_col = col

    data = {
        "files": df[file_col].values,
        "teff": df[teff_col].values if teff_col else np.zeros(len(df)),
        "logg": df[logg_col].values if logg_col else np.zeros(len(df)),
        "meta": df[meta_col].values if meta_col else np.zeros(len(df)),
        "model_dir": model_path,
    }

    return data


def load_sed(filepath):
    return np.loadtxt(filepath, unpack=True)


def prepare_sed(filepath, teff):
    wl, flux = load_sed(filepath)
    # assume it's already in F_lambda, skip unit guessing for now
    flux = renorm_to_sigmaT4(wl, flux, teff)
    return wl, flux


def normalize_to_reference(target_data, reference_data, wavelength_grid):
    """
    Normalize target SEDs to match flux scaling of closest reference model.
    """
    print(
        f"  Normalizing to reference model: {
            os.path.basename(reference_data['model_dir'])
        }"
    )

    # Check if reference data is valid
    if len(reference_data["files"]) == 0:
        print("    Error: No reference files found!")
        return np.ones(len(target_data["files"]))

    # Build KDTree for reference models
    reference_params = np.column_stack(
        [reference_data["teff"], reference_data["logg"], reference_data["meta"]]
    )

    # Check for NaN/inf values and clean them
    valid_mask = np.isfinite(reference_params).all(axis=1)
    if not valid_mask.all():
        print(
            f"    Warning: Found {
                np.sum(~valid_mask)
            } invalid reference models, removing them"
        )
        reference_params = reference_params[valid_mask]
        # Also filter the corresponding data
        reference_files = np.array(reference_data["files"])[valid_mask]
        reference_teff = np.array(reference_data["teff"])[valid_mask]
        np.array(reference_data["logg"])[valid_mask]
        np.array(reference_data["meta"])[valid_mask]
    else:
        reference_files = reference_data["files"]
        reference_teff = reference_data["teff"]

    if len(reference_params) == 0:
        print("    Error: No valid reference models found!")
        return np.ones(len(target_data["files"]))

    print("    Reference model parameter ranges:")
    print(
        f"      Teff: {reference_params[:, 0].min():.1f} - {reference_params[:, 0].max():.1f} K"
    )
    print(
        f"      logg: {reference_params[:, 1].min():.2f} - {reference_params[:, 1].max():.2f}"
    )
    print(
        f"      [M/H]: {reference_params[:, 2].min():.2f} - {reference_params[:, 2].max():.2f}"
    )

    # Normalize parameters for better distance calculation
    # Handle case where range is zero
    teff_range = reference_params[:, 0].max() - reference_params[:, 0].min()
    logg_range = reference_params[:, 1].max() - reference_params[:, 1].min()
    meta_range = reference_params[:, 2].max() - reference_params[:, 2].min()

    print(
        f"    Parameter ranges: Teff={teff_range:.1f}, logg={logg_range:.2f}, meta={meta_range:.2f}"
    )

    if teff_range > 0:
        teff_norm = (reference_params[:, 0] - reference_params[:, 0].min()) / teff_range
    else:
        teff_norm = np.zeros(len(reference_params))
        print("    Warning: Teff range is zero, using constant normalization")

    if logg_range > 0:
        logg_norm = (reference_params[:, 1] - reference_params[:, 1].min()) / logg_range
    else:
        logg_norm = np.zeros(len(reference_params))
        print("    Warning: logg range is zero, using constant normalization")

    if meta_range > 0:
        meta_norm = (reference_params[:, 2] - reference_params[:, 2].min()) / meta_range
    else:
        meta_norm = np.zeros(len(reference_params))
        print("    Warning: [M/H] range is zero, using constant normalization")

    reference_params_norm = np.column_stack([teff_norm, logg_norm, meta_norm])

    # Final check for finite values
    if not np.isfinite(reference_params_norm).all():
        print("    Error: Normalized parameters contain NaN/inf values!")
        print(
            f"    Teff range: {reference_params[:, 0].min():.1f} - {reference_params[:, 0].max():.1f}"
        )
        print(
            f"    logg range: {reference_params[:, 1].min():.2f} - {reference_params[:, 1].max():.2f}"
        )
        print(
            f"    meta range: {reference_params[:, 2].min():.2f} - {reference_params[:, 2].max():.2f}"
        )
        return np.ones(len(target_data["files"]))

    # If all parameter ranges are zero, we can't use KDTree effectively
    if teff_range == 0 and logg_range == 0 and meta_range == 0:
        print("    Warning: All parameter ranges are zero, using simple matching")
        # Just use the first reference model for all targets
        reference_tree = None
    else:
        reference_tree = KDTree(reference_params_norm)

    norm_factors = []

    for i, (file, teff, logg, meta) in enumerate(
        zip(
            target_data["files"],
            target_data["teff"],
            target_data["logg"],
            target_data["meta"],
        )
    ):
        # Check for invalid target parameters
        if not (np.isfinite(teff) and np.isfinite(logg) and np.isfinite(meta)):
            print(
                f"    Warning: Invalid target parameters for {file}: Teff={teff}, logg={
                    logg
                }, meta={meta}"
            )
            norm_factors.append(1.0)
            continue

        # Normalize target parameters the same way
        if teff_range > 0:
            target_teff_norm = (teff - reference_params[:, 0].min()) / teff_range
        else:
            target_teff_norm = 0.0

        if logg_range > 0:
            target_logg_norm = (logg - reference_params[:, 1].min()) / logg_range
        else:
            target_logg_norm = 0.0

        if meta_range > 0:
            target_meta_norm = (meta - reference_params[:, 2].min()) / meta_range
        else:
            target_meta_norm = 0.0

        query = np.array([target_teff_norm, target_logg_norm, target_meta_norm])

        # Check if query is finite
        if not np.isfinite(query).all():
            print(f"    Warning: Invalid normalized query for {file}")
            norm_factors.append(1.0)
            continue

        # Find closest reference model
        if reference_tree is None:
            # Use first reference model if we can't do proper matching
            idx = 0
        else:
            dist, idx = reference_tree.query(query)

        # Load both SEDs
        file_target = os.path.join(target_data["model_dir"], file)
        file_ref = os.path.join(reference_data["model_dir"], reference_files[idx])

        try:
            wl_tgt, flux_tgt = prepare_sed(file_target, teff)
            wl_ref, flux_ref = prepare_sed(file_ref, reference_teff[idx])

            # Interpolate both to common wavelength grid
            # Use wavelength range that both models cover
            wl_min = max(wl_tgt.min(), wl_ref.min(), wavelength_grid.min())
            wl_max = min(wl_tgt.max(), wl_ref.max(), wavelength_grid.max())
            mask = (wavelength_grid >= wl_min) & (wavelength_grid <= wl_max)

            if mask.sum() < 100:  # Need enough overlap
                norm_factors.append(1.0)
                continue

            wl_common = wavelength_grid[mask]

            interp_tgt = interp1d(wl_tgt, flux_tgt, bounds_error=False, fill_value=0)
            interp_ref = interp1d(wl_ref, flux_ref, bounds_error=False, fill_value=0)

            flux_tgt_interp = interp_tgt(wl_common)
            flux_ref_interp = interp_ref(wl_common)

            # Compute integrated fluxes
            F_tgt = np.trapz(flux_tgt_interp, wl_common)
            F_ref = np.trapz(flux_ref_interp, wl_common)

            # Normalization factor
            factor = F_ref / F_tgt if F_tgt > 0 else 1.0

            # Sanity check - don't allow extreme normalizations
            if 0.01 < factor < 100:
                norm_factors.append(factor)
            else:
                norm_factors.append(1.0)

        except Exception as e:
            print(f"    Warning: could not normalize {file} → {e}")
            norm_factors.append(1.0)

    norm_factors = np.array(norm_factors)
    print(
        f"    Normalization factors: median={np.median(norm_factors):.2f}, range=[{np.min(norm_factors):.2f}, {np.max(norm_factors):.2f}]"
    )

    return norm_factors


def create_unified_grid(all_models_data):
    """Create unified parameter grids from all models."""
    # Collect all parameter values
    all_teff = []
    all_logg = []
    all_meta = []

    for data in all_models_data:
        all_teff.extend(data["teff"])
        all_logg.extend(data["logg"])
        all_meta.extend(data["meta"])

    # Get unique sorted values
    teff_grid = np.unique(np.sort(all_teff))
    logg_grid = np.unique(np.sort(all_logg))
    meta_grid = np.unique(np.sort(all_meta))

    # Remove values that are too close (within tolerance)
    def clean_grid(grid, tol=1e-6):
        if len(grid) <= 1:
            return grid
        cleaned = [grid[0]]
        for val in grid[1:]:
            if abs(val - cleaned[-1]) > tol:
                cleaned.append(val)
        return np.array(cleaned)

    teff_grid = clean_grid(teff_grid, tol=1.0)  # 1K tolerance for Teff
    logg_grid = clean_grid(logg_grid, tol=0.01)  # 0.01 dex for log g
    meta_grid = clean_grid(meta_grid, tol=0.01)  # 0.01 dex for [M/H]

    return teff_grid, logg_grid, meta_grid


def create_common_wavelength_grid(all_models_data, sample_size=20):
    """Create a common wavelength grid by sampling models."""
    print("\nAnalyzing wavelength coverage across models...")

    min_wave = float("inf")
    max_wave = 0
    resolutions = []

    for model_data in all_models_data:
        # Sample a few SEDs from this model
        n_sample = min(sample_size, len(model_data["files"]))
        indices = np.random.choice(len(model_data["files"]), n_sample, replace=False)

        for idx in indices:
            filepath = os.path.join(model_data["model_dir"], model_data["files"][idx])
            try:
                wavelengths, _ = load_sed(filepath)

                if len(wavelengths) > 10:
                    # Focus on optical/near-IR
                    mask = (wavelengths >= 3000) & (wavelengths <= 25000)
                    if mask.sum() > 10:
                        wl_subset = wavelengths[mask]
                        min_wave = min(min_wave, wl_subset.min())
                        max_wave = max(max_wave, wl_subset.max())

                        # Estimate resolution
                        resolution = np.median(np.diff(wl_subset))
                        resolutions.append(resolution)

            except (ValueError, IndexError, TypeError):
                continue

    # Create wavelength grid
    typical_resolution = np.median(resolutions) if resolutions else 50.0
    grid_resolution = max(50.0, typical_resolution * 2)

    n_points = int((max_wave - min_wave) / grid_resolution) + 1
    n_points = min(n_points, 5000)  # Cap at 5000 points

    wavelength_grid = np.linspace(min_wave, max_wave, n_points)

    print(f"  Wavelength range: {min_wave:.0f} - {max_wave:.0f} Å")
    print(f"  Grid points: {len(wavelength_grid)}")
    print(f"  Resolution: {grid_resolution:.1f} Å")

    return wavelength_grid


def identify_reference_model(all_models_data):
    """Identify which model to use as reference (prefer Kurucz)."""
    # Look for Kurucz model
    for i, model_data in enumerate(all_models_data):
        model_name = os.path.basename(model_data["model_dir"]).lower()
        if "kurucz" in model_name:
            print(
                f"Using {os.path.basename(model_data['model_dir'])} as reference model"
            )
            return i

    # If no Kurucz, use the first model
    print(
        f"No Kurucz model found, using {os.path.basename(all_models_data[0]['model_dir'])} as reference"
    )
    return 0


def build_combined_flux_cube(
    all_models_data, teff_grid, logg_grid, meta_grid, wavelength_grid
):
    """Build the combined flux cube from all models."""
    n_teff = len(teff_grid)
    n_logg = len(logg_grid)
    n_meta = len(meta_grid)
    n_lambda = len(wavelength_grid)

    # Initialize flux cube and tracking arrays
    flux_cube = np.zeros((n_teff, n_logg, n_meta, n_lambda))
    filled_map = np.zeros((n_teff, n_logg, n_meta), dtype=bool)
    source_map = (
        np.zeros((n_teff, n_logg, n_meta), dtype=int) - 1
    )  # Which model filled each point

    print(f"\nBuilding combined flux cube: {flux_cube.shape}")
    print(f"Memory requirement: {flux_cube.nbytes / (1024**2):.1f} MB")

    # Identify reference model
    reference_idx = identify_reference_model(all_models_data)
    reference_data = all_models_data[reference_idx]

    # Calculate normalization factors for each model
    normalization_factors = {}
    for model_idx, model_data in enumerate(all_models_data):
        if model_idx == reference_idx:
            # Reference model doesn't need normalization
            normalization_factors[model_idx] = np.ones(len(model_data["files"]))
        else:
            # Normalize to reference
            normalization_factors[model_idx] = normalize_to_reference(
                model_data, reference_data, wavelength_grid
            )

    # Build KDTree for fast nearest neighbor searches
    # Normalize parameters for distance calculation
    teff_norm = (teff_grid - teff_grid.min()) / (teff_grid.max() - teff_grid.min())
    logg_norm = (logg_grid - logg_grid.min()) / (logg_grid.max() - logg_grid.min())
    meta_norm = (meta_grid - meta_grid.min()) / (meta_grid.max() - meta_grid.min())

    # Process each model
    for model_idx, model_data in enumerate(all_models_data):
        model_name = os.path.basename(model_data["model_dir"])
        print(f"\nProcessing {model_name}...")

        norm_factors = normalization_factors[model_idx]

        for i, (file, teff, logg, meta) in enumerate(
            tqdm(
                zip(
                    model_data["files"],
                    model_data["teff"],
                    model_data["logg"],
                    model_data["meta"],
                ),
                total=len(model_data["files"]),
                desc=f"Model {model_idx + 1}",
            )
        ):
            # Find grid indices
            i_teff = np.searchsorted(teff_grid, teff)
            i_logg = np.searchsorted(logg_grid, logg)
            i_meta = np.searchsorted(meta_grid, meta)

            # Clip to valid range
            i_teff = np.clip(i_teff, 0, n_teff - 1)
            i_logg = np.clip(i_logg, 0, n_logg - 1)
            i_meta = np.clip(i_meta, 0, n_meta - 1)

            # Load and interpolate SED
            filepath = os.path.join(model_data["model_dir"], file)
            try:
                model_wavelengths, model_fluxes = prepare_sed(filepath, teff)

                # Apply normalization factor
                model_fluxes *= norm_factors[i]

                # Interpolate to common grid (in log space for flux)
                log_fluxes = np.log10(np.maximum(model_fluxes, 1e-50))
                log_interpolated = np.interp(
                    wavelength_grid,
                    model_wavelengths,
                    log_fluxes,
                    left=log_fluxes[0],
                    right=log_fluxes[-1],
                )
                interpolated_flux = 10**log_interpolated

                # Store in cube
                flux_cube[i_teff, i_logg, i_meta, :] = interpolated_flux
                filled_map[i_teff, i_logg, i_meta] = True
                source_map[i_teff, i_logg, i_meta] = model_idx

            except Exception:
                continue

    # Fill gaps using nearest neighbor interpolation
    empty_points = np.sum(~filled_map)
    if empty_points > 0:
        print(f"\nFilling {empty_points} empty grid points...")

        # Get filled points
        filled_indices = np.where(filled_map)
        filled_points = np.column_stack(
            [
                teff_norm[filled_indices[0]],
                logg_norm[filled_indices[1]],
                meta_norm[filled_indices[2]],
            ]
        )

        # Build KDTree
        tree = KDTree(filled_points)

        # Fill empty points
        for i_teff in range(n_teff):
            for i_logg in range(n_logg):
                for i_meta in range(n_meta):
                    if not filled_map[i_teff, i_logg, i_meta]:
                        # Find nearest filled point
                        query_point = [
                            teff_norm[i_teff],
                            logg_norm[i_logg],
                            meta_norm[i_meta],
                        ]
                        dist, idx = tree.query(query_point, k=1)

                        # Copy flux from nearest neighbor
                        src_i = filled_indices[0][idx]
                        src_j = filled_indices[1][idx]
                        src_k = filled_indices[2][idx]

                        flux_cube[i_teff, i_logg, i_meta, :] = flux_cube[
                            src_i, src_j, src_k, :
                        ]
                        filled_map[i_teff, i_logg, i_meta] = True
                        source_map[i_teff, i_logg, i_meta] = source_map[
                            src_i, src_j, src_k
                        ]

    return flux_cube, source_map


def save_combined_data(
    output_dir,
    teff_grid,
    logg_grid,
    meta_grid,
    wavelength_grid,
    flux_cube,
    all_models_data,
):
    """Save the combined data and create unified lookup table."""
    os.makedirs(output_dir, exist_ok=True)

    # Save binary flux cube
    binary_file = os.path.join(output_dir, "flux_cube.bin")
    with open(binary_file, "wb") as f:
        # Write dimensions
        f.write(
            struct.pack(
                "4i",
                len(teff_grid),
                len(logg_grid),
                len(meta_grid),
                len(wavelength_grid),
            )
        )

        # Write grid arrays
        teff_grid.astype(np.float64).tofile(f)
        logg_grid.astype(np.float64).tofile(f)
        meta_grid.astype(np.float64).tofile(f)
        wavelength_grid.astype(np.float64).tofile(f)

        # Write flux cube
        flux_cube.astype(np.float64).tofile(f)

    print(f"\nSaved binary flux cube to: {binary_file}")

    # Create combined lookup table
    lookup_data = []
    file_counter = 0
    copied_files = 0

    print("\nCopying SED files to combined directory...")

    for model_data in tqdm(all_models_data, desc="Copying models"):
        model_name = os.path.basename(model_data["model_dir"])
        for i, (orig_file, teff, logg, meta) in enumerate(
            zip(
                model_data["files"],
                model_data["teff"],
                model_data["logg"],
                model_data["meta"],
            )
        ):
            # Create new filename that includes source model
            new_filename = f"{model_name}_{file_counter:06d}.txt"
            lookup_data.append(
                {
                    "file_name": new_filename,
                    "teff": teff,
                    "logg": logg,
                    "meta": meta,
                    "source_model": model_name,
                    "original_file": orig_file,
                }
            )

            # Copy the actual file instead of creating a symlink
            src_path = os.path.join(model_data["model_dir"], orig_file)
            dst_path = os.path.join(output_dir, new_filename)

            if os.path.exists(src_path) and not os.path.exists(dst_path):
                try:
                    # copy2 preserves metadata
                    shutil.copy2(src_path, dst_path)
                    copied_files += 1
                except Exception as e:
                    print(f"\nWarning: Could not copy {src_path}: {e}")

            file_counter += 1

    print(f"Copied {copied_files} SED files to combined directory")

    # Save lookup table
    lookup_df = pd.DataFrame(lookup_data)
    lookup_file = os.path.join(output_dir, "lookup_table.csv")

    with open(lookup_file, "w") as f:
        f.write("#file_name, teff, logg, meta, source_model, original_file\n")
        lookup_df.to_csv(f, index=False, header=False)

    print(f"Saved combined lookup table to: {lookup_file}")
    print(f"Total models in combined set: {len(lookup_df)}")

    # Calculate and display disk usage
    total_size = 0
    for f in os.listdir(output_dir):
        if f.endswith(".txt"):
            total_size += os.path.getsize(os.path.join(output_dir, f))

    print(f"Total disk space used: {total_size / (1024**3):.2f} GB")

    return lookup_df


def visualize_parameter_space(
    teff_grid, logg_grid, meta_grid, source_map, all_models_data, output_dir
):
    """Create visualizations of the parameter space coverage."""
    print("\nCreating parameter space visualizations...")

    # Get model names
    model_names = [os.path.basename(data["model_dir"]) for data in all_models_data]

    # Create color map for models
    colors = plt.cm.tab10(np.linspace(0, 1, len(model_names)))

    fig = plt.figure(figsize=(18, 12))

    # 1. 3D scatter plot of filled points
    ax1 = fig.add_subplot(2, 3, 1, projection="3d")

    # Get coordinates of filled points for each model
    for model_idx, model_name in enumerate(model_names):
        mask = source_map == model_idx
        if np.any(mask):
            indices = np.where(mask)
            ax1.scatter(
                teff_grid[indices[0]],
                logg_grid[indices[1]],
                meta_grid[indices[2]],
                c=[colors[model_idx]],
                label=model_name,
                alpha=0.6,
                s=20,
            )

    ax1.set_xlabel("Teff (K)")
    ax1.set_ylabel("log g")
    ax1.set_zlabel("[M/H]")
    ax1.set_title("3D Parameter Space Coverage")
    ax1.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

    # 2. Teff vs log g projection
    ax2 = fig.add_subplot(2, 3, 2)
    for model_idx, model_name in enumerate(model_names):
        data = all_models_data[model_idx]
        ax2.scatter(
            data["teff"],
            data["logg"],
            c=[colors[model_idx]],
            label=model_name,
            alpha=0.6,
            s=10,
        )
    ax2.set_xlabel("Teff (K)")
    ax2.set_ylabel("log g")
    ax2.set_title("Teff vs log g")
    ax2.invert_yaxis()
    ax2.grid(True, alpha=0.3)

    # 3. Teff vs [M/H] projection
    ax3 = fig.add_subplot(2, 3, 3)
    for model_idx, model_name in enumerate(model_names):
        data = all_models_data[model_idx]
        ax3.scatter(
            data["teff"],
            data["meta"],
            c=[colors[model_idx]],
            label=model_name,
            alpha=0.6,
            s=10,
        )
    ax3.set_xlabel("Teff (K)")
    ax3.set_ylabel("[M/H]")
    ax3.set_title("Teff vs [M/H]")
    ax3.grid(True, alpha=0.3)

    # 4. log g vs [M/H] projection
    ax4 = fig.add_subplot(2, 3, 4)
    for model_idx, model_name in enumerate(model_names):
        data = all_models_data[model_idx]
        ax4.scatter(
            data["logg"],
            data["meta"],
            c=[colors[model_idx]],
            label=model_name,
            alpha=0.6,
            s=10,
        )
    ax4.set_xlabel("log g")
    ax4.set_ylabel("[M/H]")
    ax4.set_title("log g vs [M/H]")
    ax4.grid(True, alpha=0.3)

    # 5. Normalisation-verification plot: λ²Fλ (should overlay after normalization)
    ax5 = fig.add_subplot(2, 3, 5)
    target = (5777, 4.44, 0.0)  # solar-like reference

    for idx, model_name in enumerate(model_names):
        data = all_models_data[idx]

        # Find closest SED to solar point
        valid_mask = (
            np.isfinite(data["teff"])
            & np.isfinite(data["logg"])
            & np.isfinite(data["meta"])
        )
        if not valid_mask.any():
            continue

        valid_teff = np.array(data["teff"])[valid_mask]
        valid_logg = np.array(data["logg"])[valid_mask]
        valid_meta = np.array(data["meta"])[valid_mask]
        valid_files = np.array(data["files"])[valid_mask]

        dist = (
            ((valid_teff - target[0]) / 1000) ** 2
            + (valid_logg - target[1]) ** 2
            + (valid_meta - target[2]) ** 2
        )

        if len(dist) > 0:
            j = np.argmin(dist)
            fpath = os.path.join(data["model_dir"], valid_files[j])

            # Quick guard against XML/FITS/binary junk
            if not fpath.lower().endswith((".txt", ".dat", ".sed")):
                continue
            try:
                # Check if file exists and is readable
                if not os.path.exists(fpath):
                    continue

                with open(fpath, "rb") as fh:
                    # skip if the first non-whitespace byte is '<' (likely XML)
                    first = fh.read(256).lstrip()[:1]
                if first == b"<":
                    continue

                wl, fl = prepare_sed(fpath, valid_teff[j])
                mask = (wl > 3000) & (wl < 10000)
                if mask.sum() > 20:
                    ax5.plot(
                        wl[mask], wl[mask] ** 2 * fl[mask], label=model_name, alpha=0.8
                    )
            except Exception as e:
                print(f"  ⚠ Skipping {model_name}: {e}")
                continue

    ax5.set_xscale("log")
    ax5.set_yscale("log")
    ax5.set_xlabel("Wavelength (Å)")
    ax5.set_ylabel(r"$λ^2 F_λ$ (arb. units)")
    ax5.set_title("Normalisation Check (solar-like)")
    ax5.grid(True, alpha=0.3)
    ax5.legend(fontsize=8)

    # 6. Grid density heatmap
    ax6 = fig.add_subplot(2, 3, 6)

    # Count models per grid cell (marginalized over metallicity)
    density_map = np.zeros((len(teff_grid), len(logg_grid)))
    for i in range(len(teff_grid)):
        for j in range(len(logg_grid)):
            # Count unique models at this Teff, log g
            models_here = source_map[i, j, :]
            unique_models = len(np.unique(models_here[models_here >= 0]))
            density_map[i, j] = unique_models

    im = ax6.imshow(
        density_map.T,
        origin="lower",
        aspect="auto",
        extent=[teff_grid.min(), teff_grid.max(), logg_grid.min(), logg_grid.max()],
        cmap="YlOrRd",
    )
    ax6.set_xlabel("Teff (K)")
    ax6.set_ylabel("log g")
    ax6.set_title("Model Density (# of models per grid point)")
    plt.colorbar(im, ax=ax6)

    plt.tight_layout()
    plot_file = os.path.join(output_dir, "parameter_space_visualization.png")
    plt.savefig(plot_file, dpi=150, bbox_inches="tight")
    print(f"Saved visualization to: {plot_file}")

    # Print summary statistics
    print("\n" + "=" * 60)
    print("COMBINED MODEL STATISTICS")
    print("=" * 60)
    print(f"Total number of source models: {len(all_models_data)}")
    print(
        f"Total grid points: {len(teff_grid)} × {len(logg_grid)} × {len(meta_grid)} = "
        f"{len(teff_grid) * len(logg_grid) * len(meta_grid):,}"
    )
    print("\nParameter ranges:")
    print(f"  Teff: {teff_grid.min():.0f} - {teff_grid.max():.0f} K")
    print(f"  log g: {logg_grid.min():.2f} - {logg_grid.max():.2f}")
    print(f"  [M/H]: {meta_grid.min():.2f} - {meta_grid.max():.2f}")

    print("\nPer-model contributions:")
    for model_idx, model_name in enumerate(model_names):
        n_points = np.sum(source_map == model_idx)
        pct = 100 * n_points / source_map.size
        print(f"  {model_name}: {n_points:,} grid points ({pct:.1f}%)")


def main():
    parser = argparse.ArgumentParser(
        description="Combine multiple stellar atmosphere models"
    )
    parser.add_argument(
        "--data_dir",
        type=str,
        default="../data/stellar_models/",
        help="Base directory containing stellar models",
    )
    parser.add_argument(
        "--output", type=str, default="combined_models", help="Output directory name"
    )
    parser.add_argument(
        "--interactive",
        action="store_true",
        default=True,
        help="Interactive model selection",
    )

    args = parser.parse_args()

    # Find available models
    model_dirs = find_stellar_models(args.data_dir)

    if not model_dirs:
        print(f"No stellar models found in {args.data_dir}")
        return

    # Select models to combine
    if args.interactive:
        selected_models = select_models_interactive(model_dirs)
    else:
        selected_models = model_dirs

    if not selected_models:
        print("No models selected.")
        return

    print(f"\nSelected {len(selected_models)} models to combine:")
    for name, path in selected_models:
        print(f"  - {name}")

    # Load all model data
    print("\nLoading model data...")
    all_models_data = []
    for name, path in selected_models:
        print(f"  Loading {name}...")
        data = load_model_data(path)
        all_models_data.append(data)

    # Create unified grids
    print("\nCreating unified parameter grids...")
    teff_grid, logg_grid, meta_grid = create_unified_grid(all_models_data)
    wavelength_grid = create_common_wavelength_grid(all_models_data)

    # Build combined flux cube
    flux_cube, source_map = build_combined_flux_cube(
        all_models_data, teff_grid, logg_grid, meta_grid, wavelength_grid
    )

    # Save combined data
    output_dir = os.path.join(args.data_dir, args.output)
    save_combined_data(
        output_dir,
        teff_grid,
        logg_grid,
        meta_grid,
        wavelength_grid,
        flux_cube,
        all_models_data,
    )

    # Create visualizations
    visualize_parameter_space(
        teff_grid, logg_grid, meta_grid, source_map, all_models_data, output_dir
    )

    print(f"\nSuccessfully combined {len(selected_models)} stellar atmosphere models!")
    print(f"Output saved to: {output_dir}")
    print("You can now use this combined model in MESA by setting:")
    print(f"  stellar_atm = '{output_dir}/'")


if __name__ == "__main__":
    main()
