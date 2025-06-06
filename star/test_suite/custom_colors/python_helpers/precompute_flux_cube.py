#!/usr/bin/env python3
import os
import numpy as np
from tqdm import tqdm
import argparse
import glob
import struct

def load_sed(filepath, index):
    """Load a spectral energy distribution file."""
    wavelengths = []
    fluxes = []
    
    with open(filepath, 'r') as f:
        # Skip header lines (comments)
        for line in f:
            if line.strip() and not line.startswith('#'):
                break
        
        # Read wavelength and flux values
        try:
            values = line.strip().split()
            wavelengths.append(float(values[0]))
            fluxes.append(float(values[1]))
        except:
            pass  # Skip any malformed lines
            
        for line in f:
            try:
                values = line.strip().split()
                wavelengths.append(float(values[0]))
                fluxes.append(float(values[1]))
            except:
                pass  # Skip any malformed lines
    
    return np.array(wavelengths), np.array(fluxes)

def load_lookup_table(lookup_file):
    """Load the lookup table with model parameters."""
    file_names = []
    teff_values = []
    logg_values = []
    meta_values = []
    
    with open(lookup_file, 'r') as f:
        # Skip header line
        header = f.readline().strip().split(',')
        
        # Find column indices
        file_col = 0  # Assume first column is filename
        teff_col = None
        logg_col = None
        meta_col = None
        
        for i, col in enumerate(header):
            col = col.strip().lower()
            if col == 'teff':
                teff_col = i
            elif col == 'logg':
                logg_col = i
            elif col in ['meta', 'feh']:
                meta_col = i
                
        # Read data rows
        for line in f:
            if line.strip():
                values = line.strip().split(',')
                if len(values) <= max(file_col, teff_col or 0, logg_col or 0, meta_col or 0):
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
    
    print(f"Column indices found - teff: {teff_col}, logg: {logg_col}, meta: {meta_col}")
    return file_names, np.array(teff_values), np.array(logg_values), np.array(meta_values)


def find_grid_index(grid, value):
    """Find the closest grid index for a value."""
    idx = np.searchsorted(grid, value)
    if idx == len(grid):
        return idx - 1
    elif idx == 0:
        return 0
    else:
        # Choose the closer of the two
        if abs(grid[idx-1] - value) < abs(grid[idx] - value):
            return idx - 1
        else:
            return idx



def get_unique_sorted(values, tolerance=1e-8):
    """Get sorted unique values with tolerance."""
    sorted_vals = np.sort(values)
    unique_indices = [0]
    
    for i in range(1, len(sorted_vals)):
        if abs(sorted_vals[i] - sorted_vals[unique_indices[-1]]) > tolerance:
            unique_indices.append(i)
            
    return sorted_vals[unique_indices]

def save_binary_file(output_file, wavelengths, teff_grid, logg_grid, meta_grid, flux_cube):
    with open(output_file, 'wb') as f:
        # Write dimensions
        f.write(struct.pack('4i', len(teff_grid), len(logg_grid), len(meta_grid), len(wavelengths)))
        
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
        with open(debug_file, 'w') as df:
            df.write("# Teff grid\n")
            for val in teff_grid:
                df.write(f"{val}\n")
            df.write("\n# Logg grid\n")
            for val in logg_grid:
                df.write(f"{val}\n")
            df.write("\n# Metallicity grid\n")
            for val in meta_grid:
                df.write(f"{val}\n")

def create_reasonable_wavelength_grid(stellar_model_dir, file_names):
    """Create a reasonable wavelength grid focused on optical/near-IR where filters typically work."""
    
    # Most photometric filters work in these ranges:
    # UV: 200-400 nm (2000-4000 Å)
    # Optical: 400-1000 nm (4000-10000 Å)  
    # Near-IR: 1-5 μm (10000-50000 Å)
    
    # Sample a few models to see their actual resolution
    print("Analyzing resolution of stellar models...")
    
    resolutions = []
    min_wavelength = float('inf')
    max_wavelength = 0
    
    for i, file_name in enumerate(file_names[:20]):  # Sample first 20
        model_path = os.path.join(stellar_model_dir, file_name)
        try:
            wavelengths, _ = load_sed(model_path, i)
            
            # Focus on optical/near-IR range where most filters work
            mask = (wavelengths >= 3000) & (wavelengths <= 25000)  # 300 nm to 2.5 μm
            if mask.sum() > 10:
                wl_subset = wavelengths[mask]
                resolution = np.median(np.diff(wl_subset))
                resolutions.append(resolution)
                
                min_wavelength = min(min_wavelength, wl_subset.min())
                max_wavelength = max(max_wavelength, wl_subset.max())
                
                if i < 5:
                    print(f"  Model {file_name}: resolution ≈ {resolution:.1f} Å in optical range")
                    
        except Exception as e:
            print(f"Error reading {model_path}: {e}")
            continue
    
    if not resolutions:
        print("Warning: No models found with optical data. Using defaults.")
        min_wavelength = 3000
        max_wavelength = 25000
        typical_resolution = 50.0
    else:
        typical_resolution = np.median(resolutions)
    
    print(f"Typical model resolution: {typical_resolution:.1f} Å")
    print(f"Optical wavelength range: [{min_wavelength:.0f}, {max_wavelength:.0f}] Å")
    
    # Create a grid with resolution similar to the models, but not finer
    # Use 2x the typical resolution to avoid over-sampling
    grid_resolution = max(50.0, typical_resolution * 2)  # At least 50 Å spacing
    
    n_points = int((max_wavelength - min_wavelength) / grid_resolution) + 1
    n_points = min(n_points, 5000)  # Cap at 5000 points maximum
    
    common_wavelengths = np.linspace(min_wavelength, max_wavelength, n_points)
    
    print(f"Created wavelength grid:")
    print(f"  Range: {min_wavelength:.0f} - {max_wavelength:.0f} Å")
    print(f"  Points: {len(common_wavelengths)}")
    print(f"  Resolution: {grid_resolution:.1f} Å")
    
    return common_wavelengths


def precompute_flux_cube_fixed(stellar_model_dir, output_file):
    """Precompute flux cube with reasonable wavelength sampling."""
    # Load the lookup table
    lookup_file = os.path.join(stellar_model_dir, 'lookup_table.csv')
    file_names, teff_values, logg_values, meta_values = load_lookup_table(lookup_file)
    
    print(f"Loaded {len(file_names)} model files from lookup table")
    
    # Get unique parameter values
    teff_grid = get_unique_sorted(teff_values)
    logg_grid = get_unique_sorted(logg_values)
    meta_grid = get_unique_sorted(meta_values)
    
    print(f"Parameter grid sizes:")
    print(f"  Teff: {len(teff_grid)} points ({teff_grid[0]:.0f} - {teff_grid[-1]:.0f} K)")
    print(f"  log g: {len(logg_grid)} points ({logg_grid[0]:.1f} - {logg_grid[-1]:.1f})")
    print(f"  [M/H]: {len(meta_grid)} points ({meta_grid[0]:.1f} - {meta_grid[-1]:.1f})")
    
    # Create reasonable wavelength grid (focused on optical/near-IR)
    wavelengths = create_reasonable_wavelength_grid(stellar_model_dir, file_names)
    n_lambda = len(wavelengths)
    
    # Create the flux cube and tracking map
    flux_cube = np.zeros((len(teff_grid), len(logg_grid), len(meta_grid), n_lambda))
    filled_map = np.zeros((len(teff_grid), len(logg_grid), len(meta_grid)), dtype=bool)
    
    print(f"Flux cube shape: {flux_cube.shape}")
    print(f"Memory usage: {flux_cube.nbytes / (1024*1024):.1f} MB")
    
    # Process all model files
    missing_points = 0
    for file_name, teff, logg, meta in tqdm(zip(file_names, teff_values, logg_values, meta_values), 
                                          total=len(file_names), desc="Processing models"):
        # Find grid indices
        i_teff = find_grid_index(teff_grid, teff)
        i_logg = find_grid_index(logg_grid, logg) 
        i_meta = find_grid_index(meta_grid, meta)
        
        # Load the SED
        model_path = os.path.join(stellar_model_dir, file_name)
        try:
            model_wavelengths, model_fluxes = load_sed(model_path, 0)
            
            # Interpolate to common wavelength grid
            # Use log interpolation for flux since it spans many orders of magnitude
            log_fluxes = np.log10(np.maximum(model_fluxes, 1e-50))  # Avoid log(0)
            
            log_interpolated = np.interp(
                wavelengths, model_wavelengths, log_fluxes,
                left=log_fluxes[0], right=log_fluxes[-1]
            )
            
            interpolated_flux = 10**log_interpolated
            
            flux_cube[i_teff, i_logg, i_meta, :] = interpolated_flux
            filled_map[i_teff, i_logg, i_meta] = True
            
        except Exception as e:
            missing_points += 1
            if missing_points < 10:  # Only print first few errors
                print(f"Error processing {model_path}: {str(e)}")
            continue
    
    # Fill empty grid points using nearest neighbor
    empty_points = np.sum(~filled_map)
    if empty_points > 0:
        print(f"Filling {empty_points} empty grid points...")
        
        for i_teff in range(len(teff_grid)):
            for i_logg in range(len(logg_grid)):
                for i_meta in range(len(meta_grid)):
                    if not filled_map[i_teff, i_logg, i_meta]:
                        # Find nearest filled neighbor
                        best_dist = float('inf')
                        best_i, best_j, best_k = -1, -1, -1
                        
                        for ii in range(len(teff_grid)):
                            for jj in range(len(logg_grid)):
                                for kk in range(len(meta_grid)):
                                    if filled_map[ii, jj, kk]:
                                        # Normalized distance
                                        d_teff = (teff_grid[i_teff] - teff_grid[ii])/1000
                                        d_logg = logg_grid[i_logg] - logg_grid[jj]
                                        d_meta = meta_grid[i_meta] - meta_grid[kk]
                                        dist = d_teff**2 + d_logg**2 + d_meta**2
                                        
                                        if dist < best_dist:
                                            best_dist = dist
                                            best_i, best_j, best_k = ii, jj, kk
                        
                        if best_i >= 0:
                            flux_cube[i_teff, i_logg, i_meta, :] = flux_cube[best_i, best_j, best_k, :]
                            filled_map[i_teff, i_logg, i_meta] = True
    
    # Save the precomputed data
    save_binary_file(output_file, wavelengths, teff_grid, logg_grid, meta_grid, flux_cube)
    
    print(f"\nPrecomputed data saved to {output_file}")
    print(f"Statistics:")
    print(f"  Models processed: {len(file_names) - missing_points}/{len(file_names)}")
    print(f"  Empty points filled: {empty_points}")
    print(f"  Final cube shape: {flux_cube.shape}")
    print(f"  File size: {flux_cube.nbytes / (1024*1024):.1f} MB")
    
    # Quick sanity check
    print(f"  Flux range: {flux_cube.min():.2e} - {flux_cube.max():.2e}")
    print(f"  Zero values: {np.sum(flux_cube == 0)}/{flux_cube.size}")
    print(f"  NaN values: {np.sum(np.isnan(flux_cube))}")


# Usage - replace the function call in your script:
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Precompute flux cube for stellar atmosphere models')
    parser.add_argument('--model_dir', type=str, required=True, help='Directory containing stellar model files')
    parser.add_argument('--output', type=str, default='flux_cube.bin', help='Output binary file')
    
    args = parser.parse_args()
    precompute_flux_cube_fixed(args.model_dir, args.output)  # Use the fixed version