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
    
    return file_names, np.array(teff_values), np.array(logg_values), np.array(meta_values)

def get_unique_sorted(values, tolerance=1e-8):
    """Get sorted unique values with tolerance."""
    sorted_vals = np.sort(values)
    unique_indices = [0]
    
    for i in range(1, len(sorted_vals)):
        if abs(sorted_vals[i] - sorted_vals[unique_indices[-1]]) > tolerance:
            unique_indices.append(i)
            
    return sorted_vals[unique_indices]

def save_binary_file(output_file, wavelengths, teff_grid, logg_grid, meta_grid, flux_cube):
    """Save the precomputed data in a simple binary format that Fortran can read."""
    with open(output_file, 'wb') as f:
        # Write dimensions
        f.write(struct.pack('4i', len(teff_grid), len(logg_grid), len(meta_grid), len(wavelengths)))
        
        # Write grid arrays
        teff_grid.astype(np.float64).tofile(f)
        logg_grid.astype(np.float64).tofile(f)
        meta_grid.astype(np.float64).tofile(f)
        wavelengths.astype(np.float64).tofile(f)
        
        # Write flux cube
        flux_cube.astype(np.float64).tofile(f)

def precompute_flux_cube(stellar_model_dir, output_file):
    """Precompute the 3D flux cube for all wavelengths."""
    # Load the lookup table
    lookup_file = os.path.join(stellar_model_dir, 'lookup_table.csv')
    file_names, teff_values, logg_values, meta_values = load_lookup_table(lookup_file)
    
    print(f"Loaded {len(file_names)} model files from lookup table")
    
    # Get unique parameter values
    teff_grid = get_unique_sorted(teff_values)
    logg_grid = get_unique_sorted(logg_values)
    meta_grid = get_unique_sorted(meta_values)
    
    print(f"Unique Teff values: {len(teff_grid)}")
    print(f"Unique logg values: {len(logg_grid)}")
    print(f"Unique metallicity values: {len(meta_grid)}")
    
    # Load first SED to get wavelength grid
    first_file = os.path.join(stellar_model_dir, file_names[0])
    wavelengths, _ = load_sed(first_file, 0)
    n_lambda = len(wavelengths)
    
    print(f"Wavelength grid has {n_lambda} points")
    
    # Create the flux cube
    flux_cube = np.zeros((len(teff_grid), len(logg_grid), len(meta_grid), n_lambda))
    
    # Process all model files
    model_index = 0
    for file_name, teff, logg, meta in tqdm(zip(file_names, teff_values, logg_values, meta_values), 
                                          total=len(file_names), desc="Processing models"):
        # Find the indices in the grids
        i_teff = np.searchsorted(teff_grid, teff)
        if i_teff == len(teff_grid):
            i_teff -= 1
            
        i_logg = np.searchsorted(logg_grid, logg)
        if i_logg == len(logg_grid):
            i_logg -= 1
            
        i_meta = np.searchsorted(meta_grid, meta)
        if i_meta == len(meta_grid):
            i_meta -= 1
        
        # Load the SED
        model_path = os.path.join(stellar_model_dir, file_name)
        try:
            model_wavelengths, model_fluxes = load_sed(model_path, model_index)
            
            # Check if wavelength grids match
            if len(model_wavelengths) == n_lambda and np.allclose(model_wavelengths, wavelengths):
                # Same grid, directly store values
                flux_cube[i_teff, i_logg, i_meta, :] = model_fluxes
            else:
                # Interpolate to the common wavelength grid
                flux_cube[i_teff, i_logg, i_meta, :] = np.interp(
                    wavelengths, model_wavelengths, model_fluxes, 
                    left=0.0, right=0.0
                )
        except Exception as e:
            print(f"Error processing {model_path}: {str(e)}")
        
        model_index += 1
    
    # Save the precomputed data
    save_binary_file(output_file, wavelengths, teff_grid, logg_grid, meta_grid, flux_cube)
    
    print(f"Precomputed data saved to {output_file}")
    
    # Print some stats for verification
    print("\nData summary:")
    print(f"Wavelength range: {wavelengths[0]} to {wavelengths[-1]}")
    print(f"Teff range: {teff_grid[0]} to {teff_grid[-1]}")
    print(f"logg range: {logg_grid[0]} to {logg_grid[-1]}")
    print(f"Metallicity range: {meta_grid[0]} to {meta_grid[-1]}")
    print(f"Flux cube shape: {flux_cube.shape}")
    print(f"Flux cube size: {flux_cube.nbytes / (1024*1024):.2f} MB")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Precompute flux cube for stellar atmosphere models')
    parser.add_argument('--model_dir', type=str, required=True, help='Directory containing stellar model files')
    parser.add_argument('--output', type=str, default='flux_cube.bin', help='Output binary file')
    
    args = parser.parse_args()
    precompute_flux_cube(args.model_dir, args.output)