#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse
from tqdm import tqdm
import struct
from scipy.spatial import KDTree
import seaborn as sns

def find_stellar_models(base_dir='../../data/stellar_models/'):
    """Find all stellar model directories containing lookup tables."""
    model_dirs = []
    
    for item in os.listdir(base_dir):
        item_path = os.path.join(base_dir, item)
        if os.path.isdir(item_path):
            lookup_file = os.path.join(item_path, 'lookup_table.csv')
            if os.path.exists(lookup_file):
                model_dirs.append((item, item_path))
    
    return model_dirs

def select_models_interactive(model_dirs):
    """Present model options to user and get selection."""
    print("\nAvailable stellar atmosphere models:")
    print("-" * 50)
    for idx, (name, path) in enumerate(model_dirs, start=1):
        # Count models in directory
        lookup_file = os.path.join(path, 'lookup_table.csv')
        try:
            df = pd.read_csv(lookup_file, comment='#')
            n_models = len(df)
        except:
            n_models = "?"
        print(f"{idx}. {name} ({n_models} models)")
    
    print("\nEnter the numbers of models to combine (comma-separated):")
    print("Example: 1,3,5 or 'all' for all models")
    
    user_input = input("> ").strip()
    
    if user_input.lower() == 'all':
        return model_dirs
    
    try:
        indices = [int(x.strip()) - 1 for x in user_input.split(',')]
        selected = [model_dirs[i] for i in indices if 0 <= i < len(model_dirs)]
        return selected
    except:
        print("Invalid input. Using all models.")
        return model_dirs

def load_model_data(model_path):
    """Load lookup table and extract parameter information."""
    lookup_file = os.path.join(model_path, 'lookup_table.csv')
    
    # Read the CSV, handling the # comment character
    with open(lookup_file, 'r') as f:
        # Read header line
        header = f.readline().strip()
        if header.startswith('#'):
            header = header[1:].strip()
        columns = [col.strip() for col in header.split(',')]
    
    # Read the data
    df = pd.read_csv(lookup_file, comment='#', names=columns, skiprows=1)
    
    # Extract parameters
    file_col = columns[0]  # Usually 'file_name' or 'filename'
    
    # Find parameter columns (case-insensitive)
    teff_col = None
    logg_col = None
    meta_col = None
    
    for col in columns:
        col_lower = col.lower()
        if 'teff' in col_lower:
            teff_col = col
        elif 'logg' in col_lower or 'log(g)' in col_lower:
            logg_col = col
        elif 'meta' in col_lower or 'feh' in col_lower or 'm/h' in col_lower:
            meta_col = col
    
    data = {
        'files': df[file_col].values,
        'teff': df[teff_col].values if teff_col else np.zeros(len(df)),
        'logg': df[logg_col].values if logg_col else np.zeros(len(df)),
        'meta': df[meta_col].values if meta_col else np.zeros(len(df)),
        'model_dir': model_path
    }
    
    return data

def load_sed(filepath):
    """Load a spectral energy distribution file."""
    wavelengths = []
    fluxes = []
    
    with open(filepath, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                try:
                    values = line.strip().split()
                    wavelengths.append(float(values[0]))
                    fluxes.append(float(values[1]))
                except:
                    pass
    
    return np.array(wavelengths), np.array(fluxes)

def create_unified_grid(all_models_data):
    """Create unified parameter grids from all models."""
    # Collect all parameter values
    all_teff = []
    all_logg = []
    all_meta = []
    
    for data in all_models_data:
        all_teff.extend(data['teff'])
        all_logg.extend(data['logg'])
        all_meta.extend(data['meta'])
    
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
    
    min_wave = float('inf')
    max_wave = 0
    resolutions = []
    
    for model_data in all_models_data:
        # Sample a few SEDs from this model
        n_sample = min(sample_size, len(model_data['files']))
        indices = np.random.choice(len(model_data['files']), n_sample, replace=False)
        
        for idx in indices:
            filepath = os.path.join(model_data['model_dir'], model_data['files'][idx])
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
            except:
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

def build_combined_flux_cube(all_models_data, teff_grid, logg_grid, meta_grid, wavelength_grid):
    """Build the combined flux cube from all models."""
    n_teff = len(teff_grid)
    n_logg = len(logg_grid)
    n_meta = len(meta_grid)
    n_lambda = len(wavelength_grid)
    
    # Initialize flux cube and tracking arrays
    flux_cube = np.zeros((n_teff, n_logg, n_meta, n_lambda))
    filled_map = np.zeros((n_teff, n_logg, n_meta), dtype=bool)
    source_map = np.zeros((n_teff, n_logg, n_meta), dtype=int) - 1  # Which model filled each point
    
    print(f"\nBuilding combined flux cube: {flux_cube.shape}")
    print(f"Memory requirement: {flux_cube.nbytes / (1024**2):.1f} MB")
    
    # Build KDTree for fast nearest neighbor searches
    # Normalize parameters for distance calculation
    teff_norm = (teff_grid - teff_grid.min()) / (teff_grid.max() - teff_grid.min())
    logg_norm = (logg_grid - logg_grid.min()) / (logg_grid.max() - logg_grid.min()) 
    meta_norm = (meta_grid - meta_grid.min()) / (meta_grid.max() - meta_grid.min())
    
    # Process each model
    for model_idx, model_data in enumerate(all_models_data):
        model_name = os.path.basename(model_data['model_dir'])
        print(f"\nProcessing {model_name}...")
        
        for i, (file, teff, logg, meta) in enumerate(tqdm(
            zip(model_data['files'], model_data['teff'], 
                model_data['logg'], model_data['meta']),
            total=len(model_data['files']), desc=f"Model {model_idx+1}")):
            
            # Find grid indices
            i_teff = np.searchsorted(teff_grid, teff)
            i_logg = np.searchsorted(logg_grid, logg)
            i_meta = np.searchsorted(meta_grid, meta)
            
            # Clip to valid range
            i_teff = np.clip(i_teff, 0, n_teff - 1)
            i_logg = np.clip(i_logg, 0, n_logg - 1)
            i_meta = np.clip(i_meta, 0, n_meta - 1)
            
            # Load and interpolate SED
            filepath = os.path.join(model_data['model_dir'], file)
            try:
                model_wavelengths, model_fluxes = load_sed(filepath)
                
                # Interpolate to common grid (in log space for flux)
                log_fluxes = np.log10(np.maximum(model_fluxes, 1e-50))
                log_interpolated = np.interp(
                    wavelength_grid, model_wavelengths, log_fluxes,
                    left=log_fluxes[0], right=log_fluxes[-1]
                )
                interpolated_flux = 10**log_interpolated
                
                # Store in cube
                flux_cube[i_teff, i_logg, i_meta, :] = interpolated_flux
                filled_map[i_teff, i_logg, i_meta] = True
                source_map[i_teff, i_logg, i_meta] = model_idx
                
            except Exception as e:
                continue
    
    # Fill gaps using nearest neighbor interpolation
    empty_points = np.sum(~filled_map)
    if empty_points > 0:
        print(f"\nFilling {empty_points} empty grid points...")
        
        # Get filled points
        filled_indices = np.where(filled_map)
        filled_points = np.column_stack([
            teff_norm[filled_indices[0]],
            logg_norm[filled_indices[1]], 
            meta_norm[filled_indices[2]]
        ])
        
        # Build KDTree
        tree = KDTree(filled_points)
        
        # Fill empty points
        for i_teff in range(n_teff):
            for i_logg in range(n_logg):
                for i_meta in range(n_meta):
                    if not filled_map[i_teff, i_logg, i_meta]:
                        # Find nearest filled point
                        query_point = [teff_norm[i_teff], logg_norm[i_logg], meta_norm[i_meta]]
                        dist, idx = tree.query(query_point, k=1)
                        
                        # Copy flux from nearest neighbor
                        src_i = filled_indices[0][idx]
                        src_j = filled_indices[1][idx]
                        src_k = filled_indices[2][idx]
                        
                        flux_cube[i_teff, i_logg, i_meta, :] = flux_cube[src_i, src_j, src_k, :]
                        filled_map[i_teff, i_logg, i_meta] = True
                        source_map[i_teff, i_logg, i_meta] = source_map[src_i, src_j, src_k]
    
    return flux_cube, source_map

def save_combined_data(output_dir, teff_grid, logg_grid, meta_grid, wavelength_grid, 
                       flux_cube, all_models_data):
    """Save the combined data and create unified lookup table."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Save binary flux cube
    binary_file = os.path.join(output_dir, 'flux_cube.bin')
    with open(binary_file, 'wb') as f:
        # Write dimensions
        f.write(struct.pack('4i', len(teff_grid), len(logg_grid), 
                           len(meta_grid), len(wavelength_grid)))
        
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
    
    for model_data in all_models_data:
        model_name = os.path.basename(model_data['model_dir'])
        for i, (orig_file, teff, logg, meta) in enumerate(
            zip(model_data['files'], model_data['teff'], 
                model_data['logg'], model_data['meta'])):
            
            # Create new filename that includes source model
            new_filename = f"{model_name}_{file_counter:06d}.txt"
            lookup_data.append({
                'file_name': new_filename,
                'teff': teff,
                'logg': logg,
                'meta': meta,
                'source_model': model_name,
                'original_file': orig_file
            })
            
            # Create symbolic link to original file
            src_path = os.path.join(model_data['model_dir'], orig_file)
            dst_path = os.path.join(output_dir, new_filename)
            if os.path.exists(src_path) and not os.path.exists(dst_path):
                try:
                    os.symlink(src_path, dst_path)
                except:
                    pass  # Symlinks might not work on all systems
            
            file_counter += 1
    
    # Save lookup table
    lookup_df = pd.DataFrame(lookup_data)
    lookup_file = os.path.join(output_dir, 'lookup_table.csv')
    
    with open(lookup_file, 'w') as f:
        f.write('#file_name, teff, logg, meta, source_model, original_file\n')
        lookup_df.to_csv(f, index=False, header=False)
    
    print(f"Saved combined lookup table to: {lookup_file}")
    print(f"Total models in combined set: {len(lookup_df)}")
    
    return lookup_df

def visualize_parameter_space(teff_grid, logg_grid, meta_grid, source_map, 
                            all_models_data, output_dir):
    """Create visualizations of the parameter space coverage."""
    print("\nCreating parameter space visualizations...")
    
    # Get model names
    model_names = [os.path.basename(data['model_dir']) for data in all_models_data]
    
    # Create color map for models
    colors = plt.cm.tab10(np.linspace(0, 1, len(model_names)))
    
    fig = plt.figure(figsize=(18, 12))
    
    # 1. 3D scatter plot of filled points
    ax1 = fig.add_subplot(2, 3, 1, projection='3d')
    
    # Get coordinates of filled points for each model
    for model_idx, model_name in enumerate(model_names):
        mask = source_map == model_idx
        if np.any(mask):
            indices = np.where(mask)
            ax1.scatter(teff_grid[indices[0]], logg_grid[indices[1]], 
                       meta_grid[indices[2]], c=[colors[model_idx]], 
                       label=model_name, alpha=0.6, s=20)
    
    ax1.set_xlabel('Teff (K)')
    ax1.set_ylabel('log g')
    ax1.set_zlabel('[M/H]')
    ax1.set_title('3D Parameter Space Coverage')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # 2. Teff vs log g projection
    ax2 = fig.add_subplot(2, 3, 2)
    for model_idx, model_name in enumerate(model_names):
        data = all_models_data[model_idx]
        ax2.scatter(data['teff'], data['logg'], c=[colors[model_idx]], 
                   label=model_name, alpha=0.6, s=10)
    ax2.set_xlabel('Teff (K)')
    ax2.set_ylabel('log g')
    ax2.set_title('Teff vs log g')
    ax2.invert_yaxis()
    ax2.grid(True, alpha=0.3)
    
    # 3. Teff vs [M/H] projection
    ax3 = fig.add_subplot(2, 3, 3)
    for model_idx, model_name in enumerate(model_names):
        data = all_models_data[model_idx]
        ax3.scatter(data['teff'], data['meta'], c=[colors[model_idx]], 
                   label=model_name, alpha=0.6, s=10)
    ax3.set_xlabel('Teff (K)')
    ax3.set_ylabel('[M/H]')
    ax3.set_title('Teff vs [M/H]')
    ax3.grid(True, alpha=0.3)
    
    # 4. log g vs [M/H] projection
    ax4 = fig.add_subplot(2, 3, 4)
    for model_idx, model_name in enumerate(model_names):
        data = all_models_data[model_idx]
        ax4.scatter(data['logg'], data['meta'], c=[colors[model_idx]], 
                   label=model_name, alpha=0.6, s=10)
    ax4.set_xlabel('log g')
    ax4.set_ylabel('[M/H]')
    ax4.set_title('log g vs [M/H]')
    ax4.grid(True, alpha=0.3)
    
    # 5. Parameter ranges bar chart
    ax5 = fig.add_subplot(2, 3, 5)
    param_ranges = []
    for data in all_models_data:
        param_ranges.append({
            'model': os.path.basename(data['model_dir']),
            'Teff_min': data['teff'].min(),
            'Teff_max': data['teff'].max(),
            'logg_min': data['logg'].min(),
            'logg_max': data['logg'].max(),
            'meta_min': data['meta'].min(),
            'meta_max': data['meta'].max(),
        })
    
    # Plot parameter ranges
    y_pos = np.arange(len(model_names))
    for i, (model_name, prange) in enumerate(zip(model_names, param_ranges)):
        # Normalize Teff to 0-1 range for plotting
        teff_norm_min = (prange['Teff_min'] - 2000) / 50000
        teff_norm_max = (prange['Teff_max'] - 2000) / 50000
        ax5.barh(i, teff_norm_max - teff_norm_min, left=teff_norm_min, 
                height=0.8, color=colors[i], alpha=0.7)
    
    ax5.set_yticks(y_pos)
    ax5.set_yticklabels(model_names)
    ax5.set_xlabel('Normalized Teff Range')
    ax5.set_title('Temperature Coverage by Model')
    ax5.set_xlim(0, 1)
    
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
    
    im = ax6.imshow(density_map.T, origin='lower', aspect='auto', 
                    extent=[teff_grid.min(), teff_grid.max(), 
                           logg_grid.min(), logg_grid.max()],
                    cmap='YlOrRd')
    ax6.set_xlabel('Teff (K)')
    ax6.set_ylabel('log g')
    ax6.set_title('Model Density (# of models per grid point)')
    plt.colorbar(im, ax=ax6)
    
    plt.tight_layout()
    plot_file = os.path.join(output_dir, 'parameter_space_visualization.png')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"Saved visualization to: {plot_file}")
    
    # Print summary statistics
    print("\n" + "="*60)
    print("COMBINED MODEL STATISTICS")
    print("="*60)
    print(f"Total number of source models: {len(all_models_data)}")
    print(f"Total grid points: {len(teff_grid)} × {len(logg_grid)} × {len(meta_grid)} = "
          f"{len(teff_grid) * len(logg_grid) * len(meta_grid):,}")
    print(f"\nParameter ranges:")
    print(f"  Teff: {teff_grid.min():.0f} - {teff_grid.max():.0f} K")
    print(f"  log g: {logg_grid.min():.2f} - {logg_grid.max():.2f}")
    print(f"  [M/H]: {meta_grid.min():.2f} - {meta_grid.max():.2f}")
    
    print(f"\nPer-model contributions:")
    for model_idx, model_name in enumerate(model_names):
        n_points = np.sum(source_map == model_idx)
        pct = 100 * n_points / source_map.size
        print(f"  {model_name}: {n_points:,} grid points ({pct:.1f}%)")

def main():
    parser = argparse.ArgumentParser(description='Combine multiple stellar atmosphere models')
    parser.add_argument('--data_dir', type=str, default='../../data/stellar_models/',
                       help='Base directory containing stellar models')
    parser.add_argument('--output', type=str, default='combined_models',
                       help='Output directory name')
    parser.add_argument('--interactive', action='store_true', default=True,
                       help='Interactive model selection')
    
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
    lookup_df = save_combined_data(
        output_dir, teff_grid, logg_grid, meta_grid, 
        wavelength_grid, flux_cube, all_models_data
    )
    
    # Create visualizations
    visualize_parameter_space(
        teff_grid, logg_grid, meta_grid, source_map, 
        all_models_data, output_dir
    )
    
    print(f"\nSuccessfully combined {len(selected_models)} stellar atmosphere models!")
    print(f"Output saved to: {output_dir}")
    print(f"You can now use this combined model in MESA by setting:")
    print(f"  stellar_atm = '{output_dir}/'")

if __name__ == "__main__":
    main()