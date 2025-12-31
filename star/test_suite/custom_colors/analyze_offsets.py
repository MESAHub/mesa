#!/usr/bin/env python3
"""
analyze_offsets.py

Analyzes the offsets between Vega, AB, and ST magnitude systems
from MESA Colors runs and compares to expected theoretical values.

Expected offsets (AB - Vega) and (ST - Vega) for Johnson filters:
These depend on Vega's SED convolved with filter transmission curves.

Niall Miller (2025)
"""

import numpy as np
import os

# Known approximate offsets for Johnson filters (from literature)
# AB - Vega offsets (Bessell & Murphy 2012, Blanton & Roweis 2007)
EXPECTED_AB_VEGA = {
    'U': 0.79,
    'B': -0.09,
    'V': 0.02,
    'R': 0.21,
    'I': 0.45,
}

# ST - Vega offsets (calculated from ST zeropoint definition)
# ST uses f_λ = 3.631e-9 erg/s/cm²/Å at all wavelengths
# These vary with filter effective wavelength
EXPECTED_ST_VEGA = {
    'U': -0.90,
    'B': -0.64,
    'V': -0.01,
    'R': 0.43,
    'I': 0.89,
}


def read_mesa_history(logs_dir):
    """Read MESA history.data file."""
    history_file = os.path.join(logs_dir, 'history.data')
    
    if not os.path.exists(history_file):
        print(f"Warning: {history_file} not found")
        return None
    
    with open(history_file, 'r') as f:
        lines = f.readlines()
    
    # Find header line
    header_line = 5
    for i, line in enumerate(lines):
        if 'model_number' in line:
            header_line = i
            break
    
    column_names = lines[header_line].split()
    data = np.loadtxt(history_file, skiprows=header_line + 1)
    
    result = {}
    for i, name in enumerate(column_names):
        if i < data.shape[1]:
            result[name] = data[:, i]
    
    return result


def analyze_offsets():
    """Compute and compare magnitude offsets."""
    
    # Load data from each run
    vega_data = read_mesa_history('LOGS_Vega')
    ab_data = read_mesa_history('LOGS_AB')
    st_data = read_mesa_history('LOGS_ST')
    
    if vega_data is None or ab_data is None or st_data is None:
        print("Error: Could not load all history files")
        return
    
    filters = ['U', 'B', 'V', 'R', 'I']
    available_filters = [f for f in filters if f in vega_data]
    
    print("=" * 70)
    print("MAGNITUDE SYSTEM OFFSET ANALYSIS")
    print("=" * 70)
    
    # Use matching model numbers or same indices
    # For simplicity, use median offset (should be constant throughout evolution)
    
    print(f"\n{'Filter':<8} {'AB-Vega (measured)':<20} {'AB-Vega (expected)':<20} {'Δ':<10}")
    print("-" * 58)
    
    ab_vega_offsets = {}
    for filt in available_filters:
        if filt in ab_data and filt in vega_data:
            # Use median to avoid edge effects
            offset = np.median(ab_data[filt] - vega_data[filt])
            ab_vega_offsets[filt] = offset
            expected = EXPECTED_AB_VEGA.get(filt, np.nan)
            diff = offset - expected if not np.isnan(expected) else np.nan
            print(f"{filt:<8} {offset:<20.3f} {expected:<20.3f} {diff:<+10.3f}")
    
    print(f"\n{'Filter':<8} {'ST-Vega (measured)':<20} {'ST-Vega (expected)':<20} {'Δ':<10}")
    print("-" * 58)
    
    st_vega_offsets = {}
    for filt in available_filters:
        if filt in st_data and filt in vega_data:
            offset = np.median(st_data[filt] - vega_data[filt])
            st_vega_offsets[filt] = offset
            expected = EXPECTED_ST_VEGA.get(filt, np.nan)
            diff = offset - expected if not np.isnan(expected) else np.nan
            print(f"{filt:<8} {offset:<20.3f} {expected:<20.3f} {diff:<+10.3f}")
    
    # Also check AB - ST (should equal AB-Vega minus ST-Vega)
    print(f"\n{'Filter':<8} {'AB-ST (measured)':<20} {'AB-ST (expected)':<20} {'Δ':<10}")
    print("-" * 58)
    
    for filt in available_filters:
        if filt in ab_data and filt in st_data:
            offset = np.median(ab_data[filt] - st_data[filt])
            expected_ab = EXPECTED_AB_VEGA.get(filt, 0)
            expected_st = EXPECTED_ST_VEGA.get(filt, 0)
            expected = expected_ab - expected_st
            diff = offset - expected
            print(f"{filt:<8} {offset:<20.3f} {expected:<20.3f} {diff:<+10.3f}")
    
    print("\n" + "=" * 70)
    print("VERIFICATION SUMMARY")
    print("=" * 70)
    
    # Check if offsets are constant throughout evolution (they should be)
    print("\nOffset stability (std dev across evolution - should be ~0):")
    for filt in available_filters:
        if filt in ab_data and filt in vega_data:
            std = np.std(ab_data[filt] - vega_data[filt])
            print(f"  {filt} AB-Vega: σ = {std:.4f} mag")
    
    print("\nNote: Expected values from Bessell & Murphy (2012) and standard")
    print("ST/AB definitions. Small differences may arise from different Vega")
    print("spectra or filter transmission curves used.")


if __name__ == '__main__':
    analyze_offsets()
