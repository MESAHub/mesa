#!/usr/bin/env python3
"""
compare_mag_systems.py

Runs MESA 3 times with different magnitude zero-point systems (Vega, AB, ST)
and overlays the resulting CMDs for comparison.

Usage: python compare_mag_systems.py
       (Run from the custom_colors test suite directory)

Niall Miller (2025)
"""

import os
import re
import shutil
import subprocess

import matplotlib.pyplot as plt
import numpy as np

# Configuration
MAG_SYSTEMS = ["Vega", "AB", "ST"]
COLORS = {"Vega": "blue", "AB": "red", "ST": "green"}
INLIST_COLORS = "inlist_colors"
INLIST_HEADER = "inlist_colors_header"

os.chdir("../")

# Fast run parameters - override these in the inlist for quick comparison
FAST_RUN_OVERRIDES = {
    "initial_mass": "3.0d0",  # Lower mass = faster evolution
    "xa_central_lower_limit(1)": "0.5",  # Stop at 50% He depletion (faster)
    "varcontrol_target": "1d-2",  # Looser tolerance
    "mesh_delta_coeff": "1.0",  # Coarser mesh
}


def read_inlist(filename):
    """Read inlist file content."""
    with open(filename, "r") as f:
        return f.read()


def write_inlist(filename, content):
    """Write inlist file content."""
    with open(filename, "w") as f:
        f.write(content)


def modify_mag_system(content, mag_system):
    """Change the mag_system parameter in inlist content."""
    pattern = r"(mag_system\s*=\s*')[^']*(')"
    return re.sub(pattern, rf"\g<1>{mag_system}\g<2>", content)


def modify_for_fast_run(content):
    """Modify inlist for faster execution."""
    # Disable pgstar for batch runs
    content = re.sub(r"pgstar_flag\s*=\s*\.true\.", "pgstar_flag = .false.", content)

    # Apply fast run overrides
    for param, value in FAST_RUN_OVERRIDES.items():
        pattern = rf"({param}\s*=\s*)[^\s!]+"
        if re.search(pattern, content):
            content = re.sub(pattern, rf"\g<1>{value}", content)

    return content


def run_mesa(mag_system, run_dir):
    """Run MESA with specified mag_system, saving LOGS to run_dir."""
    print(f"\n{'='*60}")
    print(f"Running MESA with mag_system = {mag_system}")
    print(f"{'='*60}")

    # Backup original inlist
    original_content = read_inlist(INLIST_COLORS)

    # Modify inlist
    modified_content = modify_mag_system(original_content, mag_system)
    modified_content = modify_for_fast_run(modified_content)
    write_inlist(INLIST_COLORS, modified_content)

    # Clean LOGS directory
    if os.path.exists("LOGS"):
        shutil.rmtree("LOGS")
    os.makedirs("LOGS", exist_ok=True)

    # Run MESA
    try:
        result = subprocess.run(
            ["./rn"], capture_output=True, text=True, timeout=3600  # 1 hour timeout
        )
        if result.returncode != 0:
            print(f"MESA run failed for {mag_system}")
            print(result.stderr)
    except subprocess.TimeoutExpired:
        print(f"MESA run timed out for {mag_system}")
    except Exception as e:
        print(f"Error running MESA: {e}")
    finally:
        # Restore original inlist
        write_inlist(INLIST_COLORS, original_content)

    # Copy LOGS to run-specific directory
    if os.path.exists("LOGS"):
        if os.path.exists(run_dir):
            shutil.rmtree(run_dir)
        shutil.copytree("LOGS", run_dir)
        print(f"Results saved to {run_dir}")


def read_mesa_history(logs_dir):
    """Read MESA history.data file and return as dict of arrays."""
    history_file = os.path.join(logs_dir, "history.data")

    if not os.path.exists(history_file):
        print(f"Warning: {history_file} not found")
        return None

    # Read header to get column names
    with open(history_file, "r") as f:
        lines = f.readlines()

    # Find header line (line with column names, typically line 6)
    header_line = 5  # 0-indexed, so line 6
    for i, line in enumerate(lines):
        if line.strip().startswith("model_number") or "model_number" in line:
            header_line = i
            break

    column_names = lines[header_line].split()

    # Read data (skip header lines)
    data = np.loadtxt(history_file, skiprows=header_line + 1)

    # Create dictionary
    result = {}
    for i, name in enumerate(column_names):
        if i < data.shape[1]:
            result[name] = data[:, i]

    return result


def plot_cmd_comparison(results):
    """Plot CMD comparison for all mag systems."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # CMD 1: B-V vs V
    ax1 = axes[0]
    for mag_system, data in results.items():
        if data is None:
            continue
        if "B" in data and "V" in data:
            color = data["B"] - data["V"]
            mag = data["V"]
            ax1.plot(
                color,
                mag,
                "-",
                color=COLORS[mag_system],
                label=mag_system,
                alpha=0.8,
                linewidth=1.5,
            )

    ax1.set_xlabel("B - V", fontsize=12)
    ax1.set_ylabel("V (mag)", fontsize=12)
    ax1.set_title("CMD: B-V vs V", fontsize=14)
    ax1.invert_yaxis()
    ax1.legend(loc="best")
    ax1.grid(True, alpha=0.3)

    # CMD 2: U-B vs B
    ax2 = axes[1]
    for mag_system, data in results.items():
        if data is None:
            continue
        if "U" in data and "B" in data:
            color = data["U"] - data["B"]
            mag = data["B"]
            ax2.plot(
                color,
                mag,
                "-",
                color=COLORS[mag_system],
                label=mag_system,
                alpha=0.8,
                linewidth=1.5,
            )

    ax2.set_xlabel("U - B", fontsize=12)
    ax2.set_ylabel("B (mag)", fontsize=12)
    ax2.set_title("CMD: U-B vs B", fontsize=14)
    ax2.invert_yaxis()
    ax2.legend(loc="best")
    ax2.grid(True, alpha=0.3)

    plt.suptitle("Magnitude System Comparison", fontsize=16, y=1.02)
    plt.tight_layout()

    output_file = "mag_system_comparison.png"
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    print(f"\nPlot saved to {output_file}")
    #plt.show()


def main():
    """Main execution."""
    # Check we're in the right directory
    if not os.path.exists(INLIST_COLORS):
        print(f"Error: {INLIST_COLORS} not found. Run from custom_colors directory.")
        return

    if not os.path.exists("./rn"):
        print("Error: ./rn not found. Run from custom_colors directory.")
        return

    results = {}

    # Run MESA for each mag system
    for mag_system in MAG_SYSTEMS:
        run_dir = f"LOGS_{mag_system}"
        run_mesa(mag_system, run_dir)
        results[mag_system] = read_mesa_history(run_dir)

    # Plot comparison
    plot_cmd_comparison(results)

    print("\n" + "=" * 60)
    print("Comparison complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
