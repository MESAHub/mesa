#!/usr/bin/env python3

import glob
import textwrap

import matplotlib.pyplot as plt
import mesa_reader as mr
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from static_HISTORY_check import (MesaView, read_header_columns,
                                  setup_hr_diagram_params)


def get_z_axis_selection(available_columns, default="Interp_rad"):
    """
    Prints available columns and prompts the user to select one for the Z-axis.
    """
    print("\n" + "=" * 60)
    print("AVAILABLE COLUMNS FROM HISTORY DATA")
    print("=" * 60)

    # Sort and format the list of columns to be readable
    sorted_cols = sorted(available_columns)
    col_list_str = ", ".join(sorted_cols)

    # Wrap text so it doesn't overflow the terminal
    print(textwrap.fill(col_list_str, width=80))
    print("-" * 60)

    while True:
        user_input = input(
            f"\nEnter the Z-axis column name (default: {default}): "
        ).strip()

        # Handle Default (Enter key)
        if user_input == "":
            if default in available_columns:
                print(f"Selected default: {default}")
                return default
            else:
                print(
                    f"Warning: Default '{default}' not found in file. Please select manually."
                )
                continue

        # Handle Valid Selection
        if user_input in available_columns:
            print(f"Selected: {user_input}")
            return user_input

        # Handle Invalid Selection
        print(f"Error: '{user_input}' is not a valid column. Please try again.")


def make_3d_cmd(history_file="../LOGS/history.data"):
    md = mr.MesaData(history_file)
    md = MesaView(md, 5)

    all_cols, filter_columns = read_header_columns(history_file)
    hr_color, hr_mag, hr_xlabel, hr_ylabel, color_index = setup_hr_diagram_params(
        md, filter_columns
    )

    # --- NEW INTERACTIVE SECTION ---
    # Prompt user for Z-axis, using Interp_rad as default
    z_col_name = get_z_axis_selection(all_cols, default="Interp_rad")

    # Dynamically retrieve the data column from the MesaData object
    # getattr(md, "name") is equivalent to md.name
    try:
        z_data = getattr(md, z_col_name)
    except AttributeError:
        print(f"Error: Could not retrieve data for {z_col_name}.")
        return
    # -------------------------------

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Use the dynamically selected z_data
    ax.scatter(hr_color, hr_mag, z_data, s=5)

    ax.set_xlabel(hr_xlabel)
    ax.set_ylabel(hr_ylabel)
    ax.set_zlabel(z_col_name)  # Update label to match selection

    ax.invert_yaxis()

    plt.tight_layout()
    plt.show()


def main():
    # Helper to find file if specific path not strictly enforced
    files = glob.glob("../LOGS/history.data")
    if not files:
        print("Error: ../LOGS/history.data not found.")
        return

    history_file = files[0]
    make_3d_cmd(history_file)


if __name__ == "__main__":
    main()
