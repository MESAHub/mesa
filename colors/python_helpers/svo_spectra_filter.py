import os

import pandas as pd


def filter_lookup_table(input_file):
    """
    Load a lookup table CSV, filter it based on user-specified ranges, and save the filtered table.
    """
    # Load the lookup table
    print(input_file)

    try:
        df = pd.read_csv("OG_" + input_file)
    except FileNotFoundError:
        df = pd.read_csv(input_file)

    # Ensure 'file_name' is preserved as a string
    df["#file_name"] = df["#file_name"].astype(str)

    # Clean up column headers and data
    df.columns = df.columns.str.strip()  # Strip spaces from column headers
    for column in df.columns:
        if column != "#file_name":  # Skip 'file_name' for numeric conversion
            try:
                df[column] = pd.to_numeric(
                    df[column], errors="coerce"
                )  # Convert columns to numeric if possible
            except Exception as e:
                print(f"Could not convert column {column} to numeric: {e}")

    print("\nAvailable parameters in the lookup table:")
    print([col for col in df.columns if col != "#file_name"])

    filters = {}
    for column in df.columns:
        if column == "#file_name":  # Skip filtering on file names
            continue

        # Gather column details
        unique_values = df[column].unique()
        unique_count = len(unique_values)
        print(f"\nColumn: {column}")

        if unique_count <= 5:
            print(f"Unique values: {unique_values}")
        else:
            column_range = (unique_values.min(), unique_values.max())
            median_val = df[column].median()
            print(f"Number of unique values: {unique_count}")
            print(
                f"Range: [{column_range[0]}, {column_range[1]}], Median: {median_val}"
            )

        # Prompt the user for filter input
        print(
            f"Set range or value for '{
                column
            }' (e.g., 1.0-10.0 for range, 5.0 for single value, or leave blank for no filter):"
        )
        user_input = input(f"Filter for {column}: ").strip()

        if user_input:
            try:
                if "-" in user_input:  # Check for range
                    min_val, max_val = map(float, user_input.split("-"))
                    filters[column] = ("range", min_val, max_val)
                else:  # Single value
                    value = float(user_input)
                    filters[column] = ("value", value)
            except ValueError:
                print(f"Invalid input for {column}. Skipping filter.")
                continue

    # Apply filters to the DataFrame
    for column, filter_info in filters.items():
        if filter_info[0] == "range":
            min_val, max_val = filter_info[1], filter_info[2]
            print(f"Filtering {column} between {min_val} and {max_val}...")
            df = df[(df[column] >= min_val) & (df[column] <= max_val)]
        elif filter_info[0] == "value":
            value = filter_info[1]
            print(f"Filtering {column} for exact value {value}...")
            df = df[df[column] == value]

        if df.empty:
            print("No rows match the specified filters. No changes made.")
            return

    # Rename the original table
    original_backup = input_file.replace("lookup_table.csv", "OG_lookup_table.csv")
    os.rename(input_file, original_backup)
    print(f"Original lookup table backed up as: {original_backup}")

    # Save the filtered DataFrame to the same file, keeping 'file_name' as a string
    df.to_csv(input_file, index=False)

    print(f"Filtered lookup table saved and replaced: {input_file}")


def list_lookup_tables(base_dir):
    """
    List all lookup tables in subdirectories under the base directory.
    """
    lookup_tables = []
    for root, _, files in os.walk(base_dir):
        for file in files:
            if file == "lookup_table.csv":
                lookup_tables.append(os.path.join(root, file))
    return lookup_tables


def clean_lookup_table(df):
    """
    Clean the lookup table by removing rows with missing or invalid data.
    """
    # Check for missing or incomplete rows
    if df.isnull().values.any():
        print(
            "\nWarning: The lookup table contains missing values. Cleaning the table..."
        )
        df = df.dropna(how="all").fillna(
            0
        )  # Drop fully empty rows, fill missing cells with 0
    return df


def main():
    print("Welcome to the Lookup Table Filter Tool!")

    # Base directory containing lookup tables
    base_dir = input(
        "Enter the path to the base directory containing lookup tables (Default : ../../data/stellar_models/): "
    ).strip()
    if not os.path.exists(base_dir):
        print("Error: Base directory does not exist.")
        base_dir = "../../data/stellar_models/"
        # return

    # List all lookup tables
    lookup_tables = list_lookup_tables(base_dir)
    if not lookup_tables:
        print("No lookup tables found in the specified directory.")
        return

    # Display available tables for selection
    print("\nAvailable lookup tables:")
    for idx, table in enumerate(lookup_tables, start=1):
        print(f"{idx}. {table}")

    # Allow user to select a table
    user_input = input(
        "\nEnter the number of the lookup table you want to filter: "
    ).strip()
    try:
        selected_idx = int(user_input) - 1
        if 0 <= selected_idx < len(lookup_tables):
            selected_table = lookup_tables[selected_idx]
            print(f"\nSelected lookup table: {selected_table}")

            # Filter the selected table
            filter_lookup_table(selected_table)
        else:
            print("Invalid selection. Exiting.")
    except ValueError:
        print("Invalid input. Please enter a number. Exiting.")


if __name__ == "__main__":
    main()
