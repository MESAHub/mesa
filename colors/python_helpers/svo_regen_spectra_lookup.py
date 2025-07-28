import csv
import os
import re


def clean_metadata_values(metadata):
    """
    Clean metadata values to retain only numerical parts, except for the file_name key.
    Replace missing values with 999.9.
    """
    cleaned_metadata = {}
    for key, value in metadata.items():
        if key == "file_name":
            cleaned_metadata[key] = value
        else:
            # Matches floats or integers
            match = re.search(r"[-+]?\d*\.\d+|\d+", value)
            cleaned_metadata[key] = (
                match.group() if match else "999.9"
            )  # Default to 999.9 if no match
    return cleaned_metadata


def parse_metadata(file_path):
    """
    Parse metadata from a file where metadata lines start with '#' and contain '='.
    """
    metadata = {}
    try:
        with open(file_path, "r") as file:
            for line in file:
                line = line.strip()
                if line.startswith("#") and "=" in line:
                    try:
                        key, value = line.split("=", 1)
                        key = key.strip("#").strip()
                        value = value.split("(")[0].strip()
                        metadata[key] = value
                    except ValueError:
                        continue
    except Exception as e:
        print(f"Error parsing metadata in {file_path}: {e}")
    return metadata


def scan_existing_files(base_dir):
    """
    Scan the directory structure to find existing models and spectra files.

    Returns:
        dict: A mapping of model names to their respective spectra files.
    """
    existing_files = {}
    for root, _, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".txt"):
                model_name = os.path.basename(root)
                if model_name not in existing_files:
                    existing_files[model_name] = []
                existing_files[model_name].append(os.path.join(root, file))
    return existing_files


def regenerate_lookup_table(model, files, output_dir):
    """
    Regenerate the lookup table from existing files.

    Args:
        model (str): The model name.
        files (list): List of file paths for the model.
        output_dir (str): Directory to save the lookup table.
    """
    lookup_table_path = os.path.join(output_dir, "lookup_table.csv")
    all_keys = set()
    metadata_rows = []

    for file_path in files:
        metadata = parse_metadata(file_path)
        metadata["file_name"] = os.path.basename(file_path)
        metadata = clean_metadata_values(metadata)
        metadata_rows.append(metadata)
        all_keys.update(metadata.keys())

    if metadata_rows:
        with open(lookup_table_path, mode="w", newline="") as csv_file:
            # Generate a header dynamically with #
            header = ["file_name"] + sorted(all_keys - {"file_name"})
            csv_file.write(
                "#" + ", ".join(header) + "\n"
            )  # Write the header prefixed with #
            writer = csv.DictWriter(csv_file, fieldnames=header)
            writer.writerows(metadata_rows)

        print(f"Lookup table regenerated: {lookup_table_path}")


def main():
    # Define the base directory
    base_dir = "../../data/stellar_models/"
    os.makedirs(base_dir, exist_ok=True)

    # Scan for existing models and files
    existing_files = scan_existing_files(base_dir)

    # Regenerate lookup tables for each model
    for model, files in existing_files.items():
        print(f"Processing model: {model}")
        output_dir = os.path.join(base_dir, model)
        regenerate_lookup_table(model, files, output_dir)


if __name__ == "__main__":
    main()
