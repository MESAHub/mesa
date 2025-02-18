import requests
from bs4 import BeautifulSoup
import os
import csv
import time
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
            # Use regex to extract the numerical part (e.g., remove units like 'log' or 'K')
            match = re.search(r"[-+]?\d*\.\d+|\d+", value)  # Matches floats or integers
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


def generate_lookup_header(file_path):
    """
    Generate a lookup table header based on metadata keys.

    Args:
        file_path (str): Path to a file with metadata to parse.

    Returns:
        str: Header string prefixed by a hash (e.g., '#filename, log(g), teff, me/h').
    """

    keys = []
    try:
        with open(file_path, "r") as file:
            for line in file:
                line = line.strip()
                if line.startswith("#") and "=" in line:
                    try:
                        key, value = line.split("=", 1)
                        key = key.strip("#").strip()
                        keys.append(key)
                    except ValueError:
                        continue
    except Exception as e:
        print(f"Error parsing metadata in {file_path}: {e}")
    header = "#filename, " + ", ".join(keys)
    return header


def fetch_model_links(base_url):
    """
    Fetch model links from the base URL and allow the user to select specific models.
    """
    response = requests.get(base_url)
    if response.status_code != 200:
        print("Failed to fetch model index page.")
        return []

    soup = BeautifulSoup(response.text, "html.parser")
    model_links = [
        a["href"].split("=")[1]
        for a in soup.find_all("a", href=True)
        if "models=" in a["href"]
    ]

    if not model_links:
        print("No models found.")
        return []

    # Display models for user selection
    print("Available models:")
    for idx, model in enumerate(model_links, start=1):
        print(f"{idx}. {model}")

    # Ask the user to select models
    user_input = input(
        "Enter the numbers of the models you want to select, separated by commas: "
    )

    try:
        selected_indices = [int(num.strip()) - 1 for num in user_input.split(",")]
        selected_models = [
            model_links[i] for i in selected_indices if 0 <= i < len(model_links)
        ]
        return selected_models
    except (ValueError, IndexError):
        print("Invalid input. Please enter valid numbers separated by commas.")
        return []


def download_spectrum(session, spectra_url, params, output_fp):
    """
    Download a spectrum file and save it to the specified path.
    """
    try:
        response = session.get(spectra_url, params=params, stream=True)
        if (
            response.status_code == 200 and len(response.content) > 1024
        ):  # Ensure valid content
            with open(output_fp, "wb") as file:
                file.write(response.content)
            return True
    except Exception as e:
        print(f"Error downloading spectrum: {e}")
    return False


def main():
    # Base URLs
    model_index_url = "http://svo2.cab.inta-csic.es/theory/newov2/index.php"
    spectra_base_url = "http://svo2.cab.inta-csic.es/theory/newov2/ssap.php"

    base_dir = "../data/stellar_models/"
    os.makedirs(base_dir, exist_ok=True)

    # Fetch model names
    models = fetch_model_links(model_index_url)
    if not models:
        print("No models found.")
        return
    print(f"Found models: {models}")

    session = requests.Session()  # Reuse session for efficiency

    for model in models:
        fid = 0
        output_dir = os.path.join(base_dir, model)
        os.makedirs(output_dir, exist_ok=True)
        print(f"Processing model: {model}")

        lookup_table_path = os.path.join(output_dir, "lookup_table.csv")
        all_keys = set()
        metadata_rows = []
        found_spectra = 0
        while True:
            filename = f"{model}_fid{fid}.txt"
            output_fp = os.path.join(output_dir, filename)

            if os.path.exists(output_fp):
                print(f"File already exists: {output_fp}")
                metadata = parse_metadata(output_fp)
                metadata["file_name"] = filename
                metadata = clean_metadata_values(metadata)  # Clean metadata values
                metadata_rows.append(metadata)
                all_keys.update(metadata.keys())
                fid += 1
                continue

            fid += 1
            # Query and download spectrum
            params = {"model": model, "fid": fid, "format": "ascii"}
            if download_spectrum(session, spectra_base_url, params, output_fp):
                print(f"Downloaded: {filename}")
                metadata = parse_metadata(output_fp)
                metadata["file_name"] = filename
                metadata = clean_metadata_values(metadata)  # Clean metadata values
                metadata_rows.append(metadata)
                all_keys.update(metadata.keys())
                found_spectra = 1
            else:
                print(f"No more spectra for model {model}. Last fid: {fid - 1}")
                if found_spectra > 0:

                    found_spectra = found_spectra + 1
                    if found_spectra == 20:
                        break

        # Write lookup table to CSV
        if metadata_rows:
            with open(lookup_table_path, mode="w", newline="") as csv_file:
                # Generate a header dynamically with #
                header = ["file_name"] + sorted(all_keys - {"file_name"})
                csv_file.write(
                    "#" + ", ".join(header) + "\n"
                )  # Write the header prefixed with #
                writer = csv.DictWriter(csv_file, fieldnames=header)
                writer.writerows(metadata_rows)

            print(f"Lookup table saved: {lookup_table_path}")


if __name__ == "__main__":
    main()
