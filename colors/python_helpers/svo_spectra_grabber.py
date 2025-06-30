#!/usr/bin/env python3
"""
Improved SVO Stellar Spectra Grabber
Uses proper API queries to discover and download all available spectra
"""

import csv
import json
import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests
from bs4 import BeautifulSoup
from tqdm import tqdm


class SVOSpectraGrabber:
    def __init__(self, base_dir="../data/stellar_models/", max_workers=5):
        self.base_dir = base_dir
        self.max_workers = max_workers
        self.session = requests.Session()
        self.session.headers.update(
            {"User-Agent": "Mozilla/5.0 (MESA Colors Module)"})

        # SVO endpoints
        self.model_index_url = "http://svo2.cab.inta-csic.es/theory/newov2/index.php"
        self.spectra_base_url = "http://svo2.cab.inta-csic.es/theory/newov2/ssap.php"
        self.metadata_url = "http://svo2.cab.inta-csic.es/theory/newov2/getmeta.php"

        os.makedirs(base_dir, exist_ok=True)

    def discover_models(self):
        """Discover all available stellar atmosphere models."""
        print("Discovering available models from SVO...")

        try:
            response = self.session.get(self.model_index_url, timeout=30)
            response.raise_for_status()
        except requests.RequestException as e:
            print(f"Error fetching model index: {e}")
            return []

        soup = BeautifulSoup(response.text, features="xml")
        models = []

        # Look for model links
        for link in soup.find_all("a", href=True):
            href = link["href"]
            if "models=" in href:
                model_name = href.split("models=")[1].split("&")[0]
                if model_name not in models:
                    models.append(model_name)

        return sorted(models)

    def get_model_metadata(self, model_name):
        """Get metadata about available spectra for a specific model."""
        print(f"  Fetching metadata for {model_name}...")

        # Try different approaches to get spectrum list
        spectra_info = []

        # Method 1: Direct metadata query
        try:
            params = {"model": model_name, "format": "json"}
            response = self.session.get(
                self.metadata_url, params=params, timeout=30)
            if response.status_code == 200:
                try:
                    metadata = response.json()
                    if isinstance(metadata, list):
                        spectra_info.extend(metadata)
                except json.JSONDecodeError:
                    pass
        except requests.RequestException:
            pass

        # Method 2: SSAP query for spectrum discovery
        if not spectra_info:
            spectra_info = self._discover_spectra_ssap(model_name)

        # Method 3: Brute force discovery (fallback)
        if not spectra_info:
            print(f"    Using brute force discovery for {model_name}")
            spectra_info = self._discover_spectra_brute_force(model_name)

        print(f"    Found {len(spectra_info)} spectra for {model_name}")
        return spectra_info

    def _discover_spectra_ssap(self, model_name):
        """Use SSAP protocol to discover available spectra."""
        spectra_info = []

        try:
            # Query for all spectra in this model
            params = {"model": model_name,
                      "REQUEST": "queryData", "FORMAT": "metadata"}

            response = self.session.get(
                self.spectra_base_url, params=params, timeout=30
            )
            if response.status_code == 200:
                # Parse response for spectrum IDs
                soup = BeautifulSoup(response.text, features="xml")
                # Look for spectrum identifiers in the response
                for link in soup.find_all("a", href=True):
                    href = link["href"]
                    if "fid=" in href:
                        fid = href.split("fid=")[1].split("&")[0]
                        if fid.isdigit():
                            spectra_info.append({"fid": int(fid)})

        except requests.RequestException:
            pass

        return spectra_info

    def _discover_spectra_brute_force(self, model_name, max_fid=20000, batch_size=50):
        """Brute force discovery with intelligent searching."""
        print("    Scanning for spectra (this may take a while)...")

        spectra_info = []
        consecutive_failures = 0
        max_consecutive_failures = 100  # More generous than before

        # Use batch requests to speed up discovery
        for start_fid in tqdm(range(1, max_fid, batch_size), desc="Scanning FIDs"):
            batch_fids = list(
                range(start_fid, min(start_fid + batch_size, max_fid)))

            # Test batch in parallel
            with ThreadPoolExecutor(max_workers=min(10, batch_size)) as executor:
                future_to_fid = {
                    executor.submit(self._test_spectrum_exists, model_name, fid): fid
                    for fid in batch_fids
                }

                for future in as_completed(future_to_fid):
                    fid = future_to_fid[future]
                    try:
                        exists = future.result()
                        if exists:
                            spectra_info.append({"fid": fid})
                            consecutive_failures = 0
                        else:
                            consecutive_failures += 1
                    except Exception:
                        consecutive_failures += 1

            # Early termination if we've gone too long without finding anything
            if (
                consecutive_failures > max_consecutive_failures
                and len(spectra_info) > 0
            ):
                print(
                    f"    Stopping search after {
                        consecutive_failures} consecutive failures"
                )
                break

        return spectra_info

    def _test_spectrum_exists(self, model_name, fid):
        """Test if a spectrum exists without downloading it."""
        try:
            params = {"model": model_name, "fid": fid, "format": "ascii"}
            response = self.session.head(
                self.spectra_base_url, params=params, timeout=10
            )
            return (
                response.status_code == 200
                and response.headers.get("content-length", "0") != "0"
            )
        except (requests.RequestException, ValueError, KeyError):
            return False

    def download_spectrum(self, model_name, fid, output_path):
        """Download a single spectrum."""
        try:
            params = {"model": model_name, "fid": fid, "format": "ascii"}
            response = self.session.get(
                self.spectra_base_url, params=params, timeout=30, stream=True
            )

            if response.status_code == 200 and len(response.content) > 1024:
                with open(output_path, "wb") as file:
                    file.write(response.content)
                return True
            else:
                return False

        except Exception as e:
            print(f"    Error downloading FID {fid}: {e}")
            return False

    def parse_metadata(self, file_path):
        """Parse metadata from a downloaded spectrum file."""
        metadata = {}
        try:
            with open(file_path, "r", encoding="utf-8", errors="ignore") as file:
                for line in file:
                    line = line.strip()
                    if line.startswith("#") and "=" in line:
                        try:
                            key, value = line.split("=", 1)
                            key = key.strip("#").strip()
                            value = value.split("(")[0].strip()

                            # Clean numerical values
                            if key != "file_name":
                                match = re.search(r"[-+]?\d*\.?\d+", value)
                                value = match.group() if match else "999.9"

                            metadata[key] = value
                        except ValueError:
                            continue
                    elif not line.startswith("#"):
                        break  # Stop at first data line
        except Exception as e:
            print(f"    Error parsing metadata in {file_path}: {e}")

        return metadata

    def download_model_spectra(self, model_name, spectra_info):
        """Download all spectra for a given model."""
        output_dir = os.path.join(self.base_dir, model_name)
        os.makedirs(output_dir, exist_ok=True)

        print(f"Downloading {len(spectra_info)} spectra for {model_name}...")

        metadata_rows = []
        successful_downloads = 0

        # Download spectra in parallel
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_fid = {}

            for spectrum in spectra_info:
                fid = spectrum["fid"]
                filename = f"{model_name}_fid{fid}.txt"
                output_path = os.path.join(output_dir, filename)

                # Skip if already exists
                if os.path.exists(output_path):
                    metadata = self.parse_metadata(output_path)
                    metadata["file_name"] = filename
                    metadata_rows.append(metadata)
                    successful_downloads += 1
                    continue

                future = executor.submit(
                    self.download_spectrum, model_name, fid, output_path
                )
                future_to_fid[future] = (fid, filename, output_path)

            # Process completed downloads
            for future in tqdm(
                as_completed(future_to_fid),
                total=len(future_to_fid),
                desc=f"Downloading {model_name}",
            ):
                fid, filename, output_path = future_to_fid[future]

                try:
                    success = future.result()
                    if success:
                        metadata = self.parse_metadata(output_path)
                        metadata["file_name"] = filename
                        metadata_rows.append(metadata)
                        successful_downloads += 1
                    else:
                        # Clean up failed download
                        if os.path.exists(output_path):
                            os.remove(output_path)
                except Exception as e:
                    print(f"    Error processing FID {fid}: {e}")

        # Create lookup table
        if metadata_rows:
            self._create_lookup_table(output_dir, metadata_rows)
            print(f"  Successfully downloaded {successful_downloads} spectra")
        else:
            print(
                f"  No spectra downloaded for {
                    model_name}. Consider using --force-brute if the model lacks metadata."
            )

        return successful_downloads

    def _create_lookup_table(self, output_dir, metadata_rows):
        """Create lookup table CSV for the model."""
        lookup_table_path = os.path.join(output_dir, "lookup_table.csv")

        # Get all unique metadata keys
        all_keys = set()
        for row in metadata_rows:
            all_keys.update(row.keys())

        # Define column order
        header = ["file_name"] + sorted(all_keys - {"file_name"})

        # Write CSV
        with open(
            lookup_table_path, mode="w", newline="", encoding="utf-8"
        ) as csv_file:
            csv_file.write("#" + ", ".join(header) + "\n")
            writer = csv.DictWriter(
                csv_file, fieldnames=header, extrasaction="ignore")
            writer.writerows(metadata_rows)

        print(f"    Lookup table saved: {lookup_table_path}")

    def select_models_interactive(self, available_models):
        """Allow user to select which models to download."""
        print("\nAvailable stellar atmosphere models:")
        print("-" * 50)
        for idx, model in enumerate(available_models, start=1):
            print(f"{idx:2d}. {model}")

        print("\nEnter model numbers to download (comma-separated):")
        print("Example: 1,3,5 or 'all' for all models")

        user_input = input("> ").strip()

        if user_input.lower() == "all":
            return available_models

        try:
            indices = [int(x.strip()) - 1 for x in user_input.split(",")]
            selected = [
                available_models[i] for i in indices if 0 <= i < len(available_models)
            ]
            return selected
        except (ValueError, IndexError):
            print("Invalid input. Please try again.")
            return self.select_models_interactive(available_models)

    def run(self, selected_models=None):
        """Main execution function."""
        print("SVO Stellar Atmosphere Model Downloader")
        print("=" * 50)

        # Discover available models
        available_models = self.discover_models()
        if not available_models:
            print("No models found on SVO!")
            return

        # Select models to download
        if selected_models is None:
            selected_models = self.select_models_interactive(available_models)

        if not selected_models:
            print("No models selected.")
            return

        print(f"\nSelected models: {', '.join(selected_models)}")
        print("=" * 50)

        total_spectra = 0

        # Process each model
        for model_name in selected_models:
            print(f"\nProcessing model: {model_name}")
            print("-" * 30)

            # Get metadata about available spectra
            spectra_info = self.get_model_metadata(model_name)

            if not spectra_info:
                print(f"  No spectra found for {model_name}")
                continue

            # Download spectra
            downloaded = self.download_model_spectra(model_name, spectra_info)
            total_spectra += downloaded

        print("\n" + "=" * 50)
        print("Download complete!")
        print(f"Total spectra downloaded: {total_spectra}")
        print(f"Models processed: {len(selected_models)}")
        print(f"Output directory: {self.base_dir}")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Download stellar spectra from SVO")
    parser.add_argument(
        "--output",
        type=str,
        default="../data/stellar_models/",
        help="Output directory for downloaded spectra",
    )
    parser.add_argument(
        "--workers", type=int, default=5, help="Number of parallel download workers"
    )
    parser.add_argument(
        "--models",
        type=str,
        nargs="*",
        help="Specific models to download (if not provided, interactive selection)",
    )

    args = parser.parse_args()

    # Create downloader
    downloader = SVOSpectraGrabber(
        base_dir=args.output, max_workers=args.workers)

    # Run downloader
    downloader.run(selected_models=args.models)


if __name__ == "__main__":
    main()
