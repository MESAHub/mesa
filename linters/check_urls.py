#!/usr/bin/env python

import os
import sys
import re
import requests

MESA_DIR = os.environ.get("MESA_DIR", "../")

# Look for all urls in the *.rst files in $MESA_DIR/docs/source/ directory and subdirectories,
# and check if they are valid by making a HEAD request to each URL.
# If the URL is not valid, print out the URL and the file it was found in.


def check_urls_in_file(file_path):
    with open(file_path, "r") as file:
        content = file.read()

    # Find all URLs in the file
    urls = [url.split(">`")[0] for url in re.findall(r"https?://[^\s]+", content)]
    invalid_urls = []

    print(f"Checking URLs in file: {file_path}")

    for url in urls:
        # skip urls that end in 'download' or 'download/'
        if url.endswith("download") or url.endswith("download/"):
            print(f"Skipping {url} in file {file_path}")
            continue

        try:
            response = requests.get(url, verify=True, timeout=3, stream=True)
            if response.status_code not in (200, 403, 429):
                invalid_urls.append(url)
                print(
                    f"{url} is invalid in file {file_path} with status code {response.status_code}"
                )
            else:
                print("ok!")
        except Exception:
            invalid_urls.append(url)
            print(f"{url} is invalid in file {file_path}")

    print()

    return invalid_urls


def main():
    if not os.path.exists(MESA_DIR):
        print(f"MESA_DIR '{MESA_DIR}' does not exist.")
        sys.exit(1)

    docs_dir = os.path.join(MESA_DIR, "docs", "source")
    if not os.path.isdir(docs_dir):
        print(f"Docs directory '{docs_dir}' does not exist.")
        sys.exit(1)

    invalid_urls = {}

    for root, _, files in os.walk(docs_dir):
        for file in files:
            if file.endswith(".rst"):
                file_path = os.path.join(root, file)
                file_invalid_urls = check_urls_in_file(file_path)
                if file_invalid_urls:
                    invalid_urls[file_path] = file_invalid_urls

    if invalid_urls:
        print("\n\n====================")
        print("SUMMARY:\n")
        print("Invalid URLs found:")
        for file_path, urls in invalid_urls.items():
            print(f"In {file_path}:")
            for url in urls:
                print(f"  {url}")
        sys.exit(1)
    else:
        print("All URLs are valid.")


if __name__ == "__main__":
    main()
