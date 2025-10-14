#!/usr/bin/python
import sys

TABSPACE = 3


def convert_tabs_to_spaces(file_path):
    try:
        with open(file_path, "r") as file:
            content = file.read()

        # Replace tabs with spaces
        num_tabs = content.count("\t")
        updated_content = content.replace("\t", " " * TABSPACE)

        with open(file_path, "w") as file:
            file.write(updated_content)

        print(f"Converted {num_tabs} tabs to {TABSPACE} spaces in: {file_path}")
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <file1> <file2> ...")
        sys.exit(1)

    for file_path in sys.argv[1:]:
        convert_tabs_to_spaces(file_path)
