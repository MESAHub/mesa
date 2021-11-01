#!/usr/bin/env python

import os
import re
import sys
from pathlib import Path

MESA_DIR = os.environ["MESA_DIR"]

# Search files for instances on an empty writing usin* as the format i.e write(*,*) 
# and replace with calls to write(*,'(A)') 
# This makes writes more portable to ifort

# Files or folders to skip
skip_folders = [
]

skip_files = [
]


def check_skip(path):
    for s in skip_folders:
        if path.is_relative_to(s):
            return True
    for s in skip_files:
        if path.name == s:
            return True
    return False


write = re.compile("^ *write\([*a-zA-Z0-9]+\,[*]\)$")

def replace(match):
    return match.string.split(',')[0] + ",'(A)')"

if len(sys.argv) > 1:
    files = sys.argv[1:]
else:
    files = Path("./").rglob("*.f90")


for file in Path("./").rglob("*.f90"):
    if check_skip(file):
        continue

    with open(file, "r") as f:
        lines = f.readlines()

    modified = False
    for ldx, line in enumerate(lines):
        if "write(" in line:
            if write.search(line):
                # Preserve unit number
                line = re.sub(write, replace, line)
                lines[ldx] = line
                modified = True

    if modified:
        with open(file, "w") as f:
            f.writelines(lines)
        print("Updated ", file)
