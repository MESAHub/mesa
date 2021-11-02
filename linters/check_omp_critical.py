#!/usr/bin/env python

import os
import sys
import re
from pathlib import Path

MESA_DIR = os.environ["MESA_DIR"]

# Search files for unnamed omp critical blocks

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


omp_crit_start = re.compile("omp critical *$",re.IGNORECASE)
omp_crit_end = re.compile("omp end critical *$",re.IGNORECASE)

if len(sys.argv) > 1:
    files = sys.argv[1:]
else:
    files = Path("./").rglob("*.f90")

def replace_start(line,file, num):
    whitespace = ' '*line.index('!')
    return whitespace + '!$omp critical (' + file.stem + '_' + str(num) + ')\n' 

def replace_end(line,file, num):
    whitespace = ' '*line.index('!')
    return whitespace + '!$omp end critical (' + file.stem + '_' + str(num) + ')\n' 

for file in files:
    num = 1
    if check_skip(file):
        continue

    with open(file, "r") as f:
        lines = f.readlines()

    modified = False
    for ldx, line in enumerate(lines):
        if omp_crit_start.search(line):
            line = replace_start(line,file,num)
            lines[ldx] = line
            modified = True
            # Find the end statement
            for ldx2,line2 in enumerate(lines[ldx:]):
                if omp_crit_end.search(line2):
                    line2 = replace_end(line2, file, num)
                    lines[ldx+ldx2] = line2
                    break
            num = num+1

    if modified:
        with open(file, "w") as f:
            f.writelines(lines)
        print("Updated ", file)
