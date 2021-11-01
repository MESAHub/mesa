#!/usr/bin/env python

import os
import sys
import re
from pathlib import Path

MESA_DIR = os.environ["MESA_DIR"]

# Search files for instances of stop 1 or stop 'str'
# and replace with calls to mesa_error. This way all error messages
# can be made unique

# Files or folders to skip
skip_folders = [
    "utils/",
    'eos/eosCMS_builder',
    'eos/eosFreeEOS_builder',
]

skip_files = [
    "run.f90",
]


def check_skip(path):
    for s in skip_folders:
        if path.is_relative_to(s):
            return True
    for s in skip_files:
        if path.name == s:
            return True
    return False


numeric_stop = re.compile("stop [0-9]+")
str_stop = re.compile("stop '.+'")


def mesa_error(message=None):
    if message is None:
        return "call mesa_error(__FILE__,__LINE__)"
    else:
        return "call mesa_error(__FILE__,__LINE__," + message + ")"

if len(sys.argv) > 1:
    files = sys.argv[1:]
else:
    files = Path("./").rglob("*.f90")


for file in files:
    if check_skip(file):
        continue

    with open(file, "r") as f:
        lines = f.readlines()

    modified = False
    for ldx, line in enumerate(lines):
        if "stop " in line:
            if numeric_stop.search(line):
                # Numeric case
                line = re.sub(numeric_stop, mesa_error(), line)
            elif str_stop.search(line):
                # String type, preserve error message
                match = str_stop.search(line).group()
                match = match[len("stop ") :]
                line = re.sub(str_stop, mesa_error(match), line)
            lines[ldx] = line
            modified = True

    if modified:
        with open(file, "w") as f:
            f.writelines(lines)
        print("Updated ", file)
