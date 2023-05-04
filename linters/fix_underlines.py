#!/usr/bin/env python

import os

MESA_DIR = "../"

# files that should be checked
DEFAULTS_FILES = (
    "astero/defaults/astero_search.defaults",
    "astero/defaults/astero_pgstar.defaults",
    "star/defaults/star_job.defaults",
    "star/defaults/controls.defaults",
    "star/defaults/star_job_dev.defaults",
    "star/defaults/controls_dev.defaults",
    "star/defaults/pgstar.defaults",
    "binary/defaults/binary_job.defaults",
    "binary/defaults/binary_controls.defaults",
    "binary/defaults/pgbinary.defaults",
    "eos/defaults/eos.defaults",
    "kap/defaults/kap.defaults",
)


def fix_underlines(filename):

    path = os.path.join(MESA_DIR, filename)

    with open(path, "r") as f:
        lines = f.readlines()

    # go through lines
    prev_line = None

    with open(path, "w") as f:
        for line in lines:
            if line.endswith("~\n"):
                # in an underline line
                len_prev_line = len(prev_line)
                len_line = len(line)
                if len_prev_line > len_line:
                    # label too short
                    print("fixing label (underline too short)\n")
                    print(prev_line.rstrip())
                    print(line)
                    line = line[:-1] + "~" * (len_prev_line - len_line) + "\n"
                if len_prev_line < len_line:
                    # label too long
                    print("fixing label (underline too long)\n")
                    print(prev_line.rstrip())
                    print(line)
                    diff = (len_line - len_prev_line) + 1
                    line = line[:-diff] + "\n"
            f.write(line)
            prev_line = line


if __name__ == "__main__":
    for f in DEFAULTS_FILES:
        fix_underlines(f)
