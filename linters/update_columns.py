#!/usr/bin/env python

import os
import re
import glob
import sys
from collections.abc import MutableSet
from pathlib import Path

MESA_DIR = os.environ["MESA_DIR"]

# This syncs up the history and profile files in each test case with the default versions
# leaving options enabled.


# inspiration from https://stackoverflow.com/a/27531275
class CaseInsensitiveSet(MutableSet):
    def __init__(self, iterable):
        self._values = {}
        self._fold = str.casefold
        for v in iterable:
            self.add(v)

    def __repr__(self):
        return repr(self._values.values())

    def __contains__(self, value):
        return self._fold(value) in self._values

    def __iter__(self):
        return iter(self._values.values())

    def __len__(self):
        return len(self._values)

    def add(self, value):
        self._values[self._fold(value)] = value

    def discard(self, value):
        v = self._fold(value)
        if v in self._values:
            del self._values[v]


def get_options(filename, regexp):
    """Return a set of MESA option names"""
    r = regexp
    with open(os.path.join(MESA_DIR, filename)) as f:
        matches = regexp.finditer(f.read())
    return CaseInsensitiveSet(m.group(1) for m in matches)


def get_columns(filename, regexp):
    """Return a set of MESA column names"""
    with open(os.path.join(MESA_DIR, filename)) as f:
        lines = f.readlines()
    matches = []
    for line in lines:
        m = regexp.match(line)
        if m is not None:
            matches.append(m.group(1))
    return CaseInsensitiveSet(matches)


def print_section(header):
    """Display output section header"""
    print(f"\n\n*** {header} ***\n")


def print_options(options):
    """Print a set of options"""
    for o in sorted(options):
        print(f"   {o}")

match_comments = re.compile("^[ \t]*!?[ ]?(\w+[ 0-9]*?)[ ^t]*(!.*)?$")
match_no_comment = re.compile("^[ \t]*(\w+[ 0-9]*?)[ ^t]*(!.*)?$")
match_uncomment = re.compile('(^[\s]*?)([!])(.*)')

def match_columns(filename, comments=True):
    # extract column names from history_columns.list

    # these lines look like:
    #  ! star_mass ! in Msun units
    #  ? ^^^^^^^^^ ???????????????
    # that is, they may or may not be commented out
    # and may or may not have a trailing comment

    if comments:
        regexp = match_comments
    else:
        regexp = match_no_comment 

    return get_columns(filename, regexp)


def update(default, test_suite, special_cases = None, debug=False):

    # Load default file:
    with open(default) as f:
        lines = f.readlines()

    # See match_columns for explanation of regex
    r = re.compile(match_comments)

    for file in test_suite:
        # Load enabled options for test cases columns file
        test_case = match_columns(file, False)
        with open(file, "w") as f:
            for line in lines:
                # For each line in the default file we check if it was enabled in the test case version
                match = r.match(line)
                if match:
                    if match.group(1) in test_case:
                        # If enabled make sure it stays enabled when we write the file back out
                        groups = match_uncomment.match(line)
                        # We want the first ! before the option but not ! after the option
                        if groups is not None:
                            if groups[2] == '!':
                                line = ''.join([groups[1],groups[3],'\n'])

                        test_case.discard(match.group(1))

                print(line, file=f, end="")

            # Add back in things we may miss like burning_regions 40
            if special_cases is not None and test_case:
                print(file=f)
                print('   !## Extras',file=f)
                for case in special_cases:
                    for c in test_case:
                        if c.startswith(case):
                            print(f'      {c}',file=f)

        # Debug anything leftover
        if debug:
            if test_case:
                print(file,test_case)
                print()



def update_history():
    """Run updates on history columns"""

    files = glob.glob(
        os.path.join(MESA_DIR, "star", "test_suite", "*", "history_columns.list")
    )

    default = os.path.join(MESA_DIR, "star", "defaults", "history_columns.list")

    # Special case things that are name number
    special_cases = set(['burning_regions','mixing_regions','mix_relr_regions',
                        'burn_relr_regions'])


    update(default, files, special_cases)


def update_profile():
    """Run updates on profile columns"""

    files = glob.glob(
        os.path.join(MESA_DIR, "star", "test_suite", "*", "profile_columns.list")
    )

    default = os.path.join(MESA_DIR, "star", "defaults", "profile_columns.list")

    update(default, files)


if __name__ == "__main__":
    update_history()
    update_profile()
