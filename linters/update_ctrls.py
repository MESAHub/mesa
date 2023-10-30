#!/usr/bin/env python

# This file is designed to handle moving controls from s% xx to s% ctrl %xx
# While this replicates alot of check_defaults.py we want this to be standalone
# and only depend on the python stdlib. This way users can easily use it to port their
# run_star_extras.f90 files without needing to worry about python packaging.

# Usage: python update_ctrls.up file1.f90 file2.f90 ....

# Note this only works for s% (or s %) it does not work if you renamed the star_type varaible 
# to something other than s, for instance in the binary module.

import os
import re
from collections.abc import MutableSet
import functools
import operator
import sys


MESA_DIR = os.environ["MESA_DIR"]

ctrls_files = [ os.path.join("star_data","private","star_controls.inc"),
                os.path.join("star_data","private","star_controls_dev.inc")
            ]

CRTL_NAME = 's% ctrl% '

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

    def items(self):
        return self._values.items()

    def add(self, value):
        if isinstance(value, CaseInsensitiveSet):
            for k,v in value.items():
                self._values[self._fold(k)] = v
        else:
            self._values[self._fold(value)] = value

        
    def discard(self, value):
        v = self._fold(value)
        if v in self._values:
            del self._values[v]


def get_options(filename, regexp):
    """Return a set of MESA option names"""
    r = re.compile(regexp)
    with open(os.path.join(MESA_DIR, filename)) as f:
        matches = r.finditer(f.read())
    return CaseInsensitiveSet(m.group(1) for m in matches)


def get_columns(filename, regexp):
    """Return a set of MESA column names"""
    r = re.compile(regexp)
    with open(os.path.join(MESA_DIR, filename)) as f:
        lines = f.readlines()
    matches = []
    for line in lines:
        m = r.match(line)
        if m is not None:
            matches.append(m.group(1))
    return CaseInsensitiveSet(matches)

def get_defaults(filename):
    # extract column names from defaults file

    # these lines look like:
    #  ! initial_mass = 1
    #  ? ^^^^^^^^^
    # that is, they may or may not be commented out
    # and may or may not have a( )
    # and may or may not have space before a =

    regexp = "^[ \t]*[ ]?(\w+)(\(.*\))*[ ^t]*="

    return get_columns(filename, regexp)

def load_file(filename):
    with open(os.path.join(MESA_DIR, filename),"r") as f:
        lines = f.readlines()

    return lines


def get_inc(filename):
    # extract options from a inc file
    lines = load_file(filename)

    # Remove line continutaion characters
    lines = [i.replace("&","").strip() for i in lines if i]

    # Remove type defintion (i.e real(dp) :: x) leaves just x
    # as well as anything that starstwith a comment or has a comment embeded in it
    for idl,line in enumerate(lines):
        if "::" in line:
            lines[idl] = line.split("::")[1].strip()

    lines = [i.split(",") for i in lines  if i]

    # Flatten list of lists
    lines = functools.reduce(operator.iconcat, lines, [])

    # Remove array sizes from variables
    lines  = [line.split("(")[0] for line in lines if line]
    
    # Remove comments
    lines = [line.split("!")[0] for line in lines if line]

    # Remove = x 
    lines  = [line.split("=")[0] for line in lines if line]

    # Remove remaining empty strings
    lines  = [line.strip() for line in lines if line]

    return CaseInsensitiveSet(lines)

# Load controls names
cinc = get_inc(ctrls_files[0])

for f in ctrls_files[1:]:
    cinc.add(get_inc(f))


def update(filename):
    try:
        lines = load_file(filename)
    except (UnicodeDecodeError, IsADirectoryError):
        return
        
    " s[0 or more space] % [0 or more space] [1 or more character or number or _]"
    # This wont match when s has been renamed 
    regex_all = "(s[ \t]?[a-zA-Z0-9_]?%[ \t]?[a-zA-Z0-9_]*)"

    " s [0 or more space] % [0 or more space] "
    regex_s = 's[ \ta-zA-Z0-9_]?%[ \t]?'

    r_all = re.compile(regex_all)
    r_s = re.compile(regex_s)

    for idl,line in enumerate(lines):
        # Split on s% something
        line_split = re.split(regex_all,line)
        for idm,match in enumerate(line_split):
            # Remove the s% so we can check if the variable is a control
            var = match.replace('s%','').strip()
            if var in cinc:
                # If it is a control then replace s% with CRTL_NAME
                line_split[idm] = re.sub(regex_s,CRTL_NAME,match)
        lines[idl] = ''.join(line_split)

    with open(filename,'w') as f:
        f.writelines(lines)

if __name__ == "__main__":
    for i in sys.argv[1:]:
        update(i)

# Run over MESA_DIR

# python3 linters/update_ctrls.py star/test/src/*
# python3 linters/update_ctrls.py star/work/src/*
# python3 linters/update_ctrls.py star/job/*
# python3 linters/update_ctrls.py star/p*/*.f90 
# python3 linters/update_ctrls.py star/test_suite/*/src/*.f90
# python3 linters/update_ctrls.py star/other/*.f90
# python3 linters/update_ctrls.py */test_suite/*/src/*.inc
# python3 linters/update_ctrls.py */test_suite/*/src/*/*.inc

# python3 linters/update_ctrls.py binary/test/src/*
# python3 linters/update_ctrls.py binary/work/src/*
# python3 linters/update_ctrls.py binary/job/*
# python3 linters/update_ctrls.py binary/p*/*.f90 
# python3 linters/update_ctrls.py binary/test_suite/*/src/*.f90
# python3 linters/update_ctrls.py binary/other/*.f90

# python3 linters/update_ctrls.py astero/test/src/*
# python3 linters/update_ctrls.py astero/work/src/*
# python3 linters/update_ctrls.py astero/job/*
# python3 linters/update_ctrls.py astero/p*/*.f90 
# python3 linters/update_ctrls.py astero/test_suite/*/src/*.f90
# python3 linters/update_ctrls.py astero/other/*.f90
