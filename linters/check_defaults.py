#!/usr/bin/env python

import functools
import operator
import os
import re
from collections.abc import MutableSet

MESA_DIR = "../"
ENABLE_TEST_SUITE_HIST_CHECKS = True
ENABLE_TEST_SUITE_PROF_CHECKS = True


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


def print_section(header):
    """Display output section header"""
    print(f"\n\n*** {header} ***\n")


def print_options(options):
    """Print a set of options"""
    for o in sorted(options):
        print(f"   {o}")


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
    with open(os.path.join(MESA_DIR, filename), "r") as f:
        lines = f.readlines()
    
    return lines


def get_inc(filename):
    # extract options from a inc file
    lines = load_file(filename)
    
    # Remove line continutaion characters
    lines = [i.replace("&", "").strip() for i in lines if i]
    
    # Remove type defintion (i.e real(dp) :: x) leaves just x
    # as well as anything that starstwith a comment or has a comment embeded
    # in it
    for idl, line in enumerate(lines):
        if "::" in line:
            lines[idl] = line.split("::")[1].strip()
    
    lines = [i.split(",") for i in lines if i]
    
    # Flatten list of lists
    lines = functools.reduce(operator.iconcat, lines, [])
    
    # Remove array sizes from variables
    lines = [line.split("(")[0] for line in lines if line]
    
    # Remove comments
    lines = [line.split("!")[0] for line in lines if line]
    
    # Remove = x 
    lines = [line.split("=")[0] for line in lines if line]
    
    # Remove remaining empty strings
    lines = [line.strip() for line in lines if line]
    
    return CaseInsensitiveSet(lines)


def check_io(filename, dt, var):
    # Checks that we have both dt% variable = variable and variable = dt%
    # variable
    # Where dt is the derived type (s for star_job, b for binary_job etc)
    
    lines = load_file(filename)
    
    m1 = f"{dt}% {var} = {var}"
    m2 = f"{var} = {dt}% {var}"
    
    regexp_1 = f"^[ \t]*[ ]?{dt}[ ]*%[ ]*{var}[ ]*=[ ]*{var}"
    regexp_2 = f"^[ \t]*[ ]?{var}[ ]*=[ ]*{dt}[ ]*%[ ]*{var}"
    
    rc1 = re.compile(regexp_1, flags=re.IGNORECASE)
    rc2 = re.compile(regexp_2, flags=re.IGNORECASE)
    
    match_m1 = False
    match_m2 = False
    r1 = []
    r2 = []
    
    for line in lines:
        if rc1.match(line):
            match_m1 = True
        if rc2.match(line):
            match_m2 = True
    
    if not match_m1:
        r1.append(m1)
    
    if not match_m2:
        r2.append(m2)
    
    return r1, r2


def run_checks(inc_file, defaults_file, io_file, dt, module):
    print_section(module)
    cinc = get_inc(inc_file)
    
    cdef = get_defaults(defaults_file)
    
    false_positives = (
        f"extra_{module}_inlist1_name",
        f"extra_{module}_inlist2_name",
        f"extra_{module}_inlist3_name",
        f"extra_{module}_inlist4_name",
        f"extra_{module}_inlist5_name",
        f"read_extra_{module}_inlist1",
        f"read_extra_{module}_inlist2",
        f"read_extra_{module}_inlist3",
        f"read_extra_{module}_inlist4",
        f"read_extra_{module}_inlist5",
        f"save_{module}_namelist",
        f"{module}_namelist_name",
    )
    
    print_section("Things in include but not defaults")
    
    print_options(cinc - cdef - false_positives)
    
    print_section("Things in defaults but not include")
    
    print_options(cdef - cinc - false_positives)
    
    print_section("Things in defaults but not io")
    
    m1 = []
    m2 = []
    for i in cdef - false_positives:
        r1, r2 = check_io(io_file, dt, i)
        m1.extend(r1)
        m2.extend(r2)
    
    for i in m1:
        print(i)
    
    print()
    
    for i in m2:
        print(i)
    
    print()


if __name__ == "__main__":
    run_checks("star_data/private/star_controls.inc",
               "star/defaults/controls.defaults", "star/private/ctrls_io.f90",
               "s", "controls")
    run_checks("star_data/private/star_controls_dev.inc",
               "star/defaults/controls_dev.defaults",
               "star/private/ctrls_io.f90", "s", "controls")
    
    run_checks("star_data/private/star_job_controls.inc",
               "star/defaults/star_job.defaults",
               "star/private/star_job_ctrls_io.f90", "s% job", "star_job")
    run_checks("star_data/private/star_job_controls_dev.inc",
               "star/defaults/star_job_dev.defaults",
               "star/private/star_job_ctrls_io.f90", "s% job", "star_job")
