#!/usr/bin/env python

import os
import re
from collections.abc import MutableSet

MESA_DIR = "../"


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


known_false_positives = {
    "max_years_for_timestep",
    "history_names_dict",
    "report_ierr",
    "Tsurf_factor",
    "other_photo_read",
    "job",
    "tau_factor",
    "opacity_factor",
    "binary_id",
    "chem_id",
    "pgstar_hist",
    "profile_column_spec",
    "history_values",
    "nameofequ",
    "nameofvar",
    "net_iso",
    "model_profile_filename",
    "include_binary_history_in_log_file",
    "i_dj_rot_dt",  # Shadows s% i_j_rot
    "i_equ_w_div_wc",  # Shadows s% i_w_div_wc
    "have_j_rot",  # Set from s% rotation_flag
    "nz_old",
}


def print_section(header):
    """Display output section header"""
    print(f"\n\n*** {header} ***\n")


def print_variables(variables):
    """Print a set of variables"""
    for v in sorted(variables):
        print(f"   {v}")


def get_star_info_variables(filename):
    """Extract names of variables living on s%"""

    regexp = "s% ?([A-Za-z0-9_]+)"

    with open(os.path.join(MESA_DIR, filename)) as f:
        source = f.read()

    return CaseInsensitiveSet(re.findall(regexp, source))


def get_photo_in_variables():
    """Extract variables that are (probably) read in photos"""

    filename = "star/private/photo_in.f90"
    v = get_star_info_variables(filename)

    return v - known_false_positives


def get_photo_out_variables():
    """Extract variables that are (probably) written in photos"""

    filename = "star/private/photo_out.f90"
    v = get_star_info_variables(filename)

    known_false_positives = {
        "other_photo_write",
    }

    return v - known_false_positives


def get_star_data_set_input_variables():
    """Check that the star_data input variables are read by photo_in"""

    filename = "star_data/public/star_data_step_input.inc"

    regexp = "([A-Za-z0-9_]+)"

    with open(os.path.join(MESA_DIR, filename)) as f:
        lines = f.readlines()

    # world's hackiest "parser"

    input_vars = []

    for line in lines:

        # first, lose white space
        line = line.strip()

        # if empty, skip
        if not line:
            continue

        # if a comment line, skip
        if line.startswith("!"):
            continue

        # if line has a comment, discard
        if "!" in line:
            line, _ = line.split("!", maxsplit=1)

        # if line has a type definition, discard
        if "::" in line:
            _, line = line.split("::", maxsplit=1)

        # split on commas
        if "," in line:
            input_vars.extend([v.strip() for v in line.split(",")])
        else:
            input_vars.append(line.strip())

    # throw away stuff that doesn't look like vars and return
    final_list = []
    for v in input_vars:
        m = re.match(regexp, v)
        # mostly won't match &
        if m is not None:
            final_list.append(m.group(0))

    return CaseInsensitiveSet(final_list) - known_false_positives


def check_in_out():
    """Check that the same variables appear in photo_{in,out}"""

    v_in = get_photo_in_variables()
    v_out = get_photo_out_variables()

    print_section("In, not out")
    print_variables(v_in - v_out)

    print_section("Out, not in")
    print_variables(v_out - v_in)


def check_step_input():

    sd_vars = get_star_data_set_input_variables()
    pi_vars = get_photo_in_variables()

    print_section("step input vars not read in photo")
    print_variables(sd_vars - pi_vars)


if __name__ == "__main__":

    print("=== Photo In/Out Comparison ===")
    check_in_out()
    print()

    print("=== Step Input / Photo In Comparison ===")
    check_step_input()
