#!/usr/bin/env python

import os
import re
import glob
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


def delete_command(options, filename):
    print("for i in", " ".join(options), end=";")
    print("do ", end="")
    print('sed -i "/^\s*\!$i/d" */' + filename, end=";")
    print('sed -i "/^\s*\!\s*$i/d" */' + filename, end=";")
    print('sed -i "/^\s*$i/d" */' + filename, end=";")
    print("done")


def get_history_columns(filename="star/defaults/history_columns.list"):
    # extract column names from history_columns.list

    # these lines look like:
    #  ! star_mass ! in Msun units
    #  ? ^^^^^^^^^ ???????????????
    # that is, they may or may not be commented out
    # and may or may not have a trailing comment

    regexp = "^[ \t]*!?[ ]?(\w+)[ ^t]*(!.*)?$"

    return get_columns(filename, regexp)


def get_star_history_def(filename="star/private/star_history_def.f90"):
    # extract column names from star_history_def.f90

    # these lines look like:
    #   history_column_name(h_star_mass) = 'star_mass'
    #                                       ^^^^^^^^^
    regexp = "history_column_name\\(h_\w+\\)[ ]*=[&\s]*'(\w+)'"

    return get_options(filename, regexp)


def get_star_history(filename="star/private/history.f90"):
    # extract column names from history.f90

    # these lines look like:
    #   case(h_star_mass)
    #          ^^^^^^^^^
    regexp = "case[ ]*\\(h_(\w+)\\)"

    return get_options(filename, regexp)


def check_history():
    """Run checks on history columns"""

    vals_history = get_star_history()

    vals_star_history_def = get_star_history_def()

    vals_history_list = get_history_columns()

    # make reports

    # Values that are in history.f90 but not in star_history_def.f90
    print_section("Values that are in history.f90 but not in star_history_def.f90")
    print_options(vals_history - vals_star_history_def)

    # Values that are in star_history_def.f90 but not in history.f90
    print_section("Values that are in star_history_def.f90 but not in history.f90")

    # define false positives
    known_false_positives = {
        "burn_relr_regions",
        "burning_regions",
        "mix_relr_regions",
        "mixing_regions",
        "rti_regions",
    }

    print_options(vals_star_history_def - vals_history - known_false_positives)

    # Values that are in are in history.f90 but not in history_columns.list
    print_section(
        "Values that are in are in history.f90 but not in history_columns.list"
    )

    known_false_positives = {
        "C_cntr",
        "Fe_core",
        "H_cntr",
        "H_rich",
        "He_cntr",
        "He_core",
        "Mass",
        "N_cntr",
        "Ne_cntr",
        "O_cntr",
        "Si_cntr",
        "ejecta_M",
        "eta_cntr",
        "gam_cntr",
        "lg_Dsurf",
        "lg_Lneu",
        "lg_Lnuc",
        "lg_Lphoto",
        "lg_Mdot",
        "remnant_M",
        "retries",
        "star_age_yr",
        "v_div_cs",
        "zones",
    }

    print_options(vals_history - vals_history_list - known_false_positives)

    # Values that are in are in history_columns.list but not in history.f90
    print_section(
        "Values that are in are in history_columns.list but not in history.f90"
    )

    known_false_positives = {
        "add_abs_mag",
        "add_average_abundances",
        "add_bc",
        "add_center_abundances",
        "add_log_average_abundances",
        "add_log_center_abundances",
        "add_log_lum_band",
        "add_log_surface_abundances",
        "add_lum_band",
        "add_surface_abundances",
        "add_log_total_mass",
        "add_total_mass",
        "pp",
        "cno",
        "tri_alfa",
        "c_alpha",
        "n_alpha",
        "o_alpha",
        "ne_alpha",
        "na_alpha",
        "mg_alpha",
        "si_alpha",
        "s_alpha",
        "ar_alpha",
        "ca_alpha",
        "ti_alpha",
        "fe_co_ni",
        "c12_c12",
        "c12_o16",
        "o16_o16",
        "photo",
        "pnhe4",
        "other",
        "C_cntr",
        "Fe_core",
        "H_cntr",
        "H_rich",
        "He_cntr",
        "He_core",
        "Mass",
        "N_cntr",
        "Ne_cntr",
        "O_cntr",
        "Si_cntr",
        "ejecta_M",
        "lg_Dsurf",
        "lg_Lneu",
        "lg_Lnuc",
        "lg_Lphoto",
        "lg_Mdot",
        "remnant_M",
        "retries",
        "shock",
        "v_div_cs",
        "zones",
        "misc",
        "timescales",
        "asteroseismology",
    }

    print_options(vals_history_list - vals_history - known_false_positives)

    if ENABLE_TEST_SUITE_HIST_CHECKS:
        # Value in each test case's history_columns.list but not in star/default/history_columns.list
        for i in glob.glob(
            os.path.join(MESA_DIR, "star", "test_suite", "*", "history_columns.list")
        ):
            test_case = get_history_columns(i.removeprefix(MESA_DIR))
            result = test_case - vals_history_list - known_false_positives
            if len(result):
                print_section(
                    "Values that are in are in "
                    + i
                    + " but not in history_columns.list"
                )
                print_options(result)
                delete_command(result, "history_columns.list")


def get_profile_getval(filename="star/private/profile_getval.f90"):
    # extract column names from profile_getval.f90

    # these lines look like:
    #   case(p_zone)
    #          ^^^^
    regexp = "case[ ]*\\(p_(\w+)\\)"

    return get_options(filename, regexp)


def get_profile_def(filename="star/private/star_profile_def.f90"):

    # extract column names from star_profile_def.f90

    # these lines look like:
    #   profile_column_name(p_zone) = 'zone'
    #                                  ^^^^^
    regexp = "profile_column_name\\(p_\w+\\)[ ]*=[&\s]*'(\w+)'"

    return get_options(filename, regexp)


def get_profile_columns(filename="star/defaults/profile_columns.list"):
    # extract column names from profile_columns.list

    # these lines look like:
    #  ! star_mass ! in Msun units
    #  ? ^^^^^^^^^ ???????????????
    # that is, they may or may not be commented out
    # and may or may not have a trailing comment

    regexp = "^[ \t]*!?[ ]?(\w+)[ ^t]*(!.*)?$"

    return get_columns(filename, regexp)


def check_profile():
    """Run checks on profile columns"""

    vals_profile = get_profile_getval()

    vals_star_profile_def = get_profile_def()

    # the following general info is included in a profile file
    # it looks like profile columns, but isn't, so exclude it
    general_info = {
        "initial_mass",
        "initial_z",
        "model_number",
        "num_zones",
        "star_age",
        "time_step",
        "Teff",
        "photosphere_L",
        "photosphere_r",
        "log_surface_L",
        "log_surface_radius",
        "log_surface_temp",
        "log_center_temp",
        "log_center_density",
        "log_center_P",
        "center_eta",
        "center_h1",
        "center_he3",
        "center_he4",
        "center_c12",
        "center_n14",
        "center_o16",
        "center_ne20",
        "star_mass",
        "star_mdot",
        "star_mass_h1",
        "star_mass_he3",
        "star_mass_he4",
        "star_mass_c12",
        "star_mass_n14",
        "star_mass_o16",
        "star_mass_ne20",
        "he_core_mass",
        "c_core_mass",
        "o_core_mass",
        "si_core_mass",
        "fe_core_mass",
        "tau10_mass",
        "tau10_radius",
        "tau100_mass",
        "tau100_radius",
        "dynamic_time",
        "kh_timescale",
        "nuc_timescale",
        "power_nuc_burn",
        "power_h_burn",
        "power_he_burn",
        "power_neu",
        "h1_boundary_limit",
        "he4_boundary_limit",
        "c12_boundary_limit",
        "burn_min1",
        "burn_min2",
        "ar_alpha",
        "c12_c12",
        "c12_o16",
        "c_alpha",
        "ca_alpha",
        "fe_co_ni",
        "mg_alpha",
        "n_alpha",
        "na_alpha",
        "ne_alpha",
        "o16_o16",
        "o_alpha",
        "other",
        "photo",
        "pnhe4",
        "s_alpha",
        "si_alpha",
        "ti_alpha",
    }

    vals_profile_list = get_profile_columns() - general_info

    # make reports

    # Values that are in star_profile_def.f90 but not in profile_getval.f90
    print_section(
        "Values that are in star_profile_def.f90 but not in profile_getval.f90"
    )
    print_options(vals_star_profile_def - vals_profile)

    # Values that are in profile_getval.f90 but not in star_profile_def.f90
    print_section(
        "Values that are in profile_getval.f90 but not in star_profile_def.f90"
    )
    print_options(vals_profile - vals_star_profile_def)

    # Values that are in are in profile_columns.list but not in profile_getval.f90
    print_section(
        "Values that are in are in profile_columns.list but not in profile_getval.f90"
    )

    known_false_positives = {
        "add_abundances",
        "add_log_abundances",
        "add_reaction_categories",
        "h1",
        "he3",
        "he4",
        "c12",
        "n14",
        "o16",
        "pp",
        "cno",
        "tri_alfa",
    }

    print_options(vals_profile_list - vals_profile - known_false_positives)

    # Values that are in are in profile_getval.f90 but not in profile_columns.list
    print_section(
        "Values that are in are in profile_getval.f90 but not in profile_columns.list"
    )
    print_options(vals_profile - vals_profile_list)

    if ENABLE_TEST_SUITE_PROF_CHECKS:
        # Value in each test case's profile_columns.list but not in star/default/profile_columns.list
        for i in glob.glob(
            os.path.join(MESA_DIR, "star", "test_suite", "*", "profile_columns.list")
        ):
            test_case = get_profile_columns(i.removeprefix(MESA_DIR))
            result = (
                test_case - vals_profile_list - known_false_positives - general_info
            )
            if len(result):
                print_section(
                    "Values that are in are in "
                    + i
                    + " but not in profile_columns.list"
                )
                print_options(result)
                delete_command(result, "profile_columns.list")


if __name__ == "__main__":
    check_history()
    check_profile()
