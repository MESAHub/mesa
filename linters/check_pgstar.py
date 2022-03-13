import os
import glob

import check_columns as cc

MESA_DIR = os.environ["MESA_DIR"]

profile_options = cc.CaseInsensitiveSet(
    [
        "Abundance_xaxis_name",
        "Power_xaxis_name",
        "Mixing_xaxis_name",
        "Dynamo_xaxis_name",
        "Mode_Prop_xaxis_name",
        "Summary_Burn_xaxis_name",
        "Summary_Profile_xaxis_name",
        "Summary_Profile_name",
        "Profile_Panels1_xaxis_name",
        "Profile_Panels2_xaxis_name",
        "Profile_Panels3_xaxis_name",
        "Profile_Panels4_xaxis_name",
        "Profile_Panels5_xaxis_name",
        "Profile_Panels6_xaxis_name",
        "Profile_Panels7_xaxis_name",
        "Profile_Panels8_xaxis_name",
        "Profile_Panels9_xaxis_name",
        "Profile_Panels1_yaxis_name",
        "Profile_Panels1_other_yaxis_name",
        "Profile_Panels2_yaxis_name",
        "Profile_Panels2_other_yaxis_name",
        "Profile_Panels3_yaxis_name",
        "Profile_Panels3_other_yaxis_name",
        "Profile_Panels4_yaxis_name",
        "Profile_Panels4_other_yaxis_name",
        "Profile_Panels5_yaxis_name",
        "Profile_Panels5_other_yaxis_name",
        "Profile_Panels6_yaxis_name",
        "Profile_Panels6_other_yaxis_name",
        "Profile_Panels7_yaxis_name",
        "Profile_Panels7_other_yaxis_name",
        "Profile_Panels8_yaxis_name",
        "Profile_Panels8_other_yaxis_name",
        "Profile_Panels9_yaxis_name",
        "Profile_Panels9_other_yaxis_name",
    ]
)

# Ignores color_magnitude plots as its hard to work out if the names are in the history output
history_options = cc.CaseInsensitiveSet(
    [
        "Summary_History_name",
        "Kipp_xaxis_name",
        "rti_xaxis_name",
        "History_Track1_xname",
        "History_Track2_xname",
        "History_Track3_xname",
        "History_Track4_xname",
        "History_Track5_xname",
        "History_Track6_xname",
        "History_Track7_xname",
        "History_Track8_xname",
        "History_Track9_xname",
        "History_Track1_yname",
        "History_Track2_yname",
        "History_Track3_yname",
        "History_Track4_yname",
        "History_Track5_yname",
        "History_Track6_yname",
        "History_Track7_yname",
        "History_Track8_yname",
        "History_Track9_yname",
        "History_Panels1_xaxis_name",
        "History_Panels2_xaxis_name",
        "History_Panels3_xaxis_name",
        "History_Panels4_xaxis_name",
        "History_Panels5_xaxis_name",
        "History_Panels6_xaxis_name",
        "History_Panels7_xaxis_name",
        "History_Panels8_xaxis_name",
        "History_Panels9_xaxis_name",
        "History_Panels1_yaxis_name",
        "History_Panels1_other_yaxis_name",
        "History_Panels1_points_name",
        "History_Panels2_yaxis_name",
        "History_Panels2_other_yaxis_name",
        "History_Panels2_points_name",
        "History_Panels3_yaxis_name",
        "History_Panels3_other_yaxis_name",
        "History_Panels3_points_name",
        "History_Panels4_yaxis_name",
        "History_Panels4_other_yaxis_name",
        "History_Panels4_points_name",
        "History_Panels5_yaxis_name",
        "History_Panels5_other_yaxis_name",
        "History_Panels5_points_name",
        "History_Panels6_yaxis_name",
        "History_Panels6_other_yaxis_name",
        "History_Panels6_points_name",
        "History_Panels7_yaxis_name",
        "History_Panels7_other_yaxis_name",
        "History_Panels7_points_name",
        "History_Panels8_yaxis_name",
        "History_Panels8_other_yaxis_name",
        "History_Panels8_points_name",
        "History_Panels9_yaxis_name",
        "History_Panels9_other_yaxis_name",
        "History_Panels9_points_name",
        "Text_Summary1_name",
        "Text_Summary2_name",
        "Text_Summary3_name",
        "Text_Summary4_name",
        "Text_Summary5_name",
        "Text_Summary6_name",
        "Text_Summary7_name",
        "Text_Summary8_name",
        "Text_Summary9_name",
    ]
)


profile_false_positives = cc.CaseInsensitiveSet(["Abundance", "Mixing", "Power"])

history_false_positives = cc.CaseInsensitiveSet(
    [
        "Mass",
        "lg_Mdot",
        "H_cntr",
        "H_rich",
        "He_cntr",
        "He_core",
        "C_cntr",
        "C_core",
        "O_cntr",
        "O_core",
        "Ne_cntr",
        "Ne_core",
        "Si_cntr",
        "Si_core",
        "Fe_cntr",
        "Fe_core",
        "eta_cntr",
        "lg_Lnuc",
        "lg_LNeu",
        "lg_Lphoto",
        "zones",
        "retries",
        "v_div_cs",
        "center c12",
        "center n14",
        "center o16",
        "surface_c12",
        "surface_o16",
    ]
)


profile_defaults = cc.get_profile_columns()
history_defaults = cc.get_history_columns()


def check_pgstar(filename, options, defaults, false_positives):
    if not os.path.isfile(filename):
        return

    with open(filename, "r") as f:
        lines = f.readlines()

    result = []
    for line in lines:
        line = line.strip()
        if (
            not len(line)
            or line.startswith("!")
            or line.startswith("&")
            or line.startswith("/")
        ):
            continue

        line = line.split("!")[0]

        column, *value = line.split("=")  # Things like x = 'mass'
        value = "".join(
            value
        )  # Handles things with several = signs i.e semiconvection_option = 'Langer_85 mixing; gradT = gradr'

        column_check = column.split("(")[0]  # Things like x(1) = 'mass'

        value = value.strip().replace("'", "").replace('"', "")

        if not len(value):
            continue

        if "(:)" in column:
            continue

        if column_check in options:
            if value not in defaults and value not in false_positives:
                result.append((column, value))

    return result


def check_all_pgstars(options, defualts, false_postives):
    for filename in glob.glob(
        os.path.join(MESA_DIR, "star", "test_suite", "*", "inlist*")
    ):
        values = check_pgstar(filename, options, defualts, false_postives)
        if values is None:
            continue
        if len(values):
            print(f"\n\n*** {filename} ***\n")
            for i in values:
                print(*i)
            print()
            print()


def check_all_history_pgstars():
    check_all_pgstars(history_options, history_defaults, history_false_positives)


def check_all_profile_pgstars():
    check_all_pgstars(profile_options, profile_defaults, profile_false_positives)


if __name__ == "__main__":
    check_all_history_pgstars()
    check_all_profile_pgstars()
