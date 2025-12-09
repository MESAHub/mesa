#!/usr/bin/env python3

import glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import mesa_reader as mr

from static_HISTORY_check import MesaView, read_header_columns, setup_hr_diagram_params


def make_3d_cmd(history_file="../LOGS/history.data"):
    md = mr.MesaData(history_file)
    md = MesaView(md, 5)

    all_cols, filter_columns = read_header_columns(history_file)
    hr_color, hr_mag, hr_xlabel, hr_ylabel, color_index = setup_hr_diagram_params(
        md, filter_columns
    )

    interp_rad = md.Interp_rad

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    ax.scatter(hr_color, hr_mag, interp_rad, s=5)

    ax.set_xlabel(hr_xlabel)
    ax.set_ylabel(hr_ylabel)
    ax.set_zlabel("Interp_rad")

    ax.invert_yaxis()

    plt.tight_layout()
    plt.show()


def main():
    history_file = glob.glob("../LOGS/history.data")[0]
    make_3d_cmd(history_file)


if __name__ == "__main__":
    main()
