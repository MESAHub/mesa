"""This python file was made exlusively for calulcating bi-cubic splines of LANL ATOMIC OPLIB opacity tables, not OPAL or OP (which contain 9.999 values in different table locations)"""

"""lasted updated by EbF 12-2023"""

from pathlib import Path
import numpy as np
import scipy as sp
import sys


def smooth_opac(table_path: Path):
    # importing files from history.data in each model folder
    table = table_path.read_text().splitlines()
    header = table[0:5]
    header[2] = header[2].replace(
        " 74    3.764000    9.065000", "213    3.750000    9.050000"
    )
    logr_vals = np.loadtxt(table[5:6])
    data = np.loadtxt(table[6:])
    logt_vals = data[:, 0]
    op_values = data[:, 1:]

    """now we remove the 9.9999e+0 values from opacity files"""
    for i in range(0, len(logr_vals)):
        interp_number = 0
        for j in range(0, len(logt_vals)):
            if (op_values[j, i] >= 9.99e0) and (interp_number == 0):
                interp_number = 1

                dx_to_interp = logt_vals[j - 1] - logt_vals[j - 2]
                dy_to_interp = op_values[j - 1, i] - op_values[j - 2, i]
                op_values[j, i] = (dy_to_interp / dx_to_interp) * (
                    logt_vals[j] - logt_vals[j - interp_number]
                ) + op_values[j - interp_number, i]
            elif (op_values[j, i] >= 9.99e0) and (interp_number >= 1):
                interp_number = interp_number + 1  # add 1 to interpolation number
                op_values[j, i] = (dy_to_interp / dx_to_interp) * (
                    logt_vals[j] - logt_vals[j - interp_number]
                ) + op_values[j - interp_number, i]

    """Generate a spline"""
    spline = sp.interpolate.RectBivariateSpline(
        logt_vals, logr_vals, op_values, kx=3, bbox=[3.750, None, None, None], ky=3, s=0
    )

    """spline dimensions"""
    tarr = np.linspace(3.75, 9.05, 213)
    rarr = np.linspace(-8, 1.5, 39)
    tgrid, rgrid = np.meshgrid(tarr, rarr, indexing="ij")

    kap_data_interp = spline(tgrid, rgrid, grid=False)

    """output to new spline file"""

    output_name = output_file_dir / table_path.name
    with output_name.open("w") as output:
        # write header
        output.writelines(header)

        # write log R header
        index = " ".join((f"{logr: 8.2f}" for logr in rarr))
        index = "       " + index
        output.write(index)
        output.write("\n")
        output.write("\n")
        # write logt + kap_vals
        for i in range(0, len(tarr)):
            kap_row_val = kap_data_interp[i, :]
            row = " ".join((f"{kap: 8.4f}" for kap in kap_row_val))
            row = f"{tarr[i]:.3f} {row}"
            output.write(row)
            output.write("\n")

opacity_file_dir = Path(sys.argv[1])
output_file_dir = Path(sys.argv[2])

for data_file in opacity_file_dir.glob("*.data"):
    smooth_opac(data_file)
