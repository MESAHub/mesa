import numpy as np
from scipy.sparse import lil_matrix, csc_matrix
from scipy.sparse.linalg import splu, lsqr, lsmr, eigsh

###############################################################################
# 1) Parsing functions
###############################################################################
def parse_dblk_file(filename, nvar, nz):
    dblk = np.zeros((nz, nvar, nvar), dtype=np.float64)
    with open(filename, 'r') as f:
        lines = f.readlines()
    idx = 1
    zone_count = 0
    while idx < len(lines) and zone_count < nz:
        zone_line = lines[idx].strip()
        idx += 1
        if not zone_line.startswith("Zone:"):
            raise ValueError(f"Expected 'Zone:' line, got: {zone_line}")
        for i in range(nvar):
            float_line = lines[idx].strip()
            idx += 1
            row_vals = float_line.split()
            if len(row_vals) < nvar:
                raise ValueError(f"Not enough floats in line: {float_line}")
            dblk[zone_count, i, :] = np.array(row_vals[:nvar], dtype=np.float64)
        zone_count += 1
    if zone_count != nz:
        raise ValueError(f"Finished parsing but did not find all {nz} zones in {filename}.")
    return dblk

def parse_ublk_file(filename, nvar, nz):
    ublk = np.zeros((nz, nvar, nvar), dtype=np.float64)
    with open(filename, 'r') as f:
        lines = f.readlines()
    idx = 1
    zone_count = 0
    while idx < len(lines) and zone_count < (nz - 1):
        zone_line = lines[idx].strip()
        idx += 1
        if not zone_line.startswith("Zone:"):
            raise ValueError(f"Expected 'Zone:' line, got: {zone_line}")
        for i in range(nvar):
            float_line = lines[idx].strip()
            idx += 1
            row_vals = float_line.split()
            if len(row_vals) < nvar:
                raise ValueError(f"Not enough floats in line: {float_line}")
            ublk[zone_count, i, :] = np.array(row_vals[:nvar], dtype=np.float64)
        zone_count += 1
    if zone_count != (nz - 1):
        raise ValueError(f"Didn't read all upper blocks (expected {nz-1}).")
    return ublk

def parse_lblk_file(filename, nvar, nz):
    lblk = np.zeros((nz, nvar, nvar), dtype=np.float64)
    with open(filename, 'r') as f:
        lines = f.readlines()
    idx = 1
    zone_count = 1
    while idx < len(lines) and zone_count < nz:
        zone_line = lines[idx].strip()
        idx += 1
        if not zone_line.startswith("Zone:"):
            raise ValueError(f"Expected 'Zone:' line, got: {zone_line}")
        for i in range(nvar):
            float_line = lines[idx].strip()
            idx += 1
            row_vals = float_line.split()
            if len(row_vals) < nvar:
                raise ValueError(f"Not enough floats in line: {float_line}")
            lblk[zone_count, i, :] = np.array(row_vals[:nvar], dtype=np.float64)
        zone_count += 1
    if zone_count != nz:
        raise ValueError(f"Didn't read all lower blocks (expected zones up to {nz}).")
    return lblk

def parse_residuals_file(filename, nvar, nz):
    b = np.zeros(nz * nvar, dtype=np.float64)
    with open(filename, 'r') as f:
        lines = f.readlines()
    idx = 1
    zone_count = 0
    while idx < len(lines) and zone_count < nz:
        zone_line = lines[idx].strip()
        idx += 1
        if not zone_line.startswith("Zone:"):
            raise ValueError(f"Expected 'Zone:' line, got: {zone_line}")
        float_line = lines[idx].strip()
        idx += 1
        row_vals = float_line.split()
        if len(row_vals) < nvar:
            raise ValueError(f"Not enough floats in line: {float_line}")
        offset = zone_count * nvar
        b[offset:offset + nvar] = np.array(row_vals[:nvar], dtype=np.float64)
        zone_count += 1
    if zone_count != nz:
        raise ValueError(f"Did not read all {nz} zones in residuals file.")
    return b

###############################################################################
# 2) Build the global sparse matrix
###############################################################################
def build_sparse_matrix_from_blocks(dblk, ublk, lblk, nvar, nz):
    size = nvar * nz
    J_lil = lil_matrix((size, size), dtype=np.float64)
    for k in range(nz):
        row_block = k * nvar
        col_block = k * nvar
        J_lil[row_block:row_block + nvar, col_block:col_block + nvar] = dblk[k]
        if k < nz - 1:
            J_lil[row_block:row_block + nvar,
                  col_block + nvar:col_block + 2 * nvar] = ublk[k]
        if k > 0:
            J_lil[row_block:row_block + nvar,
                  col_block - nvar:col_block] = lblk[k]
    return J_lil.tocsc()

###############################################################################
# 3) Main routine / driver
###############################################################################
def main():
    file_dblk = 'dblk_output.txt'
    file_ublk = 'ublk_output.txt'
    file_lblk = 'lblk_output.txt'
    file_res  = 'residuals_output.txt'
    file_sol  = 'B_output.txt'
    
    nvar = 27
    nz   = 1842
    
    print("Parsing block files...")
    dblk = parse_dblk_file(file_dblk, nvar, nz)
    ublk = parse_ublk_file(file_ublk, nvar, nz)
    lblk = parse_lblk_file(file_lblk, nvar, nz)
    
    print("Parsing residuals...")
    b = parse_residuals_file(file_res, nvar, nz)
    sol = parse_residuals_file(file_sol, nvar, nz)
    
    print("Building sparse matrix J...")
    J = build_sparse_matrix_from_blocks(dblk, ublk, lblk, nvar, nz)
    
    print("Matrix shape:", J.shape)
    print("Non-zero entries:", J.nnz)
    
    # Compute matrix condition number
    print("\nComputing condition number...")
    ew1, ev = eigsh(J, which='LM')
    ew2, ev = eigsh(J, sigma=1e-8)   #<--- takes a long time
    ew1 = abs(ew1)
    ew2 = abs(ew2)

    cond = ew1.max()/ew2.min()
    print("Matrix condition number:", cond )
    
    # Calculate MESA residual norm
    print("\nComputing MESA residual norm...")
    r_mesa = J @ sol - b
    rnorm_mesa = np.linalg.norm(r_mesa)
    print(f"MESA residual norm = {rnorm_mesa:e}")

    # Print residuals for zone 1817
    zone_of_interest = 1817
    zone_idx = zone_of_interest - 1  # Convert to 0-based indexing
    offset = zone_idx * nvar
    r_zone = r_mesa[offset:offset + nvar]
    print(f"\nResiduals for all variables in zone {zone_of_interest}:")
    for i, val in enumerate(r_zone, start=1):
        print(f"  Variable {i:2d}: {val:24.16e}")

if __name__ == "__main__":
    main()


