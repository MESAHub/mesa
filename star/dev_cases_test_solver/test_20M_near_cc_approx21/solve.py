import numpy as np
from scipy.sparse import lil_matrix, csc_matrix
from scipy.sparse.linalg import splu, lsqr, lsmr

###############################################################################
# 1) Parsing functions
###############################################################################
def parse_dblk_file(filename, nvar, nz):
    """
    Parse the dblk (diagonal blocks) text file.
    
    Returns
    -------
    dblk : numpy.ndarray of shape (nz, nvar, nvar)
        dblk[k,:,:] is the diagonal block for zone k (0-based in Python).
    """
    dblk = np.zeros((nz, nvar, nvar), dtype=np.float64)
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # The first line is a header, e.g. "Diagonal Block (dblk):"
    # Then for each zone we have 1 line "Zone: k" + nvar lines of data.
    idx = 1  # start reading after the header line
    zone_count = 0
    while idx < len(lines) and zone_count < nz:
        zone_line = lines[idx].strip()
        idx += 1
        
        if not zone_line.startswith("Zone:"):
            raise ValueError(f"Expected 'Zone:' line, got: {zone_line}")
        
        # Now read nvar lines, each containing nvar floats.
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
    """
    Parse the ublk (upper diagonal blocks) text file.
    
    The file has blocks for zones 1..nz-1 in Fortran terms,
    which is 0..(nz-2) in Python 0-based indexing.
    """
    ublk = np.zeros((nz, nvar, nvar), dtype=np.float64)
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    idx = 1  # skip the first line "Upper Diagonal Block (ublk):"
    zone_count = 0
    while idx < len(lines) and zone_count < (nz - 1):
        zone_line = lines[idx].strip()
        idx += 1
        
        if not zone_line.startswith("Zone:"):
            raise ValueError(f"Expected 'Zone:' line, got: {zone_line}")
        
        # read nvar lines for that block
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
    """
    Parse the lblk (lower diagonal blocks) text file.
    
    The file has blocks for zones 2..nz in Fortran terms,
    i.e. 1..(nz-1) in Python 0-based indexing.
    """
    lblk = np.zeros((nz, nvar, nvar), dtype=np.float64)
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    idx = 1  # skip the first line "Lower Diagonal Block (lblk):"
    zone_count = 1  # zone=2 in Fortran => index=1 in Python
    while idx < len(lines) and zone_count < nz:
        zone_line = lines[idx].strip()
        idx += 1
        
        if not zone_line.startswith("Zone:"):
            raise ValueError(f"Expected 'Zone:' line, got: {zone_line}")
        
        # read nvar lines
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
    """
    Parse the residuals file. The Fortran code writes out -s%equ, i.e. negative of the MESA residual.
    
    We're reading that as b here, so we solve J x = b.
    
    Returns
    -------
    b : numpy.ndarray of shape (nz*nvar,)
        Flattened right-hand side vector for the entire star.
    """
    b = np.zeros(nz * nvar, dtype=np.float64)
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    idx = 1  # skip header 'Residuals (-s%equ):'
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
from scipy.sparse import lil_matrix

def build_sparse_matrix_from_blocks(dblk, ublk, lblk, nvar, nz):
    """
    Build the global J matrix (nz*nvar x nz*nvar) in a sparse format
    from the diagonal, upper, and lower blocks.
    """
    size = nvar * nz
    J_lil = lil_matrix((size, size), dtype=np.float64)
    
    for k in range(nz):
        # place diagonal block
        row_block = k * nvar
        col_block = k * nvar
        J_lil[row_block:row_block + nvar, col_block:col_block + nvar] = dblk[k]
        
        # place upper block if k < nz - 1
        if k < nz - 1:
            J_lil[row_block:row_block + nvar,
                  col_block + nvar:col_block + 2*nvar] = ublk[k]
        
        # place lower block if k > 0
        if k > 0:
            J_lil[row_block:row_block + nvar,
                  col_block - nvar:col_block] = lblk[k]
    
    return J_lil.tocsc()

###############################################################################
# 3) Solve the system with LU, LSQR, and LSMR
###############################################################################
from scipy.sparse.linalg import splu

def solve_via_lu(J, b):
    """
    Solve J x = b using a sparse LU factorization (direct).
    Returns (x, r, residual_norm).
    """
    lu_obj = splu(J)
    x = lu_obj.solve(b)
    # next residual = Jx - b
    r = J @ x - b
    residual_norm = np.linalg.norm(r)
    return x, r, residual_norm


from scipy.sparse.linalg import lsqr, lsmr

def solve_via_lsqr(J, b):
    """
    Solve J x = b using LSQR (iterative).
    Returns (x, r, residual_norm).
    """
    sol = lsqr(J, b, atol=1e-12, btol=1e-12, iter_lim=20000)
    x = sol[0]
    r = J @ x - b
    residual_norm = np.linalg.norm(r)
    return x, r, residual_norm

def solve_via_lsmr(J, b):
    """
    Solve J x = b using LSMR (iterative).
    Returns (x, r, residual_norm).
    """
    sol = lsmr(J, b, atol=1e-12, btol=1e-12, maxiter=20000)
    x = sol[0]
    r = J @ x - b
    residual_norm = np.linalg.norm(r)
    return x, r, residual_norm


###############################################################################
# 4) Utility to write the residual array to file in zone-by-zone format
###############################################################################
def write_residual_array_to_file(r_flat, nvar, nz, filename):
    """
    r_flat: shape (nz*nvar,)
    Write to file in a zone-based approach:
      Zone: 1
        val1 val2 ... val(nvar)
      Zone: 2
        val1 val2 ... val(nvar)
      ...
    """
    with open(filename, 'w') as f:
        for zone_index in range(nz):
            f.write(f"Zone: {zone_index+1}\n")
            offset = zone_index * nvar
            row_vals = r_flat[offset:offset + nvar]
            line_str = " ".join(f"{val:24.16e}" for val in row_vals)
            f.write(line_str + "\n")

###############################################################################
# 5) Main routine / driver
###############################################################################
def main():
    file_dblk = 'dblk_output.txt'
    file_ublk = 'ublk_output.txt'
    file_lblk = 'lblk_output.txt'
    file_res  = 'residuals_output.txt'
    
    # As described:
    nvar = 27
    nz   = 1842
    
    print("Parsing block files...")
    dblk = parse_dblk_file(file_dblk, nvar, nz)
    ublk = parse_ublk_file(file_ublk, nvar, nz)
    lblk = parse_lblk_file(file_lblk, nvar, nz)
    
    print("Parsing residuals...")
    # The Fortran file is -s%equ, so we interpret that as b here
    b = parse_residuals_file(file_res, nvar, nz)
    
    print("Building sparse matrix J...")
    J = build_sparse_matrix_from_blocks(dblk, ublk, lblk, nvar, nz)
    
    print("Matrix shape:", J.shape)
    print("Non-zero entries:", J.nnz)
    
    # ---------------- LU ----------------
    print("\nSolving via LU factorization (splu)...")
    x_lu, r_lu, rnorm_lu = solve_via_lu(J, b)
    print(f"  LU residual norm = {rnorm_lu:e}")
    
    # ---------------- LSQR ----------------
    print("\nSolving via LSQR (iterative QR-like)...")
    x_qr, r_qr, rnorm_qr = solve_via_lsqr(J, b)
    print(f"  LSQR residual norm = {rnorm_qr:e}")
    
    # ---------------- LSMR ----------------
    print("\nSolving via LSMR (iterative SVD-like)...")
    x_svd, r_svd, rnorm_svd = solve_via_lsmr(J, b)
    print(f"  LSMR residual norm = {rnorm_svd:e}")
    
    # Comparison
    print("\nComparison of final residual norms:")
    print(f"  LU   : {rnorm_lu:e}")
    print(f"  LSQR : {rnorm_qr:e}")
    print(f"  LSMR : {rnorm_svd:e}")
    
    best_method = min([('LU', rnorm_lu), ('LSQR', rnorm_qr), ('LSMR', rnorm_svd)],
                      key=lambda x: x[1])
    print(f"Best method by final residual norm: {best_method[0]}, with residual={best_method[1]:e}")
    
    # ----------------------------------------------------------------
    #  Output the next residual arrays for each method.
    #  r_lu, r_qr, r_svd each has shape (nz*nvar,).
    #
    #  a) Print zone=1817's next residual (for all 27 variables) to terminal
    #  b) Write all zones' residual arrays to file
    # ----------------------------------------------------------------
    zone_of_interest = 1817
    zone_idx = zone_of_interest - 1
    if not (0 <= zone_idx < nz):
        raise IndexError(f"Requested zone {zone_of_interest} out of range [1..{nz}].")
    
    offset = zone_idx * nvar
    
    # Residual for LU
    lu_res_zone_1817 = r_lu[offset:offset + nvar]
    print(f"\nLU next-residual for zone {zone_of_interest}:")
    for i, val in enumerate(lu_res_zone_1817, start=1):
        print(f"  var#{i:2d}: {val:24.16e}")
    
    # Residual for LSQR
    qr_res_zone_1817 = r_qr[offset:offset + nvar]
    print(f"\nLSQR next-residual for zone {zone_of_interest}:")
    for i, val in enumerate(qr_res_zone_1817, start=1):
        print(f"  var#{i:2d}: {val:24.16e}")
    
    # Residual for LSMR
    svd_res_zone_1817 = r_svd[offset:offset + nvar]
    print(f"\nLSMR next-residual for zone {zone_of_interest}:")
    for i, val in enumerate(svd_res_zone_1817, start=1):
        print(f"  var#{i:2d}: {val:24.16e}")
    
    # b) Write all zones' residual arrays to file
    write_residual_array_to_file(r_lu,  nvar, nz, "residuals_next_LU.txt")
    write_residual_array_to_file(r_qr,  nvar, nz, "residuals_next_LSQR.txt")
    write_residual_array_to_file(r_svd, nvar, nz, "residuals_next_LSMR.txt")
    
    print("\nWrote next-residuals for each method to:")
    print("  residuals_next_LU.txt")
    print("  residuals_next_LSQR.txt")
    print("  residuals_next_LSMR.txt\n")


if __name__ == "__main__":
    main()
