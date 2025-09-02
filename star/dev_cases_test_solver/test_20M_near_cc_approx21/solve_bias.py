import numpy as np
from scipy.sparse import lil_matrix, csc_matrix, diags
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

    idx = 1  # start reading after the header line
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
def build_sparse_matrix_from_blocks(dblk, ublk, lblk, nvar, nz):
    """
    Build the global J matrix (nz*nvar x nz*nvar) in a sparse format
    from the diagonal, upper, and lower blocks.
    """
    size = nvar * nz
    J_lil = lil_matrix((size, size), dtype=np.float64)

    for k in range(nz):
        row_block = k * nvar
        col_block = k * nvar
        J_lil[row_block:row_block + nvar, col_block:col_block + nvar] = dblk[k]

        if k < nz - 1:
            J_lil[row_block:row_block + nvar,
                  col_block + nvar:col_block + 2*nvar] = ublk[k]

        if k > 0:
            J_lil[row_block:row_block + nvar,
                  col_block - nvar:col_block] = lblk[k]

    return J_lil.tocsc()


###############################################################################
# 3) LM-style diagonal bias + solvers
###############################################################################
def build_lm_biased_matrix(J_csc, lambda_rel=1e-3, power=1.0, min_floor=1e-16):
    """
    J' = J + lambda_rel * diag( |diag(J)|**power ) + lambda_rel*min_floor * I

    Parameters
    ----------
    J_csc : csc_matrix
        Original sparse Jacobian in CSC format.
    lambda_rel : float
        Relative strength of the diagonal bias.
    power : float
        Exponent applied to |diag(J)|; 1.0 is standard LM-like. 0.0 ~ Tikhonov I-bias.
    min_floor : float
        Ensures strictly positive bias on zero/tiny diagonals.

    Returns
    -------
    J_biased : csc_matrix
    stats : dict
    """
    if not isinstance(J_csc, csc_matrix):
        J = J_csc.tocsc()
    else:
        J = J_csc

    d = np.abs(J.diagonal())
    d_eff = d.copy()
    if power != 0.0:
        d_eff = np.power(d_eff, power)
    d_eff = d_eff + min_floor

    bias_diag = lambda_rel * d_eff
    J_biased = J + diags(bias_diag, format='csc')

    stats = {
        "lambda_rel": lambda_rel,
        "power": power,
        "diag_abs_min": float(d.min()) if d.size else 0.0,
        "diag_abs_max": float(d.max()) if d.size else 0.0,
        "bias_min": float(bias_diag.min()) if bias_diag.size else 0.0,
        "bias_max": float(bias_diag.max()) if bias_diag.size else 0.0,
    }
    return J_biased, stats


def solve_via_lu_on_matrix(J_for_factor, b):
    """
    Factorize J_for_factor with sparse LU and solve.
    Returns (x, r_J_for_factor, norm_J_for_factor, lu_obj)
    """
    lu_obj = splu(J_for_factor)
    x = lu_obj.solve(b)
    r = J_for_factor @ x - b
    residual_norm = np.linalg.norm(r)
    return x, r, residual_norm, lu_obj


def solve_via_lsqr(J, b):
    sol = lsqr(J, b, atol=1e-12, btol=1e-12, iter_lim=20000)
    x = sol[0]
    r = J @ x - b
    residual_norm = np.linalg.norm(r)
    return x, r, residual_norm


def solve_via_lsmr(J, b):
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

    nvar = 27
    nz   = 1842

    print("Parsing block files...")
    dblk = parse_dblk_file(file_dblk, nvar, nz)
    ublk = parse_ublk_file(file_ublk, nvar, nz)
    lblk = parse_lblk_file(file_lblk, nvar, nz)

    print("Parsing residuals...")
    b = parse_residuals_file(file_res, nvar, nz)

    print("Building sparse matrix J...")
    J = build_sparse_matrix_from_blocks(dblk, ublk, lblk, nvar, nz)

    print("Matrix shape:", J.shape)
    print("Non-zero entries:", J.nnz)

    # ---------------- LU with LM-style diagonal bias sweep ----------------
    print("\nSolving via LU with LM-style diagonal bias sweep...")
    lm_lambdas = [0.0, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8]
    lm_power   = 1.0

    lu_results = []
    for lam in lm_lambdas:
        if lam == 0.0:
            J_use = J
            stats = {"lambda_rel": 0.0}
        else:
            J_use, stats = build_lm_biased_matrix(J, lambda_rel=lam, power=lm_power)

        try:
            x_lu, rJuse_lu, rJuse_norm, _ = solve_via_lu_on_matrix(J_use, b)
            # Residual vs original J:
            r_orig = J @ x_lu - b
            r_orig_norm = np.linalg.norm(r_orig)

            diag_min = stats.get("diag_abs_min", np.nan)
            diag_max = stats.get("diag_abs_max", np.nan)
            bias_min = stats.get("bias_min", np.nan)
            bias_max = stats.get("bias_max", np.nan)

            print(f"  λ={lam:>9.1e} | ||J'use x - b||={rJuse_norm:8.2e} | ||Jx - b|| (orig)={r_orig_norm:8.2e} "
                  f"| diag|min,max|=({diag_min:.2e},{diag_max:.2e}) | bias[min,max]=({bias_min:.2e},{bias_max:.2e})")

            lu_results.append((f"LU(λ={lam:.1e})", x_lu, r_orig, r_orig_norm))
        except Exception as e:
            print(f"  λ={lam:>9.1e} | LU failed: {e}")

    # ---------------- LSQR (original J) ----------------
    print("\nSolving via LSQR (iterative QR-like)...")
    x_qr, r_qr, rnorm_qr = solve_via_lsqr(J, b)
    print(f"  LSQR residual norm = {rnorm_qr:e}")

    # ---------------- LSMR (original J) ----------------
    print("\nSolving via LSMR (iterative SVD-like)...")
    x_svd, r_svd, rnorm_svd = solve_via_lsmr(J, b)
    print(f"  LSMR residual norm = {rnorm_svd:e}")

    # Pick the best LU(λ) by original-J residual
    if lu_results:
        best_lu_label, best_lu_x, best_lu_r_orig, best_lu_norm = min(lu_results, key=lambda t: t[3])
        print(f"\nBest LU(λ) by original-J residual: {best_lu_label} with ||Jx-b||={best_lu_norm:e}")
        rnorm_lu_best = best_lu_norm
    else:
        best_lu_label = "LU(λ sweep failed)"
        best_lu_r_orig = None
        rnorm_lu_best = np.inf

    print("\nComparison of final residual norms (vs original J):")
    print(f"  {best_lu_label:>18s} : {rnorm_lu_best:e}")
    print(f"  {'LSQR':>18s} : {rnorm_qr:e}")
    print(f"  {'LSMR':>18s} : {rnorm_svd:e}")

    # ----------------------------------------------------------------
    #  Output the next residual arrays for each method.
    #  Use best-λ LU residual vs original J, plus LSQR/LSMR as before.
    # ----------------------------------------------------------------
    zone_of_interest = 1817
    zone_idx = zone_of_interest - 1
    if not (0 <= zone_idx < nz):
        raise IndexError(f"Requested zone {zone_of_interest} out of range [1..{nz}].")

    offset = zone_idx * nvar

    if best_lu_r_orig is not None:
        lu_res_zone = best_lu_r_orig[offset:offset + nvar]
        print(f"\nBest {best_lu_label} next-residual for zone {zone_of_interest}:")
        for i, val in enumerate(lu_res_zone, start=1):
            print(f"  var#{i:2d}: {val:24.16e}")

    qr_res_zone = r_qr[offset:offset + nvar]
    print(f"\nLSQR next-residual for zone {zone_of_interest}:")
    for i, val in enumerate(qr_res_zone, start=1):
        print(f"  var#{i:2d}: {val:24.16e}")

    svd_res_zone = r_svd[offset:offset + nvar]
    print(f"\nLSMR next-residual for zone {zone_of_interest}:")
    for i, val in enumerate(svd_res_zone, start=1):
        print(f"  var#{i:2d}: {val:24.16e}")

    # Write per-zone residual arrays
    if best_lu_r_orig is not None:
        write_residual_array_to_file(best_lu_r_orig, nvar, nz, "residuals_next_LU_LMbias.txt")
    write_residual_array_to_file(r_qr,  nvar, nz, "residuals_next_LSQR.txt")
    write_residual_array_to_file(r_svd, nvar, nz, "residuals_next_LSMR.txt")

    print("\nWrote next-residuals to:")
    if best_lu_r_orig is not None:
        print("  residuals_next_LU_LMbias.txt")
    print("  residuals_next_LSQR.txt")
    print("  residuals_next_LSMR.txt\n")


if __name__ == "__main__":
    main()

