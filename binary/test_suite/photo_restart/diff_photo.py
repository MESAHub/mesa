#!/usr/bin/env python3
"""Decode and compare MESA photos (b_*, 1_*, 2_*) field by field.

# author: Mathieu Renzo w/ Claude Opus 5.5

Photos are Fortran sequential unformatted files: every `write(iounit)` in
star/private/photo_out.f90 (stars) and binary/private/binary_photos.f90
(binary) is one record, framed by 4-byte record-length markers (gfortran).
This tool rebuilds the record layout from those write statements, decodes
each field and reports which fields differ, with values.

  photo_diff.py A B [--ignore PATTERN ...] [--keep-paths] [--max N] [--kind binary|star]

To compare pairwise photos of a binary run differing only in number of leading
zeros calling this script:

for p in 1 2 b; do for f in photos/${p}_x*; do n=${f##*_x}; n=$((10#$n)); for g in photos/${p}_x*; do m=${g##*_x}; m=$((10#$m)); ((n==m)) && [[ "$f" < "$g" ]] && python photo_diff.py $f $g && break; done; done; done

Fields holding run directories (most_recent_photo_name, *_filename) are
ignored unless --keep-paths is given.
The kind (binary or star) comes from the file name (b_ / 1_ / 2_), or --kind.
Exit status: 0 identical (after ignores), 1 different, 2 error.
Assumes gfortran, little-endian, 4-byte record markers, real(dp) = 8 bytes,
default integer and logical = 4 bytes. Records whose length does not match
the layout (code changed) are compared as raw bytes and flagged.
"""

import argparse
import fnmatch
import os
import struct
import sys

STRLEN = 256  # const_def: strlen; also net_name_len
HISTLEN = 80  # maxlen_history_column_name, maxlen_binary_history_column_name
RAND_U_LEN = 97  # star_data_def: rand_u_len
XTRA = 30  # star_num_xtra_vals, binary_num_xtra_vals
SIZE = {"r": 8, "i": 4, "l": 4, "c": STRLEN, "h": HISTLEN}
PATH_FIELDS = ["most_recent_photo_name", "*_filename"]  # contain run directories


# ---------------------------------------------------------------- records
def read_records(path):
    with open(path, "rb") as f:
        data = f.read()
    recs, pos = [], 0
    while pos < len(data):
        chunk = b""
        while True:
            (m,) = struct.unpack_from("<i", data, pos)
            n = abs(m)
            chunk += data[pos + 4 : pos + 4 + n]
            (tail,) = struct.unpack_from("<i", data, pos + 4 + n)
            if abs(tail) != n:
                raise ValueError(f"{path}: bad record marker at byte {pos}")
            pos += 8 + n
            if m >= 0:  # negative leading marker: record continues (subrecords)
                break
        recs.append(chunk)
    return recs


def unpack(t, buf, off, count):
    size = SIZE[t]
    if t == "r":
        vals = list(struct.unpack_from(f"<{count}d", buf, off))
    elif t in "il":
        vals = list(struct.unpack_from(f"<{count}i", buf, off))
        if t == "l":
            vals = [v != 0 for v in vals]
    else:
        vals = [
            buf[off + k * size : off + (k + 1) * size].decode("latin-1").rstrip()
            for k in range(count)
        ]
    return vals, off + size * count


def decode_record(buf, spec, ctx):
    """spec: [(name, type, count)]; count int, ctx key, or '?' (from length)."""
    counts = []
    for _, t, c in spec:
        counts.append(None if c == "?" else (ctx[c] if isinstance(c, str) else c))
    fixed = sum(SIZE[t] * n for (_, t, _), n in zip(spec, counts) if n is not None)
    if None in counts:  # unknown sizes, all equal (e.g. numTlim): from the length
        unit = sum(SIZE[t] for (_, t, _), n in zip(spec, counts) if n is None)
        rest = len(buf) - fixed
        if rest < 0 or rest % unit:
            return None
        counts = [rest // unit if n is None else n for n in counts]
    elif fixed != len(buf):
        return None
    fields, off = {}, 0
    for (name, t, _), n in zip(spec, counts):
        vals, off = unpack(t, buf, off, n)
        key = name
        while key in fields:  # same variable written twice
            key += "'"
        fields[key] = vals
    return fields


# ---------------------------------------------------------------- layouts
def S(names, t):
    return [(n, t, 1) for n in names.split()]


def star_layout():
    """Generator: yields record specs, receives decoded fields (ctx updated)."""
    part = [("part_number", "i", 1)]
    yield "version", [("star_def_version", "i", 1)]
    yield "part", part
    yield (
        "header",
        (
            S("initial_z", "r")
            + S("total_num_solver_iterations nz nvar_hydro nvar_chem nvar_total", "i")
            + S(
                "v_flag u_flag rotation_flag RSP2_flag RSP_flag RTI_flag w_div_wc_flag "
                "j_rot_flag D_omega_flag am_nu_rot_flag have_mlt_vc",
                "l",
            )
            + S("species num_reactions model_number", "i")
            + S(
                "star_mass mstar xmstar M_center v_center R_center L_center time dt",
                "r",
            )
            + S(
                "have_previous_conv_vel was_in_implicit_wind_limit using_revised_net_name",
                "l",
            )
            + S("revised_net_name", "c")
            + S("using_revised_max_yr_dt", "l")
            + S("revised_max_yr_dt", "r")
            + S("astero_using_revised_max_yr_dt", "l")
            + S(
                "astero_revised_max_yr_dt cumulative_energy_error cumulative_extra_heating",
                "r",
            )
            + S("have_initial_energy_integrals", "l")
            + S(
                "total_energy_initial force_tau_factor force_opacity_factor crystal_core_boundary_mass",
                "r",
            )
        ),
    )
    yield "net_name", S("net_name", "c")
    yield "part", part
    yield (
        "profiles",
        (
            [("dq", "r", "nz"), ("xa", "r", "species*nz"), ("xh", "r", "nvar_hydro*nz")]
            + [
                (n, "r", "nz")
                for n in [
                    "omega",
                    "j_rot",
                    "mlt_vc",
                    "conv_vel",
                    "D_ST_start",
                    "nu_ST_start",
                ]
            ]
            + S("have_ST_start_info", "l")
            + S("mstar_old", "r")
        ),
    )
    yield "part", part
    yield (
        "rsp",
        (
            S("rsp_num_periods", "i")
            + S("rsp_dt rsp_period", "r")
            + S("RSP_have_set_velocities", "l")
            + S("dt_limit_ratio tau_base", "r")
        ),
    )
    yield (
        "var_indices",
        S(
            "i_lnd i_lnT i_lnR i_lum i_Et_RSP i_erad_RSP i_Fr_RSP i_v i_u i_alpha_RTI i_w i_Hp "
            "i_w_div_wc i_j_rot i_dv_dt i_equL i_dlnd_dt i_dlnE_dt i_dEt_RSP_dt i_derad_RSP_dt "
            "i_dFr_RSP_dt i_du_dt i_dlnR_dt i_dalpha_RTI_dt i_detrb_dt i_equ_Hp",
            "i",
        ),
    )
    yield (
        "filenames",
        S(
            "model_controls_filename model_data_filename most_recent_profile_filename "
            "most_recent_controls_filename most_recent_model_data_filename",
            "c",
        ),
    )
    yield "part", part
    yield (
        "timestep",
        S("recent_log_header phase_of_evolution", "i")
        + S("dt_next dt_next_unclipped", "r"),
    )
    yield "part", part
    yield (
        "counters",
        (
            S(
                "num_skipped_setvars num_retries num_setvars total_num_solver_iterations "
                "total_num_solver_relax_iterations total_num_solver_calls_made "
                "total_num_solver_relax_calls_made total_num_solver_calls_converged "
                "total_num_solver_relax_calls_converged total_step_attempts total_relax_step_attempts "
                "total_step_retries total_relax_step_retries total_step_redos total_relax_step_redos "
                "total_steps_finished total_relax_steps_finished num_hydro_merges num_hydro_splits "
                "num_solver_setvars mesh_call_number solver_call_number diffusion_call_number",
                "i",
            )
            + S(
                "gradT_excess_alpha gradT_excess_alpha_old Teff mstar_dot power_nuc_burn power_h_burn "
                "power_he_burn power_z_burn power_photo",
                "r",
            )
            + S("why_Tlim", "i")
            + [("dt_why_count", "i", "?"), ("dt_why_retry_count", "i", "?")]
            + S(
                "timestep_hold model_number_for_last_retry model_number_for_last_retry_old init_model_number",
                "i",
            )
            + S("most_recent_photo_name", "c")
            + S("rand_i97 rand_j97", "i")
            + [("rand_u", "r", RAND_U_LEN)]
            + S("rand_c rand_cd rand_cm", "r")
        ),
    )
    yield "part", part
    ctx = yield "extra_lens", S("len_extra_iwork len_extra_work", "i")
    if ctx["len_extra_iwork"] > 0:
        yield "extra_iwork", [("extra_iwork", "i", "len_extra_iwork")]
    if ctx["len_extra_work"] > 0:
        yield "extra_work", [("extra_work", "r", "len_extra_work")]
    yield from history_layout("", "xtra", XTRA)


def history_layout(prefix, xtra, nxtra):
    if xtra:
        yield "ixtra", [("ixtra", "i", nxtra)]
        yield "xtra", [("xtra", "r", nxtra)]
        yield "lxtra", [("lxtra", "l", nxtra)]
    ctx = yield "history_spec_len", [("len_history_col_spec", "i", 1)]
    if ctx["len_history_col_spec"] > 0:
        yield (
            "history_spec",
            [(f"{prefix}history_column_spec", "i", "len_history_col_spec")],
        )
    ctx = yield (
        "history_info",
        S(
            f"number_of_{prefix}history_columns model_number_of_{prefix}history_values",
            "i",
        )
        + S(f"need_to_set_{prefix}history_names_etc", "l"),
    )
    n = ctx[f"number_of_{prefix}history_columns"]
    if n > 0:
        ctx = yield (
            "history_is_integer",
            [(f"{prefix}history_value_is_integer", "l", n)],
        )
        for k in range(1, n + 1):
            ctx = yield "history_name", [(f"{prefix}history_names({k})", "h", 1)]
    return ctx
    # then: rsp_photo_out (if RSP_flag), part numbers, other_photo_write: raw


CE_STATE = (
    S("CE_flag CE_flag_old CE_init CE_init_old", "l")
    + S("CE_nz", "i")
    + S(
        "CE_initial_radius CE_initial_separation CE_initial_Mdonor CE_initial_Maccretor "
        "CE_initial_age CE_initial_model_number CE_b_initial_age CE_b_initial_model_number",
        "r",
    )
    + S("CE_num1 CE_num1_old CE_num2 CE_num2_old", "i")
    + S(
        "CE_lambda1 CE_lambda1_old CE_lambda2 CE_lambda2_old CE_Ebind1 CE_Ebind1_old "
        "CE_Ebind2 CE_Ebind2_old CE_years_detached CE_years_detached_old",
        "r",
    )
)


def binary_layout(with_ce):
    """with_ce: photos with the CE variables (main); without (mr/fix_bin_photos)."""
    yield "version", [("star_def_version", "i", 1)]

    def pairs(names, t):
        return [(m, t, 1) for n in names.split() for m in (n, n + "_old")]

    yield (
        "state",
        (
            pairs("binary_age", "r")
            + pairs("model_number", "i")
            + pairs("mtransfer_rate angular_momentum_j separation eccentricity", "r")
            + [
                ("rl_relative_gap(1)", "r", 1),
                ("rl_relative_gap_old(1)", "r", 1),
                ("rl_relative_gap(2)", "r", 1),
                ("rl_relative_gap_old(2)", "r", 1),
            ]
            + [
                (f"{n}{o}({k})", "r", 1)
                for n in ("r", "rl", "m")
                for k in (1, 2)
                for o in ("", "_old")
            ]
            + pairs("dt", "r")
            + [
                ("env(1)", "r", 1),
                ("env_old(1)", "r", 1),
                ("env(2)", "r", 1),
                ("env_old(2)", "r", 1),
            ]
            + S("eq_initial_bh_mass", "r")
            + pairs("period max_timestep change_factor", "r")
            + S("min_binary_separation", "r")
            + [
                ("using_jdot_mb(1)", "l", 1),
                ("using_jdot_mb_old(1)", "l", 1),
                ("using_jdot_mb(2)", "l", 1),
                ("using_jdot_mb_old(2)", "l", 1),
            ]
            + pairs("d_i a_i point_mass_i", "i")
            + pairs("ignore_rlof_flag model_twins_flag", "l")
            + pairs("dt_why_reason", "i")
            + S("have_star_1 have_star_2", "l")
            + (CE_STATE if with_ce else [])
            + S("generations", "i")
            + [
                ("ixtra", "i", XTRA),
                ("ixtra_old", "i", XTRA),
                ("xtra", "r", XTRA),
                ("xtra_old", "r", XTRA),
                ("lxtra", "l", XTRA),
                ("lxtra_old", "l", XTRA),
            ]
        ),
    )
    ctx = yield from history_layout("binary_", None, 0)
    if with_ce and ctx["CE_init"]:
        yield (
            "CE_profiles",
            [("CE_m", "r", "CE_nz")]
            + [
                (n, "r", "4*CE_nz")
                for n in [
                    "CE_entropy",
                    "CE_U_in",
                    "CE_U_out",
                    "CE_Omega_in",
                    "CE_Omega_out",
                ]
            ],
        )
    # then: other_binary_photo_write: raw


def decode(path, kind):
    recs = read_records(path)
    if kind == "star":
        return decode_with(recs, star_layout())
    # binary: try the photo layout without and with the CE variables
    for with_ce in (False, True):
        out = decode_with(recs, binary_layout(with_ce))
        if len(out) > 1 and out[1][2] is not None:
            return out
    return out


def decode_with(recs, gen):
    out, ctx, i = [], {}, 0
    try:
        label, spec = next(gen)
        while i < len(recs):
            fields = decode_record(recs[i], spec, ctx)
            if fields is None:
                out.append((i + 1, label + " [LAYOUT MISMATCH: raw]", None, recs[i]))
                i += 1
                break
            for k, v in fields.items():
                if len(v) == 1:
                    ctx[k] = v[0]
            if "nz" in ctx and "species" in ctx:
                ctx["species*nz"] = ctx["species"] * ctx["nz"]
            if "CE_nz" in ctx:
                ctx["4*CE_nz"] = 4 * ctx["CE_nz"]
            if "nz" in ctx and "nvar_hydro" in ctx:
                ctx["nvar_hydro*nz"] = ctx["nvar_hydro"] * ctx["nz"]
            out.append((i + 1, label, fields, recs[i]))
            i += 1
            label, spec = gen.send(ctx)
    except StopIteration:
        pass
    for j in range(i, len(recs)):  # rsp, part numbers, other_*_photo_write
        out.append(
            (j + 1, "tail (rsp / part number / other_photo_write)", None, recs[j])
        )
    return out


def kind_of(path, kind):
    if kind:
        return kind
    base = os.path.basename(path)
    if base.startswith("b_"):
        return "binary"
    if base[:2] in ("1_", "2_") or base.startswith("x"):
        return "star"
    raise ValueError(f"cannot tell the photo kind of {path}; use --kind binary|star")


# ---------------------------------------------------------------- output
def fmt(v):
    return repr(v)


# 2D arrays written as x(:,1:nz): name -> scalar holding the first dimension
FIRST_DIM = {"xa": "species", "xh": "nvar_hydro"}


def same(x, y):
    """bit-for-bit for reals (so -0.0 != 0.0 and NaN == NaN with the same bits)"""
    if isinstance(x, float):
        return struct.pack("<d", x) == struct.pack("<d", y)
    return x == y


def where(name, k, scalars):
    n1 = scalars.get(FIRST_DIM.get(name))
    return f"({k % n1 + 1},{k // n1 + 1})" if n1 else f"({k + 1})"


def diff_values(name, a, b, scalars):
    if len(a) != len(b):
        return f"{name}: length {len(a)} vs {len(b)}"
    bad = [k for k in range(len(a)) if not same(a[k], b[k])]
    if not bad:
        return None
    if len(a) == 1:
        extra = ""
        if isinstance(a[0], float) and a[0] and b[0] is not None:
            extra = (
                f"  (rel {abs(b[0] - a[0]) / abs(a[0]):.3e})" if a[0] == a[0] else ""
            )
        return f"{name}: {fmt(a[0])} -> {fmt(b[0])}{extra}"
    k = bad[0]
    msg = (
        f"{name}: {len(bad)}/{len(a)} elements differ, first at {name}{where(name, k, scalars)}: "
        f"{fmt(a[k])} -> {fmt(b[k])}"
    )
    if isinstance(a[0], float):
        rel = max(abs(b[j] - a[j]) / abs(a[j]) if a[j] else abs(b[j]) for j in bad)
        msg += f", max rel {rel:.3e}"
    return msg


def cmd_diff(args):
    ka = kind_of(args.a, args.kind)
    da, db = decode(args.a, ka), decode(args.b, kind_of(args.b, args.kind))
    ignore = list(args.ignore) + ([] if args.keep_paths else PATH_FIELDS)
    scalars = {k: v[0] for _, _, f, _ in da if f for k, v in f.items() if len(v) == 1}
    lines = []
    for ra, rb in zip(da, db):
        (ia, label, fa, bufa), (_, labelb, fb, bufb) = ra, rb
        if fa is not None and fb is not None and label == labelb:
            for name in fa:
                if any(fnmatch.fnmatchcase(name, p) for p in ignore):
                    continue
                d = diff_values(name, fa[name], fb.get(name, []), scalars)
                if d:
                    lines.append(f"  record {ia:3d} {label:14s} {d}")
        elif bufa != bufb:
            first = next(
                (k for k in range(min(len(bufa), len(bufb))) if bufa[k] != bufb[k]),
                min(len(bufa), len(bufb)),
            )
            lines.append(
                f"  record {ia:3d} {label}: raw bytes differ (lengths {len(bufa)}/{len(bufb)}, "
                f"first differing byte {first})"
            )
    if len(da) != len(db):
        lines.append(f"  number of records: {len(da)} vs {len(db)}")
    if not lines:
        print(
            f"{args.a} == {args.b}"
            + (" (ignoring " + ", ".join(ignore) + ")" if ignore else "")
        )
        return 0
    print(f"{args.a} != {args.b}: {len(lines)} differences")
    for line in lines[: args.max]:
        print(line)
    if len(lines) > args.max:
        print(f"  ... {len(lines) - args.max} more (--max)")
    return 1


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("a", help="first photo")
    p.add_argument("b", help="second photo")
    p.add_argument(
        "--ignore",
        action="append",
        default=[],
        metavar="PATTERN",
        help="field name pattern to skip (fnmatch), repeatable",
    )
    p.add_argument(
        "--keep-paths",
        action="store_true",
        help="also compare fields holding run directories: " + ", ".join(PATH_FIELDS),
    )
    p.add_argument("--max", type=int, default=50, help="max differences to print")
    p.add_argument("--kind", choices=("binary", "star"))
    args = p.parse_args()
    try:
        return cmd_diff(args)
    except (OSError, ValueError, struct.error) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
