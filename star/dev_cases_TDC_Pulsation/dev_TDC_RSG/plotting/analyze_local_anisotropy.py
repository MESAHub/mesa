#!/usr/bin/env python3
"""Diagnose the local Kovacs et al. anisotropy closure on saved RSG profiles.

This is post-processing only.  It does not compile or run MESA and it does
not modify a model.  The published local closure is evaluated with its
nonlocal turbulent-energy diffusion term D_omega set to zero.  It is still
time dependent because every quantity is evaluated from the instantaneous
TDC and hydrodynamic state in each saved profile.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mesa_reader as mr
import numpy as np


CASE = Path(__file__).resolve().parents[1]
LOGS = CASE / "LOGS"
OUT = Path(__file__).resolve().parent

RSUN = 6.957e10
LSUN = 3.828e33
MSUN = 1.98847e33
DAY = 86400.0
SQRT_2_DIV_3 = np.sqrt(2.0 / 3.0)

# Canuto et al. coefficients adopted by Kovacs, Szabo & Nuspl (2026).
ALPHA1 = 1.08
BETA5 = 0.7

# This run changed alpha_M without restarting LOGS.  Keep the reconstruction
# of the current eddy stress consistent with the control active in each model.
# The anisotropy closure itself and the saved profile Eq do not use this table.
TDC_ALPHA_M_EPOCHS = (
    (0, 0.5),
    (51000, 1.0),
    (53000, 0.7),
)


def profile_number(path: Path) -> int:
    match = re.search(r"profile(\d+)\.data$", path.name)
    if match is None:
        raise ValueError(f"cannot read profile number from {path}")
    return int(match.group(1))


def read_real_control(name: str) -> float:
    """Read the last active assignment of a real control in inlist_pulses."""
    pattern = re.compile(rf"^\s*{re.escape(name)}\s*=\s*([^!\s]+)", re.I)
    value = None
    for line in (CASE / "inlist_pulses").read_text().splitlines():
        if line.lstrip().startswith("!"):
            continue
        match = pattern.match(line)
        if match is not None:
            value = float(re.sub(r"[dD]", "e", match.group(1)))
    if value is None:
        raise RuntimeError(f"missing {name} in inlist_pulses")
    return value


def face_average(values: np.ndarray, dq: np.ndarray) -> np.ndarray:
    """Reproduce the case's TDC mass interpolation from cells to faces."""
    result = np.asarray(values, dtype=float).copy()
    if result.size < 2:
        return result
    # MESA profiles run from the surface (k=1) inward.  At face k>1,
    # value = alfa*cell(k) + beta*cell(k-1), with
    # alfa = dq(k-1)/(dq(k-1)+dq(k)).
    denom = dq[:-1] + dq[1:]
    alfa = np.divide(dq[:-1], denom, out=np.full_like(denom, 0.5), where=denom > 0)
    result[1:] = alfa * values[1:] + (1.0 - alfa) * values[:-1]
    return result


def radial_derivative(values: np.ndarray, radius: np.ndarray) -> np.ndarray:
    """Differentiate on an increasing-radius copy, then restore MESA order."""
    return np.gradient(values[::-1], radius[::-1], edge_order=2)[::-1]


def weighted_quantile(
    values: np.ndarray, weights: np.ndarray, quantiles: tuple[float, ...]
) -> np.ndarray:
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(good):
        return np.full(len(quantiles), np.nan)
    x = values[good]
    w = weights[good]
    order = np.argsort(x)
    x = x[order]
    w = w[order]
    cdf = (np.cumsum(w) - 0.5 * w) / np.sum(w)
    return np.interp(quantiles, cdf, x, left=x[0], right=x[-1])


def get_header_value(profile: mr.MesaData, name: str, default: float = np.nan) -> float:
    try:
        value = getattr(profile, name)
    except (AttributeError, KeyError):
        return default
    array = np.asarray(value)
    return float(array.ravel()[0])


def alpha_m_for_model(model_number: int) -> float:
    value = TDC_ALPHA_M_EPOCHS[0][1]
    for first_model, epoch_value in TDC_ALPHA_M_EPOCHS:
        if model_number < first_model:
            break
        value = epoch_value
    return value


def compute_profile(path: Path, controls: dict[str, float]) -> dict[str, object]:
    p = mr.MesaData(str(path))
    model_number = int(round(get_header_value(p, "model_number")))
    alpha_m = alpha_m_for_model(model_number)
    r = np.asarray(p.radius) * RSUN
    v = np.asarray(p.velocity)
    dq = np.asarray(p.dq)
    rho = face_average(10.0 ** np.asarray(p.logRho), dq)
    hp = face_average(np.asarray(p.pressure_scale_height) * RSUN, dq)
    grada = face_average(np.asarray(p.grada), dq)

    dvdr = radial_derivative(v, r)
    v_div_r = np.divide(v, r, out=np.zeros_like(v), where=r > 0)

    conv_vel = np.asarray(p.mlt_vc)
    amplitude = conv_vel / SQRT_2_DIV_3
    omega = amplitude**2

    lambda0 = controls["mixing_length_alpha"] * hp
    beta_r = controls["harmonic_dissipation_length_beta"] * r
    if controls["harmonic_dissipation_length_beta"] > 0:
        mixing_length = np.divide(
            1.0,
            1.0 / np.maximum(lambda0, 1e-99) + 1.0 / np.maximum(beta_r, 1e-99),
        )
    else:
        mixing_length = lambda0

    alpha_d = (8.0 / 3.0) * SQRT_2_DIV_3 * controls["TDC_alpha_D"]
    epsilon = alpha_d * omega**1.5 / np.maximum(mixing_length, 1e-99)

    # The TDC source and convective luminosity contain the same A*Y*Cp*T
    # factor.  Their ratio therefore removes Cp and Y:
    # S = (alpha_S/alpha_C) grada Lconv/(4*pi*r^2*rho*Hp).
    lum_conv = np.asarray(p.lum_conv) * LSUN
    source = (
        controls["TDC_alpha_S"]
        / controls["TDC_alpha_C"]
        * grada
        * lum_conv
        / np.maximum(4.0 * np.pi * r**2 * rho * hp, 1e-99)
    )

    lambda_xi = (
        (2.0 / 3.0) * (ALPHA1 * v_div_r + (2.0 * ALPHA1 - 1.0) * dvdr)
        - 2.5 * np.divide(epsilon, omega, out=np.zeros_like(omega), where=omega > 0)
    )
    b_term = -(2.0 * BETA5 + 1.0) * source + epsilon  # D_omega = 0.
    c_term = (ALPHA1 - 0.2) * dvdr - 2.0 * v_div_r

    active_floor = max(1.0, 1e-3 * float(np.nanmax(amplitude)))
    active = (
        np.isfinite(omega)
        & np.isfinite(lambda_xi)
        & np.isfinite(source)
        & (amplitude > active_floor)
        & (omega > 0)
    )

    xi_raw = np.full_like(omega, np.nan)
    safe = active & (lambda_xi != 0)
    xi_raw[safe] = (
        1.0 / 3.0
        + b_term[safe] / (3.0 * omega[safe] * lambda_xi[safe])
        - 2.0 * c_term[safe] / (9.0 * lambda_xi[safe])
    )

    # This hard clip is deliberately diagnostic only.  A coupled Newton solve
    # would need a smooth realizability limiter or a bounded stress variable.
    xi_clip = np.clip(xi_raw, 0.0, 1.0)
    q_aniso_raw = 2.0 * omega * (xi_raw - 1.0 / 3.0)
    e_aniso_raw = q_aniso_raw * (dvdr - 2.0 * v_div_r)
    q_aniso = 2.0 * omega * (xi_clip - 1.0 / 3.0)
    e_aniso = q_aniso * (dvdr - 2.0 * v_div_r)

    # Continuum estimate of the current Kuhfuss eddy-viscous stress.  The
    # profile Eq column remains the authoritative discrete MESA work term.
    shear = dvdr - v_div_r
    q_kuhfuss = (
        -(4.0 / 3.0)
        * alpha_m
        * mixing_length
        * amplitude
        * shear
    )
    e_kuhfuss_cont = -q_kuhfuss * shear
    eq_mesa = np.asarray(p.Eq)

    turnover_rate = np.divide(
        amplitude, mixing_length, out=np.zeros_like(amplitude), where=mixing_length > 0
    )
    scale_rate = turnover_rate + np.abs(dvdr) + np.abs(v_div_r)
    conditioning = np.divide(
        np.abs(lambda_xi), scale_rate, out=np.full_like(lambda_xi, np.nan), where=scale_rate > 0
    )
    tau_xi_days = np.divide(
        1.0, np.abs(lambda_xi) * DAY, out=np.full_like(lambda_xi, np.nan), where=lambda_xi != 0
    )

    star_mass = get_header_value(p, "star_mass")
    dm = dq * star_mass * MSUN
    tke_weight = dm * omega

    return {
        "path": path,
        "model": model_number,
        "TDC_alpha_M": alpha_m,
        "time_days": get_header_value(p, "time_seconds") / DAY,
        "logT": np.asarray(p.logT),
        "radius_rsun": np.asarray(p.radius),
        "mass": np.asarray(p.mass),
        "dm": dm,
        "active": active,
        "omega": omega,
        "source": source,
        "epsilon": epsilon,
        "lambda_xi": lambda_xi,
        "tau_xi_days": tau_xi_days,
        "conditioning": conditioning,
        "xi_raw": xi_raw,
        "xi_clip": xi_clip,
        "q_aniso_raw": q_aniso_raw,
        "e_aniso_raw": e_aniso_raw,
        "q_aniso": q_aniso,
        "e_aniso": e_aniso,
        "q_kuhfuss": q_kuhfuss,
        "e_kuhfuss_cont": e_kuhfuss_cont,
        "eq_mesa": eq_mesa,
        "tke_weight": tke_weight,
    }


def summarize(data: dict[str, object], period_days: float) -> dict[str, float]:
    active = np.asarray(data["active"], dtype=bool)
    xi_raw = np.asarray(data["xi_raw"])
    xi_clip = np.asarray(data["xi_clip"])
    dm = np.asarray(data["dm"])
    tke_weight = np.asarray(data["tke_weight"])
    q_aniso = np.asarray(data["q_aniso"])
    e_aniso_raw = np.asarray(data["e_aniso_raw"])
    q_kuhfuss = np.asarray(data["q_kuhfuss"])
    e_aniso = np.asarray(data["e_aniso"])
    e_kuhfuss_cont = np.asarray(data["e_kuhfuss_cont"])
    eq_mesa = np.asarray(data["eq_mesa"])
    conditioning = np.asarray(data["conditioning"])
    tau = np.asarray(data["tau_xi_days"])

    finite = active & np.isfinite(xi_raw)
    weights = np.where(finite, tke_weight, 0.0)
    xi_q = weighted_quantile(xi_raw[finite], weights[finite], (0.05, 0.5, 0.95))
    clipped_mean = np.sum(weights[finite] * xi_clip[finite]) / max(np.sum(weights[finite]), 1e-99)
    oob = finite & ((xi_raw < 0.0) | (xi_raw > 1.0))
    near_singular = finite & (conditioning < 1e-2)
    backscatter = finite & (e_aniso < 0.0)

    p_aniso = float(np.sum(dm[finite] * e_aniso[finite]))
    p_aniso_raw = float(np.sum(dm[finite] * e_aniso_raw[finite]))
    p_eq = float(np.sum(dm[finite] * eq_mesa[finite]))
    p_cont = float(np.sum(dm[finite] * e_kuhfuss_cont[finite]))

    stress_floor = 1e-6 * float(np.nanmax(np.abs(q_kuhfuss[finite]))) if np.any(finite) else np.inf
    stress_good = finite & (np.abs(q_kuhfuss) > stress_floor)
    stress_ratio = weighted_quantile(
        np.abs(q_aniso[stress_good] / q_kuhfuss[stress_good]),
        tke_weight[stress_good],
        (0.5,),
    )[0]
    tau_ratio = weighted_quantile(
        tau[finite] / period_days, weights[finite], (0.5,)
    )[0]

    return {
        "model_number": float(data["model"]),
        "TDC_alpha_M": float(data["TDC_alpha_M"]),
        "time_days": float(data["time_days"]),
        "active_zones": float(np.count_nonzero(finite)),
        "xi_raw_p05": xi_q[0],
        "xi_raw_median": xi_q[1],
        "xi_raw_p95": xi_q[2],
        "xi_clipped_tke_mean": clipped_mean,
        "raw_out_of_bounds_fraction": np.count_nonzero(oob) / max(np.count_nonzero(finite), 1),
        "near_singular_fraction": np.count_nonzero(near_singular) / max(np.count_nonzero(finite), 1),
        "backscatter_zone_fraction": np.count_nonzero(backscatter) / max(np.count_nonzero(finite), 1),
        "anisotropy_power_erg_s": p_aniso,
        "raw_anisotropy_power_erg_s": p_aniso_raw,
        "mesa_Eq_power_erg_s": p_eq,
        "continuum_eddy_power_erg_s": p_cont,
        "anisotropy_to_Eq_power": p_aniso / p_eq if p_eq != 0 else np.nan,
        "raw_anisotropy_to_Eq_power": p_aniso_raw / p_eq if p_eq != 0 else np.nan,
        "continuum_to_Eq_power": p_cont / p_eq if p_eq != 0 else np.nan,
        "median_abs_stress_ratio": stress_ratio,
        "median_tau_xi_div_period": tau_ratio,
    }


def read_period_series() -> tuple[np.ndarray, np.ndarray]:
    # Reuse the restart-safe history reader maintained for this test case.
    from analyze_pulsation import cycle_end_indices, read_history_tail

    history = read_history_tail(LOGS / "history.data")
    ends = cycle_end_indices(history)
    models = history["model_number"][ends]
    periods = history["period"][ends]
    good = np.isfinite(models) & np.isfinite(periods) & (periods > 0)
    models = models[good]
    periods = periods[good]
    if periods.size == 0:
        raise RuntimeError("no completed-cycle periods in history.data")
    return models, periods


def period_at_model(
    model_number: int, period_models: np.ndarray, periods: np.ndarray
) -> float:
    return float(np.interp(model_number, period_models, periods))


def plot_latest(data: dict[str, object], period_days: float) -> None:
    active = np.asarray(data["active"], dtype=bool)
    logt = np.asarray(data["logT"])
    finite = active & np.isfinite(np.asarray(data["xi_raw"]))
    oob = finite & (
        (np.asarray(data["xi_raw"]) < 0.0) | (np.asarray(data["xi_raw"]) > 1.0)
    )

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)

    ax = axes[0, 0]
    ax.plot(logt[finite], np.asarray(data["xi_clip"])[finite], lw=1.5, label=r"$\xi$ (diagnostic clip)")
    if np.any(oob):
        raw = np.asarray(data["xi_raw"])
        ax.scatter(logt[oob], np.clip(raw[oob], 0.0, 1.0), s=10, color="tab:red", label="raw outside [0,1]")
    ax.axhline(1.0 / 3.0, color="0.35", ls="--", lw=1, label="isotropic")
    ax.set_ylim(-0.05, 1.05)
    ax.set_ylabel(r"radial TKE fraction $\xi$")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[0, 1]
    ax.plot(logt[finite], np.sqrt(3.0 * np.asarray(data["xi_clip"])[finite]), label=r"$v_{r,\rm rms}/v_{r,\rm iso}$")
    ax.axhline(1.0, color="0.35", ls="--", lw=1)
    ax.set_ylabel("radial rms velocity factor")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    ax.plot(logt[finite], np.asarray(data["q_aniso"])[finite], label="local anisotropy")
    ax.plot(logt[finite], np.asarray(data["q_kuhfuss"])[finite], label="current eddy stress (continuum)")
    qscale = np.nanpercentile(
        np.abs(np.r_[np.asarray(data["q_aniso"])[finite], np.asarray(data["q_kuhfuss"])[finite]]),
        15,
    )
    ax.set_yscale("symlog", linthresh=max(qscale, 1.0))
    ax.axhline(0, color="0.35", lw=0.8)
    ax.set_xlabel(r"$\log T$")
    ax.set_ylabel(r"specific stress $q$ (cm$^2$ s$^{-2}$)")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1, 1]
    ax.plot(logt[finite], np.asarray(data["e_aniso"])[finite], label=r"anisotropy $E_\nu$")
    ax.plot(logt[finite], np.asarray(data["eq_mesa"])[finite], label="MESA Eq")
    escale = np.nanpercentile(
        np.abs(np.r_[np.asarray(data["e_aniso"])[finite], np.asarray(data["eq_mesa"])[finite]]),
        15,
    )
    ax.set_yscale("symlog", linthresh=max(escale, 1e-12))
    ax.axhline(0, color="0.35", lw=0.8)
    ax.set_xlabel(r"$\log T$")
    ax.set_ylabel(r"TKE transfer (erg g$^{-1}$ s$^{-1}$)")
    ax.legend(frameon=False, fontsize=9)

    fig.suptitle(
        f"Local algebraic anisotropy on RSG model {data['model']}\n"
        rf"$D_\omega=0$; instantaneous TDC state; mode period $\simeq${period_days:.1f} d"
    )
    fig.savefig(OUT / "local_anisotropy_profiles.png", dpi=180)
    plt.close(fig)


def plot_evolution(summaries: list[dict[str, float]]) -> None:
    model = np.asarray([row["model_number"] for row in summaries])
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)

    ax = axes[0, 0]
    ax.plot(model, [row["xi_clipped_tke_mean"] for row in summaries], "o-", ms=3)
    ax.axhline(1.0 / 3.0, color="0.35", ls="--", lw=1)
    ax.set_ylabel(r"TKE-weighted $\xi$ (clipped)")

    ax = axes[0, 1]
    ax.plot(model, [row["raw_out_of_bounds_fraction"] for row in summaries], "o-", ms=3, label="raw outside [0,1]")
    ax.plot(model, [row["near_singular_fraction"] for row in summaries], "o-", ms=3, label=r"$|\lambda_\xi|/$local rate $<10^{-2}$")
    ax.set_ylabel("active-zone fraction")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1, 0]
    ax.plot(model, [row["anisotropy_to_Eq_power"] for row in summaries], "o-", ms=3)
    ax.axhline(0, color="0.35", lw=0.8)
    ax.axhline(1, color="0.35", ls="--", lw=0.8)
    ax.set_xlabel("model number")
    ax.set_ylabel(r"$\int E_\nu dm / \int Eq\,dm$")

    ax = axes[1, 1]
    ax.plot(model, [row["backscatter_zone_fraction"] for row in summaries], "o-", ms=3, label=r"$E_\nu<0$")
    ax.plot(model, [row["median_tau_xi_div_period"] for row in summaries], "o-", ms=3, label=r"median $|\tau_\xi|/P$")
    ax.set_xlabel("model number")
    ax.set_ylabel("fraction or timescale ratio")
    ax.set_yscale("log")
    ax.legend(frameon=False, fontsize=9)

    for ax in axes.flat:
        ax.axvline(51000, color="tab:red", ls=":", lw=1)
        ax.axvline(53000, color="tab:red", ls=":", lw=1)

    fig.suptitle("Saved-profile variation of the local algebraic anisotropy estimate")
    fig.savefig(OUT / "local_anisotropy_evolution.png", dpi=180)
    plt.close(fig)


def plot_fixed_epoch(summaries: list[dict[str, float]]) -> None:
    rows = [row for row in summaries if 41000 <= row["model_number"] <= 50000]
    if not rows:
        return
    model = np.asarray([row["model_number"] for row in rows])
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)

    ax = axes[0, 0]
    ax.plot(model, [row["xi_clipped_tke_mean"] for row in rows], "o-")
    ax.axhline(1.0 / 3.0, color="0.35", ls="--", lw=1)
    ax.set_ylabel(r"TKE-weighted $\xi$ (clipped)")

    ax = axes[0, 1]
    ax.plot(model, [row["anisotropy_to_Eq_power"] for row in rows], "o-", label="diagnostic clip")
    ax.plot(model, [row["raw_anisotropy_to_Eq_power"] for row in rows], "--", label="raw")
    ax.axhline(0, color="0.35", lw=0.8)
    ax.set_ylabel(r"$\int E_\nu dm / \int Eq\,dm$")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    ax.plot(model, [row["backscatter_zone_fraction"] for row in rows], "o-", label=r"$E_\nu<0$")
    ax.plot(model, [row["raw_out_of_bounds_fraction"] for row in rows], "o-", label=r"raw $\xi\notin[0,1]$")
    ax.set_xlabel("model number")
    ax.set_ylabel("active-zone fraction")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    ax.plot(model, [row["median_abs_stress_ratio"] for row in rows], "o-", label=r"median $|q_\nu/q_{\rm EV}|$")
    ax.plot(model, [row["continuum_to_Eq_power"] for row in rows], "o-", label="continuum work / Eq")
    ax.axhline(1, color="0.35", ls="--", lw=0.8)
    ax.set_xlabel("model number")
    ax.set_ylabel("comparison ratio")
    ax.legend(frameon=False)

    fig.suptitle(r"Fixed epoch: models 41000--50000, $\alpha_M=0.5$")
    fig.savefig(OUT / "local_anisotropy_fixed_epoch.png", dpi=180)
    plt.close(fig)


def write_outputs(
    profiles: list[dict[str, object]],
    summaries: list[dict[str, float]],
    controls: dict[str, float],
    period_days: float,
) -> None:
    csv_path = OUT / "local_anisotropy_summary.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)

    times = np.asarray([float(profile["time_days"]) for profile in profiles])
    cadence = np.median(np.diff(times)) if times.size > 1 else np.nan
    last = summaries[-1]
    text = [
        "Local algebraic anisotropy diagnostic",
        "=====================================",
        "",
        "This is a post-processing estimate, not a coupled MESA result.",
        "D_omega is the nonlocal TKE diffusion term in Kovacs et al. (2026); it is set to zero.",
        "The closure is nevertheless time dependent through the instantaneous saved TDC state.",
        "The hard 0 <= xi <= 1 clip is diagnostic only.",
        "",
        f"profiles: {len(profiles)}",
        f"median profile cadence: {cadence:.6g} day",
        f"median recent pulsation period: {period_days:.6g} day",
        f"cadence / period: {cadence / period_days:.6g}",
        "The profiles therefore locate radial effects and snapshot-to-snapshot variation,",
        "but do not resolve a smooth phase curve when cadence/period is greater than one.",
        "",
        "Controls read from inlist_pulses:",
    ]
    text.extend(f"  {name} = {value:.16g}" for name, value in controls.items())
    text.extend(
        [
            "",
            "Run-specific TDC_alpha_M history used only for the reconstructed eddy stress:",
            "  models < 51000: TDC_alpha_M = 0.5",
            "  models 51000--52000: TDC_alpha_M = 1.0",
            "  models >= 53000: TDC_alpha_M = 0.7",
            "The saved Eq column already contains the value used during each model.",
        ]
    )
    old_epoch = [
        row for row in summaries if 41000 <= row["model_number"] <= 50000
    ]
    if old_epoch:
        def values(name: str) -> np.ndarray:
            return np.asarray([row[name] for row in old_epoch])

        text.extend(
            [
                "",
                "Fixed-alpha_M comparison (models 41000--50000, alpha_M = 0.5):",
                f"  number of profiles: {len(old_epoch)}",
                f"  median clipped TKE-weighted xi: {np.median(values('xi_clipped_tke_mean')):.6g}",
                f"  median raw-bound violation fraction: {np.median(values('raw_out_of_bounds_fraction')):.6g}",
                f"  median backscatter-zone fraction: {np.median(values('backscatter_zone_fraction')):.6g}",
                f"  anisotropy / Eq power range: {np.min(values('anisotropy_to_Eq_power')):.6g} to {np.max(values('anisotropy_to_Eq_power')):.6g}",
                f"  median anisotropy / Eq power: {np.median(values('anisotropy_to_Eq_power')):.6g}",
                f"  median |q_aniso/q_eddy|: {np.median(values('median_abs_stress_ratio')):.6g}",
                f"  median continuum-current / Eq power check: {np.median(values('continuum_to_Eq_power')):.6g}",
                f"  median |tau_xi| / period: {np.median(values('median_tau_xi_div_period')):.6g}",
            ]
        )
    text.extend(
        [
            "",
            f"Latest saved profile: model {int(last['model_number'])}",
            f"  active turbulent zones: {int(last['active_zones'])}",
            f"  raw xi TKE-weighted p05/median/p95: {last['xi_raw_p05']:.6g} / {last['xi_raw_median']:.6g} / {last['xi_raw_p95']:.6g}",
            f"  clipped xi TKE-weighted mean: {last['xi_clipped_tke_mean']:.6g}",
            f"  raw xi outside [0,1] zone fraction: {last['raw_out_of_bounds_fraction']:.6g}",
            f"  near-singular lambda_xi zone fraction: {last['near_singular_fraction']:.6g}",
            f"  E_nu < 0 (backscatter) zone fraction: {last['backscatter_zone_fraction']:.6g}",
            f"  integrated anisotropy transfer: {last['anisotropy_power_erg_s']:.6e} erg/s",
            f"  integrated raw anisotropy transfer: {last['raw_anisotropy_power_erg_s']:.6e} erg/s",
            f"  integrated profile Eq: {last['mesa_Eq_power_erg_s']:.6e} erg/s",
            f"  anisotropy / Eq integrated power: {last['anisotropy_to_Eq_power']:.6g}",
            f"  raw anisotropy / Eq integrated power: {last['raw_anisotropy_to_Eq_power']:.6g}",
            f"  continuum current-stress / Eq power check: {last['continuum_to_Eq_power']:.6g}",
            f"  median |q_aniso/q_eddy|: {last['median_abs_stress_ratio']:.6g}",
            f"  median |tau_xi| / period: {last['median_tau_xi_div_period']:.6g}",
            "",
            "Interpretation rules:",
            "  E_nu > 0 transfers coherent kinetic energy into TKE; E_nu < 0 is backscatter.",
            "  The current MESA Eq term is positive definite, so it cannot backscatter.",
            "  Frequent raw bound violations or lambda_xi ~= 0 require smooth regularization",
            "  before the algebraic closure can be coupled to the nonlinear equations.",
        ]
    )
    (OUT / "local_anisotropy_diagnostics.txt").write_text("\n".join(text) + "\n")


def main() -> None:
    controls = {
        name: read_real_control(name)
        for name in (
            "mixing_length_alpha",
            "harmonic_dissipation_length_beta",
            "TDC_alpha_C",
            "TDC_alpha_S",
            "TDC_alpha_D",
            "TDC_alpha_M",
        )
    }
    paths = sorted(LOGS.glob("profile*.data"), key=profile_number)
    if not paths:
        raise RuntimeError(f"no profiles in {LOGS}")
    profiles = [compute_profile(path, controls) for path in paths]
    period_models, periods = read_period_series()
    profile_periods = [
        period_at_model(int(profile["model"]), period_models, periods)
        for profile in profiles
    ]
    summaries = [
        summarize(profile, period_days)
        for profile, period_days in zip(profiles, profile_periods)
    ]
    period_days = profile_periods[-1]
    plot_latest(profiles[-1], period_days)
    plot_evolution(summaries)
    plot_fixed_epoch(summaries)
    write_outputs(profiles, summaries, controls, period_days)
    print((OUT / "local_anisotropy_diagnostics.txt").read_text())


if __name__ == "__main__":
    main()
