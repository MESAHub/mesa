#!/usr/bin/env python3
"""Standalone algebra checks for the proposed RSP2 moment discretization.

This does not compile, import, or execute MESA. It checks the design formulas,
not the incomplete Fortran implementation or stellar-model convergence.
"""

import cmath
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np

mp.mp.dps = 70


def thermal_difference(form, y, grad_l, grad_ad, log_p_jump, coefficient, weight):
    """Return predicted ln(T_inner/T_outer) - grad_ad*ln(P_inner/P_outer).

    coefficient = -dm_bar*dlnPdm_qhse for the standard temperature row,
    or dm_bar*kappa*Lrad_per_gradT/(c*area**2*lambda*Prad_outer) for radiation.
    The coefficient includes the radiation opacity floor and flux factor.
    """
    log_t_ad = grad_ad * log_p_jump
    if form == "logarithmic":
        return log_p_jump * (y + (grad_l - grad_ad))
    if form == "standard":
        relative_t_ad = math.expm1(log_t_ad)
        jump_ad = relative_t_ad / (1 + weight * relative_t_ad)
        jump_difference = coefficient * y + (coefficient * grad_l - jump_ad)
        return math.log1p((1 - weight) * jump_difference / (1 + (1 - weight) * jump_ad)) - math.log1p(
            -weight * jump_difference / (1 - weight * jump_ad)
        )
    if form == "radiation":
        relative_prad_ad = math.expm1(4 * log_t_ad)
        jump_difference = coefficient * y + (coefficient * grad_l - relative_prad_ad)
        return 0.25 * math.log1p(jump_difference / (1 + relative_prad_ad))
    raise ValueError(form)


def reference(form, y, grad_l, grad_ad, log_p_jump, coefficient, weight):
    y, grad_l, grad_ad, log_p_jump, coefficient, weight = map(
        mp.mpf, (y, grad_l, grad_ad, log_p_jump, coefficient, weight)
    )
    if form == "logarithmic":
        log_t = (grad_l + y) * log_p_jump
    elif form == "standard":
        jump = coefficient * (grad_l + y)
        log_t = mp.log1p((1 - weight) * jump) - mp.log1p(-weight * jump)
    else:
        log_t = mp.log1p(coefficient * (grad_l + y)) / 4
    return log_t - grad_ad * log_p_jump


def main():
    report = {}
    largest_error = 0.0
    count = 0
    # Different face weights, pressure jumps, gradient offsets and radiation
    # coefficients include the effect of choosing a different opacity/lambda.
    for form in ("logarithmic", "standard", "radiation"):
        for weight in (0.1, 0.5, 0.93):
            for log_p in (1e-9, 1e-5, 0.01, 0.3):
                for grad_ad in (0.15, 0.4, 0.6):
                    for factor in (0.65, 1.0, 1.6):
                        coefficient = log_p * factor
                        for grad_l in (0.2, 0.4, 0.65):
                            for y in (-0.1, -1e-8, 0.0, 1e-12, 0.04):
                                args = (form, y, grad_l, grad_ad, log_p, coefficient, weight)
                                actual = thermal_difference(*args)
                                expected = reference(*args)
                                error = float(abs(mp.mpf(actual) - expected)) / max(
                                    abs(float(expected)), abs(grad_ad * log_p), abs(y * log_p)
                                )
                                largest_error = max(largest_error, error)
                                assert error < 3e-14, args
                                count += 1
    report["temperature_row_equivalence"] = {
        "cases": count, "largest_error_scaled_by_thermal_jump": largest_error
    }

    # Resolve a small change in Y at the same floating-point neutral baseline.
    # This is not a claim to recover input EOS/pressure digits already lost.
    weak = {}
    y, grad_ad, log_p, weight = 1e-20, 0.4, 0.02, 0.5
    for form in ("logarithmic", "standard", "radiation"):
        if form == "logarithmic":
            grad_l, coefficient = grad_ad, 1.0
            derivative = log_p
        elif form == "standard":
            q = math.expm1(grad_ad * log_p)
            grad_l, coefficient = q / (1 + weight * q), 1.0
            derivative = (1 - weight) / (1 + (1 - weight) * grad_l) + weight / (1 - weight * grad_l)
        else:
            grad_l, coefficient = math.expm1(4 * grad_ad * log_p), 1.0
            derivative = 0.25 / (1 + grad_l)
        resolved = thermal_difference(form, y, grad_l, grad_ad, log_p, coefficient, weight)
        assert math.isclose(resolved, y * derivative, rel_tol=5e-15)
        weak[form] = {"Y": y, "resolved_thermal_difference": resolved,
                      "Y_survives_direct_addition_to_gradL": (grad_l + y != grad_l)}
    report["small_Y_at_fixed_neutral_baseline"] = weak

    # Neutral-state check using the standard row's actual Tpoint/Ppoint
    # quadrature. Assume constant grad_ad, uniform composition, no time
    # centering/reconstruction, and exact discrete hydrostatic balance.
    # Then -dm_bar*dlnPdm_qhse = (P_inner-P_outer)/Ppoint. At Y=0 the
    # differential temperature row is neutral by construction, but converting
    # its temperature jump to logarithms does not preserve that neutrality.
    # High precision below shows a truncation mismatch, not roundoff noise.
    neutral_cases = []
    for weight in (0.1, 0.5, 0.9):
        for log_p in (0.001, 0.01, 0.1):
            p, a, g = map(mp.mpf, (str(log_p), str(weight), "0.4"))
            pressure_jump = mp.expm1(p)/(1+a*mp.expm1(p))
            for y in (0.0, 1e-8):
                temperature_jump = pressure_jump*(g+mp.mpf(str(y)))
                inferred_log_t = (mp.log1p((1-a)*temperature_jump) -
                                  mp.log1p(-a*temperature_jump))
                inferred_superadiabaticity = inferred_log_t/p-g
                neutral_cases.append({
                    "weight": weight, "log_pressure_jump": log_p, "Y": y,
                    "log_inversion_superadiabaticity": float(inferred_superadiabaticity),
                    "bias_relative_to_Y": float(inferred_superadiabaticity-y),
                    "scope": "Manufactured standard-row HSE state; not a MESA run."
                })
    # Equal zones: bias = grad_ad*(grad_ad**2-1)*log_pressure_jump**2/12
    # to leading order. Unequal weights generally also give a linear term.
    for weight in (0.1, 0.5, 0.9):
        p, a, g = mp.mpf("1e-7"), mp.mpf(str(weight)), mp.mpf("0.4")
        pressure_jump = mp.expm1(p)/(1+a*mp.expm1(p))
        temperature_jump = g*pressure_jump
        bias = (mp.log1p((1-a)*temperature_jump)-mp.log1p(-a*temperature_jump))/p-g
        leading = (g*(g*g-1)*p*p/12 if weight == 0.5 else
                   (mp.mpf("0.5")-a)*g*(1-g)*p)
        assert abs(bias/leading-1) < mp.mpf("1e-6")
    representative = next(case for case in neutral_cases if
                          case["weight"] == 0.5 and case["log_pressure_jump"] == 0.01
                          and case["Y"] == 1e-8)
    assert representative["log_inversion_superadiabaticity"] < 0
    report["standard_row_neutrality_counterexample"] = {
        "cases": neutral_cases,
        "conclusion": "Row-inversion algebra is accurate but does not preserve the intended Y=0 neutral state.",
        "representative_bias_over_Y": representative["bias_relative_to_Y"]/representative["Y"],
        "limitation": "No EOS noise, varying grad_ad, reconstruction, dynamics or radiation-row test."
    }

    # Independently verify derivatives of the row inversions, holding the
    # coefficient fixed. Full MESA AD must also differentiate that coefficient.
    derivative_error = 0.0
    for weight in (0.1, 0.5, 0.93):
        for jump in (-0.1, 1e-8, 0.05):
            h = 1e-25
            numerical = (cmath.log(1 + (1 - weight) * (jump + 1j*h)) -
                         cmath.log(1 - weight * (jump + 1j*h))).imag / h
            analytic = ((1-weight)/(1+(1-weight)*jump) + weight/(1-weight*jump))
            derivative_error = max(derivative_error, abs(numerical-analytic)/abs(analytic))
    assert derivative_error < 2e-15
    report["temperature_inverse_derivative_relative_error"] = derivative_error

    # A constant-coefficient local backward-Euler block at all three moments=0.
    dt, buoyancy, entropy_gradient, rad_rate = 0.1, 0.7, 0.8, 0.03
    energy_jacobian = np.array([
        [1, -dt*buoyancy, 0],
        [-dt*(2/3)*entropy_gradient, 1+dt*rad_rate, -dt*buoyancy],
        [0, -2*dt*entropy_gradient, 1+2*dt*rad_rate],
    ])
    w_jacobian = energy_jacobian.copy()
    w_jacobian[:, 0] = 0
    assert np.linalg.matrix_rank(w_jacobian) == 2
    assert np.linalg.matrix_rank(energy_jacobian) == 3
    report["zero_state_local_Jacobian"] = {
        "rank_with_w": 2, "rank_with_turbulent_energy": 3,
        "energy_coordinate_determinant": float(np.linalg.det(energy_jacobian)),
        "scope": "No shear or kinetic-energy transport; this is not the full hydro Jacobian."
    }

    # But the existing stress and interface diffusion do not become smooth
    # merely by substituting w=sqrt(e_t). Examine one-sided difference quotients.
    boundary = []
    for energy in (1e-4, 1e-8, 1e-12, 1e-16):
        eq_quotient = math.sqrt(energy)/energy  # Eq=coefficient*sqrt(e_t), coefficient=1
        # Neighbour e_t=1; omit fixed dimensional prefactors of existing Lt.
        def old_transport(e):
            return (0.5*math.sqrt(e)+0.5)*(1-e)
        old_lt_quotient = (old_transport(energy)-old_transport(0.0))/energy
        boundary.append({"e_t": energy, "Eq_difference_quotient": eq_quotient,
                         "existing_Lt_difference_quotient": old_lt_quotient})
    assert boundary[-1]["Eq_difference_quotient"] > 1e5*boundary[0]["Eq_difference_quotient"]
    assert boundary[-1]["existing_Lt_difference_quotient"] > 1e5*boundary[0]["existing_Lt_difference_quotient"]
    report["energy_coordinate_counterexamples"] = boundary

    # Reconstruct the cell covariance from the face entropy correlation,
    # Pi_face/sqrt(e_face), using the cell's own sqrt(e_t). This changes the
    # source interpolation, not the independent face moment or luminosity.
    rng = np.random.default_rng(19)
    largest_bound_ratio = 0.0
    factorization_error = 0.0
    for _ in range(1000):
        w_cell, w_outer, w_inner = 10**rng.uniform(-4, 2, 3)
        weight_outer, weight_inner = rng.uniform(0.05, 0.95, 2)
        w_faces = np.sqrt([weight_outer*w_cell**2 + (1-weight_outer)*w_outer**2,
                           weight_inner*w_inner**2 + (1-weight_inner)*w_cell**2])
        variance_faces = 10**rng.uniform(-3, 2, 2)
        flux_faces = rng.uniform(-1, 1, 2)*w_faces*np.sqrt((2/3)*variance_faces)
        cell_flux = w_cell*np.mean(flux_faces/w_faces)
        bound = (2/3)*w_cell**2*np.mean(variance_faces)
        largest_bound_ratio = max(largest_bound_ratio, cell_flux**2/bound)
        assert cell_flux**2 <= bound*(1+1e-14)
        buoyancies = rng.uniform(0.1, 2, 2)
        source_div_w = np.mean(buoyancies*flux_faces/w_faces)
        dt, work_coefficient, dissipation, eq_div_w = 0.2, 1.03, 0.7, 0.04
        original = work_coefficient*w_cell**2 + dt*dissipation*w_cell**3 - dt*w_cell*(source_div_w+eq_div_w)
        divided = work_coefficient*w_cell + dt*dissipation*w_cell**2 - dt*(source_div_w+eq_div_w)
        error = abs(original-w_cell*divided)/max(abs(original), abs(w_cell*divided), 1e-100)
        factorization_error = max(factorization_error, error)
        assert error < 3e-14
    report["cell_covariance_reconstruction"] = {
        "cases": 1000,
        "largest_cell_bound_ratio_given_admissible_face_moments": largest_bound_ratio,
        "largest_local_energy_factorization_relative_error": factorization_error,
        "scope": "Factorization: accepted w=0, no Lt, Eq proportional to cell w; not the general u-flag row."
    }

    # Signed buoyancy at an empty cell adjacent to a turbulent cell.
    face_w, face_flux, buoyancy = math.sqrt(0.5), -0.1, 1.0
    old_source = 0.5*buoyancy*face_flux
    reconstructed_source = 0.0*(0.5*buoyancy*face_flux/face_w)
    assert old_source < 0 and reconstructed_source == 0
    report["empty_cell_negative_buoyancy"] = {
        "direct_face_source_average": old_source,
        "reconstructed_cell_source": reconstructed_source,
        "interpretation": "Direct averaging can demand negative cell energy at zero energy and no imported Lt/Eq."
    }

    # Imported Lt/Eq terms: keep w and use a positive quadratic predictor.
    # It is an initial guess only; Newton must still solve the full residual.
    from scipy.optimize import brentq
    predictor_cases = []
    for linear in (-2.0, 0.0, 3.0):
        for imported_energy in (1e-12, 0.01, 2.0):
            work_coefficient, dissipation = 1.2, 0.7
            discriminant = math.sqrt(linear**2+4*work_coefficient*imported_energy)
            if linear >= 0:
                guess = (linear+discriminant)/(2*work_coefficient)
            else:
                guess = 2*imported_energy/(discriminant-linear)
            def residual(w):
                return work_coefficient*w*w+dissipation*w**3-linear*w-imported_energy
            root = brentq(residual, 0, max(1.0, 2*guess), xtol=1e-30, rtol=1e-14)
            error = abs(residual(root))/(work_coefficient*root*root+dissipation*root**3+abs(linear*root)+imported_energy)
            assert guess > 0 and root > 0 and error < 1e-13
            predictor_cases.append({"linear_coefficient": linear, "imported_energy": imported_energy,
                                    "positive_initial_guess": guess, "full_residual_positive_root": root})
    report["imported_energy_initial_guess"] = predictor_cases

    # Dimensional, start-state scales tied to the existing luminosity tolerance.
    area, rho, temperature, luminosity_scale = 3e24, 2e-8, 2e4, 7e36
    flux_ref = luminosity_scale/(area*rho*temperature)
    velocity_ref = (luminosity_scale/(area*rho))**(1/3)
    variance_ref = 1.5*(flux_ref/velocity_ref)**2
    assert math.isclose(area*rho*temperature*flux_ref/luminosity_scale, 1.0)
    # Change the entropy unit by a factor b. The implied temperature unit is
    # divided by b so T*s retains energy units; both normalized rows are invariant.
    for b in (1e-5, 1, 1e6):
        transformed_flux_ref = luminosity_scale/(area*rho*(temperature/b))
        transformed_variance_ref = 1.5*(transformed_flux_ref/velocity_ref)**2
        assert math.isclose(transformed_flux_ref, b*flux_ref, rel_tol=3e-15)
        assert math.isclose(transformed_variance_ref, b*b*variance_ref, rel_tol=3e-15)
    report["normalization"] = {"entropy_flux_reference": flux_ref,
                               "entropy_variance_reference": variance_ref,
                               "checks": "finite zero-moment scale; luminosity equivalence; entropy-unit covariance"}
    report["limitations"] = [
        "Design algebra only; the incomplete Fortran has not been validated.",
        "No MESA compilation or model runs.",
        "Accurate row inversions do not establish small-driving accuracy: the standard-row neutral-state counterexample fails that requirement.",
        "Zero-boundary source reconstruction and local branch selection are design checks, not a full coupled solver test.",
        "Fixed-coefficient derivative check is not a full reconstructed-EOS Jacobian check."
    ]
    path = Path(__file__).resolve().parents[1]/"output/review/rsp2_three_equation_20260919/design_checks.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
