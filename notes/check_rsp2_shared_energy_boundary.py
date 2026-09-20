#!/usr/bin/env python3
"""Standalone checks of a proposed shared RSP2 energy-row treatment.

No MESA compilation or execution. Geometry and thermodynamic coefficients
are fixed; the coupled momentum/EOS solver and LNA are not exercised.
"""

import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from check_rsp2_transport_startup import face_flux, predict, residual_jacobian, solve
from check_rsp2_three_equation_zero_boundary import solve_case


def solve_divided_transport(start, dm, dt, theta, mass_interp, guess):
    w = guess.copy()
    assert np.all(w > 0)
    for iteration in range(80):
        energy, jacobian = residual_jacobian(w, start, dm, dt, theta, mass_interp)
        divided = energy/w
        if np.max(np.abs(energy)) < 1e-13 and np.max(np.abs(divided)) < 1e-12:
            return w, iteration
        # Include the derivative of the divisor, not just a frozen row scale.
        divided_jacobian = jacobian/w[:, None]-np.diag(energy/w**2)
        correction = np.linalg.solve(divided_jacobian, -divided)
        step = 1.0
        decreasing = correction < 0
        if np.any(decreasing):
            step = min(step, .95*np.min(-w[decreasing]/correction[decreasing]))
        for _ in range(40):
            trial = w+step*correction
            trial_energy, _ = residual_jacobian(trial, start, dm, dt, theta, mass_interp)
            if np.max(np.abs(trial_energy/trial)) < np.max(np.abs(divided)):
                w = trial
                break
            step *= .5
        else:
            raise AssertionError(("line search", start, dm, dt, theta, iteration))
    raise AssertionError(("iterations", start, dm, dt, theta))


def main():
    report = {}
    count, worst, energy_error = 0, 0, 0.0
    for ratio in (1e-4, .01, .1, 1., 10., 100., 1e4):
        dm = np.array([ratio, 1.])
        for reverse in (False, True):
            for quiet in (0., 1e-14, 1e-8):
                start = np.array([quiet, 1.])
                if reverse:
                    start = start[::-1].copy()
                for mass_interp in (False, True):
                    for theta in (.5, 1.):
                        for step in (1e-6, .01, .1, .5):
                            dt = step*min(dm)
                            guess = predict(start, dm, dt, theta, mass_interp)
                            root, iterations = solve_divided_transport(start, dm, dt, theta, mass_interp, guess)
                            error = abs(np.dot(dm, root*root-start*start))/sum(dm)
                            assert error < 1e-12
                            energy_error = max(energy_error, error)
                            worst = max(worst, iterations)
                            count += 1
    report["one_equation_transport"] = {
        "cases": count, "max_iterations": worst,
        "max_mass_weighted_energy_error": energy_error,
        "scope": "Original face-w transport and theta_L; no source, viscosity, work or EOS."
    }

    # Previously demonstrated trap: a simple incoming-energy guess still
    # sends the raw-energy Newton step toward negative w in a narrow cell.
    dm, start, dt = np.array([1e-4, 1.]), np.array([0., 1.]), 1e-6
    flux, _ = face_flux(start, dm)
    guess = np.array([math.sqrt(dt*flux/dm[0]), 1.])
    assert solve(start, dm, dt, .5, guess=guess)[0] is None
    root, iterations = solve_divided_transport(start, dm, dt, .5, True, guess)
    report["narrow_cell_Newton_direction"] = {
        "raw_energy_row_stalls": True, "divided_row_iterations": iterations,
        "positive_solution": root.tolist()
    }

    # One-equation local residual, including work/radiative quadratic terms,
    # cubic dissipation and signed linear production. At zero start the
    # divided expression extends analytically to w=0.
    count = 0
    for start_w in (0., 1e-14, 1.):
        for production in (-.5, 0., .5):
            for dt in (.001, .01, .1):
                for quadratic in (.7, 1., 1.3):
                    for dissipation in (0., .7):
                        def energy(w):
                            return quadratic*w*w+dt*dissipation*w**3-dt*production*w-start_w**2
                        def divided(w):
                            return quadratic*w+dt*dissipation*w*w-dt*production-start_w**2/w
                        if start_w == 0 and production <= 0:
                            root = 0.
                            analytic_limit = -dt*production
                            assert min(root, analytic_limit) == 0 and energy(root) == 0
                        else:
                            root = brentq(divided, 1e-100, 10., xtol=1e-100, rtol=1e-14)
                            scale = quadratic*root*root+dt*dissipation*root**3+abs(dt*production*root)+start_w**2
                            assert abs(energy(root))/scale < 2e-14
                            derivative = quadratic+2*dt*dissipation*root+start_w**2/root**2
                            assert derivative > 0
                        count += 1
    report["one_equation_local_branches"] = {"cases": count,
        "scope": "Scalar frozen-coefficient branch/root checks, not a full one-equation hydro solve."}

    # A quiet u-flag block can have collective viscous driving even when
    # each diagonal production derivative is negative. A purely cell-local
    # onset test is insufficient; the predictor must consider neighbours.
    weights = np.array([1/3, 2/3])
    production_jacobian = .125*np.tile(weights, (2, 1))-.1*np.eye(2)
    assert np.all(np.diag(production_jacobian) < 0)
    growth = float(max(np.linalg.eigvals(production_jacobian)))
    dt = .1
    positive_root = np.full(2, dt*growth)
    assert growth > 0 and np.max(np.abs(positive_root**2-dt*production_jacobian@positive_root)) < 1e-19
    report["coupled_u_flag_onset_counterexample"] = {
        "own_w_production_derivatives": np.diag(production_jacobian).tolist(),
        "collective_positive_production_eigenvalue": growth,
        "positive_root": positive_root.tolist(),
        "conclusion": "An own-w quadratic guess alone is not a complete quiet-block activation rule."
    }

    moment_cases = []
    for velocity in ("v", "u"):
        for alfat in (0., .01, .2):
            for alfam in (0., .25):
                for theta in (.5, 1.):
                    for dt in (.001, .01, .1):
                        for driving in (-.5, .5):
                            flux = math.copysign(.03, driving)
                            for start in ([1, 0, flux, .02], [1, 1e-14, flux, .02],
                                          [0, 0, 0, 0], [0, 0, 0, .02], [1, 0, 0, 0]):
                                moment_cases.append(solve_case(velocity, alfat, alfam, dt, driving,
                                                              start, True, theta))
    failed = [case for case in moment_cases if not case["converged"]]
    assert not failed, failed
    report["three_equation_energy_moments"] = {
        "cases": len(moment_cases), "failed": failed,
        "max_iterations": max(case["iterations"] for case in moment_cases),
        "max_original_residual": max(case["original_energy_and_moment_residual"] for case in moment_cases),
        "scope": "Two cells, one interior face; fixed thermodynamics and mean flow; simplified positive viscous-heating coefficients."
    }
    report["limits"] = [
        "No source installation, MESA compilation, model evolution or LNA validation.",
        "No general uniqueness/convergence or moment-realizability proof.",
        "Time-centered mean-strain sign reversals and fully coupled hydro remain untested.",
        "A quiet multi-cell u-flag region may require a coupled heating predictor; local own-w derivatives are not a universal activation test.",
        "Dividing by w is valid for w>0; zero needs a justified limiting row or a positive initial guess, never a denominator floor."
    ]
    path = Path(__file__).resolve().parents[1]/"output/review/rsp2_three_equation_20260919/shared_energy_boundary_checks.json"
    path.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
