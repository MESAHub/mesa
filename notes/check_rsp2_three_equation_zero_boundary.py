#!/usr/bin/env python3
"""Two-cell moment/energy prototype, with fixed thermodynamics and mean flow.

Tests the proposed source reconstruction and zero-state branch selection.
The u/v cases retain the respective face/cell dependence of viscous heating;
this is not MESA, a complete momentum solve, or a full hydro Jacobian test.
"""

import json
import math
from pathlib import Path

import numpy as np


def solve_case(velocity, alfat, alfam, dt, entropy_gradient, start,
               divide_energy_by_w=False, theta_L=1.0):
    # Surface-side cell has mass 2, inner cell mass 1. No boundary fluxes.
    mass = np.array([2.0, 1.0])
    weight = np.array([1/3, 2/3])
    beta_flux = 6*math.sqrt(2/3)
    beta_variance = 4*math.sqrt(2/3)
    dissipation = (8/3)*math.sqrt(2/3)
    buoyancy, radiation_rate = 0.7, 0.1
    start = np.array(start, dtype=float)

    def terms(x):
        w = x[:2]
        flux, variance = x[2:]
        energy_face = np.dot(weight, w*w)
        w_face = np.sqrt(energy_face) if energy_face.real > 0 else 0*x[0]
        source = 0*w
        source_div_w = 0*x[0]
        if w_face.real > 0:
            source = 0.5*buoyancy*(w/w_face)*flux
            source_div_w = 0.5*buoyancy*flux/w_face
        mean_w = np.dot(weight, w)
        transported_luminosity = alfat*mean_w*(w[0]**2-w[1]**2)
        transport = np.array([-transported_luminosity/mass[0], transported_luminosity/mass[1]])
        if velocity == "v":
            eq_div_w = alfam*np.array([0.4, 0.7])
            eq = eq_div_w*w
        else:
            eq_div_w = np.zeros(2)
            eq = np.full(2, 0.5*alfam*mean_w, dtype=x.dtype)
        flux_rhs = ((2/3)*energy_face*entropy_gradient + buoyancy*variance -
                    (beta_flux*w_face+radiation_rate)*flux)
        variance_rhs = 2*entropy_gradient*flux-(beta_variance*w_face+2*radiation_rate)*variance
        return source, source_div_w, transport, eq, eq_div_w, flux_rhs, variance_rhs

    def residual(x, select_branch=True):
        w = x[:2]
        source, source_div_w, transport, eq, eq_div_w, flux_rhs, variance_rhs = terms(x)
        transport_start = terms(start)[2]
        weighted_transport = theta_L*transport+(1-theta_L)*transport_start
        residual_energy = w*w-start[:2]**2-dt*(source+weighted_transport+eq-dissipation*w**3)
        residual_flux = x[2]-start[2]-dt*flux_rhs
        residual_variance = x[3]-start[3]-dt*variance_rhs
        residuals = np.array([*residual_energy, residual_flux, residual_variance])
        if not select_branch:
            return residuals
        for k in range(2):
            if divide_energy_by_w and w[k].real > 0:
                residuals[k] = residual_energy[k]/w[k]
            if start[k] == 0 and alfat == 0 and (velocity == "v" or alfam == 0):
                divided = w[k]+dt*dissipation*w[k]**2-dt*(source_div_w+eq_div_w[k])
                residuals[k] = w[k] if w[k].real <= divided.real else divided
        if np.all(w.real == 0) and np.all(start[:2] == 0) and alfam == 0:
            if start[2] == 0 and start[3] == 0:
                residuals[:] = x  # Exact homogeneous dormant solution.
        return residuals

    # Initial guesses only. Accepted start values are never changed.
    x = start.copy()
    if np.all(x[:2] == 0) and x[3] > 0:
        x[:2] = dt*abs(buoyancy)*math.sqrt(0.5*x[3])
        x[2] = start[2]+dt*buoyancy*start[3]
    else:
        source, source_div_w, transport, eq, eq_div_w, flux_rhs, variance_rhs = terms(x)
        x[2] += dt*flux_rhs
        x[3] = max(0.0, start[3]+2*dt*entropy_gradient*x[2])
        source, source_div_w, transport, eq, eq_div_w, *_ = terms(x)
        if divide_energy_by_w:
            # A quadratic guess from the entire energy RHS and its own-w
            # derivative, with accepted history and neighbours fixed here.
            # The equations solved below still use all CURRENT dependencies.
            saved_guess = x.copy()
            old_transport = terms(start)[2]
            def energy_rhs(trial):
                source, _, transport, eq, *_ = terms(trial)
                return start[:2]**2+dt*(source+theta_L*transport+
                       (1-theta_L)*old_transport+eq-dissipation*trial[:2]**3)
            rhs = energy_rhs(saved_guess)
            for k in range(2):
                trial = saved_guess.astype(complex)
                trial[k] += 1e-28j
                derivative = energy_rhs(trial)[k].imag/1e-28
                available = rhs[k]-derivative*saved_guess[k]
                if available < 0:
                    continue  # Not repaired by changing the physical RHS.
                discriminant = math.sqrt(derivative**2+4*available)
                if derivative > 0:
                    predicted = 0.5*(derivative+discriminant)
                elif available > 0:
                    predicted = 2*available/(discriminant-derivative)
                else:
                    predicted = 0.0
                x[k] = max(x[k], predicted)
        for k in range(2) if not divide_energy_by_w else ():
            if start[k] != 0:
                continue
            imported = dt*(transport[k]+eq[k])
            if imported > 0:
                x[k] = math.sqrt(imported)
            elif source_div_w+eq_div_w[k] > 0:
                x[k] = dt*(source_div_w+eq_div_w[k])

    last_norm = math.inf
    converged = False
    for iteration in range(1, 81):
        r = residual(x)
        last_norm = float(np.linalg.norm(r, ord=np.inf))
        physical_norm = float(np.linalg.norm(residual(x, False), ord=np.inf))
        if last_norm < 1e-11 and physical_norm < 1e-11:
            converged = True
            break
        jacobian = np.zeros((4, 4))
        for j in range(4):
            z = x.astype(complex)
            z[j] += 1e-28j
            jacobian[:, j] = residual(z).imag/1e-28
        try:
            delta = np.linalg.solve(jacobian, -r)
        except np.linalg.LinAlgError:
            break
        # A projected trial direction mirrors a nonnegative variable domain;
        # convergence still requires the original physical residual to vanish.
        for j in (0, 1, 3):
            if x[j]+delta[j] < 0:
                local_zero_branch = (j < 2 and start[j] == 0 and alfat == 0 and
                                     (velocity == "v" or alfam == 0))
                if divide_energy_by_w and j < 2 and x[j] > 0 and not local_zero_branch:
                    delta[j] = -0.95*x[j]
                else:
                    delta[j] = -x[j]
        step = 1.0
        for _ in range(32):
            candidate = x+step*delta
            new_norm = np.linalg.norm(residual(candidate), ord=np.inf)
            if new_norm < last_norm or new_norm < 1e-11:
                x = candidate
                break
            step *= 0.5
        else:
            break
    return {
        "velocity_placement": velocity, "alfat": alfat, "alfam": alfam,
        "divide_energy_by_w": divide_energy_by_w, "theta_L": theta_L,
        "dt": dt, "entropy_gradient": entropy_gradient,
        "start": start.tolist(), "converged": converged, "iterations": iteration,
        "selected_residual": last_norm,
        "original_energy_and_moment_residual": float(np.linalg.norm(residual(x, False), ord=np.inf)),
        "solution": x.tolist(),
    }


def main():
    cases = []
    for velocity in ("v", "u"):
        for alfat in (0.0, 0.01, 0.2):
            for alfam in (0.0, 0.25):
                for dt in (1e-3, 0.01, 0.1):
                    for driving in (-0.5, 0.5):
                        flux = math.copysign(0.03, driving)
                        cases.append(solve_case(velocity, alfat, alfam, dt, driving, [1, 0, flux, 0.02]))
    for velocity in ("v", "u"):
        for alfat in (0.0, 0.2):
            cases.append(solve_case(velocity, alfat, 0, 0.01, 0.5, [0, 0, 0, 0]))
            cases.append(solve_case(velocity, alfat, 0, 0.01, 0.5, [0, 0, 0, 0.02]))
            cases.append(solve_case(velocity, alfat, 0.25, 0.01, 0.5, [1, 0, 0, 0]))
    passed = [c for c in cases if c["converged"]]
    failed = [c for c in cases if not c["converged"]]
    report = {"cases": len(cases), "converged": len(passed), "failed": len(failed),
              "max_iterations_for_converged_cases": max(c["iterations"] for c in passed),
              "max_original_residual_for_converged_cases": max(c["original_energy_and_moment_residual"] for c in passed),
              "limitations": "Fixed thermodynamics, geometry and velocities; u/v stress dependence only. Not a full MESA hydro or LNA test.",
              "failed_cases": failed, "all_cases": cases}
    path = Path(__file__).resolve().parents[1]/"output/review/rsp2_three_equation_20260919/zero_boundary_checks.json"
    path.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({k:v for k,v in report.items() if k != "all_cases"}, indent=2))


if __name__ == "__main__":
    main()
