"""Independent derivatives of the proposed undivided local residuals.

This verifies the algebra to be translated into MESA AD. It does not execute
MESA's wrappers, opacity/EOS routines, or the full stellar Jacobian.
"""
from pathlib import Path
import json
import numpy as np

OUT = Path(__file__).resolve().parent
CD = (8 / 3) * np.sqrt(2 / 3)
CPHI = 4 * np.sqrt(2 / 3)
CPI = (CD + CPHI) / 2
rng = np.random.default_rng(20260923)


def residual(x, initial, dt, mass, rho_start, theta, alfap):
    w, Pi, Phi, buoyancy, entropy_gradient, Lambda, rad, rho1, rho2, Eq = x
    rho = np.array([rho1, rho2])
    weight = alfap * mass / (3 * .5 * sum(mass))
    work_new = theta * sum(weight * (1 - rho / rho_start))
    work_start = (1 - theta) * sum(weight * (rho_start / rho - 1))
    w_start, Pi_start, Phi_start = initial
    return np.array([
        (1 + work_new) * w*w - (1 - work_start) * w_start*w_start -
        dt * (buoyancy * Pi - CD * w**3 / Lambda + Eq),
        (1 + work_new / 2) * Pi - np.sqrt(1 - work_start) * Pi_start -
        dt * ((2 / 3) * w*w * entropy_gradient + buoyancy * Phi / 3 -
              (CPI * w / Lambda + rad) * Pi),
        Phi - Phi_start - dt * (2 * entropy_gradient * Pi - (CPHI * w / Lambda + 2*rad) * Phi),
    ])


def derivatives(x, initial, dt, mass, rho_start, theta, alfap):
    w, Pi, Phi, buoyancy, entropy_gradient, Lambda, rad, rho1, rho2, Eq = x
    rho = np.array([rho1, rho2])
    weight = alfap * mass / (3 * .5 * sum(mass))
    work_new = theta * sum(weight * (1 - rho / rho_start))
    work_start = (1 - theta) * sum(weight * (rho_start / rho - 1))
    w_start, Pi_start, Phi_start = initial
    jac = np.zeros((3, 10))
    jac[:, 0] = [2 * (1 + work_new) * w + 3 * dt * CD * w*w / Lambda,
                 -dt * ((4 / 3) * w * entropy_gradient - CPI * Pi / Lambda),
                 dt * CPHI * Phi / Lambda]
    jac[:, 1] = [-dt * buoyancy, 1 + work_new / 2 + dt * (CPI * w / Lambda + rad),
                 -2 * dt * entropy_gradient]
    jac[:, 2] = [0, -dt * buoyancy / 3, 1 + dt * (CPHI * w / Lambda + 2 * rad)]
    jac[:, 3] = [-dt * Pi, -dt * Phi / 3, 0]
    jac[:, 4] = [0, -(2 / 3) * dt * w*w, -2 * dt * Pi]
    jac[:, 5] = [-dt * CD * w**3 / Lambda**2,
                 -dt * CPI * w * Pi / Lambda**2, -dt * CPHI * w * Phi / Lambda**2]
    jac[:, 6] = [0, dt * Pi, 2 * dt * Phi]
    for j in range(2):
        d_new = -theta * weight[j] / rho_start[j]
        d_start = -(1 - theta) * weight[j] * rho_start[j] / rho[j]**2
        jac[:, 7+j] = [d_new * w*w + d_start * w_start*w_start,
                       d_new * Pi / 2 + d_start * Pi_start / (2 * np.sqrt(1 - work_start)), 0]
    jac[:, 9] = [-dt, 0, 0]
    return jac


largest_error = 0.
for _ in range(500):
    w, Phi = 10 ** rng.uniform(-4, 4, 2)
    Pi = rng.uniform(-1.5, 1.5) * w * np.sqrt(Phi)
    # Trial Phi may be negative: residual and Jacobian remain algebraic.
    if rng.random() < .3:
        Phi = -Phi
    rho_start = 10 ** rng.uniform(-3, 3, 2)
    rho = rho_start * 10 ** rng.uniform(-.2, .2, 2)
    x = np.array([w, Pi, Phi, 10 ** rng.uniform(-2, 2), rng.uniform(-3, 3),
                  10 ** rng.uniform(-2, 2), 10 ** rng.uniform(-4, 1), *rho, rng.uniform(0, 1)])
    initial = 10 ** rng.uniform(-2, 2, 3)
    initial[1] *= rng.choice([-1, 1])
    dt = 10 ** rng.uniform(-4, 4)
    mass = 10 ** rng.uniform(-2, 2, 2)
    theta, alfap = rng.uniform(0, 1, 2)
    expected = derivatives(x, initial, dt, mass, rho_start, theta, alfap)
    observed = np.empty_like(expected)
    for j in range(len(x)):
        trial = x.astype(complex)
        trial[j] += 1e-30j
        observed[:, j] = residual(trial, initial, dt, mass, rho_start, theta, alfap).imag / 1e-30
    error = np.max(abs(expected - observed) / np.maximum(1., np.maximum(abs(expected), abs(observed))))
    largest_error = max(largest_error, float(error))
assert largest_error < 1e-12
report = {'cases': 500, 'derivatives_per_case': 30,
          'maximum_scaled_derivative_error': largest_error,
          'negative_Phi_trials_checked': True,
          'scope': 'Proposed residual algebra, not the native MESA Jacobian.'}
(OUT / 'residual_derivative_checks.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
