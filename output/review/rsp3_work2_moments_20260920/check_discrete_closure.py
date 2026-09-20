"""Checks of the proposed RSP3 discretization, without compiling or running MESA.

The mathematical test unknowns are radial velocity variance, Pi, and Phi.
MESA would continue to store w, Pi, and Phi. No stellar state is modified.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq

OUT = Path(__file__).resolve().parent
CD = (8 / 3) * np.sqrt(2 / 3)
CPHI = 4 * np.sqrt(2 / 3)
rng = np.random.default_rng(20260921)


def min_covariance_eigenvalue(moments):
    radial_variance, Pi, Phi = moments
    covariance = np.array([[radial_variance, Pi], [Pi, Phi]])
    # Different physical units are scaled separately before the eigenvalue test.
    scales = np.sqrt(np.maximum(np.abs([radial_variance, Phi]), 1e-280))
    covariance = covariance / scales[:, None] / scales[None, :]
    return float(np.linalg.eigvalsh(covariance)[0])


def pressure_coefficients(rho, rho_start, dm, alfap, theta):
    dm_face = 0.5 * np.sum(dm)
    new_work = alfap * theta * np.sum(dm * (1 - rho / rho_start)) / (3 * dm_face)
    old_work = alfap * (1 - theta) * np.sum(dm * (rho_start / rho - 1)) / (3 * dm_face)
    return new_work, old_work


def matrix(w, dt, buoyancy, entropy_gradient, Lambda, inverse_rad_time,
           new_work, alfa_pi=1.0):
    energy_decay = CD * w / Lambda
    variance_decay = CPHI * w / Lambda
    Pi_decay = alfa_pi * 0.5 * (energy_decay + variance_decay) + inverse_rad_time
    return np.array([
        [1 + new_work + dt * energy_decay, -(2 / 3) * dt * buoyancy, 0],
        [-dt * entropy_gradient, 1 + new_work / 2 + dt * Pi_decay, -dt * buoyancy / 3],
        [0, -2 * dt * entropy_gradient, 1 + dt * (variance_decay + 2 * inverse_rad_time)],
    ])


def implicit_local(initial, dt, buoyancy, entropy_gradient, Lambda,
                   inverse_rad_time, new_work, old_work, Eq, alfa_pi):
    if old_work >= 1:
        raise ValueError('Old pressure work removes all available old energy.')
    forcing = initial * np.array([1 - old_work, np.sqrt(1 - old_work), 1])
    forcing[0] += (2 / 3) * dt * Eq
    if np.all(forcing == 0):
        return np.zeros(3), 0.0, 0.0

    def block(w):
        return matrix(w, dt, buoyancy, entropy_gradient, Lambda,
                      inverse_rad_time, new_work, alfa_pi)

    def margin(w):
        return np.min(np.linalg.eigvals(block(w)).real)

    # Select the branch with a positive covariance resolvent. A large unstable
    # step requires enough nonlinear turnover damping; solve for it, do not
    # cap every evolutionary timestep at the turnover time.
    lower = 0.0
    upper = max(1.0, np.sqrt(1.5 * initial[0]))
    while margin(upper) <= 0:
        upper *= 2
    if margin(0) <= 0:
        threshold = brentq(margin, 0, upper, xtol=1e-14, rtol=1e-14)
        lower = threshold + max(1e-12, 1e-10 * threshold)

    def solve(w):
        return np.linalg.solve(block(w), forcing)

    def consistency(w):
        return w * w - 1.5 * solve(w)[0]

    if consistency(lower) > 0:
        raise ArithmeticError('No bracket on the admissible branch.')
    upper = max(upper, 2 * lower)
    while consistency(upper) < 0:
        upper *= 2
    w = brentq(consistency, lower, upper, xtol=1e-280, rtol=4e-15)
    result = solve(w)
    error = abs(consistency(w)) / max(w * w, abs(1.5 * result[0]), 1e-280)
    residual = block(w) @ result - forcing
    relative_residual = np.max(abs(residual)) / max(np.max(abs(block(w)) @ abs(result)),
                                                 np.max(abs(forcing)), 1e-280)
    return result, error, relative_residual


report = {}
worst_pressure_error = 0.0
for _ in range(1000):
    rho_start = 10 ** rng.uniform(-4, 4, 2)
    rho = rho_start * 10 ** rng.uniform(-0.25, 0.25, 2)
    dm = 10 ** rng.uniform(-2, 2, 2)
    alfap, theta = rng.uniform(0, 1, 2)
    energy, energy_start = 10 ** rng.uniform(-4, 4, 2)
    new_work, old_work = pressure_coefficients(rho, rho_start, dm, alfap, theta)
    native = np.sum((alfap / 3) * dm / (0.5 * np.sum(dm)) *
                    (theta * rho * energy + (1 - theta) * rho_start * energy_start) *
                    (1 / rho - 1 / rho_start))
    rewritten = new_work * energy + old_work * energy_start
    error = abs(native - rewritten) / max(abs(new_work * energy) +
                                         abs(old_work * energy_start), 1e-280)
    worst_pressure_error = max(worst_pressure_error, error)
assert worst_pressure_error < 2e-12
report['pressure_work'] = {'cases': 1000, 'max_scaled_rewrite_error': worst_pressure_error}

minimum, max_consistency, max_residual = 1.0, 0.0, 0.0
counts = {'stable': 0, 'unstable': 0, 'zero_energy': 0, 'zero_variance': 0}
for trial in range(800):
    radial_variance, Phi = 10 ** rng.uniform(-3, 3, 2)
    Pi = rng.uniform(-1, 1) * np.sqrt(radial_variance * Phi)
    if trial % 20 == 0:
        radial_variance, Pi = 0.0, 0.0
        counts['zero_energy'] += 1
    elif trial % 20 == 1:
        Phi, Pi = 0.0, 0.0
        counts['zero_variance'] += 1
    initial = np.array([radial_variance, Pi, Phi])
    sign = -1 if trial % 2 else 1
    counts['stable' if sign == -1 else 'unstable'] += 1
    dt = 10 ** rng.uniform(-4, 5)
    buoyancy = 10 ** rng.uniform(-1, 1)
    entropy_gradient = sign * 10 ** rng.uniform(-1, 1)
    Lambda = 10 ** rng.uniform(-1, 1)
    inverse_rad_time = 10 ** rng.uniform(-5, 1)
    rho_start = 10 ** rng.uniform(-4, 4, 2)
    rho = rho_start * 10 ** rng.uniform(-0.2, 0.2, 2)
    new_work, old_work = pressure_coefficients(
        rho, rho_start, 10 ** rng.uniform(-2, 2, 2), rng.uniform(0, 1), rng.uniform(0, 1))
    Eq = 10 ** rng.uniform(-8, 0)
    alfa_pi = rng.uniform(1, 3)
    new, consistency, residual = implicit_local(
        initial, dt, buoyancy, entropy_gradient, Lambda,
        inverse_rad_time, new_work, old_work, Eq, alfa_pi)
    minimum = min(minimum, min_covariance_eigenvalue(new))
    max_consistency = max(max_consistency, consistency)
    max_residual = max(max_residual, residual)
assert minimum >= -1e-10
assert max_consistency < 1e-9
assert max_residual < 1e-12
report['nonlinear_implicit'] = dict(counts, cases=800,
    min_scaled_covariance_eigenvalue=minimum,
    max_energy_consistency_error=max_consistency, max_scaled_residual=max_residual)

# Keeping the old Pi on the right while exporting old kinetic energy is
# inconsistent even with no buoyancy or entropy gradient. Pair the old terms.
initial = np.array([1., 1., 1.])
new_work, old_work = 0.2, 0.2
unpaired = initial / np.array([1 + new_work, 1 + new_work / 2, 1])
unpaired[0] *= 1 - old_work
paired = unpaired.copy()
paired[1] *= np.sqrt(1 - old_work)
assert min_covariance_eigenvalue(unpaired) < 0
assert min_covariance_eigenvalue(paired) > 0
report['pressure_pairing'] = {
    'unpaired_min_scaled_eigenvalue': min_covariance_eigenvalue(unpaired),
    'paired_min_scaled_eigenvalue': min_covariance_eigenvalue(paired),
}

# Linear interpolation in the actual solver variables need not preserve the
# covariance bound, even when both endpoint states are admissible.
w, Pi, Phi = 0.5, 0.5 * np.sqrt(2 / 3), 0.5
det = (2 / 3) * w * w * Phi - Pi * Pi
assert det < 0
report['linear_w_trial_counterexample'] = {'covariance_determinant': det}

# Conservative common implicit transport on unequal face control volumes.
mass = 10 ** rng.uniform(-2, 2, 9)
conductance = 10 ** rng.uniform(-2, 2, 8)
stiffness = np.zeros((9, 9))
for k, coefficient in enumerate(conductance):
    stiffness[k:k + 2, k:k + 2] += coefficient * np.array([[1, -1], [-1, 1]])
start = []
for _ in range(9):
    U, Phi = 10 ** rng.uniform(-2, 2, 2)
    start.append([U, rng.uniform(-1, 1) * np.sqrt(U * Phi), Phi])
start = np.array(start)
new = np.linalg.solve(np.diag(mass) + 100 * stiffness, mass[:, None] * start)
transport_min = min(min_covariance_eigenvalue(row) for row in new)
conservation = np.max(abs(mass @ (new - start))) / np.max(mass @ abs(start))
assert transport_min > 0
assert conservation < 1e-11
report['unequal_mass_transport'] = {'min_scaled_eigenvalue': transport_min,
                                  'max_relative_conservation_error': conservation}

(OUT / 'discrete_closure_checks.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
