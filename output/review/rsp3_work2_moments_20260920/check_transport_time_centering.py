"""Common moment transport with the actual theta time weighting.

The high-frequency counterexample is independent of the stellar Newton solve.
All faces have positive turbulent energy, so the diffusion coefficient is
positive even when entropy variance has compact support.
"""
from pathlib import Path
import json
import numpy as np

OUT = Path(__file__).resolve().parent
rng = np.random.default_rng(20260922)


def diffusion(mass, conductance):
    stiffness = np.zeros((len(mass), len(mass)))
    for k, coefficient in enumerate(conductance):
        stiffness[k:k+2, k:k+2] += coefficient * np.array([[1, -1], [-1, 1]])
    return stiffness / mass[:, None]


def step(old, operator, dt, theta):
    identity = np.eye(len(operator))
    return np.linalg.solve(identity + theta * dt * operator,
                           (identity - (1 - theta) * dt * operator) @ old)


operator = diffusion(np.ones(3), np.ones(2))
old = np.array([[1., 0., 0.], [1., 0., 1.], [1., 0., 0.]])
centered = step(old, operator, 10., .5)
implicit = step(old, operator, 10., 1.)
assert centered[1, 2] < 0
assert np.min(implicit[:, 2]) > 0
assert np.allclose(centered[:, 0], 1.)

minimum_weight, max_conservation = 1., 0.
for _ in range(1000):
    mass = 10 ** rng.uniform(-2, 2, 8)
    conductance = 10 ** rng.uniform(-2, 2, 7)
    operator = diffusion(mass, conductance)
    theta = rng.uniform(.5, 1)
    # Sufficient positivity condition for the explicit old-state part.
    dt = rng.uniform(.01, .99) / ((1 - theta) * np.max(np.diag(operator)))
    transfer = step(np.eye(8), operator, dt, theta)
    minimum_weight = min(minimum_weight, float(transfer.min()))
    error = np.max(abs(mass @ transfer - mass)) / np.max(mass)
    max_conservation = max(max_conservation, float(error))
assert minimum_weight >= 0
assert max_conservation < 2e-11

report = {'centered_variance_dt10': centered[:, 2].tolist(),
          'implicit_variance_dt10': implicit[:, 2].tolist(),
          'positive_theta_cases': 1000,
          'minimum_transfer_weight': minimum_weight,
          'maximum_mass_conservation_error': max_conservation,
          'scope': 'Transport only; pressure work changes the old-state positivity condition.'}
(OUT / 'transport_time_centering_checks.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
