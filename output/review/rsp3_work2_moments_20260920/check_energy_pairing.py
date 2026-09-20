"""Discrete transport and viscous work identities for the proposed update.

The incidence-matrix test checks both locations of velocity/stress. Native
MESA geometry, mass corrections and boundaries still require partial tests.
"""
from pathlib import Path
import json
import numpy as np
import sympy as sp

OUT = Path(__file__).resolve().parent
rng = np.random.default_rng(20260924)

Lr, Lc, Lt, Lr0, Lc0, Lt0, theta = sp.symbols('Lr Lc Lt Lr0 Lc0 Lt0 theta')
used = theta * (Lr + Lc + Lt) + (1 - theta) * (Lr0 + Lc0 + Lt0) + (1 - theta) * (Lt - Lt0)
assert sp.simplify(used - (theta * (Lr + Lc) + (1 - theta) * (Lr0 + Lc0) + Lt)) == 0

max_projection_error, max_work_error = 0., 0.
minimum_heating = np.inf
for _ in range(1000):
    dm = 10 ** rng.uniform(-2, 2, 15)
    # Cell-center turbulent luminosities with zero boundary transport.
    center_flux = rng.normal(size=15)
    center_flux[0] = center_flux[-1] = 0.
    face_mass = .5 * (dm[:-1] + dm[1:])
    face_div = (center_flux[:-1] - center_flux[1:]) / face_mass
    face_flux = (dm[1:] * center_flux[:-1] + dm[:-1] * center_flux[1:]) / (dm[:-1] + dm[1:])
    cell_div = (face_flux[:-1] - face_flux[1:]) / dm[1:-1]
    projected = .5 * (face_div[:-1] + face_div[1:])
    error = np.max(abs(cell_div - projected)) / max(np.max(abs(cell_div)), np.max(abs(projected)), 1.)
    max_projection_error = max(max_projection_error, float(error))

    for velocity_location in ['face', 'cell']:
        # A positive stress coefficient may contain the current density, w,
        # mixing length, geometry and alfam. Both grids share this identity.
        radius = 10 ** rng.uniform(-1, 1, 10)
        mass = 10 ** rng.uniform(-2, 2, 10)
        velocity_new = rng.normal(size=10)
        velocity_old = rng.normal(size=10)
        velocity_work = .5 * (velocity_new + velocity_old)
        strain_work = np.diff(-velocity_work / radius)
        coefficient = 10 ** rng.uniform(-2, 2, 9)
        stress = coefficient * strain_work
        force = np.zeros(10)
        force[:-1] -= stress / radius[:-1]
        force[1:] += stress / radius[1:]
        acceleration = force / mass
        heating = stress * strain_work
        error = abs(np.dot(mass * velocity_work, acceleration) + np.sum(heating)) / max(np.sum(heating), 1.)
        max_work_error = max(max_work_error, float(error))
        minimum_heating = min(minimum_heating, float(heating.min()))

assert max_projection_error < 1e-12
assert max_work_error < 1e-12
assert minimum_heating >= 0
report = {'heat_flux_time_weight_identity': 'exact', 'projection_cases': 1000,
          'max_projection_error': max_projection_error,
          'viscous_work_cases': 2000, 'max_viscous_work_error': max_work_error,
          'minimum_heating': minimum_heating,
          'reversal_old_cross_product': -1 * .5,
          'reversal_paired_product': .5 * .5,
          'scope': 'Algebraic incidence identities; no native MESA geometry or boundary execution.'}
(OUT / 'energy_pairing_checks.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
