"""Exact identities for the proposed RSP3 closure and its pressure work."""
from pathlib import Path
import json
import sympy as sp

OUT = Path(__file__).resolve().parent
energy, Pi, Phi, buoyancy, entropy_gradient = sp.symbols(
    'energy Pi Phi buoyancy entropy_gradient', real=True)
energy_decay, variance_decay, Pi_decay, Eq, production_factor = sp.symbols(
    'energy_decay variance_decay Pi_decay Eq production_factor', real=True)

determinant = sp.Rational(2, 3) * energy * Phi - Pi**2
energy_rhs = buoyancy * Pi - energy_decay * energy + Eq
Pi_rhs = sp.Rational(2, 3) * energy * entropy_gradient + production_factor * buoyancy * Phi - Pi_decay * Pi
Phi_rhs = 2 * entropy_gradient * Pi - variance_decay * Phi
determinant_rhs = sp.expand(sum(sp.diff(determinant, variable) * rhs for variable, rhs in
                              zip((energy, Pi, Phi), (energy_rhs, Pi_rhs, Phi_rhs))))
expected = (sp.Rational(2, 3) - 2 * production_factor) * buoyancy * Pi * Phi
expected -= (energy_decay + variance_decay) * determinant
expected += (2 * Pi_decay - energy_decay - variance_decay) * Pi**2
expected += sp.Rational(2, 3) * Eq * Phi
assert sp.simplify(determinant_rhs - expected) == 0
corrected = sp.factor(determinant_rhs.subs({
    production_factor: sp.Rational(1, 3),
    Pi_decay: (energy_decay + variance_decay) / 2,
}))
assert sp.simplify(corrected + (energy_decay + variance_decay) * determinant -
                   sp.Rational(2, 3) * Eq * Phi) == 0

rho, rho_start, energy_start, theta = sp.symbols('rho rho_start energy_start theta', positive=True)
native_work = (theta * rho * energy + (1 - theta) * rho_start * energy_start) * (1 / rho - 1 / rho_start)
new_work = theta * (1 - rho / rho_start)
old_work = (1 - theta) * (rho_start / rho - 1)
assert sp.simplify(native_work - new_work * energy - old_work * energy_start) == 0

# The pressure term's linearization belongs in the LNA inertia. The weights
# below are the two adjacent cell masses divided by the face control mass.
drho1, drho2, mass1, mass2, alfap = sp.symbols('drho1 drho2 mass1 mass2 alfap', real=True)
eps = sp.symbols('eps', real=True)
dm_face = (mass1 + mass2) / 2
new_work = alfap * theta / (3 * dm_face) * (-mass1 * eps * drho1 - mass2 * eps * drho2)
old_work = alfap * (1 - theta) / (3 * dm_face) * (
    mass1 * (1 / (1 + eps * drho1) - 1) + mass2 * (1 / (1 + eps * drho2) - 1))
Pi_work = new_work * Pi / 2 + (1 - sp.sqrt(1 - old_work)) * Pi
linear_work = sp.simplify(sp.diff(Pi_work, eps).subs(eps, 0))
expected_work = -alfap * Pi * (mass1 * drho1 + mass2 * drho2) / (6 * dm_face)
assert sp.simplify(linear_work - expected_work) == 0
assert sp.simplify(sp.diff(linear_work, theta)) == 0

report = {'covariance_identity': str(sp.factor(expected)),
          'candidate_covariance_identity': str(corrected),
          'pressure_rewrite': 'exact',
          'Pi_pressure_inertia': str(linear_work),
          'Pi_pressure_inertia_independent_of_theta': True}
(OUT / 'symbolic_closure_checks.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
