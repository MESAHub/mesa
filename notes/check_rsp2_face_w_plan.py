"""Check the fixed-mass algebra in rsp2_face_w_implementation.md.

These are quadrature/remap checks, not MESA or closure evolution tests.
Both mathematical boundary faces are included; imposed boundary values
and their energy exchange require separate implementation checks.
"""
from pathlib import Path
import json
import numpy as np


def face_edges(cell_edges):
    return np.r_[cell_edges[0], .5*(cell_edges[:-1]+cell_edges[1:]), cell_edges[-1]]


def overlap_weights(old_edges, new_edges):
    overlap = np.maximum(0., np.minimum(new_edges[1:, None], old_edges[None, 1:])
                         - np.maximum(new_edges[:-1, None], old_edges[None, :-1]))
    return overlap/np.diff(new_edges)[:, None]


def scaled_error(left, right):
    left, right = np.asarray(left), np.asarray(right)
    return float(np.max(abs(left-right))/max(np.max(abs(left)), np.max(abs(right)), 1.))


rng = np.random.default_rng(20260920)
dm = 10.**rng.uniform(-3, 3, 41)
edges = np.r_[0., np.cumsum(dm)]
dm = np.diff(edges)
dual_edges = face_edges(edges)
mu = np.r_[dm[0]/2, .5*(dm[:-1]+dm[1:]), dm[-1]/2]
w = rng.uniform(.1, 5., len(mu))
energy = w*w
cell_energy = .5*(energy[:-1]+energy[1:])
checks = {'energy_quadrature': scaled_error(dm@cell_energy, mu@energy)}

# Outward center luminosity and the actual outer/inner boundary fluxes.
center_Lt = rng.normal(size=len(dm))
boundary_Lt = rng.normal(size=2)
all_Lt = np.r_[boundary_Lt[0], center_Lt, boundary_Lt[1]]
face_rate = np.diff(all_Lt)/mu
gas_Lt = np.r_[boundary_Lt[0],
               (dm[1:]*center_Lt[:-1]+dm[:-1]*center_Lt[1:])/(dm[:-1]+dm[1:]),
               boundary_Lt[1]]
cell_rate = np.diff(gas_Lt)/dm
checks['Lt_local_projection'] = scaled_error(cell_rate, .5*(face_rate[:-1]+face_rate[1:]))
checks['Lt_global_boundary_budget'] = scaled_error(mu@face_rate, boundary_Lt[1]-boundary_Lt[0])

# Eq/w can have either sign; the allocation identity does not assume positivity.
Eq_div_w_cell = rng.normal(size=len(dm))
Eq_stress_cell = .5*(w[:-1]+w[1:])*Eq_div_w_cell
Eq_div_w_face = np.r_[Eq_div_w_cell[0],
                     (dm[:-1]*Eq_div_w_cell[:-1]+dm[1:]*Eq_div_w_cell[1:])/(dm[:-1]+dm[1:]),
                     Eq_div_w_cell[-1]]
Eq_face = w*Eq_div_w_face
checks['Eq_global_allocation'] = scaled_error(mu@Eq_face, dm@Eq_stress_cell)
checks['Eq_projected_cell_budget'] = scaled_error(dm@(.5*(Eq_face[:-1]+Eq_face[1:])), dm@Eq_stress_cell)

# Time-weighted pressure increments split by their face energy contributions.
rho, rho_start = rng.uniform(.3, 3., (2, len(dm)))
start_energy = rng.uniform(.1, 5., len(mu))**2
theta_P, alfap = .37, .8
dV = 1/rho-1/rho_start
cell_work_outer = alfap/3*(theta_P*rho*energy[:-1]
                            +(1-theta_P)*rho_start*start_energy[:-1])*dV
cell_work_inner = alfap/3*(theta_P*rho*energy[1:]
                            +(1-theta_P)*rho_start*start_energy[1:])*dV
cell_work = cell_work_outer+cell_work_inner
face_work = (np.r_[dm*cell_work_outer, 0.]+np.r_[0., dm*cell_work_inner])/mu
checks['pressure_global_allocation'] = scaled_error(dm@cell_work, mu@face_work)
projected_work = .5*(face_work[:-1]+face_work[1:])
checks['pressure_redistribution_zero_integral'] = scaled_error(dm@(cell_work-projected_work), 0.)
assert np.max(abs(cell_work-projected_work)) > .01  # global identity is not pointwise equality

# Common overlap of the covariance matrix preserves admissibility and integrals.
Phi = rng.uniform(.1, 3., len(mu))
Pi = rng.uniform(-.95, .95, len(mu))*np.sqrt((2/3)*energy*Phi)
moments = np.column_stack((energy, Pi, Phi))
new_edges = np.r_[0., np.sort(rng.uniform(0., edges[-1], 62)), edges[-1]]
new_dual_edges = face_edges(new_edges)
new_mu = np.diff(new_dual_edges)
weights = overlap_weights(dual_edges, new_dual_edges)
new_moments = weights@moments
old_mu_geometry = np.diff(dual_edges)
checks['dual_geometry_matches_quadrature'] = scaled_error(mu, old_mu_geometry)
checks['remap_moment_integrals'] = scaled_error(old_mu_geometry@moments, new_mu@new_moments)
checks['remap_constant'] = scaled_error(weights@np.ones(len(mu)), np.ones(len(new_mu)))
checks['remap_identity'] = scaled_error(overlap_weights(dual_edges, dual_edges), np.eye(len(mu)))
det = (2/3)*new_moments[:, 0]*new_moments[:, 2]-new_moments[:, 1]**2
assert np.min(det) >= -1e-12

# Converting from cell averages conserves total energy, but is not reversible.
old_cell_energy = rng.uniform(.1, 10., len(dm))
converted_face = overlap_weights(edges, dual_edges)@old_cell_energy
returned_cell = .5*(converted_face[:-1]+converted_face[1:])
checks['cell_to_face_energy_integral'] = scaled_error(dm@old_cell_energy, old_mu_geometry@converted_face)
checks['conversion_roundtrip_energy_integral'] = scaled_error(dm@old_cell_energy, dm@returned_cell)
assert np.max(abs(returned_cell-old_cell_energy)) > .01

# Uniform-background old source averages away an alternating Pi pattern.
alternating_Pi = (-1.)**np.arange(32)
old_source = .5*(alternating_Pi+np.roll(alternating_Pi, -1))
assert np.max(abs(old_source)) == 0.
assert np.linalg.norm(alternating_Pi) > 0.

for name, error in checks.items():
    assert error < 2e-10, (name, error)
result = dict(checks=checks, minimum_remapped_covariance_determinant=float(np.min(det)),
              mass_ratio=float(np.max(dm)/np.min(dm)),
              maximum_local_pressure_redistribution=float(np.max(abs(cell_work-projected_work))),
              maximum_roundtrip_cell_energy_change=float(np.max(abs(returned_cell-old_cell_energy))),
              scope='Fixed-mass algebra and positive overlap only. No MESA run or full hydro validation.')
root = Path(__file__).resolve().parents[1]
out = root/'output/review/rsp3_face_w_plan_20260920'
out.mkdir(parents=True, exist_ok=True)
(out/'checks.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
