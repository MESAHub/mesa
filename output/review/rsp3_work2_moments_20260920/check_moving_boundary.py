"""Nonlinear face-moment experiments with prescribed moving stratification.

These are small dimensionless systems, not MESA models. Log energy is used
only in this independent reference solver. The proposed MESA unknown is w.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import least_squares

OUT = Path(__file__).resolve().parent
CD = (8 / 3) * np.sqrt(2 / 3)
CPHI = 4 * np.sqrt(2 / 3)
CPI = (CD + CPHI) / 2
nz = 12
mass = np.linspace(.5, 1.5, nz)


def equilibrium(gradient):
    energy = np.where(gradient > 0, (3 / 16) * gradient, 1e-4)
    Pi = np.where(gradient > 0, CD * energy**1.5, 0.)
    Phi = np.where(gradient > 0, 2 * gradient * Pi / (CPHI * np.sqrt(energy)), 1e-4)
    return np.column_stack(((2 / 3) * energy, Pi, Phi))


def solve(initial, gradient, dt, alfat, guess=None, rate_scaled=False):
    steady = equilibrium(gradient)
    # A convex covariance predictor, rather than separate clipping or floors.
    fraction = dt / (1 + dt)
    if guess is None:
        guess = (1 - fraction) * initial + fraction * steady
    scale = np.maximum(np.maximum(abs(initial), abs(steady)), 1e-4)
    w_reference = np.sqrt(1.5 * scale[:, 0])
    rate_scale = np.column_stack((
        (2/3) * scale[:, 1] + CD * w_reference * scale[:, 0],
        abs(gradient) * scale[:, 0] + scale[:, 2] / 3 + CPI * w_reference * scale[:, 1],
        2 * abs(gradient) * scale[:, 1] + CPHI * w_reference * scale[:, 2],
    ))
    reference_transport = alfat * .5 * (w_reference[1:] + w_reference[:-1])[:, None] * (scale[:-1] + scale[1:])
    rate_scale[:-1] += reference_transport / mass[:-1, None]
    rate_scale[1:] += reference_transport / mass[1:, None]
    rate_scale += scale / dt

    def unpack(x):
        x = x.reshape(nz, 3)
        moments = x * scale
        moments[:, 0] = np.exp(x[:, 0])
        return moments

    def residual(x):
        moments = unpack(x)
        U, Pi, Phi = moments.T
        w = np.sqrt(1.5 * U)
        source_factor = 1. if rate_scaled else dt
        result = (moments - initial) / (dt if rate_scaled else 1.)
        result[:, 0] -= source_factor * ((2 / 3) * Pi - CD * w * U)
        result[:, 1] -= source_factor * (gradient * U + Phi / 3 - CPI * w * Pi)
        result[:, 2] -= source_factor * (2 * gradient * Pi - CPHI * w * Phi)
        # Common new-state diffusion coefficient and identical face volumes.
        if alfat > 0:
            flux = alfat * .5 * (w[1:] + w[:-1])[:, None] * (moments[:-1] - moments[1:])
            result[:-1] += source_factor * flux / mass[:-1, None]
            result[1:] -= source_factor * flux / mass[1:, None]
        return (result / (rate_scale if rate_scaled else scale)).ravel()

    x0 = guess / scale
    x0[:, 0] = np.log(guess[:, 0])
    lower = np.full((nz, 3), -np.inf)
    upper = np.full((nz, 3), np.inf)
    # Reference solver overflow guards, well outside all physical solutions.
    lower[:, 0], upper[:, 0] = -500, 500
    solution = least_squares(residual, x0.ravel(), jac='cs',
                             bounds=(lower.ravel(), upper.ravel()),
                             ftol=1e-13, xtol=1e-13, gtol=1e-13, max_nfev=600)
    moments = unpack(solution.x)
    error = np.max(abs(residual(solution.x)))
    ratio = moments[:, 1]**2 / (moments[:, 0] * moments[:, 2])
    return moments, float(error), float(np.max(ratio)), solution.nfev


def positive_predictor(initial, gradient, dt, alfat, current, steps):
    """Damped fixed-point predictor with a positive covariance resolvent.

    The added diagonal term also appears on the right and vanishes at a
    fixed point. It changes neither the target equation nor its initial data.
    """
    minimum_Phi = np.inf
    maximum_ratio = 0.
    for _ in range(steps):
        w = np.sqrt(1.5 * current[:, 0])
        operator = np.zeros((nz, nz))
        conductance = alfat * .5 * (w[1:] + w[:-1])
        for k, coefficient in enumerate(conductance):
            operator[k:k+2, k:k+2] += coefficient * np.array([[1, -1], [-1, 1]]) / mass[k:k+2, None]
        block = np.kron(dt * operator, np.eye(3))
        for k in range(nz):
            block[3*k:3*k+3, 3*k:3*k+3] += np.array([
                [1 + dt * CD * w[k], -(2 / 3) * dt, 0],
                [-dt * gradient[k], 1 + dt * CPI * w[k], -dt / 3],
                [0, -2 * dt * gradient[k], 1 + dt * CPHI * w[k]],
            ])
        damping = 1 + dt * (CPHI * np.max(w) + np.max(np.sqrt(abs(gradient))) +
                             np.max(np.diag(operator)))
        minimum_eigenvalue = np.min(np.linalg.eigvals(block).real)
        damping = max(damping, 1 - minimum_eigenvalue)
        current = np.linalg.solve(block + damping * np.eye(3*nz),
                                   (initial + damping * current).ravel()).reshape(nz, 3)
        minimum_Phi = min(minimum_Phi, float(current[:, 2].min()))
        maximum_ratio = max(maximum_ratio, float(np.max(current[:, 1]**2 / (current[:, 0] * current[:, 2]))))
    return current, minimum_Phi, maximum_ratio


initial_gradient = np.where(np.arange(nz) < 4, 1., -1.)
initial = equilibrium(initial_gradient)
results = []
for alfat in [0., .1, 1.]:
    for boundary in [1, 4, 9, 12]:
        gradient = np.where(np.arange(nz) < boundary, 1., -1.)
        for dt in [1e-3, 1., 1e3, 1e6]:
            moments, error, ratio, evaluations = solve(initial, gradient, dt, alfat)
            record = dict(alfat=alfat, boundary=boundary, dt=dt,
                          max_scaled_residual=error, max_covariance_ratio=ratio,
                          minimum_Phi=float(moments[:, 2].min()), evaluations=evaluations)
            results.append(record)

report = {'cases': results,
          'max_scaled_residual': max(r['max_scaled_residual'] for r in results),
          'max_covariance_ratio': max(r['max_covariance_ratio'] for r in results),
          'minimum_Phi': min(r['minimum_Phi'] for r in results),
          'scope': 'Prescribed gradient, fixed unequal face masses, constant Lambda, no viscosity or pressure work.'}
(OUT / 'moving_boundary_direct_checks.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({k: v for k, v in report.items() if k != 'cases'}, indent=2))
print('Cases with residual above 1e-7:', sum(r['max_scaled_residual'] > 1e-7 for r in results))
print('Cases outside covariance domain:', sum(r['minimum_Phi'] < 0 or r['max_covariance_ratio'] > 1 + 1e-8 for r in results))

continued = []
for record in results:
    if record['max_scaled_residual'] <= 1e-7 and record['minimum_Phi'] >= 0 and record['max_covariance_ratio'] <= 1 + 1e-8:
        continued.append(dict(record, continuation_stages=0))
        continue
    gradient = np.where(np.arange(nz) < record['boundary'], 1., -1.)
    guess = initial.copy()
    stages = np.geomspace(1e-3, record['dt'], max(2, int(np.ceil(np.log2(record['dt'] / 1e-3))) + 1))
    error, ratio, evaluations = np.inf, np.inf, 0
    failed = False
    for dt in stages:
        guess, error, ratio, evaluations = solve(initial, gradient, dt, record['alfat'], guess)
        if error > 1e-7 or np.min(guess[:, 2]) < 0 or ratio > 1 + 1e-8:
            failed = True
            break
    continued.append(dict(record, max_scaled_residual=error, max_covariance_ratio=ratio,
                          minimum_Phi=float(guess[:, 2].min()), continuation_stages=len(stages),
                          continuation_failed=failed, last_trial_dt=float(dt)))
report = {'cases': continued,
          'max_scaled_residual': max(r['max_scaled_residual'] for r in continued),
          'max_covariance_ratio': max(r['max_covariance_ratio'] for r in continued),
          'minimum_Phi': min(r['minimum_Phi'] for r in continued),
          'failed_cases': sum(r.get('continuation_failed', False) for r in continued),
          'scope': 'Algebraic continuation holds the original starting moments fixed.'}
(OUT / 'moving_boundary_checks.json').write_text(json.dumps(report, indent=2) + '\n')
print('After algebraic continuation:')
print(json.dumps({k: v for k, v in report.items() if k != 'cases'}, indent=2))

predicted = []
for record in continued:
    if not record.get('continuation_failed', False):
        predicted.append(dict(record, positive_predictor_steps=0))
        continue
    gradient = np.where(np.arange(nz) < record['boundary'], 1., -1.)
    predictor = initial.copy()
    minimum_Phi, maximum_ratio = np.inf, 0.
    steps = 0
    failed = True
    for batch in [10, 40, 150, 800]:
        predictor, trial_min, trial_ratio = positive_predictor(
            initial, gradient, record['dt'], record['alfat'], predictor, batch)
        steps += batch
        minimum_Phi = min(minimum_Phi, trial_min)
        maximum_ratio = max(maximum_ratio, trial_ratio)
        moments, error, ratio, evaluations = solve(initial, gradient, record['dt'], record['alfat'], predictor)
        if error <= 1e-7 and np.min(moments[:, 2]) >= 0 and ratio <= 1 + 1e-8:
            failed = False
            break
    predicted.append(dict(record, max_scaled_residual=error, max_covariance_ratio=ratio,
                          minimum_Phi=float(moments[:, 2].min()), positive_predictor_steps=steps,
                          predictor_minimum_Phi=minimum_Phi, predictor_maximum_ratio=maximum_ratio,
                          predictor_failed=failed, evaluations=evaluations))
report = {'cases': predicted,
          'max_scaled_residual': max(r['max_scaled_residual'] for r in predicted),
          'max_covariance_ratio': max(r['max_covariance_ratio'] for r in predicted),
          'minimum_Phi': min(r['minimum_Phi'] for r in predicted),
          'failed_cases': sum(r.get('predictor_failed', False) for r in predicted),
          'scope': 'Physical equations and initial moments stay fixed during predictor and Newton solve.'}
(OUT / 'moving_boundary_predictor_checks.json').write_text(json.dumps(report, indent=2) + '\n')
print('After positive predictor:')
print(json.dumps({k: v for k, v in report.items() if k != 'cases'}, indent=2))

rate_results = []
for record in results:
    gradient = np.where(np.arange(nz) < record['boundary'], 1., -1.)
    moments, error, ratio, evaluations = solve(initial, gradient, record['dt'], record['alfat'], rate_scaled=True)
    predictor = initial.copy()
    steps = 0
    for batch in [10, 40, 150, 800]:
        if error <= 1e-10 and np.min(moments[:, 2]) >= 0 and ratio <= 1 + 1e-8:
            break
        predictor, _, _ = positive_predictor(initial, gradient, record['dt'], record['alfat'], predictor, batch)
        steps += batch
        moments, error, ratio, evaluations = solve(initial, gradient, record['dt'], record['alfat'],
                                                   predictor, rate_scaled=True)
    rate_results.append(dict(alfat=record['alfat'], boundary=record['boundary'], dt=record['dt'],
                             max_scaled_rate_residual=error, max_covariance_ratio=ratio,
                             minimum_Phi=float(moments[:, 2].min()), positive_predictor_steps=steps,
                             evaluations=evaluations))
rate_report = {'cases': rate_results,
               'max_scaled_rate_residual': max(r['max_scaled_rate_residual'] for r in rate_results),
               'max_covariance_ratio': max(r['max_covariance_ratio'] for r in rate_results),
               'minimum_Phi': min(r['minimum_Phi'] for r in rate_results),
               'failed_cases': sum(r['max_scaled_rate_residual'] > 1e-10 or r['minimum_Phi'] < 0 or
                                   r['max_covariance_ratio'] > 1 + 1e-8 for r in rate_results),
               'scope': 'Fixed normalization combines state-change and source-rate scales; physical residual unchanged.'}
(OUT / 'moving_boundary_rate_checks.json').write_text(json.dumps(rate_report, indent=2) + '\n')
print('With rate normalization and positive predictor:')
print(json.dumps({k: v for k, v in rate_report.items() if k != 'cases'}, indent=2))
