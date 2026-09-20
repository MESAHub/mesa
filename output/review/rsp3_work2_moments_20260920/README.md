# RSP3 closure investigation

These scripts read copied MESA output or solve independent mathematical
models. They do not compile, relink or run MESA, and do not change a test case.

The active derivation is in `notes/rsp3_covariance_discretization.md`.
The evidence and literature normalization audit are in
`notes/rsp3_covariance_closure.md`.

## Saved stellar evidence

- `analyze.py`: supplied work2 trace and copied version-22 photos.
- `analyze_pulsation.py`: supplied pulsation trace and copied photos.
- `inputs/` and `pulsation/`: the input snapshots used by those analyses.

## Independent checks

Run a script with `/Users/owner/opt/anaconda3/bin/python SCRIPT.py` from this
directory or use its full path. NumPy, SciPy and SymPy are required.

| Script | Result |
| --- | --- |
| `check_closure.py` | Original ODE counterexample, candidate local covariance and implicit blocks |
| `check_symbolic_closure.py` | Exact covariance, pressure-work and LNA inertia identities |
| `check_discrete_closure.py` | Pressure work, local nonlinear branches, zero moments and common transport |
| `check_residual_derivatives.py` | Hand-derived residual derivatives against complex-step differentiation |
| `check_transport_time_centering.py` | Centered-diffusion counterexample and positivity condition |
| `check_energy_pairing.py` | Heat-flux time weighting, face/cell energy projection and viscous work |
| `check_moving_boundary.py` | 48 nonlinear local/nonlocal tests with a prescribed moving convective boundary |

For `check_moving_boundary.py`, `OPENBLAS_NUM_THREADS=1` avoids threading
overhead in its small matrix solves. The reference solver uses log energy
internally; this is not a proposal to replace MESA's w variable.

JSON reports retain failed generic Newton and timestep-continuation attempts.
They are evidence of the solver issue, not passing solutions. The final
`moving_boundary_rate_checks.json` reports the positive predictor and rate
normalization: all 48 cases pass, maximum scaled rate residual 3.93e-13.
The physical residual, starting model and requested timestep remain fixed
during the predictor. The predictor is an iterative algebraic solve.

The tests include frozen gas coefficients or prescribed stratification.
They do not validate the coupled stellar thermal/mechanical solution, native
MESA derivatives, physical closure calibration, or stellar evolution.
Enhanced stable-layer dissipation has not been added by these scripts.
