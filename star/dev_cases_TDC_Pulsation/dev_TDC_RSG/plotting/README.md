# Star LNA Kick Diagnostics

Run from the test case directory:

```sh
conda run -n base python plotting/plot_star_lna_kick.py
```

The script reads the existing `LOGS` and `LNA` output and writes:

- `kick_history.png`: timestep, surface velocity, and retry history;
- `kick_profile.png`: imposed velocity shape and profile evolution;
- `mode_quality.png`: displacement nodes and the kinematic eigenvector check;
- `kick_diagnostics.txt`: numerical values used in the figures.

It does not run or compile MESA.

# Energy-error diagnostics

Run from the test case directory:

```sh
conda run -n base python plotting/analyze_energy_error.py 301 314
```

The script compares the global history energy error with the sum of the saved
cell energy-equation residuals and checks whether the Lagrangian `dq` grid
changed between the selected profiles. For an excised-envelope model it uses
`dm = dq * (Mstar - M_center)`, matching MESA's `xmstar` normalization rather
than treating `dq` as a fraction of the complete stellar mass.

# Local algebraic anisotropy diagnostic

Run from the test case directory:

```sh
conda run -n base python plotting/analyze_local_anisotropy.py
```

The script evaluates the local anisotropy prescription of Kovacs, Szabo &
Nuspl (2026) on every existing `LOGS/profile*.data`.  It sets the paper's
nonlocal turbulent-energy diffusion term `D_omega` to zero, but uses the
instantaneous TDC amplitude, buoyant source, dissipation, and velocity
gradients from each saved model.  Thus the diagnostic is algebraic in the
current time-dependent state; it does not evolve an additional anisotropy
variable.

It writes:

- `local_anisotropy_profiles.png`: the radial closure, stress, and work in the
  latest profile;
- `local_anisotropy_evolution.png`: variation over all saved profiles;
- `local_anisotropy_fixed_epoch.png`: the ten profiles at
  `TDC_alpha_M = 0.5` immediately before the control change;
- `local_anisotropy_summary.csv`: one numerical summary row per profile;
- `local_anisotropy_diagnostics.txt`: assumptions and latest-profile values.

The `0 <= xi <= 1` clip in the post-processing is only a diagnostic.  A
coupled implementation would need a smooth realizability treatment.  The
profiled `Eq` is used as the authoritative current MESA eddy-viscous work;
the plotted continuum eddy stress is only a reconstruction because the saved
velocity does not retain the exact time-centered discrete strain.

This particular `LOGS` changed `TDC_alpha_M` during the run. The script uses
0.5 for models below 51000, 1.0 for models 51000--52000, and 0.7 beginning at
model 53000 when reconstructing the comparison eddy stress. This epoch table
does not affect the anisotropy closure or the saved `Eq` column.

The script does not run, compile, or modify MESA.
