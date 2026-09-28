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
